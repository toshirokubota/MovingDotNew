// IPCameraCaptureOpenCV.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <mex.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/video/background_segm.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <CoreParticle.h>
#include <CoreParticleMakeGraph.h>
#include <szConnectedComponentC.h>
#include <MovingDot.h>
//using namespace cv;
using namespace std;

struct CoreBranch
{
	static int _id;
	CoreBranch(vector<TK::MovingDot*>& dots)
	{
		this->dots = dots;
		vector<CoreParticle*> particles;
		for (int i = 0; i < dots.size(); ++i)
		{
			particles.push_back(dots[i]->p);
		}
		vector<float> center0 = centroid(particles);
		center[0] = center0[0];
		center[1] = center0[1];
		float xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
		for (int i = 0; i < dots.size(); ++i)
		{
			float x0 = dots[i]->x - center[0];
			float y0 = dots[i]->y - center[1];
			xx += x0*x0;
			xy += x0*y0;
			yy += y0*y0;
		}
		cv::Size sz(2, 2);
		cv::Mat m = cv::Mat::zeros(sz, cv::DataType<float>::type);
		m.at<float>(0, 0) = xx / dots.size();
		m.at<float>(0, 1) = xy / dots.size();
		m.at<float>(1, 0) = xy / dots.size();
		m.at<float>(1, 1) = yy / dots.size();
		cv::Mat evals;
		cv::Mat evecs;
		cv::eigen(m, evals, evecs);
		laxis = evals.at<float>(0);
		saxis = evals.at<float>(1);
		theta = atan2(evecs.at<float>(0, 1), evecs.at<float>(0, 0));
		label = 0;
		float sdx = 0, sdy = 0;
		for (int i = 0; i < dots.size(); ++i)
		{
			sdx += dots[i]->dx;
			sdy += dots[i]->dy;
		}
		velocity[0] = sdx / dots.size(); 
		velocity[1] = sdy / dots.size();
		id = _id++;
	}
	void estimateVelocity() {
		if (dots.empty()) return;

		float sx = 0, sy = 0;
		int count = 0;
		for (set<CoreBranch*>::iterator it=ascendants.begin(); it != ascendants.end(); ++it)
		{
			CoreBranch* br = *it;
			sx += br->center[0] * br->dots.size();
			sy += br->center[1] * br->dots.size();
			count += br->dots.size();
		}
		if (count > 0)
		{
			velocity[0] = center[0] - sx / count;
			velocity[1] = center[1] - sy / count;
		}
		else
		{
			velocity[0] = 0;
			velocity[1] = 0;
		}
	}

	vector<TK::MovingDot*> dots;
	float center[2];
	float laxis;
	float saxis;
	float theta;
	int id;
	int label;
	float velocity[2];
	set<CoreBranch*> decendants;
	set<CoreBranch*> ascendants;
};

int CoreBranch::_id = 0;

/*
Convert a 2D Mat to a std::vector*/
template<class T>
void
MAT2StdVector(vector<T>& F, cv::Mat& m, int& ndim, int dims[])
{
	ndim = 2;
	dims[0] = m.cols, dims[1] = m.rows;

	F.resize(dims[0] * dims[1]);
	if (m.type() == CV_8U)
	{
		for (int j = 0; j < dims[0] * dims[1]; ++j)
		{
			F[j] = (T)m.at<unsigned char>(j);
		}
	}
	else if (m.type() == CV_32F)
	{
		for (int j = 0; j < dims[0] * dims[1]; ++j)
		{
			F[j] = (T)m.at<float>(j);
		}
	}
	else if (m.type() == CV_32S)
	{
		for (int j = 0; j < dims[0] * dims[1]; ++j)
		{
			F[j] = (T)m.at<int>(j);
		}
	}
	else {
		//un-supported format...
	}
}

template<class T>
void
StdVector2Mat(vector<T>& F, cv::Mat& m, int frame, int ndim, const int* dims)
{
	//m.resize(dims[0] * dims[1]);
	//m.reshape(dims[0], dims[1]);
	int ofs = frame * dims[0] * dims[1];
	if (m.type() == CV_8U)
	{
		for (int i = 0; i < m.rows * m.cols; ++i)
		{
			m.at<unsigned char>(i) = F[ofs + i];
		}
	}
	else if (m.type() == CV_32F)
	{
		for (int i = 0; i <  m.rows * m.cols; ++i)
		{
			m.at<float>(i) = F[ofs + i];
		}
	}
	else if (m.type() == CV_32S)
	{
		for (int i = 0; i <  m.rows * m.cols; ++i)
		{
			m.at<int>(i) = F[ofs + i];
		}
	}
	else {
		//un-supported format...
	}
}

void
drawEdges(cv::Mat& M,
	vector<vector<TK::Edge<CoreParticle*>*>>& mst, int f1, int f2)
{
	for (int f = f1; f < f2; ++f)
	{
		for (int i = 0; i < mst[f].size(); ++i)
		{
			CoreParticle* p = mst[f][i]->u->key;
			CoreParticle* q = mst[f][i]->v->key;
			//cv::line(M, cv::Point2i(p->x, p->y), cv::Point2i(q->x, q->y), cv::Scalar(128));
			cv::circle(M, cv::Point2i(p->x, p->y), 0, cv::Scalar(128));
		}
	}
	for (int i = 0; i < mst[f2].size(); ++i)
	{
		CoreParticle* p = mst[f2][i]->u->key;
		CoreParticle* q = mst[f2][i]->v->key;
		//cv::line(M, cv::Point2i(p->x, p->y), cv::Point2i(q->x, q->y), cv::Scalar(255));
		cv::circle(M, cv::Point2i(p->x, p->y), 0, cv::Scalar(255));
	}
}

void
drawMovingDots(cv::Mat& M,
	vector<TK::MovingDot*>& dots)
{
	for (int i = 0; i < dots.size(); ++i)
	{
		CoreParticle* p = dots[i]->p;
		cv::circle(M, cv::Point2i(p->x, p->y), 0, cv::Scalar(255));
		if (dots[i]->match != NULL)
		{
			CoreParticle* q = dots[i]->match->p;
			cv::circle(M, cv::Point2i(q->x, q->y), 0, cv::Scalar(128));
			cv::line(M, cv::Point2i(p->x, p->y), cv::Point2i(q->x, q->y), cv::Scalar(128));
		}
	}
}

/*
Break a graph into smaller branches.
A branch is used to estimate the motion at a semi-macro level.
*/
vector<CoreBranch*>
cutTrees(vector<TK::MovingDot*>& dots, float thres)
{
	vector<TK::Node<TK::MovingDot*>*> nodes;
	map<TK::MovingDot*, int> vmap;
	for (int i = 0; i < dots.size(); ++i)
	{
		nodes.push_back(TK::makeset(dots[i]));
		vmap[dots[i]] = i;
	}

	for (int i = 0; i < dots.size(); ++i)
	{
		TK::MovingDot* p = dots[i];
		for (int j = i+ 1; j < dots.size(); ++j)
		{
			TK::MovingDot* q = dots[j];
			float len = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
			if (len <= thres)
			{
				TK::merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<TK::Node<TK::MovingDot*>*> rep = TK::clusters(nodes);
	map<TK::Node<TK::MovingDot*>*, int> rmap;
	for (int i = 0; i < rep.size(); ++i)
	{
		rmap[rep[i]] = i;
	}
	vector<vector<TK::MovingDot*>> grouped(rep.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = rmap[findset(nodes[i])];
		grouped[k].push_back(nodes[i]->key);
	}
	vector<CoreBranch*> branches;
	for (int i = 0; i < grouped.size(); ++i)
	{
		//if (grouped[i].size() >= 20)
		{
			branches.push_back(new CoreBranch(grouped[i]));
		}
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}

	return branches;
}

TK::MovingDot*
closestPoint(float x, float y, vector<TK::MovingDot*>& dots, float& minD)
{
	minD = std::numeric_limits<float>::infinity();
	TK::MovingDot* closest = NULL;
	for (int i = 0; i < dots.size(); ++i)
	{
		float d = length(x - dots[i]->x, y- dots[i]->y, 0.0f, 0.0f);
		if (d < minD)
		{
			minD = d;
			closest = dots[i];
		}
	}
	return closest;
}

/*
Establish correspondence between a branch in the current frame and those in the previous frame.
The routine handles a case where multiple branches merge into one in the current frame.
When found, the merged branch is splitted into multiple.
*/
vector<CoreBranch*> 
establishCorresponedence(vector<CoreBranch*>& current,
	vector<vector<CoreBranch*>>& past,
	float thres)
{
	vector<CoreBranch*> result;
	if(!past.empty())
	{
		vector<CoreBranch*> prev = past[past.size() - 1];
		for (int i = 0; i < current.size(); ++i)
		{
			CoreBranch* br = current[i];
			vector<TK::MovingDot*> dots = br->dots;
			vector<CoreBranch*> matches;
			for (int j = 0; j < prev.size(); ++j)
			{
				CoreBranch* bs = prev[j];
				vector<TK::MovingDot*> dots2 = bs->dots;
				float maxD = 0;
				for (int k = 0; k < dots2.size(); ++k)
				{
					float mind;
					closestPoint(dots2[k]->x+bs->velocity[0], dots2[k]->y+bs->velocity[1], dots, mind);
					if (mind > maxD)
					{
						maxD = mind;
					}
				}
				if (maxD <= thres)
				{
					matches.push_back(bs);
				}
			}
			if (matches.size() == 0)
			{
				result.push_back(br);
			}
			else if (matches.size() == 1)
			{
				br->ascendants.insert(matches[0]);
				matches[0]->decendants.insert(br);
				result.push_back(br);
			}
			else //if (matches.size() > 1)
			{
				//for each dot in matched branch, find the closest one in the current branch.
				map<TK::MovingDot*, int> dmap;
				for (int j = 0; j < dots.size(); ++j)
				{
					float mind = std::numeric_limits<float>::infinity();
					int idx = -1;
					for (int k = 0; k < matches.size(); ++k)
					{
						float d;
						closestPoint(dots[j]->x-matches[k]->velocity[0], dots[j]->y - matches[k]->velocity[1], matches[k]->dots, d);
						if (d < mind)
						{
							mind = d;
							idx = k;
						}
					}
					dmap[dots[j]] = idx;
				}
				vector<vector<TK::MovingDot*>> splitted(matches.size());
				for (int j = 0; j < dots.size(); ++j)
				{
					splitted[dmap[dots[j]]].push_back(dots[j]);
				}
				for (int j = 0; j < matches.size(); ++j)
				{
					if (!splitted[j].empty())
					{
						CoreBranch* newbr = new CoreBranch(splitted[j]);
						newbr->ascendants.insert(matches[j]);
						matches[j]->decendants.insert(newbr);
						result.push_back(newbr);
					}
				}
			}
		}
	}
	else
	{
		result = current;
	}
	for (int i = 0; i < result.size(); ++i)
	{
		result[i]->estimateVelocity();
	}

	return result;
}

float
distance(CoreBranch* a, CoreBranch* b)
{
	float d1 = length(a->center[0] - b->center[0], a->center[1] - b->center[1], 0.0f, 0.0f);
	float d2;
	if (a->ascendants.empty() && b->ascendants.empty())
	{ //velocity is not available.
		d2 = 1.0;
	}
	else
	{
		d2 = length(a->velocity[0] - b->velocity[0], a->velocity[1] - b->velocity[1], 0.0f, 0.0f);
	}
	return sqrt(d1 * d2);
}

float differenceBranches(CoreBranch* a, CoreBranch* b, int ngen)
{
	float sumD = 0;
	int count = 0;
	set<CoreBranch*> S;
	set<CoreBranch*> T;
	S.insert(a);
	T.insert(b);
	for (int i = 0; i < ngen; ++i)
	{
		for (set<CoreBranch*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			for (set<CoreBranch*>::iterator jt = T.begin(); jt != T.end(); ++jt)
			{
				float dval = distance(*it, *jt);
				sumD += dval;
				count++;
			}
		}
		set<CoreBranch*> S2;
		for (set<CoreBranch*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			S2.insert((*it)->ascendants.begin(), (*it)->ascendants.end());
		}
		set<CoreBranch*> T2;
		for (set<CoreBranch*>::iterator it = T.begin(); it != T.end(); ++it)
		{
			T2.insert((*it)->ascendants.begin(), (*it)->ascendants.end());
		}
		S = S2;
		T = T2;
	}

	return sumD / count;
}

/*
Merge branches in the current frame if they are likely to be moving in the same direction.
*/
vector<CoreBranch*> 
mergeBranches(vector<CoreBranch*>& branches, float thres)
{
	vector<TK::Node<CoreBranch*>*> nodes;
	for (int i = 0; i < branches.size(); ++i)
	{
		nodes.push_back(TK::makeset(branches[i]));
	}
	for (int i = 0; i < branches.size(); ++i)
	{
		for (int j = i + 1; j < branches.size(); ++j)
		{
			float dval = differenceBranches(branches[i], branches[j], 5);
			if (dval < thres)
			{
				TK::merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<TK::Node<CoreBranch*>*> reps = TK::clusters(nodes);
	vector<vector<CoreBranch*>> mergedBranches(reps.size());
	for (int i = 0; i < branches.size(); ++i)
	{
		int k = distance(reps.begin(), find(reps.begin(), reps.end(), TK::findset(nodes[i])));
		mergedBranches[k].push_back(branches[i]);
	}
	vector<CoreBranch*> result(reps.size());
	for (int i = 0; i < mergedBranches.size(); ++i)
	{
		vector<TK::MovingDot*> vb;
		for (int j = 0; j < mergedBranches[i].size(); ++j)
		{
			vb.insert(vb.end(), mergedBranches[i][j]->dots.begin(), mergedBranches[i][j]->dots.end());
		}
		result[i] = new CoreBranch(vb);
		for (int j = 0; j < mergedBranches[i].size(); ++j)
		{
			CoreBranch* br = mergedBranches[i][j];
			for (set<CoreBranch*>::iterator it = br->ascendants.begin(); it != br->ascendants.end(); ++it)
			{
				result[i]->ascendants.insert(*it);
				(*it)->decendants.insert(result[i]);
				(*it)->decendants.erase(br);
			}
		}
	}
	for (int i = 0; i < result.size(); ++i)
	{
		result[i]->estimateVelocity();
	}

	//clean up
	for (int i = 0; i < branches.size(); ++i)
	{
		delete branches[i];
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}

	return result;
}

int TK::MovingDot::_id = 0;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: P = MotionDetect(video_file)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classV;
	vector<float> P;
	LoadData(P, prhs[0], classV, ndim, &dims);

	int nvoxels = numberOfElements(ndim, dims);

	float merge_thres = 0.01f;
	float est_thres = 0.05f;
	float dist_thres = 5.0f; //a stringent threshold to cut tree into small branches.
	bool bDisplay = true;
	int sleepTime = 30;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(merge_thres, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(est_thres, prhs[2], classMode);
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(dist_thres, prhs[3], classMode);
	}

	int end_frame = 0;
	for (int i = 0; i < dims[0]; ++i)
	{
		int k = (int)GetData2(P, i, 2, dims[0], dims[1], 0.0f);
		if (end_frame < k) end_frame = k;
	}
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<vector<CoreParticle*>> particles(end_frame+1);
	vector<vector<TK::MovingDot*>> cores(end_frame + 1);
	for (int i = 0; i < dims[0]; ++i)
	{
		float x = GetData2(P, i, 0, dims[0], dims[1], 0.0f);
		float y = GetData2(P, i, 1, dims[0], dims[1], 0.0f);
		float t = GetData2(P, i, 2, dims[0], dims[1], 0.0f);
		CoreParticle* p = factory.makeParticle(x, y, 0, t);
		particles[(int)t].push_back(factory.makeParticle(x, y, 0, t));
		cores[(int)t].push_back(new TK::MovingDot(p));
	}
	
	vector<vector<CoreBranch*>> branches;
	vector<TK::MovingDot*> dots;
	cv::Size destSize(256, 256);

	int lval = 0;
	float tscale = 1.0f; //temporal unit copmared to the pixel size.
	cv::Mat prevFrame;
	int npixels = dims[0] * dims[1];

	for(int frcount=0; frcount<particles.size(); ++frcount)
	{
		vector<TK::MovingDot*> cores0 = cores[frcount];
		set<CoreParticle*> noncore;
		dots.insert(dots.end(), cores0.begin(), cores0.end());
		cores.push_back(cores0);

		vector<CoreBranch*> branch0 = cutTrees(cores0, dist_thres);
		branch0 = establishCorresponedence(branch0, branches, dist_thres);
		branch0 = mergeBranches(branch0, merge_thres);
		branch0 = establishCorresponedence(branch0, branches, dist_thres);
		//branch0 = mergeBranches(branch0, merge_thres);
		//branch0 = establishCorresponedence(branch0, branches, dist_thres);
		branches.push_back(branch0);
	}

	cv::destroyAllWindows(); 
	if (nlhs >= 1)
	{
		int count = 0;
		for (int i = 0; i < branches.size(); ++i)
		{
			for (int j = 0; j < branches[i].size(); ++j)
			{
				count += branches[i][j]->dots.size();
			}
		}
		int dimsF[2] = { count, 7 };
		vector<int> F(dimsF[0] * dimsF[1], 0);
		for (int i = 0, m = 0; i < branches.size(); ++i)
		{
			for (int j = 0; j < branches[i].size(); ++j)
			{
				for (int k = 0; k < branches[i][j]->dots.size(); ++k, m++)
				{
					TK::MovingDot* p = branches[i][j]->dots[k];
					if (p)
					{
						SetData2(F, m, 0, dimsF[0], dimsF[1], p->x);
						SetData2(F, m, 1, dimsF[0], dimsF[1], p->y);
						SetData2(F, m, 2, dimsF[0], dimsF[1], p->t);
						SetData2(F, m, 3, dimsF[0], dimsF[1], p->id+1);
						if (p->match)
						{
							SetData2(F, m, 4, dimsF[0], dimsF[1], p->match->id+1);
						}
						SetData2(F, m, 5, dimsF[0], dimsF[1], branches[i][j]->id + 1);
						SetData2(F, m, 6, dimsF[0], dimsF[1], (int)branches[i][j]->decendants.size());
					}
				}
			}
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dimsF);
	}
	if (nlhs >= 2)
	{
		int count = 0;
		for (int i = 0; i < branches.size(); ++i)
		{
			count += branches[i].size();
		}
		int dimsF[2] = { count, 9 };
		vector<float> F(dimsF[0] * dimsF[1], 0);
		for (int i = 0, m = 0; i < branches.size(); ++i)
		{
			for (int j = 0; j < branches[i].size(); ++j, ++m)
			{
				CoreBranch* br = branches[i][j];
				SetData2(F, m, 0, dimsF[0], dimsF[1], br->center[0]);
				SetData2(F, m, 1, dimsF[0], dimsF[1], br->center[1]);
				SetData2(F, m, 2, dimsF[0], dimsF[1], br->velocity[0]);
				SetData2(F, m, 3, dimsF[0], dimsF[1], br->velocity[1]);
				SetData2(F, m, 4, dimsF[0], dimsF[1], (float)br->dots.size());
				SetData2(F, m, 5, dimsF[0], dimsF[1], (float)br->decendants.size());
				SetData2(F, m, 6, dimsF[0], dimsF[1], (float)br->id + 1);
				SetData2(F, m, 7, dimsF[0], dimsF[1], (float)(br->ascendants.size()));
				SetData2(F, m, 8, dimsF[0], dimsF[1], (float)i);
			}
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dimsF);
	}

	TK::GraphFactory<CoreParticle*>::GetInstance().Clean();
	CoreParticleFactory::getInstance().clean();
	for (int i = 0; i < dots.size(); ++i)
	{
		delete dots[i];
	}
	mexUnlock();
}
