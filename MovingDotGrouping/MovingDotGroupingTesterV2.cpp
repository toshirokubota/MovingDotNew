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
		ascendant = NULL;
		id = _id++;
	}
	float intraSimilarityMeasure(const CoreBranch* br) const
	{
		float alpha = 0.5f;
		float sd = 0;
		float minMin = std::numeric_limits<float>::infinity();
		for (int i = 0; i < dots.size(); ++i)
		{
			TK::MovingDot* p = dots[i];
			TK::MovingDot* sel = NULL;
			float mind = std::numeric_limits<float>::infinity();
			for (int j = 0; j < br->dots.size(); ++j)
			{
				TK::MovingDot* q = br->dots[j];
				float d = sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y));
				if (d < mind)
				{
					mind = d;
					sel = q;
				}
			}
			sd += mind*mind;
			minMin = Min(minMin, mind);
		}
		sd = (sd / dots.size() + minMin*minMin) / 2.0; //average min-difference plus min-min
		float dvx = velocity[0] - br->velocity[0];
		float dvy = velocity[1] - br->velocity[1];
		float d = sqrt(alpha*(sd) + 100*(1.0 - alpha)*(dvx*dvx + dvy*dvy));
		return 1.0 / (1.0 + d);
	}
	float interSimilarityMeasure(const CoreBranch* br) const
	{
		float alpha = 0.5f;
		float svx = 0, svy = 0, sd=0;
		for (int i = 0; i < dots.size(); ++i)
		{
			TK::MovingDot* p = dots[i];
			TK::MovingDot* sel = NULL;
			float mind = std::numeric_limits<float>::infinity();
			for (int j = 0; j < br->dots.size(); ++j)
			{
				TK::MovingDot* q = br->dots[j];
				float d = sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y));
				if (d < mind)
				{
					mind = d;
					sel = q;
				}
			}
			svx += p->x - sel->x;
			svy += p->y - sel->y;
			sd += mind*mind;
		}
		float vx = svx / br->dots.size();
		float vy = svy / br->dots.size();
		float dvx = vx - br->velocity[0];
		float dvy = vy - br->velocity[1];
		float d = sqrt(alpha*(sd/dots.size()) + (1.0-alpha)*(dvx*dvx + dvy*dvy));
		return 1.0 / (1.0 + d);
	}

	void estimateVelocity() {
		if (dots.empty()) return;

		float sx = 0, sy = 0;
		for (int i = 0; i < dots.size(); ++i)
		{
			sx += dots[i]->dx;
			sy += dots[i]->dy;
		}
		velocity[0] = sx / dots.size();
		velocity[1] = sy / dots.size();
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
	CoreBranch* ascendant;
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

/*
For each new branch, find the best correspondence in the previous frames (going back to
at most 5 frames).
Using the corresponding branches, find matching pairs of moving dots and use them to 
estimate the velocity.
Finally, refine the velocity estimate of the new branch using dot pairs.
*/
void
estimateMotions(vector<CoreBranch*>& branch0, vector<vector<CoreBranch*>>& branches,
	float thres,
	int ndim, const int* dims)
{
	for (int i = 0; i < branch0.size(); ++i)
	{
		CoreBranch* br = branch0[i];
		for (int j = 0; j < br->dots.size(); ++j)
		{
			br->dots[j]->match = NULL;
		}

		float maxs = 0;
		int idx = -1;
		int jx = (int)branches.size() - 1;
		while (jx >= Max(0, (int)branches.size() - 5))
		{
			for (int j = 0; j < branches[jx].size(); ++j)
			{
				float s = branch0[i]->interSimilarityMeasure(branches[jx][j]);
				if (s > maxs)
				{
					maxs = s;
					idx = j;
				}
			}
			if (maxs >= thres) break; //found a good match
			jx--;
		}

		if (maxs < thres) continue; //not good enough
		br->ascendant = branches[jx][idx];
		br->ascendant->decendants.insert(br);

		float dx = br->center[0] - br->ascendant->center[0];
		float dy = br->center[1] - br->ascendant->center[1];
		for (int j = 0; j < br->dots.size(); ++j)
		{
			float mind = std::numeric_limits<float>::infinity();
			TK::MovingDot* dot = br->dots[j];
			for (int k = 0; k < br->ascendant->dots.size(); k++)
			{
				TK::MovingDot* dot2 = br->ascendant->dots[k];
				float dx2 = dot->x - dot2->x;
				float dy2 = dot->y - dot2->y;
				float d = sqrt(dx2*dx2 + dy2*dy2);
				if (d < mind)
				{
					dot->match = dot2;
					dot->dx = dot->x - dot2->x;
					dot->dy = dot->y - dot2->y;
					mind = d;
				}
			}
		}
		br->estimateVelocity();
	}
}

#include <SimilarityMatrix.h>

/*
Merge branches if they are close enough in both spatial and velocity.
Use centroid based hierarchical clustering.
*/
/*vector<CoreBranch*>
mergeCoreBranches0(vector<CoreBranch*> branches, float thres, float rate, int tm)
{
	//int tm = branches[0]->dots[0]->t;
	SimilarityMatrix sm(branches.size());
	for (int i = 0; i < branches.size(); ++i)
	{
		for (int j = i + 1; j < branches.size(); ++j)
		{
			float sval = branches[i]->similarityMeasure(branches[j]);
			sm.set(i, j, sval);
		}
	}
	if (tm == 49)
	{
		for (int i = 0; i < sm.size(); ++i)
		{
			for (int j = 0; j < sm.size(); ++j)
			{
				printf("%3.3f ", sm.m.at<float>(i, j));
			}
			printf("\n");
		}
	}
	for (int iter = 0; iter < 5; ++iter)
	{
		bool bchanged = false;
		while (true)
		{
			int row, col;
			float maxval = sm.maxVal(row, col);
			if (maxval > thres)
			{
				//printf("merging %d and %d\n", row, col);
				vector<TK::MovingDot*> dots = branches[row]->dots;
				dots.insert(dots.end(), branches[col]->dots.begin(), branches[col]->dots.end());
				delete branches[row];
				delete branches[col];
				branches[row] = new CoreBranch(dots);
				branches.erase(branches.begin() + col);
				sm.remove(col);
				bchanged = true;

				for (int i = 0; i < sm.size(); ++i)
				{
					if (i != row)
					{
						float sval = branches[row]->similarityMeasure(branches[i]);
						sm.set(row, i, sval);
					}
				}
			}
			else
			{
				break;
			}
		}
		if (!bchanged) break;
		thres *= rate;
	}
	if (tm == 49)
	{
		for (int i = 0; i < sm.size(); ++i)
		{
			for (int j = 0; j < sm.size(); ++j)
			{
				printf("%3.3f ", sm.m.at<float>(i, j));
			}
			printf("\n");
		}
	}
	return branches;
}*/

vector<CoreBranch*>
mergeCoreBranches(vector<CoreBranch*> branches, float thres, float rate, int tm)
{
	for (int iter = 0; iter < 5; ++iter)
	{
		bool bchanged = false;
		while (true)
		{
			float maxval = 0;
			pair<CoreBranch*, CoreBranch*> pr(NULL, NULL);
			int row, col;
			for (int i = 0; i < branches.size(); ++i)
			{
				for (int j = 0; j < branches.size(); ++j)
				{
					if (i == j) continue;
					float sval = branches[i]->intraSimilarityMeasure(branches[j]);
					if (sval > maxval)
					{
						pr.first = branches[i];
						pr.second = branches[j];
						row = i;
						col = j;
						maxval = sval;
					}
				}
			}
			if (maxval > thres)
			{
				//printf("merging %d and %d\n", row, col);
				vector<TK::MovingDot*> dots = pr.first->dots;
				dots.insert(dots.end(), pr.second->dots.begin(), pr.second->dots.end());
				delete branches[row];
				delete branches[col];
				branches[row] = new CoreBranch(dots);
				branches.erase(branches.begin() + col);
				bchanged = true;
			}
			else
			{
				break;
			}
		}
		if (!bchanged) break;
		thres *= rate;
	}
	if (tm == 49)
	{
		for (int i = 0; i < branches.size(); ++i)
		{
			for (int j = 0; j < branches.size(); ++j)
			{
				//float sval = (branches[i]->similarityMeasure(branches[j]) + branches[j]->similarityMeasure(branches[i])) / 2.0;
				float sval = branches[i]->intraSimilarityMeasure(branches[j]);
				printf("%3.3f ", sval);
			}
			printf("\n");
		}
	}
	return branches;
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


	if (bDisplay)
	{
		cv::namedWindow("input", CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
	}

	int end_frame = 0;
	for (int i = 0; i < dims[0]; ++i)
	{
		int k = (int)GetData2(P, i, 2, dims[0], dims[1], 0.0f);
		if (end_frame < k) end_frame = k;
	}
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<vector<CoreParticle*>> particles(end_frame+1);
	for (int i = 0; i < dims[0]; ++i)
	{
		float x = GetData2(P, i, 0, dims[0], dims[1], 0.0f);
		float y = GetData2(P, i, 1, dims[0], dims[1], 0.0f);
		float t = GetData2(P, i, 2, dims[0], dims[1], 0.0f);
		particles[(int)t].push_back(factory.makeParticle(x, y, 0, t));
	}
	
	vector<vector<TK::MovingDot*>> cores;
	vector<vector<CoreBranch*>> branches;
	vector<TK::MovingDot*> dots;
	cv::Size destSize(256, 256);

	int lval = 0;
	float tscale = 1.0f; //temporal unit copmared to the pixel size.
	cv::Mat prevFrame;
	int npixels = dims[0] * dims[1];
	int label = 0;

	for(int frcount=0; frcount<particles.size(); ++frcount)
	{
		vector<TK::MovingDot*> cores0;
		set<CoreParticle*> noncore;
		map<CoreParticle*, TK::MovingDot*> pmap;
		for (int i = 0; i < particles[frcount].size(); ++i)
		{
			TK::MovingDot* dot = new TK::MovingDot(particles[frcount][i]);
			dot->t = particles[frcount][i]->t;
			cores0.push_back(dot);
			pmap[particles[frcount][i]] = dot;
		}
		dots.insert(dots.end(), cores0.begin(), cores0.end());
		cores.push_back(cores0);
		vector<CoreBranch*> branch0 = cutTrees(cores0, dist_thres);
		if (!branches.empty())
		{
			estimateMotions(branch0, branches, est_thres, 2, dims);
		}
		int numbr = branch0.size();
		map<TK::MovingDot*, TK::MovingDot*> idmap;
		for (int i = 0; i < cores0.size(); ++i)
		{
			idmap[cores0[i]] = cores0[i]->match;
		}
		branch0 = mergeCoreBranches(branch0, merge_thres, 0.9, frcount);
		for (int i = 0; i < cores0.size(); ++i)
		{
			if (idmap[cores0[i]] != cores0[i]->match)
			{
				printf("Inconsistent match: %d %d\n", frcount, cores0[i]->id);
			}
		}
		//printf("%d: # branches = %d, after merge %d\n", frcount, numbr, branch0.size());
		for (int i = 0; i < branch0.size(); ++i)
		{
			branch0[i]->label = lval++;
		}
		branches.push_back(branch0);
		if (bDisplay)
		{
			cv::Mat M = cv::Mat::zeros(destSize, cv::DataType<int>::type); //blank frame
			drawMovingDots(M, cores0);
			if (cv::waitKey(sleepTime) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
			{
				cout << "esc key is pressed by user" << endl;
				break;
			}
		}
		//if(frcount == 2) break;
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
				SetData2(F, m, 7, dimsF[0], dimsF[1], (float)(br->ascendant ? br->ascendant->id + 1 : 0));
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
