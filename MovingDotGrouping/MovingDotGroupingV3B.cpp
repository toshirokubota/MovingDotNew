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
		center = centroid(particles);
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
		velocity[0] = velocity[1] = 0.0f;
		ascendant = NULL;
		id = _id++;
	}
	float similarityMeasure0(CoreBranch* br)
	{
		float alpha = 1.0f;
		float dx = center[0] - br->center[0];
		float dy = center[1] - br->center[1];
		float dvx = velocity[0] - br->velocity[0];
		float dvy = velocity[1] - br->velocity[1];
		float d = sqrt(dx*dx + dy*dy + alpha*(dvx*dvx + dvy*dvy));
		return 1.0 / (1.0 + d);
	}

	float similarityMeasure(CoreBranch* br)
	{
		int ngen = 5;
		set<CoreBranch*> S;
		S.insert(this);
		set<CoreBranch*> T;
		T.insert(br);
		float maxVal = 0.0f;
		for (int i = 0; i < ngen; ++i)
		{
			float minVal = std::numeric_limits<float>::infinity();
			for (set<CoreBranch*>::iterator it = S.begin(); it != S.end(); ++it)
			{
				for (set<CoreBranch*>::iterator jt = T.begin(); jt != T.end(); ++jt)
				{
					float val = (*it)->similarityMeasure0(*jt);
					if (val < minVal) minVal = val;
				}
			}
			if (minVal > maxVal) maxVal = minVal;
		}
		return maxVal;
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
	vector<float> center;
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
			if (len < thres)
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
				float s = branch0[i]->similarityMeasure(branches[jx][j]);
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
				float dx2 = dot->x - dx - dot2->x;
				float dy2 = dot->y - dy - dot2->y;
				float d = sqrt(dx2*dx2 + dy2*dy2);
				if (d < mind)
				{
					dot->match = dot2;
					dot->dx = dot->x - dot2->x;
					dot->dy = dot->y - dot2->y;
				}
			}
		}
		br->estimateVelocity();
	}
}

float
weight(TK::MovingDot* p, TK::MovingDot* q)
{
	float d = sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) +
		(p->dx - q->dx)*(p->dx - q->dx) + (p->dy - q->dy)*(p->dy - q->dy));
	return d;
}


float 
similarity(TK::MovingDot* p, TK::MovingDot* q)
{
	float val = 0;
	for (int gen = 0; gen < 3; ++gen)
	{
		if (p != NULL && q != NULL)
		{
			float dx = p->x - q->x;
			float dy = p->y - q->y;
			float dvx = p->dx - q->dx;
			float dvy = p->dy - q->dy;
			float d = sqrt(dx*dx + dy*dy + dvx*dvx + dvy*dvy);
			if (d > val)
			{
				val = d;
			}
			p = p->match;
			q = q->match;
		}
		else
		{
			val = std::numeric_limits<float>::infinity();
			break;
		}
	}
	return val;
}

float
similarity(CoreBranch* a, CoreBranch* b, float wgt = 1.0f)
{
	float mind = std::numeric_limits<float>::infinity();
	pair<TK::MovingDot*, TK::MovingDot*> selected(NULL, NULL);
	for (int i = 0; i < a->dots.size(); ++i)
	{
		for (int j = 0; j < b->dots.size(); ++j)
		{
			//float d = length(a.dots[i]->x - b.dots[j]->x, a.dots[i]->y - b.dots[j]->y, 0, 0);
			float d = similarity(a->dots[i], b->dots[j]);
			if (d < mind)
			{
				mind = d;
				selected.first = a->dots[i];
				selected.second = b->dots[j];
			}
		}
	}
	if (mind < std::numeric_limits<float>::infinity())
	{
		float dvec = length(selected.first->dx - selected.second->dx, selected.first->dy - selected.second->dy, 0.0f, 0.0f);
		return 1.0 / (1.0 + mind + wgt * dvec);
	}
	else
	{
		return 0.0f;
	}
}

#include <SimilarityMatrix.h>

vector<CoreBranch*>
mergeCoreBranches(vector<CoreBranch*> branches, float thres, float rate)
{
	SimilarityMatrix sm(branches.size());
	for (int i = 0; i < branches.size(); ++i)
	{
		for (int j = i + 1; j < branches.size(); ++j)
		{
			float sval = branches[i]->similarityMeasure(branches[j]);
			sm.set(i, j, sval);
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
				vector<TK::MovingDot*> dots = branches[row]->dots;
				dots.insert(dots.end(), branches[col]->dots.begin(), branches[col]->dots.end());
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
	return branches;
}

void
groupCoreBranches(vector<CoreBranch*>& branches, float thres, int& label)
{
	vector<TK::Node<int>*> nodes;
	for (int i = 0; i < branches.size(); ++i)
	{
		nodes.push_back(TK::makeset(i));
	}
	for (int i = 0; i < branches.size(); ++i)
	{
		for (int j = i + 1; j < branches.size(); ++j)
		{
			if (similarity(branches[i], branches[j], 5.0f) > thres)
			{
				TK::merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<TK::Node<int>*> reps = TK::clusters(nodes);
	map<TK::Node<int>*, int> rmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		rmap[reps[i]] = label+i;
	}
	for (int i = 0; i < branches.size(); ++i)
	{
		branches[i]->label = rmap[TK::findset(nodes[i])];
	}
	label += reps.size();
}

void
labelMovingDots(vector<TK::MovingDot*>& dots)
{
	vector<TK::Node<TK::MovingDot*>*> nodes(dots.size());
	map < TK::MovingDot*, TK::Node<TK::MovingDot*>*> nmap;
	for (int i = 0; i < dots.size(); ++i)
	{
		nodes[i] = TK::makeset(dots[i]);
		nmap[dots[i]] = nodes[i];
	}
	for (int i = 0; i < dots.size(); ++i)
	{
		if (dots[i]->match)
		{
			merge(nodes[i], nmap[dots[i]->match]);
		}
	}
	vector<TK::Node<TK::MovingDot*>*> reps = TK::clusters(nodes);
	map<TK::Node<TK::MovingDot*>*, int> rmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		rmap[reps[i]] = i + 1;
	}
	for (int i = 0; i < dots.size(); ++i)
	{
		dots[i]->label = rmap[findset(nodes[i])];
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
}

/*
The code is taken from:
http://felix.abecassis.me/2011/09/opencv-morphological-skeleton/
*/
cv::Mat
skeleton(cv::Mat& img)
{
	cv::Mat skel(img.size(), CV_8UC1, cv::Scalar(0));
	cv::Mat temp(img.size(), CV_8UC1);
	cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
	bool done;
	do
	{
		cv::morphologyEx(img, temp, cv::MORPH_OPEN, element);
		cv::bitwise_not(temp, temp);
		cv::bitwise_and(img, temp, temp);
		cv::bitwise_or(skel, temp, skel);
		cv::erode(img, img, element);

		double max;
		cv::minMaxLoc(img, 0, &max);
		done = (max == 0);
	} while (!done);

	return skel;
}

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
	vector<unsigned char> V;
	LoadData(V, prhs[0], classV, ndim, &dims);

	int nvoxels = numberOfElements(ndim, dims);

	int gray_thres = 30;
	float decRate = 1;
	int minSize = 3;
	float group_thres = 0.1f;
	bool bDisplay = true;
	int sleepTime = 30;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(gray_thres, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(group_thres, prhs[2], classMode);
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(sleepTime, prhs[3], classMode);
	}
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(bDisplay, prhs[4], classMode);
	}

	if (bDisplay)
	{
		cv::namedWindow("input", CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
	}

	vector<cv::Mat> detected;
	vector<vector<TK::MovingDot*>> cores;
	vector<vector<CoreBranch*>> branches;
	vector<TK::MovingDot*> dots;
	cv::Size destSize(dims[1], dims[0]);

	int lval = 0;
	float tscale = 1.0f; //temporal unit copmared to the pixel size.
	float dist_thres = 2.0f; //a stringent threshold to cut tree into small branches.
	cv::Mat prevFrame;
	int npixels = dims[0] * dims[1];
	int label = 0;

	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	for(int frcount=0; frcount<dims[2]; ++frcount)
	{
		cv::Mat gray(dims[1], dims[0], CV_8U);
		StdVector2Mat(V, gray, frcount, ndim, dims);
		if (frcount <= 1)
		{
			prevFrame = gray;
			continue; //skip the first frame since there is no previous frame, yet
		}
		cv::Mat diff = gray - prevFrame;
		prevFrame = gray;
		cv::Mat fg = cv::Mat::zeros(destSize, cv::DataType<int>::type); //blank frame
		cv::threshold(cv::abs(diff), fg, (double)gray_thres, 255.0, 0);

		cv::Mat label, stats, centroids;
		int nc = connectedComponentsWithStats(fg, label, stats, centroids);
		for (int i = 0; i < fg.rows; ++i)
		{
			unsigned char* ptr = fg.ptr(i);
			for (int j = 0; j < fg.cols; ++j)
			{
				int lb = label.at<int>(i, j);
				if (lb > 0 && stats.at<int>(lb, 4) <= minSize)
				{
					ptr[j] = 0;
				}
			}
		}
		cv::Mat skel = skeleton(fg.clone());

		detected.push_back(skel);

		vector<TK::MovingDot*> cores0;
		set<CoreParticle*> noncore;
		map<CoreParticle*, TK::MovingDot*> pmap;
		for (int i = 0; i < skel.rows; ++i)
		{
			unsigned char* ptr = skel.ptr(i);
			for (int j = 0; j < skel.cols; ++j)
			{
				if (ptr[j])
				{
					CoreParticle* p = factory.makeParticle(j, i);
					TK::MovingDot* dot = new TK::MovingDot(p);
					dot->t = frcount;
					cores0.push_back(dot);
					pmap[p] = dot;
				}
			}
		}
		dots.insert(dots.end(), cores0.begin(), cores0.end());
		cores.push_back(cores0);
		vector<CoreBranch*> branch0 = cutTrees(cores0, dist_thres);
		if (!branches.empty())
		{
			estimateMotions(branch0, branches, decRate / 20.0f, 2, dims);
		}
		branch0 = mergeCoreBranches(branch0, 0.01, 0.9);
		for (int i = 0; i < branch0.size(); ++i)
		{
			branch0[i]->label = lval++;
		}
		branches.push_back(branch0);
		if (bDisplay)
		{
			cv::imshow("input", skel); //show the frame in "MyVideo" window
			if (cv::waitKey(sleepTime) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
			{
				cout << "esc key is pressed by user" << endl;
				break;
			}
		}
		prevFrame = gray;

		CoreParticleFactory::getInstance().deleteParticle(noncore);
	}
	//labelMovingDots(dots);

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
		int dimsF[2] = { count, 4 };
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
						SetData2(F, m, 3, dimsF[0], dimsF[1], branches[i][j]->label);
					}
				}
			}
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dimsF);
	}
	
	TK::GraphFactory<CoreParticle*>::GetInstance().Clean();
	CoreParticleFactory::getInstance().clean();
	for (int i = 0; i < dots.size(); ++i)
	{
		delete dots[i];
	}
	mexUnlock();
}
