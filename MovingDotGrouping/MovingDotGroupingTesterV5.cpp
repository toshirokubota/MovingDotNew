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

TK::MovingDot*
closestPoint(float x, float y, vector<TK::MovingDot*>& dots, float& minD)
{
	minD = std::numeric_limits<float>::infinity();
	TK::MovingDot* closest = NULL;
	for (int i = 0; i < dots.size(); ++i)
	{
		float d = length(x - dots[i]->x, y - dots[i]->y, 0.0f, 0.0f);
		if (d < minD)
		{
			minD = d;
			closest = dots[i];
		}
	}
	return closest;
}

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

	/*
	This version compares the centroids of this branch and the ascendent branches.
	It is important that centroid of the ascendents is computed using only points matched to the
	points in the current branch. Else, the centroid displacement can be drastically incorrect when 
	only a fraction of points are matched.
	*/
	void estimateVelocity(float dist_thres) {
		if (dots.empty()) return;

		vector<TK::MovingDot*> ascDots;
		for (set<CoreBranch*>::iterator it = ascendants.begin(); it != ascendants.end(); ++it)
		{
			CoreBranch* br = *it;
			ascDots.insert(ascDots.end(), br->dots.begin(), br->dots.end());
		}
		set<TK::MovingDot*> matched;
		for (int i = 0; i < dots.size(); ++i)
		{
			float d;
			TK::MovingDot* p = closestPoint(dots[i]->x, dots[i]->y, ascDots, d);
			if (d < dist_thres)
			{
				matched.insert(p);
			}
		}
		float sx = 0, sy = 0;
		for (set<TK::MovingDot*>::iterator it = matched.begin(); it != matched.end(); ++it)
		{
			TK::MovingDot* p = *it;
			sx += p->x;
			sy += p->y;
		}

		if (matched.size() > 0)
		{
			velocity[0] = center[0] - sx / matched.size();
			velocity[1] = center[1] - sy / matched.size();
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

float matchValue(CoreBranch* ref, CoreBranch* branch, bool bBackward)
{
	float maxD = 0;
	float sumD = 0;
	for (int j = 0; j < ref->dots.size(); ++j)
	{
		TK::MovingDot* p = ref->dots[j];
		float d;
		if (bBackward)
		{
			closestPoint(p->x - ref->velocity[0], p->y - ref->velocity[1], branch->dots, d);
			//closestPoint(p->x, p->y, branch->dots, d);
		}
		else
		{
			closestPoint(p->x + ref->velocity[0], p->y + ref->velocity[1], branch->dots, d);
			//closestPoint(p->x, p->y, branch->dots, d);
		}
		/*if (ref->id == 4205 - 1)// && branch->id == 4266 - 1)
		{
			TK::MovingDot* q = closestPoint(p->x + ref->velocity[0], p->y + ref->velocity[1], branch->dots, d);
			printf("Closest [%d,%d] -> [%d,%d], %f\n",
				p->x, p->y, q->x, q->y, d);
		}*/
		if (d > maxD)
		{
			maxD = d;
		}
		sumD += d;
	}
	//return maxD;
	return sumD / ref->dots.size();
}
/*
return true if all points in REF are a part of BRANCH (within the margin of THRES).
If bCausal, the velocity of REF is reversed.
*/
bool
isMatch(CoreBranch* ref, CoreBranch* branch, float thres, bool bBackward)
{
	float val = matchValue(ref, branch, bBackward);
	/*if (ref->id == 4195 - 1 && branch->id == 4252 - 1)
	{
		printf("%d, %d: matchValue = %f, thres = %f\n",
			ref->id+1, branch->id+1, val, thres);
	}*/
	return val < thres;
}

/*
collect all branches in BRANCHES that are good match to REF -- all points in REF are 
a part of the match.
*/
vector<CoreBranch*>
findMatches(CoreBranch* ref, vector<CoreBranch*>& branches, float thres, bool bBackward)
{
	vector<CoreBranch*> matches;
	for (int i = 0; i < branches.size(); ++i)
	{
		CoreBranch* br = branches[i];
		if (isMatch(ref, br, thres, bBackward))
		{
			matches.push_back(br);
		}
	}
	return matches;
}

#include <Hungarian.h>

/*
Establish correspondence between a branch in the current frame and those in the previous frame.
The routine handles a case where multiple branches merge into one in the current frame.
When found, the merged branch is splitted into multiple.
*/
vector<CoreBranch*> 
establishCorresponedence(vector<CoreBranch*>& current,
	vector<CoreBranch*>& prev,
	float thres)
{
	float penaltyCost = 1.0e4;
	int n = Max(prev.size(), current.size());
	int dims[] = { n, n };
	vector<float> C(dims[0] * dims[1], penaltyCost);
	map<CoreBranch*, int> bmap;
	for (int i = 0; i < current.size(); ++i)
	{
		bmap[current[i]] = i;
	}
	for (int i = 0; i < prev.size(); ++i)
	{
		for (int j = 0; j < current.size(); ++j)
		{
			float val = matchValue(prev[i], current[j], false);
			/*if (current[j]->id == 4266 - 1)
			{
				printf("%d %d: %f %f\n", prev[i]->id+1, current[j]->id+1, val, thres);
			}*/
			//if (val < thres)
			{
				SetData2(C, i, j, dims[0], dims[1], val);
			}
		}
	}
	Hungarian hungarian(C, dims);
	hungarian.solve(false);
	for (int i = 0; i < current.size(); ++i)
	{
		/*if (current[i]->id == 4266 - 1)
		{
			for (int j = 0; j < prev.size(); ++j)
			{
				printf("cost %d - %d = %f %f\n", prev[j]->id+1, current[i]->id+1, 
					hungarian.getCost(j, i),
					GetData2(C, j, i, dims[0], dims[1], -1000.0f));
				int k = hungarian.getWorker(j);
				if (k >= 0 && k < current.size())
				{
					printf("assignment %d -> %d with %f %f\n", 
						prev[j]->id + 1, current[k]->id + 1, 
						hungarian.getCost(j, k),
						GetData2(C, j, k, dims[0], dims[1], -1000.0f));
				}
				else
				{
					printf("assignment %d -> NA with %f\n", 
						prev[j]->id + 1, hungarian.getCost(j, k));
				}
			}
		}*/
	}

	for (int i = 0; i < prev.size(); ++i)
	{
		int k = hungarian.getWorker(i);
		if (k >= 0 && k < current.size())
		{
			if (hungarian.getCost(i, k) < thres)
			{
				current[k]->ascendants.insert(prev[i]);
				prev[i]->decendants.insert(current[k]);
			}
		}
	}
	for (int i = 0; i < current.size(); ++i)
	{
		current[i]->estimateVelocity(thres);
	}
	return current;
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
	//return sqrt(d1 * d2);
	return d2;
}

vector<float>
centroid2(set<CoreBranch*>& S)
{
	vector<float> center(2, 0.0f);
	int count = 0;
	for (set<CoreBranch*>::iterator it = S.begin(); it != S.end(); ++it)
	{
		for (int j = 0; j < (*it)->dots.size(); ++j)
		{
			TK::MovingDot* p = (*it)->dots[j];
			center[0] += p->x;
			center[1] += p->y;
			count++;
		}
	}
	if (count > 0)
	{
		center[0] /= count;
		center[1] /= count;
	}
	return center;
}

void
TraceBranchTree(CoreBranch* b, int maxdepth)
{
	set<CoreBranch*> S;
	S.insert(b);
	int iter = 0;
	while (!S.empty() && iter < maxdepth)
	{
		set<CoreBranch*> S2;
		for (set<CoreBranch*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreBranch* br = *it;
			printf("Branch: %d\n", br->id);
			for (int i = 0; i < br->dots.size(); ++i)
			{
				printf("%d %d %d\n", br->dots[i]->x, br->dots[i]->y, br->dots[i]->t);
			}
			S2.insert(br->ascendants.begin(), br->ascendants.end());
		}
		S = S2;
		iter++;
	}
}

float differenceBranches(CoreBranch* a, CoreBranch* b, int ngen)
{
	set<CoreBranch*> S;
	set<CoreBranch*> T;
	S.insert(a);
	T.insert(b);
	vector<vector<float>> vS;
	vector<vector<float>> vT;
	for (int i = 0; i < ngen; ++i)
	{
		vector<float> vcS = centroid2(S);
		vector<float> vcT = centroid2(T);
		vS.push_back(vcS);
		vT.push_back(vcT);
		/*if (a->id == 5829-1 && b->id == 5835-1)
		{
			printf("Gen = %d\n", i);
			for (set<CoreBranch*>::iterator it = S.begin(); it != S.end(); ++it)
			{
				CoreBranch* br = *it;
				printf("S: (%3.3f, %3.3f) (%3.3f, %3.3f) %d\n", 
					br->center[0], br->center[1], vcS[0], vcS[1], br->dots.size());
			}
			set<CoreBranch*> T2;
			for (set<CoreBranch*>::iterator it = T.begin(); it != T.end(); ++it)
			{
				CoreBranch* br = *it;
				printf("T: (%3.3f, %3.3f) (%3.3f, %3.3f) %d\n",
					br->center[0], br->center[1], vcT[0], vcT[1], br->dots.size());
			}
		}*/

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
	float val = 0;
	float wgt = 1.0;
	float swgt = 0.0;
	float rate = 0.8;
	for (int i = 0; i < vS.size(); ++i)
	{
		float dx1 = vT[i][0] - vS[i][0];
		float dy1 = vT[i][1] - vS[i][1];
		float sx = 0, sy = 0;
		int count = 0;
		for (int j = 0; j < vS.size(); ++j)
		{
			if (i == j) continue;
			float dx = vT[j][0] - vS[j][0];
			float dy = vT[j][1] - vS[j][1];
			sx += dx;
			sy += dy;
			count++;
		}
		float dfx = dx1 - sx / count;
		float dfy = dy1 - sy / count;
		if (a->id == 4015 - 1 && b->id == 4021 - 1)
		{
			printf("%d %d: (%3.3f,%3.3f) - (%3.3f,%3.3f)=(%3.3f,%3.3f)\n",
				a->id+1, b->id+1, dx1, dy1, sx/count, sy/count, dfx, dfy);
		}
		val += wgt * (dfx*dfx + dfy*dfy);
		swgt += wgt;
		wgt *= rate;
	}
	if (a->id == 4015 - 1 && b->id == 4021 - 1)
	{
		for (int i = 0; i < vS.size(); ++i)
		{
			printf("Centroid %d %d: (%3.3f,%3.3f), (%3.3f,%3.3f)\n",
				a->id+1, b->id+1,
				vS[i][0], vS[i][1], vT[i][0], vT[i][1]);
		}
		printf("result: %d (%3.3f,%3.3f), %d (%3.3f,%3.3f), %f\n", 
			a->id+1, a->center[0], a->center[1],
			b->id+1, b->center[0], b->center[1], sqrt(val / swgt));
	}

	return sqrt(val / swgt);
}

void
labelBranches(vector<CoreBranch*> branches, float merge_thres)
{
	bool bChanged = false;
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
			if (branches[i]->id == 4015-1 && branches[j]->id == 4021-1)
			{
				TraceBranchTree(branches[i], 5);
				TraceBranchTree(branches[j], 5);
				printf("%d (%d) -- %d (%d) => %f\n",
				branches[i]->id, branches[i]->dots.size(),
				branches[j]->id, branches[j]->dots.size(), dval);
			}
			if (dval < merge_thres)
			{
				TK::merge(nodes[i], nodes[j]);
				bChanged = true;
			}
		}
	}
	vector<TK::Node<CoreBranch*>*> reps = TK::clusters(nodes);
	for (int i = 0; i < branches.size(); ++i)
	{
		int k = distance(reps.begin(), find(reps.begin(), reps.end(), TK::findset(nodes[i])));
		branches[i]->label = k + 1;
	}
	//clean up
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
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
		if (branches.size() > 1)
		{
			branch0 = establishCorresponedence(branch0, branches[branches.size()-1], est_thres);
			labelBranches(branch0, merge_thres);
		}
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
						SetData2(F, m, 4, dimsF[0], dimsF[1], branches[i][j]->label);
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
		int dimsF[2] = { count, 15 };
		vector<float> F(dimsF[0] * dimsF[1], 0);
		int maxNAsc = 0;
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
				SetData2(F, m, 6, dimsF[0], dimsF[1], (float)br->ascendants.size());
				SetData2(F, m, 7, dimsF[0], dimsF[1], (float)br->id + 1);
				SetData2(F, m, 8, dimsF[0], dimsF[1], (float)i);
				SetData2(F, m, 9, dimsF[0], dimsF[1], (float)br->label);
				float d;
				TK::MovingDot* p = closestPoint(br->center[0], br->center[1], br->dots, d);
				SetData2(F, m, 10, dimsF[0], dimsF[1], (float)p->x);
				SetData2(F, m, 11, dimsF[0], dimsF[1], (float)p->y);
				int cnt = 0;
				for (set<CoreBranch*>::iterator it = br->ascendants.begin(); it != br->ascendants.end() && cnt < 3; it++, cnt++)
				{
					SetData2(F, m, 12+cnt, dimsF[0], dimsF[1], (float)((*it)->id + 1));
				}
				maxNAsc = Max(maxNAsc, br->ascendants.size());
			}
		}
		printf("The # of ascendnets are at most %d\n", maxNAsc);
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
