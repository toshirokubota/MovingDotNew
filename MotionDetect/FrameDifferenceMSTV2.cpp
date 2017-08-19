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
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <CoreParticle.h>
#include <CoreParticleMakeGraph.h>
#include <Kruskal.h>
#include <MovingDot.h>
//using namespace cv;
using namespace std;

struct CoreBranch
{
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
		cv::Size sz(2,2);
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
	}
	float similarityMeasure(CoreBranch& br)
	{
		float d = length(center[0] - br.center[0], center[1] - br.center[1], 0, 0);
		return 1.0 / (1.0 + d); 
	}

	vector<TK::MovingDot*> dots;
	vector<float> center;
	float laxis;
	float saxis;
	float theta;
};

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
StdVector2Mat(vector<T>& F, cv::Mat& m, int ndim, int dims[])
{
	m.resize(dims[0] * dims[1]);
	m.reshape(dims[0], dims[1]);
	if (m.type() == CV_8U)
	{
		for (int i = 0; i < F.size(); ++i)
		{
			m.at<unsigned char>(i) = F[i];
		}
	}
	else if (m.type() == CV_32F)
	{
		for (int i = 0; i < F.size(); ++i)
		{
			m.at<float>(i) = F[i];
		}
	}
	else if (m.type() == CV_32S)
	{
		for (int i = 0; i < F.size(); ++i)
		{
			m.at<int>(i) = F[i];
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
vector<CoreBranch>
cutTrees(vector<TK::MovingDot*>& dots, float thres)
{
	vector<TK::Node<int>*> nodes;
	map<TK::Vertex<CoreParticle*>*, int> vmap;
	for (int i = 0; i < dots.size(); ++i)
	{
		nodes.push_back(TK::makeset(i));
		vmap[dots[i]->p->vertex] = i;
	}

	for (int i = 0; i < dots.size(); ++i)
	{
		TK::Vertex<CoreParticle*>* u = dots[i]->p->vertex;
		CoreParticle* p = u->key;
		for (int j = 0; j < u->aList.size(); ++j)
		{
			TK::Vertex<CoreParticle*>* v = u->aList[j]->v;
			CoreParticle* q = v->key;
			float len = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
			if (len < thres)
			{
				TK::merge(nodes[i], nodes[vmap[v]]);
			}
		}
	}
	vector<TK::Node<int>*> rep = TK::clusters(nodes);
	map<TK::Node<int>*, int> rmap;
	for (int i = 0; i < rep.size(); ++i)
	{
		rmap[rep[i]] = i;
	}
	vector<vector<TK::MovingDot*>> grouped(rep.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = rmap[findset(nodes[i])];
		grouped[k].push_back(dots[nodes[i]->key]);
	}
	vector<CoreBranch> branches;
	for (int i = 0; i < grouped.size(); ++i)
	{
		//if (grouped[i].size() >= 20)
		{
			branches.push_back(CoreBranch(grouped[i]));
		}
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}

	return branches;
}

void
estimateMotions(vector<CoreBranch>& branch0, vector<vector<CoreBranch>>& branches, 
			float thres,
			int ndim, const int* dims)
{
	for (int i = 0; i < branch0.size(); ++i)
	{
		for (int j = 0; j < branch0[i].dots.size(); ++j)
		{
			branch0[i].dots[j]->match = NULL;
		}

		float maxs = 0;
		int idx = -1;
		int jx = (int)branches.size() - 1;
		while (jx >= Max(0, (int)branches.size() - 5))
		{
			for (int j = 0; j < branches[jx].size(); ++j)
			{
				float s = branch0[i].similarityMeasure(branches[jx][j]);
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

		float dx = branch0[i].center[0] - branches[jx][idx].center[0];
		float dy = branch0[i].center[1] - branches[jx][idx].center[1];
		for (int j = 0; j < branch0[i].dots.size(); ++j)
		{
			float mind = std::numeric_limits<float>::infinity();
			TK::MovingDot* dot = branch0[i].dots[j];
			for (int k = 0; k < branches[jx][idx].dots.size(); k++)
			{
				TK::MovingDot* dot2 = branches[jx][idx].dots[k];
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
	}
}

float
weight(TK::MovingDot* p, TK::MovingDot* q)
{
	float d = sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) +
		(p->dx - q->dx)*(p->dx - q->dx) + (p->dy - q->dy)*(p->dy - q->dy));
	return d;
}

void
pruneEdges(vector<TK::Vertex<CoreParticle*>*>& vertices,
	vector<TK::Edge<CoreParticle*>*>& mst,
	float thres)
{
	for (int i = mst.size() - 1; i >= 0; --i)
	{
		if (mst[i]->w > thres)
		{
			mst.erase(mst.begin() + i);
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: P = MotionDetect(video_file)");
		return;
	}
	cv::VideoCapture cap;
	char ch[256];
	int len = mxGetString(prhs[0], ch, 256);
	cap.open(ch);

	//cap.open(0); //open a local camera (web-cam)
	if (!cap.isOpened())  // if not success, exit program
	{
		mexErrMsgTxt("Cannot open the video.");
		return;
	}
	int gray_thres = 30;
	float decRate = 1;
	int minSize = 3;
	float weight_thres = 10.0f;
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
		ReadScalar(decRate, prhs[2], classMode);
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[3], classMode);
	}
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(weight_thres, prhs[4], classMode);
	}	
	if (nrhs >= 6)
	{
		mxClassID classMode;
		ReadScalar(sleepTime, prhs[5], classMode);
	}

	double dWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH); //get the width of frames of the video
	double dHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT); //get the height of frames of the video
	//cout << "Frame size : " << dWidth << " x " << dHeight << endl;
	double fps = cap.get(CV_CAP_PROP_FPS); //get the frames per seconds of the video

	cout << "Frame per seconds : " << fps << endl;

	if (bDisplay)
	{
		cv::namedWindow("input", CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
		cv::namedWindow("motion", CV_WINDOW_AUTOSIZE);
		cv::namedWindow("skeleton", CV_WINDOW_AUTOSIZE);
	}

	cv::Mat frame; //after resizing frame
	cv::Mat frame0; //original frame
	cv::Size destSize((int)(dWidth/decRate), int(dHeight/decRate));
	cv::Mat prevFrame = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //previous frame
	int frameCount = 0;
	vector<vector<TK::MovingDot*>> cores;
	vector<vector<TK::Edge<CoreParticle*>*>> mst;
	vector<vector<CoreBranch>> branches;
	vector<TK::MovingDot*> dots;
	vector<cv::Mat> frames;

	int nframes = 2; //this many previous frames are considered for MST
	float tscale = 1.0f; //temporal unit copmared to the pixel size.
	float dist_thres = 2.0f; //a stringent threshold to cut tree into small branches.
	while (true)
	{
		bool bSuccess = cap.read(frame0); // read a new frame from video
		frameCount++;

		if (!bSuccess) //if not success, break loop
		{
			cout << "Cannot read a frame from video stream" << endl;
			break;
		}
		cv::resize(frame0, frame, destSize, 0.0, 0.0, cv::INTER_AREA);

		cv::Mat gray;
		cv::cvtColor(frame, gray, CV_BGR2GRAY);

		if (frameCount <= 1)
		{
			prevFrame = gray;
			continue; //skip the first frame since there is no previous frame, yet
		}
		frames.push_back(gray);

		cv::Mat diff = gray - prevFrame;
		prevFrame = gray;
		cv::Mat fg = cv::Mat::zeros(destSize, cv::DataType<int>::type); //blank frame
		cv::threshold(abs(diff), fg, (double)gray_thres, 255.0, 0);

		cv::Mat label, stats, centroids;
		int nc = connectedComponentsWithStats(fg, label, stats, centroids);
		for (int i = 0; i < fg.rows; ++i)
		{
			unsigned char* ptr = fg.ptr(i);
			for (int j = 0; j < fg.cols; ++j)
			{
				int lb = label.at<int>(i, j);
				if (lb > 0 && stats.at<int>(lb, 4) < minSize)
				{
					ptr[j] = 0;
				}
			}
		}

		int ndim;
		int dims[2];
		vector<unsigned char> L;
		MAT2StdVector(L, fg, ndim, dims);

		int nvoxels = numberOfElements(ndim, dims);
		vector<CoreParticle*> mp = generateParticleMap(L, ndim, dims);
		vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndim, dims);
		for (int i = 0; i < particles.size(); ++i)
		{
			particles[i]->t = frameCount; //used in constructing MST
		}
		vector<int> S(nvoxels, 0);
		propagateParticles(particles, mp, S, ndim, dims);
		vector<int> S2(nvoxels, 0);
		vector<TK::Vertex<CoreParticle*>*> vertices = makeGraphStructure(mp, S2, ndim, dims);
		vector<TK::MovingDot*> cores0;
		set<CoreParticle*> noncore;
		map<CoreParticle*, TK::MovingDot*> pmap;
		for (int i = 0; i < particles.size(); ++i)
		{
			if (particles[i]->vertex != NULL)
			{
				TK::MovingDot* dot = new TK::MovingDot(particles[i]);
				cores0.push_back(dot);
				pmap[particles[i]] = dot;
			}
			else
			{
				noncore.insert(particles[i]);
			}
		}
		dots.insert(dots.end(), cores0.begin(), cores0.end());
		cores.push_back(cores0);
		vector<TK::Edge<CoreParticle*>*> mst0;
		vector<CoreBranch> branch0 = cutTrees(cores0, dist_thres);
		if (branches.empty())
		{
			mst0 = Kruskal(vertices);
			mst.push_back(mst0);
		}
		else
		{
			estimateMotions(branch0, branches, decRate / 20.0f, ndim, dims);
			for (int i = 0; i < vertices.size(); ++i)
			{
				TK::Vertex<CoreParticle*>* u = vertices[i];
				TK::MovingDot* up = pmap[u->key];
				for (int j = 0; j < u->aList.size(); ++j)
				{
					TK::Vertex<CoreParticle*>* v = u->aList[j]->v;
					TK::MovingDot* vp = pmap[v->key];
					u->aList[j]->w = weight(up, vp);
				}
			}
			mst0 = Kruskal(vertices);
		}

		pruneEdges(vertices, mst0, weight_thres);
		mst.push_back(mst0);
		branches.push_back(branch0);
		if (bDisplay)
		{
			cv::Mat blank = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //blank frame
			if (mst.size() >= 3)
			{
				drawEdges(blank, mst, Max(0, mst.size() - 3), mst.size() - 1);
			}
			cv::Mat edg = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //blank frame
			//cv::Canny(frame, edg, 100, 200);
			drawMovingDots(edg, cores0);

			cv::imshow("input", frame); //show the frame in "MyVideo" window
			cv::imshow("motion", edg);
			cv::imshow("skeleton", blank);
		}
		if (cv::waitKey(sleepTime) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
			cout << "esc key is pressed by user" << endl;
			break;
		}
		CoreParticleFactory::getInstance().deleteParticle(noncore);
	}
	cap.release(); 

	cv::destroyAllWindows();
	if (nlhs >= 1)
	{
		const int dims[] = { frame.cols, frame.rows, frames.size() };
		vector<unsigned char> F(dims[0] * dims[1] * dims[2], (unsigned char)0);
		for (int i = 0; i < frames.size(); ++i)
		{
			for (int j = 0; j < frame.rows; ++j)
			{
				unsigned char* ptr = frames[i].ptr(j);
				for (int k = 0; k < frame.cols; ++k)
				{
					SetData3(F, k, j, i, dims[0], dims[1], dims[2], ptr[k]);
				}
			}
		}
		plhs[0] = StoreData(F, mxUINT8_CLASS, 3, dims);
	}

	TK::GraphFactory<CoreParticle*>::GetInstance().Clean();
	CoreParticleFactory::getInstance().clean();
	for (int i = 0; i < dots.size(); ++i)
	{
		delete dots[i];
	}
	mexUnlock();
}

