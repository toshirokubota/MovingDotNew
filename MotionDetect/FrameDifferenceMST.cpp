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
#include <MovingDot.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <CoreParticle.h>
#include <CoreParticleMakeGraph.h>
#include <szDistanceTransformNonIsotropic.h>
#include <Kruskal.h>

//using namespace cv;
using namespace std;

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

/*float weight2(CoreParticle* p, CoreParticle* q, float tscale)
{
	float len = length(p->src->x - q->src->x, p->src->y - q->src->y, p->src->z - q->src->z, tscale * (p->src->t - q->src->t));
	return len*len;
}

vector<Edge<CoreParticle*>*> 
incrementalMST(vector<Vertex<CoreParticle*>*>& cores0,
	vector<vector<Vertex<CoreParticle*>*>>& cores,
	int nframes, float tscale, float thres,
	int ndim, const int* dims)
{
	GraphFactory<CoreParticle*>& factory = GraphFactory<CoreParticle*>::GetInstance();

	vector<pair<float, Edge<CoreParticle*>*>> pairs;
	for (int i = 0; i < cores0.size(); ++i)
	{
		for (int j = 0; j < cores0[i]->aList.size(); ++j)
		{
			Edge<CoreParticle*>* ed = cores0[i]->aList[j];
			pair<float, Edge<CoreParticle*>*> pr(ed->w, ed);
			pairs.push_back(pr);
		}
	}
	for (int k = Max(cores.size() - nframes + 1, 0); k < cores.size(); ++k)
	{
		for (int i = 0; i < cores[k].size(); ++i)
		{
			Vertex<CoreParticle*>* v = cores[k][i];
			for (int j = 0; j < cores0.size(); ++j)
			{
				Vertex<CoreParticle*>* u = cores0[j];
				float w = weight2(u->key, v->key, tscale);
				if (w < thres)
				{
					Edge<CoreParticle*>* ed = factory.makeEdge(v, u, w);
					pair<float, Edge<CoreParticle*>*> pr(ed->w, ed);
					pairs.push_back(pr);
				}
			}
		}
	}
	vector<Edge<CoreParticle*>*> mst;
	set<Vertex<CoreParticle*>*> S;
}*/

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
	int thres = 30;
	float decRate = 1;
	int minSize = 3;
	bool bDisplay = true;
	int sleepTime = 30;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
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
		unsigned int val;
		mxClassID classMode;
		ReadScalar(val, prhs[4], classMode);
		bDisplay = val == 0 ? false : true;
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
	vector<vector<TK::Vertex<CoreParticle*>*>> cores;
	vector<vector<TK::Edge<CoreParticle*>*>> mst;
	int nframes = 2; //this many previous frames are considered for MST
	float tscale = 1.0f; //temporal unit copmared to the pixel size.
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

		cv::Mat diff = gray - prevFrame;
		prevFrame = gray;
		cv::Mat fg = cv::Mat::zeros(destSize, cv::DataType<int>::type); //blank frame
		cv::threshold(abs(diff), fg, (double)thres, 255.0, 0);

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
		vector<unsigned char> iL(L.size(), 0);
		for (int i = 0; i < L.size(); ++i)
		{
			iL[i] = L[i] ? 0 : 1;
		}

		int nvoxels = numberOfElements(ndim, dims);
		vector<CoreParticle*> mp = generateParticleMap(L, ndim, dims);
		vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndim, dims);
		for (int i = 0; i < particles.size(); ++i)
		{
			particles[i]->t = frameCount; //used in constructing MST
		}
		vector<float> D(L.size(), 0.0f);
		{
			vector<float> vs(ndim, 1.0f);
			DistanceTransformEuclidF(D, iL, vs, ndim, dims);
			for (int i = 0; i < D.size(); ++i)
			{
				if (mp[i] != NULL)
				{
					mp[i]->dval = D[i];
				}
			}
			//iL.clear(); //to save memory 
			D.clear(); //to save memory 
		}
		vector<int> S(nvoxels, 0);
		propagateParticles(particles, mp, S, ndim, dims);
		vector<int> S2(nvoxels, 0);
		vector<TK::Vertex<CoreParticle*>*> vertices = makeGraphStructure(mp, S2, ndim, dims);
		vector<CoreParticle*> cores0;
		set<CoreParticle*> noncore;
		for (int i = 0; i < particles.size(); ++i)
		{
			if (strongMedialParticle(particles[i], ndim, dims))
			{
				cores0.push_back(particles[i]);
			}
			else
			{
				noncore.insert(particles[i]);
			}
		}

		//vector<Edge<CoreParticle*>*> mst0 = incrementalMST(vertices, cores, nframes, tscale, ndim, dims);
		vector<TK::Edge<CoreParticle*>*> mst0 = Kruskal(vertices);
		cores.push_back(vertices);
		mst.push_back(mst0);

		cv::Mat blank = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //blank frame
		if (mst.size() >= 3)
		{
			drawEdges(blank, mst, Max(0, mst.size() - 3), mst.size() - 1);
		}
		cv::Mat edg = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //blank frame
		cv::Canny(frame, edg, 100, 200);
		if (bDisplay)
		{
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
	mexUnlock();
}

