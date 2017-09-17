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
#include <MovingEdgeDetector.h>
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
cv::Mat
StdVector2Mat(vector<T>& F, int ndim, const int* dims)
{
	//m.resize(dims[0] * dims[1]);
	//m.reshape(dims[0], dims[1]);
	cv::Size destSize(dims[0], dims[1]);
	cv::Mat m = cv::Mat::zeros(destSize, cv::DataType<T>::type);
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
	return m;
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
Update PREV by incorporating a new image in CUR.
*/
void
temporalSmoothing(cv::Mat& prev, cv::Mat cur)
{
	for (int i = 0; i < cur.rows; ++i)
	{
		unsigned char* ptr = cur.ptr<unsigned char>(i);
		unsigned char* ptr2 = prev.ptr<unsigned char>(i);
		for (int j = 0; j < cur.cols; ++j)
		{
			ptr2[j] = ptr2[j] / 2 + ptr[j] / 2;
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 1)
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
	float threshold = 0.5;
	float percentage = 0.7;
	float ratio = 0.7;
	int minSize = 5;
	bool bDisplay = true;
	float srate = 1.0f; //subsample rate
	int sleepTime = 30;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(threshold, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[2], classMode);
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(percentage, prhs[3], classMode);
	}
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(ratio, prhs[4], classMode);
	}
	if (nrhs >= 6)
	{
		mxClassID classMode;
		ReadScalar(srate, prhs[5], classMode);
	}
	if (nrhs >= 7)
	{
		mxClassID classMode;
		ReadScalar(sleepTime, prhs[6], classMode);
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
		cv::namedWindow("edge", CV_WINDOW_AUTOSIZE);
	}

	cv::Size destSize(dWidth/srate, dHeight/srate);
	cv::Mat smoothFrame; // = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //previous frame

	TK::MovingEdgeDetector detector(threshold, percentage, ratio, minSize);

	int frameCount = 0;
	vector<cv::Mat> frames;
	vector<cv::Mat> motions;
	vector<cv::Mat> edges;
	vector<vector<int>> points;
	while (true)
	{
		cv::Mat frame0; //original frame
		bool bSuccess = cap.read(frame0); // read a new frame from video
		frameCount++;

		if (!bSuccess) //if not success, break loop
		{
			cout << "Cannot read a frame from video stream" << endl;
			break;
		}
		cv::Mat frame; //after resizing frame
		resize(frame0, frame, destSize);

		cv::Mat gray;
		cv::cvtColor(frame, gray, CV_BGR2GRAY);

		//cv::Mat edge = cv::Mat::zeros(destSize, cv::DataType<unsigned char>::type); //blank frame
		int dims[] = { destSize.width, destSize.height };
		int ndim;
		vector<float> fim;
		MAT2StdVector(fim, gray, ndim, dims);
		detector.run(fim, ndim, dims);
		//cv::Canny(smoothFrame, edge, low, high); // , 7, true);
		vector<unsigned char> _edge = detector.retrieveEdges();
		cv::Mat edge = StdVector2Mat(_edge, ndim, dims);

		if (frameCount <= 1)
		{
			continue; //skip the first frame since there is no previous frame, yet
		}
		/*printf("%d edge retained out of %d.\n",
			detector.edges.size(),
			detector.edgeDetector[(detector.next + 1) % 2].edges.size());*/
		vector<unsigned char> _dif = detector.retrieveMovingEdges();
		cv::Mat dif = StdVector2Mat(_dif, ndim, dims);

		frames.push_back(gray);
		edges.push_back(edge);
		motions.push_back(dif);

		if (bDisplay)
		{
			imshow("input", gray); //show the frame in "MyVideo" window
			imshow("motion", dif);
			imshow("edge", edge);
		}
		if (cv::waitKey(sleepTime) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
			cout << "esc key is pressed by user" << endl;
			break;
		}
	}
	cap.release();

	if (nlhs >= 1)
	{
		const int dims[] = { points.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < points.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i][0]);
			SetData2(F, i, 1, dims[0], dims[1], points[i][1]);
			SetData2(F, i, 2, dims[0], dims[1], points[i][2]);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { (int)destSize.width, (int)destSize.height, motions.size() };
		vector<unsigned char> F(dims[0] * dims[1] * dims[2], 0);
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				unsigned char* ptr = motions[i].ptr<unsigned char>(j);
				for (int k = 0; k < dims[0]; ++k)
				{
					SetData3(F, k, j, i, dims[0], dims[1], dims[2], (unsigned char)ptr[k]);
				}
			}
		}
		plhs[1] = StoreData(F, mxUINT8_CLASS, 3, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { (int)destSize.width, (int)destSize.height, edges.size() };
		vector<unsigned char> F(dims[0] * dims[1] * dims[2], (unsigned char)0);
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				unsigned char* ptr = edges[i].ptr(j);
				for (int k = 0; k < dims[0]; ++k)
				{
					SetData3(F, k, j, i, dims[0], dims[1], dims[2], ptr[k]);
				}
			}
		}
		plhs[2] = StoreData(F, mxUINT8_CLASS, 3, dims);
	}
	if (nlhs >= 4)
	{
		const int dims[] = { (int)destSize.width, (int)destSize.height, frames.size() };
		vector<unsigned char> F(dims[0] * dims[1] * dims[2], (unsigned char)0);
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				unsigned char* ptr = frames[i].ptr(j);
				for (int k = 0; k < dims[0]; ++k)
				{
					SetData3(F, k, j, i, dims[0], dims[1], dims[2], ptr[k]);
				}
			}
		}
		plhs[3] = StoreData(F, mxUINT8_CLASS, 3, dims);
	}
	cv::destroyAllWindows();
	mexUnlock();
}
