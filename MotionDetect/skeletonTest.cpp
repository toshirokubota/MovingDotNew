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
StdVector2Mat(vector<T>& F, cv::Mat& m, int ndim, const int* dims)
{
	//m.resize(dims[0] * dims[1]);
	//m.reshape(dims[0], dims[1]);
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
	int thres = 30;
	int blockSize = 5;
	int minSize = 3;
	bool bDisplay = false;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(blockSize, prhs[2], classMode);
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
	cv::Size destSize(dWidth, dHeight);
	cv::Mat prevFrame = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //previous frame
	int frameCount = 0;
	vector<cv::Mat> detected;
	while (true)
	{
		bool bSuccess = cap.read(frame0); // read a new frame from video
		frameCount++;

		if (!bSuccess) //if not success, break loop
		{
			cout << "Cannot read a frame from video stream" << endl;
			break;
		}
		resize(frame0, frame, destSize);

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
		cv::threshold(cv::abs(diff), fg, (double)thres, 255.0, 0);

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
		cv::Mat skel = skeleton(fg.clone());
		detected.push_back(skel);

		if (bDisplay)
		{
			imshow("input", frame); //show the frame in "MyVideo" window
			imshow("motion", fg);
			imshow("Skeleton", skel);
		}
		if (cv::waitKey(30) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
			cout << "esc key is pressed by user" << endl;
			break;
		}
	}
	cap.release();

	if (nlhs >= 1)
	{
		const int dims[] = { frame.cols, frame.rows, detected.size() };
		vector<unsigned char> F(dims[0] * dims[1] * dims[2], (unsigned char)0);
		for (int i = 0; i < detected.size(); ++i)
		{
			for (int j = 0; j < frame.rows; ++j)
			{
				unsigned char* ptr = detected[i].ptr(j);
				for (int k = 0; k < frame.cols; ++k)
				{
					SetData3(F, k, j, i, dims[0], dims[1], dims[2], ptr[k]);
				}
			}
		}
		plhs[0] = StoreData(F, mxUINT8_CLASS, 3, dims);
	}

	cv::destroyAllWindows();
	mexUnlock();
}
