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
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: M = Video2Mat(video_file)");
		return;
	}
	cv::VideoCapture cap;
	char ch[256];
	int len = mxGetString(prhs[0], ch, 256);
	cap.open(ch);
	double decRate = 1.0;
	int sleepTime = 30;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(decRate, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(sleepTime, prhs[2], classMode);
	}

	//cap.open(0); //open a local camera (web-cam)
	if (!cap.isOpened())  // if not success, exit program
	{
		mexErrMsgTxt("Cannot open the video.");
		return;
	}

	double dWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH); //get the width of frames of the video
	double dHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT); //get the height of frames of the video
	double fps = cap.get(CV_CAP_PROP_FPS); //get the frames per seconds of the video
	cv::Size destSize((int)(dWidth / decRate), int(dHeight / decRate));

	cv::namedWindow("input", CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"

	cv::Mat frame; //after resizing frame
	cv::Mat frame0; //original frame
	vector<cv::Mat> frames;

	while (true)
	{
		bool bSuccess = cap.read(frame0); // read a new frame from video

		if (!bSuccess) //if not success, break loop
		{
			cout << "Cannot read a frame from video stream" << endl;
			break;
		}

		cv::resize(frame0, frame, destSize, 0.0, 0.0, cv::INTER_AREA);
		cv::Mat gray;
		cv::cvtColor(frame, gray, CV_BGR2GRAY);

		cv::imshow("input", gray); //show the frame in "MyVideo" window
		frames.push_back(gray);
		if (cv::waitKey(sleepTime) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
			cout << "esc key is pressed by user" << endl;
			break;
		}
	}
	cap.release();

	cv::destroyAllWindows();
	if (nlhs >= 1)
	{
		const int dims[] = { (int)destSize.width, (int)destSize.height, frames.size() };
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
	mexUnlock();
}

