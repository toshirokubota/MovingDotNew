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

using namespace cv;
using namespace std;

vector<Point2i>
collectMotions(Mat& frame, int blockSize, int thres, int minCount)
{
	vector<Point2i> motions;
	int height = frame.size().height;
	int width = frame.size().width;
	for (int i = 0; i < height; i += blockSize)
	{
		for (int j = 0; j < width; j += blockSize)
		{
			int count = 0;
			for (int y = i; y < min(i + blockSize, height); ++y)
			{
				for (int x = j; x < min(j + blockSize, width); ++x)
				{
					if ((*frame.ptr(y, x)) > thres)
					{
						count++;
					}
				}
			}
			if (count >= minCount)
			{
				motions.push_back(Point2i(j + blockSize / 2, i + blockSize / 2));
			}
		}
	}
	return motions;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: P = MotionDetect(video_file)");
		return;
	}
	VideoCapture cap; 
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
		namedWindow("MyVideo", CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
		namedWindow("MyMotionVideo", CV_WINDOW_AUTOSIZE);
	}

	Mat frame; //after resizing frame
	Mat frame0; //original frame
	Size destSize(dWidth, dHeight);
	Mat prevFrame = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //previous frame
	int frameCount = 0;
	vector<Mat> detected;
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

		Mat gray;
		cv::cvtColor(frame, gray, CV_BGR2GRAY);

		if (frameCount <= 1)
		{
			prevFrame = gray;
			continue; //skip the first frame since there is no previous frame, yet
		}

		Mat diff = gray - prevFrame;
		prevFrame = gray;
		Mat fg = cv::Mat::zeros(destSize, cv::DataType<int>::type); //blank frame
		threshold(abs(diff), fg, (double)thres, 255.0, 0);

		Mat label, stats, centroids;
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
		detected.push_back(fg);

		if (bDisplay)
		{
			imshow("MyVideo", frame); //show the frame in "MyVideo" window
			imshow("MyMotionVideo", fg);
		}
		if (waitKey(30) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
			cout << "esc key is pressed by user" << endl;
			break;
		}
	}
	cap.release();

	if (nlhs >= 1)
	{
		const int dims[] = { frame.cols, frame.rows, detected.size()};
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

	destroyAllWindows();
	mexUnlock();
}

