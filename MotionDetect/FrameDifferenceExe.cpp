// IPCameraCaptureOpenCV.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/video/background_segm.hpp>
#include <iostream>
#include <vector>
#include <fstream>

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

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <videofile> <logfile>\n", argv[0]);
		exit(1);
	}
	VideoCapture cap; //

	string filename = argv[1];
	string logfile = argv[2];
	//path = path + filename + ".mod";
	cap.open(filename.c_str());
	ofstream os(logfile.c_str());

	//cap.open(0); //open a local camera (web-cam)
	if (!cap.isOpened())  // if not success, exit program
	{
		cout << "Cannot open the video cam" << endl;
		return -1;
	}

	double dWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH); //get the width of frames of the video
	double dHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT); //get the height of frames of the video
	int decRate = 1;
	cout << "Frame size : " << dWidth << " x " << dHeight << endl;

	double fps = cap.get(CV_CAP_PROP_FPS); //get the frames per seconds of the video

	cout << "Frame per seconds : " << fps << endl;

	namedWindow("MyVideo", CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
	namedWindow("MyMotionVideo", CV_WINDOW_AUTOSIZE);

	Mat frame; //after resizing frame
	Mat frame0; //original frame
	Size destSize(dWidth / decRate, dHeight / decRate);
	Mat prevFrame = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //previous frame
	int N = 1;
	int frameCount = 0;
	int blockSize = 10;
	if (argc >= 4)
	{
		blockSize = atoi(argv[3]); 
	}
	int thres = 30;
	int minCount = blockSize*blockSize / 5 + 1;

	os << filename << endl;
	os << dWidth << " " << dHeight << " " << decRate << " " << blockSize << " " << thres << " " << minCount << endl;
	cout << dWidth << " " << dHeight << " " << decRate << " " << blockSize << " " << thres << " " << minCount << endl;

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
		imshow("MyVideo", frame); //show the frame in "MyVideo" window

		Mat gray;
		cv::cvtColor(frame, gray, CV_BGR2GRAY);

		if (frameCount == 1) continue; //skip the first frame since there is no previous frame, yet

		Mat diff = gray - prevFrame;
		Mat blank = cv::Mat::zeros(destSize, cv::DataType<uchar>::type); //blank frame

		//imshow("MyMotionVideo", diff);
		if (frameCount % N == 0)
		{
			vector<Point2i> motions = collectMotions(diff, blockSize, thres, minCount);
			//cout << "#motions: " << motions.size() << endl;

			for (int k = 0; k < motions.size(); ++k)
			{
				circle(blank, motions[k], 2, Scalar(255, 255, 255));
				os << motions[k].x << " " << motions[k].y << " " << frameCount << endl;
			}
			imshow("MyMotionVideo", blank);
		}
		prevFrame = gray;

		if (waitKey(30) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
		{
			cout << "esc key is pressed by user" << endl;
			break;
		}
	}
	os.close();
	return 0;
}

