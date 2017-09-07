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
#include <CannyStruct.h>
#include <szConvolutionTemplate.h>

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: P = MotionDetect(video_file)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classV;
	vector<float> im;
	LoadData(im, prhs[0], classV, ndim, &dims);

	float percent = 0.7f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(percent, prhs[1], classMode);
	}
	float ratio = 0.4;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(ratio, prhs[2], classMode);
	}
	int minLength = 5;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minLength, prhs[3], classMode);
	}

	//vector<float> filter = makeGaussianFilter(0.0f, sgm, (int)(2 * sgm) * 2 + 1, true);
	vector<float> filter = makeGaussianFilter(0.0f, 1.0f, 5, true);
	vector<float> sm(numberOfElements(ndim, dims), 0.0f);
	//separableConvolution(sm, im, filter, CBE_MIRROR, ndim, dims);
	sm = im;
	/*vector<float> gr(numberOfElements(ndim, dims), 0.0f);
	vector<float> dx(numberOfElements(ndim, dims), 0.0f);
	vector<float> dy(numberOfElements(ndim, dims), 0.0f);
	for (int i = 1; i < dims[1] - 1; ++i)
	{
		for (int j = 1; j < dims[0] - 1; ++j)
		{
			float nval = GetData2(sm, j, i - 1, dims[0], dims[1], 0.0f);
			float sval = GetData2(sm, j, i + 1, dims[0], dims[1], 0.0f);
			float wval = GetData2(sm, j - 1, i, dims[0], dims[1], 0.0f);
			float eval = GetData2(sm, j + 1, i, dims[0], dims[1], 0.0f);
			float dfx = (eval - wval) / 2;
			float dfy = (sval - nval) / 2;
			SetData2(dx, j, i, dims[0], dims[1], (eval - wval) / 2);
			SetData2(dy, j, i, dims[0], dims[1], (sval - nval) / 2);
			SetData2(gr, j, i, dims[0], dims[1], sqrt(dfx*dfx + dfy*dfy));
		}
	}
	float percent = 0.7f;
	float ratio = 0.4f;
	vector<float> emag = CannyNonMaximumSuppression(gr, dy, dx, ndim, dims);
	high = ThresholdSelection(gr, percent, ndim, dims);
	low = high*ratio;
	emag = EdgeTrace(emag, dy, dx, high, low, ndim, dims);

	//CannyEdge(edge, im, low, high, 0.5f, ndim, dims);
	CannyEdge(edge, im, ndim, dims);*/

	TK::Canny canny(percent, ratio, minLength);
	canny.run(im, ndim, dims);
	vector<float> emag = canny.retrieve();
	vector<unsigned char> edge(im.size(), 0);
	for (int i = 0; i < emag.size(); ++i)
	{
		if (emag[i] > 0) edge[i] = 255;
	}
	printf("Threshold values: %f, %f\n", canny.low, canny.high);
	printf("There are %d edge groups\n", canny.edges.size());

	cv::Size imsize(dims[0], dims[1]);
	vector<unsigned char> cim(im.size());
	for (int i = 0; i < im.size(); ++i)
	{
		cim[i] = (unsigned char)im[i];
	}
	cv::Mat im2 = cv::Mat::zeros(imsize, cv::DataType<unsigned char>::type); //blank frame
	StdVector2Mat(cim, im2, ndim, dims);
	cv::Mat edge2;
	cv::Canny(im2, edge2, canny.low, canny.high); // , 7, true);
	vector<unsigned char> edge3;
	int dims2[2];
	MAT2StdVector(edge3, edge2, ndim, dims2);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(edge, mxUINT8_CLASS, ndim, dims);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreData(edge3, mxUINT8_CLASS, ndim, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(emag, mxSINGLE_CLASS, ndim, dims);
	}
	mexUnlock();
}
