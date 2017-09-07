#pragma once
#include <vector>
using namespace std;

vector<float>
CannyNonMaximumSuppression(const vector<float>& dmag,
	const vector<float>& dy,
	const vector<float>& dx,
	int ndim, const int* dims);

vector<float>
EdgeTrace(const vector<float>& mag,
	const vector<float>& gy,
	const vector<float>& gx,
	float high, float low,
	int ndim, const int* dims);

float
ThresholdSelection(const vector<float> grd,
	float percent,
	int ndim, const int* dims);

void
CannyEdge(vector<unsigned char>& edge,
	vector<float>& im,
	float lowThres, float highThres, float sigma,
	int ndim, const int* dims);

void
CannyEdge(vector<unsigned char>& edge,
	vector<float>& im,
	int ndim, const int* dims);
