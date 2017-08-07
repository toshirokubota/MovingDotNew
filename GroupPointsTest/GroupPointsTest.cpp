/*
This routine implements grouping of moving points.
*/
#include <mex.h>

#include <iostream>
#include <fstream>
#include <set>
using namespace std;
#include <stdlib.h>
#include <MotionPoint.h>
#include <MotionPointProcessing.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>

int MotionPoint::_id = 0;
int Factor = 4;
int Width=256;
int Height=256;
float FrameRate=1.0;
float Rate=0.25;
//int NumIter=20;
//float NoiseVar=1.0;
int NumGroup=3;
int PointsPerGroup=8;
int NoisePoints=20;
int NumPoints=NumGroup*PointsPerGroup+NoisePoints;
float Epsilon=0.00001;
float Threshold1=.01;
float Threshold2=.01;
int Seed=54;
//float Sigma1=sqrt(2.0)*NoiseVar;
//float Sigma2=1.0;
//float Sigma3=1.0;

const int MaxNumGroups=10;
//const int MaxCandidates=10;
float MaxRange=25.0;
const int MaxNumPoints = 256;
const int NumMaxCandidates = 10;

float VelocityY[]={20.0, 20.0, -0.0};
float VelocityX[]={-20.0, 20.0, 20.0};

char infile[256];
char outfile[256];



vector<MotionPoint*>
MovePoints(const vector<MotionPoint*>& points, float rate, float var)
{
	int nump=points.size();
	vector<MotionPoint*> res;

	int i;
	for(i=0; i<nump; ++i) 
	{
		MotionPointSynth* p = (MotionPointSynth*)points[i];
		MotionPointSynth* q = new MotionPointSynth(p->x, p->y, p->ux, p->uy, p->groupId);
		q->Move(rate, var);
		res.push_back(q);
	}
	return res;
}

vector<MotionPoint*>
ShufflePoints(const vector<MotionPoint*>& points)
{
	int nump=points.size();
	vector<MotionPoint*> res(nump);

	int i;
	for(i=0; i<nump; ++i)
	{ 
		res[nump-i-1]=points[i];
	}
	return res;
}

vector<MotionPoint*>
PlacePoints(int ng, int nppg, int nrp, int height, int width) 
{
	int nump=ng * nppg + nrp;
	int offy=30;
	int offx=30;
	vector<MotionPoint*> points(nump);
	for(int i=0,m=0; i<ng; ++i) //ng = number of gruops
	{
		for(int j=0; j<nppg; ++j,m++) //nppg = number of points per group
		{
			float y = (height - 2 * offy)*rndm(0) + offy;
			float x = (width - 2 * offx)*rndm(0) + offx;
			MotionPointSynth* p = new MotionPointSynth(x, y, VelocityX[i], VelocityY[i], i);
			points[m] = p;
		}
	}
	for(int m=ng * nppg; m<nump; ++m)
	{
		float x = rndm(0) * width;
		float y = rndm(0) * height;
		float vx = rndm(0) * 50;
		float vy = rndm(0) * 50;
		MotionPointSynth* p = new MotionPointSynth(x, y, vx, vy, ng);
		points[m] = p;
	}
	return points;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	MotionPoint::_id = 0;
	rndm(Seed);

	if (nrhs < 0 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = GroupPoints(numiter, noise_var)");
		return;
	}
	int numIter = 10;
	if (nrhs >= 1)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[0], classMode);
	}
	int numNoise = 20;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(numNoise, prhs[1], classMode);
	}
	float noiseVar = 1.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(noiseVar, prhs[2], classMode);
	}
	float sigma = Max(1.0f, sqrt(2.0) * noiseVar);

	vector<MotionPoint*> frame1 = PlacePoints(NumGroup, PointsPerGroup, numNoise, Height, Width);
	vector<MotionPoint*> frame2 = MovePoints(frame1, FrameRate, noiseVar);
	//frame2=ShufflePoints(frame2); //to make sure the array index is not helping
	EstimateTransform(frame1, frame2, NumMaxCandidates, MaxRange);
	EstimateTransform(frame2, frame1, NumMaxCandidates, MaxRange);

	for (int i = 0; i < frame1.size(); ++i)
	{
		frame1[i]->initializeLinks(frame1);
	}
	for (int i = 0; i < frame2.size(); ++i)
	{
		frame2[i]->initializeLinks(frame2);
	}

	for (int i = 0; i<numIter; ++i) {
		printf("Iteration %d: \n", i);
		UpdateLinkWeights(frame1,Threshold1,sigma);
		UpdateLinkWeights(frame2, Threshold1, sigma);
		UpdateTransformEstimate(frame1, sigma, Rate);
		UpdateTransformEstimate(frame2, sigma, Rate);
		UpdateProbMeasureTmp(frame1, frame2, Threshold2, sigma);
		UpdateProbMeasureTmp(frame2, frame1, Threshold2, sigma);
		UpdateProbMeasure(frame1);
		UpdateProbMeasure(frame2);
	}
	vector<int> vlabel1 = labelIntraFramePoints(frame1);
	vector<int> vlabel2 = labelIntraFramePoints(frame2);
	if (nlhs >= 1)
	{
		const int dims[] = { frame1.size(), 6 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			MotionPointSynth* p = (MotionPointSynth*) frame1[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], (float)vlabel1[i]);
			SetData2(F, i, 3, dims[0], dims[1], p->ux);
			SetData2(F, i, 4, dims[0], dims[1], p->uy);
			SetData2(F, i, 5, dims[0], dims[1], (float)p->groupId);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		const int dims[] = { frame2.size(), 6 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			MotionPointSynth* p = (MotionPointSynth*)frame2[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], (float)vlabel2[i]);
			SetData2(F, i, 3, dims[0], dims[1], p->ux);
			SetData2(F, i, 4, dims[0], dims[1], p->uy);
			SetData2(F, i, 5, dims[0], dims[1], (float)p->groupId);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 3)
	{
		vector<pair<int, int>> matches = findIntraframeMatches(frame1);
		const int dims[] = { matches.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first + 1);
			SetData2(F, i, 1, dims[0], dims[1], matches[i].second + 1);
			SetData2(F, i, 2, dims[0], dims[1], vlabel1[matches[i].first]);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 4)
	{
		vector<pair<int, int>> matches = findIntraframeMatches(frame2);
		const int dims[] = { matches.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first+1);
			SetData2(F, i, 1, dims[0], dims[1], matches[i].second+1);
			SetData2(F, i, 2, dims[0], dims[1], vlabel2[matches[i].first]);
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 5)
	{
		vector<pair<int, int>> matches = findInterframeMatches(frame1, frame2);
		const int dims[] = { matches.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first+1);
			SetData2(F, i, 1, dims[0], dims[1], matches[i].second + 1);
			SetData2(F, i, 2, dims[0], dims[1], vlabel1[matches[i].first]);
			SetData2(F, i, 3, dims[0], dims[1], vlabel2[matches[i].second]);
		}
		plhs[4] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	for (int i = 0; i < frame1.size(); ++i)
	{
		delete frame1[i];
	}
	for (int i = 0; i < frame2.size(); ++i)
	{
		delete frame2[i];
	}
	mexUnlock();
}

