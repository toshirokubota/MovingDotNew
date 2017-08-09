/*
This routine implements grouping of moving points.
*/
#include <mex.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>
using namespace std;
#include <stdlib.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <szParticle.h>
#include <szConnectedComponent.h>
#include <szMyNeighborOp.h>
#include <CoreParticle.h>
#include <CoreParticleUtility.h>

struct MotionRegion {
	vector<CoreParticle*> points;
};

vector<MotionRegion>
initialRegions(vector<unsigned char>& L, int ndim, const int* dims)
{
	vector<int> C(L.size());
	int nc = ConnectedComponentAnalysisBigger(C, L, NeighborhoodEight, (unsigned char)0, ndim, dims);
	vector<MotionRegion> regions(nc);
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			int val = GetData2(C, j, i, dims[0], dims[1], 0);
			if (val > 0)
			{
				int k = val - 1;
				//regions[k].points.push_back(CParticle(j, i));
			}
		}
	}
	return regions;
}

void
paintRegions(vector<int>& C, vector<MotionRegion>& regions, int ndim, const int* dims)
{
}

void
propagate(vector<CoreParticle*>& particles, int ndim, const int* dims)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->src = particles[i];
		particles[i]->gen = 1;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = MotionGrouping2D(P)");
		return;
	}

	vector<unsigned char> L;
	int ndimL;
	const int* dimsL;
	{
		mxClassID classIdP;
		LoadData(L, prhs[0], classIdP, ndimL, &dimsL);
	}
	if (ndimL != 2)
	{
		mexErrMsgTxt("MotionGrouping2D: Input has to be 2D.");
		return;
	}
	int numIter = 1;

	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<unsigned char> ones(L.size(), (unsigned char)1);
	vector<CoreParticle*> mp = generateParticleMap(ones, ndimL, dimsL);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndimL, dimsL);
	vector<CoreParticle*> fg;
	for (int i = 0; i < nvoxels; ++i)
	{
		mp[i]->value = L[i];
		if (L[i] > 0)
		{
			fg.push_back(mp[i]);
		}
	}
	propagate(fg, ndimL, dimsL);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(R, mxINT32_CLASS, ndimL, dimsL);
	}

	mexUnlock();
}

