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
#include <CoreParticle.h>
//using namespace cv;
using namespace std;
int xoff[] = { 0, -1, 1, 0, -1, 1, -1, 1};
int yoff[] = { -1, 0, 0, 1, -1, -1, 1, 1};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: seg = MovingDotGrouping2D(dot_image)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classV;
	vector<unsigned char> im;
	LoadData(im, prhs[0], classV, ndim, &dims);
	int numNeighbors = 4;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		int val = 0;
		ReadScalar(val, prhs[1], classMode);
		if (val) numNeighbors = 8;
	}

	int nvoxels = numberOfElements(ndim, dims);

	vector<int> S(nvoxels, 0);
	vector<CoreParticle*> D(nvoxels, NULL);
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			if (GetData2(im, j, i, dims[0], dims[1], (unsigned char)0) != 0)
			{
				SetData2(D, j, i, dims[0], dims[1], factory.makeParticle(j, i));
			}
		}
	}
	set<CoreParticle*> front;
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			CoreParticle* p = GetData2(D, j, i, dims[0], dims[1], (CoreParticle*)NULL);
			if (p != NULL)
			{
				if (i == 0 || j == 0 || i == dims[1] - 1 || j == dims[0] - 1) continue;

				for (int k = 0; k < numNeighbors; ++k)
				{
					int i2 = p->y + yoff[k];
					int j2 = p->x + xoff[k];
					CoreParticle* q = GetData2(D, j2, i2, dims[0], dims[1], (CoreParticle*)NULL);
					if(q == NULL)
					{
						front.insert(p);
						break;
					}
				}
			}
		}
	}
	int iter = 0;
	while (!front.empty())
	{
		set<CoreParticle*> next;
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			SetData2(S, p->x, p->y, dims[0], dims[1], iter);

			for (int k = 0; k < numNeighbors; ++k)
			{
				int i2 = p->y + yoff[k];
				int j2 = p->x + xoff[k];
				if (i2 < 0 || j2 < 0 || i2 >= dims[1] || j2 >= dims[0]) continue;

				CoreParticle* q = GetData2(D, j2, i2, dims[0], dims[1], (CoreParticle*)NULL);
				if (q == NULL)
				{
					q = factory.makeParticle(j2, i2);
					SetData2(D, j2, i2, dims[0], dims[1], q);
					next.insert(q);
				}
			}
		}
		front = next;
		iter++;
	}
	if (nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, 2, dims);
	}

	factory.clean();
	mexUnlock();
}
