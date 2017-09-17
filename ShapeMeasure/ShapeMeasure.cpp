// ShapeMeasure.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <mex.h>
#include <szMexUtility.h>
#include <szmexutilitytemplate.h>
#include <mexFileIO.h>
#include <szConvexHull2D.h>
#include <szParticleF.h>
#include <szMiscOperations.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: O = ShapeMeasure(ponints)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classV;
	vector<float> P;
	LoadData(P, prhs[0], classV, ndim, &dims);
	vector<CParticleF> points;
	for (int i = 0; i < dims[0]; ++i)
	{
		float x = GetData2(P, i, 0, dims[0], dims[1], 0.0f);
		float y = GetData2(P, i, 1, dims[0], dims[1], 0.0f);
		points.push_back(CParticleF(x, y));
		//printf("Points: %d %d\n", (int)x, (int)y);
	}
	vector<CParticleF> hull = ConvexHull2D(points);
	/*for (int i = 0; i < hull.size(); ++i)
	{
		printf("Hull: %d %d\n", (int)hull[i].m_X, (int)hull[i].m_Y);
	}*/
	float area = polygonArea(hull);
	vector<float> measure(3);
	measure[0] = points.size() / area;
	measure[1] = area;
	measure[2] = (float)hull.size();
	//printf("measure=%f, area=%f, hull size=%d, #points=%d\n",
	//	measure[0], area, hull.size(), points.size());


	if (nlhs >= 1)
	{
		const int dims2[] = { measure.size(), 1 };
		int ndim2 = 2;
		plhs[0] = StoreData(measure, mxSINGLE_CLASS, ndim2, dims2);
	}
	if (nlhs >= 2)
	{
		const int dims2[] = { hull.size(), 2 };
		int ndim2 = 2;
		vector<float> F(hull.size() * 2);
		for (int i = 0; i < dims2[0]; ++i)
		{
			SetData2(F, i, 0, dims2[0], dims2[1], hull[i].m_X);
			SetData2(F, i, 1, dims2[0], dims2[1], hull[i].m_Y);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, ndim2, dims2);
	}
	mexUnlock();
}


