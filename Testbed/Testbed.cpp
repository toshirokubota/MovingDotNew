// Testbed.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <mex.h>
#include <vector>
using namespace std;
#include <szMexUtility.h>
#include <szmexutilitytemplate.h>
#include <mexFileIO.h>

/*
separable convolution in semi-optimized way.
*/
void
_separableConvolution2D(float out[], float im[], float g[],
	int n, int ndim, const int* dims)
{
	float* out1 = new float[dims[0] * dims[1]];
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			float s = 0;
			for (int k = 0, k2 = j - n / 2; k < n; ++k, ++k2)
			{
				int j2 = j + k - n / 2;
				j2 = j2 < 0 ? 0 : (j2 >= dims[0] ? dims[0] - 1 : j2);
				int idx = i*dims[0] + j2;
				s += im[idx] * g[k];
			}
			out1[i*dims[0] + j] = s;
		}
	}
	for (int i = 0; i < dims[0]; ++i)
	{
		for (int j = 0; j < dims[1]; ++j)
		{
			float s = 0;
			for (int k = 0; k < n; ++k)
			{
				int j2 = j + k - n / 2;
				j2 = j2 < 0 ? 0 : (j2 >= dims[1] ? dims[1] - 1 : j2);
				int idx = j2*dims[0] + i;
				s += out1[idx] * g[k];
			}
			out[j*dims[0] + i] = s;
		}
	}
	delete out1;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 2 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: O = Testbed(image, filter)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classV;
	vector<float> im;
	LoadData(im, prhs[0], classV, ndim, &dims);

	int ndimF;
	const int* dimsF;
	mxClassID classF;
	vector<float> filter;
	LoadData(filter, prhs[1], classV, ndimF, &dimsF);

	vector<float> out(numberOfElements(ndim, dims));
	
	/*float* _out = new float[out.size()];
	float* _im = new float[im.size()];
	for (int i = 0; i < im.size(); ++i)
	{
		_im[i] = im[i];
	}
	float* _g = new float[filter.size()];
	for (int i = 0; i < filter.size(); ++i)
	{
		_g[i] = filter[i];
	}

	_separableConvolution2D(_out, _im, _g, filter.size(), ndim, dims);

	for (int i = 0; i < out.size(); ++i)
	{
		out[i] = _out[i];
	}
	delete[] _out;
	delete[] _im;
	delete[] _g;*/

	_separableConvolution2D((float*)(&out[0]), (float*)(&im[0]), (float*)(&filter[0]), 
		filter.size(), ndim, dims);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(out, mxSINGLE_CLASS, ndim, dims);
	}
	mexUnlock();
}

