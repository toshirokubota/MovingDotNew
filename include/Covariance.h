#pragma once
#include <Feature.h>
#include <cmatrix>
#include <vector>
using namespace std;

class Covariance
{
public:
	static Covariance computeCovariance(const vector<Feature>& x, Feature m);
	static Covariance computeCovariance(const vector<Feature>& x)
	{
		return computeCovariance(x, Feature::computeMean(x));
	}

	Covariance(int n=0)
	{
		ndim = n;
	}

	int ndim;
	techsoft::matrix<float> mat;
	techsoft::matrix<float> imat;
};
