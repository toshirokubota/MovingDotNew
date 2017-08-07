#pragma once
#include <Feature.h>
#include <Covariance.h>
#include <vector>
#include <szMiscOperations.h>
using namespace std;

class GaussianDistribution
{
public:
	GaussianDistribution()
	{
	}
	GaussianDistribution(vector<Feature>& x)
	{
		mu = Feature::computeMean(x);
		sgm = Covariance::computeCovariance(x, mu);
	}
	float eval(Feature x0)
	{
		int k = x0.length();
		techsoft::matrix<float> m = Feature::toColumnVector(mu);
		techsoft::matrix<float> x = Feature::toColumnVector(x0);
		techsoft::matrix<float> d = (~(m - x))*sgm.imat*(m - x);
		return exp(-d(0, 0) / 2.0) / (pow(2 * PI, k / 2.0) * sqrt(sgm.mat.det()));
	}

	Feature mu;
	Covariance sgm;
};
