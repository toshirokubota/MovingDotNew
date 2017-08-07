#pragma once

#include <vector>
#include <Feature.h>
#include <GaussianDistribution.h>
using namespace std;

int
selectGaussianCluster(Feature& p, vector<GaussianDistribution>& gaussians);

vector<GaussianDistribution>
EMGaussianClustering(vector<Feature>& points, int K, int maxiter = 10, float tolerance = 1.0e-1);
