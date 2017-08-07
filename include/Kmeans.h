#pragma once

#include <vector>
#include <Feature.h>
using namespace std;

int
selectKmeansCluster(Feature& p, vector<Feature>& centers);

vector<Feature>
KmeansClustering(vector<Feature>& points, int K, float tolerance = 1.0e-1);
