#ifndef ___DISTANCE_TRANSFORM_NON_ISOTROPIC_H___
#define ___DISTANCE_TRANSFORM_NON_ISOTROPIC_H___

#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

#include "szDistanceTransform.h"

void
DistanceTransformEuclid(vector<double>& D, 
                        const vector<unsigned char>& L, 
						const vector<double>& v,
                        int ndim, 
                        const int* dims);

void
DistanceTransformEuclidF(vector<float>& D, 
                        const vector<unsigned char>& L, 
						const vector<float>& v,
                        int ndim, 
                        const int* dims);

#endif /* ___DISTANCE_TRANSFORM_NON_ISOTROPIC_H___ */