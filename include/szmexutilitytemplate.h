#ifndef ___MY_MEX_UTILITY_TEMPLATE_H____
#define ___MY_MEX_UTILITY_TEMPLATE_H____

#include <fstream>
#include <string>
#include <cmath>
#include <vector>
using namespace std;

#include "szMexUtility.h"

//Generic data assess
template<class Item>
Item
GetData(const vector<Item>& A, int i, bool& success);

template<class Item>
Item
GetData(const vector<Item>& A, int i, Item defval);

template<class Item>
bool
SetData(vector<Item>& A, int i, const Item val);

//4D data access
template<class Item>
Item
GetData4(const vector<Item>& A, int x, int y, int z, int t, int xD, int yD, int zD, int tD, bool& success);

template<class Item>
Item
GetData4(const vector<Item>& A, int x, int y, int z, int t, int xD, int yD, int zD, int tD, Item defval);

template<class Item>
bool
SetData4(vector<Item>& A, int x, int y, int z, int t, int xD, int yD, int zD, int tD, const Item val);

//3D data access
template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, bool& success);

template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, Item defval);

template<class Item>
bool
SetData3(vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, const Item val);

//2D data access
template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, bool& success);

template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, Item defval);

template<class Item>
bool
SetData2(vector<Item>& A, int x, int y, int xD, int yD, const Item val);

//N-D data access
template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, Item& defval);

template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, bool& success);

template<class Item>
bool
SetDataN(vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, const Item& val);

template<class Item>
bool
getMaximum(const vector<Item>& A, Item& val);

template<class Item>
bool
getMinimum(const vector<Item>& A, Item& val);

template<class Item>
bool
getRange(const vector<Item>& A, Item& minVal, Item& maxVal);

template<class Item>
int
numMoreThan(const vector<Item>& A, const Item& val);

template<class Item>
int
numLessThan(const vector<Item>& A, const Item& val);

template<class Item>
vector<Item>
isotropicResamplingZ(const vector<Item>& A, const vector<float>& v, const int* dims);

template<class Item>
Item
linearInterpolate(const vector<Item>& A,
                  const vector<float>& c,
                  int offset,
                  int ndim,
                  const int* dims);

template<class Item>
bool
doResampling(vector<Item>& B, //interpolated data
			 const vector<Item>& A, //original data
			 const vector<float>& vspace, //spacing of each pixel dimension
			 const vector<float>& vcenterOld, //center of the interpolation
			 const vector<float>& vcenterNew, //center of the interpolation
			 const int* dims_new,	//the dimension of the interpolated data
			 int ndim,	//the number of dimensions
			 const int* dims //the dimension of the original data.
			 );

template<class Item>
vector<Item>
adjustCoordinates(const vector<Item>& vcoord, 
				  const vector<float>& vnew, 
				  const vector<float>& vold, 
				  int ndim);

template<class Item>
vector<Item>
doThreshold(const vector<Item>& A, 
			const Item& thres, 
			const Item& valueOne,
			const Item& valueZero,
			int ndim, 
			const int* dims);

template<class Item>
bool
getMinMaxInteisity(const vector<Item>& A,
				   Item& minV,
				   Item& maxV);

template<class Item>
vector<Item>
adjustIntensity(const vector<Item>& A,
				double offset,
				double scale);

template<class Item>
vector<double>
computeMoments(const vector<Item>& A,
			   int order, 
			   bool bCentral = true,
			   int begin = 0,
			   int end = -1);

template<class Item>
int
binarySearch(const vector<Item>& A,
			 const Item& val,
			 int bind = 0,
			 int eind = MAX_INTEGER);

template<class Item>
bool
writeVolumeToFile(const char* filename, 
				  const vector<Item>& A, 
				  int ndim, 
				  const int* dims);

template<class Item>
vector<int>
IndexedSort(vector<Item>& A);

template<class Item>
bool
IndexVector(vector<int>& S, 
			const vector<Item>& B, 
			const vector<Item>& vedges);

template<class Item>
void
gradientMagnitude(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims);

template<class Item>
void
laplacianResponse(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims);

template<class Item>
void
gradientAngle(vector<float>& G,
			  const vector<Item>& A,
			  const vector<unsigned char>& S,
			  int x, int y, int z,
			  const int* dims);


#include "../src/szMexUtilityTemplate.cpp"

#endif /* ___MY_MEX_UTILITY_TEMPLATE_H____ */