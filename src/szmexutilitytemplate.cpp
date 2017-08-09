//#include <mex.h>
//#include "MexUtilityTemplate.h"
#include <szIndexedData.h>
#include <algorithm>
using namespace std;

//Generic data access
template<class Item>
Item
GetData(const vector<Item>& A, int i, bool& success) {
	if (i>=0 && i<A.size()) {
		success=true;
		return A[i];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData(const vector<Item>& A, int i, Item defval) {
	if (i>=0 && i<A.size()) {
		return A[i];
	}
	else {
		return (Item)defval;
	}
}

template<class Item>
bool
SetData(vector<Item>& A, int i, const Item val) {
	if (i>=0 && i<A.size()) {
		A[i]=val;
		return true;
	}
	else {
		return false;
	}
}


//4D data access
template<class Item>
Item
GetData4(const vector<Item>& A, int x, int y, int z, int t, int xD, int yD, int zD, int tD, bool& success) {
	if (x >= 0 && x<xD && y >= 0 && y<yD && z >= 0 && z<zD && t >= 0 && t<tD) {
		success = true;
		return A[t*xD*yD*zD + z*xD*yD + y*xD + x];
	}
	else {
		success = false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData4(const vector<Item>& A, int x, int y, int z, int t, int xD, int yD, int zD, int tD, Item defval) {
	if (x >= 0 && x<xD && y >= 0 && y<yD && z >= 0 && z<zD && t >= 0 && t<tD) {
		return A[t*xD*yD*zD + z*xD*yD + y*xD + x];
	}
	else {
		return defval;
	}
}

//3D data access
template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, bool& success) {
	if (x>=0 && x<xD && y>=0 && y<yD && z>=0 && z<zD) {
		success=true;
		return A[z*xD*yD+y*xD+x];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, Item defval) {
	if (x>=0 && x<xD && y>=0 && y<yD && z>=0 && z<zD) {
		return A[z*xD*yD+y*xD+x];
	}
	else {
		return defval;
	}
}


template<class Item>
bool
SetData4(vector<Item>& A, int x, int y, int z, int t, int xD, int yD, int zD, int tD, Item val) {
	if (x >= 0 && x<xD && y >= 0 && y<yD && z >= 0 && z<zD && t>=0 && t<tD) {
		A[t*xD*yD*zD + z*xD*yD + y*xD + x] = val;
		return true;
	}
	else {
		return false;
	}
}

template<class Item>
bool
SetData3(vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, Item val) {
	if (x>=0 && x<xD && y>=0 && y<yD && z>=0 && z<zD) {
		A[z*xD*yD+y*xD+x]=val;
		return true;
	}
	else {
		return false;
	}
}

//2D data access
template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, bool& success) {
	if (x>=0 && x<xD && y>=0 && y<yD) {
		success=true;
		return A[y*xD+x];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, Item defval) {
	if (x>=0 && x<xD && y>=0 && y<yD) {
		return A[y*xD+x];
	}
	else {
		return defval;
	}
}

template<class Item>
bool
SetData2(vector<Item>& A, int x, int y, int xD, int yD, Item val) {
	if (x>=0 && x<xD && y>=0 && y<yD) {
		A[y*xD+x]=val;
		return true;
	}
	else {
		return false;
	}
}

//N-D data access
template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, Item defval) {
	int id=0;
	int skip=1;
	bool within=true;
	for(int i=0; i<ndim; ++i) {
		id+=vsub[i]*skip;
		skip*=dims[i];
		if(vsub[i]<0 || vsub[i]>=dims[i]) {
			within=false;
			break;
		}
	}

	if (within && id>=0 && id<skip) {
		return A[id];
	}
	else {
		return defval;
	}
}

template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, bool& success) {
	int id=0;
	int skip=1;
	bool within=true;
	for(int i=0; i<ndim; ++i) {
		id+=vsub[i]*skip;
		skip*=dims[i];
		if(vsub[i]<0 || vsub[i]>=dims[i]) {
			within=false;
			break;
		}
	}

	if (within && id>=0 && id<skip) {
		success=true;
		return A[id];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
bool
SetDataN(vector<Item>& A, vector<int> vsub, const int* dims, int ndim, Item val) 
{
	int id=0;
	int skip=1;
	bool within=true;
	for(int i=0; i<ndim; ++i) 
	{
		id+=vsub[i]*skip;
		skip*=dims[i];
		if(vsub[i]<0 || vsub[i]>=dims[i]) 
		{
			within=false;
			break;
		}
	}

	if (within && id>=0 && id<skip) 
	{
		A[id]=val;
		return true;
	}
	else 
	{
		return false;
	}
}

template<class Item>
bool
getMaximum(const vector<Item>& A, Item& val)
{
	if(A.empty())
		return false;

	val = A[0];
	for(int i=1; i<A.size(); ++i)
	{
		val = Max(A[i], val);
	}
	return true;
}

template<class Item>
bool
getMinimum(const vector<Item>& A, Item& val)
{
	if(A.empty())
		return false;

	val = A[0];
	for(int i=1; i<A.size(); ++i)
	{
		val = Min(A[i], val);
	}
	return true;
}

template<class Item>
bool
getRange(const vector<Item>& A, Item& minVal, Item& maxVal)
{
	if(A.empty())
		return false;

	minVal = A[0];
	maxVal = A[0];
	for(int i=1; i<A.size(); ++i)
	{
		minVal = Min(A[i], minVal);
		maxVal = Max(A[i], maxVal);
	}
	return true;
}

template<class Item>
int
numMoreThan(const vector<Item>& A, const Item& val)
{
	int count = 0;
	for(int i=1; i<A.size(); ++i)
	{
		if(A[i] > val)
		{
			count++;
		}
	}
	return count;
}

template<class Item>
int
numLessThan(const vector<Item>& A, const Item& val)
{
	int count = 0;
	for(int i=1; i<A.size(); ++i)
	{
		if(A[i] < val)
		{
			count++;
		}
	}
	return count;
}

/*
adjust coordinates to comply with an interpolation operation.
*/
template<class Item>
vector<Item>
adjustCoordinates(const vector<Item>& vcoord, 
				  const vector<float>& vnew, 
				  const vector<float>& vold, 
				  int ndim)
{
	int nelm = vcoord.size() / ndim;
	vector<Item> vres(vcoord.size());
	int i, j;
	for(i=0; i<nelm; ++i)
	{
		for(j=0; j<ndim; ++j)
		{
			int k = ndim*i + j;
			Item x = vcoord[k];
			vres[k] = (Item)(x * vold[j] / vnew[j]);
		}
	}
	return vres;
}

/*
reinterpolate a sub-volume so that Z-space is as small as XY spacing.
*/
template<class Item>
vector<Item>
isotropicResamplingZ(const vector<Item>& A, const vector<float>& v, const int* dims) {
	vector<Item> B(A.size());

	float space=Min(v[0],v[1]);
	if(space>=v[2]) {
		B=A;
		return B;
	}
	else {
		int i,j,k;
		int midIdx=(dims[2]-1)/2;
		float inc=space/v[2];
		float z = midIdx - inc*midIdx;
		//printf("isotropicResamplingZ: first Z = %f, space = %f\n", z,space);
		for(i=0; i<dims[2]; ++i, z+=inc) {
			int slice1=(int)z;
			int slice2=slice1+1;
			float t=z-(float)slice1;
			//printf("isotropicResamplingZ: %d: s1=%d, s2=%d, t=%f\n", i, slice1, slice2, t);
			for(j=0; j<dims[1]; ++j) {
				for(k=0; k<dims[0]; ++k) {
					bool bRet;
					Item v1=GetData3(A,k,j,slice1,dims[0],dims[1],dims[2],bRet);
					if(bRet==false) {
						//this cannot happen!!!
						continue;
					}
					Item v2=GetData3(A,k,j,slice2,dims[0],dims[1],dims[2],v1);
					Item val = (Item)((float)v1*(1.0-t)+(float)v2*t);
					//if(j==10 && k==10)
					//  printf("t=%f, v1=%d, v2=%d, val=%d\n", t, (int)v1, (int)v2, (int)val);
					SetData3(B,k,j,i,dims[0],dims[1],dims[2],val);
				}
			}
		}
	}

	return B;
}

template<class Item>
Item
linearInterpolate(const vector<Item>& A,
				  const vector<float>& c,
				  int offset,
				  int ndim,
				  const int* dims)
{
	vector<Item>::const_iterator p = A.begin() + offset;
	if(ndim <= 0)
	{
		return Item(0);
	}
	else if(ndim == 1)
	{
		if(c[0] < 0)
		{
			return *p;
		}
		else if(c[0] >= dims[0]-1)
		{
			return *(p + dims[0] - 1);
		}
		else
		{
			int x1 = (int)c[0];
			int x2 = x1 + 1;
			Item v1 = *(p + x1);
			Item v2 = *(p + x2);
			double t = c[0] - x1;
			Item val = (Item)((double)v1 * (1.0-t) + (double)v2 * t);
			return val;
		}
	}
	else
	{
		int n = numberOfElements(ndim-1, dims);
		int m = ndim - 1;
		int z1, z2;
		if(c[m] <= 0)
		{
			z1 = z2 = 0;
		}
		else if(c[m] >= dims[m]-1)
		{
			z1 = z2 = dims[m]-1;
		}
		else
		{
			z1 = (int)c[m];
			z2 = z1 + 1;
		}

		Item v1 = linearInterpolate(A, c, offset + n * z1, m, dims);
		Item v2 = linearInterpolate(A, c, offset + n * z2, m, dims);
		double t = c[m] - (int)c[m];
		Item val = (Item)((double)v1 * (1.0-t) + (double)v2 * t);
		return val;
	}
}


/*
reinterpolate a sub-volume so that Z-space is as small as XY spacing.
*/
/*template<class Item>
vector<Item>
doResampling(const vector<Item>& A, 
			 const vector<float>& vnew, 
			 const vector<float>& vold,
			 int ndim,
			 const int* dims,
			 int* newdims) 
{
	int* dims2 = new int[ndim];
	int i, j;
	for(i=0; i<ndim; ++i)
	{
		float r = dims[i] * vold[i] / vnew[i];
		dims2[i] = (int)(r + 0.5);
	}
	//printf("doResampling: (%d %d %d) -> (%d %d %d)\n", 
	//	dims[0], dims[1], dims[2], dims2[0], dims2[1], dims2[2]);
	int nvoxels = numberOfElements(ndim, dims2);

	vector<Item> B(nvoxels);
	for(i = 0; i<nvoxels; ++i)
	{
		vector<int> vsub = Ind2Sub(i, ndim, dims2);
		vector<float> vc(ndim);
		for(j=0; j<ndim; ++j)
		{
			vc[j] = vsub[j] * vnew[j] / vold[j];
		}

		Item val = linearInterpolate(A, vc, 0, ndim, dims);
		B[i] = val;
	}

	if(newdims)
	{
		memcpy(newdims, dims2, ndim * sizeof(int));
	}

	delete [] dims2;

	return B;
}*/

/*
linearly interpolate an N-dim array around the center of the array
*/
template<class Item>
bool
doResampling(vector<Item>& B, //interpolated data
			 const vector<Item>& A, //original data
			 const vector<float>& vspace, //spacing of each pixel dimension
			 const vector<float>& vcenterOld, //center of the interpolation
			 const vector<float>& vcenterNew, //center of the interpolation
			 const int* dims_new,	//the dimension of the interpolated data
			 int ndim,	//the number of dimensions
			 const int* dims) //the dimension of the original data.
{
	int i, j;
	int nvoxels = numberOfElements(ndim, dims_new);
	if(B.size() < nvoxels)
	{
		return false;
	}

	vector<float> vc(ndim);
	for(i = 0; i<nvoxels; ++i)
	{
		vector<int> vsub = Ind2Sub(i, ndim, dims_new);
		//bool bOutside = false;
		for(j=0; j<ndim; ++j)
		{
			vc[j] = (vsub[j] - vcenterNew[j]) * vspace[j] + vcenterOld[j];
			//if(vc[j] < 0 || vc[j] > dims[j]-1) {
				//bOutside = true;
				//break;
			//}
		}
		/*if(i < 3 || i >= nvoxels - 3)
		{
			for(j=0; j<ndim; ++j)
			{
				printf("%d, %d: vsub = %d, vc = %f, vsp = %f, vcN = %d, vcO = %d, b=%d\n",
					i, j, vsub[j], vc[j], vspace[j], vcenterNew[j], vcenterOld[j], bOutside);
			}
		}*/
		//if(!bOutside)
		{
			Item val = linearInterpolate(A, vc, 0, ndim, dims);
			B[i] = val;
		}
	}

	return true;
}

template<class Item>
vector<Item>
doThreshold(const vector<Item>& A, 
			const Item& thres, 
			const Item& valueOne,
			const Item& valueZero,
			int ndim, 
			const int* dims)
{
	vector<Item> B(A.size());
	int i;
	for(i=0; i<A.size(); ++i)
	{
		if(A[i] > thres)
		{
			B[i] = valueOne;
		}
		else
		{
			B[i] = valueZero;
		}
	}
	return B;
}

template<class Item>
bool
getMinMaxInteisity(const vector<Item>& A,
				   Item& minV,
				   Item& maxV)
{
	if(A.empty())
	{
		return false;
	}
	minV=A[0];
	maxV=A[0];
	for(int i=1; i<A.size(); ++i)
	{
		if(A[i]<minV)
		{
			minV=A[i];
		}
		if(A[i]>maxV)
		{
			maxV=A[i];
		}
	}
	return true;
}

template<class Item>
vector<Item>
adjustIntensity(const vector<Item>& A,
				double offset,
				double scale)
{
	vector<Item> B(A.size());
	for(int i=1; i<A.size(); ++i)
	{
		B[i] = (Item)(((double)A[i]-offset)*scale);
	}
	return B;
}

template<class Item>
vector<double>
computeMoments(const vector<Item>& values,
			   int maxOrder,
			   bool bCentral,
			   int begin,
			   int end)
{
	vector<double> moments(maxOrder, 0);
	if(values.empty())
		return moments;

	if(end<0)
		end = values.size() - 1;

	int i, j;
	double mean =0;
	for(i=begin; i<=end; ++i)
	{
		mean += values[i];
	}
	mean /= (end-begin + 1);
	moments[0] = mean;

	if(maxOrder > 1)
	{
		if(!bCentral)
			mean = 0; 

		for(i=begin; i<=end; ++i)
		{
			double df = values[i] - mean;
			double prod = df;
			for(j=1; j<maxOrder; ++j)
			{
				prod *= df;
				moments[j] += prod;
			}
		}
		for(j=1; j<maxOrder; ++j)
		{
			moments[j] /= (end-begin + 1);
		}
	}
	return moments;
}

/*
Find closest an index of a sorted array A that is closest to 
the given value (val).  The array is sorted in the ascending order.
bind and eind are optional indices that can restrict the range
of the array.
*/
template<class Item>
int
binarySearch(const vector<Item>& A,
			 const Item& val,
			 int bind,
			 int eind)
{
	bind = Max(0, bind);
	eind = Min(eind, A.size() - 1);
	if(val <= A[bind])
	{
		return bind;
	}
	else if(val>=A[eind])
	{
		return eind;
	}

	bool bFound = false;
	int ind;
	while(!bFound)
	{
		if(bind >= eind)
		{
			bFound = true;
			ind = bind;
		}
		else if(bind == eind - 1)
		{
			bFound = true;
			if(Abs(A[bind]-val) <= Abs(A[eind]-val))
				ind = bind;
			else
				ind = eind;
		}
		else
		{
			int mind = (eind + bind)/2;
			if(val > A[mind])
			{
				bind = mind;
			}
			else if(val < A[mind])
			{
				eind = mind;
			}
			else
			{
				bFound = true;
				ind = mind;
			}
		}
	}

	return ind;
}

template<class Item>
bool
writeVolumeToFile(const char* filename, 
				  const vector<Item>& A, 
				  int ndim, 
				  const int* dims)
{
	ofstream out(filename);
	if(out.fail())
	{
		return false;
	}
	int nvoxels = numberOfElements(ndim, dims);
	Item* buffer = new Item[nvoxels];
	if(!buffer)
	{
		out.close();
		return false;
	}

	for(int i=0; i<nvoxels; ++i)
	{
		buffer[i] = A[i];
	}
	for(int i=0; i<nvoxels; ++i)
	{
		out << buffer[i] << endl;
	}
	
	bool bSuccess = true; //out.write((char*) buffer, nvoxels * sizeof(Item));
	out.close();
	delete [] buffer;

	return bSuccess;
}
		
template<class Item>
vector<int>
IndexedSort(vector<Item>& A)
{
	int i;
	vector<indexedData> vid(A.size());
	for(i=0; i<A.size(); ++i)
	{
		vid[i].data = (double)A[i];
		vid[i].index = i;
	}
	sort(vid.begin(), vid.end());
	
	vector<int> vindex(A.size());
	for(i=0; i<A.size(); ++i)
	{
		A[i] = (Item)vid[i].data;
		vindex[i] = vid[i].index;
	}

	return vindex;
}

/*
for each voxel in B, find the bin number according to vedges and 
store the number in S.
*/
template<class Item>
bool
IndexVector(vector<int>& S, 
			const vector<Item>& B, 
			const vector<Item>& vedges)
{
	int i;
	for(i=0; B.size(); ++i)
	{
		int id = binarySearch(vedges, B[i]);
		S[i] = id;
	}

	return true;
}


template<class Item>
void
gradientMagnitude(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims)
{
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					Item c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					Item w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					Item e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					Item n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					Item s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					Item f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					Item b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float gr = sqrt((float)(w-e)*(w-e) + (n-s)*(n-s) + (f-b)*(f-b));
					SetData3(G, k, j, i, dims[0], dims[1], dims[2], gr);
				}
			}
		}
	}
	return;
}

template<class Item>
void
laplacianResponse(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims)
{
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					Item c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					Item w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					Item e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					Item n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					Item s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					Item f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					Item b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float lp = (float) c - ((float)w + e + s + n + f + b) / 6.0;
					SetData3(G, k, j, i, dims[0], dims[1], dims[2], lp);
				}
			}
		}
	}
	return;
}

template<class Item>
void
gradientAngle(vector<float>& G,
			  const vector<Item>& A,
			  const vector<unsigned char>& S,
			  int x, int y, int z,
			  const int* dims)
{
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		int dz = i-z;
		for(j=0; j<dims[1]; ++j)
		{
			int dy = j-y;
			for(k=0; k<dims[0]; ++k)
			{
				int dx = k-x;
				if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					Item c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					Item w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					Item e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					Item n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					Item s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					Item f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					Item b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float ang = 0;
					if(dz == 0 && dy == 0 && dx==0)
					{
						ang = 0;
					}
					else if(w==e && s==n && f==b)
					{
						ang = 0;
					}
					else
					{
						float l1 = sqrt((float)dx*dx + dy*dy + dz*dz);
						float l2 = sqrt((float)(w-e)*(w-e) + (s-n)*(s-n) + (f-b)*(f-b));
						ang = (dx*(e-w)+dy*(s-n)+dz*(b-f))/(l1*l2);
					}
					SetData3(G, k, j, i, dims[0], dims[1], dims[2], ang);
				}
			}
		}
	}
	return;
}

