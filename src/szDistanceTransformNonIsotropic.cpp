#include <szDistanceTransformNonIsotropic.h>

#define dt_square(a) ((a)*(a))

vector<double>  
dt_Felzenszwalb(const vector<double>& f, double sp, int n) 
{
	vector<double> d(n);
	vector<int> v(n);
	vector<double> z(n+1);
	int k = 0;
	v[0] = 0;
	z[0] = -INFINITY;
	z[1] = +INFINITY;
	for (int i = 1; i <= n-1; i++) 
	{
		double s  = ((f[i]+dt_square(i*sp))-(f[v[k]]+dt_square(sp*v[k])))/(2*i*sp-2*v[k]*sp);
		while (s <= z[k]) 
		{
			k--;
			s  = ((f[i]+dt_square(i*sp))-(f[v[k]]+dt_square(v[k]*sp)))/(2*i*sp-2*v[k]*sp);
		}
		k++;
		v[k] = i;
		z[k] = s;
		z[k+1] = +INFINITY;
	}

	k = 0;
	for (int i = 0; i <= n-1; i++) 
	{
		while (z[k+1] < i*sp)
		{
			k++;
		}
		d[i] = dt_square(i*sp-v[k]*sp) + f[v[k]];
	}
	return d;
}


void
DistanceTransformEuclid(vector<double>& D, 
                        const vector<unsigned char>& L, 
						const vector<double>& v,
                        int ndim, 
                        const int* dims) 
{
	int i,n;

	for(i=0; i<D.size(); ++i) 
	{
		if(L[i])
			D[i]=0;
		else
			D[i]=INFINITY;
	}

	int nvoxels=numberOfElements(ndim,dims);
	int stride=1;
	vector<int> vsub(ndim); //subscript buffer
	for(n=0; n<ndim; ++n) 
	{
		vector<double> buffer(dims[n]);
		int offset=0;
		for(int m=0; m<nvoxels/dims[n]; m++) 
		{
			//compute the offset
			int m2=m;
			int k;
			for(k=0; k<ndim; ++k) 
			{
				if(k!=n) 
				{
					vsub[k]=m2 % dims[k];
					m2/=dims[k];
				}
				else
					vsub[k]=0;
			}
			int offset=0;
			int stride2=1;
			for(k=0; k<ndim; ++k) 
			{
				offset+=vsub[k]*stride2;
				stride2*=dims[k];
			}

			//copy relevant line of data
			for(k=0; k<dims[n]; ++k)
				buffer[k]=GetData(D,offset+k*stride,(double)0);

			//do the computation
			vector<double> dst = dt_Felzenszwalb(buffer, v[n], dims[n]);

			//update the distance
			for(k=0; k<dims[n]; ++k)
				SetData(D,offset+k*stride,(double)(dst[k]));
		}
		stride*=dims[n];
	}

	//take the square root of the distance square
	for(n=0; n<D.size(); ++n)
		D[n]=sqrt(D[n]); //no array limit checking...
}

void
DistanceTransformEuclidF(vector<float>& D, 
                        const vector<unsigned char>& L, 
						const vector<float>& v,
                        int ndim, 
                        const int* dims) 
{
	int i,n;

	for(i=0; i<D.size(); ++i) 
	{
		if(L[i])
			D[i]=0;
		else
			D[i]=INFINITYF;
	}

	int nvoxels=numberOfElements(ndim,dims);
	int stride=1;
	vector<int> vsub(ndim); //subscript buffer
	for(n=0; n<ndim; ++n) 
	{
		vector<double> buffer(dims[n]);
		int offset=0;
		for(int m=0; m<nvoxels/dims[n]; m++) 
		{
			//compute the offset
			int m2=m;
			int k;
			for(k=0; k<ndim; ++k) 
			{
				if(k!=n) 
				{
					vsub[k]=m2 % dims[k];
					m2/=dims[k];
				}
				else
					vsub[k]=0;
			}
			int offset=0;
			int stride2=1;
			for(k=0; k<ndim; ++k) 
			{
				offset+=vsub[k]*stride2;
				stride2*=dims[k];
			}

			//copy relevant line of data
			for(k=0; k<dims[n]; ++k)
				buffer[k]=GetData(D,offset+k*stride,(float)0);

			//do the computation
			vector<double> dst = dt_Felzenszwalb(buffer, (double)v[n], dims[n]);

			//update the distance
			for(k=0; k<dims[n]; ++k)
				SetData(D,offset+k*stride,(float)(dst[k]));
		}
		stride*=dims[n];
	}

	//take the square root of the distance square
	for(n=0; n<D.size(); ++n)
		D[n]=sqrt(D[n]); //no array limit checking...
}
