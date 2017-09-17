#pragma once
#include <vector>
using namespace std;
#include <szParticle.h>

namespace TK
{
	struct Canny
	{
		static const int SearchDistance = 1;
		Canny(float pc=0.7f, float r=0.4f, int ml=5, int mg=2)
		{
			percent = pc;
			ratio = r;
			minLength = ml;
			margin = mg;
			ndim = 0;
			dims[0] = dims[1] = 0;
			high = low = 0;
		}
		~Canny()
		{
			for (int i = 0; i < mp.size(); ++i)
			{
				if(mp[i] != NULL) delete mp[i];
			}
		}
		void run(vector<float>& im, int ndim, const int* dims);
		vector<float> retrieve();
		void _EdgeTrace(const vector<float>& mag,
			const vector<float>& gy,
			const vector<float>& gx,
			float high, float low);

		float ratio;
		float percent;
		int minLength;
		int margin; //ignore edges that are constrained within this margin from the image boundary
		vector<vector<CParticle*>> edges;
		vector<CParticle*> mp;
		int ndim;
		int dims[2];
		float low, high; //selected threshold values.
	};
};

vector<float>
CannyNonMaximumSuppression(const vector<float>& dmag,
	const vector<float>& dy,
	const vector<float>& dx,
	int ndim, const int* dims);

vector<float>
EdgeTrace(const vector<float>& mag,
	const vector<float>& gy,
	const vector<float>& gx,
	float high, float low,
	int ndim, const int* dims);

float
ThresholdSelection(const vector<float> grd,
	float percent,
	int ndim, const int* dims);

void
CannyEdge(vector<unsigned char>& edge,
	vector<float>& im,
	float lowThres, float highThres, float sigma,
	int ndim, const int* dims);

void
CannyEdge(vector<unsigned char>& edge,
	vector<float>& im,
	int ndim, const int* dims);
