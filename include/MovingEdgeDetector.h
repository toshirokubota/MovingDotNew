#pragma once
#include <CannyStruct.h>

namespace TK
{
	struct MovingEdgeDetector
	{
		static const int NumDetectors = 2;

		MovingEdgeDetector(float thres, float perc, float r, float ml)
		{
			threshold = thres;
			percent = perc;
			ratio = r;
			minLength = ml;
			for (int i = 0; i < NumDetectors; ++i)
			{
				edgeDetector[i] = Canny(percent, ratio, minLength);
				valid[i] = false;
			}
			next = 0;
		}
		~MovingEdgeDetector()
		{
			for (int i = 0; i < mp.size(); ++i)
			{
				if (mp[i] != NULL) delete mp[i];
			}
		}
		void run(vector<float>& im, int ndim, const int* dims);
		vector<unsigned char> retrieveMovingEdges();
		vector<unsigned char> retrieveEdges();
		vector<float> retrieveMovingEdgeStrength();
		vector<float> retrieveEdgeStrength();

		Canny edgeDetector[NumDetectors];
		bool valid[NumDetectors];
		int next;
		float threshold; //between 0 and 1 - filter out moving or not moving.
		float ratio;
		float percent;
		int minLength;
		vector<vector<CParticle*>> edges;
		vector<CParticle*> mp;
		int ndim;
		int dims[2];
	};
};