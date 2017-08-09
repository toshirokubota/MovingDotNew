#pragma once
#include <vector>
#include <algorithm>
using namespace std;
#include <szMyNeighborOp.h>
#include <szParticle4D.h>

float
length(float x, float y, float z, float t)
{
	return sqrt(x*x + y*y + z*z + t*t);
}

vector<vector<int>> MakeEightNeighborhoodHere(int n)
{
	if (n == 2)
	{
		vector<vector<int>> idx(8);
		int xoff[] = { 0, -1, 1, 0, -1, 1, -1, 1 };
		int yoff[] = { -1, 0, 0, 1, -1, -1, 1, 1 };
		for (int i = 0; i < idx.size(); ++i)
		{
			vector<int> jdx(2);
			jdx[0] = xoff[i];
			jdx[1] = yoff[i];
			idx[i] = jdx;
		}
		return idx;
	}
	else if (n == 3)
	{
		vector<pair<float, CParticle4D>> vals;
		for (int x = -1; x <= 1; ++x)
		{
			for (int y = -1; y <= 1; ++y)
			{
				for (int z = -1; z <= 1; ++z)
				{
					if (!(x == 0 && y == 0 && z == 0))
					{
						CParticle4D p(x, y, z, 0);
						float len = length(x, y, z, 0);
						vals.push_back(pair<float, CParticle4D>(len, p));
					}
				}
			}
		}
		sort(vals.begin(), vals.end());
		vector<vector<int>> idx;
		for (int i = 0; i < vals.size(); ++i)
		{
			vector<int> jdx(3);
			jdx[0] = vals[i].second.m_X;
			jdx[1] = vals[i].second.m_Y;
			jdx[2] = vals[i].second.m_Z;
			idx.push_back(jdx);
		}
		return idx;
	}
	else
	{
		return MakeEightNeighborhood(n);
	}
}

struct NeighborhoodFactory
{
public:
	static NeighborhoodFactory& getInstance(int n = 0)
	{
		static NeighborhoodFactory instance(n);
		if (n > 0 && n != instance.ndim)
		{
			instance.neighbor4 = MakeFourNeighborhood(n);
			instance.neighbor8 = MakeEightNeighborhoodHere(n);
			instance.ndim = n;
		}
		return instance;
	}
	vector<vector<int>> neighbor4;
	vector<vector<int>> neighbor8;
private:
	int ndim;
	NeighborhoodFactory(int n)
	{
	}
	~NeighborhoodFactory()
	{
	}
	NeighborhoodFactory(NeighborhoodFactory& f) {}
	NeighborhoodFactory operator=(NeighborhoodFactory& f) {}
};
