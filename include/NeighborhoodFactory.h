#pragma once
#include <vector>
#include <algorithm>
using namespace std;
#include <szMyNeighborOp.h>
#include <szParticle4D.h>


vector<vector<int>> MakeEightNeighborhoodHere(int n);

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
