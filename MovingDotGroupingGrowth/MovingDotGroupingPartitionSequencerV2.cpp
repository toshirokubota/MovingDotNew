// IPCameraCaptureOpenCV.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <mex.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/video/background_segm.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <iterator>
#include <set>
#include <map>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <CoreParticle.h>
#include <CoreParticleMakeGraph.h>
#include <szConnectedComponent.h>
#include <CoreParticle.h>
#include <CoreParticleUtility.h>
#include <NeighborhoodFactory.h>
#include <CoreParticleMakeGraph.h>
#include <szConvexHull2D.h>
//using namespace cv;
using namespace std;

struct RegionInfo
{
	RegionInfo(vector<CoreParticle*>& seq, set<CoreParticle*>& pnts0, int g)
	{
		contour = seq;
		pnts.insert(pnts0.begin(), pnts0.end());
		if (contour.size() < 8) // need at least 8 to surround a region.
		{
			area = 0.0f;
			perimeter = 0.0f;
			measure = 0.0f;
		}
		else
		{
			vector<CParticleF> vp;
			for (vector<CoreParticle*>::iterator it = contour.begin(); it != contour.end(); ++it)
			{
				CoreParticle* p = *it;
				vp.push_back(CParticleF(p->x, p->y, p->z));
				cores.insert(p->src);
			}
			area = polygonArea(vp);
			perimeter = 0;
			for (int i = 0; i < vp.size(); ++i)
			{
				int j = i == vp.size() - 1 ? 0 : i + 1;
				perimeter += length(vp[i].m_X - vp[j].m_X, vp[i].m_Y - vp[j].m_Y, vp[i].m_Z - vp[j].m_Z, 0.0f);
			}
			measure = area / perimeter;
		}
		gen = g;
		error = 0;
		bKeep = false;
	}

	set<CoreParticle*> pnts;
	vector<CoreParticle*> contour;
	set<CoreParticle*> cores;
	float area;
	float perimeter;
	float measure;
	int gen;
	int error; //0 - no error, 1 - indefinite loop, 2 - reaching null, +10 - missing points
	bool bKeep;
};

bool
region_info_less(RegionInfo* a, RegionInfo* b)
{
	return a->measure < b->measure;
}

float
shapeMeasure2(vector<CoreParticle*>& boundary, set<CoreParticle*>& region)
{
	return sqrt(2.0*PI*(float)region.size()) / (float)boundary.size();
}
bool
less_than(CoreParticle* p, CoreParticle* q)
{
	return p->y < q->y || (p->y == q->y && p->x < q->x);
}

/*
Given a region in S (with the given label as its value), construct particles along the outer shell of the region.
*/
set<CoreParticle*> 
traceBoundary(vector<int>& S, vector<CoreParticle*>& prev, 
	int label, int gen,
	int ndim, const int* dims)
{
	int dims2[] = { dims[0] + 2, dims[1] + 2 };
	vector<CoreParticle*> mp(dims2[0] * dims2[1], NULL);
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();

	set<CoreParticle*> bdry;
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			int val = GetData2(S, j, i, dims[0], dims[1], 0);
			if (val == label)
			{
				for (int i2 = i - 1; i2 <= i + 1; ++i2)
				{
					for (int j2 = j - 1; j2 <= j + 1; ++j2)
					{
						int val2 = GetData2(S, j2, i2, dims[0], dims[1], 0);
						if (val2 == 0)
						{
							CoreParticle* q = GetData2(mp, j2 + 1, i2 + 1, dims2[0], dims2[1], (CoreParticle*)NULL);
							if (q == NULL)
							{
								q = factory.makeParticle(j2, i2);
								q->gen = gen;
								q->label = label;
								SetData2(mp, j2 + 1, i2 + 1, dims2[0], dims2[1], q);
								bdry.insert(q);
							}
						}
					}
				}
			}
		}
	}
	//make 4-neighbors
	int xoff[] = { -1, 0, 1, 0 };
	int yoff[] = { 0, -1, 0, 1 };
	for (set<CoreParticle*>::iterator it = bdry.begin(); it != bdry.end(); ++it)
	{
		CoreParticle* p = *it;
		for (int n = 0; n < 4; ++n)
		{
			int x = p->x + xoff[n];
			int y = p->y + yoff[n];
			CoreParticle* q = GetData2(mp, x + 1, y + 1, dims2[0], dims2[1], (CoreParticle*)NULL);
			if (q && q->label == p->label)
			{
				p->neighbors.push_back(q);
			}
		}
	}
	//set ancestors
	for (set<CoreParticle*>::iterator it = bdry.begin(); it != bdry.end(); ++it)
	{
		CoreParticle* p = *it;
		for (int y = p->y - 1; y <= p->y + 1; ++y)
		{
			for (int x = p->x - 1; x <= p->x + 1; ++x)
			{
				if (x == p->x && y == p->y) continue;
				int val = GetData2(S, x, y, dims[0], dims[1], 0);
				if (val == label)
				{
					CoreParticle* q = GetData2(prev, x, y, dims[0], dims[1], (CoreParticle*)NULL);
					if (q)
					{
						p->ascendents.insert(q);
						q->descendents.insert(p);
					}
				}
			}
		}
	}
	return bdry;
}

/*
Cluster particles into spatially disjoint sets.
*/
vector<set<CoreParticle*>>
clusterFront(set<CoreParticle*>& front, bool bEightNeighbor, int ndim, const int* dims)
{
	vector<TK::Node<CoreParticle*>*> nodes;
	for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
	{
		TK::Node<CoreParticle*>* n = TK::makeset(*it);
		nodes.push_back(n);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			CoreParticle* q = nodes[j]->key;
			if (bEightNeighbor)
			{
				if (Abs(p->x - q->x) <= 1 && Abs(p->y - q->y) <= 1 && Abs(p->z - q->z) <= 1) //check for 8-neighbors
				{
					TK::merge(nodes[i], nodes[j]);
				}
			}
			else
			{
				if (Abs(p->x - q->x) + Abs(p->y - q->y) + Abs(p->z - q->z) == 1) //check for 4-neighbors
				{
					TK::merge(nodes[i], nodes[j]);
				}
			}
		}
	}
	vector<vector<TK::Node<CoreParticle*>*>> dsets = TK::grouping(nodes);
	vector<set<CoreParticle*>> res(dsets.size());
	for (int i = 0; i < dsets.size(); ++i)
	{
		set<CoreParticle*> sp;
		for (int j = 0; j < dsets[i].size(); ++j)
		{
			sp.insert(dsets[i][j]->key);
		}
		res[i] = sp;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return res;
}

/*
Among U, find those that are closest to particles in V.
*/
set<CoreParticle*>
closestParticle(set<CoreParticle*> U, set<CoreParticle*>& V)
{
	set<CoreParticle*> choice;
	float mind = std::numeric_limits<float>::infinity();
	int sizeZ = 0;
	for (set<CoreParticle*>::iterator it = U.begin(); it != U.end(); ++it)
	{
		CoreParticle* p = *it;
		set<CoreParticle*> Z;
		float mind0 = std::numeric_limits<float>::infinity();
		for (set<CoreParticle*>::iterator jt = V.begin(); jt != V.end(); ++jt)
		{
			CoreParticle* q = *jt;
			float d = length(p->x - q->x, p->y - q->y, p->z - q->z, 0.0f);
			if (d < mind0)
			{
				Z.clear();
				Z.insert(p);
				mind0 = d;
			}
			else if (d == mind0)
			{
				Z.insert(p);
			}
		}
		if (mind0 < mind)
		{
			mind = mind0;
			choice.clear();
			choice.insert(p);
			sizeZ = Z.size();
		}
		else if (mind0 == mind)
		{
			if (sizeZ < Z.size())
			{
				choice.clear();
				choice.insert(p);
				sizeZ = Z.size();
			}
			else
			{
				choice.insert(p);
			}
		}
	}
	return choice;
}

/*
for a complete set of boundary particles, sequence them.
*/
vector<CoreParticle*>
sequence(set<CoreParticle*>& bdry,
	int& error_code,
	int ndim, const int* dims)
{
	//find the most top-left point and use it as a starting point
	CoreParticle* tl = NULL;
	for (set<CoreParticle*>::iterator it = bdry.begin(); it != bdry.end(); ++it)
	{
		CoreParticle* p = *it;
		p->value = 0; //indidate unvisited site. change it to 1 once visited.
		if (tl == NULL)
		{
			tl = p; //initialize tl
		}
		else if (less_than(p, tl))
		{
			tl = p;
		}
	}
	//search 4-neighbors in this order.
	int xoff[] = { -1, 0, 1, 0 };
	int yoff[] = { 0, -1, 0, 1 };

	vector<CoreParticle*> vp;
	CoreParticle* p = tl;
	CoreParticle* prev = NULL; //previously visited site
	error_code = 0;
	do
	{
		if (p->x + 1 == 118 && p->y + 1 == 73)
		{
			p->x += 0;
		}
		CoreParticle* q = NULL;
		set<CoreParticle*> candidates;
		float minVal = std::numeric_limits<float>::infinity();
		for (vector<CoreParticle*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
		{
			CoreParticle* r = *it;
			if (r)
			{
				if (r == prev) continue;

				if (r->gen == p->gen)
				{
					if (r->value < minVal)
					{
						minVal = r->value;
						candidates.clear();
						candidates.insert(r);
					}
					else if (r->value == minVal)
					{
						candidates.insert(r);
					}
				}
			}
		}
		if (candidates.empty())
		{
			q = prev;
		}
		else if (candidates.size() == 1)
		{
			q = *(candidates.begin());
		}
		else
		{
			CoreParticle* x = p;
			CoreParticle* q0 = NULL;  //backup
			int vpidx = vp.size() - 1;
			do
			{
				set<CoreParticle*> filtered = closestParticle(candidates, x->ascendents);
				if (filtered.size() == 1)
				{
					q = *(filtered.begin());
					break;
				}
				else if (filtered.size() > 1)
				{
					if (vpidx >= 0)
					{
						x = vp[vpidx--];
						q0 = *(filtered.begin());
					}
					else
					{
						if (!vp.empty())
						{
							printf("Cannot resolve tie-breaker at (%d,%d) gen=%d with boundary size=%d.\n",
								p->x, p->y, p->gen, bdry.size());
						}
						q = *(filtered.begin());
						break;
					}
				}
				else 
				{
					if (q0) q = q0;
					else q = *(candidates.begin());
					break;
				}
			} while (true);
		}
		vp.push_back(p);
		p->value += 1;
		prev = p;
		p = q;
		if (p == NULL)
		{
			printf("Sequence: premature exit. NULL resulted from (%d, %d) gen=%d.\n",
				prev->x, prev->y, prev->gen);
			error_code = 2;
			break; //something is wrong
		}
		if (vp.size() >= bdry.size() * 2)
		{
			printf("Sequence: premature exit. Stuck in an end-less loop. %d original pnts, now %d contour pnts. Generation=%d.\n",
				bdry.size(), vp.size(), p->gen);
			error_code = 1;
			break; //something is wrong
		}
		if (p == tl)
		{
			bool bDone = true;
			for (set<CoreParticle*>::iterator it = bdry.begin(); it != bdry.end(); ++it)
			{
				CoreParticle* a = *it;
				if (a->value == 0)
				{
					bDone = false;
					break;
				}
			}
			if (bDone) break; //truely done!
		}
	} while (true);
	if (vp.size() < bdry.size())
	{
		printf("Warning: a sequence with %d points from %d points. Some points missing...\n",
			vp.size(), bdry.size());
		error_code += 10;
	}
	return vp;
}

/*
Inner regions are a bit tricky.
We only concern those that are touching the image boundary, as the other ones (enclosed by the current region) would be
taken care of by the outer disconnected regions. Thus, we remove them to avoid further computation.
The front needs to move inward to comply with the traceBoundary, which adds a shell layer. Thus, we erode the region, except
the image boundary where we need to keep it intact.

Finally, we consider a region that is significantly larger than the others as the background and removes it.
We consider 'the significance' as 5-times as large as the second largest one.
*/
vector<int>
adjustInnerRegions(vector<int>& S,
	int numC, //number of clusters in S
	int ndim, const int* dims)
{
	vector<bool> bKeep(numC + 1, false);
	vector<int> vcount(numC + 1, 0);
	for (int i = 0; i < dims[1]; ++i)
	{
		int val1 = GetData2(S, 0, i, dims[0], dims[1], 0);
		int val2 = GetData2(S, dims[0] - 1, i, dims[0], dims[1], 0);
		bKeep[val1] = true;
		bKeep[val2] = true;
	}
	for (int i = 0; i < dims[0]; ++i)
	{
		int val1 = GetData2(S, i, 0, dims[0], dims[1], 0);
		int val2 = GetData2(S, i, dims[1] - 1, dims[0], dims[1], 0);
		bKeep[val1] = true;
		bKeep[val2] = true;
	}
	vector<int> S2 = S;
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			int val = GetData2(S, j, i, dims[0], dims[1], 0);
			vcount[val] ++;
			if (val && bKeep[val])
			{
				//erode the region
				bool bOuter = false;
				for (int i2 = i - 1; i2 <= i + 1 && !bOuter; ++i2)
				{
					for (int j2 = j - 1; j2 <= j + 1 && !bOuter; ++j2)
					{
						if (GetData2(S, j2, i2, dims[0], dims[1], 1) == 0)
						{
							bOuter = true;
						}
					}
				}
				if (bOuter)
				{
					SetData2(S2, j, i, dims[0], dims[1], 0);
				}
			}
			else if (val && !bKeep[val])
			{
				SetData2(S2, j, i, dims[0], dims[1], 0); //not touching the image boundary. Remove the region
			}
		}
	}

	//try removing the background
	int first = 0, second = 0;
	int first_idx = -1, second_idx = -1;
	for (int i = 1; i < vcount.size(); ++i)
	{
		if (vcount[i] > first)
		{
			second = first;
			first = vcount[i];
			second_idx = first_idx;
			first_idx = i;
		}
		else if (vcount[i] > second)
		{
			second = vcount[i];
			second_idx = i;
		}
	}
	if (first == 0); //no more inner region. do nothing.
	else if (second == 0 ||
		first > 5 * second)
	{
		//remove the largest region.
		for (int i = 0; i < S.size(); ++i)
		{
			if (S[i] == first_idx)
			{
				S2[i] = 0;
			}
		}
	}
	return S2;
}

/*
This routine grows regions from a set of dots (in DOTS) for a fixed number of iterations.
MP gives the spatial arrangement of particles.

The routine also establish a hierarchical structure as initial dot (connected) components are merged.
*/
vector<RegionInfo*>
grow(vector<int>& S,
	vector<unsigned char>& L, 
	int niter,
	bool bEightNeighorhood,
	int ndim, const int* dims)
{
	NeighborhoodFactory& nfactory = NeighborhoodFactory::getInstance(ndim);
	vector<vector<int>> nbh4 = nfactory.neighbor4;
	vector<vector<int>> nbh8 = nfactory.neighbor8;

	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> mp = generateParticleMap(L, ndim, dims);
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			if (GetData2(L, j, i, dims[0], dims[1], (unsigned char)0))
			{
				CoreParticle* p = factory.makeParticle(j, i);
				p->gen = 1;
				SetData2(S, j, i, dims[0], dims[1], 1);
				SetData2(mp, j, i, dims[0], dims[1], p);
			}
		}
	}

	vector<RegionInfo*> regions;
	int iter = 1;
	int gen = 2; //initial front gets 1.
	while (iter <= niter)
	{
		vector<CoreParticle*> particles;
		vector<int> Outer(S.size(), 0);
		vector<int> Inner(S.size(), 0);
		int nc1 = ConnectedComponentAnalysisBigger(Outer, S, NeighborhoodEight, 0, ndim, dims);
		int nc2 = ConnectedComponentAnalysisEqual(Inner, S, NeighborhoodEight, 0, ndim, dims);
		Inner = adjustInnerRegions(Inner, nc2, ndim, dims);

		for (int m = 1; m <= nc1; ++m)
		{
			set<CoreParticle*> front = traceBoundary(Outer, mp, m, gen, ndim, dims);
			particles.insert(particles.end(), front.begin(), front.end());
			vector<set<CoreParticle*>> boundaries = clusterFront(front, false, ndim, dims);
			for (int n = 0; n < boundaries.size(); ++n)
			{
				if (boundaries[n].size() <= 2) continue; //ignore small ones
				int error_code = 0;
				vector<CoreParticle*> seq = sequence(boundaries[n], error_code, ndim, dims);
				RegionInfo* r = new RegionInfo(seq, boundaries[n], gen);
				r->error = error_code;
				regions.push_back(r);
			}
		}
		for (int m = 1; m <= nc2; ++m)
		{
			set<CoreParticle*> front = traceBoundary(Inner, mp, m, gen, ndim, dims);
			particles.insert(particles.end(), front.begin(), front.end());
			vector<set<CoreParticle*>> boundaries = clusterFront(front, false, ndim, dims);
			for (int n = 0; n < boundaries.size(); ++n)
			{
				if (boundaries[n].size() <= 2) continue; //ignore small ones
				int error_code = 0;
				vector<CoreParticle*> seq = sequence(boundaries[n], error_code, ndim, dims);
				RegionInfo* r = new RegionInfo(seq, boundaries[n], gen);
				r->error = error_code;
				regions.push_back(r);
			}
		}
		for (int i = 0; i < particles.size(); ++i)
		{
			SetData2(mp, particles[i]->x, particles[i]->y, dims[0], dims[1], particles[i]);
			SetData2(S, particles[i]->x, particles[i]->y, dims[0], dims[1], gen);
		}
		//saved.insert(saved.end(), regions.begin(), regions.end());
		printf("iter %d: %d regions.\n", iter, regions.size());
		iter++;
		gen++;
	}
	//sort them based on the fitness measure.
	vector<pair<float, RegionInfo*>> pairs;
	for (int i = 0; i < regions.size(); ++i)
	{
		pairs.push_back(pair<float, RegionInfo*>(regions[i]->measure, regions[i]));
	}
	sort(pairs.begin(), pairs.end());
	for (int i = 0; i < regions.size(); ++i)
	{
		regions[i] = pairs[regions.size() - i - 1].second;
	}
	//final clean up
	printf("There are %d regions.\n", regions.size());

	return regions;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: seg = MovingDotGrouping2D(dot_image)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classV;
	vector<unsigned char> im;
	LoadData(im, prhs[0], classV, ndim, &dims);
	bool bEightNeighborhood = true; //this version has to use 8-neighborhood.
	int numIter = 30;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[1], classMode);
	}
	int numKeep = 5;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(numKeep, prhs[2], classMode);
	}
	int nvoxels = numberOfElements(ndim, dims);

	CoreParticleFactory::getInstance().clean(); //initialize the factory before using it.

	/*vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	set<CoreParticle*> dots;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) dots.insert(mp[i]);
	}
	vector<CoreParticle*> dotsv = setupParticleNeighbors(mp, ndim, dims);*/
	vector<int> S(nvoxels, 0);
	vector<RegionInfo*> regions = grow(S, im, numIter, bEightNeighborhood, ndim, dims);
	//vector<int> S2(nvoxels, 0);
	//partitionParticles(dots, S2, numIter, bEightNeighborhood, ndim, dims);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, ndim, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { regions.size(), 1 };
		vector<vector<int>> F; // (dims[0] * dims[1]);
		int ncol = 4;
		for (int i = 0; i < regions.size(); ++i)
		{
			RegionInfo* r = regions[i];
			const int dims2[] = { r->pnts.size(), ncol };
			vector<int> f(dims2[0] * dims2[1]);
			int j = 0;
			for (set<CoreParticle*>::iterator it = r->pnts.begin(); it != r->pnts.end(); ++it, j++)
			{
				CoreParticle* p = *it;
				SetData2(f, j, 0, dims2[0], dims2[1], p->x);
				SetData2(f, j, 1, dims2[0], dims2[1], p->y);
				SetData2(f, j, 2, dims2[0], dims2[1], regions[i]->gen);
				SetData2(f, j, 3, dims2[0], dims2[1], (int)(regions[i]->measure*1000.0));
			}
			//delete regions[i];
			F.push_back(f);
		}
		plhs[1] = StoreDataCell(F, mxINT32_CLASS, ndim, dims, ncol);
	}
	if (nlhs >= 3)
	{
		//sort regions first
		//sort(regions.begin(), regions.end(), region_info_less);
		const int dims[] = { regions.size(), 1 };
		vector<vector<int>> F; // (dims[0] * dims[1]);
		int ncol = 5;
		for (int i = 0; i < regions.size(); ++i)
		{
			RegionInfo* r = regions[i];
			printf("%d Region: %d %d %d %f %f %f %d\n",
				i, r->gen, r->pnts.size(), r->contour.size(), r->measure, r->area, r->perimeter, r->error);
			const int dims2[] = { r->contour.size(), ncol };
			vector<int> f(dims2[0] * dims2[1]);
			for (int j = 0; j < r->contour.size(); ++j)
			{
				CoreParticle* p = r->contour[j];
				SetData2(f, j, 0, dims2[0], dims2[1], p->x);
				SetData2(f, j, 1, dims2[0], dims2[1], p->y);
				SetData2(f, j, 2, dims2[0], dims2[1], regions[i]->gen);
				SetData2(f, j, 3, dims2[0], dims2[1], (int)(regions[i]->measure*1000.0));
				SetData2(f, j, 4, dims2[0], dims2[1], p->id);
			}
			//delete regions[i];
			F.push_back(f);
		}
		plhs[2] = StoreDataCell(F, mxINT32_CLASS, ndim, dims, ncol);
	}
	for (int i = 0; i < regions.size(); ++i)
	{
		delete regions[i];
	}
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}
