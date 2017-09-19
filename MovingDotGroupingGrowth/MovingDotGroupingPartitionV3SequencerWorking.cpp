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
	RegionInfo(vector<CoreParticle*>& seq, vector<CoreParticle*>& pnts0, int g)
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

/*float
shapeMeasure(set<CoreParticle*>& particles)
{
	vector<CParticleF> pnts;
	for (set<CoreParticle*>::iterator it = particles.begin(); it != particles.end(); ++it)
	{
		CoreParticle* p = *it;
		pnts.push_back(CParticleF(p->x, p->y));
	}
	vector<CParticleF> hull = ConvexHull2D(pnts);
	float area = polygonArea(hull);
	return pnts.size() / area;
}*/

float
shapeMeasure2(vector<CoreParticle*>& boundary, set<CoreParticle*>& region)
{
	return sqrt(2.0*PI*(float)region.size()) / (float)boundary.size();
}

/*
Cluster particles into spatially disjoint sets.
*/
vector<vector<CoreParticle*>> 
clusterFront(set<CoreParticle*>& front, bool bEightNeighbor, int ndim, const int* dims)
{
	vector<TK::Node<CoreParticle*>*> nodes;
	for (set<CoreParticle*>::iterator it=front.begin(); it != front.end(); ++it)
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
				if (Abs(p->x - q->x) <=1 && Abs(p->y - q->y) <= 1 && Abs(p->z - q->z) <= 1) //check for 8-neighbors
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
	vector<vector<CoreParticle*>> res(dsets.size());
	for (int i = 0; i < dsets.size(); ++i)
	{
		vector<CoreParticle*> vp;
		for (int j = 0; j < dsets[i].size(); ++j)
		{
			vp.push_back(dsets[i][j]->key);
		}
		res[i] = vp;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return res;
}

bool
less_than (CoreParticle* p, CoreParticle* q)
{
	return p->y < q->y || (p->y == q->y && p->x < q->x);
}

/*
Among U, find those that are closest to particles in V.
*/
set<CoreParticle*>
closestParticle(set<CoreParticle*> U, set<CoreParticle*>& V)
{
	set<CoreParticle*> Z;
	float mind = std::numeric_limits<float>::infinity();
	for (set<CoreParticle*>::iterator it = U.begin(); it != U.end(); ++it)
	{
		CoreParticle* p = *it;
		for (set<CoreParticle*>::iterator jt = V.begin(); jt != V.end(); ++jt)
		{
			CoreParticle* q = *jt;
			float d = length(p->x - q->x, p->y - q->y, p->z - q->z, 0.0f);
			if (d < mind)
			{
				Z.clear();
				Z.insert(p);
				mind = d;
			}
			else if (d == mind)
			{
				Z.insert(p);
			}
		}
	}
	return Z;
}


/*
for a complete set of boundary particles, sequence them.
*/
vector<CoreParticle*> 
sequence(vector<CoreParticle*>& bdry, 
	int& error_code,
	int ndim, const int* dims)
{
	//find the most top-left point and use it as a starting point
	CoreParticle* tl = NULL;
	int gen;
	for (int i = 0; i < bdry.size(); ++i)
	{
		CoreParticle* p = bdry[i];
		if (i == 0)
		{
			gen = p->gen; //remember the generation number
			tl = p; //initialize tl
		}
		p->gen = -1; //to distinguish quickly from particles with the same gen but not in this cluster
		p->value = 0; //indidate unvisited site. change it to 1 once visited.
		if (less_than(p, tl))
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
		if (p->x == 162 && p->y == -1)
		{
			p->value += 0;
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
								p->x, p->y, gen, bdry.size());
						}
						q = *(filtered.begin());
						break;
					}
				}
				else //no choice...
				{
					q = q0;
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
				prev->x, prev->y, gen);
			error_code = 2;
			break; //something is wrong
		}
		if (vp.size() >= bdry.size() * 2)
		{
			printf("Sequence: premature exit. Stuck in an end-less loop. %d original pnts, now %d contour pnts. Generation=%d.\n",
				bdry.size(), vp.size(), gen);
			error_code = 1;
			break; //something is wrong
		}
		if(p == tl )
		{
			bool bDone = true;
			for (vector<CoreParticle*>::iterator it = bdry.begin(); it != bdry.end(); ++it)
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
	for (int i = 0; i < bdry.size(); ++i)
	{
		bdry[i]->gen = gen; //restore the generation number
	}
	return vp;
}

/*
Given a sequenced path around a region, project the path to initial shape.
*/
vector<CoreParticle*>
projectToSources(vector<CoreParticle*>& vp, int ndim, const int* dims)
{
	vector<CoreParticle*> proj;
	set<CoreParticle*> sp;
	for (int i = 0; i < vp.size(); ++i)
	{
		CoreParticle* p = vp[i];
		for (set<CoreParticle*>::iterator it = p->sources.begin(); it != p->sources.end(); ++it)
		{
			CoreParticle* q = *it;
			//if (sp.find(q) == sp.end())
			{
				proj.push_back(q);
				sp.insert(q);
			}
		}
	}
	return proj;
}

/*
This routine grows regions from a set of dots (in DOTS) for a fixed number of iterations.
MP gives the spatial arrangement of particles.

The routine also establish a hierarchical structure as initial dot (connected) components are merged.
*/
vector<RegionInfo*>
grow(vector<int>& S,
	vector<CoreParticle*>& mp,
	set<CoreParticle*>& particles,
	int numKeep,
	int niter,
	bool bEightNeighorhood,
	int ndim, const int* dims)
{
	NeighborhoodFactory& nfactory = NeighborhoodFactory::getInstance(ndim);
	vector<vector<int>> nbh4 = nfactory.neighbor4;
	vector<vector<int>> nbh8 = nfactory.neighbor8;
	vector<vector<int>> nbh = nbh4;
	if (bEightNeighorhood)
	{
		nbh = nbh8;
	}

	set<CoreParticle*> front;
	front.insert(particles.begin(), particles.end());
	for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
	{
		CoreParticle* p = *it;
		p->gen = 1;
		p->src = p;
		p->sources.insert(p);
		TK::Node<CoreParticle*>* n = TK::makeset(p);
		SetVoxel(S, p, 1, ndim, dims);
	}
	//vector<int> L0;
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<RegionInfo*> allRegions;
	vector<RegionInfo*> regions;
	vector<RegionInfo*> saved;
	vector<RegionInfo*> errornous;
	int nc0 = 0;
	int iter = 1;
	int gen = 2; //initial front gets 1.
	while (!front.empty() && iter <= niter)
	{
		vector<RegionInfo*> regions0 = regions;
		regions.clear();
		vector<int> loc(ndim);
		set<CoreParticle*> next;
		vector<vector<CoreParticle*>> groups = clusterFront(front, true, ndim, dims);
		vector<set<CoreParticle*>> expanded;
		//for each cluster, get the neighbor shell - there may be overap with another cluster.
		for (int m = 0; m < groups.size(); ++m)
		{
			int dims2[] = { dims[0] + 2, dims[1] + 2 };
			vector<CoreParticle*> mp2(dims2[0]*dims2[1], NULL);
			vector<CoreParticle*> vg = groups[m];
			set<CoreParticle*> sp;
			for (vector<CoreParticle*>::iterator it = vg.begin(); it != vg.end(); ++it)
			{
				CoreParticle* p = *it;
				for (int n = 0; n < nbh.size(); ++n)
				{
					int x2 = p->x + nbh[n][0];
					int y2 = p->y + nbh[n][1];
					if (GetData2(S, x2, y2, dims[0], dims[1], 0) == 0)
					{
						CoreParticle* q = GetData2(mp2, x2 + 1, y2 + 1, dims2[0], dims2[1], (CoreParticle*)NULL);
						if (q == NULL)
						{
							q = factory.makeParticle(x2, y2);
							SetData2(mp2, q->x + 1, q->y + 1, dims2[0], dims2[1], q);
							q->src = p->src; ///TK!! we need to consider multipe sources.
						}
						if (Abs(p->x - q->x) + Abs(p->y - q->y) + Abs(p->z - q->z) == 1)
						{
							p->neighbors.push_back(q);
							q->neighbors.push_back(p);
						}
						if (Abs(p->x - q->x) <= 1 && Abs(p->y - q->y) <= 1 && Abs(p->z - q->z) <= 1)
						{
							p->neighbors8.push_back(q);
							q->neighbors8.push_back(p);
						}
						sp.insert(q);
						q->gen = gen;
					}
				}
			}
			for (set<CoreParticle*>::iterator it = sp.begin(); it != sp.end(); ++it)
			{
				CoreParticle* p = *it;
				set<CoreParticle*>::iterator jt = it;
				jt++;
				for (; jt != sp.end(); ++jt)
				{
					CoreParticle* q = *jt;
					if (Abs(p->x - q->x) + Abs(p->y - q->y) + Abs(p->z - q->z) == 1)
					{
						p->neighbors.push_back(q);
						q->neighbors.push_back(p);
					}
				}
			}
			for (vector<CoreParticle*>::iterator it = vg.begin(); it != vg.end(); ++it)
			{
				CoreParticle* p = *it;
				vector<CoreParticle*> vq = closestDescendents(p, bEightNeighorhood);
				for (int j = 0; j < vq.size(); ++j)
				{
					p->descendents.insert(vq[j]);
					vq[j]->ascendents.insert(p);
					vq[j]->sources.insert(p->sources.begin(), p->sources.end());
				}
			}
			for (set<CoreParticle*>::iterator it = sp.begin(); it != sp.end(); ++it)
			{
				CoreParticle* p = *it;
				vector<CoreParticle*> vq = closestAscendents(p, bEightNeighorhood);
				for (int j = 0; j < vq.size(); ++j)
				{
					p->ascendents.insert(vq[j]);
					vq[j]->descendents.insert(p);
					p->sources.insert(vq[j]->sources.begin(), vq[j]->sources.end());
				}
			}
			vector<vector<CoreParticle*>> bdries = clusterFront(sp, false, ndim, dims);
			for (int i = 0; i < bdries.size(); ++i)
			{
				if (bdries[i].size() > 1) //can't trace with a single point.
				{
					int error_code;
					vector<CoreParticle*> seq = sequence(bdries[i], error_code, ndim, dims);
					//RegionInfo* r = new RegionInfo(projectToSources(seq, ndim, dims), bdries[i], iter);
					RegionInfo* r = new RegionInfo(seq, bdries[i], iter);
					r->error = error_code;
					regions.push_back(r);
					allRegions.push_back(r);
					if (r->error)
					{
						errornous.push_back(r);
					}
				}
			}
			expanded.push_back(sp);
		}
		//saved.insert(saved.end(), regions.begin(), regions.end());
		printf("iter %d: %d regions.\n", iter, regions.size());
		//for each region, check if it has been merged.
		for (int i = 0; i < regions0.size(); ++i)
		{
			RegionInfo* r0 = regions0[i];
			if (r0->error) continue;
			for (int j = 0; j < regions.size(); ++j)
			{
				RegionInfo* r = regions[j];
				if (r->error) continue;

				//no error. Check if merged. Then fitness.
				set<CoreParticle*> diff1;
				set_difference(r0->cores.begin(), r0->cores.end(), r->cores.begin(), r->cores.end(), std::inserter(diff1, diff1.begin()));
				if (diff1.empty()) continue; 
				set<CoreParticle*> diff2;
				set_difference(r->cores.begin(), r->cores.end(), r0->cores.begin(), r0->cores.end(), std::inserter(diff2, diff2.begin()));
				if (diff2.empty()) continue;
				//if (r0->measure > r->measure &&
				//	(saved.size() < numKeep || saved[saved.size() - 1]->measure < r->measure))
				if(r0->measure > r->measure)
				{
					int k = saved.size() - 1;
					while (k>=0 && saved[k]->measure < r->measure)
					{
						k--;
					}
					saved.insert(saved.begin() + (k + 1), r0);
					if (saved.size() > numKeep)
					{
						saved.resize(numKeep);
					}
				}
				break;
			}
		}
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			strongMedialParticle(p, ndim, dims);
		}
		for (int i = 0; i < expanded.size(); ++i)
		{
			for (set<CoreParticle*>::iterator it = expanded[i].begin(); it != expanded[i].end(); ++it)
			{
				CoreParticle* p = *it;
				if (p->x < 0 || p->x >= dims[0] || p->y < 0 || p->y >= dims[1]) continue; //ignore outside the image
				CoreParticle* q = GetData2(mp, p->x, p->y, dims[0], dims[1], (CoreParticle*)NULL);
				if (q == NULL)
				{
					SetData2(mp, p->x, p->y, dims[0], dims[1], p);
					next.insert(p);
				}
				else
				{
					q->ascendents.insert(p->ascendents.begin(), p->ascendents.end());
					q->sources.insert(p->sources.begin(), p->sources.end());
					//next.insert(q);
				}
			}
		}
		for (set<CoreParticle*>::iterator it = next.begin(); it != next.end(); ++it)
		{
			CoreParticle* p = *it;
			SetData2(S, p->x, p->y, dims[0], dims[1], gen);
		}
		front = next;
		iter++;
		gen++;
	}
	//final clean up
	printf("There are %d best regions and %d errornous regions.\n", saved.size(), errornous.size());
	saved.insert(saved.end(), errornous.begin(), errornous.end());

	for (int i = 0; i < saved.size(); ++i)
	{
		saved[i]->bKeep = true;
	}
	for (int i = 0; i < allRegions.size(); ++i)
	{
		if (allRegions[i]->bKeep == false)
		{
			delete allRegions[i];
		}
	}

	return saved;
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

	vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	set<CoreParticle*> dots;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) dots.insert(mp[i]);
	}
	vector<CoreParticle*> dotsv = setupParticleNeighbors(mp, ndim, dims);
	vector<int> S(nvoxels, 0);
	vector<RegionInfo*> regions = grow(S, mp, dots, numKeep, numIter, bEightNeighborhood, ndim, dims);
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
			for (set<CoreParticle*>::iterator it=r->pnts.begin(); it != r->pnts.end(); ++it, j++)
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
			printf("%d Region: %d %d %d %f %f %f\n",
				i, r->gen, r->pnts.size(), r->contour.size(), r->measure, r->area, r->perimeter);
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
