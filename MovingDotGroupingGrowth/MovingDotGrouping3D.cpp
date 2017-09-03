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
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <CoreParticle.h>
#include <CoreParticleMakeGraph.h>
#include <szConnectedComponentC.h>
#include <CoreParticle.h>
#include <CoreParticleUtility.h>

//using namespace cv;
using namespace std;
int xoff[] = { 0, -1, 1, 0, 0, 0, -1, 1, -1, 1, -1, 0, 1, -1, 1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 0, 1 };
int yoff[] = { -1, 0, 0, 1, 0, 0, -1, -1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1 };
int zoff[] = { 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1 };

void
grow2(vector<int>& S,
	vector<CoreParticle*>& D,
	int numNeighbors, int maxIter,
	int ndim, const int* dims)
{
	int nvoxels = numberOfElements(ndim, dims);
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();

	set<CoreParticle*> front;
	for (int i = 0; i < D.size(); ++i)
	{
		CoreParticle* p = D[i];
		if (p == NULL) continue;
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1)
			continue;

		for (int k = 0; k < numNeighbors; ++k)
		{
			int i2 = p->x + zoff[k];
			int j2 = p->y + yoff[k];
			CoreParticle* q = NULL;
			q = GetData2(D, i2, j2, dims[0], dims[1], (CoreParticle*)NULL);
			if (q == NULL)
			{
				front.insert(p);
				break;
			}
		}
	}
	int iter = 0;
	while (!front.empty() && iter <= maxIter)
	{
		set<CoreParticle*> next;
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			SetData2(S, p->x, p->y, dims[0], dims[1], iter);

			for (int k = 0; k < numNeighbors; ++k)
			{
				int j2 = p->y + yoff[k];
				int k2 = p->x + xoff[k];
				if (j2 < 0 || k2 < 0 || j2 >= dims[1] || k2 >= dims[0]) continue;

				CoreParticle* q = GetData2(D, k2, j2, dims[0], dims[1], (CoreParticle*)NULL);
				if (q == NULL)
				{
					q = factory.makeParticle(k2, j2);
					q->gen = iter + 1; //save the state in the particle as well.
					SetData2(D, k2, j2, dims[0], dims[1], q);
					next.insert(q);
				}
			}
		}
		front = next;
		iter++;
	}
}

void
grow3(vector<int>& S, 
	vector<CoreParticle*>& D, 
	int numNeighbors, int maxIter, 
	int ndim, const int* dims)
{
	int nvoxels = numberOfElements(ndim, dims);
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();

	set<CoreParticle*> front;
	for (int i = 0; i < D.size(); ++i)
	{
		CoreParticle* p = D[i];
		if (p == NULL) continue;
		if (p->x == 0 || p->y == 0 || p->x == dims[0] - 1 || p->y == dims[1] - 1 || 
			p->z == 0 || p->z == dims[2] - 1)
			continue;

		for (int k = 0; k < numNeighbors; ++k)
		{
			int i2 = p->z + zoff[k];
			int j2 = p->y + yoff[k];
			int k2 = p->x + xoff[k];
			CoreParticle* q = NULL;
			q = GetData3(D, k2, j2, i2, dims[0], dims[1], dims[2], (CoreParticle*)NULL);
			if (q == NULL)
			{
				front.insert(p);
				break;
			}
		}
	}
	int iter = 0;
	while (!front.empty() && iter <= maxIter)
	{
		set<CoreParticle*> next;
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			SetData3(S, p->x, p->y, p->z, dims[0], dims[1], dims[2], iter);
			p->gen = iter; //save the state in the particle as well.

			for (int k = 0; k < numNeighbors; ++k)
			{
				int i2 = p->z + zoff[k];
				int j2 = p->y + yoff[k];
				int k2 = p->x + xoff[k];
				if (i2 < 0 || j2 < 0 || k2 < 0 || i2 >= dims[2] || j2 >= dims[1] || k2 >= dims[0]) continue;

				CoreParticle* q = GetData3(D, k2, j2, i2, dims[0], dims[1], dims[2], (CoreParticle*)NULL);
				if (q == NULL)
				{
					q = factory.makeParticle(k2, j2, i2);
					SetData3(D, k2, j2, i2, dims[0], dims[1], dims[2], q);
					next.insert(q);
				}
			}
		}
		front = next;
		iter++;
	}
}

void
grow(vector<int>& S,
	vector<CoreParticle*>& D,
	int numNeighbors, int maxIter,
	int ndim, const int* dims)
{
	if (ndim == 2)
	{
		grow2(S, D, numNeighbors, maxIter, ndim, dims);
	}
	else if (ndim == 3)
	{
		grow3(S, D, numNeighbors, maxIter, ndim, dims);
	}
}

/*
Establish ascendent/decendent relationships among particles based on the
region growth pattern. The iteration number of the growth is stored in gen instance variable
of CoreParticle.
*/
void
linkParticles(vector<CoreParticle*>& particles,
	bool bEightNeighborhood,
	int ndim, const int* dims)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		float minDLow = std::numeric_limits<float>::infinity();
		float minDHigh = std::numeric_limits<float>::infinity();
		vector<CoreParticle*> vHigh;
		vector<CoreParticle*> vLow;
		vector<CoreParticle*> nbrs = bEightNeighborhood ? p->neighbors8 : p->neighbors;
		for (int j = 0; j < nbrs.size(); ++j)
		{
			CoreParticle* q = nbrs[j];
			if (p->gen < q->gen)
			{
				float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
				if (d < minDHigh)
				{
					minDHigh = d;
					vHigh.clear();
					vHigh.push_back(q);
				}
				else if (d == minDHigh)
				{
					vHigh.push_back(q);
				}
			}
			else if (p->gen > q->gen)
			{
				float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
				if (d < minDLow)
				{
					minDLow = d;
					vLow.clear();
					vLow.push_back(q);
				}
				else if (d == minDLow)
				{
					vLow.push_back(q);
				}
			}
		}
		for(int j=0; j<vHigh.size(); ++j)
		{
			p->descendents.insert(vHigh[j]);
			vHigh[j]->ascendents.insert(p);
		}
		for (int j = 0; j<vLow.size(); ++j)
		{
			p->ascendents.insert(vLow[j]);
			vLow[j]->descendents.insert(p);
		}
	}
}

/*
Apply connected component on particles with the same growth number.
*/
vector<vector<CoreParticle*>>
connectedGrowthComponents(vector<CoreParticle*>& particles,
	bool bEightNeighborhood,
	int ndim, const int* dims)
{
	vector<TK::Node<CoreParticle*>*> nodes;
	map<CoreParticle*, TK::Node<CoreParticle*>*> pmap;
	for (int i = 0; i < particles.size(); ++i)
	{
		TK::Node<CoreParticle*>* n = TK::makeset(particles[i]);
		nodes.push_back(n);
		pmap[n->key] = n;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		vector<CoreParticle*> nbrs = bEightNeighborhood ? p->neighbors8 : p->neighbors;
		for (int j = 0; j < nbrs.size(); ++j)
		{
			CoreParticle* q = nbrs[j];
			if (p->gen == q->gen)
			{
				TK::merge(nodes[i], pmap[q]);
			}
		}
	}
	vector<TK::Node<CoreParticle*>*> reps = TK::clusters(nodes);
	map<TK::Node<CoreParticle*>*, int> nrmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		nrmap[reps[i]] = i;
	}
	vector<vector<CoreParticle*>> groups(reps.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = nrmap[TK::findset(nodes[i])];
		groups[k].push_back(nodes[i]->key);
		nodes[i]->key->label = k + 1;  //label it with the cluster index
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return groups;
}

vector<vector<CoreParticle*>>
localMaximumComponents(vector<CoreParticle*>& particles,
	bool bEightNeighborhood,
	int ndim, const int* dims)
{
	vector<vector<CoreParticle*>> groups = connectedGrowthComponents(particles, bEightNeighborhood, ndim, dims);
	vector<int> lmaxIdx(groups.size(), true);
	vector<vector<CoreParticle*>> lmaxComponents;
	for (int i = 0; i < groups.size(); ++i)
	{
		for (int j = 0; j < groups[i].size() && lmaxIdx[i]; ++j)
		{
			CoreParticle* p = groups[i][j];
			vector<CoreParticle*> nbrs = bEightNeighborhood ? p->neighbors8 : p->neighbors;
			for (int k = 0; k < nbrs.size(); ++k)
			{
				CoreParticle* q = nbrs[k];
				if (p->gen < q->gen)
				{
					lmaxIdx[i] = false;
					break;
				}
			}
		}
		if (lmaxIdx[i])
		{
			vector<CoreParticle*> comp;
			for (int j = 0; j < groups[i].size(); ++j)
			{
				comp.push_back(groups[i][j]);
			}
			lmaxComponents.push_back(comp);
		}
	}
	return lmaxComponents;
}

/*
Based on the growth pattern (captured in ascendent/descendent relations), pickup core particles
and origins of the growth. 
Use the following numeric encoding:
core = 0: regular particle
core = 1: origin particle
core = 2: core particle
core = 3: strong core particle -- those in a locally maximum component
core = 4: growth boundary (where the growth stopped as it reached the maximum iteration).
*/
void
classifyParticles(vector<CoreParticle*>& particles,
				bool bEightNeighborhood,
				int maxIter,
				int ndim, const int* dims)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		p->core = 0;
		if (p->descendents.empty())
		{
			//check if no descendent because of premature ending (by the maximum iteration setting).
			if (p->gen < maxIter && p->neighbors8.size() == (int)pow(3,ndim)-1)
			{
				p->core = 2;
			}
		}
		if (p->ascendents.empty())
		{
			p->core = 1;
		}
	}
	vector<vector<CoreParticle*>> lmaxC = localMaximumComponents(particles, bEightNeighborhood, ndim, dims);
	for (int i = 0; i < lmaxC.size(); ++i)
	{
		for (int j = 0; j < lmaxC[i].size(); ++j)
		{
			CoreParticle* p = lmaxC[i][j];
			if (p->gen < maxIter && p->neighbors8.size() == (int)pow(3, ndim) - 1)
			{
				p->core = 3;
			}
		}
	}
}

/*
Track the particle DAG from strong cores, and color the trace in T.
*/
void
reverseGrow2(vector<int>& T,
	vector<CoreParticle*>& particles,
	bool bEightNeighborhood,
	int ndim, const int* dims)
{
	set<int> sgen;
	int count = 0;
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
		if (particles[i]->core == 2 || particles[i]->core == 3)
		//if (particles[i]->core == 3)
		{
			sgen.insert(particles[i]->gen);
			count++;
		}
	}
	printf("There are %d strong cores.\n", count);
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (int i = 0; i<particles.size(); ++i)
		{
			if (particles[i]->gen == gval && (particles[i]->core == 2 || particles[i]->core == 3))
			//if (particles[i]->gen == gval && particles[i]->core == 3)
			{
				Q.push_back(particles[i]);
			}
		}
		while (Q.empty() == false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				SetVoxel(T, p, gval, ndim, dims);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
					{
						q->label = p->label;
						Q2.insert(q);
					}
				}
			}
			Q.clear();
			Q.insert(Q.end(), Q2.begin(), Q2.end());
		}
	}
}

void
reverseGrow(vector<int>& T,
	vector<CoreParticle*>& particles,
	bool bEightNeighborhood,
	int ndim, const int* dims)
{
	vector<vector<CoreParticle*>> comps = localMaximumComponents(particles, bEightNeighborhood, ndim, dims);
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
	}
	set<int> sgen;
	for (int i = 0; i<comps.size(); ++i)
	{
		CoreParticle* p = comps[i][0];
		sgen.insert(p->gen);
		for (int j = 0; j < comps[i].size(); ++j)
		{
			comps[i][j]->label = i + 1; //unique label per component
		}
	}
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (int i = 0; i<comps.size(); ++i)
		{
			if (comps[i][0]->gen == gval)
			{
				for (int j = 0; j < comps[i].size(); ++j)
				{
					//if (comps[i][j]->core == 3)
					{
						Q.push_back(comps[i][j]);
					}
				}
				//Q.insert(Q.end(), comps[i].begin(), comps[i].end());
			}
		}
		while (Q.empty() == false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				SetVoxel(T, p, p->label, ndim, dims);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
					{
						q->label = p->label;
						Q2.insert(q);
					}
				}
			}
			Q.clear();
			Q.insert(Q.end(), Q2.begin(), Q2.end());
		}
	}
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
	int numNeighbors = 6;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		int val = 0;
		ReadScalar(val, prhs[1], classMode);
		if (val) numNeighbors = 26;
	}
	int maxIter = 64;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(maxIter, prhs[2], classMode);
	}
	int nvoxels = numberOfElements(ndim, dims);

	vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	vector<int> S(nvoxels, 0);
	grow(S, mp, numNeighbors, maxIter, ndim, dims);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndim, dims);
	linkParticles(particles, numNeighbors > 6, ndim, dims);
	classifyParticles(particles, numNeighbors > 6, maxIter, ndim, dims);
	vector<int> S2(nvoxels, 0);
	reverseGrow2(S2, particles, numNeighbors > 6, ndim, dims);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, ndim, dims);
	}
	if (nlhs >= 2)
	{
		vector<unsigned char> C(S.size(), 0);
		for (int i = 0; i < particles.size(); ++i)
		{
			if (particles[i]->core)
			{
				SetVoxel(C, particles[i], (unsigned char)particles[i]->core, ndim, dims);
			}
		}
		plhs[1] = StoreData(C, mxUINT8_CLASS, ndim, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(S2, mxINT32_CLASS, ndim, dims);
	}

	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}
