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
#include <szConnectedComponent.h>
#include <CoreParticle.h>
#include <CoreParticleUtility.h>
#include <NeighborhoodFactory.h>
#include <CoreParticleMakeGraph.h>

//using namespace cv;
using namespace std;

struct DotComponent
{
	static int _id;
	DotComponent()
	{
		id = _id++;
	}
	DotComponent(vector<CoreParticle*> dots) : DotComponent()
	{
		this->dots = dots;
	}
	vector<CoreParticle*> dots;
	int id;
};
int DotComponent::_id = 0;

struct DotComponentPair
{
	DotComponentPair(DotComponent* a, DotComponent* b, int gen)
	{
		this->a = a;
		this->b = b;
		this->gen = gen;
		this->merged = NULL;
		value = 0.0f;
	}
	DotComponent* a;
	DotComponent* b;
	DotComponent* merged;
	int gen; //generation where the merge happened.
	float value; 
};

vector<DotComponent*>
initialComponents(
	vector<CoreParticle*>& mp,
	vector<CoreParticle*>& dots,
	bool bEightNeighorhood,
	int ndim, const int* dims)
{
	vector<int> S(numberOfElements(ndim, dims));
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) S[i] = 1;
	}
	vector<int> C(S.size());
	int nc;
	if (bEightNeighorhood)
	{
		nc = ConnectedComponentAnalysisBigger(C, S, NeighborhoodEight, 0, ndim, dims);
	}
	else
	{
		nc = ConnectedComponentAnalysisBigger(C, S, NeighborhoodFour, 0, ndim, dims);
	}
	vector<DotComponent*> comps(nc);
	for (int n = 0; n < nc; ++n)
	{
		comps[n] = new DotComponent();
	}
	for (int i = 0; i < dots.size(); ++i)
	{
		int val = GetVoxel(C, dots[i], 0, ndim, dims);
		comps[val - 1]->dots.push_back(dots[i]);
	}
	return comps;
}


/*
This routine grows regions from a set of dots (in DOTS) until all the dots are connected.
MP gives the spatial arrangement of particles.

The routine also establish a hierarchical structure as initial dot (connected) components are merged.
*/
vector<DotComponent*> 
grow(vector<int>& S,
	vector<CoreParticle*>& mp,
	vector<DotComponent*>& comps0,
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
	vector<DotComponent*> comps = comps0;
	set<CoreParticle*> front;
	for (int i = 0; i < comps.size(); ++i)
	{
		front.insert(comps[i]->dots.begin(), comps[i]->dots.end());
	}
	int iter = 1;
	vector<DotComponentPair> pairs;
	while (!front.empty())
	{
		vector<int> loc(ndim);
		set<CoreParticle*> next;
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			SetVoxel(S, p, iter, ndim, dims);
			p->gen = iter;
		}
		vector<int> C(S.size());
		int nc;
		if (bEightNeighorhood)
		{
			nc = ConnectedComponentAnalysisBigger(C, S, NeighborhoodEight, 0, ndim, dims);
		}
		else
		{
			nc = ConnectedComponentAnalysisBigger(C, S, NeighborhoodFour, 0, ndim, dims);
		}
		if (nc != comps.size())
		{
			//check if any DotComponents are merged.
			vector<TK::Node<DotComponent*>*> nodes;
			map<DotComponent*, TK::Node<DotComponent*>*> nmap;
			for (int i = 0; i < comps.size(); ++i)
			{
				TK::Node<DotComponent*>* n = TK::makeset(comps[i]);
				nodes.push_back(n);
				nmap[comps[i]] = n;
			}
			int nc0 = comps.size();
			vector<DotComponentPair> newpairs;
			for (int i = 0; i < nc0; ++i)
			{
				int lb1 = GetVoxel(C, *(comps[i]->dots.begin()), 0, ndim, dims);
				for (int j = i + 1; j < nc0; ++j)
				{
					int lb2 = GetVoxel(C, *(comps[j]->dots.begin()), 0, ndim, dims);
					if (lb1 == lb2)
					{
						DotComponent* a = TK::findset(nodes[i])->key;
						DotComponent* b = TK::findset(nodes[j])->key;
						if (a != b)
						{
							DotComponent* c = new DotComponent();
							c->dots = a->dots;
							c->dots.insert(c->dots.end(), b->dots.begin(), b->dots.end());
							DotComponentPair pr(a, b, iter);
							pr.merged = c;
							newpairs.push_back(pr);
							TK::Node<DotComponent*>* n = TK::makeset(c);
							TK::mergeFirst(n, nmap[comps[i]]);
							TK::mergeFirst(n, nmap[comps[j]]);
							nodes.push_back(n);
							nmap[c] = n;
							/*printf("Generation %d: %d(%d)=>%d(%d) + %d(%d)=>%d(%d) -> %d (%d)\n",
								pr.gen,
								comps[i]->id, comps[i]->dots.size(),
								pr.a->id, pr.a->dots.size(),
								comps[j]->id, comps[j]->dots.size(),
								pr.b->id, pr.b->dots.size(),
								pr.merged ? pr.merged->id : -1, pr.merged ? pr.merged->dots.size() : 0);
								*/
						}
					}
				}
			}
			vector<TK::Node<DotComponent*>*> reps = TK::clusters(nodes);
			comps.clear();
			for (int i = 0; i < reps.size(); ++i)
			{
				comps.push_back(reps[i]->key);
			}
			pairs.insert(pairs.end(), newpairs.begin(), newpairs.end());
			for (int i = 0; i < nodes.size(); ++i)
			{
				if (TK::findset(nodes[i]) != nodes[i]) //not a representative
				{
					delete nodes[i]->key;
				}
			}
			for (int i = 0; i < nodes.size(); ++i)
			{
				delete nodes[i];
			}
		}
		if (nc <= 1) break;

		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			vector<int> sub = coreParticle2Index(p, ndim);
			for (int n = 0; n < nbh.size(); ++n)
			{
				if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
				{
					int idx = Sub2Ind(sub, nbh[n], ndim, dims);
					if (mp[idx] == NULL)
					{
						for (int k = 0; k < ndim; ++k)
						{
							loc[k] = sub[k] + nbh[n][k];
						}
						CoreParticle* q = coreParticleNdim(loc, ndim);
						SetVoxel(mp, q, q, ndim, dims);
						next.insert(q);
						q->gen = iter;
					}
				}
			}
		}
		front = next;
		iter++;
	}
	/*for (int i = 0; i < pairs.size(); ++i)
	{
		printf("Generation %d: %d(%d) + %d(%d) -> %d(%d)\n",
			pairs[i].gen, 
			pairs[i].a->id, pairs[i].a->dots.size(), 
			pairs[i].b->id, pairs[i].b->dots.size(),
			pairs[i].merged->id, pairs[i].merged->dots.size());
	}
	for (int i = 0; i < allComps.size(); ++i)
	{
		delete allComps[i];
	}*/
	return comps;
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
	bool bEightNeighborhood;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		int val = 0;
		ReadScalar(val, prhs[1], classMode);
		if (val) bEightNeighborhood = true;
	}
	int maxIter = 64;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(maxIter, prhs[2], classMode);
	}
	int nvoxels = numberOfElements(ndim, dims);

	vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	vector<CoreParticle*> dots;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) dots.push_back(mp[i]);
	}
	vector<int> S(nvoxels, 0);
	vector<DotComponent*> comps = initialComponents(mp, dots, bEightNeighborhood, ndim, dims);
	comps = grow(S, mp, comps, bEightNeighborhood, ndim, dims);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndim, dims);
	vector<int> S2(nvoxels, 0);
	vector<CoreParticle*> cores = propagateParticles(particles, mp, S2, ndim, dims);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, ndim, dims);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreData(S2, mxINT32_CLASS, ndim, dims);
	}	
	
	/*if (nlhs >= 2)
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
	*/
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}
