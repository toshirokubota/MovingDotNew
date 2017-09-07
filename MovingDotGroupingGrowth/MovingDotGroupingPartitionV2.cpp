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
This routine grows regions from a set of dots (in DOTS) for a fixed number of iterations.
MP gives the spatial arrangement of particles.

The routine also establish a hierarchical structure as initial dot (connected) components are merged.
*/
set<CoreParticle*>
grow(vector<int>& S,
	vector<CoreParticle*>& mp,
	set<CoreParticle*>& particles,
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
		p->gen = 0;
		p->src = p;
	}

	int iter = 1;
	while (!front.empty() && iter <= niter)
	{
		vector<int> loc(ndim);
		set<CoreParticle*> next;
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			SetVoxel(S, p, iter, ndim, dims);
			p->gen = iter;
		}

		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			vector<int> sub = coreParticle2Index(p, ndim);
			for (int n = 0; n < nbh.size(); ++n)
			{
				if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
				{
					int idx = Sub2Ind(sub, nbh[n], ndim, dims);
					if (S[idx] == 0)
					{
						CoreParticle* q = mp[idx];
						if (q == NULL)
						{
							for (int k = 0; k < ndim; ++k)
							{
								loc[k] = sub[k] + nbh[n][k];
							}
							q = coreParticleNdim(loc, ndim);
							SetVoxel(mp, q, q, ndim, dims);
							next.insert(q);
							q->gen = iter;
							q->src = p->src;
							particles.insert(q);
						}
						q->ascendents.insert(p);
						p->descendents.insert(q);
					}
				}
			}
		}
		front = next;
		iter++;
	}
	return particles;
}

/*
Label particles based on asendency-descendency from cores.
*/
void
partitionParticles(set<CoreParticle*>& dots,
	vector<int>& S,
	int maxIter,
	bool bEightNeighborhood,
	int ndim, const int* dims)
{
	vector<TK::Node<CoreParticle*>*> nodes;
	map<CoreParticle*, int> imap;
	for (set<CoreParticle*>::iterator it = dots.begin(); it != dots.end(); ++it)
	{
		CoreParticle* p = *it;
		if (p->descendents.empty() && p->gen < maxIter)
		{
			TK::Node<CoreParticle*>* n = TK::makeset(p);
			imap[p] = nodes.size();
			nodes.push_back(n);
			p->core = 1;
		}
		else
		{
			p->core = 0;
		}
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		/*set<CoreParticle*> nbh;
		if (bEightNeighborhood)
		{
			nbh.insert(p->neighbors8.begin(), p->neighbors8.end());
		}
		else
		{
			nbh.insert(p->neighbors.begin(), p->neighbors.end());
		}*/
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			CoreParticle* q = nodes[j]->key;
			if(length(p->x-q->x, p->y-q->y, p->z-q->z, p->t-q->t) < 3.0f)
			//if (nbh.find(q) != nbh.end())
			{
				TK::merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<vector<TK::Node<CoreParticle*>*>> groups = TK::grouping(nodes);
	for (int i = 0; i < groups.size(); ++i)
	{
		for (int j = 0; j < groups[i].size(); ++j)
		{
			groups[i][j]->key->label = i + 1;
		}
	}
	printf("There are %d cores into %d groups.\n", nodes.size(), groups.size());
	set<CoreParticle*> Q;
	for (int i = 0; i < nodes.size(); ++i)
	{
		Q.insert(nodes[i]->key);
	}
	while (!Q.empty())
	{
		set<CoreParticle*> Q2;
		for (set<CoreParticle*>::iterator it = Q.begin(); it != Q.end(); ++it)
		{
			CoreParticle* p = *it;
			SetVoxel(S, p, p->label, ndim, dims);
			for (set<CoreParticle*>::iterator jt = p->ascendents.begin(); jt != p->ascendents.end(); ++jt)
			{
				CoreParticle* q = *jt;
				q->label = p->label;
				Q2.insert(q);
			}
		}
		Q = Q2;
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
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
	bool bEightNeighborhood;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		int val = 0;
		ReadScalar(val, prhs[1], classMode);
		if (val) bEightNeighborhood = true;
	}
	int numIter = 30;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[2], classMode);
	}
	int nvoxels = numberOfElements(ndim, dims);

	vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	set<CoreParticle*> dots;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) dots.insert(mp[i]);
	}
	vector<int> S(nvoxels, 0);
	//vector<DotComponent*> comps = initialComponents(mp, dots, bEightNeighborhood, ndim, dims);
	printf("There are %d particles before growth.\n", dots.size());
	dots = grow(S, mp, dots, numIter, bEightNeighborhood, ndim, dims);
	printf("There are %d particles after growth.\n", dots.size());
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndim, dims);
	vector<int> S2(nvoxels, 0);
	partitionParticles(dots, S2, numIter, bEightNeighborhood, ndim, dims);

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
