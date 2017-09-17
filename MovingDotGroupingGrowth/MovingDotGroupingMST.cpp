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
#include <Kruskal.h>

//using namespace cv;
using namespace std;

vector<vector<TK::Edge<CoreParticle*>*>>
constructMSTForrest(vector<CoreParticle*>& particles, float maxGap)
{
	TK::GraphFactory<CoreParticle*>& factory = TK::GraphFactory<CoreParticle*>::GetInstance();
	vector<TK::Vertex<CoreParticle*>*> vertices;
	vector<TK::Node<CoreParticle*>*> nodes;
	map<CoreParticle*,int> imap;
	for (int i = 0; i < particles.size(); ++i)
	{
		vertices.push_back(factory.makeVertex(particles[i]));
		nodes.push_back(TK::makeset(particles[i]));
		imap[particles[i]] = i;
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		CoreParticle* p = vertices[i]->key;
		for (int j = i + 1; j < vertices.size(); ++j)
		{
			CoreParticle* q = vertices[j]->key;
			float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
			if (d <= maxGap)
			{
				vertices[i]->Add(factory.makeEdge(vertices[i], vertices[j], d*d));
				vertices[j]->Add(factory.makeEdge(vertices[j], vertices[i], d*d));
				TK::merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<TK::Edge<CoreParticle*>*> mst = Kruskal(vertices);
	vector<TK::Node<CoreParticle*>*> reps = TK::clusters(nodes);
	map<TK::Node<CoreParticle*>*, int> rmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		rmap[reps[i]] = i;
	}
	vector<vector<TK::Edge<CoreParticle*>*>> forrest(reps.size());
	for (int i = 0; i < mst.size(); ++i)
	{
		int j = imap[mst[i]->u->key];
		int k = rmap[nodes[j]];
		forrest[k].push_back(mst[i]);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return forrest;
}

vector<CoreParticle*> collectAscendents(CoreParticle* p)
{
	set<CoreParticle*> S;
	set<CoreParticle*> Q;
	Q.insert(p);
	while (!Q.empty())
	{
		set<CoreParticle*> Q2;
		for (set<CoreParticle*>::iterator it = Q.begin(); it != Q.end(); ++it)
		{
			CoreParticle* q = *it;
			if (q->ascendents.empty())
			{
				S.insert(q);
			}
			Q2.insert(q->ascendents.begin(), q->ascendents.end());
		}
		Q = Q2;
	}
	vector<CoreParticle*> vp;
	vp.insert(vp.end(), S.begin(), S.end());
	return vp;
}

/*
From each MST, perform region-growth.
Then, select cores and from each, derive a cycle connecting two disjoint particles whose
descendents meet at the core.
*/
vector<vector<CoreParticle*>> 
grow(vector<int>& S, 
	vector<vector<TK::Edge<CoreParticle*>*>>& forrest, 
	bool bEightNeighborhood, int maxIter, int minSize,
	int ndim, const int* dims)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	int xoff[] = { -1, 0, 1, -1, 1, -1, 0, 1 }; //only applicable to 2D for now.
	int yoff[] = { -1, -1, -1, 0, 0, 1, 1, 1 }; //only applicable to 2D for now.
	vector<vector<CoreParticle*>> paths;
	for (int fi = 0; fi < forrest.size(); ++fi)
	{
		vector<CoreParticle*> particles;

		if (forrest[fi].size() < minSize) continue;
		vector<int> S0(numberOfElements(ndim, dims), 0);
		vector<CoreParticle*> mp(numberOfElements(ndim, dims), NULL);
		set<CoreParticle*> Q;
		for (int i = 0; i < forrest[fi].size(); ++i)
		{
			CoreParticle* p = forrest[fi][i]->u->key;
			CoreParticle* q = forrest[fi][i]->v->key;
			Q.insert(p);
			Q.insert(q);
			SetVoxel(mp, p, p, ndim, dims);
			SetVoxel(mp, q, q, ndim, dims);
			particles.push_back(p);
			particles.push_back(q);
		}
		for (int iter = 1; iter <= maxIter && !Q.empty(); ++iter)
		{
			set<CoreParticle*> Q2;
			for (set<CoreParticle*>::iterator it = Q.begin(); it != Q.end(); ++it)
			{
				CoreParticle* p = *it;
				SetVoxel(S0, p, iter, ndim, dims);
			}
			for (set<CoreParticle*>::iterator it = Q.begin(); it != Q.end(); ++it)
			{
				CoreParticle* p = *it;
				for (int k = 0; k < sizeof(xoff) / sizeof(xoff[0]); ++k)
				{
					int x2 = p->x + xoff[k];
					int y2 = p->y + yoff[k];
					if (GetData2(mp, x2, y2, dims[0], dims[1], (CoreParticle*)NULL) == NULL)
					{
						CoreParticle* q = factory.makeParticle(x2, y2);
						SetVoxel(mp, q, q, ndim, dims);
						particles.push_back(q);
						Q2.insert(q);
					}
					CoreParticle* q = GetData2(mp, x2, y2, dims[0], dims[1], (CoreParticle*)NULL);
				}
				for (set<CoreParticle*>::iterator it = Q.begin(); it != Q.end(); ++it)
				{
					CoreParticle* p = *it;
					float minD = std::numeric_limits<float>::infinity();
					vector<CoreParticle*> vp;
					for (set<CoreParticle*>::iterator jt = Q2.begin(); jt != Q2.end(); ++jt)
					{
						CoreParticle* q = *jt;
						float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
						if (d < minD)
						{
							minD = d;
							vp.clear();
							vp.push_back(q);
						}
						else if (d == minD)
						{
							vp.push_back(q);
						}
					}
					for (int i = 0; i < vp.size(); ++i)
					{
						vp[i]->ascendents.insert(p);
						p->descendents.insert(vp[i]);
					}
				}
				for (set<CoreParticle*>::iterator it = Q2.begin(); it != Q2.end(); ++it)
				{
					CoreParticle* p = *it;
					float minD = std::numeric_limits<float>::infinity();
					vector<CoreParticle*> vp;
					for (set<CoreParticle*>::iterator jt = Q.begin(); jt != Q.end(); ++jt)
					{
						CoreParticle* q = *jt;
						float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
						if (d < minD)
						{
							minD = d;
							vp.clear();
							vp.push_back(q);
						}
						else if (d == minD)
						{
							vp.push_back(q);
						}
					}
					for (int i = 0; i < vp.size(); ++i)
					{
						vp[i]->descendents.insert(p);
						p->ascendents.insert(vp[i]);
					}
				}
			}
			Q = Q2;
		}
		vector<pair<CoreParticle*, CoreParticle*>> pairs;
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			if (p->descendents.empty() && p->ascendents.size()>=2)
			{
				vector<CoreParticle*> sa = collectAscendents(p);
				for (int j = 0; j < sa.size(); ++j)
				{
					for (int k = j + 1; k < sa.size(); ++k)
					{
						pairs.push_back(pair<CoreParticle*, CoreParticle*>(sa[j], sa[k]));
					}
				}
			}
		}
		vector<vector<CoreParticle*>> path0;
	}
	return paths;
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
	int maxIter = 16;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(maxIter, prhs[2], classMode);
	}
	float maxGap = 10;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(maxGap, prhs[3], classMode);
	}
	int minSize = 10;
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[4], classMode);
	}

	int nvoxels = numberOfElements(ndim, dims);

	vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	vector<CoreParticle*> dots;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) dots.push_back(mp[i]);
	}
	setupParticleNeighbors(mp, ndim, dims);

	vector<vector<TK::Edge<CoreParticle*>*>> forrest = constructMSTForrest(dots, maxGap);
	
	vector<int> S(nvoxels, 0);
	vector<vector<CoreParticle*>> paths = grow(S, forrest, bEightNeighborhood, maxIter, minSize, ndim, dims);

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
