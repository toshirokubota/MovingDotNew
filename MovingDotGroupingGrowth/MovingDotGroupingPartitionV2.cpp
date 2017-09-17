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
#include <szConvexHull2D.h>
//using namespace cv;
using namespace std;

struct RegionInfo
{
	RegionInfo(vector<CoreParticle*>& pnts, float m, int g, int n)
	{
		contour = pnts;
		measure = m;
		gen = g;
		size = n;
	}
	vector<CoreParticle*> contour;
	float measure;
	int gen;
	int size;
};

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
This routine grows regions from a set of dots (in DOTS) for a fixed number of iterations.
MP gives the spatial arrangement of particles.

The routine also establish a hierarchical structure as initial dot (connected) components are merged.
*/
vector<RegionInfo>
grow(vector<int>& S,
	vector<CoreParticle*>& mp,
	set<CoreParticle*>& particles,
	float thres,
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
		TK::Node<CoreParticle*>* n = TK::makeset(p);
	}
	vector<int> L0;
	vector<RegionInfo> regions;
	vector<RegionInfo> saved;
	int nc0 = 0;
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
		vector<int> L(numberOfElements(ndim, dims), 0);
		int nc = ConnectedComponentAnalysisBigger(L, S, NeighborhoodEight, 0, ndim, dims);
		vector<RegionInfo> regions0 = regions;
		regions = vector<RegionInfo>();
		for (int k = 1; k <= nc; ++k)
		{
			set<CoreParticle*> sp;
			vector<CoreParticle*> vp;
			for (int m = 0; m < L.size(); ++m)
			{
				if (L[m] == k)
				{
					//if (mp[m]->gen == 1)
					{
						sp.insert(mp[m]);
					}
					if (mp[m]->gen == iter)
					{
						vp.push_back(mp[m]);
					}
				}
			}
			regions.push_back(RegionInfo(vp, shapeMeasure2(vp,sp), iter, sp.size()));
		}

		if (nc != nc0)
		{
			//check which regions have been merged.
			for (int i = 0; i < regions0.size(); ++i)
			{
				int label1 = GetVoxel(L, *(regions0[i].contour.begin()), 0, ndim, dims);
				for (int j = 0; j < regions0.size(); ++j)
				{
					if (i == j) continue;
					int label2 = GetVoxel(L, *(regions0[j].contour.begin()), 0, ndim, dims);
					if (label1 == label2)
					{
						printf("Merge %d %d: %d,%d,%f  => %d,%d,%f, %c\n",
							iter, i,
							regions0[i].size, regions0[i].contour.size(), regions0[i].measure,
							regions[label1 - 1].size, regions[label1 - 1].contour.size(), 
							regions[label1 - 1].measure,
							regions0[i].measure>regions[label1-1].measure ? '*' : 'x');
						if (regions0[i].measure > thres && 
							regions0[i].measure > regions[label1 - 1].measure)
						{
							saved.push_back(regions0[i]);
						}
						break;
					}
				}
			}
			/*
			//compare the measure.
			for (int i = 0; i < pairs.size(); ++i)
			{
				float m1 = regions0[pairs[i].first].measure;
				float m2 = regions0[pairs[i].second].measure;
				int label = GetVoxel(L, *(regions0[pairs[i].first].contour.begin()), 0, ndim, dims);
				float m3 = regions[label - 1].measure;
				printf("Merge %d %d: %d + %d => %d, %f * %f => %f, %c, %c\n",
					iter, i,
					regions0[pairs[i].first].contour.size(),
					regions0[pairs[i].second].contour.size(),
					regions[label - 1].contour.size(), 
					m1, m2, m3,
					m1*m2>m3*m3 && m1 > thres ? '*': 'x',
					m1*m2>m3*m3 && m2 > thres ? '*': 'x');
				if (m1*m2 > m3*m3)
				{
					//saved.push_back(regions[label - 1]);
					if (m1 > thres)
					{
						saved.push_back(regions0[pairs[i].first]);
					}
					if (m2 > thres)
					{
						saved.push_back(regions0[pairs[i].second]);
					}
				}
			}*/
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
		for (set<CoreParticle*>::iterator it = front.begin(); it != front.end(); ++it)
		{
			CoreParticle* p = *it;
			strongMedialParticle(p, ndim, dims);
		}
		front = next;
		iter++;
		nc0 = nc;
	}
	return saved;
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
	float thres =  0.5;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[3], classMode);
	}
	int nvoxels = numberOfElements(ndim, dims);

	vector<CoreParticle*> mp = generateParticleMap(im, ndim, dims);
	set<CoreParticle*> dots;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i]) dots.insert(mp[i]);
	}
	vector<CoreParticle*> dotsv = setupParticleNeighbors(mp, ndim, dims);
	vector<int> S(nvoxels, 0);
	vector<RegionInfo> regions = grow(S, mp, dots, thres, numIter, bEightNeighborhood, ndim, dims);
	vector<int> S2(nvoxels, 0);
	//partitionParticles(dots, S2, numIter, bEightNeighborhood, ndim, dims);

	if (nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, ndim, dims);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreData(S2, mxINT32_CLASS, ndim, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { regions.size(), 1 };
		vector<vector<int>> F; // (dims[0] * dims[1]);
		int ncol = 4;
		for(int i=0; i<regions.size(); ++i)
		{
			const int dims2[] = { regions[i].contour.size(), ncol };
			vector<int> f(dims2[0] * dims2[1]);
			for (int j = 0; j < regions[i].contour.size(); ++j)
			{
				CoreParticle* p = regions[i].contour[j];
				SetData2(f, j, 0, dims2[0], dims2[1], p->x);
				SetData2(f, j, 1, dims2[0], dims2[1], p->y);
				SetData2(f, j, 2, dims2[0], dims2[1], regions[i].gen);
				SetData2(f, j, 3, dims2[0], dims2[1], (int)(regions[i].measure*1000.0));
			}
			F.push_back(f);
		}
		plhs[2] = StoreDataCell(F, mxINT32_CLASS, ndim, dims, ncol);
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
