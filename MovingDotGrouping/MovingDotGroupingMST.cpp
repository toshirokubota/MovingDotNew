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
#include <Kruskal.h>
//using namespace cv;
using namespace std;

vector<TK::Edge<CoreParticle*>*>
constructMST(vector<CoreParticle*>& particles, float dthres)
{
	TK::GraphFactory<CoreParticle*>& gfactory = TK::GraphFactory<CoreParticle*>::GetInstance();
	vector<TK::Node<CoreParticle*>*> nodes;
	map<CoreParticle*, int> imap;
	vector<TK::Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < particles.size(); ++i)
	{
		TK::Node<CoreParticle*>* n = TK::makeset(particles[i]);
		imap[particles[i]] = i;
		nodes.push_back(n);
		vertices.push_back(gfactory.makeVertex(particles[i]));
	}
	vector<TK::Edge<CoreParticle*>*> edges;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		for (int j = 0; j < p->neighbors.size(); ++j)
		{
			CoreParticle* q = p->neighbors[j];
			int k = imap[q];
			TK::merge(nodes[i], nodes[k]);
			float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
			TK::Edge<CoreParticle*>* ed = gfactory.makeEdge(vertices[i], vertices[k], d);
			vertices[i]->aList.push_back(ed);
			edges.push_back(ed);
		}
	}
	printf("There are %d neighbor edges.\n", edges.size());
	vector<vector<TK::Node<CoreParticle*>*>> groups = TK::grouping(nodes);
	printf("There are %d disjoint regions.\n", groups.size());
	for (int i = 0; i < groups.size(); ++i)
	{
		for (int j = i + 1; j < groups.size(); ++j)
		{
			pair<TK::Node<CoreParticle*>*, TK::Node<CoreParticle*>*> pr(NULL,NULL);
			float mind = dthres;
			for (int i2 = 0; i2 < groups[i].size(); ++i2)
			{
				CoreParticle* p = groups[i][i2]->key;
				for (int j2 = 0; j2 < groups[j].size(); ++j2)
				{
					CoreParticle* q = groups[j][j2]->key;
					float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
					if (d < mind)
					{
						mind = d;
						pr.first = groups[i][i2];
						pr.second = groups[j][j2];
					}
				}
			}
			if (pr.first && pr.second)
			{
				int k1 = imap[pr.first->key];
				int k2 = imap[pr.second->key];
				TK::Edge<CoreParticle*>* ed1 = gfactory.makeEdge(vertices[k1], vertices[k2], mind);
				TK::Edge<CoreParticle*>* ed2 = gfactory.makeEdge(vertices[k2], vertices[k1], mind);
				vertices[k1]->Add(ed1);
				vertices[k2]->Add(ed2);
			}
		}
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	vector<TK::Edge<CoreParticle*>*> mst = Kruskal(vertices);

	return mst;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: T = MovingDotGrouping(P)");
		return;
	}
	int ndim;
	const int* dims;
	mxClassID classL;
	vector<unsigned char> L;
	LoadData(L, prhs[0], classL, ndim, &dims);

	float dist_thres = 10.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(dist_thres, prhs[1], classMode);
	}

	TK::GraphFactory<CoreParticle*>& gfactory = TK::GraphFactory<CoreParticle*>::GetInstance();
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> mp = generateParticleMap(L, ndim, dims);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndim, dims);
	printf("There are %d particles.\n", particles.size());
	vector<TK::Edge<CoreParticle*>*> mst = constructMST(particles, dist_thres);

	if (nlhs >= 1)
	{
		int dims[] = { particles.size(), 5 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->z);
			SetData2(F, i, 3, dims[0], dims[1], p->t);
			SetData2(F, i, 4, dims[0], dims[1], p->id);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { mst.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < mst.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], mst[i]->u->key->id);
			SetData2(F, i, 1, dims[0], dims[1], mst[i]->v->key->id);
			SetData2(F, i, 2, dims[0], dims[1], (int)(mst[i]->w));
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	gfactory.Clean();
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}
