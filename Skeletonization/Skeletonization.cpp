#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>
#include <map>
#include <set>
#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szParticleF.h>
#include <Feature.h>
#include <DisjointSet.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Dijkstra.h>

vector<int>
initialClustering(vector<Feature>& P, float space)
{
	vector<Node<int>*> nodes;
	for (int i = 0; i < P.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			if (Feature::distance(P[i], P[j]) <= space)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<int>*> c = clusters(nodes);
	vector<int> labels(nodes.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = distance(c.begin(), find(c.begin(), c.end(), findset(nodes[i])));
		labels[i] = k;
	}
	return labels;
}

/*
Divide the feature set (P) into multiple bins, according to the labeling (LABEL).
Small clusters less than (MINSIZE) are not kept.
*/
vector<vector<Feature>> divideClusters(vector<Feature>& P, vector<int>& labels, int minSize)
{
	int numC = 0;
	for (int i = 0; i < labels.size(); ++i)
	{
		numC = Max(labels[i], numC);
	}
	vector<int> counts(numC+1, 0);
	for (int i = 0; i < labels.size(); ++i)
	{
		counts[labels[i]]++;
	}
	vector<vector<Feature>> bins(numC + 1);
	for (int lb = 0; lb <= numC; ++lb)
	{
		if (counts[lb] >= minSize)
		{
			for (int j = 0; j < labels.size(); ++j)
			{
				if (labels[j] == lb)
				{
					bins[lb].push_back(P[j]);
				}
			}
		}
	}
	for (int i = bins.size() - 1; i >= 0; i--)
	{
		if (bins[i].empty())
		{
			bins.erase(bins.begin() + i);
		}
	}
	return bins;
}

/*
Find a skeleton as a sequence of points that are on a longest shortest path between a pair of points.
Collect skeleton from a collection of points.
*/
vector<vector<Feature>>
skeleton(vector<Feature>& P, float space, float minDistance, float prop)
{
	GraphFactory<int>& factory = GraphFactory<int>::GetInstance();
	vector<Vertex<int>*> vertices;
	for (int i = 0; i < P.size(); ++i)
	{
		vertices.push_back(factory.makeVertex(i));
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		for (int j = i; j < vertices.size(); ++j)
		{
			float d = Feature::distance(P[i], P[j]);
			if (d <= space)
			{
				Edge<int>* e1 = factory.makeEdge(vertices[i], vertices[j], d);
				vertices[i]->Add(e1);
				Edge<int>* e2 = factory.makeEdge(vertices[j], vertices[i], d);
				vertices[j]->Add(e2);
			}
		}
	}
	vector<Vertex<int>*> uncovered = vertices;
	set<Vertex<int>*> covered;
	set<Vertex<int>*> traced;

	vector<vector<Feature>> paths;
	while (uncovered.empty() == false)
	{
		float dmax = 0;
		Vertex<int>* src = NULL;
		Vertex<int>* dst = NULL;
		//find the longest shortest path
		for (int i = 0; i < uncovered.size(); ++i)
		{
			Dijkstra(vertices, uncovered[i]);
			for (int j = 0; j < vertices.size(); ++j)
			{
				if (uncovered[i] == vertices[j]) continue;

				if (vertices[j]->d < std::numeric_limits<float>::infinity() && vertices[j]->d > dmax)
				{
					src = uncovered[i];
					dst = vertices[j];
					dmax = vertices[j]->d;
				}
			}
		}
		if (src != NULL)
		{
			//perform dijkstra one more time from the chosen src.
			Dijkstra(vertices, src);
			vector<Vertex<int>*> tr = tracePath(dst);
			vector<Feature> Q;
			for (int j = 0; j < tr.size(); ++j)
			{
				Q.push_back(P[tr[j]->key]);
				if (traced.find(tr[j]) != traced.end()) break;
				traced.insert(tr[j]);
			}
			paths.push_back(Q);

			//update covered and uncovered
			for (int j = uncovered.size() - 1; j >= 0; --j)
			{
				float dmin = std::numeric_limits<float>::infinity();
				Feature fj = P[uncovered[j]->key];
				bool bdiscovered = false;
				for (int k = 0; k < tr.size(); ++k)
				{
					Feature fk = P[tr[k]->key];
					if (Feature::distance(fj, fk) <= Max(minDistance, dmax*prop))
					{
						bdiscovered = true;
						break;
					}
				}
				if (bdiscovered)
				{
					covered.insert(uncovered[j]);
					uncovered.erase(uncovered.begin() + j);
				}
			}
		}
	}
	return paths;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: L = EMGClustering(P, K)");
		return;
	}
	//Points
	vector<Feature> P;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			Feature f(2);
			f.vals[0] = x;
			f.vals[1] = y;
			P.push_back(f);
		}
	}
	float space = 10.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(space, prhs[1], classMode);
	}
	int minSize = 2;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[2], classMode);
	}
	float minDistance = 10.0f;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(minDistance, prhs[3], classMode);
	}

	vector<int> labels = initialClustering(P, space);
	vector<vector<Feature>> clusters = divideClusters(P, labels, minSize);

	vector<vector<Feature>> clusters2;
	for (int i = 0; i < clusters.size(); ++i)
	{
		vector<vector<Feature>> t = skeleton(clusters[i], space, minDistance, 0.2f);
		clusters2.insert(clusters2.end(), t.begin(), t.end());
	}

	if (nlhs >= 1)
	{
		vector<Feature> A;
		vector<Feature> B;
		vector<int> segments;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size()-1; ++j)
			{
				A.push_back(clusters2[i][j]);
				B.push_back(clusters2[i][j+1]);
				segments.push_back(i);
			}
		}
		const int dims[] = { A.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], A[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], A[i].vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], B[i].vals[0]);
			SetData2(F, i, 3, dims[0], dims[1], B[i].vals[1]);
			SetData2(F, i, 4, dims[0], dims[1], (float) segments[i]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		vector<Feature> Q;
		vector<int> segments;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size(); ++j)
			{
				Q.push_back(clusters2[i][j]);
				segments.push_back(i);
			}
		}
		const int dims[] = { Q.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], Q[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], Q[i].vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], (float)segments[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		vector<int> labels2;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size(); ++j)
			{
				labels2.push_back(i);
			}
		}
		const int dims[] = { labels2.size(), 1 };
		plhs[2] = StoreData(labels2, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 4) //initial clusters
	{
		vector<Feature> A;
		vector<int> L;
		for (int i = 0; i < clusters.size(); ++i)
		{
			for (int j = 0; j < clusters[i].size() - 1; ++j)
			{
				A.push_back(clusters[i][j]);
				L.push_back(i);
			}
		}
		const int dims[] = { A.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], A[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], A[i].vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], (float)L[i]);
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

