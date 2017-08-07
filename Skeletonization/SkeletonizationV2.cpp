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
#include <szConvexHull2D.h>
#include <SimplePoint.h>
int SimplePoint::_id = 0;

vector<int>
initialClustering(vector<SimplePoint*>& P, float space)
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
			if (Distance(P[i]->x, P[i]->y, P[j]->x, P[j]->y) <= space)
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
vector<vector<SimplePoint*>> divideClusters(vector<SimplePoint*>& P, vector<int>& labels, int minSize)
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
	vector<vector<SimplePoint*>> bins(numC + 1);
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
Derive a 2D convex hull from the clustre of features.
*/
vector<SimplePoint*>
convexify(vector<SimplePoint*>& cluster)
{
	vector<CParticleF> points;
	for (int i = 0; i < cluster.size(); ++i)
	{
		points.push_back(CParticleF(cluster[i]->x, cluster[i]->y));
	}
	vector<CParticleF> h = ConvexHull2D(points);
	vector<SimplePoint*> res;
	for (int i = 0; i < h.size(); ++i)
	{
		for (int j = 0; j < cluster.size(); ++j)
		{
			if (cluster[j]->x == h[i].m_X && cluster[j]->y == h[i].m_Y)
			{
				res.push_back(cluster[j]);
				break;
			}
		}
	}
	return res;
}
 
/*
From a convex hull, pick a pair that is farthest apart. 
*/
pair<SimplePoint*, SimplePoint*>
selectExtermity(vector<SimplePoint*>& hull)
{
	float dmax = 0;
	pair<SimplePoint*, SimplePoint*> extremity;
	for (int i = 0; i < hull.size(); ++i)
	{
		for (int j = i + 1; j < hull.size(); ++j)
		{
			float d = Distance(hull[i]->x, hull[i]->y, hull[j]->x, hull[j]->y);
			if (d > dmax)
			{
				dmax = d;
				extremity.first = hull[i];
				extremity.second = hull[j];
			}
		}
	}
	return extremity;
}



/*
Find shortest paths from each convex vertex to the source.
*/
vector<Edge<SimplePoint*>*>
skeletonTree(vector<SimplePoint*>& P, vector<SimplePoint*>& hull, pair<SimplePoint*, SimplePoint*>& extremity, 
			float space, 
			float minprop = 0.25 //the length of branche has to be at least this proportion of the spine.
			)
{
	GraphFactory<SimplePoint*>& factory = GraphFactory<SimplePoint*>::GetInstance();
	vector<Vertex<SimplePoint*>*> vertices;
	Vertex<SimplePoint*>* src = NULL;
	Vertex<SimplePoint*>* sink = NULL;
	vector<Vertex<SimplePoint*>*> vhull;
	
	//first find a shortest path between the two extremity points.
	for (int i = 0; i < P.size(); ++i)
	{
		Vertex<SimplePoint*>* v = factory.makeVertex(P[i]);
		if (P[i] == extremity.first)
		{
			src = v;
		}
		else if (P[i] == extremity.second)
		{
			sink = v;
		}
		else
		{
			if (find(hull.begin(), hull.end(), P[i]) != hull.end())
			{
				vhull.push_back(v);
			}
		}
		vertices.push_back(v);
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		for (int j = i; j < vertices.size(); ++j)
		{
			float d = Distance(P[i]->x, P[i]->y, P[j]->x, P[j]->y);
			if (d <= space)
			{
				Edge<SimplePoint*>* e1 = factory.makeEdge(vertices[i], vertices[j], d*d);
				vertices[i]->Add(e1);
				Edge<SimplePoint*>* e2 = factory.makeEdge(vertices[j], vertices[i], d*d);
				vertices[j]->Add(e2);
			}
		}
	}
	Dijkstra(vertices, src);
	vector<Vertex<SimplePoint*>*> spine = tracePath(sink, src);
	//now set weights of the path between the source and sink to 0
	//also add the path to the skeleton edge set, and to the skeleton vertex set.
	set<Vertex<SimplePoint*>*> vset;
	for (int i = 0; i < spine.size(); ++i)
	{
		vset.insert(spine[i]);
	}
	set<Edge<SimplePoint*>*> eset;
	for (int i = 0; i < spine.size() - 1; ++i)
	{
		vset.insert(spine[i]);
		Edge<SimplePoint*>* ed = spine[i]->findEdge(spine[i + 1]);
		ed->w = 0.0f;
		eset.insert(ed);
		Edge<SimplePoint*>* ed2 = spine[i + 1]->findEdge(spine[i]);
		ed2->w = 0.0f;
	}
	//now run dijkstra from each convex hull vertex
	vector<vector<Vertex<SimplePoint*>*>> traces;
	vector<pair<float, int>> pairs;
	for (int i = 0; i < vhull.size(); ++i)
	{
		Dijkstra(vertices, vhull[i]);
		//from vhull[i], keep the trace until it meets the long spine of the skeleton
		vector<Vertex<SimplePoint*>*> tr = tracePath(sink, vhull[i]);
		for (int j = 0; j < tr.size(); ++j)
		{
			if (vset.find(tr[j]) != vset.end())
			{
				pairs.push_back(pair<float, int>(tr[j]->d, i));
				break;
			}
		}
		traces.push_back(tr);
	}
	sort(pairs.begin(), pairs.end());
	//go from the longest to shortest.
	//since it is too tedious to store weight sum for all intermediate nodes in each trace,
	//we only use the number of vertices as approximation to the trace length
	for (int i = pairs.size() - 1; i >= 0; i--)
	{
		int id = pairs[i].second;
		vector<Vertex<SimplePoint*>*> tr;
		for (int j = 0; j < traces[id].size(); ++j)
		{
			if (vset.find(traces[id][j]) != vset.end())
			{
				break;
			}
			tr.push_back(traces[id][j]);
		}
		if ((float)tr.size() >(float)spine.size() * minprop)
		{
			for (int j = 0; j < tr.size(); ++j)
			{
				vset.insert(tr[j]);
			}
			for (int j = 0; j < tr.size() - 1; ++j)
			{
				eset.insert(tr[j]->findEdge(tr[j + 1]));
			}
		}
	}

	vector<Edge<SimplePoint*>*> skeleton;
	for (set<Edge<SimplePoint*>*>::iterator it = eset.begin(); it != eset.end(); ++it)
	{
		skeleton.push_back(*it);
	}
	return skeleton;
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
	vector<SimplePoint*> P;
	SimplePointFactory& factory = SimplePointFactory::getInstance();
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
			P.push_back(factory.makePont(x, y));
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
	float minProportion = 0.25f;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(minProportion, prhs[3], classMode);
	}

	vector<int> labels = initialClustering(P, space);
	vector<vector<SimplePoint*>> clusters = divideClusters(P, labels, minSize);

	vector<vector<pair<SimplePoint*,SimplePoint*>>> skeletons;
	for (int i = 0; i < clusters.size(); ++i)
	{
		vector<SimplePoint*> hull = convexify(clusters[i]);
		pair<SimplePoint*, SimplePoint*> extremity = selectExtermity(hull);
		vector<Edge<SimplePoint*>*> t = skeletonTree(clusters[i], hull, extremity, space, minProportion);

		vector<pair<SimplePoint*,SimplePoint*>> vf;
		for (int k = 0; k < t.size(); ++k)
		{
			vf.push_back(pair<SimplePoint*, SimplePoint*>(t[k]->u->key, t[k]->v->key));
		}
		skeletons.push_back(vf);
	}

	if (nlhs >= 1)
	{
		const int dims[] = { skeletons.size(), 1 };
		vector<vector<float>> vvfl;
		for (int i = 0; i < skeletons.size(); ++i)
		{
			const int dims2[] = { skeletons[i].size(), 4 };
			vector<float> vfl(dims2[0] * dims2[1]);
			for (int j = 0; j < dims2[0]; ++j)
			{
				SetData2(vfl, j, 0, dims2[0], dims2[1], skeletons[i][j].first->x);
				SetData2(vfl, j, 1, dims2[0], dims2[1], skeletons[i][j].first->y);
				SetData2(vfl, j, 2, dims2[0], dims2[1], skeletons[i][j].second->x);
				SetData2(vfl, j, 3, dims2[0], dims2[1], skeletons[i][j].second->y);
			}
			vvfl.push_back(vfl);
		}
		plhs[0] = StoreDataCell(vvfl, mxSINGLE_CLASS, 2, dims, 4);
	}
	if (nlhs >= 2) //initial clusters
	{
		const int dims[] = { P.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], P[i]->x);
			SetData2(F, i, 1, dims[0], dims[1], P[i]->y);
			SetData2(F, i, 2, dims[0], dims[1], (float)labels[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	GraphFactory<SimplePoint*>::GetInstance().Clean();
	SimplePointFactory::getInstance().clean();

	mexUnlock();
}

