// Kruskal.cpp : Defines the entry point for the console application.
//

#include <vector>
#include <algorithm>
using namespace std;
#include <Graph.h>
#include <DisjointSet.h>

template<class T>
vector<Edge<T>*>
Kruskal(vector<Edge<T>*>& edges, 
		vector<Vertex<T>*>& vertices)
{
	vector<Edge<T>*> mst;
	//sortEdges(edges, 0, edges.size()-1);
	vector<pair<float, int>> pairs(edges.size());
	for (int i = 0; i < edges.size(); ++i)
	{
		pairs[i] = pair<float, int>(edges[i]->w, i);
	}
	sort(pairs.begin(), pairs.end());

	vector<Node<Vertex<T>*>*> nodes;
	map<Vertex<T>*, int> vmap;
	for(int i=0; i<vertices.size(); ++i)
	{
		nodes.push_back(makeset(vertices[i]));
		vmap[vertices[i]] = i;
	}

	for(int i=0; i<pairs.size(); ++i)
	{
		int j = pairs[i].second;
		int ui = vmap[edges[j]->u]; // distance(vertices.begin(), find(vertices.begin(), vertices.end(), edges[j]->u));
		int uj = vmap[edges[j]->v]; // distance(vertices.begin(), find(vertices.begin(), vertices.end(), edges[j]->v));
		if(findset(nodes[ui]) != findset(nodes[uj]))
		{
			merge(nodes[ui], nodes[uj]);
			mst.push_back(edges[j]);
		}
	}
	for(int i=0; i<vertices.size(); ++i)
	{
		delete nodes[i];
	}
	return mst;
}

#include <set>
template<class T>
vector<Edge<T>*>
Kruskal(vector<Edge<T>*>& edges)
{
	set<Vertex<T>*> V;
	for(int i=0; i<edges.size(); ++i)
	{
		V.insert(edges[i]->u);
		V.insert(edges[i]->v);
	}
	vector<Vertex<T>*> vertices;
	for(set<Vertex<T>*>::iterator p = V.begin(); p!=V.end(); p++)
	{
		vertices.push_back(*p);
	}
	return Kruskal(edges, vertices);
}

template<class T>
vector<Edge<T>*>
Kruskal(vector<Vertex<T>*>& vertices)
{
	vector<Edge<T>*> edges;
	for (int i = 0; i < vertices.size(); ++i)
	{
		edges.insert(edges.end(), vertices[i]->aList.begin(), vertices[i]->aList.end());
	}
	return Kruskal(edges, vertices);
}
