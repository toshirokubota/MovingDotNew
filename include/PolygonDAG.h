#pragma once
#include <Graph.h>
#include <GraphFactory.h>
#include <Polygon.h>
#include <vector>
#include <set>
#include <map>
using namespace std;

class PolygonDAG
{
public:
	PolygonDAG();
	void add(Polygon* u, Polygon* v);
	bool bSorted;
	void sort(); //topological sort
	vector<Polygon*> getPolygons() const;
	vector<Polygon*> predecessor(Polygon* poly); //

	vector<Edge<Polygon*>*> edges;
	map<Polygon*,Vertex<Polygon*>*> pmap;
	vector<Vertex<Polygon*>*> vertices;
};
