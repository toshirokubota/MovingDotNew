#include <PolygonDAG.h>

PolygonDAG::PolygonDAG()
{
	bSorted = false;
}

void 
PolygonDAG::add(Polygon* u, Polygon* v)
{
	GraphFactory<Polygon*>& factory = GraphFactory<Polygon*>::GetInstance();
	Vertex<Polygon*>* from = NULL;
	Vertex<Polygon*>* to = NULL;
	if (pmap.find(u) == pmap.end())
	{
		from = factory.makeVertex(u);
		pmap[u] = from;
		vertices.push_back(from);
	}
	else
	{
		from = pmap[u];
	}
	if (pmap.find(v) == pmap.end())
	{
		to = factory.makeVertex(v);
		pmap[v] = to;
		vertices.push_back(to);
	}
	else
	{
		to = pmap[v];
	}
	bool bfound = false;
	for (int i = 0; i < from->aList.size(); ++i)
	{
		if (from->aList[i]->v == to)
		{
			bfound = true;
			break;
		}
	}
	if (bfound == false)
	{
		Edge<Polygon*>* edge = factory.makeEdge(from, to);
		from->aList.push_back(edge);
		edges.push_back(edge);
		bSorted = false;
	}
}

vector<Polygon*>
PolygonDAG::getPolygons() const
{
	vector<Polygon*> polygons;
	for (int i = 0; i < vertices.size(); ++i)
	{
		polygons.push_back(vertices[i]->key);
	}
	return polygons;
}

vector<Polygon*>
PolygonDAG::predecessor(Polygon* poly)
{
	vector<Polygon*> polygons;
	Vertex<Polygon*>* u = pmap[poly];
	if (u != NULL)
	{
		for (int i = 0; i < u->aList.size(); ++i)
		{
			polygons.push_back(u->aList[i]->v->key);
		}
	}
	return polygons;
}

void
_DFSVisit(Vertex<Polygon*>* node, int& time)
{
	time++;
	node->color = Gray;
	node->d = time;
	for (int i = 0; i<node->aList.size(); ++i)
	{
		Vertex<Polygon*>* v = node->aList[i]->v;
		if (v->color == White)
		{
			v->pi = node;
			_DFSVisit(v, time);
		}
	}
	node->color = Black;
	node->f = ++time;
}

void
_DFS(vector<Vertex<Polygon*>*> vnodes)
{
	for (int i = 0; i<vnodes.size(); ++i)
	{
		vnodes[i]->Reset();
	}

	int time = 0;
	for (int i = 0; i<vnodes.size(); ++i)
	{
		if (vnodes[i]->color == White)
		{
			_DFSVisit(vnodes[i], time);
		}
	}
}


void
PolygonDAG::sort()
{
	_DFS(vertices);
	vector<pair<int, Vertex<Polygon*>*>> pairs;
	for (int i = 0; i<vertices.size(); ++i)
	{
		pair<int, Vertex<Polygon*>*> p(vertices[i]->f, vertices[i]);
		pairs.push_back(p);
	}
	std::sort(pairs.begin(), pairs.end());
	vertices.clear();
	for (int i = pairs.size() - 1; i >= 0; --i)
	{
		Vertex<Polygon*>* p = pairs[i].second;
		vertices.push_back(p);
	}
	bSorted = true;
}
