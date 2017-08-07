#ifndef ___GRAPH_FACTORY_H___
#define ___GRAPH_FACTORY_H___

#include <Graph.h>
#include <Tree.h>

template <class T>
class GraphFactory
{
public:
	static GraphFactory& GetInstance()
	{
		static GraphFactory* _instance = NULL;
		if (_instance == NULL)
		{
			_instance = new GraphFactory();
		}
		return *_instance;
	}
	Vertex<T>* makeVertex(T key)
	{
		Vertex<T>* u = new Vertex<T>(key);
		vertices.push_back(u);
		return u;
	}
	TreeNode<T>* makeTreeNode(T key)
	{
		TreeNode<T>* u = new TreeNode<T>(key);
		vertices.push_back(u);
		return u;
	}
	Edge<T>* makeEdge(Vertex<T>* u, Vertex<T>* v, float w=0, EdgeType type = Tr)
	{
		Edge<T>* e = new Edge<T>(u, v, w, type);
		edges.push_back(e);
		return e;
	}
	FlowEdge<T>* makeFlowEdge(Vertex<T>* u, Vertex<T>* v, float w=0, EdgeType type = Tr)
	{
		FlowEdge<T>* e = new FlowEdge<T>(u, v, w, type);
		edges.push_back(e);
		return e;
	}
	bool RemoveVertex(Vertex<T>* u)
	{
		std::vector<Vertex<T>*>::iterator p = find(vertices.begin(), vertices.end(), u);
		if(p == vertices.end())
		{
			return false;
		}
		else
		{
			vertices.erase(p);
			return true;
		}
	}
	bool RemoveEdge(Edge<T>* e)
	{
		std::vector<Edge<T>*>::iterator p = find(edges.begin(), edges.end(), u);
		if(p == edges.end())
		{
			return false;
		}
		else
		{
			edges.erase(p);
			return true;
		}
	}
	void Clean()
	{
		for(int i=0; i<edges.size(); ++i)
		{
			delete edges[i];
			edges[i] = 0;
		}
		edges.clear();
		for(int i=0;  i<vertices.size(); ++i)
		{
			delete vertices[i];
			vertices[i] = 0;
		}
		vertices.clear();
	}
	void Clean(Vertex<T>* u) //remove this vertex from the graph
	{
		//remove all edges incident from u
		for (int i = 0; i < u->aList.size(); ++i)
		{
			Edge<T>* e = u->aList[i];
			std::vector<Edge<T>*>::iterator it = find(edges.begin(), edges.end(), e);
			edges.erase(it);
			delete e;
		}
		//remove all edges incident to u
		for (int i = edges.size() - 1; i >= 0; --i)
		{
			Edge<T>* e = edges[i];
			if (e->v == u)
			{
				e->u->aList.erase(find(e->u->aList.begin(), e->u->aList.end(), e));
				edges.erase(edges.begin() + i);
				delete e;
			}
		}
		//finally delete the vertex
		std::vector<Vertex<T>*>::iterator it = find(vertices.begin(), vertices.end(), u);
		vertices.erase(it);
		delete u;
	}
	std::vector<Vertex<T>*> vertices;
	std::vector<Edge<T>*> edges;
private:
	GraphFactory() {};
	~GraphFactory()
	{
		Clean();
	}
	GraphFactory(GraphFactory& f)
	{
	}
	GraphFactory operator=(GraphFactory& f)
	{
	}
};

#endif /*  ___GRAPH_FACTORY_H___ */