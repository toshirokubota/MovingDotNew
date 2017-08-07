// TopologicalSort.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <limits>
#include <algorithm>
using namespace std;
#include <Graph.h>


/*
A modified DFSVisit.
As vertices are visited, they are stored in a stack.
The topological sort can simply pop them out of the stack to
order them in the decreasing order of the finish time.
*/
template <class T>
void
__DFSVisitTS(Vertex<T>* node,
		stack<Vertex<T>*>& st,
		 int& time)
{
	time++;
	node->color = Gray;
	node->d = time;
	for(int i=0; i<node->aList.size(); ++i)
	{
		Vertex<T>* v = node->aList[i]->v;
		if(v->color == White)
		{
			v->pi = node;
			__DFSVisitTS(v, st, time);
		}
	}
	node->color = Black;
	node->f = ++time;
	st.push(node);
}

/*
The function takes a set of vertices in a directed graph, 
sort them in a topological order, and returned the same set of vertices
in the sorted order.
*/
template <class T>
vector<Vertex<T>*>
TopoSort(vector<Vertex<T>*> vnodes)
{
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i]->Reset();
	}
	vector<Vertex<T>*> q;
	int time = 0;
	stack<Vertex<T>*> st;
	for(int i=0; i<vnodes.size(); ++i)
	{
		if(vnodes[i]->color == White)
		{
			__DFSVisitTS(vnodes[i], st, time);
		}
	}
	while(st.empty() == false)
	{
		q.push_back(st.top());
		st.pop();
	}

	return q;
}

