#pragma once;
#include <MovingParticle.h>
#include <mex.h>
#include <Snapshot.h>
#include <StationaryParticle.h>
#include <SkeletonDAG.h>

class OffsetPolygonDAGBuilder
{
public:
	OffsetPolygonDAGBuilder(const mxArray* points, const mxArray* edges)
	{
		dag = SkeletonDAG::fromMxArray(points, edges);
		for (int i = 0; i < dag.vertices.size(); ++i)
		{
			aMap[dag.vertices[i]->key] = true;
		}
	}
	virtual bool Polygonify();
	virtual bool Expand(float endtime);

	vector<vector<MovingParticle*>> forrest;
	map<StationaryParticle*,bool> aMap;
	SkeletonDAG dag;
protected:
	bool initializePolygon(vector<MovingParticle*>& particles);
};

