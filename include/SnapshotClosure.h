#pragma once
#include <Snapshot.h>

class SnapshotTrace : public Snapshot
{
public:
	SnapshotTrace(float t, float s, vector<MovingParticle*>& particles) : Snapshot(t, s, particles)
	{
		srcPolygon = NULL;
	}
	Polygon* srcPolygon;
};
