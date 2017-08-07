#include <Snapshot.h>
#include <szmexutilitytemplate.h>
#include <mexFileIO.h>

vector<Snapshot>
Snapshot::TakeSnapshot(float time)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<Snapshot> polygons;
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory.activeSet.begin(); it != factory.activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = MovingParticle::vectorize(p);
				Snapshot shot(time, time, vp);
				bdone = false;
				polygons.push_back(shot);
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
				}
			}
		}
		if (bdone) break;
	}
	return polygons;
}

mxArray*
Snapshot::StoreSnapshot(Snapshot& snapshot)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<CParticleF> shape = snapshot.polygon->project(snapshot.getProjectionTime());
	const int dims[] = { shape.size(), 4 };
	vector<float> F(dims[0] * dims[1]);
	for (int j = 0; j < dims[0]; ++j)
	{
		SetData2(F, j, 0, dims[0], dims[1], shape[j].m_X);
		SetData2(F, j, 1, dims[0], dims[1], shape[j].m_Y);
		SetData2(F, j, 2, dims[0], dims[1], snapshot.projection_time);
		SetData2(F, j, 3, dims[0], dims[1], snapshot.created_time);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

/*
Store stationary particles.
*/
mxArray*
Snapshot::StoreSnapshot0(Snapshot& snapshot)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<MovingParticle*> vp = snapshot.polygon->getParticles();
	const int dims[] = { vp.size(), 3 };
	vector<float> F(dims[0] * dims[1]);
	for (int j = 0; j < dims[0]; ++j)
	{
		SetData2(F, j, 0, dims[0], dims[1], vp[j]->getInitParticle()->getX());
		SetData2(F, j, 1, dims[0], dims[1], vp[j]->getInitParticle()->getY());
		SetData2(F, j, 2, dims[0], dims[1], (float)vp[j]->getInitParticle()->getId());
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

mxArray*
Snapshot::StoreSnapshots(vector<Snapshot>& snapshots)
{
	mxArray* plhs;
	const int dims[] = { snapshots.size(), 1 };
	plhs = mxCreateCellArray((mwSize)2, (const mwSize*)dims);
	for (int i = 0; i<snapshots.size(); ++i)
	{
		mxArray* a = StoreSnapshot(snapshots[i]);
		mxSetCell(plhs, i, a);
	}
	return plhs;
}

mxArray*
Snapshot::StoreSnapshots0(vector<Snapshot>& snapshots)
{
	mxArray* plhs;
	const int dims[] = { snapshots.size(), 1 };
	plhs = mxCreateCellArray((mwSize)2, (const mwSize*)dims);
	for (int i = 0; i<snapshots.size(); ++i)
	{
		mxArray* a = StoreSnapshot0(snapshots[i]);
		mxSetCell(plhs, i, a);
	}
	return plhs;
}

