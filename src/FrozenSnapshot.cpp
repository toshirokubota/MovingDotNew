#include <FrozenSnapshot.h>
#include <szmexutilitytemplate.h>
#include <mexFileIO.h>


mxArray*
FrozenSnapshot::StoreSnapshot(FrozenSnapshot& snapshot)
{
	const int dims[] = { snapshot.particles.size(), 4 };
	vector<float> F(dims[0] * dims[1]);
	for (int j = 0; j < dims[0]; ++j)
	{
		SetData2(F, j, 0, dims[0], dims[1], snapshot.particles[j].m_X);
		SetData2(F, j, 1, dims[0], dims[1], snapshot.particles[j].m_Y);
		SetData2(F, j, 2, dims[0], dims[1], snapshot.particles[j].m_Z);
		SetData2(F, j, 3, dims[0], dims[1], snapshot.particles[j].m_Life);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

mxArray*
FrozenSnapshot::StoreSnapshots(vector<FrozenSnapshot>& snapshots)
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

