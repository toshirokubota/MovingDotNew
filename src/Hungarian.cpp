#include <Hungarian.h>
//#include <mex.h>

/*
This version of Hungarian algorithm is by Andrei Lopatin. See http://e-maxx.ru/algo/assignment_hungary.
*/
float Hungarian::solve(bool bMaximize) {
	float INF = std::numeric_limits<float>::infinity();
	int m = nJobs;
	int n = nWorkers;
	//store the cost matrix as 1 based index.
	vector<vector<float>> a(m + 1);
	a[0] = vector<float>(n + 1, 0);
	for (int i = 1; i <= m; ++i)
	{
		vector<float> a0(n + 1, 0);
		for (int j = 1; j <= n; ++j)
		{
			a0[j] = C[(i - 1)* n + j - 1];
			if (bMaximize)
			{
				a0[j] = -a0[j];
			}
		}
		a[i] = a0;
	}

	vector<float> u(n + 1);
	vector<float> v(m + 1);
	vector<int> p(m + 1);
	vector<int> way(m + 1);
	for (int i = 1; i <= n; ++i) {
		p[0] = i;
		int j0 = 0;
		vector<float> minv(m + 1, INF);
		vector<char> used(m + 1, false);
		do {
			used[j0] = true;
			int i0 = p[j0];
			float delta = INF;
			int j1;
			for (int j = 1; j <= m; ++j)
			{
				if (!used[j]) 
				{
					float cur = a[i0][j] - u[i0] - v[j];
					if (cur < minv[j])
						minv[j] = cur, way[j] = j0;
					if (minv[j] < delta)
					{
						delta = minv[j], j1 = j;
						//printf("delta = %f\n", delta);
					}
				}
			}
			for (int j = 0; j <= m; ++j)
			{
				if (used[j])
					u[p[j]] += delta, v[j] -= delta;
				else
					minv[j] -= delta;
			}
			j0 = j1;
		} while (p[j0] != 0);
		do {
			int j1 = way[j0];
			p[j0] = p[j1];
			j0 = j1;
		} while (j0);
	}

	for (int j = 1; j <= nJobs; ++j)
		jobs[p[j] - 1] = j - 1;
	for (int j = 1; j <= nWorkers; ++j)
		workers[j-1] = p[j]-1;

	return bMaximize ?  v[0]: -v[0];
}
