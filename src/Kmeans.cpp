#include <Kmeans.h>
#include <szMiscOperations.h>

int
selectKmeansCluster(Feature& p, vector<Feature>& centers)
{
	if (centers.size() == 0) return -1;

	float mind = Feature::distance(p, centers[0]);
	int sel = 0;
	for (int i = 1; i < centers.size(); ++i)
	{
		float d = Feature::distance(p, centers[i]);
		if (d < mind)
		{
			mind = d;
			sel = i;
		}
	}
	return sel;
}

vector<Feature>
KmeansClustering(vector<Feature>& points, int K, float tolerance)
{
	if (K >= points.size())
	{
		return points;
	}
	int N = points[0].length();

	Feature zero(N);
	vector<Feature> centers(K, zero);
	vector<int> ri = randomIndices(points.size(), K);
	for (int i = 0; i < K; ++i)
	{
		centers[i] = points[ri[i]];
	}
	int maxIter = 1000;
	int iter = 0;
	while (true)
	{
		iter++;
		if (iter >= maxIter)
		{
			printf("KMeans: Maximum iteration reached before convergence.\n");
			break;
		}
		vector<Feature> centers2(K, zero);
		vector<int> counts(K, 0);
		for (int i = 0; i < points.size(); ++i)
		{
			int c = selectKmeansCluster(points[i], centers);
			for (int j = 0; j < N; ++j)
			{
				centers2[c].vals[j] += points[i].vals[j];
			}
			counts[c]++;
		}
		bool bDone = true;
		//remove empty cluster
		for (int i = centers.size() - 1; i >= 0; --i)
		{
			if (counts[i] == 0)
			{
				centers.erase(centers.begin() + i);
				counts.erase(counts.begin() + i);
				centers2.erase(centers2.begin() + i);
				printf("KMeans: Removed empty cluster. K is now %d.\n", centers.size());
				bDone = false;
				K--;
			}
		}
		//new centroids
		for (int i = 0; i < centers.size(); ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				centers2[i].vals[j] = centers2[i].vals[j] / counts[i];
			}
		}
		//check for convergence
		float sum = 0;
		for (int i = 0; i < centers.size(); ++i)
		{
			sum += Feature::distance(centers[i], centers2[i]);
		}
		if (sum > tolerance)
		{
			bDone = false;
		}
		centers = centers2;

		if (bDone) break;
	}

	return centers;
}
