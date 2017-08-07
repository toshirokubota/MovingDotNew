#include <EMGaussian.h>
#include <Kmeans.h>
#include <szMiscOperations.h>

/*
Select the Gaussian with the largest p.
*/
int
selectGaussianCluster(Feature& p, vector<GaussianDistribution>& gaussians)
{
	if (gaussians.size() == 0) return -1;

	float maxp = gaussians[0].eval(p);
	int sel = 0;
	for (int i = 1; i < gaussians.size(); ++i)
	{
		float d = gaussians[i].eval(p);
		if (d > maxp)
		{
			maxp = d;
			sel = i;
		}
	}
	return sel;
}

vector<GaussianDistribution>
EMGaussianClustering(vector<Feature>& points, int K, int maxIter, float tolerance)
{
	if (K >= points.size())
	{
		vector<GaussianDistribution> dist(points.size());
		for (int i = 0; i < points.size(); ++i)
		{
			dist[i].mu = points[i];
			vector<Feature> f(1, points[i]);
			dist[i].sgm = Covariance::computeCovariance(f);
		}
		return dist;
	}

	int N = points[0].length();
	Feature zero(N);

	//initialize the Gaussian distributions from K-means
	vector<Feature> centers = KmeansClustering(points, K, tolerance); //use K-means for initialization
	K = centers.size();
	Feature zero2(K);

	vector<int> labels(points.size());
	for (int i = 0; i < points.size(); ++i)
	{
		labels[i] = selectKmeansCluster(points[i], centers);
	}
	vector<GaussianDistribution> dist(centers.size());
	for (int i = 0; i < centers.size(); ++i)
	{
		vector<Feature> x;
		for (int j = 0; j < labels.size(); ++j)
		{
			if (labels[j] == i)
			{
				x.push_back(points[j]);
			}
		}
		dist[i] = GaussianDistribution(x);
	}

	int iter = 0;
	while (true)
	{
		iter++;
		if (iter >= maxIter)
		{
			//printf("EMGaussian: Maximum iteration reached before convergence.\n");
			break;
		}
		//E-step: compute the weights
		vector<Feature> w(points.size(), zero2);
		for (int i = 0; i < points.size(); ++i)
		{
			float sum = 0;
			for (int j = 0; j < K; ++j)
			{
				float pval = dist[j].eval(points[i]);
				w[i].vals[j] = pval;
				sum += pval;
			}
			for (int j = 0; j < K; ++j)
			{
				w[i].vals[j] /= sum;
			}
		}
		//M-step:
		//Compute the means
		vector<GaussianDistribution> dist2(K);
		vector<float> vsum(K);
		for (int i = 0; i < K; ++i)
		{
			Feature m = zero;
			float sum = 0;
			for (int j = 0; j < points.size(); ++j)
			{
				for (int k = 0; k < m.length(); ++k)
				{
					m.vals[k] += w[j].vals[i] * points[j].vals[k];
				}
				sum += w[j].vals[i];
			}
			for (int k = 0; k < m.length(); ++k)
			{
				m.vals[k] /= sum;
			}
			dist2[i].mu = m;
			vsum[i] = sum;
		}
		//Compute the covariances
		for (int i = 0; i < K; ++i)
		{
			techsoft::matrix<float> cov(N, N, 0.0f);
			techsoft::matrix<float> m = Feature::toColumnVector(dist2[i].mu);
			float sum = 0;
			for (int j = 0; j < points.size(); ++j)
			{
				techsoft::matrix<float> df(N, 1);
				for (int k = 0; k < N; ++k)
				{
					df(k,0) = points[j].vals[k] - m(k, 0);
				}
				cov = cov + w[j].vals[i] * (df * (~df));
			}
			cov /= vsum[i];
			dist2[i].sgm = Covariance(K);
			dist2[i].sgm.mat = cov;
			dist2[i].sgm.imat = !cov;
		}
		dist = dist2;
	}

	return dist;
}
