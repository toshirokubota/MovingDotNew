#pragma once
#include <vector>
#include <limits>
#include <cmatrix>
using namespace std;
//a generic feature vector

class Feature
{
public:
	Feature(int n=0, float val=0.0f)
	{
		vals = vector<float>(n, val);
	}
	Feature(vector<float> v)
	{
		vals = v;
	}

	static float distance(const Feature& a, const Feature& b)
	{
		if (a.length() != b.length())
		{
			return std::numeric_limits<float>::quiet_NaN();
		}
		float sum = 0;
		for (int i = 0; i < a.length(); ++i)
		{
			float df = (a.vals[i] - b.vals[i]);
			sum += df * df;
		}
		return sqrt(sum);
	}

	static Feature computeMean(const vector<Feature>& a)
	{
		if (a.size()==0)
		{
			return Feature(0);
		}
		int K = a[0].length();
		Feature mu(K);
		for (int i = 0; i < a.size(); ++i)
		{
			for (int j = 0; j < K; ++j)
			{
				mu.vals[j] += a[i].vals[j];
			}
		}
		for (int i = 0; i < K; ++i)
		{
			mu.vals[i] /= a.size();
		}
		return mu;
	}

	bool operator == (Feature& f)
	{
		if (f.vals.size() != vals.size()) return false;
		for (int k = 0; k < vals.size(); ++k)
		{
			if (vals[k] != f.vals[k]) return false;
		}
		return true;
	}

	static techsoft::matrix<float> toColumnVector(const Feature& f)
	{
		techsoft::matrix<float> v(f.length(), 1);
		for (int i = 0; i < f.length(); ++i)
		{
			v[i][0] = f.vals[i];
		}
		return v;
	}

	int length() const { return vals.size(); }
	vector<float> vals;
};

