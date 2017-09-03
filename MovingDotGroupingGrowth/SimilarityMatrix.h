#ifndef ___SIMILARITY_MATRIX_H____
#define ___SIMILARITY_MATRIX_H____
#include <opencv2/core/core.hpp>
//#include <cstdlib>

struct SimilarityMatrix {
	SimilarityMatrix(int n)
	{
		m = cv::Mat::zeros(n, n, CV_32F);
	}
	int size() {
		return m.rows;
	}
	bool set(int row, int col, float val)
	{
		if (row >= m.rows || col >= m.cols) return false;
		m.at<float>(row, col) = val; //make it symmetric
		m.at<float>(col, row) = val;
	}
	float maxVal(int& row, int& col)
	{
		float maxval = 0;
		for (int i = 0; i < m.rows; ++i)
		{
			for (int j = i + 1; j < m.rows; ++j)
			{
				if (maxval < m.at<float>(i, j))
				{
					maxval = m.at<float>(i, j);
					row = i;
					col = j;
				}
			}
		}
		return maxval;
	}
	float minVal(int& row, int& col)
	{
		float minval = std::numeric_limits<float>::infinity();
		for (int i = 0; i < m.rows; ++i)
		{
			for (int j = i + 1; j < m.rows; ++j)
			{
				if (minval > m.at<float>(i, j))
				{
					minval = m.at<float>(i, j);
					row = i;
					col = j;
				}
			}
		}
		return minval;
	}
	void remove(int idx)
	{
		int n = m.rows;
		cv::Mat m2 = m.clone();
		m = cv::Mat(n - 1, n - 1, m.type());
		for (int i = 0, i2 = 0; i < n; ++i)
		{
			if (i == idx) continue;
			for (int j = 0, j2 = 0; j < n; ++j)
			{
				if (j == idx) continue;

				m.at<float>(i2, j2) = m2.at<float>(i, j);
				j2++;
			}
			i2++;
		}
	}

	cv::Mat m;
};

#endif /* ___SIMILARITY_MATRIX_H____ */
