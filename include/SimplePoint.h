#ifndef ___SIMPLE_POINT_H___
#define ___SIMPLE_POINT_H___
#include <vector>
using namespace std;
struct SimplePoint {
	SimplePoint(float x0 = 0, float y0 = 0, float z0 = 0, float f=0, int g=0)
	{
		x = x0;
		y = y0;
		z = z0;
		frame = f;
		gid = g;
		id = _id++;
	}
	float x;
	float y;
	float z;
	int id;
	float frame; //for spatio-temporal point.
	int gid;  //group id for debugging.
	static int _id;
};

struct SimplePointFactory
{
	static SimplePointFactory& getInstance()
	{
		static SimplePointFactory* _instance = NULL;
		if (_instance == NULL)
		{
			_instance = new SimplePointFactory();
		}
		return *_instance;

	}
	SimplePoint* makePont(float x=0, float y=0, float z=0, float f=0, int g=0)
	{
		SimplePoint* p = new SimplePoint(x, y, z, f, g);
		points.push_back(p);
		return p;
	}
	void clean()
	{
		for (int i = 0; i<points.size(); ++i)
		{
			delete points[i];
		}
		points.clear();
	}
	vector<SimplePoint*> points;
private:
	SimplePointFactory()
	{
		SimplePoint::_id = 0;
	}
	SimplePointFactory(SimplePointFactory& f)
	{
	}
	SimplePointFactory operator=(SimplePointFactory& f)
	{
	}
	~SimplePointFactory()
	{
		clean();
	}
};

#endif /* ___SIMPLE_POINT_H___ */