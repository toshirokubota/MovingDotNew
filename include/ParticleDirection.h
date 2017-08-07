#pragma once

struct ParticleDirection {
	float x;
	float y;
	float z;
	float length;
	ParticleDirection(float x = 0, float y = 0, float z = 0,  bool norm=false) 
	{
		_init(x, y, z, norm);
	}
	ParticleDirection(CParticleF& from, CParticleF& to, bool norm = false){
		_init(to.m_X - from.m_X, to.m_Y - from.m_Y, to.m_Z - from.m_Z, norm);
	}
private:
	void _init(float x, float y, float z, bool norm)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		length = sqrt(x*x + y*y + z*z);
		if (norm)
		{
			this->x /= length;
			this->y /= length;
			this->z /= length;
		}
	}
};
