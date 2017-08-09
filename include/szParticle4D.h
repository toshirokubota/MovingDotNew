#ifndef ___SZ_CPARTICLE_4D_H___
#define ___SZ_CPARTICLE_4D_H___

class CParticle4D
{
public:
	CParticle4D(int x=0, int y=0, int z=0, int t=0)
	{
		m_X = x;
		m_Y = y;
		m_Z = z;
		m_T = t;
	}
	bool operator < (const CParticle4D& p) const
	{
		return m_T < p.m_T || (m_T == p.m_T && m_Z < p.m_Z) || (m_T == p.m_T && m_Z == p.m_Z && m_Y < p.m_Y) || (m_T == p.m_T && m_Z == p.m_Z && m_Y == p.m_Y && m_X < p.m_X);
	}

public:
	int m_X;
	int m_Y;
	int m_Z;
	int m_T;
};

#endif /* ___SZ_CPARTICLE_4D_H___ */