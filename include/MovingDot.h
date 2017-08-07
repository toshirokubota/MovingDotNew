#pragma once

class MovingDot {
public:
	MovingDot(int x, int y, int t)
	{
		this->x = x;
		this->y = y;
		this->t = t;
	}
	int x;
	int y;
	int t;
};