#include <cmath>
#include "interface.h"

void Interface::roe()
// Roe's Riemann solver for 4-element vectors
// Left and right values of the Riemann problem are defined by U1_[4] and U2_[4]
// Reference: Toro Chapter 11
{
	double h1, h2, h, p1, p2, u, v, a, theta, delta[4], dv[4];
	int i;
	
	p1 = (GAMMA - 1.0) * (U1_[3] - 0.5 * (U1_[1] * U1_[1] + U1_[2] * U1_[2]) / U1_[0]);
	p2 = (GAMMA - 1.0) * (U2_[3] - 0.5 * (U2_[1] * U2_[1] + U2_[2] * U2_[2]) / U2_[0]);
	h1 = (U1_[3] + p1) / U1_[0];
	h2 = (U2_[3] + p2) / U2_[0];
	
	theta = sqrt(U1_[0]) / (sqrt(U1_[0]) + sqrt(U2_[0]));
	u = theta * U1_[1] / U1_[0] + (1.0 - theta) * U2_[1] / U2_[0];
	v = theta * U1_[2] / U1_[0] + (1.0 - theta) * U2_[2] / U2_[0];
	h = theta * h1 + (1.0 - theta) * h2;
	a = sqrt((GAMMA - 1.0) * (h - 0.5 * u * u - 0.5 * v * v));
	
	for (i = 0; i < 4; i++)
		delta[i] = U2_[i] - U1_[i];
	
	dv[2] = delta[2] - v * delta[0];
	dv[1] = (GAMMA - 1.0) / (a * a) * (delta[0] * (h - u * u) + u * delta[1] - delta[3] + v * (delta[2] - v * delta[0]));
	dv[0] = ((u + a) * delta[0] - delta[1] - a * dv[1]) / (2.0 * a);
	dv[3] = delta[0] - dv[0] - dv[1];

	if (u - a > 0)
	{
		F_[0] = U1_[1];
		F_[1] = U1_[1] * U1_[1] / U1_[0] + p1;
		F_[2] = U1_[1] * U1_[2] / U1_[0];
		F_[3] = U1_[1] * h1;
	}
	else if (u > 0)
	{
		F_[0] = U1_[1] + (u - a) * dv[0];
		F_[1] = U1_[1] * U1_[1] / U1_[0] + p1 + (u - a) * (u - a) * dv[0];
		F_[2] = U1_[1] * U1_[2] / U1_[0] + (u - a) * v * dv[0];
		F_[3] = U1_[1] * h1 + (u - a) * (h - a * u) * dv[0];
	}
	else if (u + a > 0)
	{
		F_[0] = U2_[1] - (u + a) * dv[3];
		F_[1] = U2_[1] * U2_[1] / U2_[0] + p2 - (u + a) * (u + a) * dv[3];
		F_[2] = U2_[1] * U2_[2] / U2_[0] - (u + a) * v * dv[3];
		F_[3] = U2_[1] * h2 - (u + a) * (h + a * u) * dv[3];
	}
	else
	{
		F_[0] = U2_[1];
		F_[1] = U2_[1] * U2_[1] / U2_[0] + p2;
		F_[2] = U2_[1] * U2_[2] / U2_[0];
		F_[3] = U2_[1] * h2;
	}
}

