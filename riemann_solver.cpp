#include <cmath>
#include "interface.h"
#include "cell.h"

void Interface::roe(char direction)
// Roe's Riemann solver for 4-element vectors
// Left and right values of the Riemann problem are defined by U1_[4] and U2_[4]
// The calculated numerical flux is stored in F_[4]
// Reference: Toro, section 11.2
{
	double h1, h2, h, p1, p2, u, v, a, theta, delta[4], dv[4], temp;
	int i;
	
	if (direction == 'y')
	{
			temp = U1_[1];
			U1_[1] = U1_[2];
			U1_[2] = temp;
			
			temp = U2_[1];
			U2_[1] = U2_[2];
			U2_[2] = temp;			
	}
	
	// calculate pressure and enthalpy in left and right region
	p1 = (GAMMA - 1.0) * (U1_[3] - 0.5 * (U1_[1] * U1_[1] + U1_[2] * U1_[2]) / U1_[0]);
	p2 = (GAMMA - 1.0) * (U2_[3] - 0.5 * (U2_[1] * U2_[1] + U2_[2] * U2_[2]) / U2_[0]);
	h1 = (U1_[3] + p1) / U1_[0];
	h2 = (U2_[3] + p2) / U2_[0];
	
	// calculate the Roe average value of u, v, h and a
	theta = sqrt(U1_[0]) / (sqrt(U1_[0]) + sqrt(U2_[0]));
	u = theta * U1_[1] / U1_[0] + (1.0 - theta) * U2_[1] / U2_[0];
	v = theta * U1_[2] / U1_[0] + (1.0 - theta) * U2_[2] / U2_[0];
	h = theta * h1 + (1.0 - theta) * h2;
	a = sqrt((GAMMA - 1.0) * (h - 0.5 * u * u - 0.5 * v * v));
	
	// calculate the jump at the interface
	for (i = 0; i < 4; i++)
		delta[i] = U2_[i] - U1_[i];
	
	// decompose the jump at the interface into the four characteristics
	dv[2] = delta[2] - v * delta[0];
	dv[1] = (GAMMA - 1.0) / (a * a) * (delta[0] * (h - u * u) + u * delta[1] - delta[3] + v * (delta[2] - v * delta[0]));
	dv[0] = ((u + a) * delta[0] - delta[1] - a * dv[1]) / (2.0 * a);
	dv[3] = delta[0] - dv[0] - dv[1];
	
	// calculate the numerical flux
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
	
	if (direction == 'y')
	{
		temp = F_[1];
		F_[1] = F_[2];
		F_[2] = temp;
		
		temp = U1_[1];
		U1_[1] = U1_[2];
		U1_[2] = temp;
		
		temp = U2_[1];
		U2_[1] = U2_[2];
		U2_[2] = temp;
	}
}

void Interface::hlle(char direction)
// HLLE Riemann solver for 4-element vectors
// Left and right values of the Riemann problem are defined by U1_[4] and U2_[4]
// The calculated numerical flux is stored in F_[4]
// Reference: Toro, section 10.3
{
	double h1, h2, h, p1, p2, u1, u2, u, v1, v2, v, a1, a2, a, theta, s1, s2, F1[4], F2[4], temp;
	int i;

	if (direction == 'y')
	{
			temp = U1_[1];
			U1_[1] = U1_[2];
			U1_[2] = temp;
			
			temp = U2_[1];
			U2_[1] = U2_[2];
			U2_[2] = temp;			
	}
	
	// calculate the pressure, velocity and enthalpy in the left and right region	
	p1 = (GAMMA - 1.0) * (U1_[3] - 0.5 * (U1_[1] * U1_[1] + U1_[2] * U1_[2]) / U1_[0]);
	p2 = (GAMMA - 1.0) * (U2_[3] - 0.5 * (U2_[1] * U2_[1] + U2_[2] * U2_[2]) / U2_[0]);
	h1 = (U1_[3] + p1) / U1_[0];
	h2 = (U2_[3] + p2) / U2_[0];
	u1 = U1_[1] / U1_[0];
	u2 = U2_[1] / U2_[0];
	v1 = U1_[2] / U1_[0];
	v2 = U2_[2] / U2_[0];
	a1 = sqrt(GAMMA * p1 / U1_[0]);
	a2 = sqrt(GAMMA * p2 / U2_[0]);

	// calculate the numerical flux in the left and right region
	F1[0] = U1_[1];
	F1[1] = U1_[1] * U1_[1] / U1_[0] + p1;
	F1[2] = U1_[1] * U1_[2] / U1_[0];
	F1[3] = U1_[1] * h1;
	F2[0] = U2_[1];
	F2[1] = U2_[1] * U2_[1] / U2_[0] + p2;
	F2[2] = U2_[1] * U2_[2] / U2_[0];
	F2[3] = U2_[1] * h2;

	// calculate the Roe average value of u, v, a for wave speed estimation
	theta = sqrt(U1_[0]) / (sqrt(U1_[0]) + sqrt(U2_[0]));
	u = theta * u1 + (1.0 - theta) * u2;
	v = theta * v1 + (1.0 - theta) * v2;	
	h = theta * h1 + (1.0 - theta) * h2;
	a = sqrt((GAMMA - 1.0) * (h - 0.5 * u * u - 0.5 * v * v));

	// estimate the minimum and maximum wave speed	
	s1 = fmin(u - a, u1 - a1);
	s2 = fmax(u + a, u2 + a2);

	// calculate the numerical flux
	if (s1 > 0)	// use numerical flux in the left region
		for (i = 0; i < 4; i++)
			F_[i] = F1[i];
	else if (s2 > 0)	// use numerical flux in the star region
		for (i = 0; i < 4; i++)
			F_[i] = (s2 * F1[i] - s1 * F2[i] + s1 * s2 * (U2_[i] - U1_[i])) / (s2 - s1);
	else	// use numerical flux in the right region
		for (i = 0; i < 4; i++)
			F_[i] = F2[i];
			
	if (direction == 'y')
	{
		temp = F_[1];
		F_[1] = F_[2];
		F_[2] = temp;
		
		temp = U1_[1];
		U1_[1] = U1_[2];
		U1_[2] = temp;
		
		temp = U2_[1];
		U2_[1] = U2_[2];
		U2_[2] = temp;
	}
}

void Interface::hllc(char direction)
// HLLC Riemann solver for 4-element vectors
// Left and right values of the Riemann problem are defined by U1_[4] and U2_[4]
// The calculated numerical flux is stored in F_[4]
// Reference: Toro, section 10.6
{
	double h1, h2, p1, p2, pstar, u1, u2, a1, a2, s1, s2, sstar, Ustar[4], F1[4], F2[4], temp;
	int i;
	
	if (direction == 'y')
	{
			temp = U1_[1];
			U1_[1] = U1_[2];
			U1_[2] = temp;
			
			temp = U2_[1];
			U2_[1] = U2_[2];
			U2_[2] = temp;			
	}
	
	// calculate the pressure, velocity and enthalpy in the left and right region	
	p1 = (GAMMA - 1.0) * (U1_[3] - 0.5 * (U1_[1] * U1_[1] + U1_[2] * U1_[2]) / U1_[0]);
	p2 = (GAMMA - 1.0) * (U2_[3] - 0.5 * (U2_[1] * U2_[1] + U2_[2] * U2_[2]) / U2_[0]);
	h1 = (U1_[3] + p1) / U1_[0];
	h2 = (U2_[3] + p2) / U2_[0];
	u1 = U1_[1] / U1_[0];
	u2 = U2_[1] / U2_[0];
	a1 = sqrt(GAMMA * p1 / U1_[0]);
	a2 = sqrt(GAMMA * p2 / U2_[0]);

	// calculate the numerical flux in the left and right region
	F1[0] = U1_[1];
	F1[1] = U1_[1] * U1_[1] / U1_[0] + p1;
	F1[2] = U1_[1] * U1_[2] / U1_[0];
	F1[3] = U1_[1] * h1;
	F2[0] = U2_[1];
	F2[1] = U2_[1] * U2_[1] / U2_[0] + p2;
	F2[2] = U2_[1] * U2_[2] / U2_[0];
	F2[3] = U2_[1] * h2;

	// pressure-based wave speed estimation	
	pstar = fmax(0.0, 0.5 * (p1 + p2) - 0.5 * (u2 - u1) * 0.5 * (U1_[0] + U2_[0]) * 0.5 * (a1 + a2));
	if (pstar < p1)
		s1 = u1 - a1;
	else
		s1 = u1 - a1 * sqrt(1.0 + (GAMMA + 1.0) / (2.0 * GAMMA) * (pstar / p1 - 1.0));
	if (pstar < p2)
		s2 = u2 + a2;
	else
		s2 = u2 + a2 * sqrt(1.0 + (GAMMA + 1.0) / (2.0 * GAMMA) * (pstar / p2 - 1.0));
	// calculate the estimated wave speed in the star region
	sstar = (p2 - p1 + U1_[1] * (s1 - u1) - U2_[1] * (s2 - u2)) / (U1_[0] * (s1 - u1) - U2_[0] * (s2 - u2));
	
	// calculate the numerical flux
	if (s1 > 0)	// use numerical flux in the left region
		for (i = 0; i < 4; i++)
			F_[i] = F1[i];
	else if (sstar > 0)	// use numerical flux in the left star region
	{
		Ustar[0] = U1_[0] * (s1 - u1) / (s1 - sstar);
		Ustar[1] = Ustar[0] * sstar;
		Ustar[2] = Ustar[0] * U1_[2] / U1_[0];
		Ustar[3] = Ustar[0] * (U1_[3] / U1_[0] + (sstar - u1) * (sstar + p1 / (U1_[0] * (s1 - u1))));
		for (i = 0; i < 4; i++)
			F_[i] = F1[i] + s1 * (Ustar[i] - U1_[i]);
	}
	else if (s2 > 0)	// use numerical flux in the right star region
	{
		Ustar[0] = U2_[0] * (s2 - u2) / (s2 - sstar);
		Ustar[1] = Ustar[0] * sstar;
		Ustar[2] = Ustar[0] * U2_[2] / U2_[0];
		Ustar[3] = Ustar[0] * (U2_[3] / U2_[0] + (sstar - u2) * (sstar + p2 / (U2_[0] * (s2 - u2))));
		for (i = 0; i < 4; i++)
			F_[i] = F2[i] + s2 * (Ustar[i] - U2_[i]);
	}
	else	// use numerical flux in the right region
		for (i = 0; i < 4; i++)
			F_[i] = F2[i];
			
	if (direction == 'y')
	{
		temp = F_[1];
		F_[1] = F_[2];
		F_[2] = temp;
		
		temp = U1_[1];
		U1_[1] = U1_[2];
		U1_[2] = temp;
		
		temp = U2_[1];
		U2_[1] = U2_[2];
		U2_[2] = temp;
	}
}
