#include <iostream>
#include <cmath>
#include <fstream>
#include "model.h"

void grid(double radius, int division, const char *filename)
// generate a file contain the grid information
{
	int Nx, Nx1, Nx2, Ny, Ny1, Ny2;
	int interfaceid[4], cellid[2];
	double dx1, dx2, dy1, dy2, startx, starty, x, y, dx, dy, x0, y0;
	int i, j, k;
	ofstream fout;
	
	Ny2 = division * 2;
	dy2 = 2 * radius / Ny2;
	dy1 = dy2 * 2.0;
	dx1 = dy1;
	dx2 = dy2;
	Ny1 = Ny2 / 2;
	Nx1 = Ny1 * 3 / 4;
	Nx2 = Ny2;
	Nx = Nx1 + Nx2;
	Ny = Ny1 + Ny2;

// calculate all the coordinate information that is needed	
	int saw_position[Ny2], v_num, h_num[Ny2], c_num[Ny2 + 1];
	
	startx = - Nx2 * dx2 - Nx1 * dx1 + 0.5 * dx1;
	starty = - 0.5 * Ny2 * dy2 - 0.5 * Ny1 * dy1 + 0.5 * dy1;
	
	x0 = - Nx2 * dx2 + 0.5 * dx2;
	y0 = - 0.5 * Ny2 * dy2 + 0.5 * dy2;

// saw_position defines the boundary of the cynlinder	
	for (j = 0; j < Ny2; j++)
	{
		i = 0;
		while (pow(x0 + i * dx2, 2.0) + pow(y0 + j * dy2, 2.0) > pow(0.5 * Ny2 * dy2, 2.0))
		{
			i = i + 1;
		}
		if (i == Nx2)
		{
			saw_position[j] = Nx - 1;
		}
		else
		{
			saw_position[j] = Nx1 + i;
		}
	}
	
	v_num = (Nx + 1) * Ny1;
	for (i = 0; i < Ny2; i++)
	{
		v_num = v_num + saw_position[i] + 1;
	}
	
	for (i = 0; i < Ny2; i++)
	{
		h_num[i] = 0;
		for (j = 0; j < i; j++)
		{
			if (saw_position[j] > saw_position[j + 1])
			{
				h_num[i] = h_num[i] + saw_position[j];
			}
			else
			{
				h_num[i] = h_num[i] + saw_position[j + 1];
			}
		}
	}
	
	for (i = 0; i < Ny2 + 1; i++)
	{
		c_num[i] = 0;
		for (j = 0; j < i; j++)
		{
			c_num[i] = c_num[i] + saw_position[j];
		}
	}

// start writing	
	fout.open(filename);
	
	fout << "# Nvertical, Nx, Ny" << endl;
	fout << v_num << '\t' << Nx << '\t' << Ny << endl;
	
// write cell information
	fout << "# cells" << endl;

// below the cylinder	
	for (j = 0; j < Ny1 / 2; j++)
	{
		for (i = 0; i < Nx; i++)
		{
			interfaceid[0] = j * (Nx + 1) + i;
			interfaceid[1] = interfaceid[0] + 1;
			interfaceid[2] = v_num + j * Nx + i;
			interfaceid[3] = interfaceid[2] + Nx;
			cellid[0] = j * Nx + i;
			y = starty + j * dy1;
			dy = dy1;
			if (i < Nx1)
			{
				x = startx + i * dx1;
				dx = dx1;
			}
			else
			{
				x = startx + (Nx1 - 1) * dx1 + (dx1 + dx2) / 2.0 + (i - Nx1) * dx2;
				dx = dx2;
			}
			for (k = 0; k < 4; k++)
			{
				fout << interfaceid[k] << '\t';
			}
			fout << x << '\t' << y << '\t' << dx << '\t' << dy << '\t' << cellid[0] << endl;
		}
	}
	
// first line of the cylinder (j = Ny1 / 2)
	for (i = 0; i < saw_position[0]; i++)
	{
		interfaceid[0] = Ny1 / 2 * (Nx + 1) + i;
		interfaceid[1] = interfaceid[0] + 1;
		interfaceid[2] = v_num + Ny1 / 2 * Nx + i;
		interfaceid[3] = interfaceid[2] + Nx;
		cellid[0] = Ny1 / 2 * Nx + i;
		y = y0;
		dy = dy2;
		if (i < Nx1)
		{
			x = startx + i * dx1;
			dx = dx1;
		}
		else
		{
			x = x0 + (i - Nx1) * dx2;
			dx = dx2;
		}
		for (k = 0; k < 4; k++)
		{
			fout << interfaceid[k] << '\t';
		}
		fout << x << '\t' << y << '\t' << dx << '\t' << dy << '\t' << cellid[0] << endl;
	}

// cross the cylinder	
	for (j = Ny1 / 2 + 1; j < Ny2 + Ny1 /2; j++)
	{
		for (i = 0; i < saw_position[j - Ny1 / 2]; i++)
		{
			interfaceid[0] = Ny1 / 2 * (Nx + 1) + c_num[j - Ny1 / 2] + (j - Ny1 / 2) + i;
			interfaceid[1] = interfaceid[0] + 1;
			if (saw_position[j - Ny1 / 2 - 1] > saw_position[j - Ny1 / 2])
			{
				interfaceid[2] = v_num + Ny1 / 2 * Nx + Nx + h_num[j - Ny1 / 2] - saw_position[j - Ny1 / 2 - 1] + i;
			}
			else
			{
				interfaceid[2] = v_num + Ny1 / 2 * Nx + Nx + h_num[j - Ny1 / 2] - saw_position[j - Ny1 / 2] + i;
			}
			interfaceid[3] = v_num + Ny1 / 2 * Nx + Nx + h_num[j - Ny1 / 2] + i;
			cellid[0] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2] + i;
			y = y0 + (j - Ny1 / 2) * dy2;
			dy = dy2;
			if (i < Nx1)
			{
				x = startx + i * dx1;
				dx = dx1;
			}
			else
			{
				x = x0 + (i - Nx1) * dx2;
				dx = dx2;
			}
			for (k = 0; k < 4; k++)
			{
				fout << interfaceid[k] << '\t';
			}
			fout << x << '\t' << y << '\t' << dx << '\t' << dy << '\t' << cellid[0] << endl;
		}
	}

// above the cylinder	
	for (j = Ny2 + Ny1 / 2; j < Ny; j++)
	{
		for (i = 0; i < Nx; i++)
		{
			interfaceid[0] = Ny1 / 2 * (Nx + 1) + c_num[Ny2] + Ny2 + (j - Ny1 / 2 - Ny2) * (Nx + 1) + i;
			interfaceid[1] = interfaceid[0] + 1;
			interfaceid[2] = v_num + (Ny1 / 2 + 1) * Nx + h_num[Ny2 - 1] + (j - Ny1 / 2 - Ny2) * Nx + i;
			interfaceid[3] = interfaceid[2] + Nx;
			cellid[0] = Ny1 / 2 * Nx + c_num[Ny2] + (j - Ny1 / 2 - Ny2) * Nx + i;
			y = starty + (Ny1 / 2 - 1) * dy1  + (Ny2 - 1) * dy2 + (dy1 + dy2) + (j - Ny1 / 2 - Ny2) * dy1;
			dy = dy1;
			if (i < Nx1)
			{
				x = startx + i * dx1;
				dx = dx1;
			}
			else
			{
				x = x0 + (i - Nx1) * dx2;
				dx = dx2;
			}
			for (k = 0; k < 4; k++)
			{
				fout << interfaceid[k] << '\t';
			}
			fout << x << '\t' << y << '\t' << dx << '\t' << dy << '\t' << cellid[0] << endl;
		}
	}

// write interface information	
	fout << "# interfaces" << endl;

// vertical interfaces
// below the cylinder	
	for (j = 0; j < Ny1 /2; j++)
	{
		for (i = 0; i < Nx + 1; i++)
		{
			if (i == 0)
			{
				cellid[0] = -1;
				cellid[1] = j * Nx;
			}
			else if (i == Nx)
			{
				cellid[0] = j * Nx + Nx - 1;
				cellid[1] = -2;
			}
			else
			{
				cellid[0] = j * Nx + i - 1;
				cellid[1] = cellid[0] + 1;
			}
			interfaceid[0] = j * (Nx + 1) + i;
			fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
		}
	}

// across the cylinder	
	for (j = Ny1 / 2; j < Ny1 / 2 + Ny2; j++)
	{
		for (i = 0; i <= saw_position[j - Ny1 / 2]; i++)
		{
			if (i == 0)
			{
				cellid[0] = -1;
				cellid[1] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2 + 1] - saw_position[j - Ny1 / 2];
			}
			else if (i == saw_position[j - Ny1 / 2])
			{
				cellid[0] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2 + 1] - 1;
				cellid[1] = -3;
			}
			else
			{
				cellid[0] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2 + 1] - saw_position[j - Ny1 / 2] + i - 1;
				cellid[1] = cellid[0] + 1;
			}
			interfaceid[0] = Ny1 / 2 * (Nx + 1) + c_num[j - Ny1 / 2 + 1] - saw_position[j - Ny1 / 2] + j - Ny1 / 2 + i;
			fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
		}
	}

// above the cylinder	
	for (j = Ny1 / 2 + Ny2; j < Ny; j++)
	{
		for (i = 0; i < Nx + 1; i++)
		{
			if (i == 0)
			{
				cellid[0] = -1;
				cellid[1] = Ny1 / 2 * Nx + c_num[Ny2] + (j - Ny1 / 2 - Ny2) * Nx;
			}
			else if (i == Nx)
			{
				cellid[0] = Ny1 / 2 * Nx + c_num[Ny2] + (j - Ny1 / 2 - Ny2) * Nx + Nx - 1;
				cellid[1] = -2;
			}
			else
			{
				cellid[0] = Ny1 / 2 * Nx + c_num[Ny2] + (j - Ny1 / 2 - Ny2) * Nx + i - 1;
				cellid[1] = cellid[0] + 1;
			}
			interfaceid[0] = Ny1 / 2 * (Nx + 1) + c_num[Ny2] + Ny2 + (j - Ny1 / 2 - Ny2) * (Nx + 1) + i;
			fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
		}
	}

// horizontal intefaces
// below the cylinder	
	for (j = 0; j < Ny1 / 2; j++)
	{
		for (i = 0; i < Nx; i++)
		{
			if (j == 0)
			{
				cellid[0] = -2;
				cellid[1] = i;
			}
			else
			{
				cellid[0] = (j - 1) * Nx + i;
				cellid[1] = cellid[0] + Nx;
			}
			interfaceid[0] = v_num + j * Nx + i;
			fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
		}
	}
	
// first line of the cylinder (j = Ny1 / 2)
	for (i = 0; i < Nx; i++)
	{
		cellid[0] = (Ny1 / 2 - 1) * Nx + i;
		if (i < saw_position[0])
		{
			cellid[1] = cellid[0] + Nx;
		}
		else
		{
			cellid[1] = -3;
		}
		interfaceid[0] = v_num + Ny1 / 2 * Nx + i;
		fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
	}

// across the cylinder	
	for (j = Ny1 / 2 + 1; j < Ny1 /2 + Ny2; j++)
	{
		for (i = 0; i <= saw_position[j - Ny1 / 2] - 1 || i <= saw_position[j - Ny1 / 2 - 1] - 1; i++)
		{
			if (saw_position[j - Ny1 / 2] > saw_position[j - Ny1 / 2 - 1])
			{
				if (i < saw_position[j - Ny1 / 2 - 1])
				{
					cellid[0] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2] - saw_position[j - Ny1 / 2 - 1] + i;
				}
				else
				{
					cellid[0] = -3;
				}
				cellid[1] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2] + i;
				interfaceid[0] = v_num + (Ny1 / 2 + 1) * Nx + h_num[j - Ny1 / 2] - saw_position[j - Ny1 / 2] + i;
			}
			else
			{
				cellid[0] = Ny1 /2 * Nx + c_num[j - Ny1 / 2] - saw_position[j - Ny1 / 2 - 1] + i;
				if (i < saw_position[j - Ny1 / 2])
				{
					cellid[1] = Ny1 / 2 * Nx + c_num[j - Ny1 / 2] + i;
				}
				else
				{
					cellid[1] = -3;
				}
				interfaceid[0] = v_num + (Ny1 / 2 + 1) * Nx + h_num[j - Ny1 / 2] - saw_position[j - Ny1 / 2 - 1] + i;
			}
			fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
		}
	}
	
// last line of the cylinder (j = Ny1 / 2 + Ny2)
	for (i = 0; i < Nx; i++)
	{
		if (i < saw_position[Ny2 - 1])
		{
			cellid[0] = Ny1 / 2 * Nx + c_num[Ny2] - saw_position[Ny2 - 1] + i;
		}
		else
		{
			cellid[0] = -3;
		}
		cellid[1] = Ny1 / 2 * Nx + c_num[Ny2] + i;
		interfaceid[0] = v_num + (Ny1 / 2 + 1) * Nx + h_num[Ny2 - 1] + i;
		fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
	}

// above the cylinder	
	for (j = Ny1 / 2 + Ny2 + 1; j < Ny + 1; j++)
	{
		for (i = 0; i < Nx; i++)
		{
			if (j == Ny)
			{
				cellid[0] = Ny1 / 2 * Nx + c_num[Ny2] + (Ny1 / 2 - 1) * Nx + i;
				cellid[1] = -2;
			}
			else
			{
				cellid[0] = Ny1 / 2 * Nx + c_num[Ny2] + (j - Ny1 / 2 - Ny2 - 1) * Nx + i;
				cellid[1] = cellid[0] + Nx;
			}
			interfaceid[0] = v_num + (Ny1 / 2 + 1) * Nx + h_num[Ny2 - 1] + (j - Ny1 / 2 - Ny2) * Nx + i;
			fout << cellid[0] << '\t' << cellid[1] << '\t' << interfaceid[0] << endl;
		}
	}
	
	fout.close();
}
