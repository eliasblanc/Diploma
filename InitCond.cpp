#include "Instances.h"

void InitConditions(model* model)
{
	Type angles[9][3] = {
		{1,      0,     0},    {0.866, 0.5,   0},
		{0.940,  0.342, 0},    {0.985, 0.174, 0},
		{0.866, -0.5,   0},    {0.940,-0.342, 0},
		{0.985, -0.174, 0},    {0.707, 0,     0.707},
		{0.707,  0,    -0.707}
	};
	size_t j(0);
	for (size_t i = 0; i < Nparts; i++)
	{
		model[i].x  = L;
		model[i].y  = 0;
		model[i].z  = 0;
		model[i].Vx = angles[j][0] * tVelo;
		model[i].Vy = angles[j][1] * tVelo;
		model[i].Vz = angles[j][2] * tVelo;
		model[i].De = dencity;
		model[i].En = energy;
		model[i].Pr = pressure;
		j == 8 ? j = 0 : j++;
	}
}