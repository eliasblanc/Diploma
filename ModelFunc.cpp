#include "Instances.h"

void ModelFunc(model& t_model, model vec, Type De, Type En)
{
	t_model.De = De;
	t_model.En = En;
	t_model.Pr = 0.666 * t_model.De * t_model.En;

	Type p1, p2;
	p1 = pow((vec.x + mu) * (vec.x + mu) + vec.y * vec.y + vec.z * vec.z, 1.5);
	p2 = pow((vec.x - (1 - mu)) * (vec.x - (1 - mu)) + vec.y * vec.y + vec.z * vec.z, 1.5);

	t_model.x  = vec.Vx;
	t_model.Vx =  2 * vec.Vy + vec.x - (1 - mu) * (vec.x + mu) / p1 - mu * (vec.x - (1 - mu)) / p2;
	t_model.y  = vec.Vy;
	t_model.Vy = -2 * vec.Vx + vec.y - (1 - mu) * vec.y / p1 - mu * vec.y / p2;
	t_model.z  = vec.Vz;
	t_model.Vz = -4 * PI * PI * (1 - mu) * vec.z / p1 - 4 * PI * PI * mu * vec.z / p2;
}