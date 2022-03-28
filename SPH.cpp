#include "Instances.h"
//#define crit1 20
Type W1(Type u, Type Sh)
{
	Type temp = 0.23873 / (Sh * Sh * Sh * Sh);
	if ((u >= 0) && (u <= 1))
	{
		return temp * (3 * u * u - 4 * u);
	}
	if ((u >= 1) && (u <= 2))
	{
		return -temp * (2 - u) * (2 - u);
	}
	return 0;
}
Type W2(Type u, Type Sh, Type r)
{
	Type temp = 0.23874 / (Sh * Sh * Sh * Sh);
	if ((u >= 0) && (u <= 1))
	{
		return temp * (3 * u - 4) / Sh;
	}
	if ((u >= 1) && (u <= 2))
	{
		return temp * (2 - u) * (2 - u) / r;
	}
	return 0.0;
}
Type sign(Type num)
{
	if (num > 0) return 1;
	return -1;
}
Type BigP(Type temp, Type Rab, Type r, Type Sh)
{
	if (temp <= 0)
	{
		Type mu = Sh * temp / (r * r + NUsph * NUsph * Sh * Sh);
		return (-Asph * CS * mu + Bsph * mu * mu) / Rab;
	}
	return 0;
}
Type BigM(Type u)
{
	if ((u >= 0) && (u <= 1))
	{
		return (1.33 - 1.2 * u * u + 0.5 * pow(u, 3)) * pow(u, 3);
	}
	if ((u >= 1) && (u < 2))
	{
		return 0.066 - 2.67 * pow(u, 3) + 3 * pow(u, 4) - 1.2 * pow(u, 5) + 0.16 * pow(u, 6);
	}
	return 1;
}

void VelocityFunction_S(model& res, model vec, int index)
{
	model t_model;
	Type oVx, oVy, oVz;

	Type rx, ry, rz, rVx, rVy, rVz, r, H, W, Fx, 
	Fy, Fz, Fgx, Fgy, Fgz, u, temp, Rab, PP, tDe, tEn;
	model tmp;

	tmp = vec;

	Fx = Fy = Fz = 0;
	Fgx = Fgy = Fgz = 0;
	tDe = tEn = t_model.De = t_model.En = t_model.Pr = 0;

	Fgrav0 = 5e16;
	if (tmp.x < -0.45) Fgrav0 = 5e14;
	//if ((tmp.x > -0.45) && ( tmp.x < -0.3)) Fgrav0 = 1e16;
	if ((tmp.x > -0.3) && (tmp.x < 0)) Fgrav0 = 5e16;


	omp_set_num_threads(4);
#pragma omp parallel for private(rx, ry, rz, rVx, rVy, rVz, r, Sh, u, W, H, temp, Rab, PP) \
reduction(+:Fx) reduction(+:Fy) reduction(+:Fz) reduction(+:Fgx) reduction(+:Fgy) reduction(+:Fgz) \
reduction(+:tDe) reduction(+:tEn)
	for (int i = 0; i < SPHsize; i++)
	{
		if (i != index)
		{
			rx = tmp.x - M[i].x; rVx = tmp.Vx - M[i].Vx;
			ry = tmp.y - M[i].y; rVy = tmp.Vy - M[i].Vy;
			rz = tmp.z - M[i].z; rVz = tmp.Vz - M[i].Vz;
			r = sqrt(rx * rx + ry * ry + rz * rz);

			Sh = 2 * pow(m / tmp.De, 0.33);
			u = r / Sh;

			H = tmp.Pr / (tmp.De * tmp.De) + M[i].Pr / (M[i].De * M[i].De);
			temp = rVx * rx + rVy * ry + rVz * rz;
			Rab = (tmp.De + M[i].De) / 2;
			PP = BigP(temp, Rab, r, Sh);
			H += PP;                 

			W = W2(u, Sh, r);

			/// F (begin)
			Fx -= m * H * rx * W;
			Fy -= m * H * ry * W;
			Fz -= m * H * rz * W;
			/// F (end)

			// Hydro (begin)
			tDe += m * temp * W;
			tEn += -0.5 * m * H * temp * W;
			// Hydro (end)

			//
			temp = BigM(u);
			if (r != 0)
			{
				Fgx += -Fgrav0 * temp * m * rx / pow(r, 3);
				Fgy += -Fgrav0 * temp * m * ry / pow(r, 3);
				Fgz += -Fgrav0 * temp * m * rz / pow(r, 3);
			}
			//
		}
	}

	t_model.De = tDe;
	t_model.En = tEn;
	t_model.Pr = 0.66 * t_model.De * t_model.En;
	
	Type p1, p2;
	p1 = pow((tmp.x + mu) * (tmp.x + mu) + tmp.y * tmp.y + tmp.z * tmp.z, 1.5);
	p2 = pow((tmp.x - (1 - mu)) * (tmp.x - (1 - mu)) + tmp.y * tmp.y + tmp.z * tmp.z, 1.5);

	t_model.x = tmp.Vx;
	t_model.Vx = 2 * tmp.Vy + tmp.x - (1 - mu) * (tmp.x + mu) / p1 - mu * (tmp.x - (1 - mu)) / p2;
	t_model.y = tmp.Vy;
	t_model.Vy = -2 * tmp.Vx + tmp.y - (1 - mu) * tmp.y / p1 - mu * tmp.y / p2;
	t_model.z = tmp.Vz;
	t_model.Vz = -4 * PI * PI * (1 - mu) * tmp.z / p1 - 4 * PI * PI * mu * tmp.z / p2;

	t_model.Vx += Fx;
	t_model.Vy += Fy;
	t_model.Vz += Fz;

	t_model.Vx += Fgx;
	t_model.Vy += Fgy;
	t_model.Vz += Fgz;

	res = t_model;
}
