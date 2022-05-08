#include "Instances.h"

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
Type sign(Type num)
{
	if (num > 0) return 1;
	return -1;
}

void SPH(model& res, model vec, int index)
{
	Type Fx,  Fy,  Fz, 
		 Fgx, Fgy, Fgz, 
		 tDe, tEn;

	Fx  = Fy  = Fz  = \
	Fgx = Fgy = Fgz = \
	tDe = tEn = 0;

// 5e17 - right 
// 5e16 - too much
// 5e15 - good
// 3e15 - +/-
// 7e15 - hmm
// 9e15 - 5k try
	Fgrav0 = 9e15;// 5e16 5e17
	if (vec.x < -0.45) Fgrav0 = 5e14;
	if ((vec.x > -0.3) && (vec.x < 0)) Fgrav0 = 5e16;

#pragma omp parallel for reduction(+:Fx)  reduction(+:Fy)  reduction(+:Fz)  \
	                     reduction(+:Fgx) reduction(+:Fgy) reduction(+:Fgz) \
                         reduction(+:tDe) reduction(+:tEn)
	for (int i = 0; i < SPHsize; i++)
	{
		if (i != index)
		{
			Type rx, ry, rz, r,
				 Vx, Vy, Vz,
				 Sh,
				 u,
				 W,
				 H,
				 temp,
				 Rab,
				 PP;

			rx = vec.x - M[i].x; Vx = vec.Vx - M[i].Vx;
			ry = vec.y - M[i].y; Vy = vec.Vy - M[i].Vy;
			rz = vec.z - M[i].z; Vz = vec.Vz - M[i].Vz;
			r = sqrt(rx * rx + ry * ry + rz * rz);

			Sh = 2 * pow(m / vec.De, 0.33);
			u = r / Sh;

			H = vec.Pr / (vec.De * vec.De) + M[i].Pr / (M[i].De * M[i].De);
			temp = Vx * rx + Vy * ry + Vz * rz;
			Rab = (vec.De + M[i].De) / 2;
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
			tEn -= 0.5 * m * H * temp * W;
			// Hydro (end)
			
			//
			temp = BigM(u);
			if (r != 0)
			{
				Fgx -= Fgrav0 * temp * m * rx / pow(r, 3);
				Fgy -= Fgrav0 * temp * m * ry / pow(r, 3);
				Fgz -= Fgrav0 * temp * m * rz / pow(r, 3);
			}
			//
		}
	}

	////

	model t_model;
	t_model.De = t_model.En = t_model.Pr = 0;
	ModelFunc(t_model, vec, tDe, tEn);

	////



	if (fabs(Fx) > KSTL * fabs(t_model.Vx)) Fx = 0;//sign(Fx) * 0.1 * fabs(t_model.Vx);
	if (fabs(Fy) > KSTL * fabs(t_model.Vy)) Fy = 0;//sign(Fy) * 0.1 * fabs(t_model.Vy);
	if (fabs(Fz) > KSTL * fabs(t_model.Vz)) Fz = 0;//sign(Fz) * 0.1 * fabs(t_model.Vz);

	if (fabs(Fgx) > KSTL2 * fabs(t_model.Vx)) Fgx = 0; //Fgx = sign(Fgx) * 0.2 * fabs(t_model.Vx);
	if (fabs(Fgy) > KSTL2 * fabs(t_model.Vy)) Fgy = 0; //Fgy = sign(Fgy) * 0.2 * fabs(t_model.Vy);
	if (fabs(Fgz) > KSTL2 * fabs(t_model.Vz)) Fgz = 0; //Fgz = sign(Fgz) * 0.2 * fabs(t_model.Vz);

	t_model.Vx += Fx;
  t_model.Vy += Fy;
	t_model.Vz += Fz;

	t_model.Vx += Fgx;
	t_model.Vy += Fgy;
	t_model.Vz += Fgz;

	res = t_model;
}