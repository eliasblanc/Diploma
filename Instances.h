#pragma once
#include "iostream"
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <omp.h>

typedef long double Type;

class model
{
public:
	Type x, y, z;        // Coord
	Type Vx, Vy, Vz;     // Vel
	Type De, En, Pr;     // Hydro

	model& operator += (model model)
	{
		x += model.x;   y += model.y;   z += model.z;
		Vx += model.Vx; Vy += model.Vy; Vz += model.Vz;
		De += model.De; En += model.En; Pr += model.Pr;
		return *this;
	}
	model& operator *= (Type value)
	{
		x *= value;  y *= value;  z *= value;
		Vx *= value; Vy *= value; Vz *= value;
		De *= value; En *= value; Pr *= value;
		return *this;
	}
};


extern std::ofstream fout;
extern std::ofstream cool;
extern std::ofstream notcool;
extern std::ofstream keenetic;
extern model* M;
extern Type* lnD;
extern Type* lnT;
extern Type* lnL;

extern Type mu;
extern Type L;
extern const int dim;
extern Type h;
extern Type PI;
extern Type Sh;
extern Type KSTL;
extern Type KSTL2;

extern int falled;
extern int spdlim;
extern int area;

extern int SPHsize;
extern const int Nparts;
extern Type tVelo;

extern const Type Asph;
extern const Type Bsph;
extern const Type NUsph;
extern const Type CS;

extern Type m;
extern Type dencity;
extern Type pressure;
extern Type energy;

extern Type G;
extern Type Fgrav0;

void InitConditions_S(model* model);
void MatrPlusVecParall_S(model* A, model vec, int index);
void LastSum_S(model& tk, model k1, model k2, model k3, model k4, Type h);
void ReadFile();

Type sign(Type num);
Type W1(Type u, Type Sh);
Type W2(Type u, Type Sh, Type r);
Type BigP(Type temp, Type Rab, Type r, Type Sh);
Type BigM(Type u);
void VelocityFunction_S(model& res, model vec, int index);
void RK4_S();
