#include "Instances.h"

void MatrPlusVecParall(model* A, model vec, int index)
{
	A[index].x  += vec.x;
	A[index].y  += vec.y;
	A[index].z  += vec.z;
	A[index].Vx += vec.Vx;
	A[index].Vy += vec.Vy;
	A[index].Vz += vec.Vz;
	A[index].De += vec.De;
	A[index].En += vec.En;
	A[index].Pr += vec.Pr;
}

Type LastSumP(Type v1, Type v2, Type v3, Type v4, Type h)
{
	return h * (v1 * 0.166 + v2 * 0.333 + v3 * 0.333 + v4 * 0.166);
}

void LastSum(model& tk, model k1, model k2, model k3, model k4, Type h)
{
	tk.x  = LastSumP(k1.x, k2.x, k3.x, k4.x, h);
	tk.y  = LastSumP(k1.y, k2.y, k3.y, k4.y, h);
	tk.z  = LastSumP(k1.z, k2.z, k3.z, k4.z, h);
	tk.Vx = LastSumP(k1.Vx, k2.Vx, k3.Vx, k4.Vx, h);
	tk.Vy = LastSumP(k1.Vy, k2.Vy, k3.Vy, k4.Vy, h);
	tk.Vz = LastSumP(k1.Vz, k2.Vz, k3.Vz, k4.Vz, h);
	tk.De = LastSumP(k1.De, k2.De, k3.De, k4.De, h);
	tk.En = LastSumP(k1.En, k2.En, k3.En, k4.En, h);
	tk.Pr = LastSumP(k1.Pr, k2.Pr, k3.Pr, k4.Pr, h);
}

void ReadFile()
{
	std::ifstream rcf;
	rcf.open("rcf_H.dat");

	std::string str_temp;

	double temp;

	for (int i = 0; i < 4; i++)
		std::getline(rcf, str_temp);

	int D_index = 0, T_index = 0;
	for (int i = 0; i < 32829; i++)
	{
		rcf >> temp;
		if (i % 353 == 0)
		{
			
			lnD[D_index] = exp(temp) / 28223.0;
			++D_index;
		}
		
		rcf >> temp;
		if (T_index < 353)
		{
			//lnT[T_index] = exp(temp) / 3000;
			lnT[T_index] = exp(temp);
			++T_index;
		}
		
		rcf >> temp;
		lnL[i] = 10 * exp(temp) / (4.30339);
	}

	rcf.close();
}

void falled(model& model)
{
	model.x = 10; model.y = model.z = 0;
	model.Vx = model.Vy = model.Vz = tVelo * 0.577;
	model.De = dencity;
	model.En = energy;
	model.Pr = pressure;
}

void speedlim(model& model)
{
	model.x = 10; model.y = model.z = 0;
	model.Vx /= 100; model.Vy /= 100; model.Vz /= 100;
	model.De = dencity;
	model.En = energy;
	model.Pr = pressure;
}
