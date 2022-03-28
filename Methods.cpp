#include "Instances.h"

void InitConditions_S(model* model)
{
	Type angles[9][3] = {
		{1, 0, 0},
        {0.866,0.5,0},
		{0.940,0.342,0},
		{0.985,0.174,0},
		{0.866,-0.5,0},
		{0.940,-0.342,0},
		{0.985,-0.174,0},
		{0.707, 0, 0.707},
		{0.707, 0, -0.707}
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

void MatrPlusVecParall_S(model* A, model vec, int index)
{
	A[index].x += vec.x;
	A[index].y += vec.y;
	A[index].z += vec.z;
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

void LastSum_S(model& tk, model k1, model k2, model k3, model k4, Type h)
{
	tk.x  = LastSumP( k1.x,  k2.x,  k3.x,  k4.x, h);
	tk.y  = LastSumP( k1.y,  k2.y,  k3.y,  k4.y, h);
	tk.z  = LastSumP( k1.z,  k2.z,  k3.z,  k4.z, h);
	tk.Vx = LastSumP(k1.Vx, k2.Vx, k3.Vx, k4.Vx, h);
	tk.Vy = LastSumP(k1.Vy, k2.Vy, k3.Vy, k4.Vy, h);
	tk.Vz = LastSumP(k1.Vz, k2.Vz, k3.Vz, k4.Vz, h);
	tk.De = LastSumP(k1.De, k2.De, k3.De, k4.De, h);
	tk.En = LastSumP(k1.En, k2.En, k3.En, k4.En, h);
	tk.Pr = LastSumP(k1.Pr, k2.Pr, k3.Pr, k4.Pr, h);
}

void ReadFile()
{
	std::ifstream file;
	file.open("rcf_H.dat");
	std::string str_temp;
	double temp;

	std::getline(file, str_temp);
	std::getline(file, str_temp);
	std::getline(file, str_temp);
	std::getline(file, str_temp);

	int D_index, T_index; D_index = T_index = 0;
	for (int i = 0; i < 32829; i++)
	{
		file >> temp;
		//printf("D = %e ", exp(temp) / 28223.0);
		//std::cout << "i = " << i << " temp = " << temp << "\n";
		if (i % 353 == 0)
		{
			//lnD[D_index] = 1e-3 * exp(temp) / 28223.0;
			lnD[D_index] = exp(temp) / 28223.0;
			++D_index;
		}
		file >> temp;
		//printf(" T = %f ", exp(temp));
		if (T_index < 353)
		{
			//lnT[T_index] = exp(temp) / 3000;
			lnT[T_index] = exp(temp);
			++T_index;
		}
		file >> temp;
		//printf(" L = %e ", 10 * exp(temp) / 4.30339);
		//printf("\n");
		//lnL[i] = 10 * exp(temp) / (4.30339 * 1e11);
		lnL[i] = 10 * exp(temp) / (4.30339);
	}

	/*for (int i = 0; i < 32829; i++)
	{
		std::cout << lnL[i] << "\n";
	}*/

	file.close();
}

void RK4_S()
{
	model* resM = new model[Nparts];

	std::memcpy(resM, M, Nparts * sizeof(model));

	omp_set_num_threads(4);
#pragma omp parallel for
	for (int i = 0; i < SPHsize; i++)
	{
		model tk, k1, k2, k3, k4, vecc, temp;
		temp = M[i];

		Type th(0.003);

		vecc = temp;
		VelocityFunction_S(k1, vecc, i);

		tk = k1;
		tk *= th / 2;
		vecc = temp;
		vecc += tk;
		VelocityFunction_S(k2, vecc, i);

		tk = k2;
		tk *= th / 2;
		vecc = temp;
		vecc += tk;
		VelocityFunction_S(k3, vecc, i);

		tk = k3;
		tk *= th;
		vecc = temp;
		vecc += tk;
		VelocityFunction_S(k4, vecc, i);

		LastSum_S(tk, k1, k2, k3, k4, h);

		double lmtr(1.01);
		if ((fabs(resM[i].En / (resM[i].En + tk.En)) > lmtr) || (fabs((resM[i].En + tk.En) / resM[i].En) > lmtr))
		{
			tk.En = 0;
		}
		if ((fabs(resM[i].Pr / (resM[i].Pr + tk.Pr)) > lmtr) || (fabs((resM[i].Pr + tk.Pr) / resM[i].Pr) > lmtr))
		{
			tk.Pr = 0;
		}
		if ((fabs(resM[i].De / (resM[i].De + tk.De)) > lmtr) || (fabs((resM[i].De + tk.De) / resM[i].De) > lmtr))
		{
			tk.De = 0;
		}

		bool flag = true;
		Type S = sqrt((temp.x - 0.0806045) * (temp.x - 0.0806045) + temp.y * temp.y + temp.z * temp.z);
		if (S < 0.07)
		{
			flag = false;
			temp.x = 10; temp.y = temp.z = 0;
			temp.Vx = temp.Vy = temp.Vz = tVelo * 0.577;
			temp.De = dencity;
			temp.En = energy;
			temp.Pr = pressure;
			resM[i] = temp;
			falled++;
		}
		if (sqrt(temp.Vx * temp.Vx + temp.Vy * temp.Vy + temp.Vz * temp.Vz) > 100 * tVelo)
		{
			flag = false;
			temp.x = 10; temp.y = temp.z = 0;
			temp.Vx /= 100;
			temp.Vy /= 100;
			temp.Vz /= 100;
			temp.De = dencity;
			temp.En = energy;
			temp.Pr = pressure;
			resM[i] = temp;
			spdlim++;
		}

		// Out of area

		//if (fabs(temp.x) > 0.7 || fabs(temp.y) > 0.7)
		//{
		//	flag = false;
		//	area++;
		//	//std::cout << "\nArea";
		//}  

		if (flag)
		{
			MatrPlusVecParall_S(resM, tk, i);

			const Type epsD = 1e-15;
			const Type epsT = 50;
			Type L(0);
			for (int j = 0; j < 93; j++)
			{
				/*printf("resM[i].De = %e", resM[i].De);
				printf("lnD[j] = %e", lnD[j]);
				printf("\n");*/
				if (fabs(lnD[j] - resM[i].De) < epsD)
				{
					for (int k = 0; k < 353; k++)
					{
						if (fabs(lnT[k] - 168e4 * resM[i].Pr / resM[i].De) < epsT)
						{
							/*printf("resM[i].Pr / ... = %f", 168e4 * resM[i].Pr / resM[i].De);
							printf(" lnT[k] = %f", lnT[k]);
							printf("\n");*/
							L = 0.1 * lnL[j * 353 + k];
							break;
						}
					}
				}
			}

			/*printf("En = %f", resM[i].En);
			printf("L = %e", L);
			printf("\n");*/

			double lmtr2(1.02);//1.05
			if ((fabs(resM[i].En / (resM[i].En - L)) > lmtr2) \
				|| (fabs((resM[i].En - L) / resM[i].En) > lmtr2) || (resM[i].En - L < 0))
			{
				//printf("YES \n");
				L = 0;
			}
			/*else
			{
				printf("NO \n");
			}*/
			resM[i].En -= L;
		}
	}

	for (int i = 0; i < SPHsize; i++)
	{
		fout << resM[i].x << " " << resM[i].y << std::endl;
	}

	int indx(7);

	//cool << resM[indx].En << std::endl;
	keenetic << (resM[indx].Vx * resM[indx].Vx + resM[indx].Vy * resM[indx].Vy + resM[indx].Vz * resM[indx].Vz) * m / 2 << std::endl;

	/*printf("D = %e\t", resM[indx].De);
	printf("En = %f\t", resM[indx].En);
	printf("Pr = %e\t", resM[indx].Pr);
	printf("T = %f\t", 168e4 * resM[indx].Pr / resM[indx].De);
	printf("(x, y, z) = (%f, %f, %f)\t", resM[indx].x, resM[indx].y, resM[indx].z);
	printf("\n");*/

	std::memcpy(M, resM, Nparts * sizeof(model));
	delete[] resM;
}