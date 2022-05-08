#include "Instances.h"

void RK4()
{
	model* resM = new model[Nparts];

	std::memcpy(resM, M, Nparts * sizeof(model));

#pragma omp parallel for
	for (int i = 0; i < SPHsize; i++)
	{
		model tk, k1, k2, k3, k4, vecc, temp;
		temp = M[i];

		Type th(0.003);

		vecc = temp;
		SPH(k1, vecc, i);

		tk = k1;
		tk *= th / 2;
		vecc = temp;
		vecc += tk;
		SPH(k2, vecc, i);

		tk = k2;
		tk *= th / 2;
		vecc = temp;
		vecc += tk;
		SPH(k3, vecc, i);

		tk = k3;
		tk *= th;
		vecc = temp;
		vecc += tk;
		SPH(k4, vecc, i);

		LastSum(tk, k1, k2, k3, k4, h);

		double lmtr(1.01);
		if ((fabs(resM[i].De / (resM[i].De + tk.De)) > lmtr) || (fabs((resM[i].De + tk.De) / resM[i].De) > lmtr)) tk.De = 0;
		if ((fabs(resM[i].En / (resM[i].En + tk.En)) > lmtr) || (fabs((resM[i].En + tk.En) / resM[i].En) > lmtr)) tk.En = 0;
		if ((fabs(resM[i].Pr / (resM[i].Pr + tk.Pr)) > lmtr) || (fabs((resM[i].Pr + tk.Pr) / resM[i].Pr) > lmtr)) tk.Pr = 0;

		bool flag = true;
		Type S = sqrt((temp.x - 0.0806045) * (temp.x - 0.0806045) + temp.y * temp.y + temp.z * temp.z);
		if (S < 0.05) // 0.07
		{
			flag = false;
			falled(temp);
			resM[i] = temp;
			Nfall++;
		}
		if (sqrt(temp.Vx * temp.Vx + temp.Vy * temp.Vy + temp.Vz * temp.Vz) > 100 * tVelo)
		{
			flag = false;
			speedlim(temp);
			resM[i] = temp;
			Nspdlim++;
		}

		if (flag)
		{
			MatrPlusVecParall(resM, tk, i);

      /// comm
      
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
			///*else
			{
				printf("NO \n");
			}
			resM[i].En -= L;
      
      /// comm
		}
	}


if ((SPHsize >= 450 * 9) && (SPHsize <= 455 * 9)) ///
	{
		for (int i = 0; i < SPHsize; i++)
		{
			double T = 168e4 * resM[i].Pr / resM[i].De;
			if (T > 3200)             color1 << "gold" << std::endl;
			if (T > 3000 && T < 3200) color1 << "orange" << std::endl;
			if (T > 2700 && T < 3000) color1 << "coral" << std::endl;
			if (T > 2300 && T < 2700) color1 << "orangered" << std::endl;
			if (T > 2100 && T < 2300) color1 << "red" << std::endl;
			if (T > 1600 && T < 2100) color1 << "firebrick" << std::endl;
			if (T < 1600)             color1 << "maroon" << std::endl;
      fout1 << resM[i].x << " " << resM[i].y << std::endl;
		}
	}
 /*
if ((SPHsize >= 600 * 9) && (SPHsize <= 605 * 9)) ///
	{
		for (int i = 0; i < SPHsize; i++)
		{
			double T = 168e4 * resM[i].Pr / resM[i].De;
			if (T > 3200)             color2 << "gold" << std::endl;
			if (T > 3000 && T < 3200) color2 << "orange" << std::endl;
			if (T > 2700 && T < 3000) color2 << "coral" << std::endl;
			if (T > 2300 && T < 2700) color2 << "orangered" << std::endl;
			if (T > 2100 && T < 2300) color2 << "red" << std::endl;
			if (T > 1600 && T < 2100) color2 << "firebrick" << std::endl;
			if (T < 1600)             color2 << "maroon" << std::endl;
      fout2 << resM[i].x << " " << resM[i].y << std::endl;
		}
	}
 
 if ((SPHsize >= 900 * 9) && (SPHsize <= 901 * 9)) ///
	{
		for (int i = 0; i < SPHsize; i++)
		{
			double T = 168e4 * resM[i].Pr / resM[i].De;
			if (T > 3200)             color3 << "gold" << std::endl;
			if (T > 3000 && T < 3200) color3 << "orange" << std::endl;
			if (T > 2700 && T < 3000) color3 << "coral" << std::endl;
			if (T > 2300 && T < 2700) color3 << "orangered" << std::endl;
			if (T > 2100 && T < 2300) color3 << "red" << std::endl;
			if (T > 1600 && T < 2100) color3 << "firebrick" << std::endl;
			if (T < 1600)             color3 << "maroon" << std::endl;
      fout3 << resM[i].x << " " << resM[i].y << std::endl;
		}
	}
 
 if (SPHsize <= 300 * 9) ///
	{
		for (int i = 0; i < SPHsize; i++)
		{
			fout << resM[i].x << " " << resM[i].y << std::endl;
		}
	}
 */
 

	int indx(7);
std::cout << resM[indx].z << "\n";
	//cool << resM[indx].En << std::endl;
	//keenetic << (resM[indx].Vx * resM[indx].Vx + resM[indx].Vy * resM[indx].Vy + resM[indx].Vz * resM[indx].Vz) * m / 2 << std::endl;

	// print
  /*
	printf("D = %e\t", resM[indx].De);
	printf("En = %f\t", resM[indx].En);
	printf("Pr = %e\t", resM[indx].Pr);
	printf("T = %f\t", 168e4 * resM[indx].Pr / resM[indx].De);
	printf("(x, y, z) = (%f, %f, %f)\t", resM[indx].x, resM[indx].y, resM[indx].z);
	printf("\n");
  */


	std::memcpy(M, resM, Nparts * sizeof(model));
	delete[] resM;
}