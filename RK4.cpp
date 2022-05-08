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

		bool flag = true;
		Type S = sqrt((temp.x - 0.0806045) * (temp.x - 0.0806045) + temp.y * temp.y + temp.z * temp.z);
		}
	}


if ((SPHsize >= 450 * 9) && (SPHsize <= 455 * 9))
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

	int indx(7);
        std::cout << resM[indx].z << "\n";

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
