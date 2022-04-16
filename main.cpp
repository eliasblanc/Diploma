#include "Instances.h"

//#define __CRTDBG_MAP_ALLOC
//#include <crtdbg.h>
//#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#define new DEBUG_NEW

std::ofstream fout;
std::ofstream color;
std::ofstream cool;
std::ofstream notcool;
std::ofstream keenetic;

Type mu    = 0.919395;
Type L     = -0.647241;
Type h     = 0.003;
Type PI    = 3.1416;
Type Sh    = 0.015;
Type KSTL  = 1.5;//1.1
Type KSTL2 = 1.2;//1

const int dim = 9;

int SPHsize;
const int Nparts(1200);
Type tVelo(0.1);

const Type Asph(1.0);
const Type Bsph(2.0);
const Type NUsph(0.01);
const Type CS(0.0153);

int Nfall(0);
int Nspdlim(0);
int area(0);

Type m = 18e-21;                          // СИ
Type dencity = 176549e-19;                // Старый темп
Type pressure = 315e-19;                  //
Type energy = 0.00267986;                 //

//Type m = 18e-21;                         // СГС
//Type dencity = 176549e-19;               // Старый темп
//Type pressure = 158e-19;                 //
//Type energy = 0.00133912;                //

//Type m = 18e-18;                         // СГС
//Type dencity = 1765e-14;                 // Новый темп
//Type pressure = 1576e-17;                //
//Type energy = 0.00133912;                //

Type G(667e-13);
Type Fgrav0(1e14);

model* M  = new model[Nparts];
Type* lnD = new Type[93];
Type* lnT = new Type[353];
Type* lnL = new Type[32829];

int main()
{
	fout.open("royal.txt", std::ios_base::out | std::ios_base::trunc);
	color.open("color.txt", std::ios_base::out | std::ios_base::trunc);
	cool.open("cool.txt", std::ios_base::out | std::ios_base::trunc);
	notcool.open("notcool.txt", std::ios_base::out | std::ios_base::trunc);
	keenetic.open("keenetic.txt", std::ios_base::out | std::ios_base::trunc);

	ReadFile();

	SPHsize = 9;

	InitConditions_S(M);

	for (size_t i = 0; i < 9; i++)
	{
		fout << L << " " << 0 << std::endl;
	}

	Type t = -omp_get_wtime();
	for (size_t i = 0; i < Nparts; i++)
	{
		if (i % 10 == 0 && i != 0)
		{
			SPHsize += 9;
			std::cout << i << std::endl;
		}
		RK4_S();
	};

	std::cout << "\n----------------\ntime = " << t + omp_get_wtime() \
		<< "\nFalls = " << Nfall << "\nSpd lmts = " << Nspdlim << "\n Area lost = " << area;

	keenetic.close();
	keenetic.clear();
	notcool.close();
	notcool.clear();
	cool.close();
	cool.clear();
	fout.close();
	fout.clear();
	color.close();
	color.clear();

	delete[] M;
	delete[] lnD;
	delete[] lnT;
	delete[] lnL;

	//_CrtDumpMemoryLeaks();

	return 0;
} 
