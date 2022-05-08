#include "Instances.h"

//#define __CRTDBG_MAP_ALLOC
//#include <crtdbg.h>
//#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#define new DEBUG_NEW

std::ofstream fout;
std::ofstream fout1;
std::ofstream fout2;
std::ofstream fout3;
std::ofstream color;
std::ofstream color1;
std::ofstream color2;
std::ofstream color3;
std::ofstream cool;
std::ofstream notcool;
std::ofstream keenetic;

Type mu    = 0.919395;
Type L     = -0.647241;
Type h     = 0.003;
Type PI    = 3.1416;
Type Sh    = 0.015;
Type KSTL  = 1;//1.1
Type KSTL2 = 0.45;//1

const int dim = 9;

int SPHsize;
const int Nparts(4700);
Type tVelo(0.1);

const Type Asph(1.0);
const Type Bsph(2.0);
const Type NUsph(0.01);
const Type CS(0.0153);

int Nfall(0);
int Nspdlim(0);
int area(0);

Type m = 18e-21;                          // ��
Type dencity = 176549e-19;                // ������ ����
Type pressure = 315e-19;                  //
Type energy = 0.00267986;                 //

//Type m = 18e-21;                         // ���
//Type dencity = 176549e-19;               // ������ ����
//Type pressure = 158e-19;                 //
//Type energy = 0.00133912;                //

//Type m = 18e-18;                         // ���
//Type dencity = 1765e-14;                 // ����� ����
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
  fout1.open("royal1.txt", std::ios_base::out | std::ios_base::trunc);
  fout2.open("royal2.txt", std::ios_base::out | std::ios_base::trunc);
  fout3.open("royal3.txt", std::ios_base::out | std::ios_base::trunc);
  
	color.open("color.txt", std::ios_base::out | std::ios_base::trunc);
  color1.open("color1.txt", std::ios_base::out | std::ios_base::trunc);
  color2.open("color2.txt", std::ios_base::out | std::ios_base::trunc);
  color3.open("color3.txt", std::ios_base::out | std::ios_base::trunc);
 
	//cool.open("cool.txt", std::ios_base::out | std::ios_base::trunc);
	//notcool.open("notcool.txt", std::ios_base::out | std::ios_base::trunc);
	//keenetic.open("keenetic.txt", std::ios_base::out | std::ios_base::trunc);

	//ReadFile();

	SPHsize = 9;

	InitConditions(M);

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
			//std::cout << i << " " << t + omp_get_wtime() << std::endl;
		}
		RK4();
	};

	std::cout << "\n-------\ntime = " << t + omp_get_wtime() \
		<< "\nFalls = " << Nfall << "\nSpd lmts = " << Nspdlim \
		<< "\nNparts = " << Nparts;

	//keenetic.close();
	//keenetic.clear();
	//notcool.close();
	//notcool.clear();
	//cool.close();
	//cool.clear();
 
	fout.close(); 
	fout.clear();
  fout1.close();
	fout1.clear();
  fout2.close();
	fout2.clear();
  fout3.close();
	fout3.clear();
	color.close();
	color.clear();
  color1.close();
	color1.clear();
  color2.close();
	color2.clear();
  color3.close();
	color3.clear();

	delete[] M;
	delete[] lnD;
	delete[] lnT;
	delete[] lnL;

	//_CrtDumpMemoryLeaks();

	return 0;
} 
