#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

typedef double real;

class LOS {
   int N;
   int maxiter;
   double eps, nev;
   real* ggl, * ggu, * di, * pr;
   int* ig, * jg;
   real* r, * z, * p, * x;
   real alpha, beta;
public:
   void input(ifstream& kuslau, ifstream& ig, ifstream& jg, ifstream& ggl, ifstream& ggu, ifstream& di, ifstream& pr);
   void output(ofstream& output);
   void los(ofstream& iteration);
   void los_d(ofstream& iteration);
   void los_LU(ofstream& iteration);
private:
   real LUdirect(real* y);
   real LUreverse(real* res);
   void LU();
   void r0();
   void xk();
   void rk();
   void zk();
   void pk(real* Ak);
   void calcAlpha();
   void calcBeta(real* Ar);
   void copyVec(real* a, real* b, int size);
   void matMul(real* vec, real* res);
   void vecDivD(real* res);
   real norm(real* vec);
   real scalar(real* a, real* b, int size);
};

