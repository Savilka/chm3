#pragma once
#include <iostream>
#include <fstream>

using namespace std;

typedef double real;

class LOS {
   int N;
   int maxiter;
   double eps, nev;
   real* ggl, * ggu, * di, * pr;
   int* ig, * jg;
   real* r, * z, * p, * x, * res;
   real alpha, beta;
public:
   void input(ifstream& kuslau, ifstream& ig, ifstream& jg, ifstream& ggl, ifstream& ggu, ifstream& di, ifstream& pr);
   void output(ofstream& output);
   void los();
   void los_LU();
private:
   void LU();
   void r0();
   void xk();
   void rk();
   void zk();
   void pk();
   void calcAlpha();
   void calcBeta();
   void copyVec(real* a, real* b, int size);
   void matMul(real* vec);
   real norm(real* vec);
   real scalar(real* a, real* b, int size);
};

