#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include<functional>
#include <iomanip>

using namespace std;
typedef double real;
real eps = 1e-13;

struct Data
{
	// количество узлов
	real nodes;
	vector<double> k;
	// количество точек для разбиения
	vector<double>n;
} x, y;

class lab {


	real lymbda, gamma;
	int tmp, maxiter ;
	real** view;
	real* result;
	real* F;
	int* check;
	real* n1;
	real* n2;
	real* n3;
	real* n4;
	real* di;
	int* mx;
	int* my;
public:
	double FuncF(double x, double y);
	double FuncU(double x, double y);
	void Read(Data& area, string file, int& m);
	void MatrixIns(real* Ox, real* Oy, int i, int j);
	void MatrixUnIns(real* Ox, real* Oy, int i, int j);
	void MatrixBound(real* Ox, real* Oy, int i, int j);
	void BuildGrid(Data S, real* res);
	void Iteration(real* xk, real* xknext, double w, int nx, int ny);
	void Multiply(real* x, int nx, int ny);
	void Zeidel(double w, real* x, int nx, int ny);
	double Norm(real* vec, int n);
	double Addition(int i, real* x, int kx);
	void GaussSeidel(real* x, real* f, double w, int kx);
	void Output(real* x0, string file);
	void AreaUn(double a, double b, double k, int n, real* res, int& q);
	void BuildGridUn(Data S, real* res);
	void CheckPoint(Data x, Data y, int Ox_size, int Oy_size);
	void CheckArea(real* Ox, real* Oy, int c);


};