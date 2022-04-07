#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include<functional>
#include <iomanip>


using namespace std;
typedef double real;
struct Data
{
	// количество узлов
	vector<real> nodes;
	vector<real> k;
	// количество точек для разбиения
	vector<real> n;
};

class lab {

public:
	vector<vector<real>> view;
	vector<real> result;
	vector<real> F;
	vector<int> check;
	vector<real> n1;
	vector<real> n2;
	vector<real> n3;
	vector<real> n4;
	vector<real> di;
    vector<int> mxy;
	int mx;
	int my;

	real eps = 1e-13;
	real lymbda = 1;
	real gamma = 1;
	int maxiter = 10000;

	Data x, y;



public:
	real FuncF(real x, real y);
	real FuncU(real x, real y);
	void Read(Data& area, string file, int& m);

	void MatrixIns(vector<real> Ox, vector<real>Oy, int i, int j);

	void MatrixUnIns(vector<real> Ox, vector<real> Oy, int i, int j);

	void MatrixBound(vector<real> Ox, vector<real> Oy, int i, int j);

	void BuildGrid(Data& S, vector<real>& res);

	void Iteration(vector<real> xk, vector<real>& xknext, real w, int nx, int ny);

	void Multiply(vector<real>x, int nx, int ny, int tmp);

	void Zeidel(real w, vector<real>& x, int nx, int ny, vector<real> xk, vector<real>& xknext, int tmp);

	real Norm(vector <real>& vec, int n);

	real Addition(int i, vector<real>& x, int kx, int tmp);

	void GaussSeidel(vector<real>& x, vector<real>& f, real w, int kx, int tmp);
	void Output(vector<real> x0, string file);
	void AreaUn(real a, real b, real k, int n, vector<real>& res, int& q);
	void BuildGridUn(Data& S, vector<real>& res);
	void CheckPoint(Data& x, Data& y, int Ox_size, int Oy_size, int &tmp);
	void CheckArea(Data& x, Data& y,vector<real> Ox, vector<real> Oy, int c, int tmp);


};

