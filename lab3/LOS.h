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
}x, y;
real eps = 1e-13;
real lymbda = 1;
real gamma = 1;
int maxiter = 10000;
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
	int* mx;
	int* my;



public:
	double FuncF(double x, double y);
	double FuncU(double x, double y);
	void Read(Data* area, string file, int& m);

	void MatrixIns(vector<real> Ox, vector<real>Oy, int i, int j, vector<double> n1, vector<double> n2, vector<double> di, vector<double>n3,
		vector<double> n4, vector<double> F);

	void MatrixUnIns(vector<real> Ox, vector<real> Oy, int i, int j, vector<double> n1, vector<double> n2, vector<double> di, vector<double>n3,
		vector<double> n4, vector<double> F);

	void MatrixBound(vector<real> Ox, vector<real> Oy, int i, int j);

	void BuildGrid(Data* S, vector<real> res);

	void Iteration(vector<real> xk, vector<real> xknext, double w, int nx, int ny);

	void Multiply(vector<real> x, int nx, int ny, int tmp);

	void Zeidel(double w, vector<real> x, int nx, int ny, int tmp);

	double Norm(vector<real> vec, int n);

	double Addition(int i, vector<real> x, int kx, int tmp);

	void GaussSeidel(vector<real> x, vector<real> f, double w, int kx, int tmp);
	void Output(vector<real> x0, string file);
	void AreaUn(double a, double b, double k, int n, vector<real> res, int& q, vector<real>& mxy);
	void BuildGridUn(Data* S, vector<real> res);
	void CheckPoint(Data* x, Data* y, int Ox_size, int Oy_size, vector<int> check, vector<vector<real>> view,int &tmp);
	void CheckArea(Data x, Data y,vector<real> Ox, vector<real> Oy, int c, vector<double> n1, vector<double> n2, vector<double> di, vector<double>n3,
		vector<double> n4, vector<double> F, vector<int> check, int tmp);


};

