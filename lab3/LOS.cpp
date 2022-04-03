#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include<functional>
#include <iomanip>
#include"LOS.h"



double FuncF(double x, double y)
{
	return 2 * x + 3;
}

double FuncU(double x, double y)
{
	return 2 * x + 3;
}
// считывание из входного файла
void Read(Data& area, string file, int& m)

{
	//area.nodes = new vector<real>;

	ifstream fin(file + ".txt");
	fin >> m;
	area.nodes.resize(m, 0);
	area.k.resize(m - 1, 0);
	// количество точек для разбиений вместе с границами на отрезке
	area.n.resize(m - 1, 0);
	// данные в двух разных структурах x,y
	for (int i = 0; i < m; i++)
		fin >> area.nodes[i];
	for (int i = 0; i < m - 1; i++)
		fin >> area.k[i];
	for (int i = 0; i < m - 1; i++)
		fin >> area.n[i];
	fin.close();
}
// внутренние узлы в равномерной матрице
void MatrixIns(vector<double> Ox, vector<double> Oy, int i, int j, vector<double> n1, vector<double> n2, vector<double> di, vector<double>n3,	
	vector<double> n4, vector<double> F)
{
	double hx = Ox[i + 1] - Ox[i];
	double hy = Oy[j + 1] - Oy[j];
	int row = j * Ox.size() + i;
	n1[row] = -lymbda / (hy * hy);
	n2[row] = -lymbda / (hx * hx);
	di[row] = 2 * lymbda * (1 / (hx * hx) + 1 / (hy * hy)) + gamma;
	n3[row - 1] = -lymbda / (hx * hx);
	n4[row - Ox.size()] = -lymbda / (hy * hy);
	F[row] = FuncF(Ox[i], Oy[j]); //есть ещё funcU
}
// внутренние узлы в неравномерной матрице
void MatrixUnIns(vector<double> Ox, vector<double> Oy, int i, int j, vector<double> n1, vector<double> n2, vector<double> di, vector<double>n3,
	vector<double> n4, vector<double> F)
{
	double hx = Ox[i + 1] - Ox[i];
	double hy = Oy[j + 1] - Oy[j];
	double hx_pred = Ox[i] - Ox[i - 1];
	double hy_pred = Oy[j] - Oy[j - 1];
	int row = j * Ox.size() + i;
	n1[row] = -2 * lymbda / (hy * (hy + hy_pred));
	n2[row] = -2 * lymbda / (hx * (hx + hx_pred));
	di[row] = 2 * lymbda * (1 / (hx * hx_pred) + 1 / (hy * hy_pred)) + gamma;
	n3[row - 1] = -2 * lymbda / (hx_pred * (hx + hx_pred));
	n4[row - Ox.size()] = -2 * lymbda / (hy_pred * (hy + hy_pred));
	F[row] = FuncF(Ox[i], Oy[j]);
}
// краевые условия первого рода
void MatrixBound(vector<double>Ox, vector<double>Oy, int i, int j, vector<double>& F)
{
	int glob = j * Ox.size() + i;
	F[glob] = FuncU(Ox[i], Oy[j]);
}

void BuildGrid(Data* S, vector<double>& res, vector<double>& mxy)
{
	double h, coord;

	int q = 0, i;
	for (int j = 0; j < S.n.size(); j++)
	{
		double b = S.nodes[j + 1];
		coord = S.nodes[j];
		//делим на количество областей
		h = (b - coord) / (S.n[j] - 1);
		for (i = 0; i < S.n[j] - 1; i++)
		{
			res.push_back(coord);
			coord += h;
		}
		q += 1;
		if (q == 1)
			mxy.push_back(i);
	}
	res.push_back(S.nodes[S.nodes.size() - 1]);
}
// подсчет итерация
void Iteration(vector<double> xk, vector<double>& xknext, double w, int nx, int ny, vector<double>& n1, vector<double>& n2, vector<double>& di, vector<double>& n3,
	vector<double>& n4, vector<double>& F)
{
	int i = 0;
	int m = nx;
	double sum, w1;
	int n = nx * ny;
	for (i; i < n; i++)
	{
		int j = 0;
		sum = 0;
		if (i > 0) sum += n3[i - 1] * xk[i - 1];
		if (i >= m) sum += n4[i - m] * xk[i - m];
		sum += di[i] * xk[i];
		if (i < n - 1) sum += n2[i] * xk[i + 1];
		if (i < n - m) sum += n1[i] * xk[i + m];
		sum = F[i] - sum;
		w1 = w / di[i];
		xknext[i] = xk[i] + w1 * sum;
	}
}
// нахождение размерности
void Multiply(vector<double>x, int nx, int ny, vector<double>& n1, vector<double>& n2, vector<double>& di, vector<double>& n3,
	vector<double>& n4, vector<int>& result)
{
	int i = 0;
	double sum;
	int m = nx;
	int n = nx * ny;
	for (i; i < tmp; i++)
	{
		sum = 0;
		if (i > 0) sum += n3[i - 1] * x[i - 1];
		if (i >= m) sum += n4[i - m] * x[i - m];
		sum += di[i] * x[i];
		if (i < n - 1) sum += n2[i] * x[i + 1];
		if (i < n - m) sum += n1[i] * x[i + m];
		result[i] = sum;
	}
}
//  метод Зейделя
void Zeidel(double w, vector<double>& x, int nx, int ny, vector<double> xk, vector<double>& xknext, vector<double>& n1, vector<double>& n2, vector<double>& di, vector<double>& n3,
	vector<double>& n4, vector<double>& F, vector<int>& result)
{
	int q = 0;
	int i;
	x.resize(tmp);
	result.resize(tmp);
	double nev, norm_b = 0, norm, d;
	for (i = 0; i < tmp; i++)
	{
		x[i] = 0;
		norm_b += F[i] * F[i];
	}
	nev = 1;
	while (q < maxiter && nev < eps)
	{
		Iteration(x, x, w, nx, ny ,xk ,xknext,n1,n2,n4,F);
		Multiply(x, nx, ny, n1, n2, di, n3, n4, result);
		norm = 0;
		for (i = 0; i < tmp; i++)
		{
			d = F[i] - result[i];
			norm += d * d;
		}
		nev = sqrt(norm / norm_b);
		q++;
	}
}
// нахождение нормы
double Norm(vector <double>& vec, int n)
{
	double norma = 0;
	for (int i = 0; i < n; i++)
		norma += vec[i] * vec[i];
	return sqrt(norma);
}

double Addition(int i, vector<double>& x, int kx, vector<double>& n1, vector<double>& n2, vector<double>& di, vector<double>& n3,
	vector<double>& n4)
{
	double sum = di[i] * x[i];
	if (i < tmp - 1) sum += n2[i] * x[i + 1];
	if (i < tmp - kx) sum += n1[i] * x[i + kx];
	if (i >= 1) sum += n3[i - 1] * x[i - 1];
	if (i >= kx) sum += n4[i - kx] * x[i - kx];
	return sum;
}

void GaussSeidel(vector<double>& x, vector<double>& f, double w, int kx, vector<double>& n1, vector<double>& n2, vector<double>& di, vector<double>& n3,
	vector<double>& n4)
{
	x.resize(tmp);
	double nev = 0;
	int n = tmp;
	double norma_b = 0, error = 0, errorx = 0;
	for (int i = 0; i < n; i++)
	{
		norma_b += f[i] * f[i];
		double sum = f[i] - Addition(i, x, kx, n1,n2,di,n3,n4);
		nev += sum * sum;
	}
	norma_b = sqrt(norma_b);
	nev = sqrt(nev) / norma_b;
	for (int itr = 0; nev > eps && itr < maxiter; itr++)
	{
		nev = 0;
		for (int i = 0; i < n; i++)
		{
			double sum = f[i] - Addition(i, x, kx, n1, n2, di, n3, n4);
			x[i] = x[i] + w * sum / di[i];
			nev += sum * sum;
		}
		nev = sqrt(nev) / norma_b;
	}
}
// вывод в файл
void Output(vector<double>x0, string file)
{
	ofstream fout(file + ".txt");
	for (int i = 0; i < x0.size(); i++)
		fout << setprecision(15) << x0[i] << endl;
}

void AreaUn(double a, double b, double k, int n, vector<double>& res, int& q)
{
	double h0, h, coord;
	coord = a;
	//k^n, n-к-во разбиений
	h0 = (b - a) * (1 - k) / (1 - pow(k, n - 1));
	h = h0;
	q += 1;
	int i;
	for (i = 1; i < n - 1; i++)
	{
		coord += h;
		res.push_back(coord);
		h *= k;
	}
	if (q == 1)
		mxy.push_back(i);
}

void BuildGridUn(Data *S, vector<double>& res)
{
	int k = 0;
	for (int i = 0; i < S.n.size(); i++)
	{
		res.push_back(S.nodes[i]);
		AreaUn(S.nodes[i], S.nodes[i + 1], S.k[i], S.n[i], res, k);
	}
	res.push_back(S.nodes[S.nodes.size() - 1]);
}
// обнуление фиктивных узлов
void CheckPoint(Data *x, Data *y, int Ox_size, int Oy_size, vector<int>& check, vector<vector<real>>& view)
{
	tmp = Ox_size * Oy_size;
	// 1 - узел существует, не фиктивный
	check.resize(tmp, 1);
	for (int l = 0; l < Oy_size; l++)
	{
		for (int i = 0, q = 0; i < Ox_size; i++, q++)
		{
			if (view[l][q] == 0)
			{
				int row = l * Ox_size + i;
				check[row] = 0;
			}
		}
	}
}

void CheckArea(vector<double>Ox, vector<double>Oy, int c, vector<double>& n1, vector<double>& n2, vector<double>& di, vector<double>& n3,
	vector<double>& n4, vector<double>& F, vector<int>& check)
{
	int nx = Ox.size();
	int ny = Oy.size();
	di.resize(tmp, 1);
	n1.resize(tmp - nx + 1, 0);
	n2.resize(tmp - 1, 0);
	n3.resize(tmp - 1, 0);
	n4.resize(tmp - nx + 1, 0);
	F.resize(tmp, 0);
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			int flagx = 1;
			int flagy = 1;
			int row = j * nx + i;
			// пропуск фиктивных узлов
			if (check[row] == 0) continue;

			//обработка внутренних узлов
			for (int q = 0; q < x.nodes.size(); q++)
			{
				if ((Ox[i]) == x.nodes[q]) flagx = 0;
			}
			for (int q = 0; q < y.nodes.size(); q++)
			{
				if ((Oy[j]) == y.nodes[q]) flagy = 0;
			}
			if ((flagx == 1 || flagy == 1) && Ox[i] != Ox[0] && Ox[i] != Ox[Ox.size() - 1] && Oy[j] != Oy[0] && Oy[j] != Oy[Oy.size() - 1])
			{
				if (check[row - 1] == 1 && check[row + 1] == 1)
				{
					if (c == 1)
						// для неравномерной сетки
						MatrixUnIns(Ox, Oy, i, j, n1,  n2, di,  n3, n4, F);
					else if (c == 2)
						// для равномерной сетки 
						MatrixIns(Ox, Oy, i, j, n1, n2, di, n3, n4, F);
				}
			}
			else
			{
				// обработка первых граничных условий
				MatrixBound(Ox, Oy, i, j, F);
				continue;
			}
		}
	}
}

