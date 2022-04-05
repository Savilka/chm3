#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include<functional>
#include <iomanip>

#include "LOS.h"

using namespace std;

int main() {
	lab la;
	int tmp = 0;
	vector<real> Ox;
	vector<real> Oy;
	vector<real> x0;
	int mx;
	int my;
	vector<real> mxy;
	int m1, m2;
	la.Read(&x, "X", m1);
	la.Read(&y, "Y", m2);
	int nx, ny;
	int c;
	double w = 1.6;
	setlocale(LC_ALL, "Russian");
	cout << "Нажмите 1 для неравномерной сетки" << endl;
	cout << "Нажмите 2 для равномерной сетки" << endl;
	cin >> c;
	switch (c)
	{
	case 1:
	{
		la.BuildGrid(&x, Ox);
		// построили неравномерную прямоугольную сетку
		la.BuildGridUn(&y, Oy);
		mx = mxy[0]; my = mxy[1];
		nx = Ox.size(); ny = Oy.size();
		la.view.resize(ny);//x
		for (int i = 0; i < ny; i++)
		{
			la.view[i].resize(nx);//y
		}
		ifstream fin;
		fin.open("input.txt");
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				fin >> la.view[j][i];
		fin.close();
		// выявили фиктивные вершины
		la.CheckPoint(&x, &y, nx, ny, la.check, la.view, tmp);
		// обрабатываем граничную и внутреннюю области
		la.CheckArea(x,y,Ox, Oy, c, la.n1, la.n2,  la.di, la.n3,
			 la.n4,  la.F, la.check, tmp);
		la.GaussSeidel( x0 , la.F, w, Ox.size(), tmp);
		la.Output(x0, "output");
		la.Output(Ox, "Ox");
		la.Output(Oy, "Oy");
		break;
	}
	case 2:
	{
		// построили равномерную прямоугольную сетку
		la.BuildGrid(&x, Ox);
		la.BuildGrid(&y, Oy);
		mx = mxy[0];*la.my = mxy[1];
		nx = Ox.size(); ny = Oy.size();
		la.view.resize(ny);//x
		for (int i = 0; i < ny; i++)
		{
			la.view[i].resize(nx);//y
		}
		ifstream fin;
		fin.open("input.txt");
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				fin >> la.view[j][i];
		fin.close();
		// выявили фиктивные вершины
		la.CheckPoint(&x, &y, nx, ny, la.check, la.view,tmp);
		// обрабатываем граничную и внутреннюю области
		la.CheckArea(x,y, Ox, Oy, c, la.n1, la.n2, la.di, la.n3,
			la.n4, la.F,la.check,tmp);
		//решение слау методом Гаусса-Зейделя
		la.GaussSeidel(x0, la.F, w, Ox.size(),tmp);
		la.Output(x0, "output");
		la.Output(Ox, "Ox");
		la.Output(Oy, "Oy");
		break;
	}
	}
}