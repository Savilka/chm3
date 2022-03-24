#include <iostream>
#include "LOS.h"

using namespace std;

int main() {
	lab la;
	int m1, m2;
	la.Read(x, "X", m1);
	la.Read(y, "Y", m2);
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
		la.BuildGrid(x, la.Ox);
		// построили неравномерную прямоугольную сетку
		la.BuildGridUn(y, la.Oy);
		la.mx = mxy[0]; my = mxy[1];
		nx = Ox.size(); ny = Oy.size();
		view.resize(ny);//x
		for (int i = 0; i < ny; i++)
		{
			view[i].resize(nx);//y
		}
		ifstream fin;
		fin.open("input.txt");
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				fin >> view[j][i];
		fin.close();
		// выявили фиктивные вершины
		CheckPoint(x, y, nx, ny);
		// обрабатываем граничную и внутреннюю области
		CheckArea(Ox, Oy, c);
		GaussSeidel(x0, F, w, Ox.size());
		Output(x0, "output");
		Output(Ox, "Ox");
		Output(Oy, "Oy");
		break;
	}
	case 2:
	{
		// построили равномерную прямоугольную сетку
		BuildGrid(x, Ox);
		BuildGrid(y, Oy);
		mx = mxy[0]; my = mxy[1];
		nx = Ox.size(); ny = Oy.size();
		view.resize(ny);//x
		for (int i = 0; i < ny; i++)
		{
			view[i].resize(nx);//y
		}
		ifstream fin;
		fin.open("input.txt");
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				fin >> view[j][i];
		fin.close();
		// выявили фиктивные вершины
		CheckPoint(x, y, nx, ny);
		// обрабатываем граничную и внутреннюю области
		CheckArea(Ox, Oy, c);
		//решение слау методом Гаусса-Зейделя
		GaussSeidel(x0, F, w, Ox.size());
		Output(x0, "output");
		Output(Ox, "Ox");
		Output(Oy, "Oy");
		break;
	}
	}
}