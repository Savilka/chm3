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
	cout << "������� 1 ��� ������������� �����" << endl;
	cout << "������� 2 ��� ����������� �����" << endl;
	cin >> c;
	switch (c)
	{
	case 1:
	{
		la.BuildGrid(&x, Ox);
		// ��������� ������������� ������������� �����
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
		// ������� ��������� �������
		la.CheckPoint(&x, &y, nx, ny, la.check, la.view);
		// ������������ ��������� � ���������� �������
		la.CheckArea(Ox, Oy, c, la.n1, la.n2,  la.di, la.n3,
			 la.n4,  la.F, la.check);
		la.GaussSeidel( x0 , la.F, w, Ox.size());
		la.Output(x0, "output");
		la.Output(Ox, "Ox");
		la.Output(Oy, "Oy");
		break;
	}
	case 2:
	{
		// ��������� ����������� ������������� �����
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
		// ������� ��������� �������
		la.CheckPoint(&x, &y, nx, ny, la.check, la.view);
		// ������������ ��������� � ���������� �������
		la.CheckArea(Ox, Oy, c, la.n1, la.n2, la.di, la.n3,
			la.n4, la.F,la.check);
		//������� ���� ������� ������-�������
		la.GaussSeidel(x0, la.F, w, Ox.size());
		la.Output(x0, "output");
		la.Output(Ox, "Ox");
		la.Output(Oy, "Oy");
		break;
	}
	}
}