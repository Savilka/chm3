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
	int m1, m2;
	la.Read(la.x, "x", m1);
	la.Read(la.y, "y", m2);
	int nx, ny;
	int c;
	real w = 1.6;
	setlocale(LC_ALL, "Russian");
	cout << "������� 1 ��� ������������� �����" << endl;
	cout << "������� 2 ��� ����������� �����" << endl;
	cin >> c;
	switch (c)
	{
	case 1:
	{
		la.BuildGrid(la.x, Ox);
		// ��������� ������������� ������������� �����
		la.BuildGridUn(la.y, Oy);
		la.mx = la.mxy[0]; la.my = la.mxy[1];
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
		la.CheckPoint(la.x, la.y, nx, ny, tmp);
		// ������������ ��������� � ���������� �������
		la.CheckArea(la.x, la.y,Ox, Oy, c, tmp);
		la.GaussSeidel( x0 , la.F, w, Ox.size(), tmp);
		la.Output(x0, "output");
		la.Output(Ox, "Ox");
		la.Output(Oy, "Oy");
		break;
	}
	case 2:
	{
		// ��������� ����������� ������������� �����
		la.BuildGrid(la.x, Ox);
		la.BuildGrid(la.y, Oy);
		la.mx = la. mxy[0]; la.my = la.mxy[1];
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
		la.CheckPoint(la.x, la.y, nx, ny, tmp);
		// ������������ ��������� � ���������� �������
		la.CheckArea(la.x,la.y, Ox, Oy, c, tmp);
		//������� ���� ������� ������-�������
		la.GaussSeidel(x0, la.F, w, Ox.size(),tmp);
		la.Output(x0, "output");
		la.Output(Ox, "Ox");
		la.Output(Oy, "Oy");
		break;
	}
	}
}