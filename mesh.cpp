//#include "pch.h"
//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>

using namespace std;

void file_read(const char name[50], int& N, int& NP, int& NE1, int& NE2, int& type, vector<double>& x, vector<double>& y, vector<double>& z)
{
	N = 1;
	string str;
	ifstream TEMP(name);
	if (!TEMP.is_open())
		cerr << "Файл не открыт... ";

	TEMP >> NP; TEMP >> NP;
	getline(TEMP, str);
	istringstream inf(str);
	while (inf >> str)
		N++;

	ifstream ifs(name);
	if (!ifs.is_open())
		cerr << "Файл не открыт... ";
	ifs >> NP;

	if (N == 1)
	for (int i = 0; i < NP; i++)
	{
		double temp;
		ifs >> temp;
		x.push_back(temp);
		y.push_back(0);
		z.push_back(0);
	}
	else if (N == 2)
	for (int i = 0; i < NP; i++)
	{
		double temp;
		ifs >> temp;
		x.push_back(temp);
		ifs >> temp;
		y.push_back(temp);
		z.push_back(0);
	}
	else
	for (int i = 0; i < NP; i++)
	{
		double temp;
		ifs >> temp;
		x.push_back(temp);
		ifs >> temp;
		y.push_back(temp);
		ifs >> temp;
		z.push_back(temp);
	}

	ifs >> NE1;
	if (NP == 4)
		ifs >> NE2;
	else NE2 = 0;
	ifs >> type;

}
void read_check(int& NP, int& NE1, int& NE2, int& type, vector<double>& x, vector<double>& y, vector<double>& z)
{
	cout << "NP = " << NP << endl;
	cout << "NE1 = " << NE1 << endl;
	cout << "NE2 = " << NE2 << endl;
	cout << "type = " << type << endl;

	for (int i = 0; i < NP; i++)
	{
		cout << "x[" << i << "] = " << x[i] << " ";
		cout << "y[" << i << "] = " << y[i] << " ";
		cout << "z[" << i << "] = " << z[i] << endl;
	}
		
}
void dots_check(int NE1, int NE2, vector<double>& r)
{
	ofstream f("dots.dat");
	for (int i = 0; i < r.size()/3; i++)
		f << i << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;
	f.close();
}
void dots_check2(int NE1, vector<double>& r)
{
	ofstream f("dots2.dat");
	for (int i = 0; i < r.size() / 3; i++)
		f << i << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;
	f.close();
}
void find_middle(vector<int>& d, vector<double>& r)
{
	double x1, x2, y1, y2, z1, z2, a1, a2, b1, b2, c1, c2;
	double x, y, z;

	x1 = r[(d[0] - 1) * 3];
	y1 = r[(d[0] - 1) * 3 + 1];
	z1 = r[(d[0] - 1) * 3 + 2];

	x2 = r[(d[2] - 1) * 3];
	y2 = r[(d[2] - 1) * 3 + 1];
	z2 = r[(d[2] - 1) * 3 + 2];

	a1 = r[(d[1] - 1) * 3];
	b1 = r[(d[1] - 1) * 3 + 1];
	c1 = r[(d[1] - 1) * 3 + 2];

	a2 = r[(d[3] - 1) * 3];
	b2 = r[(d[3] - 1) * 3 + 1];
	c2 = r[(d[3] - 1) * 3 + 2];

	x = (a1 * (b2 * (x1 - x2) + x2 * y1 - x1 * y2) + a2 * (b1 * (-x1 + x2) - x2 * y1 + x1 * y2)) \
		/ (-(b1 - b2) * (x1 - x2) + (a1 - a2) * (y1 - y2));

	y = (a1 * b2 * (y1 - y2) + a2 * b1 * (-y1 + y2) + (b1 - b2) * (x2 * y1 - x1 * y2))\
		/ (-(b1 - b2) * (x1 - x2) + (a1 - a2) * (y1 - y2));

	z = ((b1 - b2) * (x2 * z1 - x1 * z2) + a1 * (b2 * z1 - y2 * z1 - b2 * z2 + y1 * z2)\
		+ a2 * (y2 * z1 - y1 * z2 + b1 * (-z1 + z2))) / (-(b1 - b2) * (x1 - x2) + (a1 - a2) * (y1 - y2));

	r.push_back(x);
	r.push_back(y);
	r.push_back(z);
}

void mesh2_1(int NE1, int type, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& r)
{
	int M = NE1 * type + 1;   //кол. точек, которые будут в векторе r

	r.push_back(x[0]); r.push_back(y[0]); r.push_back(z[0]);
	double* temp_v = new double[3];

	temp_v[0] = (x[1] - x[0]) / (M - 1);
	temp_v[1] = (y[1] - y[0]) / (M - 1);
	temp_v[2] = (z[1] - z[0]) / (M - 1);

	for (int i = 0; i < M - 2; i++)
	{
		r.push_back((x[0] + temp_v[0] * (i + 1)));
		r.push_back((y[0] + temp_v[1] * (i + 1)));
		r.push_back((z[0] + temp_v[2] * (i + 1)));
	}

	r.push_back(x[1]); r.push_back(y[1]); r.push_back(z[1]);
	dots_check2(NE1, r);

	delete[] temp_v;
}
void mesh2_2(int NE1, int type, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& r)
{
	mesh2_1(NE1, type, x, y, z, r);
	dots_check2(NE1, r);
	/*	double* temp_v = new double[3];
	temp_v[0] = (x[1] - x[0]) / NE1;
	temp_v[1] = (y[1] - y[0]) / NE1;
	temp_v[2] = (z[1] - z[0]) / NE1;

	vector <double>::iterator it;
	int g = 1;

	for (int i = 0; i < NE1 - 1; i++)
	{
	it = r.begin() + i + g;
	r.insert(it, x[0] + 1 / 2 / NE1 * (x[1] - x[0]) + temp_v[0] * i);
	cout << "r" << r[1] << " ";
	r.insert(it, y[0] + 1 / 2 / NE1 * (y[1] - y[0]) + temp_v[0] * i);
	cout << "r" << r[1] << " ";
	r.insert(it, z[0] + 1 / 2 / NE1 * (z[1] - z[0]) + temp_v[0] * i);
	g = 2;
	}

	// рассчитали координаты, а теперь посчитаем значения для файла

	o_NE = NE1;
	o_NP = 2 * NE1 + 1;
	o_NC = 2;

	delete[] temp_v;
	*/
}
void mesh4_123(int NE1, int NE2, int type, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& r)
{
	int r_size = 3 * 2 * (NE1 + NE2) - 1;
	vector<double> ri;  // наши радиус-вектора 
	
	for (int i = 0; i < 4; i++)
	{
		if ((i == 0) || (i == 2))
		{
			ri.push_back((x[i + 1] - x[i]) / NE1);            // 4 радиус-вектора контура
			ri.push_back((y[i + 1] - y[i]) / NE1);
			ri.push_back((z[i + 1] - z[i]) / NE1);
		}
	
		if (i == 1) 
		{
			ri.push_back((x[i + 1] - x[i]) / NE2);            // 4 радиус-вектора контура
			ri.push_back((y[i + 1] - y[i]) / NE2);
			ri.push_back((z[i + 1] - z[i]) / NE2);
		}

		if (i == 3)
		{
			ri.push_back((x[0] - x[3]) / NE2);            // 4 радиус-вектора контура
			ri.push_back((y[0] - y[3]) / NE2);
			ri.push_back((z[0] - z[3]) / NE2);
		}

	}

	for (int i = 0; i < NE1; i++)   // записываем узлы, которые по контуру
	{
		r.push_back(x[0] + ri[0] * i);
		r.push_back(y[0] + ri[1] * i);
		r.push_back(z[0] + ri[2] * i);
	}

	for (int i = 0; i < NE2; i++)   // записываем узлы, которые по контуру
	{
		r.push_back(x[1] + ri[3] * i);
		r.push_back(y[1] + ri[4] * i);
		r.push_back(z[1] + ri[5] * i);
	}

	for (int i = 0; i < NE1; i++)   // записываем узлы, которые по контуру
	{
		r.push_back(x[2] + ri[6] * i);
		r.push_back(y[2] + ri[7] * i);
		r.push_back(z[2] + ri[8] * i);
	}

	for (int i = 0; i < NE2; i++)   // записываем узлы, которые по контуру
	{
		r.push_back(x[3] + ri[9] * i);
		r.push_back(y[3] + ri[10] * i);
		r.push_back(z[3] + ri[11] * i);
	}


	vector<double> ri_temp1(3 * (NE2 - 1));
	vector<double> ri_temp2(3 * (NE2 - 1));

	for (int i = 0; i < NE2 - 1; i++)   // готовимся к созданию р.в. для внутренних точек
	{
		ri_temp1[3*i] = r[r_size - 3*i - 2];
		ri_temp1[3*i + 1] = r[r_size - 3*i - 1];
		ri_temp1[3*i + 2] = r[r_size  - 3*i];
	}

	for (int i = 0; i < NE2 - 1; i++)   // готовимся к созданию р.в. для внутренних точек
	{
		ri_temp2[3 * i] = r[r_size - 2 + 3*(-NE2 + 1 - NE1 - NE2 + 1) + 3*i];
		ri_temp2[3 * i + 1] = r[r_size - 1 + 3*(-NE2 + 1 - NE1 - NE2 + 1) + 3 * i];
		ri_temp2[3 * i + 2] = r[r_size + 3*(-NE2 + 1 - NE1 - NE2 + 1) + 3 * i];
	}

	for (int i = 0; i < 3*(NE2 - 1); i++)      // создаем р.в для внутренних точек
	{
		ri.push_back((ri_temp2[i] - ri_temp1[i]) / NE1);
	}

	for(int j = 0; j < NE2 - 1; j++)
	for (int i = 0; i < NE1 - 1; i++)      // создаем сами внутренние точки
	{
		r.push_back(ri_temp1[0 + j*3] + (i + 1) * ri[12 + j * 3]);
		r.push_back(ri_temp1[1 + j*3] + (i + 1) * ri[13 + j * 3]);
		r.push_back(ri_temp1[2 + j*3] + (i + 1) * ri[14 + j * 3]);
	} 
	
	dots_check(NE1, NE2, r);	
}
void mesh4_4(int NE1, int NE2, int type, vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& r) 
{
	mesh4_123(NE1, NE2, type, x, y, z, r);

	int dN = 2 * (NE1 + NE2) - 1;
	vector<int> d(4);

	for (int i = 0; i < NE1; i++) // в вектор d запишем номера узлов элемента, середину которого хотим отыскать.
	{
		if (i == (NE1 - 1))
		{
			d[0] = i + 1;
			d[1] = i + 2;
			d[2] = i + 3;
			d[3] = i + 1 + dN;
		}
		else
		{
			d[0] = i + 1;
			d[1] = i + 2;
			d[2] = i + 2 + dN;
			d[3] = i + 1 + dN;
		}

		find_middle(d, r);
	}

	for (int j = 0; j < NE2 - 2; j++)
	{
		for (int i = 0; i < NE1; i++)
		{
			if (i == 0)
			{ 
				d[0] = dN - j + 1;
				d[1] = dN + 2 + (NE1 - 1) * j;
				d[2] = dN + 2 + (NE1 - 1) * (j + 1);
				d[3] = dN - j;
			}
			else
				if (i == (NE1 - 1))
				{
					d[0] = (NE1 - 1) * (j + 1) + dN + 1;
					d[1] = NE1 + 2 + j;
					d[2] = NE1 + 3 + j;
					d[3] = (NE1 - 1) * (j + 2) + dN + 1;
				}
				else
				{
					d[0] = i + 1 + dN + (NE1 - 1) * j;
					d[1] = i + 2 + dN + (NE1 - 1) * j;
					d[2] = i + 2 + dN + (NE1 - 1) * (j + 1);
					d[3] = i + 1 + dN + (NE1 - 1) * (j + 1);
				}
			find_middle(d, r);
		}
	}

	for (int i = 0; i < NE1; i++)
	{
		if (i == 0)
		{
			d[0] = dN + 1 - NE2 + 2;
			d[1] = (NE2 - 2) * (NE1 - 1) + 2 + dN;
			d[2] = dN + 1 - NE2;
			d[3] = dN + 1 - NE2 + 1;
		}
		else
			if (i == (NE1 - 1))
			{
				d[0] = dN + NE1 + (NE1 - 1) * (NE2 - 2);
				d[1] = NE1 + NE2;
				d[2] = NE1 + NE2 + 1;
				d[3] = NE1 + NE2 + 2;
			}
			else
			{
				d[0] = i + 1 + dN + (NE1 - 1) * (NE2 - 2);
				d[1] = i + 2 + dN + (NE1 - 1) * (NE2 - 2);
				d[2] = 2 * NE1 + NE2 - i;
				d[3] = 2 * NE1 + NE2 - i + 1;
			}
		find_middle(d, r);

	}

	dots_check(NE1, NE2, r); // check
}
void output_2(const char name[50], int NP, int NE1, int type, vector<double>& r)
{
	ofstream f("output2.txt");
	int o_NE = NE1; // кол. элементов
	int o_NP = NE1 * type + 1; // кол. узлов
	int o_NC = 2; // кол-во контуров 
	int o_ENP; // кол. узлов в эл-те
	//	int o_EN номер элемента

	int M = NE1 * type + 1; //кол. точек, которые будут в векторе r

	f << o_NE << " " << o_NP << " " << o_NC << endl;  //***
	(type == 1) ? (o_ENP = 2) : (o_ENP = 3);  // это пригодится в будущем

	for (int i = 0; i < o_NE; i++)   // цикл для граничных значений
	{
		f << i + 1 << " " << o_ENP << " ";  // выводим номер и количество точек в элементе

		if ((i == 0) && (o_ENP == 2))       // граничный узел слева для  первого типа
			f << 1 << " " << 3 << " " << endl;
		else
		if ((i == 0) && (o_ENP == 3))       // граничный узел слева для  второго типа
			f << 1 << " " << 4 << " " << 3 << " " << endl;
		else
		if ((i == (o_NE - 1)) && (o_ENP == 2))     // граничный узел справа для  первого типа
			f << M << " " << 2 << " " << endl;
		else
		if ((i == (o_NE - 1)) && (o_ENP == 3)) // граничный узел справа для второго типа
			f << M - 1 << " " << 2 << " " << M << " " << endl;

		else
		{
			if (type == 1) {
				for (int j = 1; j < o_ENP + 1; j++)
					f << i * type + j + 1 << " ";
				f << endl;
			}
			else
			{
				f << i * type + 1 + 1 << " ";
				int temp = i * type + 2 + 1;
				f << i * type + 3 + 1 << " " << temp << endl;
			}
		}
	//	cout << "i = " << i << endl;
	
	}

	

	f << 1 << " " << r[0] << " " << r[1] << " " << r[2] << " " << endl;   // опять граничныые точки выводит
	f << 2 << " " << r[M * 3 - 3] << " " << r[M * 3 - 2] << " " << r[M * 3 - 1] << endl;

	for (int i = 1; i < M - 1; i++)  // и все остальные в цикле
		f << i + 2 << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;

	f << 1 << " " << 1 << " " << endl;
	f << 1 << endl;    // это тоже  надо вывести
	f << 2 << endl;
	f.close();
}
void output_4(const char name[50], int NP, int NE1, int NE2, int type, vector<double>& r)
{
	ofstream f(name);
	int xc = (NE1 + 1) * (NE2 + 1) + 1; // номер первого центральног узла
	int o_NE ; // кол. элементов
	if (type == 1) o_NE = NE1 * NE2;
	else if ((type == 2) || (type == 3)) o_NE = NE1 * NE2 * 2;
	else  o_NE = NE1 * NE2 * 4;

	int o_NP; 
	if (type == 4) o_NP = (NE1 + 1) * (NE2 + 1) + NE1 * NE2; 
	else o_NP = (NE1 + 1) * (NE2 + 1); // кол. узлов

	int o_NC = 1; // кол-во контуров 
	int o_ENP; // кол. узлов в эл-те
	//	int o_EN номер элемента

	int dN = 2 * (NE1 + NE2) - 1; // число точек в контуре минус 1 
//	int M = NE1 * type + 1; //кол. точек, которые будут в векторе r

	f << o_NE << " " << o_NP << " " << o_NC << endl;  //*** 
	if ((type == 2) || (type == 3) || (type == 4)) o_ENP = 3; else o_ENP = 4; // это пригодится в будущем

	if (type == 1)
	{
		for (int i = 0; i < NE1; i++)
		{
			if (i == (NE1 - 1))
				f << i + 1 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " " << i + 3 << " " << i + 1 + dN << endl;
			else
				f << i + 1 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " " << i + 2 + dN << " " << i + 1 + dN << endl;
		}

		cout << "dN = " << dN << endl;

		for (int j = 0; j < NE2 - 2; j++)
		{
			for (int i = 0; i < NE1; i++)
			{
				if (i == 0)
					f << (j + 1) * NE1 + (i + 1) << " " << o_ENP << " " << dN - j + 1 << " " << dN + 2 + (NE1 - 1) * j << " " \
					<< dN + 2 + (NE1 - 1) * (j + 1) << " " << dN - j << endl;
				else
					if (i == (NE1 - 1))
						f << (j + 1) * NE1 + (i + 1) << " " << o_ENP << " " << (NE1 - 1) * (j + 1) + dN + 1 << " " << NE1 + 2 + j << " " \
						<< NE1 + 3 + j << " " << (NE1 - 1) * (j + 2) + dN + 1 << endl;
					else
						f << (j + 1) * NE1 + (i + 1) << " " << o_ENP << " " << i + 1 + dN + (NE1 - 1) * j << " " \
						<< i + 2 + dN + (NE1 - 1) * j << " " << \
						i + 2 + dN + (NE1 - 1) * (j + 1) << " " << i + 1 + dN + (NE1 - 1) * (j + 1) << endl;
			}
		}

		for (int i = 0; i < NE1; i++)
		{
			if (i == 0)
				f << (NE2 - 1) * NE1 + i + 1 << " " << o_ENP << " " << dN + 1 - NE2 + 2 << " " << (NE2 - 2) * (NE1 - 1) + 2 + dN << \
				" " << dN + 1 - NE2 << " " << dN + 1 - NE2 + 1 << endl;
			else
				if (i == (NE1 - 1))
					f << (NE2 - 1) * NE1 + i + 1 << " " << o_ENP << " " << dN + NE1 + (NE1 - 1) * (NE2 - 2) << \
					" " << NE1 + NE2 << \
					" " << NE1 + NE2 + 1 << " " << NE1 + NE2 + 2 << endl;
				else
					f << (NE2 - 1) * NE1 + i + 1 << " " << o_ENP << " "\
					<< i + 1 + dN + (NE1 - 1) * (NE2 - 2) << " " << i + 2 + dN + (NE1 - 1) * (NE2 - 2) \
					<< " " << 2 * NE1 + NE2 - i << " " << 2 * NE1 + NE2 - i + 1 << endl;
		}


		// теперь вывод узлов

		for (int i = 0; i < (NE1 + 1) * (NE2 + 1); i++)  // и все остальные в цикле
			f << i + 1 << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;

		f << dN + 1 << endl;
		for (int i = 0; i < dN + 1; i++)
			f << i + 1 << endl;
	}

	if (type == 2) 
	{
		for (int i = 0; i < NE1; i++)
		{
			if (i == (NE1 - 1)) {
				f << 2 * NE1 - 1 << " " << o_ENP << " " << i + 1 << " " << i + 3 << " " << i + 1 + dN << endl;
				f << 2 * NE1 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " " << i + 3  << endl;
			}
			else {
				f << 2*i + 1 << " " << o_ENP << " " << i + 1 << " " << i + 2 + dN << " " << i + 1 + dN << endl;
				f << 2*i + 2 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " " << i + 2 + dN  << endl;
			}
		}

		for (int j = 0; j < NE2 - 2; j++)
		{
			for (int i = 0; i < NE1; i++)
			{
				if (i == 0) {
					f << 2*((j + 1) * NE1 + (i + 1)) - 1 << " " <<\
						o_ENP << " " << dN - j + 1 << " "  \
						<< dN + 2 + (NE1 - 1) * (j + 1) << " " << dN - j << endl;
					f << 2 * ((j + 1) * NE1 + (i + 1)) << " " <<\
						o_ENP << " " << dN - j + 1 << " " << dN + 2 + (NE1 - 1) * j << " " \
						<< dN + 2 + (NE1 - 1) * (j + 1)  << endl;
				}
				else
					if (i == (NE1 - 1)){
						f << 2 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << \
							o_ENP << " " << (NE1 - 1) * (j + 1) + dN + 1 << " " \
							<< NE1 + 3 + j << " " << (NE1 - 1) * (j + 2) + dN + 1 << endl;
						f << 2 * ((j + 1) * NE1 + (i + 1))  << " " << \
							o_ENP << " " << (NE1 - 1) * (j + 1) + dN + 1 << " " << NE1 + 2 + j << " " \
							<< NE1 + 3 + j << endl;
					}
					else {
						f << 2 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << \
							o_ENP << " " << i + 1 + dN + (NE1 - 1) * j << " " << \
							i + 2 + dN + (NE1 - 1) * (j + 1) << " " << i + 1 + dN + (NE1 - 1) * (j + 1) << endl;
						f << 2 * ((j + 1) * NE1 + (i + 1))  << " " << \
							o_ENP << " " << i + 1 + dN + (NE1 - 1) * j << " " \
							<< i + 2 + dN + (NE1 - 1) * j << " " << \
							i + 2 + dN + (NE1 - 1) * (j + 1) << endl;
					}
			}
		}

		for (int i = 0; i < NE1; i++)
		{
			if (i == 0) {
				f << 2 * ((NE2 - 1) * NE1) + 2 * i  + 1<< " " << o_ENP << " " \
					<< dN + 1 - NE2 + 2 << \
					" " << dN + 1 - NE2 << " " << dN + 1 - NE2 + 1 << endl;

				f << 2 * ((NE2 - 1) * NE1) + 2 * i  + 2 << " " << o_ENP << " " \
					<< dN + 1 - NE2 + 2 << " " << (NE2 - 2) * (NE1 - 1) + 2 + dN << \
					" " << dN + 1 - NE2 <<  endl;
			}
			else
				if (i == (NE1 - 1)) {
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 1<< " " << o_ENP << " " << \
						dN + NE1 + (NE1 - 1) * (NE2 - 2) << \
						" " << NE1 + NE2 + 1 << " " << NE1 + NE2 + 2 << endl;
					f << 2 * ((NE2 - 1) * NE1) + 2 * i  + 2 << " " << o_ENP << " " << \
						dN + NE1 + (NE1 - 1) * (NE2 - 2) << \
						" " << NE1 + NE2 << \
						" " << NE1 + NE2 + 1 << endl;
				}
				else {
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 1<< " " << o_ENP << " "\
						<< i + 1 + dN + (NE1 - 1) * (NE2 - 2)  \
						<< " " << 2 * NE1 + NE2 - i << " " << 2 * NE1 + NE2 - i + 1 << endl;
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 2 << " " << o_ENP << " "\
						<< i + 1 + dN + (NE1 - 1) * (NE2 - 2) << " " << i + 2 + dN + (NE1 - 1) * (NE2 - 2) \
						<< " " << 2 * NE1 + NE2 - i << endl;
				}
		}


		// теперь вывод узлов

		for (int i = 0; i < (NE1 + 1) * (NE2 + 1); i++)  // и все остальные в цикле
			f << i + 1 << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;

		f << dN + 1 << endl;
		for (int i = 0; i < dN + 1; i++)
			f << i + 1 << endl;
	}

	if (type == 3)
	{
		for (int i = 0; i < NE1; i++)
		{
			if (i == (NE1 - 1)) {
				f << 2 * NE1 - 1 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " "  << i + 1 + dN << endl;
				f << 2 * NE1 << " " << o_ENP << " " << i + 2 << " " << i + 3 << " " << i + 1 + dN << endl;
			}
			else {
				f << 2 * i + 1 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " " << i + 1 + dN << endl;
				f << 2 * i + 2 << " " << o_ENP << " " << i + 2 << " " << i + 2 + dN << " " << i + 1 + dN << endl;
			}
		}

		for (int j = 0; j < NE2 - 2; j++)
		{
			for (int i = 0; i < NE1; i++)
			{
				if (i == 0) {
					f << 2 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << o_ENP << " " << dN - j + 1 << " " << dN + 2 + (NE1 - 1) * j << " " \
						 << dN - j << endl;
					f << 2 * ((j + 1) * NE1 + (i + 1))  << " " << o_ENP << " "  << dN + 2 + (NE1 - 1) * j << " " \
						<< dN + 2 + (NE1 - 1) * (j + 1) << " " << dN - j << endl;
				}
				else
					if (i == (NE1 - 1)) {
						f << 2 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << o_ENP << " " << (NE1 - 1) * (j + 1) + dN + 1 << " " << NE1 + 2 + j << " " \
							 << (NE1 - 1) * (j + 2) + dN + 1 << endl;
						f << 2 * ((j + 1) * NE1 + (i + 1))  << " " << o_ENP <<  " " << NE1 + 2 + j << " " \
							<< NE1 + 3 + j << " " << (NE1 - 1) * (j + 2) + dN + 1 << endl;
					}
					else {
						f << 2 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << o_ENP << " " << i + 1 + dN + (NE1 - 1) * j << " " \
							<< i + 2 + dN + (NE1 - 1) * j << " " << \
							 i + 1 + dN + (NE1 - 1) * (j + 1) << endl;
						f << 2 * ((j + 1) * NE1 + (i + 1))  << " " << o_ENP << " "  \
							<< i + 2 + dN + (NE1 - 1) * j << " " << \
							i + 2 + dN + (NE1 - 1) * (j + 1) << " " << i + 1 + dN + (NE1 - 1) * (j + 1) << endl;
					}
			}
		}

		for (int i = 0; i < NE1; i++)
		{
			if (i == 0) {
				f << 2 * ((NE2 - 1) * NE1) + 2 * i + 1 << " " << o_ENP << " " << dN + 1 - NE2 + 2 << " " << (NE2 - 2) * (NE1 - 1) + 2 + dN << \
					 " " << dN + 1 - NE2 + 1 << endl;
				f << 2 * ((NE2 - 1) * NE1) + 2 * i + 2 << " " << o_ENP << " "  << (NE2 - 2) * (NE1 - 1) + 2 + dN << \
					" " << dN + 1 - NE2 << " " << dN + 1 - NE2 + 1 << endl;
			}
			else
				if (i == (NE1 - 1)) {
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 1 << " " << o_ENP << " " << dN + NE1 + (NE1 - 1) * (NE2 - 2) << \
						" " << NE1 + NE2 << \
						" "  << NE1 + NE2 + 2 << endl;
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 2 << " " << o_ENP  << \
						" " << NE1 + NE2 << \
						" " << NE1 + NE2 + 1 << " " << NE1 + NE2 + 2 << endl;
				}
				else {
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 1 << " " << o_ENP << " "\
						<< i + 1 + dN + (NE1 - 1) * (NE2 - 2) << " " << i + 2 + dN + (NE1 - 1) * (NE2 - 2) \
						<< " " << 2 * NE1 + NE2 - i + 1 << endl;
					f << 2 * ((NE2 - 1) * NE1) + 2 * i + 2 << " " << o_ENP << " "\
						 << i + 2 + dN + (NE1 - 1) * (NE2 - 2) \
						<< " " << 2 * NE1 + NE2 - i << " " << 2 * NE1 + NE2 - i + 1 << endl;
				}
		}


		// теперь вывод узлов

		for (int i = 0; i < (NE1 + 1) * (NE2 + 1); i++)  // и все остальные в цикле
			f << i + 1 << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;

		f << dN + 1 << endl;
		for (int i = 0; i < dN + 1; i++)
			f << i + 1 << endl;
	}

	if (type == 4)
	{
		for (int i = 0; i < NE1; i++)
		{
			if (i == (NE1 - 1))
			{
				f << 4 * NE1 - 3 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " " << xc << endl;
				f << 4 * NE1 - 2 << " " << o_ENP << " " << i + 2 << " " << i + 3 << xc << endl;
				f << 4 * NE1 - 1 << " " << o_ENP << " " << i + 3 << " " << i + 1 + dN << " " << xc << endl;
				f << 4 * NE1 << " " << o_ENP << " " << i + 1 + dN << " " << i + 1 << " "<< xc  << endl;
				xc++;
			}
			else
			{
				f << 4 * i + 1 << " " << o_ENP << " " << i + 1 << " " << i + 2 << " "  << xc << endl;
				f << 4 * i + 2 << " " << o_ENP << " " << i + 2 << " " << i + 2 + dN << " " << xc << endl;
				f << 4 * i + 3 << " " << o_ENP << " "  << i + 2 + dN << " " << i + 1 + dN << " " << xc <<endl;
				f << 4 * i + 4 << " " << o_ENP << " " << i + 1 + dN << " " << i + 1 << " " << xc << endl;
				xc++;
			}
				
		}

		for (int j = 0; j < NE2 - 2; j++)
		{
			for (int i = 0; i < NE1; i++)
			{
				if (i == 0)
				{
					f << 4 * ((j + 1) * NE1 + (i + 1)) - 3 << " " << o_ENP << " " << \
						dN - j + 1 << " " << dN + 2 + (NE1 - 1) * j << " " \
						<< xc << endl;
					f << 4 * ((j + 1) * NE1 + (i + 1)) - 2 << " " << o_ENP << " " << \
						dN + 2 + (NE1 - 1) * j << " " \
						<< dN + 2 + (NE1 - 1) * (j + 1) << " " << xc << endl;
					f << 4 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << o_ENP << " "  \
						<< dN + 2 + (NE1 - 1) * (j + 1) << " " << dN - j << " " << xc << endl;
					f << 4 * ((j + 1) * NE1 + (i + 1)) << " " << o_ENP << " " << \
						dN - j << "  " << dN - j + 1 << " " << xc << endl;
					xc++;
				}
					
				else
					if (i == (NE1 - 1))
					{
						f << 4 * ((j + 1) * NE1 + (i + 1)) - 3 << " " << o_ENP << " " << \
							(NE1 - 1)* (j + 1) + dN + 1 << " " << NE1 + 2 + j << " " \
							<< xc << endl;
						f << 4 * ((j + 1) * NE1 + (i + 1)) - 2 << " " << o_ENP << " " << \
							 NE1 + 2 + j << " " \
							<< NE1 + 3 + j << " " << xc << endl;
						f << 4 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << o_ENP << " " << \
							 NE1 + 3 + j << " " << (NE1 - 1) * (j + 2) + dN + 1 << " " << xc <<  endl;
						f << 4 * ((j + 1) * NE1 + (i + 1))  << " " << o_ENP << " " << \
							 (NE1 - 1) * (j + 2) + dN + 1 << " " << (NE1 - 1) * (j + 1) + dN + 1 << " " << xc << endl;
						xc++;
					}
					else
					{
						f << 4 * ((j + 1) * NE1 + (i + 1)) - 3 << " " << o_ENP << " " << \
							i + 1 + dN + (NE1 - 1) * j << " " \
							<< i + 2 + dN + (NE1 - 1) * j << " " << \
							xc << endl;
						f << 4 * ((j + 1) * NE1 + (i + 1)) - 2 << " " << o_ENP << " " << \
							i + 2 + dN + (NE1 - 1) * j << " " << \
							i + 2 + dN + (NE1 - 1) * (j + 1) << " " << \
							xc << endl;
						f << 4 * ((j + 1) * NE1 + (i + 1)) - 1 << " " << o_ENP << " " << \
							i + 2 + dN + (NE1 - 1) * (j + 1) << " " << \
							i + 1 + dN + (NE1 - 1) * (j + 1) << " " << xc << endl;
						f << 4 * ((j + 1) * NE1 + (i + 1)) << " " << o_ENP << " " << i + 1 + dN + (NE1 - 1) * (j + 1) << \
							" " << i + 1 + dN + (NE1 - 1) * j << " " << xc << endl;
						xc++;
					}
						
			}
		}

		for (int i = 0; i < NE1; i++)
		{
			if (i == 0)
			{
				f << 4 * ((NE2 - 1) * NE1) + 4 * i + 1 << " " << o_ENP << " " \
					<< dN + 1 - NE2 + 2 << " " \
					<< (NE2 - 2) * (NE1 - 1) + 2 + dN << " " \
					<< xc << endl;
				f << 4 * ((NE2 - 1) * NE1) + 4 * i + 2 << " " << o_ENP << " " \
					<< (NE2 - 2) * (NE1 - 1) + 2 + dN << " " \
					<< dN + 1 - NE2 << " " \
					<< xc << endl;
				f << 4 * ((NE2 - 1) * NE1) + 4 * i + 3 << " " << o_ENP << " " \
					<< dN + 1 - NE2 << " " \
					<< dN + 1 - NE2 + 1 << " " << xc << endl;
				f << 4 * ((NE2 - 1) * NE1) + 4 * i + 4 << " " << o_ENP << " " \
					<< dN + 1 - NE2 + 1 << " " \
					<< dN + 1 - NE2 + 2 << " " << xc <<endl;
				xc++;
			}
				
			else
				if (i == (NE1 - 1))
				{
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 1 << " " << o_ENP << " " \
						<< dN + NE1 + (NE1 - 1) * (NE2 - 2) << " " \
						<< NE1 + NE2 << " " \
						<< xc << endl;
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 2 << " " << o_ENP << " " \
						<< NE1 + NE2 << " " \
						<< NE1 + NE2 + 1 << " " \
						<< xc << endl;
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 3 << " " << o_ENP << " " \
						<< NE1 + NE2 + 1 << " " \
						<< NE1 + NE2 + 2 << " " \
						<< xc << endl;
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 4 << " " << o_ENP << " " \
						<< NE1 + NE2 + 2 <<  " " \
						<< dN + NE1 + (NE1 - 1) * (NE2 - 2) << " " \
						<< xc << endl;
					xc++;
				}
					
				else
				{
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 1 << " " << o_ENP << " "\
						<< i + 1 + dN + (NE1 - 1) * (NE2 - 2) << " " \
						<< i + 2 + dN + (NE1 - 1) * (NE2 - 2) << " " \
						<< xc << endl;
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 2 << " " << o_ENP << " "\
						<< i + 2 + dN + (NE1 - 1) * (NE2 - 2) << " " \
						<< 2 * NE1 + NE2 - i << " " \
						<< xc << endl;
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 3 << " " << o_ENP << " "\
						<< 2 * NE1 + NE2 - i << " " \
						<< 2 * NE1 + NE2 - i + 1 << " " \
						<< xc <<endl;
					f << 4 * ((NE2 - 1) * NE1) + 4 * i + 4 << " " << o_ENP << " "\
						<< 2 * NE1 + NE2 - i + 1 << " " \
						<< i + 1 + dN + (NE1 - 1) * (NE2 - 2) << " " << xc << endl;
					xc++;
				}
					
		}


		// теперь вывод узлов

		for (int i = 0; i < o_NP; i++)  // и все остальные в цикле
			f << i + 1 << " " << r[i * 3] << " " << r[i * 3 + 1] << " " << r[i * 3 + 2] << endl;

		f << dN + 1 << endl;
		for (int i = 0; i < dN + 1; i++)
			f << i + 1 << endl;
	}
	f.close();
}

int main()
{
	int NP, NE1, NE2, N, type;
	int o_NP, o_NE, o_NC, o_EN, o_ENP;
	vector<double> x;
	vector<double> y;   // исходные координаты точек конутра
	vector<double> z;
	vector<double> r;

	file_read("Np4_random.txt", N, NP, NE1, NE2, type, x, y, z);
	read_check(NP, NE1, NE2, type, x, y, z);

	if (NP == 2)
	{ 
		if (type == 1)
		{
			mesh2_1(NE1, type, x, y, z, r);
			output_2("output2.txt", NP, NE1, type, r);
		}
		else if (type == 2)
		{
			mesh2_2(NE1, type, x, y, z, r);
			output_2("output2.txt", NP, NE1, type, r);
		}
		else cout << "type variable is out of range, it must be 1 or 2" << endl;
	}

	if (NP == 4)
	{
		if (type == 1)
		{
			mesh4_123(NE1, NE2, type, x, y, z, r);
			output_4("output4.txt", NP, NE1, NE2, type, r);
		}
		else if (type == 2)
		{
			mesh4_123(NE1, NE2, type, x, y, z, r);
			output_4("output4.txt", NP, NE1, NE2, type, r);
		}
		else if (type == 3)
		{
			mesh4_123(NE1, NE2, type, x, y, z, r);
			output_4("output4.txt", NP, NE1, NE2, type, r);
		}
		else if(type == 4)
		{
			mesh4_4(NE1, NE2, type, x, y, z, r);
			output_4("output4.txt", NP, NE1, NE2, type, r);
		}
		else cout << "type variable is out of range, it must be in 1,2,3,4" << endl;
	}


//	cin.get();
	return 0;
}