/*
 * prog2.cpp
 * ��������� �����, ����� � �������� ���������� ���� ����������� ���������� �^- ��������������
 *  Created on: 3 ��� 2020 �.
 *      Author: 1
 */
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <vector>
#include "mx.hpp"
using namespace std;

const int maxLen2 = 200000;
double p[maxLen2][4], Layer1[maxLen2][4], Layer2[maxLen2][4];   // t,x,R,L

struct prn1 {
	int qp;    // ���������� �����, ������� ������� ����� � �������� ���������� ����
	double t;
	double x;
	double u;
	double c;
};
vector<prn1> netFinal;

void crossLayerTf(double l1[maxLen2][4], int count)
{
	double lambda, u0, u1, c0, c1;
	prn1 p;
	for (int i=0; i<count;i++)
		if ((l1[i][0]<=tf)and(l1[i+1][0]>=tf))
			{
			lambda = (tf-l1[i][0])/(l1[i+1][0]-l1[i][0]);
			p.t = l1[i][0] + lambda*(l1[i+1][0]-l1[i][0]);
			p.x = l1[i][1] + lambda*(l1[i+1][1]-l1[i][1]);
			u0 = (l1[i][2]+l1[i][3])/2;
			c0 = (l1[i][2]-l1[i][3])*(gamma-1)/4;
			u1 = (l1[i+1][2]+l1[i+1][3])/2;
			c1 = (l1[i+1][2]-l1[i+1][3])*(gamma-1)/4;
			p.u = u0 + lambda*(u1-u0);
			p.c = c0 + lambda*(c1-c0);
			if (p.u < 0.1) netFinal.push_back(p);  // ������� ������ ������������� �����, ����������� � (t_f, r_f)
			break;
			}
}

void L2ToL1_OutToFile(double l1[maxLen2][4], double l2[maxLen2][4], const int count)
{
	for (int i=0; i <= count; i++)
		{
		l1[i][0] = l2[i][0];
		l1[i][1] = l2[i][1];
		l1[i][2] = l2[i][2];
		l1[i][3] = l2[i][3];
//		double u = (l1[i][2]+l1[i][3])/2;
//		double c = (l1[i][2]-l1[i][3])*(gamma-1)/4;
//		output << numOfLayer << ";" << i << ";" << l1[i][0] <<  ";" << Layer1[i][1] << ";" << u << ";" <<
//				c << ";" << u+c << ";" << u-c << ";" << endl;
		}
}

//  ����� �������� ��������� ======================================================================
int main_prog2(){
//	double p[maxLen2][4], Layer1[maxLen2][4], Layer2[maxLen2][4];   // t,x,R,L
	cout << fixed << setprecision(16);
	readPiston("piston5.csv", p); // ��������� ���������� ������

	double timeEndRegular, rEndRegular;
	int nr = 25000;
	double dr = (rw-rs)/nr;
	int i;

	// �������������� ������ ����
	int numPointOnLayer = 0;
	int numOfLayer = 0;
	int maxNumPointOnLayer1;
	Layer1[numPointOnLayer][0] = ts;
	Layer1[numPointOnLayer][1] = rs;
	Layer1[numPointOnLayer][2] = 2*c0/(gamma-1);
	Layer1[numPointOnLayer][3] = -2*c0/(gamma-1);
	cout << "��������� ���� ����� = " << numOfLayer << endl;
// ������ ��������� ���� �� ����� r=rw
	cout << "������ ������� 1: ����� ���� 1, ����� �� ������ = " << numOnPiston << endl;
	while (Layer1[0][1] < rw)
		{
		numOfLayer ++;
		Layer2[0][0] = ts;
		Layer2[0][1] = rs + numOfLayer*dr;
		Layer2[0][2] = 2*c0/(gamma-1);
		Layer2[0][3] = -2*c0/(gamma-1);
		for (i=1; i < numOfLayer*2; i++)
			{
			newPoint1(Layer2[i-1], Layer1[i-1], Layer2[i]);
			}
		while (newPoint3_prog2(p[numOnPiston], p[numOnPiston-1], Layer2[numOfLayer*2-1], Layer2[numOfLayer*2]) != 1)
			{
			numOnPiston -=1;
			}

		L2ToL1_OutToFile(Layer1, Layer2, numOfLayer*2);
		}
		cout << "��������� ������� 1: ����� ���� = " << numOfLayer << " ����� �� ������ = " << numOnPiston << endl;

	maxNumPointOnLayer1 = numOfLayer*2;
// ������ ��������� ���� ����� r=rw ���� ����� �������������� ��� ���������� ���������� ������
	cout << "������ ������� 2: ����� ���� = " << numOfLayer << " ����� �� ������ = " << numOnPiston << endl;
	while (numOnPiston > 1)
		{
		newPoint2(Layer1[1], Layer2[0]);
		for (i=2; i <= maxNumPointOnLayer1; i++)
			{
			newPoint1(Layer2[i-2], Layer1[i], Layer2[i-1]);
			}
		while (newPoint3_prog2(p[numOnPiston], p[numOnPiston-1], Layer2[maxNumPointOnLayer1-1], Layer2[maxNumPointOnLayer1]) != 1)
			{
			numOnPiston -=1;
			if (numOnPiston == 0)
				{
				timeEndRegular = Layer2[maxNumPointOnLayer1-1][0];
				rEndRegular = Layer2[maxNumPointOnLayer1-1][1];
				break;
				}
			}
		L2ToL1_OutToFile(Layer1, Layer2, maxNumPointOnLayer1);
		numOfLayer ++;
//		cout << "������� 2: ��������� ���� ����� = " << numOfLayer <<
//			", ����� �� ������ = " << numOnPiston << endl;
		}
	cout << "��������� ������� 2: ����� ���� = " << numOfLayer << " ����� �� ������ = " << numOnPiston << endl;
// ������ ��������� ���� ����� ����� r=rw ����� ���������� ������ �����������
	cout << "������ ������� 3: ����� ���� = " << numOfLayer << " ����� ����� �� ���� = " << maxNumPointOnLayer1 << endl;
	while (maxNumPointOnLayer1 > 1)
		{
		newPoint2(Layer1[1], Layer2[0]);
		for (i=2; i <= maxNumPointOnLayer1; i++)
			{
			newPoint1(Layer2[i-2], Layer1[i], Layer2[i-1]);
			}
		maxNumPointOnLayer1 -= 1;
		if ((Layer2[0][0]<=tf)and(Layer2[maxNumPointOnLayer1][0]>=tf))
			{
//			cout << "���� ����� = " << numOfLayer << ", ��������� ����� t == tf" << endl;
			crossLayerTf(Layer2, maxNumPointOnLayer1);
			}
		else{
			// �� ��������� t=tf, ��������� ��� ����������� �^- �������������� �������� � ����� (tf,rf)
			double k = (Layer2[maxNumPointOnLayer1][2]+Layer2[maxNumPointOnLayer1][3])/2;
			k = k - (gamma-1)*(Layer2[maxNumPointOnLayer1][2]-Layer2[maxNumPointOnLayer1][3])/4;
//			cout << "�������������� �������� � ����� " << Layer2[maxNumPointOnLayer1][1] + k*(tf-Layer2[maxNumPointOnLayer1][0]) << endl;
			}

		L2ToL1_OutToFile(Layer1, Layer2, maxNumPointOnLayer1);
		numOfLayer ++;
//		cout << "������� 3: ��������� ���� ����� = " << numOfLayer << ", ����� ����� �� ���� = " << maxNumPointOnLayer1 << endl;
		}
	cout << "��������� ������� 3: ����� ���� = " << numOfLayer << " ����� ����� �� ���� = " << maxNumPointOnLayer1 << endl;
	cout << "������ �������� " << endl;
//	output.close();

// ================================================================================================
	ofstream outNetFinal;
	outNetFinal.open("netFinal.csv");
	outNetFinal << fixed << setprecision(16);
	outNetFinal << "t;x;u;c;rho" << endl;
	for (auto pnt:netFinal)
		{
		cout << "t = " << pnt.t << ", x = " << pnt.x << ", u = " << pnt.u << ", c = " << pnt.c << ", rho = " << pow(pnt.c, 2/(gamma-1)) << endl;
		outNetFinal << pnt.t << ";" << pnt.x << ";" << pnt.u << ";" << pnt.c << ";" << pow(pnt.c, 2/(gamma-1)) << endl;
		}
	outNetFinal.close();
// ================================================================================================
	ofstream tabl1;
	tabl1.open("tabl1.csv");
	tabl1 << fixed << setprecision(16);
	tabl1 << "r_s;t_f;r_f;maxNum;u_max;c_max;rho_max" << endl;
	tabl1 << rs << ";" << tf << ";" << rf << ";" << maxNum << ";" << u_max << ";" << c_max << ";" << rho_max << endl;
	tabl1.close();
// ================================================================================================
	ofstream tabl2;
	tabl2.open("tabl2.csv");
	tabl2 << fixed << setprecision(160);
	tabl2 << "gamma;nu;rho_;m;n;����������C����" << endl;
	tabl2 << gamma << ";" << nu << ";" << rho_ << ";" << m << ";" << nr << ";" << numOfLayer << endl;
	tabl2.close();
// ================================================================================================
	ofstream tabl3;
	tabl3.open("tabl3.csv");
	tabl3 << fixed << setprecision(16);
	tabl3 << "dm;t_E1;r_E1;r_E2;T%;d_reg%" << endl;
	tabl3 << dm(rs, rho0, netFinal[0].x, pow((netFinal[0].c+netFinal[netFinal.size()-1].c)/2, 2/(gamma-1)), rw)
			<< ";" << timeEndRegular << ";" << rEndRegular << ";" << netFinal[0].x << ";" << 100*timeEndRegular/tf
			<< ";" << 100*(rw-netFinal[0].x)/(rw-rf) << endl;
	tabl3.close();
// ================================================================================================
	cout << "=====================================================================================" << endl;
	cout << "����������� ���� �����: " << dm(rs, rho0, netFinal[0].x,
			pow((netFinal[0].c+netFinal[netFinal.size()-1].c)/2, 2/(gamma-1)), rw) << "%" <<endl;
	cout << "��������� ������ ������� ����� ����� ���� ���������� (����� �1): " << timeEndRegular << endl;
	cout << "��������� �������� ���������� r ����� ����� ���� ���������� (����� �1): " << rEndRegular << endl;
	cout << "�������� ���������� r ����� ����� ����� ����� ���������� (����� �� ����� t = tf, ����� �2): " << netFinal[0].x << endl;
	cout << 100*timeEndRegular/tf << " % ������� ������ ����� ��������� ��������� (T_reg)" << endl;
//	cout << "���������� �� ������� ��� ����� �� ��������: " << netFinal[0].x-rEndRegular << endl;
//	cout << 100*(netFinal[0].x-rEndRegular)/(rw-rf) << " %, ����������� ���������� ���������� ������������ � ���������� ����� ������� ������ � ����������� �������" << endl;
//	cout << 100*(netFinal[0].x-rEndRegular)/(rw-rs) << " %, ����������� ���������� ���������� ������������ � ���������� ����� ������� ������ � ����������� �������" << endl;
	cout << 100*(rw-netFinal[0].x)/(rw-rf) << " % ���������� � ������ ������ ��� ��������� ���������� �����, d_reg" << endl;

	cout << "ts = " << ts << endl;
	cout << "rs = " << rs << endl;
	cout << "tf = " << tf << endl;
	cout << "rf = " << rf << endl;
	cout << "rw = " << rw << endl;
	cout << "c_ = " << pow(rho_,(gamma-1)/2) << endl;

	cout << "==  ������ �������  =================================================================" << endl;
	return 0;
}


