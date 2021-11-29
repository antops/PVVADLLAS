/**************************************************************************************
* version 1.0 2019/2/26
* Designer ����
* E-mail jinmingaps@163.com
* Fuction��
* ������ת��
* ���� Eu0, Ev0���õ�������ת���Eut,Evt;
* ��Ҫ����Ĳ�����Ƶ�ʣ�Ŀ����룬�������N*M;
Ŀ������б�ǣ�theta��phi����ds��Ĭ�ϲ���/3.5��
***************************************************************************************/
#pragma once
#ifndef FFTASRotation_H
#define FFTASRotation_H

#include <cmath>
#include <complex>
#include <vector>
#include "../../Util/GEMS_Memory.h"
#include "../../Util/Constant_Var.h"
#include "../../Util/Vector3.h"
#include "../../Util/Matrix4D.h"
#include "../../Util/Vector2.h"
#include "FFT.h"

class FFTASRotation
{
public:
	FFTASRotation(double _f = 10.65e9, int _Nu = 361, int _Nv = 361);
	~FFTASRotation();
	void SetParas(double _f, double _Nu, double _Nv, double _ds);
	void SetRotationParasinRad(double _theta_rou, double _phi_rou);
	void Allocate();
	void FreeCal();
	//���Ex Ey Ez������FreeCal
	void output(std::complex<double> ** EuOut, std::complex<double> ** EvOut);
	//�������벢��0
	void SetInput(std::complex<double> ** EuIn, std::complex<double> ** EvIn);
	void PerformRotate();
	//���Ex Ey Ez Hx Hy Hz ������FreeCal

private:

	bool InterVal_PointinTriangle(const Vector2 & A,
		const Vector2 & B, const Vector2 & C, const Vector2 & P);

private:
	double freq; // Ƶ�� Ĭ�� 10.65e9
	double lambda; // ����
	double k;
	double ds;
	double theta_rou, phi_rou;//��ת
	int Nu, Nv; //����ĵ��� һ����Ϊ���� Ĭ��361
	int Nu2, Nv2; // N2 = 2 * N -1

	std::complex<double> ** Eu0, ** Ev0, ** En0;
	std::complex<double> ** Eut, ** Evt, ** Ent;
	Matrix4D RotationMatrix;
	Matrix4D InvRotationMatrix;
};



#endif // !CALCUlATION_H
