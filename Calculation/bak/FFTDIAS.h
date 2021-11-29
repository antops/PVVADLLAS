
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
#ifndef FFTDIAS_H
#define FFTDIAS_H

#include <cmath>
#include <complex>
#include "../../Util/GEMS_Memory.h"
#include "../../Util/Constant_Var.h"
#include "../../Util/Vector3.h"
#include "FFT.h"

class FFTDIAS
{
public:
	FFTDIAS(double f = 10.65e9, double z0 = 1,int N = 361, int M = 361);
	~FFTDIAS();

	void Initialization(); 
	void FreeCal();
	void Setds(double ds1); //����ds
	void SetPropagationParas(double _xp, double _yp, double _zp, double _theta_rou, double _phi_rou);

	//�������벢��0
	void SetInput(std::complex<double> ** EuIn, std::complex<double> ** EvIn);

	void StartCal();

	//���Ex Ey Ez Hx Hy Hz ������FreeCal
	void output(std::complex<double> ** Eu, std::complex<double> ** Ev, std::complex<double> ** En, 
		std::complex<double> ** Hu, std::complex<double> ** Hv, std::complex<double> ** Hn );

private:
	//���帴���ĳ�
	void MulComp(const double r1, const double i1, const double r2, const double i2, double & outR, double & outI);
	std::complex<double> MulComp(std::complex<double> & In, double r2, double i2);
	//��ֵ����
	double InserVal(const double x0, const double yUp, const double y0, const double yDown);

	//���㴫�ﺯ��
	void Calh0();

	//�����볡������ת
	void axisRotate();//������������䳡�ֲ�����������ת����б����

	//�����볡������ת
	void Rotate();

	//�����볡���в�ֵ
	//void InserExEy();

private:
	double f; // Ƶ�� Ĭ�� 10.65e9
	double lamda; // ����
	double k;
	//double theta;
	//double theta_h0;
	double ds;
	double z0; //Ŀ�����

	double xp, yp, zp;//��������
	double theta_rou, phi_rou;//��ת


	int N, M; //����ĵ��� һ����Ϊ���� Ĭ��360
	int N2, M2; // N2 = 2 * N -1

	//���ݺ���
	std::complex<double> ** HExy, ** HEz_x0, **HEz_y0;
	std::complex<double> ** HHx_x0, ** HHx_y0,
		** HHy_x0, ** HHy_y0,
		** HHz_x0, ** HHz_y0;

	//��0�������
	std::complex<double> ** Ex0, ** Ey0;

	// ��˹ά������
	double 	PowerHExy, PowerFFTHExy;
	double 	PowerHEz_x0, PowerFFTHEz_x0;
	double 	PowerHEz_y0, PowerFFTHEz_y0;
	double 	PowerHHx_x0, PowerFFTHHx_x0;
	double 	PowerHHx_y0, PowerFFTHHx_y0;
	double 	PowerHHy_x0, PowerFFTHHy_x0;
	double 	PowerHHy_y0, PowerFFTHHy_y0;
	double 	PowerHHz_x0, PowerFFTHHz_x0;
	double 	PowerHHz_y0, PowerFFTHHz_y0;

	//Ex1 Ey1 Ez1 Hx1 Hy1 Hz1
	std::complex<double> ** Ex1, ** Ey1, **Ez1;
	std::complex<double> ** Hx1, ** Hy1, **Hz1;
};



#endif // !CALCUlATION_H
