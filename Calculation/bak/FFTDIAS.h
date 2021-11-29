
/**************************************************************************************
* version 1.0 2019/2/26
* Designer 金铭
* E-mail jinmingaps@163.com
* Fuction： 
* 坐标旋转：
* 输入 Eu0, Ev0，得到坐标旋转后的Eut,Evt;
* 需要输入的参数，频率，目标距离，计算点数N*M;
  目标面倾斜角（theta和phi）、ds（默认波长/3.5）
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
	void Setds(double ds1); //设置ds
	void SetPropagationParas(double _xp, double _yp, double _zp, double _theta_rou, double _phi_rou);

	//设置输入并补0
	void SetInput(std::complex<double> ** EuIn, std::complex<double> ** EvIn);

	void StartCal();

	//输出Ex Ey Ez Hx Hy Hz 并调用FreeCal
	void output(std::complex<double> ** Eu, std::complex<double> ** Ev, std::complex<double> ** En, 
		std::complex<double> ** Hu, std::complex<double> ** Hv, std::complex<double> ** Hn );

private:
	//定义复数的乘
	void MulComp(const double r1, const double i1, const double r2, const double i2, double & outR, double & outI);
	std::complex<double> MulComp(std::complex<double> & In, double r2, double i2);
	//插值函数
	double InserVal(const double x0, const double yUp, const double y0, const double yDown);

	//计算传达函数
	void Calh0();

	//对输入场进行旋转
	void axisRotate();//这个函数将入射场分布进行坐标旋转至倾斜口面

	//对输入场进行旋转
	void Rotate();

	//对输入场进行插值
	//void InserExEy();

private:
	double f; // 频率 默认 10.65e9
	double lamda; // 波长
	double k;
	//double theta;
	//double theta_h0;
	double ds;
	double z0; //目标距离

	double xp, yp, zp;//传播距离
	double theta_rou, phi_rou;//旋转


	int N, M; //计算的点数 一般设为奇数 默认360
	int N2, M2; // N2 = 2 * N -1

	//传递函数
	std::complex<double> ** HExy, ** HEz_x0, **HEz_y0;
	std::complex<double> ** HHx_x0, ** HHx_y0,
		** HHy_x0, ** HHy_y0,
		** HHz_x0, ** HHz_y0;

	//补0后的输入
	std::complex<double> ** Ex0, ** Ey0;

	// 帕斯维尔定理
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
