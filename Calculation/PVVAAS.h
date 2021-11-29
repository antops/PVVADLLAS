#pragma once
/*
*	created by liyun 2018/1/23
	modified by Ming Jin 2019/03/07
*   function 计算反射面 输入镜子的指针 输出反射光线和交点以及是否相交
*   输入 源Ex_In,Ey_In, 源位置、传播方向，反射面
*   依赖RayTracing做光线追踪部分
*   version 1.0
*/
#ifndef PVVAAS_H
#define PVVAAS_H

#include <iostream>
#include <vector>
#include <cassert>
#include <vector>
#include <list>
#include <complex>

#include "../../Util/Vector3.h"
#include "../../Util/Vector2.h"
#include "../../Util/GEMS_Memory.h"
#include "../../Util/Constant_Var.h"
#include "../../Util/GraphTrans.h"
#include "../../Util/Fieldbase.h"
#include "../../FFTDIDLL/FFTDI.h"
#include "../../MirrorsAndRestriction/RayTracing.h"
#include "../../MirrorsAndRestriction/STLMirror.h"
#include "../../MirrorsAndRestriction/QuadricSurfaceMirror.h"
#include "../VTK/Field.h"

using namespace std;


class PVVAAS
{
public:
	PVVAAS();
	~PVVAAS();

	void Initialization(); // 数据初始化

	void setUnit(double factor); // 设置单位

	void setSource(const FieldBase* _in);

	void setOut(const FieldBase* _out);

	void setFre(double _fre);

	void setMirror(Mirror * mirror);

	void setOneSide(bool _in);

	void CalMirrorCenter();

	void CalZ0Theta();

	FieldBase getFieldBaseout();

	void CalArbitraryPlane(GraphTrans GT0, GraphTrans GTt, int mode = 0);

   // 得到源的中点坐标以及传播方向
	void getSourcePoint(Vector3 & interPoint, Vector3 & n_Point) const;

	// 源的中点坐标以及传播方向
	void setSourcePoint(const Vector3 & interPoint, const Vector3 & n_Point);

	static GraphTrans getSourceGraphTrans(const Vector3 & n_Point);

	bool CheckModle(Vector3 & interPoint, Vector3 & n_Point); // check是否能计算

	void ReflectCUDA();// 计算反射 改变了Exyz_R的相位和极化 更新了反射坐标

	void CalAmplitude(); // 计算反射前后的投影面积 更新Exyz_R幅度变换

	void InterVal(); // 对反射回来的波进行插值，使其落入标准网格中

	void Result(); // 得到结果面

	void getPlane(Vector3 ** &_org, Vector3 ** &_n) const; // 得到平面点和点的方向

	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// 注册回调函数

	void SetCATRMode(bool in){ CATRMode = in; }

	void SetThreadNum(int _threadNum) { threadNum = _threadNum; }

private:

	void Poynting(); //计算面的坡应廷矢量

	void AllocateMemory();  //分配内存

	void updateSource_n(const Vector3& new_n); // 更新源的法向量以及旋转信息

												   //复数的共轭相乘 A * B(共轭)
	complex<double> ConjugateMul(const complex<double> &A, const complex<double> &B) const;
	//复数的相乘 A * B
	complex<double> ComMul(const complex<double> &A, const complex<double> &B) const;

	double CalDistance(const Vector3 &a, const Vector3 &b) const;

	/************************************************
	* 直线与平面相交检测 orig直线起点 dir 射线方向
	* 返回t t>0 相交在直接前方 t<0 相交在直线后方
	* 平面为点法式 n*(P-P0) = 0
	* P0为平面上的一点 n为平面法线
	***************************************************/
	double IntersectPlane(const Vector3 &orig, const Vector3 &dir, const Vector3 &Plane_org,
		const Vector3 &Plane_n, Vector3 &intersection);

	double IntersectPoint(const Vector3 &orig, const Vector3 &dir, const Vector3 &v0,
		const Vector3 &E1, const Vector3 &E2, Vector3 &intersection);

	// 计算反射电场 输入法向量 入射电场ExEyEz 
	void CalReflectExyz(const Vector3 &n, const complex<double> &Ex, const complex<double> &Ey,
		const complex<double> &Ez, complex<double> &Ex_out, complex<double> &Ey_out, complex<double> &Ez_out);

	// 计算四个点围城的面积  等效于计算三角形ABC 和 三角形ADC的面积和
	double CalSquare(const Vector3 &A, const Vector3 &B, const Vector3 &C, const Vector3 &D) const;
	double CalSquare(const Vector3 &A, const Vector3 &B, const Vector3 &C) const;

	// 判断两线段是否有交点 有则求出交点
	bool get_line_intersection(const Vector2 &A, const Vector2 &B,
		const Vector2 &C, const Vector2 &D, Vector2 &O);

	// 插值用的——判断点是否在三角形内部
	bool InterVal_PointinTriangle(const Vector2 &A, const Vector2 &B,
		const Vector2 &C, const Vector2 &P);

private:
	int N, M; // 计算矩形的长宽 N*M
	double f; // 频率 默认 10.65e9
	double lamda; // 波长
				  //double k;
	double dx, dy;
	int threadNum;

	bool oneside;	//单次倾斜还是多次倾斜

	double unit_factor; // 保存模型的单位转化为米的系数

	Mirror * mirror;
	FieldBase * outFieldBasePtr;

	GraphTrans SourceGraphTrans;// 源的旋转参数
	GraphTrans OutFieldGraphTrans;//输出面的旋转参数

	complex<double> ** Px, ** Py, ** Pz;	// 坡应廷矢量 
	complex<double> ** Ex1, ** Ey1, ** Ez1;	// 经过calculation计算后的到的电磁场
	complex<double> ** Hx1, ** Hy1, ** Hz1;	// 

	complex<double> ** Ex_In, ** Ey_In;	// 源 

	complex<double> ** Ex_R, ** Ey_R, ** Ez_R;	// 经过反射计算后的到的电场

	//Org_Source平面中点 n_Source平面法向量 R_Source反射后的平面法向量
	Vector3 Org_Source, n_Source, R_Source; 

	Vector3 **Plane;    // 平面各点的坐标（绝对坐标）
	Vector3 **n_Plane;  // 平面传播各点后的法向量
	Vector3 **R_Plane;  // 反射后的各点的坐标
	Vector3 **Rn_Plane; // 反射后的各点的法向量

	bool ** Eno;


	double alfa;
	double beta;
	

	Vector3 n_Mirror;	//镜面的平均指向
	Vector3 c_Mirror;	//镜面的中心位置

	void(*returnFloat)(float, void*);
	void *user; // 回调函数的类指针
	bool CATRMode;

};




#endif
