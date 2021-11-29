#pragma once
/*
*	created by liyun 2018/1/23
	modified by Ming Jin 2019/03/07
*   function ���㷴���� ���뾵�ӵ�ָ�� ���������ߺͽ����Լ��Ƿ��ཻ
*   ���� ԴEx_In,Ey_In, Դλ�á��������򣬷�����
*   ����RayTracing������׷�ٲ���
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

	void Initialization(); // ���ݳ�ʼ��

	void setUnit(double factor); // ���õ�λ

	void setSource(const FieldBase* _in);

	void setOut(const FieldBase* _out);

	void setFre(double _fre);

	void setMirror(Mirror * mirror);

	void setOneSide(bool _in);

	void CalMirrorCenter();

	void CalZ0Theta();

	FieldBase getFieldBaseout();

	void CalArbitraryPlane(GraphTrans GT0, GraphTrans GTt, int mode = 0);

   // �õ�Դ���е������Լ���������
	void getSourcePoint(Vector3 & interPoint, Vector3 & n_Point) const;

	// Դ���е������Լ���������
	void setSourcePoint(const Vector3 & interPoint, const Vector3 & n_Point);

	static GraphTrans getSourceGraphTrans(const Vector3 & n_Point);

	bool CheckModle(Vector3 & interPoint, Vector3 & n_Point); // check�Ƿ��ܼ���

	void ReflectCUDA();// ���㷴�� �ı���Exyz_R����λ�ͼ��� �����˷�������

	void CalAmplitude(); // ���㷴��ǰ���ͶӰ��� ����Exyz_R���ȱ任

	void InterVal(); // �Է�������Ĳ����в�ֵ��ʹ�������׼������

	void Result(); // �õ������

	void getPlane(Vector3 ** &_org, Vector3 ** &_n) const; // �õ�ƽ���͵�ķ���

	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// ע��ص�����

	void SetCATRMode(bool in){ CATRMode = in; }

	void SetThreadNum(int _threadNum) { threadNum = _threadNum; }

private:

	void Poynting(); //���������Ӧ͢ʸ��

	void AllocateMemory();  //�����ڴ�

	void updateSource_n(const Vector3& new_n); // ����Դ�ķ������Լ���ת��Ϣ

												   //�����Ĺ������ A * B(����)
	complex<double> ConjugateMul(const complex<double> &A, const complex<double> &B) const;
	//��������� A * B
	complex<double> ComMul(const complex<double> &A, const complex<double> &B) const;

	double CalDistance(const Vector3 &a, const Vector3 &b) const;

	/************************************************
	* ֱ����ƽ���ཻ��� origֱ����� dir ���߷���
	* ����t t>0 �ཻ��ֱ��ǰ�� t<0 �ཻ��ֱ�ߺ�
	* ƽ��Ϊ�㷨ʽ n*(P-P0) = 0
	* P0Ϊƽ���ϵ�һ�� nΪƽ�淨��
	***************************************************/
	double IntersectPlane(const Vector3 &orig, const Vector3 &dir, const Vector3 &Plane_org,
		const Vector3 &Plane_n, Vector3 &intersection);

	double IntersectPoint(const Vector3 &orig, const Vector3 &dir, const Vector3 &v0,
		const Vector3 &E1, const Vector3 &E2, Vector3 &intersection);

	// ���㷴��糡 ���뷨���� ����糡ExEyEz 
	void CalReflectExyz(const Vector3 &n, const complex<double> &Ex, const complex<double> &Ey,
		const complex<double> &Ez, complex<double> &Ex_out, complex<double> &Ey_out, complex<double> &Ez_out);

	// �����ĸ���Χ�ǵ����  ��Ч�ڼ���������ABC �� ������ADC�������
	double CalSquare(const Vector3 &A, const Vector3 &B, const Vector3 &C, const Vector3 &D) const;
	double CalSquare(const Vector3 &A, const Vector3 &B, const Vector3 &C) const;

	// �ж����߶��Ƿ��н��� �����������
	bool get_line_intersection(const Vector2 &A, const Vector2 &B,
		const Vector2 &C, const Vector2 &D, Vector2 &O);

	// ��ֵ�õġ����жϵ��Ƿ����������ڲ�
	bool InterVal_PointinTriangle(const Vector2 &A, const Vector2 &B,
		const Vector2 &C, const Vector2 &P);

private:
	int N, M; // ������εĳ��� N*M
	double f; // Ƶ�� Ĭ�� 10.65e9
	double lamda; // ����
				  //double k;
	double dx, dy;
	int threadNum;

	bool oneside;	//������б���Ƕ����б

	double unit_factor; // ����ģ�͵ĵ�λת��Ϊ�׵�ϵ��

	Mirror * mirror;
	FieldBase * outFieldBasePtr;

	GraphTrans SourceGraphTrans;// Դ����ת����
	GraphTrans OutFieldGraphTrans;//��������ת����

	complex<double> ** Px, ** Py, ** Pz;	// ��Ӧ͢ʸ�� 
	complex<double> ** Ex1, ** Ey1, ** Ez1;	// ����calculation�����ĵ��ĵ�ų�
	complex<double> ** Hx1, ** Hy1, ** Hz1;	// 

	complex<double> ** Ex_In, ** Ey_In;	// Դ 

	complex<double> ** Ex_R, ** Ey_R, ** Ez_R;	// ������������ĵ��ĵ糡

	//Org_Sourceƽ���е� n_Sourceƽ�淨���� R_Source������ƽ�淨����
	Vector3 Org_Source, n_Source, R_Source; 

	Vector3 **Plane;    // ƽ���������꣨�������꣩
	Vector3 **n_Plane;  // ƽ�洫�������ķ�����
	Vector3 **R_Plane;  // �����ĸ��������
	Vector3 **Rn_Plane; // �����ĸ���ķ�����

	bool ** Eno;


	double alfa;
	double beta;
	

	Vector3 n_Mirror;	//�����ƽ��ָ��
	Vector3 c_Mirror;	//���������λ��

	void(*returnFloat)(float, void*);
	void *user; // �ص���������ָ��
	bool CATRMode;

};




#endif
