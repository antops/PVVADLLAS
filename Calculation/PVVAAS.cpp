#include "PVVAAS.h"
#include "../../Util/Position3D.h"
#include "../../Util/Matrix4D.h"
#include "../../Util/Vector3D.h"
#include "../../Util/Rays.h"
#include "../../FFTASRotationDLL/FFTASRotation.h"
#include "../../FFTDIDLL/FFTDI.h"
#include <cmath>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkSTLWriter.h>


using namespace calculation;


PVVAAS::PVVAAS()
{
	Initialization();
}

PVVAAS::~PVVAAS()
{
	Free_2D(Px);
	Free_2D(Py);
	Free_2D(Pz);

	Free_2D(Plane);
	Free_2D(n_Plane);
	Free_2D(R_Plane);
	
	Free_2D(Ex1);
	Free_2D(Ey1);
	Free_2D(Ez1);
	Free_2D(Hx1);
	Free_2D(Hy1);
	Free_2D(Hz1);

	Free_2D(Ex_R);
	Free_2D(Ey_R);
	Free_2D(Ez_R);
	Free_2D(Eno);


}

void PVVAAS::Initialization()
{
	N = 0; M = 0;
	f = 10.65e9;
	lamda = 0; 
	//k = 0;
	dx = 0;	dy = 0;
	unit_factor = 1;

	Ex_In = NULL;
	Ey_In = NULL;

	Px = NULL;
	Py = NULL;
	Pz = NULL;

	Plane = NULL;
	n_Plane = NULL;
	R_Plane = NULL;
	Rn_Plane = NULL;
	Eno = NULL;

	Ex1 = NULL;
	Ey1 = NULL;
	Ez1 = NULL;
	Hx1 = NULL;;
	Hy1 = NULL;
	Hz1 = NULL;

	Ex_R = NULL;
	Ey_R = NULL;
	Ez_R = NULL;

	returnFloat = NULL;
	user = NULL;

	threadNum = 1;
	CATRMode = false;	oneside = false;
}

void PVVAAS::setUnit(double factor)
{
	unit_factor = factor;
}

void PVVAAS::setFre(double _fre)
{
	f = _fre;
	lamda = C_Speed / f;

}

void PVVAAS::setSource(const FieldBase* _in)
{
	const vector<vector<complex<double>>>&Ex_temp = _in->Ex;
	const vector<vector<complex<double>>>&Ey_temp = _in->Ey;

	SourceGraphTrans = _in->graphTransField;
	SourceGraphTrans.normalization(unit_factor);
	dx = _in->ds_x*unit_factor;	dy = _in->ds_y*unit_factor;
	N = _in->N_width;			M = _in->M_depth;

	AllocateMemory();

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Ex_In[i][j] = Ex_temp[i][j];
			Ey_In[i][j] = Ey_temp[i][j];
		}

	Vector3D RotateAsix(SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrix = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D translateMatrix = Matrix4D::getTranslateMatrix(SourceGraphTrans.getTrans_x(),
		SourceGraphTrans.getTrans_y(), SourceGraphTrans.getTrans_z());


	// ����Դ�������λ��
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Position3D tempPoint((i - (N - 1) / 2) * dx, (j - (M - 1) / 2) * dy, 0);
			tempPoint = translateMatrix * rotatMatrix * tempPoint;
			Plane[i][j].set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
		}

	// ����Դ�ķ���
	Position3D tempPoint(0, 0, 1);
	tempPoint = rotatMatrix * tempPoint;
	n_Source.set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
	n_Source.Normalization();

	// ƽ����е�
	Position3D tempPoint1(0, 0, 0);
	tempPoint1 = translateMatrix * rotatMatrix * tempPoint1;
	Org_Source.set(tempPoint1.X(), tempPoint1.Y(), tempPoint1.Z());

	//	updateSource_n(Vector3(1, 1, 0));
}

void PVVAAS::setMirror(Mirror * mirror)
{
	this->mirror = mirror;
}

void PVVAAS::setOut(const FieldBase* _out) {
	OutFieldGraphTrans = _out->graphTransField;
}

void PVVAAS::setOneSide(bool in) {
	oneside = in;
}

//���㾵�����ģ�ȷ��������1��λ�ã�����ɴ�������
void PVVAAS::CalMirrorCenter()//����STL����������
{
	//�������ӵ������ҳ����ӵ�����λ��
	vtkSmartPointer<vtkPolyData> polyData = mirror->getPolyData();
	int EleNum = polyData->GetNumberOfCells();

	//��ʼ��
	n_Mirror = Vector3(0.0, 0.0, 0.0);
	c_Mirror = Vector3(0.0, 0.0, 0.0);
	GraphTrans RotProGraphTrans;

	double t;
	for (int i = 0; i < EleNum; i++)  //���뷴����Ľ���
	{
		vtkIdList * p;
		p = polyData->GetCell(i)->GetPointIds();
		double * point;
		point = polyData->GetPoint(p->GetId(0));
		Vector3 NodesXYZ1(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(1));
		Vector3 NodesXYZ2(point[0], point[1], point[2]);
		point = polyData->GetPoint(p->GetId(2));
		Vector3 NodesXYZ3(point[0], point[1], point[2]);

		Vector3 tempa = NodesXYZ1 - NodesXYZ2;
		Vector3 tempb = NodesXYZ1 - NodesXYZ3;
		Vector3 n_light = tempa.Cross(tempb);  //������

		n_Mirror = n_Mirror + n_light;
		c_Mirror = c_Mirror + (NodesXYZ1 + NodesXYZ2 + NodesXYZ3);
	}
	//���ӵ����ķ���
	n_Mirror.Normalization();
	Vector3 temp = Vector3(c_Mirror.x / EleNum / 3.0, c_Mirror.y / EleNum / 3.0, c_Mirror.z / EleNum / 3.0);
	//���ӵ�����
	c_Mirror = temp;
	c_Mirror = c_Mirror;// +n_Mirror*(-0.08);
	//ע�⴫�����Ǿ������ķ���ı���Ŷ
	GraphTrans::updateSource_n(Vector3(-n_Mirror.x,-n_Mirror.y,-n_Mirror.z), RotProGraphTrans);
	RotProGraphTrans.updateTranslate(c_Mirror);

	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(20, user);
	}
	//�����䳡������������1�� ��Ex_in ���� Ex1
	CalArbitraryPlane(SourceGraphTrans, RotProGraphTrans);


	//�������ˣ��������䳡�Ѿ����㵽��������1����
	//����Դ�������������Դ��������������1����
	Org_Source = c_Mirror;
	n_Source = Vector3(-n_Mirror.x, -n_Mirror.y, -n_Mirror.z);
	SourceGraphTrans = RotProGraphTrans;
	
	//����������1�ϸ������λ����Ϣ
	Vector3D RotateAsix(SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrix = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D translateMatrix = Matrix4D::getTranslateMatrix(SourceGraphTrans.getTrans_x(),
		SourceGraphTrans.getTrans_y(), SourceGraphTrans.getTrans_z());
	// ����Դ�������λ��-�ֲ�����ϵ��ȫ������ϵ
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Position3D tempPoint((i - (N - 1) / 2) * dx, (j - (M - 1) / 2) * dy, 0);
			tempPoint = translateMatrix * rotatMatrix * tempPoint;
			Plane[i][j].set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
		}

	//���ڸ���������2��ָ����Ϣ
	if (oneside) {//������2�������ƽ��
		R_Source = Matrix4D::getRotateMatrix(OutFieldGraphTrans.getRotate_theta(), OutFieldGraphTrans.getRotate_x(), OutFieldGraphTrans.getRotate_y(), OutFieldGraphTrans.getRotate_z())
			     * Vector3(0, 0, 1.0);
	}
	else {//������ָ���뾵��ƽ��ָ����ͬ
		R_Source = n_Mirror;
	}

	//����������1�����䳡��poyntingʸ��
	Poynting();
	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(30, user);
	}

}

void PVVAAS::CalZ0Theta()//���ڽ����� ���ڽ�����
{
	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(10, user);
	}

	Vector3 tempReflect;
	double t;

	int i = 0;
	RayTracing rayTracing(mirror);
	bool tempIsIntersect = false;


	rayTracing.calcNormalOfLine_Mirror(Org_Source, n_Source, n_Mirror, c_Mirror, tempIsIntersect, t);
	n_Mirror.Normalization();
	if (!tempIsIntersect)
		return; // û���ཻ

				//����ƽ�Ƶľ��루��ʱ��������û�иı䣩
	GraphTrans GTsource = SourceGraphTrans;
	GraphTrans RotProGraphTrans;
	RotProGraphTrans.updateTranslate(c_Mirror);
	GraphTrans::updateSource_n(n_Mirror*(-1.0), RotProGraphTrans);

	CalArbitraryPlane(GTsource, RotProGraphTrans);

	//�������ˣ��������䳡�Ѿ����㵽��������1����
	//����Դ�������������Դ��������������1����
	Org_Source = c_Mirror;
	n_Source = n_Mirror*(-1.0);
	SourceGraphTrans = RotProGraphTrans;

	//����������1�ϸ������λ����Ϣ
	Vector3D RotateAsix(SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrix = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D translateMatrix = Matrix4D::getTranslateMatrix(SourceGraphTrans.getTrans_x(),
		SourceGraphTrans.getTrans_y(), SourceGraphTrans.getTrans_z());
	// ����Դ�������λ��-�ֲ�����ϵ��ȫ������ϵ
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Position3D tempPoint((i - (N - 1) / 2) * dx, (j - (M - 1) / 2) * dy, 0);
			tempPoint = translateMatrix * rotatMatrix * tempPoint;
			Plane[i][j].set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
		}

	FILE* writefile;
	string fileAddress = "./PVVD/outasbRef.dat";
	errno_t err;
	err = fopen_s(&writefile, fileAddress.c_str(), "wb");
	if (err != 0) return;
	double tx, ty, tz, rx, ry, rz, rth;
	int flag = 0;
	tx = SourceGraphTrans.getTrans_x(); ty = SourceGraphTrans.getTrans_y(); tz = SourceGraphTrans.getTrans_z();
	rx = SourceGraphTrans.getRotate_x(); ry = SourceGraphTrans.getRotate_y(); rz = SourceGraphTrans.getRotate_z();
	rth = SourceGraphTrans.getRotate_theta();
	fwrite(&f, sizeof(double), 1, writefile);
	fwrite(&tx, sizeof(double), 1, writefile); fwrite(&ty, sizeof(double), 1, writefile); fwrite(&tz, sizeof(double), 1, writefile);
	fwrite(&rx, sizeof(double), 1, writefile); fwrite(&ry, sizeof(double), 1, writefile); fwrite(&rz, sizeof(double), 1, writefile);
	fwrite(&rth, sizeof(double), 1, writefile);
	fwrite(&N, sizeof(int), 1, writefile);
	fwrite(&M, sizeof(int), 1, writefile);
	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);
	fwrite(&flag, sizeof(int), 1, writefile);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ex1[i][j], sizeof(complex<double>), 1, writefile);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ey1[i][j], sizeof(complex<double>), 1, writefile);
		}
	fclose(writefile);

	//���ڸ���������2��ָ����Ϣ
	if (oneside) {//������2�������ƽ��
		R_Source = Matrix4D::getRotateMatrix(OutFieldGraphTrans.getRotate_theta(), OutFieldGraphTrans.getRotate_x(), OutFieldGraphTrans.getRotate_y(), OutFieldGraphTrans.getRotate_z())
			* Vector3(0, 0, 1.0);
	}
	else {//������ָ���뾵��ƽ��ָ����ͬ
		R_Source = n_Mirror;
	}

	//����������1�����䳡��poyntingʸ��
	Poynting();

	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(20, user);
	}
}

//��������������֮���б��������Ex_InΪ���룬��Ex1Ϊ���
void PVVAAS::CalArbitraryPlane(GraphTrans _GT0, GraphTrans _GTt, int mode) {
	//mode == 0 ����ת�ٴ��� �ʺϾ�ǰ����
	//mode == 1 �ȴ�������ת �ʺϾ��󴫲�
	GraphTrans GT0 = _GT0;
	GraphTrans GTt = _GTt;

	//�����䳡������ת����Ŀ��1һ�µ�����ϵ
	Vector3 U0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z())
		* Vector3(1.0, 0.0, 0.0);
	Vector3 V0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z())
		* Vector3(0.0, 1.0, 0.0);
	Vector3 N0 = Matrix4D::getRotateMatrix(GT0.getRotate_theta(), GT0.getRotate_x(), GT0.getRotate_y(), GT0.getRotate_z())
		* Vector3(0.0, 0.0, 1.0);
	Vector3 Ut = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())
		* Vector3(1.0, 0.0, 0.0);
	Vector3 Vt = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())
		* Vector3(0.0, 1.0, 0.0);
	Vector3 Nt = Matrix4D::getRotateMatrix(GTt.getRotate_theta(), GTt.getRotate_x(), GTt.getRotate_y(), GTt.getRotate_z())
		* Vector3(0.0, 0.0, 1.0);

	if (mode == 0) {//����ת�ٴ���
		if (Ut.Dot(Vector3(U0)) > 0.999 && Vt.Dot(Vector3(V0)) > 0.999) {
			//������ȫ��ͬ����ת��
		}
		else {
			FFTASRotation* RotField = new FFTASRotation;
			RotField->SetParas(this->f, this->N, this->M, this->dx, this->dy);
			RotField->SetRotationParas(GT0, GTt);
			RotField->SetInput(Ex_In, Ey_In);
			RotField->PerformRotate();
			RotField->output(Ex_In, Ey_In);
			delete RotField;
		}

		//�����䳡������������1��ע�⣬��������������1���ڵ�����ϵ��
		Vector3 dis; dis = Vector3(GTt.getTrans_x() - GT0.getTrans_x(), GTt.getTrans_y() - GT0.getTrans_y(), GTt.getTrans_z() - GT0.getTrans_z());
		//������������ת�����ȫ������ϵ�µ�disת��������1����ϵ

		double Up, Vp, Np;
		Up = Ut.Dot(dis);	Vp = Vt.Dot(dis);	Np = Nt.Dot(dis);
		FFTDI* FDI = new FFTDI(this->f, Up, Vp, Np, N, M);
		FDI->SetInput(Ex_In, Ey_In);
		FDI->Setds(dx,dy);
		FDI->StartCal();
		FDI->output(Ex1, Ey1, Ez1, Hx1, Hy1, Hz1);
		delete FDI;
	}
	else if (mode == 1) {
		//�����䳡������������1��ע�⣬��������������1���ڵ�����ϵ��
		Vector3 dis; dis = Vector3(GTt.getTrans_x() - GT0.getTrans_x(), GTt.getTrans_y() - GT0.getTrans_y(), GTt.getTrans_z() - GT0.getTrans_z());
		//������������ת�����ȫ������ϵ�µ�disת��������1����ϵ

		double Up, Vp, Np;
		Up = U0.Dot(dis);	Vp = V0.Dot(dis);	Np = N0.Dot(dis);
		FFTDI* FDI = new FFTDI(this->f, Up, Vp, Np, N, M);
		FDI->SetInput(Ex_In, Ey_In);
		FDI->Setds(dx,dy);
		FDI->StartCal();
		FDI->outsource(Ex_In, Ey_In);
		delete FDI;

		if (Ut.Dot(Vector3(U0)) > 0.999 && Vt.Dot(Vector3(V0)) > 0.999) {
			//������ȫ��ͬ����ת��
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < M; j++) {
					Ex1[i][j] = Ex_In[i][j];
					Ey1[i][j] = Ey_In[i][j];
				}
			}
		}
		else {
			FFTASRotation* RotField = new FFTASRotation;
			RotField->SetParas(this->f, this->N, this->M, this->dx,this->dy);
			RotField->SetRotationParas(GT0, GTt);
			RotField->SetInput(Ex_In, Ey_In);
			RotField->PerformRotate();
			RotField->output(Ex1, Ey1, Ez1);
			delete RotField;
		}
	}
}

void PVVAAS::Poynting()
{
	double absx, absy, absz;
	//double sum = 0;

	// Դ������ϵ��������1ת������������ϵ��ֻ����ת��
	Vector3D RotateAsixSou(SourceGraphTrans.getRotate_x(),
						   SourceGraphTrans.getRotate_y(),
		                   SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrixSou = Matrix4D::getRotateMatrix(
		SourceGraphTrans.getRotate_theta(), RotateAsixSou);

	double max = 0;
	double Power;

	Vector3 FieldsReal;
	Vector3 FieldsImag;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Px[i][j] = ConjugateMul(Ey1[i][j], Hz1[i][j]) - ConjugateMul(Ez1[i][j], Hy1[i][j]);
			Py[i][j] = ConjugateMul(Ez1[i][j], Hx1[i][j]) - ConjugateMul(Ex1[i][j], Hz1[i][j]);
			Pz[i][j] = ConjugateMul(Ex1[i][j], Hy1[i][j]) - ConjugateMul(Ey1[i][j], Hx1[i][j]);

			absx = Px[i][j].real();
			absy = Py[i][j].real();
			absz = Pz[i][j].real();

			n_Plane[i][j].set(absx, absy, absz);
			Power = n_Plane[i][j].Length();
			if (max < Power) max = Power;
			//����ӡ͢ʸ���ķ����������1�ľֲ�����ϵת����ȫ������ϵ
			n_Plane[i][j] = rotatMatrixSou * n_Plane[i][j];
			n_Plane[i][j].Normalization();
			//���糡��������1�ľֲ�����ϵת����ȫ������ϵ
			FieldsReal = Vector3(Ex1[i][j].real(), Ey1[i][j].real(), Ez1[i][j].real());
			FieldsImag = Vector3(Ex1[i][j].imag(), Ey1[i][j].imag(), Ez1[i][j].imag());

			FieldsReal = rotatMatrixSou *FieldsReal;
			FieldsImag = rotatMatrixSou *FieldsImag;

			Ex1[i][j] = complex<double>(FieldsReal.x, FieldsImag.x);
			Ey1[i][j] = complex<double>(FieldsReal.y, FieldsImag.y);
			Ez1[i][j] = complex<double>(FieldsReal.z, FieldsImag.z);

		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Px[i][j] = ConjugateMul(Ey1[i][j], Hz1[i][j]) - ConjugateMul(Ez1[i][j], Hy1[i][j]);
			Py[i][j] = ConjugateMul(Ez1[i][j], Hx1[i][j]) - ConjugateMul(Ex1[i][j], Hz1[i][j]);
			Pz[i][j] = ConjugateMul(Ex1[i][j], Hy1[i][j]) - ConjugateMul(Ey1[i][j], Hx1[i][j]);

			absx = Px[i][j].real();
			absy = Py[i][j].real();
			absz = Pz[i][j].real();

			Vector3 n_Power(absx, absy, absz);
			Power = n_Power.Length();
			if (Power>max*0.00000001)  Eno[i][j] = true;
			else Eno[i][j] = false;
		}
	}

}

void PVVAAS::ReflectCUDA() {
	//���빦�ܣ�����������������
	RayVec rayvec;

	//cout << "Calculating Reflecting surface" << endl;

	// ����ÿһ���� N * M * EleNum 
	Vector3 InterPoint;
	Vector3 Reflight, n_light; // ������� �� ������
	Vector3 n_light_Plane;  // �����ƽ��ķ�����

	// ��������ϵת����Դ������ϵ-������1��ֻ����ת��
	Vector3D RotateAsixSou(SourceGraphTrans.getRotate_x(),
		SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	double theta = SourceGraphTrans.getRotate_theta();
	Matrix4D rotatMatrixSou = Matrix4D::getRotateMatrix(
		-theta, RotateAsixSou);

	Vector3 tempa, tempb;

	double dir_t;  // ��ӡ͢ʸ�������ϵ��
	double plane_t;
	double d1, d2; //����·��
	double tempphase; // ��λ
	complex<double> tempejphase; // = exp(j*phase)
	RayTracing rayTracing(mirror);
	rayTracing.setthreadNum(threadNum);
	//rayTracing.SetupNodes();//������趨һ��Buffer ��ڵ�

	int returnNum = 0;
	int returnScale = N / 10;

	bool isInter = false;
	int num = N*M;
	vector<Vector3> startPiont(num);
	vector<Vector3> direction(num);
	vector<Vector3> nomal(num);
	vector<Vector3> intersection(num);
	vector<bool> isIntersect(num);
	vector<float> port(num);

	for (int i = 0; i < N; i++)//����
	{
		for (int j = 0; j < M; j++)//����
		{
			startPiont[i*N + j] = Plane[i][j];//������λ��-ȫ������
			direction[i*N + j] = n_Plane[i][j];//��������-ȫ������
		}
	}
	rayTracing.calcReflectBatch(startPiont, direction, nomal, intersection, isIntersect, port);

	//ȷ���迼�ǵĹ��߸���
	int RayNum = 0;
	for (int i = 0; i < N; i += 10) {
		for (int j = 0; j < N; j += 10) {
			if (isIntersect[i*N + j]) RayNum++;
		}
	}
	rayvec.Rays.resize(RayNum);
	//��ӹ���׷�����
	int countRay;
	countRay = 0;
	for (int i = 0; i < N; i += 10) {
		for (int j = 0; j < N; j += 10) {
			if (isIntersect[i*N + j]) {
				rayvec.AddPointToRay(countRay, startPiont[i*N + j]);
				countRay++;
			}
		}
	}
	//��ӵ����侵����
	countRay = 0;
	for (int i = 0; i < N; i += 10) {
		for (int j = 0; j < N; j += 10) {
			if (isIntersect[i*N + j]) {
				rayvec.AddPointToRay(countRay, intersection[i*N + j]);
				countRay++;
			}
		}
	}

	//��ɽ������
	for (int i = 0; i < N; i++)
	{
		if (returnFloat) // ���û��ע���򲻻ص�
		{
			if (i >= returnNum * returnScale)
			{
				returnFloat(20 + 5 * returnNum, user);
				returnNum++;
			}
		}
		for (int j = 0; j < M; j++)
		{
			int index = j + i*M;
			Ex_R[i][j] = 0;
			Ey_R[i][j] = 0;
			Ez_R[i][j] = 0;
			R_Plane[i][j] = 0;
			//���￪ʼ���������׷�ٽ���-���������ĵط���ֻ��Ҫ���������CUDA���ټ���
			if (!isIntersect[index])	continue;//û�ཻ�Ͳ�������
												 //����ľͶ����ཻ���ˣ�
			n_light = nomal[index];
			InterPoint = intersection[index];
			dir_t = port[index];
			n_light.Normalization();
									
			// ȫ������ϵ�����糡�ķ���
			CalReflectExyz(n_light, Ex1[i][j], Ey1[i][j],
				Ez1[i][j], Ex_R[i][j], Ey_R[i][j], Ez_R[i][j]);

			// �������-ȫ������-�������������-����㷨��-�õ������������
			Reflight = RayTracing::reflectLight(n_Plane[i][j], n_light);
			Reflight.Normalization();
			Rn_Plane[i][j] = Reflight;//ȫ������

			if (dir_t > 0.0000000001)  // ƽ���ڷ������ǰ��
			{
				plane_t = IntersectPlane(InterPoint, Reflight,
					Org_Source, R_Source, R_Plane[i][j]);
				d2 = CalDistance(InterPoint, R_Plane[i][j]);
				d1 = CalDistance(InterPoint, Plane[i][j]);

				if (plane_t > 0.0)  // ������2�ڷ������ǰ��
					tempphase = -(d1 + d2) / lamda * 2 * Pi;
				else
					tempphase = -(d1 - d2) / lamda * 2 * Pi;
				tempejphase = complex <double>(cos(tempphase), sin(tempphase));
				Ex_R[i][j] = ComMul(Ex_R[i][j], tempejphase);  // ֻ����λ�任
				Ey_R[i][j] = ComMul(Ey_R[i][j], tempejphase);
				Ez_R[i][j] = ComMul(Ez_R[i][j], tempejphase);
			}
			else if (dir_t < -0.0000000001)   // ƽ���ڷ�����ĺ���
			{
				plane_t = IntersectPlane(InterPoint, Reflight,
					Org_Source, R_Source, R_Plane[i][j]);
				d2 = CalDistance(InterPoint, R_Plane[i][j]);
				d1 = CalDistance(InterPoint, Plane[i][j]);

				if (plane_t < 0.0)  // ������2�ڷ�����ĺ���
					tempphase = (d1 + d2) / lamda * 2 * Pi;
				else
					tempphase = (d1 - d2) / lamda * 2 * Pi;
				tempejphase = complex <double>(cos(tempphase), sin(tempphase));
				Ex_R[i][j] = ComMul(Ex_R[i][j], tempejphase); // ֻ����λ�任
				Ey_R[i][j] = ComMul(Ey_R[i][j], tempejphase);
				Ez_R[i][j] = ComMul(Ez_R[i][j], tempejphase);
			}
			else  // ƽ���뷴�����ཻ
			{
				R_Plane[i][j] = InterPoint;
			}
		}
	} // endloop

	  //��ӵ����侵����
	countRay = 0;
	for (int i = 0; i < N; i += 10) {
		for (int j = 0; j < N; j += 10) {
			if (isIntersect[i*N + j]) {
				rayvec.AddPointToRay(countRay, R_Plane[i][j]);
				countRay++;
			}
		}
	}

	rayvec.writeRays2BiFile("./pvvaasrays.dat");

	//������ֵ-����ȫ�־�������ϵ�µ�
	CalAmplitude();  // ֻ�����ȱ任-�����õ���n_source ��R_source


	//Դ�Ĵ�������ı�-�Ӵ���������2����-����������1��������2�ǹ����ĵģ��������ﲻ��Ҫ����λ��
	updateSource_n(R_Source);//Դ��GT�Ѹı�
	//Ȼ�󽫳�������ȫ�־�������ϵ�任��������2�ľֲ��������ϵ
	Vector3 Ut = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(), SourceGraphTrans.getRotate_z())*Vector3(1.0,0.0,0.0);
	Vector3 Vt = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(), SourceGraphTrans.getRotate_z())*Vector3(0.0,1.0, 0.0);
	Vector3 Nt = Matrix4D::getRotateMatrix(SourceGraphTrans.getRotate_theta(), SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(), SourceGraphTrans.getRotate_z())*Vector3(0.0,0.0, 1.0);
	
	//������ת
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			Vector3 tempReal(Ex_R[i][j].real(), Ey_R[i][j].real(), Ez_R[i][j].real());
			Vector3 tempImag(Ex_R[i][j].imag(), Ey_R[i][j].imag(), Ez_R[i][j].imag());
			if (tempReal.Length() > 0.00001) {
				int check = 0;
			}
			Ex_R[i][j] = complex<double>(Ut.Dot(tempReal), Ut.Dot(tempImag));
			Ey_R[i][j] = complex<double>(Vt.Dot(tempReal), Vt.Dot(tempImag));
			Ez_R[i][j] = complex<double>(Nt.Dot(tempReal), Nt.Dot(tempImag));
		}
	}
	//�ͷ��ڴ�
	// �õ����
//	delete[] intersected; 
//	delete[] prot;

//	delete[] inter_x;	delete[] inter_y;	delete[] inter_z;
//	delete[] STLIndex;
//	delete[] point;

	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(80, user);
	}

}

void PVVAAS::CalAmplitude()
{
	Vector3 tempA1, tempA2;
	Vector3 tempB1, tempB2;
	Vector3 tempC1, tempC2;
	Vector3 tempD1, tempD2;
	double tempcos1, tempcos2;
	double PreSquare, PostSquare;
	double tempratio;
	double LAE, TAE;
	n_Source.Normalization();
	//���Ĳ����ĸ�����
	for (int i = 1; i < N - 2; i++)
	{
		for (int j = 1; j < M - 2; j++)
		{
			// ����ǰͶӰ���
			tempA1 = Plane[i - 1][j];
			tempB1 = Plane[i][j + 1];
			tempC1 = Plane[i + 1][j];
			tempD1 = Plane[i][j - 1];
			Vector3 direction = n_Plane[i][j];
			tempcos1 = direction.Dot(n_Source);//������1�ϸ�������������ķ���֮��ĵ��

			LAE = pow((Ex1[i][j].real() * Ex1[i][j].real() + Ex1[i][j].imag() * Ex1[i][j].imag()
				+ Ey1[i][j].real() * Ey1[i][j].real() + Ey1[i][j].imag() * Ey1[i][j].imag()
				+ Ez1[i][j].real() * Ez1[i][j].real() + Ez1[i][j].imag() * Ez1[i][j].imag()), 0.5);

			double Square1 = CalSquare(tempA1, tempB1, tempC1, tempD1);
			PreSquare = Square1*pow(tempcos1,1.0);


			// �����ͶӰ���
			tempA2 = R_Plane[i - 1][j];
			tempB2 = R_Plane[i][j + 1];
			tempC2 = R_Plane[i + 1][j];
			tempD2 = R_Plane[i][j - 1];
			direction = Rn_Plane[i][j];
			tempcos2 = direction.Dot(R_Source);//������2�������������ķ���֮��ĵ��

			TAE = pow((Ex_R[i][j].real() * Ex_R[i][j].real() + Ex_R[i][j].imag() * Ex_R[i][j].imag()
				+ Ey_R[i][j].real() * Ey_R[i][j].real() + Ey_R[i][j].imag() * Ey_R[i][j].imag()
				+ Ez_R[i][j].real() * Ez_R[i][j].real() + Ez_R[i][j].imag() * Ez_R[i][j].imag()), 0.5);

			double Square2 = CalSquare(tempA2, tempB2, tempC2, tempD2);
			PostSquare = Square2*pow(tempcos2,1.0);
			if (Square2 > 0) {
				int check = 0;
			}
			//PostSquare = Square2;
			/*
			if (abs(SquareTest1 - SquareTest2) < 0.00000001)
				PostSquare = tempcos2 * SquareTest1;
			else
			{
				PostSquare = PreSquare;
			}
			if (PostSquare == 0)
				PostSquare = 1;
			if (TAE == 0)
				TAE = 1;
			*/
			tempratio = pow(fabs(PreSquare / PostSquare), 0.5);// *LAE / TAE;

			Ex_R[i][j] = Ex_R[i][j] * tempratio;
			Ey_R[i][j] = Ey_R[i][j] * tempratio;
			Ez_R[i][j] = Ez_R[i][j] * tempratio;
		}
	}
}


void PVVAAS::InterVal()
{

	//��ȫ������ϵת��������ͶӰ��ľֲ�����ϵ��
	// �����
	Vector3D RotateAsix(SourceGraphTrans.getRotate_x(),
		SourceGraphTrans.getRotate_y(),
		SourceGraphTrans.getRotate_z());
	Matrix4D R_rotatMatrix = Matrix4D::getRotateMatrix(
		-SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D R_translateMatrix = Matrix4D::getTranslateMatrix(
		-SourceGraphTrans.getTrans_x(),
		-SourceGraphTrans.getTrans_y(), -SourceGraphTrans.getTrans_z());
	Matrix4D R_Matrix = R_rotatMatrix * R_translateMatrix;

	//��ȫ������ϵת�������ճ������ľֲ�����ϵ��
	Vector3D RotateAsixOut(OutFieldGraphTrans.getRotate_x(),
		OutFieldGraphTrans.getRotate_y(),
		OutFieldGraphTrans.getRotate_z());
	Matrix4D R_rotatMatrixOut = Matrix4D::getRotateMatrix(
		-OutFieldGraphTrans.getRotate_theta(), RotateAsixOut);
	Vector3 tempRPlaneOut;
	complex<double> tempc;
	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			if (R_Plane[i][j].x == 100) {}
			else {
				tempRPlaneOut = R_rotatMatrixOut * R_Plane[i][j];
				tempc = exp(complex<double>(0.0, 2 * Pi / lamda*tempRPlaneOut.z));
				Ex_R[i][j] = Ex_R[i][j] * tempc;
				Ey_R[i][j] = Ey_R[i][j] * tempc;
			}
		}


	bool ** hit = NULL;//�ж���û�ཻ
	hit = Allocate_2D(hit, N, M);
	double maxValue = 0.0;
	double Value;
	Vector3 tempRPlane;
	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			Value = sqrt(Ex_R[i][j].real()*Ex_R[i][j].real() + Ex_R[i][j].imag()*Ex_R[i][j].imag() + Ey_R[i][j].real()*Ey_R[i][j].real() + Ey_R[i][j].imag()*Ey_R[i][j].imag());
			if (Value > maxValue) maxValue = Value;

		}

	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			Value = sqrt(Ex_R[i][j].real()*Ex_R[i][j].real() + Ex_R[i][j].imag()*Ex_R[i][j].imag() + Ey_R[i][j].real()*Ey_R[i][j].real() + Ey_R[i][j].imag()*Ey_R[i][j].imag());

			if (R_Plane[i][j].x == 100)
				hit[i][j] = false;
			else if (Value > maxValue*0.00001) 
				hit[i][j] = true;
			else hit[i][j] = false;
				
			tempRPlane = R_Matrix * R_Plane[i][j];
			R_Plane[i][j] = tempRPlane;
		}

	//����ͶӰ��ֲ�����ϵ��
	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			Ex_In[i][j] = complex<double>(0, 0);
			Ey_In[i][j] = complex<double>(0, 0);
			Position3D tempPoint((i - (N - 1) / 2) * dx, (j - (M - 1) / 2) * dy, 0);
			Plane[i][j].set(tempPoint.X(), tempPoint.Y(), tempPoint.Z());
		}

	FILE* writefile;
	string fileAddress = "./PVVD/outas0.dat";
	errno_t err;
	err = fopen_s(&writefile, fileAddress.c_str(), "wb");
	if (err != 0) return;
	double tx, ty, tz, rx, ry, rz, rth;
	int flag = 0;
	tx = SourceGraphTrans.getTrans_x(); ty = SourceGraphTrans.getTrans_y(); tz = SourceGraphTrans.getTrans_z();
	rx = SourceGraphTrans.getRotate_x(); ry = SourceGraphTrans.getRotate_y(); rz = SourceGraphTrans.getRotate_z();
	rth = SourceGraphTrans.getRotate_theta();
	fwrite(&f, sizeof(double), 1, writefile);
	fwrite(&tx, sizeof(double), 1, writefile); fwrite(&ty, sizeof(double), 1, writefile); fwrite(&tz, sizeof(double), 1, writefile);
	fwrite(&rx, sizeof(double), 1, writefile); fwrite(&ry, sizeof(double), 1, writefile); fwrite(&rz, sizeof(double), 1, writefile);
	fwrite(&rth, sizeof(double), 1, writefile);
	fwrite(&N, sizeof(int), 1, writefile);
	fwrite(&M, sizeof(int), 1, writefile);
	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);
	fwrite(&flag, sizeof(int), 1, writefile);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ex_R[i][j], sizeof(complex<double>), 1, writefile);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ey_R[i][j], sizeof(complex<double>), 1, writefile);
		}
	fclose(writefile);

	if (CATRMode) {
		int NumSlice = N / 100;
		if (NumSlice == 0) NumSlice = 1;
		vector<int> NStart;
		vector<int> NNum;
		NStart.resize(NumSlice);		NNum.resize(NumSlice);

		for (int ii = 0; ii < NumSlice; ii++) {
			int nn = N / NumSlice;
			int resn = N%NumSlice;
			if ((ii < resn) && (resn != 0)) NNum[ii] = nn + 1;
			else NNum[ii] = nn;
		}
		for (int ii = 0; ii < NumSlice; ii++) {
			if (ii == 0) NStart[ii] = 0;
			else	NStart[ii] = NStart[ii - 1] + NNum[ii - 1];
		}
		//omp acceleration Preparation
		int totalNum = NumSlice*NumSlice;
		vector<int> startS;	startS.resize(threadNum);
		vector<int> numS;	numS.resize(threadNum);
		vector<int> stopS;	stopS.resize(threadNum);
		for (int i = 0; i < threadNum; i++) {
			int num = totalNum / threadNum;
			int rr = (totalNum) % threadNum;
			if (i<rr && rr != 0) num += 1;
			numS[i] = num;
			//starti
			if (i == 0) startS[i] = 0;
			else startS[i] = startS[i - 1] + numS[i - 1];
			//stopi
			stopS[i] = startS[i] + numS[i];
		}

		omp_set_num_threads(threadNum);
#pragma omp parallel
		{
			int id = omp_get_thread_num();
			for (int ip = 0; ip < numS[id]; ip++) {
				int io = ip*threadNum + id;

				int is = io / NumSlice;
				int js = io % NumSlice;

				if (id == 0) { if (returnFloat) returnFloat(((ip - startS[id] + 0.5) / numS[id] * 20 + 60), user); }
				vtkFloatArray *scalars = vtkFloatArray::New();
				vtkSmartPointer<vtkPoints> points =
					vtkSmartPointer<vtkPoints>::New();

				double MAXX = dx * (N - 1) / 2;
				double MAXY = dy * (M - 1) / 2;


				int NxStart = NStart[is];	if (is > 0) NxStart = NxStart - 10;
				int NxEnd = NStart[is] + NNum[is];	if (is < NumSlice - 1) NxEnd = NxEnd + 10;
				//int NxNum = NxEnd - NxStart;

				int NyStart = NStart[js];	if (js > 0) NyStart = NyStart - 10;
				int NyEnd = NStart[js] + NNum[js];	if (js < NumSlice - 1) NyEnd = NyEnd + 10;
				//int NyNum = NyEnd - NyStart;

				int num = 0;
				//���Ƶ㼯-�������λ�� ������������ת����ͶӰ������ϵ��
				for (int i = NxStart; i < NxEnd; i++) {
					for (int j = NyStart; j < NyEnd; j++) {
						if (!hit[i][j])
							continue;
						if (R_Plane[i][j].x < -MAXX || R_Plane[i][j].x > MAXX ||
							R_Plane[i][j].y < -MAXY || R_Plane[i][j].y > MAXY)
							continue;
						points->InsertNextPoint(R_Plane[i][j].x, R_Plane[i][j].y, 0);
						scalars->InsertTuple1(num, i*M + j + 0.01);//��������Scalar��ʵ��Ϊ�˱�ǵ��index
						num++;
					}
				}

				vtkSmartPointer<vtkPolyData> polydata =
					vtkSmartPointer<vtkPolyData>::New();
				polydata->SetPoints(points);
				polydata->GetPointData()->SetScalars(scalars);

				double scalarsRange[2];
				scalars->GetRange(scalarsRange);
				scalars->Delete();

				vtkSmartPointer<vtkDelaunay2D> delaunay =
					vtkSmartPointer<vtkDelaunay2D>::New();
				delaunay->SetInputData(polydata);
				delaunay->Update();
				polydata = delaunay->GetOutput();

				int EleNum = polydata->GetNumberOfCells();
				//CATRģʽ�˴�����
				double ds = (dx + dy)*0.5;
				float Maxlength;
				vector<int> LongTriangleId;	LongTriangleId.resize(0);
				vector<Vector2> Trip; Trip.resize(3);
				vector<double> ll;	ll.resize(3);
				for (int i = 0; i < EleNum; i++) {
					vtkIdList * p;
					p = polydata->GetCell(i)->GetPointIds();
					double * point;
					point = polydata->GetPoint(p->GetId(0));	Trip[0] = Vector2(point[0], point[1]);
					point = polydata->GetPoint(p->GetId(1));	Trip[1] = Vector2(point[0], point[1]);
					point = polydata->GetPoint(p->GetId(2));	Trip[2] = Vector2(point[0], point[1]);
					ll[0] = Vector2(Trip[0].x - Trip[1].x, Trip[0].y - Trip[1].y).Length();
					ll[1] = Vector2(Trip[1].x - Trip[2].x, Trip[1].y - Trip[2].y).Length();
					ll[2] = Vector2(Trip[2].x - Trip[0].x, Trip[2].y - Trip[0].y).Length();
					Maxlength = 0;
					for (int k = 0; k < 3; k++) {
						if (ll[k] > Maxlength) Maxlength = ll[k];
					}
					if (Maxlength > 2 * ds) {
						LongTriangleId.insert(LongTriangleId.end(), i);
					}
				}
				int indexNum = LongTriangleId.size();
				for (int i = 0; i < indexNum; i++) {
					polydata->DeleteCell(LongTriangleId[indexNum - 1 - i]);
				}
				polydata->RemoveDeletedCells();
				polydata->Modified();
				
				EleNum = polydata->GetNumberOfCells();

				for (int i = 0; i < EleNum; i++)  //���뷴����Ľ���
				{
					double maxX = -MAXX;
					double minX = MAXX;
					double maxY = -MAXY;
					double minY = MAXY;

					vtkIdList * p;
					p = polydata->GetCell(i)->GetPointIds();
					double * point;

					vector<Vector2> A(3);
					vector<int> dataIJ(3);
					Vector2 P;
					for (int k = 0; k < 3; k++)
					{
						dataIJ[k] = polydata->GetPointData()
							->GetScalars()->GetVariantValue(p->GetId(k)).ToInt();//ȡ��Scalar
						point = polydata->GetPoint(p->GetId(k));//�ҵ����λ��
						if (maxX < point[0])
							maxX = point[0];
						if (minX > point[0])
							minX = point[0];
						if (maxY < point[1])
							maxY = point[1];
						if (minY > point[1])
							minY = point[1];
						A[k].set(point[0], point[1]);
					}
					int iiMin = ceil(minX / ds) + (N - 1) / 2;
					int iiMax = floor(maxX / ds) + (N - 1) / 2;
					int jjMin = ceil(minY / ds) + (M - 1) / 2;
					int jjMax = floor(maxY / ds) + (M - 1) / 2;
					if (iiMin < 0)	iiMin = 0;
					if (jjMin < 0)	jjMin = 0;
					register complex<double> tempEx;
					register complex<double> tempEy;
					for (int ii = iiMin; ii <= iiMax&&ii < N; ii++)
						for (int jj = jjMin; jj <= jjMax&&jj < M; jj++)
						{
							P.set(Plane[ii][jj].x, Plane[ii][jj].y);
							if (InterVal_PointinTriangle(A[0], A[1], A[2], P))//������ȸ��λ���������������棬Ȼ��������������������Բ�ֵ
							{
								double a = (P - A[0]).area();
								double b = (P - A[1]).area();
								double c = (P - A[2]).area();

								tempEx = (
									b * c * Ex_R[dataIJ[0] / M][dataIJ[0] % M] +
									a * c * Ex_R[dataIJ[1] / M][dataIJ[1] % M] +
									b * a * Ex_R[dataIJ[2] / M][dataIJ[2] % M]) /
									(b * c + a * c + b * a);

								tempEy = (
									b * c * Ey_R[dataIJ[0] / M][dataIJ[0] % M] +
									a * c * Ey_R[dataIJ[1] / M][dataIJ[1] % M] +
									b * a * Ey_R[dataIJ[2] / M][dataIJ[2] % M]) /
									(b * c + a * c + b * a);
								if (abs(tempEx) + abs(tempEy) > 1e-5) {
									Ex_In[ii][jj] = tempEx;
									Ey_In[ii][jj] = tempEy;
								}
								else {
									ii = ii;
								}
							}//if
						}//two for
				}//for i  EleNum
			}
		}//omp
	}
	else {
	
		//No Omp acceleration
		int num = 0;
		vtkFloatArray *scalars = vtkFloatArray::New();
		vtkSmartPointer<vtkPoints> points =	vtkSmartPointer<vtkPoints>::New();

		double MAXX = dx * (N - 1) / 2;
		double MAXY = dy * (M - 1) / 2;
		//���Ƶ㼯-�������λ�� ������������ת����ͶӰ������ϵ��
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				if (!hit[i][j])
					continue;
				if (R_Plane[i][j].x < -MAXX || R_Plane[i][j].x > MAXX ||
					R_Plane[i][j].y < -MAXY || R_Plane[i][j].y > MAXY)
					continue;
				points->InsertNextPoint(R_Plane[i][j].x, R_Plane[i][j].y, 0);
				scalars->InsertTuple1(num, i*M + j + 0.01);//��������Scalar��ʵ��Ϊ�˱�ǵ��index
				num++;
			}
		}

		vtkSmartPointer<vtkPolyData> polydata =
			vtkSmartPointer<vtkPolyData>::New();
		polydata->SetPoints(points);
		polydata->GetPointData()->SetScalars(scalars);

		double scalarsRange[2];
		scalars->GetRange(scalarsRange);
		scalars->Delete();

		vtkSmartPointer<vtkDelaunay2D> delaunay =
			vtkSmartPointer<vtkDelaunay2D>::New();
		delaunay->SetInputData(polydata);
		delaunay->Update();
		polydata = delaunay->GetOutput();

		int EleNum = polydata->GetNumberOfCells();
		float Maxlength;
		double ds = (dx + dy)*0.5;
		vector<int> LongTriangleId;	LongTriangleId.resize(0);
		vector<Vector2> Trip; Trip.resize(3);
		vector<double> ll;	ll.resize(3);
		for (int i = 0; i < EleNum; i++) {
			vtkIdList * p;
			p = polydata->GetCell(i)->GetPointIds();
			double * point;
			point = polydata->GetPoint(p->GetId(0));	Trip[0] = Vector2(point[0], point[1]);
			point = polydata->GetPoint(p->GetId(1));	Trip[1] = Vector2(point[0], point[1]);
			point = polydata->GetPoint(p->GetId(2));	Trip[2] = Vector2(point[0], point[1]);
			ll[0] = Vector2(Trip[0].x - Trip[1].x, Trip[0].y - Trip[1].y).Length();
			ll[1] = Vector2(Trip[1].x - Trip[2].x, Trip[1].y - Trip[2].y).Length();
			ll[2] = Vector2(Trip[2].x - Trip[0].x, Trip[2].y - Trip[0].y).Length();
			Maxlength = 0;
			for (int k = 0; k < 3; k++) {
				if (ll[k] > Maxlength) Maxlength = ll[k];
			}
			if (Maxlength > 5 * ds) {
				LongTriangleId.insert(LongTriangleId.end(), i);
			}
		}
		int indexNum = LongTriangleId.size();
		for (int i = 0; i < indexNum; i++) {
			polydata->DeleteCell(LongTriangleId[indexNum - 1 - i]);
		}
		polydata->RemoveDeletedCells();
		polydata->Modified();


		vtkSmartPointer<vtkSTLWriter> writer =
			vtkSmartPointer<vtkSTLWriter>::New();
		writer->SetInputData(polydata);
		writer->SetFileName("./PVVD/testpvvaas.stl");
		writer->Update();

		EleNum = polydata->GetNumberOfCells();

		for (int i = 0; i < EleNum; i++)  //���뷴����Ľ���
		{
			double maxX = -MAXX;
			double minX = MAXX;
			double maxY = -MAXY;
			double minY = MAXY;

			vtkIdList * p;
			p = polydata->GetCell(i)->GetPointIds();
			double * point;

			vector<Vector2> A(3);
			vector<int> dataIJ(3);
			Vector2 P;
			for (int k = 0; k < 3; k++)
			{
				dataIJ[k] = polydata->GetPointData()
					->GetScalars()->GetVariantValue(p->GetId(k)).ToInt();//ȡ��Scalar
				point = polydata->GetPoint(p->GetId(k));//�ҵ����λ��
				if (maxX < point[0])
					maxX = point[0];
				if (minX > point[0])
					minX = point[0];
				if (maxY < point[1])
					maxY = point[1];
				if (minY > point[1])
					minY = point[1];
				A[k].set(point[0], point[1]);
			}
			int iiMin = ceil(minX / dx) + (N - 1) / 2;
			int iiMax = floor(maxX / dx) + (N - 1) / 2;
			int jjMin = ceil(minY / dy) + (M - 1) / 2;
			int jjMax = floor(maxY / dy) + (M - 1) / 2;
			if (iiMin < 0)	iiMin = 0;
			if (jjMin < 0)	jjMin = 0;
			register complex<double> tempEx;
			register complex<double> tempEy;
			for (int ii = iiMin; ii <= iiMax&&ii < N; ii++)
				for (int jj = jjMin; jj <= jjMax&&jj < M; jj++)
				{
					P.set(Plane[ii][jj].x, Plane[ii][jj].y);
					if (InterVal_PointinTriangle(A[0], A[1], A[2], P))//������ȸ��λ���������������棬Ȼ��������������������Բ�ֵ
					{
						double a = (P - A[0]).area();
						double b = (P - A[1]).area();
						double c = (P - A[2]).area();

						tempEx = (
							b * c * Ex_R[dataIJ[0] / M][dataIJ[0] % M] +
							a * c * Ex_R[dataIJ[1] / M][dataIJ[1] % M] +
							b * a * Ex_R[dataIJ[2] / M][dataIJ[2] % M]) /
							(b * c + a * c + b * a);

						tempEy = (
							b * c * Ey_R[dataIJ[0] / M][dataIJ[0] % M] +
							a * c * Ey_R[dataIJ[1] / M][dataIJ[1] % M] +
							b * a * Ey_R[dataIJ[2] / M][dataIJ[2] % M]) /
							(b * c + a * c + b * a);
						if (abs(tempEx) + abs(tempEy) > 1e-5) {
							Ex_In[ii][jj] = tempEx;
							Ey_In[ii][jj] = tempEy;
						}
						else {
							ii = ii;
						}
					}//if
				}//two for
		}//for i  EleNum
	}


	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(70, user);
	}

	Free_2D(hit);	hit = NULL;


	//�����󣬰��������Ͼֲ�����ϵ��λ�ñ任��ȫ�������ֵ
	//RotateAsix = Vector3D(SourceGraphTrans.getRotate_x(), SourceGraphTrans.getRotate_y(), SourceGraphTrans.getRotate_z());
	Matrix4D rotatMatrix = Matrix4D::getRotateMatrix(
		SourceGraphTrans.getRotate_theta(), RotateAsix);
	Matrix4D translateMatrix = Matrix4D::getTranslateMatrix(
		SourceGraphTrans.getTrans_x(),
		SourceGraphTrans.getTrans_y(), SourceGraphTrans.getTrans_z());
	// ����Դ�������λ�� ������2��
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			Plane[i][j] = translateMatrix * rotatMatrix * Plane[i][j];
		}

	for (int i = 0; i<N; i++)
		for (int j = 0; j < M; j++) {
			if (abs(Ex_R[i][j]) + abs(Ey_In[i][j])>0.0) {
				tempRPlaneOut = R_rotatMatrixOut * Plane[i][j];
				tempc = exp(complex<double>(0.0, -2 * Pi / lamda*tempRPlaneOut.z));
				Ex_In[i][j] = Ex_In[i][j] * tempc;
				Ey_In[i][j] = Ey_In[i][j] * tempc;
			}

		}

	fileAddress = "./PVVD/outas1.dat";
	err = fopen_s(&writefile, fileAddress.c_str(), "wb");
	if (err != 0) return;
	tx = SourceGraphTrans.getTrans_x(); ty = SourceGraphTrans.getTrans_y(); tz = SourceGraphTrans.getTrans_z();
	rx = SourceGraphTrans.getRotate_x(); ry = SourceGraphTrans.getRotate_y(); rz = SourceGraphTrans.getRotate_z();
	rth = SourceGraphTrans.getRotate_theta();
	fwrite(&f, sizeof(double), 1, writefile);
	fwrite(&tx, sizeof(double), 1, writefile); fwrite(&ty, sizeof(double), 1, writefile); fwrite(&tz, sizeof(double), 1, writefile);
	fwrite(&rx, sizeof(double), 1, writefile); fwrite(&ry, sizeof(double), 1, writefile); fwrite(&rz, sizeof(double), 1, writefile);
	fwrite(&rth, sizeof(double), 1, writefile);
	fwrite(&N, sizeof(int), 1, writefile);
	fwrite(&M, sizeof(int), 1, writefile);
	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);
	fwrite(&flag, sizeof(int), 1, writefile);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ex_In[i][j], sizeof(complex<double>), 1, writefile);
		}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			fwrite(&Ey_In[i][j], sizeof(complex<double>), 1, writefile);
		}
	fclose(writefile);


	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(80, user);
	}
	//ע�⣬���ｫ������2�ϵĳ���Ϊ��һ�������������ʼԴ
}

void PVVAAS::Result()
{
	CalArbitraryPlane(SourceGraphTrans,OutFieldGraphTrans,1);
	if (returnFloat) // ���û��ע���򲻻ص�
	{
		returnFloat(90, user);
	}
}

double PVVAAS::IntersectPlane(const Vector3 & orig, const Vector3 & dir,
	const Vector3 & Plane_org, const Vector3 & Plane_n, Vector3 & intersection)
{
	double temp1 = dir.Dot(Plane_n);
	if(temp1 == 0)
		return 0.0;
	double temp2 = Plane_org.Dot(Plane_n) - orig.Dot(Plane_n);
	double t = temp2 / temp1;
	intersection = orig + dir * t;
	if (t > 10)
	{
		double l = 0;
	}
	return t;

}

double PVVAAS::IntersectPoint(const Vector3 &orig, const Vector3 &dir,
	const Vector3 &v0, const Vector3 &E1, const Vector3 &E2, Vector3 &intersection)
{
	// P
	Vector3 P = dir.Cross(E2);

	// determinant
	double det = E1.Dot(P);

	double u, v, t;
	Vector3 T;
	T = orig - v0;

	// If determinant is near zero, ray lies in plane of triangle
	if (det < 10e-10 && det > -10e-10)
		return 0;

	u = T.Dot(P);
	double fInvDet = 1.0f / det;
	u *= fInvDet;

	// Q
	Vector3 Q = T.Cross(E1);

	v = dir.Dot(Q);
	v *= fInvDet;

	// Calculate t, scale parameters, ray intersects triangle
	t = E2.Dot(Q);
	t *= fInvDet;

	intersection = orig + dir * t;
	return t;
}

void PVVAAS::getPlane(Vector3 **& _org, Vector3 **& _n) const
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
		{
			_org[i][j] = Plane[i][j];
			_n[i][j] = n_Plane[i][j];
		}
}

void PVVAAS::SetReturnFloat(void(*returnFloat)(float, void *), void * _user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}

void PVVAAS::updateSource_n(const Vector3& new_n)
{
	//new_n.Normalization();
	if (new_n.x != 0 || new_n.y != 0 || new_n.z != 1)
	{
		Vector3 rotate_axis = Vector3(0, 0, 1).Cross(new_n); // ��ת��
		rotate_axis.Normalization();
		double rotate_theta = acos(Vector3(0, 0, 1).Dot(new_n));
		rotate_theta = rotate_theta / Pi * 180;
		SourceGraphTrans.updateRotate(rotate_axis, rotate_theta);
	}
	else
	{
		SourceGraphTrans.updateRotate(Vector3(0, 0, 1), 0);
	}
}

void PVVAAS::CalReflectExyz(const Vector3 & n, const complex<double>& Ex,
	const complex<double>& Ey, const complex<double>& Ez, complex<double>& Ex_out,
	complex<double>& Ey_out, complex<double>& Ez_out)
{
	/*
	Vector3 CrossRe;	Vector3 CrossIm;
	CrossRe = n.Cross(Vector3(Ex.real(), Ey.real(), Ez.real()));
	CrossIm = n.Cross(Vector3(Ex.imag(), Ey.imag(), Ez.imag()));

	Ex_out = Ex + 2.0*complex<double>(CrossRe.x, CrossIm.x);
	Ey_out = Ey + 2.0*complex<double>(CrossRe.y, CrossIm.y);
	Ez_out = Ez + 2.0*complex<double>(CrossRe.z, CrossIm.z);
*/
	
	complex<double> ne(n.x * Ex.real() + n.y * Ey.real() + n.z * Ez.real(),
		n.x * Ex.imag() + n.y * Ey.imag() + n.z * Ez.imag());
	complex<double> tempx = 2 * n.x * ne;
	complex<double> tempy = 2 * n.y * ne;
	complex<double> tempz = 2 * n.z * ne;

	Ex_out = tempx - Ex;
	Ey_out = tempy - Ey;
	Ez_out = tempz - Ez;
	
}

double PVVAAS::CalSquare(const Vector3 & A, const Vector3 & B, const Vector3 & C, const Vector3 & D) const
{
	double AB = (B - A).Length();
	double AC = (C - A).Length();
	double AD = (D - A).Length();
	double BC = (C - B).Length();
	double DC = (C - D).Length();
	double p1 = (AB + AC + BC) / 2;
	double p2 = (AD + AC + DC) / 2;

	double tempS1 = sqrt(p1*(p1 - AB)*(p1 - BC)*(p1 - AC));
	double tempS2 = sqrt(p2*(p2 - AC)*(p2 - AD)*(p2 - DC));

	return tempS1 + tempS2;
}

//�������������
double PVVAAS::CalSquare(const Vector3 & A, const Vector3 & B,
	const Vector3 & C) const
{
	double AB = (B - A).Length();
	double AC = (C - A).Length();
	double BC = (B - C).Length();
	double p1 = (AB + AC + BC) / 2;
	
	return sqrt(p1*(p1 - AB)*(p1 - BC)*(p1 - AC));
}

bool PVVAAS::get_line_intersection(
	const Vector2 & A, const Vector2 & B, 
	const Vector2 & C, const Vector2 & D, Vector2 & O)
{
	double s10_x = B.x - A.x;
	double s10_y = B.y - A.y;
	double s32_x = D.x - C.x;
	double s32_y = D.y - C.y;

	double denom = s10_x * s32_y - s32_x * s10_y;
	if (abs(denom) < 0.000001)//ƽ�л���
		return 0; // Collinear
	bool denomPositive = denom > 0.0;

	double s02_x = A.x - C.x;
	double s02_y = A.y - C.y;
	double s_numer = s10_x * s02_y - s10_y * s02_x;
	if ((s_numer < 0.0) == denomPositive)
		//�����Ǵ��ڵ���0��С�ڵ���1�ģ����ӷ�ĸ����ͬ���ҷ���С�ڵ��ڷ�ĸ
		return false; // No collision

	double t_numer = s32_x * s02_y - s32_y * s02_x;
	if ((t_numer < 0.0) == denomPositive)
		return false; // No collision

	if (((s_numer > denom) == denomPositive) || ((t_numer > denom) == denomPositive))
		return false; // No collision
				  // Collision detected
	double t = t_numer / denom;
	
	O.set(A.x + (t * s10_x), A.y + (t * s10_y));

	return true;
}

void PVVAAS::AllocateMemory()
{
	//��ӡ͢ʸ��
	Px = Allocate_2D(Px, N, M);
	Py = Allocate_2D(Py, N, M);
	Pz = Allocate_2D(Pz, N, M);

	//Դ
	Ex_In = Allocate_2D(Ex_In, N, M);
	Ey_In = Allocate_2D(Ey_In, N, M);

	//������ĵ�ų�
	Ex1 = Allocate_2D(Ex1, N, M);
	Ey1 = Allocate_2D(Ey1, N, M);
	Ez1 = Allocate_2D(Ez1, N, M);

	Hx1 = Allocate_2D(Hx1, N, M);
	Hy1 = Allocate_2D(Hy1, N, M);
	Hz1 = Allocate_2D(Hz1, N, M);

	//�����ĵ糡
	Ex_R = Allocate_2D(Ex_R, N, M);
	Ey_R = Allocate_2D(Ey_R, N, M);
	Ez_R = Allocate_2D(Ez_R, N, M);

	// ƽ��λ�ã������ȣ�����󣩺͸���ķ�����
	n_Plane = Allocate_2D(n_Plane, N, M);
	Plane = Allocate_2D(Plane, N, M);
	R_Plane = Allocate_2D(R_Plane, N, M);
	Rn_Plane = Allocate_2D(Rn_Plane, N, M);

	Eno = Allocate_2D(Eno, N, M);

}

void PVVAAS::getSourcePoint(Vector3 & interPoint, Vector3 & n_Point) const
{
	interPoint = Org_Source;
	n_Point = n_Source;
}

void PVVAAS::setSourcePoint(const Vector3 & interPoint, const Vector3 & n_Point)
{
	Org_Source = interPoint;
	n_Source = n_Point;
}

GraphTrans PVVAAS::getSourceGraphTrans(const Vector3 & n_Point)
{
	Vector3 a = n_Point;
	GraphTrans res;
	a.Normalization();
	if (a.x != 0 || a.y != 0 || a.z != 1)
	{
		Vector3 rotate_axis = Vector3(0, 0, 1).Cross(a); // ��ת��
		double rotate_theta = acos(Vector3(0, 0, 1).Dot(a));
		rotate_theta = rotate_theta / Pi * 180;
	
		res.updateRotate(rotate_axis, rotate_theta);
	}
	else
	{
		res.updateRotate(Vector3(0, 0, 1), 0);
	}
	return res;
}

bool PVVAAS::CheckModle(Vector3 & interPoint, Vector3 & n_Point)
{
	return false;
}

FieldBase PVVAAS::getFieldBaseout() {
	FieldBase result;
	result.ds_x = dx;	result.ds_y = dy;
	result.graphTransField = OutFieldGraphTrans;
	result.N_width = N;	result.M_depth = M;
	result.Ex.resize(N);	result.Ey.resize(N);	result.Ez.resize(N);
	result.Hx.resize(N);	result.Hy.resize(N);	result.Hz.resize(N);
	for (int i = 0; i < N; i++) {
		result.Ex[i].resize(M);	result.Ey[i].resize(M);	result.Ez[i].resize(M);
		result.Hx[i].resize(M);	result.Hy[i].resize(M);	result.Hz[i].resize(M);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			result.Ex[i][j] = Ex1[i][j];	result.Ey[i][j] = Ey1[i][j];	result.Ez[i][j] = Ez1[i][j];
			result.Hx[i][j] = Hx1[i][j];	result.Hy[i][j] = Hy1[i][j];	result.Hz[i][j] = Hz1[i][j];
		}
	}
	return result;
}

complex<double> PVVAAS::ConjugateMul(const complex<double> &A, const complex<double> &B) const
{
	return complex<double>(A.real() * B.real() + A.imag() * B.imag(),
		-A.real() * B.imag() + A.imag() * B.real());
}

complex<double> PVVAAS::ComMul(const complex<double>& A, const complex<double>& B) const
{
	return complex<double>(A.real() * B.real() - A.imag() * B.imag(),
		A.real() * B.imag() + A.imag() * B.real());
}

double PVVAAS::CalDistance(const Vector3 &a, const Vector3 &b) const
{
	return pow(pow((a.x - b.x), 2) + pow((a.y - b.y), 2)+ pow((a.z - b.z), 2), 0.5);
}

bool PVVAAS::InterVal_PointinTriangle(const Vector2 & A, 
	const Vector2 & B, const Vector2 & C, const Vector2 & P)
{
	Vector2 v0 = C - A;
	Vector2 v1 = B - A;
	Vector2 v2 = P - A;

	double dot00 = v0.Dot(v0);
	double dot01 = v0.Dot(v1);
	double dot02 = v0.Dot(v2);
	double dot11 = v1.Dot(v1);
	double dot12 = v1.Dot(v2);

	double inverDeno = 1.0 / (dot00 * dot11 - dot01 * dot01);

	double u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < -0.0001 || u > 1.00001) // if u out of range, return directly
	{
		return false;
	}

	double v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < -0.0001 || v > 1.00001) // if v out of range, return directly
	{
		return false;
	}

	return u + v <= 1.00001;
}

