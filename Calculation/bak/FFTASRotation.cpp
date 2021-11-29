#include "FFTASRotation.h"
#include "FFT.h"
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkSmartPointer.h>

FFTASRotation::FFTASRotation(double _f, int _Nu, int _Nv)
{
	freq = _f;	lambda = C_Speed / _f;	k = 2 * Pi / lambda;
	Nu = _Nu;	Nu2 = 2 * _Nu - 1;
	Nv = _Nv;	Nv2 = 2 * _Nv - 1;
	ds = lambda / 3.5;
	Allocate();
}

FFTASRotation::~FFTASRotation() {
	FreeCal();
}

void FFTASRotation::Allocate()
{
	Eu0 = Allocate_2D(Eu0, Nu2, Nv2);	Ev0 = Allocate_2D(Ev0, Nu2, Nv2);	En0 = Allocate_2D(En0, Nu2, Nv2);
	Evt = Allocate_2D(Eut, Nu2, Nv2);	Evt = Allocate_2D(Evt, Nu2, Nv2);	Ent = Allocate_2D(Ent, Nu2, Nv2);
}

void FFTASRotation::FreeCal() {
	Free_2D(Eu0);	Free_2D(Ev0);	Free_2D(En0);
	Free_2D(Eut);	Free_2D(Evt);	Free_2D(Ent);
	Eu0 = NULL;		Ev0 = NULL;		En0 = NULL;
	Eut = NULL;		Evt = NULL;		Ent = NULL;
}

void FFTASRotation::output(std::complex<double>** EuOut, std::complex<double>** EvOut)
{
	int ShiftU = (Nu - 1) / 2;
	int ShiftV = (Nv - 1) / 2;
	for (int i = 0; i < Nu; i++) {
		for (int j = 0; j < Nv; j++) {
			EuOut[i][j] = Eut[i + ShiftU][j + ShiftV];
			EvOut[i][j] = Evt[i + ShiftU][j + ShiftV];
		}
	}
}

void FFTASRotation::SetInput(std::complex<double>** EuIn, std::complex<double>** EvIn)
{
	int ShiftU = (Nu - 1) / 2;
	int ShiftV = (Nv - 1) / 2;
	for (int i = 0; i < Nu2; i++) {
		for (int j = 0; j < Nv2; j++) {
			Eu0[i][j] = 0.0;
			Ev0[i][j] = 0.0;
			En0[i][j] = 0.0;
		}
	}
	for (int i = 0; i < Nu; i++) {
		for (int j = 0; j < Nv; j++) {
			Eu0[i + ShiftU][j + ShiftV] = EuIn[i][j];
			Ev0[i + ShiftU][j + ShiftV] = EvIn[i][j];
		}
	}
}

void FFTASRotation::PerformRotate()
{
	std::vector<complex<double>> phaseU;	phaseU.resize(Nu2);
	std::vector<complex<double>> phaseV;	phaseV.resize(Nv2);
	std::vector<complex<double>> IphaseU;	IphaseU.resize(Nu2);
	std::vector<complex<double>> IphaseV;	IphaseV.resize(Nv2);
	//谱操作的相位补偿项
	for (int i = 0; i < Nu2; i++) { 
		phaseU[i] = complex<double>(cos(2 * Pi*0.5*(Nu2 - 1) / Nu2*i), sin(2 * Pi*0.5*(Nu2 - 1) / Nu2*i));
		IphaseU[i] = complex<double>(cos(-2 * Pi*0.5*(Nu2 - 1) / Nu2*i), sin(-2 * Pi*0.5*(Nu2 - 1) / Nu2*i)); 
	}
	for (int j = 0; j < Nv2; j++) {
		phaseV[j] = complex<double>(cos(2 * Pi*0.5*(Nv2 - 1) / Nv2*j), sin(2 * Pi*0.5*(Nv2 - 1) / Nv2*j));
		IphaseV[j] = complex<double>(cos(-2 * Pi*0.5*(Nv2 - 1) / Nv2*j), sin(-2 * Pi*0.5*(Nv2 - 1) / Nv2*j));
	}
	//相位补偿-空间
	for (int i = 0; i < Nu2; i++)
		for (int j = 0; j < Nv2; j++) 
		{
			Eu0[i][j] = Eu0[i][j] * IphaseU[i] * IphaseV[j];
			Ev0[i][j] = Ev0[i][j] * IphaseU[i] * IphaseV[j];
		}
	//逆傅里叶变换到平面波谱
	FFT fft;
	fft.IFFT_2D(Eu0, Eu0, Nu2, Nv2);
	fft.IFFT_2D(Ev0, Ev0, Nu2, Nv2);
	//到谱了哦
	//相位补偿-谱域
	for (int i = 0; i < Nu2; i++)
		for (int j = 0; j < Nv2; j++)
		{
			Eu0[i][j] = Eu0[i][j] * IphaseU[i] * IphaseV[j];
			Ev0[i][j] = Ev0[i][j] * IphaseU[i] * IphaseV[j];
		}
	//作为标准谱域网格坐标
	double ** cosalfaU = NULL;	double** cosbetaU = NULL;	double** cosgammaU = NULL;
	//作为非均匀谱域网格坐标
	double ** cosalfaN = NULL;	double** cosbetaN = NULL;	double** cosgammaN = NULL;
	cosalfaU = Allocate_2D(cosalfaU, Nu2, Nv2);	cosbetaU = Allocate_2D(cosbetaU, Nu2, Nv2);	cosgammaU = Allocate_2D(cosgammaU, Nu2, Nv2);
	cosalfaN = Allocate_2D(cosalfaN, Nu2, Nv2);	cosbetaN = Allocate_2D(cosbetaN, Nu2, Nv2);	cosgammaN = Allocate_2D(cosgammaN, Nu2, Nv2);

	for (int i = 0; i < Nu2; i++)
		for (int j = 0; j < Nv2; j++) {
			//均匀谱网格
			cosalfaU[i][j] = lambda*(i - 0.5*(Nu2 - 1) / ds / Nu2);
			cosbetaU[i][j] = lambda*(j - 0.5*(Nv2 - 1) / ds / Nv2);
			if (cosalfaU[i][j] * cosalfaU[i][j] + cosbetaU[i][j] * cosbetaU[i][j] < 1){
				cosgammaU[i][j] = sqrt(1 - cosalfaU[i][j] * cosalfaU[i][j] - cosbetaU[i][j] * cosbetaU[i][j]);
				if(cosgammaU[i][j] < 0.02) cosgammaU[i][j] = -2;//cosgamma过小
				else {//正常的
					En0[i][j] = (Eu0[i][j] * cosalfaU[i][j] + Ev0[i][j] * cosbetaU[i][j])
						/ (-cosgammaU[i][j]);
				}
			}
			else cosgammaU[i][j] = -2;//瞬逝波
			//非均匀谱网格
			cosalfaN[i][j] = 0;	cosbetaN[i][j] = 0;	cosgammaN[i][j] = 0;

			//极化旋转-借这个循环顺道做了
			Vector3 FieldReal;
			Vector3 FieldImag;
			FieldReal = RotationMatrix*Vector3(Eu0[i][j].real(), Ev0[i][j].real(), En0[i][j].real());
			FieldImag = RotationMatrix*Vector3(Eu0[i][j].imag(), Ev0[i][j].imag(), En0[i][j].imag());
			Eu0[i][j] = complex<double>(FieldReal.x, FieldImag.x);
			Ev0[i][j] = complex<double>(FieldReal.y, FieldImag.y);
			En0[i][j] = complex<double>(FieldReal.z, FieldImag.z);
			//极化旋转后的原谱
		}
	//重头戏，谱域插值

	//首先插出个平面三角形 非均匀网格
	vtkFloatArray *scalars = vtkFloatArray::New();
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	int num = 0;
	for (int i = 0; i < Nu2; ++i)
		for (int j = 0; j < Nv2; ++j)
		{
			if (cosgammaU[i][j] > -2) {
				Vector3 cosU(cosalfaU[i][j], cosbetaU[i][j], cosgammaU[i][j]);
				Vector3 cosN = RotationMatrix*cosU;
				//注意：这里的插值和mathcad示例里面的有所不同，
				//      这里cos0表示非均匀系，cosT表示均匀采样系
				//本质上，是将旋转前的谱坐标投影到旋转后的谱坐标上（即用非均匀网格上的数据插值得到均匀网格上的数据）
				//			 所以应该是 (kxt,kyt,kzt) = Ma(kx0,ky0,kz0)
				//        是将旋转前的谱值  投影到旋转后的谱坐标系上，将形成一堆不均匀的点
				//        这些不均匀的点用delaunay形成三角网格，然后用三角网格的坐标关系，
				//        用三角网格上的顶点值，插值均匀网格落在三角内的点上的值

				cosalfaN[i][j] = cosN.x;
				cosbetaN[i][j] = cosN.y;
				cosgammaN[i][j] = cosN.z;
				if (cosgammaN[i][j] > 0) {
					points->InsertNextPoint(cosalfaN[i][j], cosbetaN[i][j], 0);
					scalars->InsertTuple1(num++, i*Nv2 + j);//注意后面这个参量有个index的作用，在插值中会用到，也是设计的巧妙的地方
				}

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
		
	//遍历每个三角形
	for (int i = 0; i < EleNum; i++){
		double du = lambda / ds / Nu2;
		double dv = lambda / ds / Nv2;
		double maxX = lambda*(0 - 0.5*(Nu2 - 1) / ds / Nu2); 
		double minX = lambda*(Nu2 - 1 - 0.5*(Nu2 - 1) / ds / Nu2);
		double maxY = lambda*(0 - 0.5*(Nv2 - 1) / ds / Nv2);
		double minY = lambda*(Nv2 - 1 - 0.5*(Nv2 - 1) / ds / Nv2);
		vtkIdList * p;
		p = polydata->GetCell(i)->GetPointIds();

		double * point;
		vector<Vector2> A(3);//三角形存下来
		vector<int> dataIJ(3);
		Vector2 P;
		//找到三角形的每个点 并找到三角形对应的四边形
		for (int k = 0; k < 3; k++)
		{
			dataIJ[k] = polydata->GetPointData()
				->GetScalars()->GetVariantValue(p->GetId(k)).ToInt();
			point = polydata->GetPoint(p->GetId(k));

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
		int iiMin = ceil(minX / du) + (Nu2 - 1) / 2;
		int iiMax = floor(maxX / du) + (Nu2 - 1) / 2;
		int jjMin = ceil(minY / dv) + (Nv2 - 1) / 2;
		int jjMax = floor(maxY / dv) + (Nv2 - 1) / 2;
		if (iiMin < 0)
			iiMin = 0;
		if (jjMin < 0)
			jjMin = 0;
		for (int ii = iiMin; ii <= iiMax&&ii<Nu2; ii++)
			for (int jj = jjMin; jj <= jjMax&&jj<Nv2; jj++)
			{
				P.set(cosalfaU[ii][jj], cosbetaU[ii][jj]);//网格位置点
				if (InterVal_PointinTriangle(A[0], A[1], A[2], P))//判断网格位置点在不在三角形里
				{//在三角形里，插值
					double a = (P - A[0]).area();
					double b = (P - A[1]).area();
					double c = (P - A[2]).area();

					Eut[ii][jj] = (1.0/	(b * c + a * c + b * a))*
					   (b * c * Eu0[dataIJ[0] / Nv2][dataIJ[0] % Nv2] +
						a * c * Eu0[dataIJ[1] / Nv2][dataIJ[1] % Nv2] +
						b * a * Eu0[dataIJ[2] / Nv2][dataIJ[2] % Nv2]);
					Evt[ii][jj] = (1.0 / (b * c + a * c + b * a))*
					   (b * c * Ev0[dataIJ[0] / Nv2][dataIJ[0] % Nv2] +
						a * c * Ev0[dataIJ[1] / Nv2][dataIJ[1] % Nv2] +
						b * a * Ev0[dataIJ[2] / Nv2][dataIJ[2] % Nv2]);
					Ent[ii][jj] = (1.0 / (b * c + a * c + b * a))*
						(b * c * En0[dataIJ[0] / Nv2][dataIJ[0] % Nv2] +
						 a * c * En0[dataIJ[1] / Nv2][dataIJ[1] % Nv2] +
						 b * a * En0[dataIJ[2] / Nv2][dataIJ[2] % Nv2]);
				}
			}
	}
	//插值结束
	//旋转后的谱已准备好
	//变换回场前相位补偿-对应FFT-还有乘上雅克比行列式
	for (int i = 0; i < Nu2; i++)
		for (int j = 0; j < Nv2; j++) {
			Eut[i][j] = Eut[i][j] * phaseU[i] * phaseV[j];
			Evt[i][j] = Evt[i][j] * phaseU[i] * phaseV[j];
			Ent[i][j] = Ent[i][j] * phaseU[i] * phaseV[j];
			if (cosgammaU[i][j]>0.02) {
				Eut[i][j] = Eut[i][j] * abs(cosgammaU[i][j] / cosgammaN[i][j]);
				Evt[i][j] = Evt[i][j] * abs(cosgammaU[i][j] / cosgammaN[i][j]);
				Ent[i][j] = Ent[i][j] * abs(cosgammaU[i][j] / cosgammaN[i][j]);
			}
		}
	fft.FFT_2D(Eut, Eut, Nu2, Nv2);
	fft.FFT_2D(Evt, Evt, Nu2, Nv2);
	//变换回场后相位补偿
	for (int i = 0; i < Nu2; i++) 
		for (int j = 0; j < Nv2; j++) {
			Eut[i][j] = Eut[i][j] * phaseU[i] * phaseV[j];
			Evt[i][j] = Evt[i][j] * phaseU[i] * phaseV[j];
			Ent[i][j] = Ent[i][j] * phaseU[i] * phaseV[j];
		}
	//完成计算


	Free_2D(cosalfaU);	Free_2D(cosbetaU);	Free_2D(cosgammaU);
	Free_2D(cosalfaN);	Free_2D(cosbetaN);	Free_2D(cosgammaN);
	cosalfaU = NULL;	cosbetaU = NULL;	cosgammaU = NULL;
	cosalfaN = NULL;	cosbetaN = NULL;	cosgammaN = NULL;

}

void FFTASRotation::SetParas(double _f, double _Nu, double _Nv, double _ds) {
	freq = _f;
	Nu = _Nu;	Nu2 = Nu * 2 - 1;
	Nv = _Nv;	Nv2 = Nv * 2 - 1;
	ds = _ds;
	FreeCal();
	Allocate();
}

void FFTASRotation::SetRotationParasinRad(double _theta_rou, double _phi_rou) {
	theta_rou = _theta_rou;
	phi_rou = _phi_rou;
	RotationMatrix = Matrix4D::getRotateMatrixRad(theta_rou, phi_rou);
	InvRotationMatrix = Matrix4D::getInvRotateMatrixRad(theta_rou, phi_rou);
}

bool FFTASRotation::InterVal_PointinTriangle(const Vector2 & A,
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
