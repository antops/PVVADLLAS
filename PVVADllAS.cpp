#include "PVVADllAS.h"
#include "Calculation\PVVAAS.h"
#include "VTK\Field.h"
#include "../FFTASRotationDLL/FFTASRotation.h"
#include "../MirrorsAndRestriction/Mirror.h"
#include "../MirrorsAndRestriction/MirrorFactory.h"
#include "../MirrorsAndRestriction/Restriction.h"
#include "../MirrorsAndRestriction/Paraboloid.h"
#include "../MirrorsAndRestriction/Ellipsoid.h"
#include "../MirrorsAndRestriction/STLMirror.h"
#include "../Util/FieldBase.h"

PVVADllAS::PVVADllAS()
{
	this->returnFloat = NULL;
	this->user = NULL;
	CATRMode = false;
	ThreadNum = 4;
}

PVVADllAS::~PVVADllAS()
{
}

void PVVADllAS::setInField(FieldBase * _in) {
	FieldBaseIn.Ex = _in->Ex;
	FieldBaseIn.Ey = _in->Ey;
	FieldBaseIn.N_width = _in->N_width;
	cout << FieldBaseIn.N_width << endl;
	FieldBaseIn.M_depth = _in->M_depth;
	FieldBaseIn.graphTransField = _in->graphTransField;
	FieldBaseIn.ds_x = _in->ds_x;
	FieldBaseIn.ds_y = _in->ds_y;
}

void PVVADllAS::setOutField(FieldBase * _out) {
	FieldBaseOut.graphTransField = _out->graphTransField;
}

void PVVADllAS::setAnalyticalModelFile(std::string _mirrorfile) {
	Json::Reader reader;
	Json::Value js;
	ifstream file(_mirrorfile);
	reader.parse(file, js);
	file.close();
	mirrorptr = MirrorFactory::getMirrorByJson(js["Mirror"][0], ".");
	mirrorptr->updateData();
	AnalyticalMirror = true;
}

void PVVADllAS::setSTLModelFile(std::string _stlfile) {
	mirrorptr = new STLMirror;
	mirrorptr->setNameFile(_stlfile);
}

void PVVADllAS::setCATRMode(bool in) {
	CATRMode = in;
}

void PVVADllAS::setSingleTilt(bool _in)
{
	oneside = _in;
}

void PVVADllAS::calculate(double fre)
{
	Field inputField;
	inputField.setFileAddress(inputFieldFile);
	inputField.readDataBinary();

	Mirror* tempmirror = (Mirror*) mirrorptr;

	//STLMirror stlMirror;
	//stlMirror.setNameFile(stlMirrorFile);
	//stlMirror.readData();

	
	//�����롢�������ģ�ʹ��ݸ�PVVA
	PVVAAS pvva;
	// ���õ�λ
	pvva.SetReturnFloat(returnFloat, user);
	pvva.setUnit(1);
	// ����Ƶ��
	pvva.setFre(fre);

	pvva.setOneSide(false);

	pvva.SetThreadNum(ThreadNum);
	
	// ����Դ�������ڴ�
	pvva.setSource(&FieldBaseIn);//���볡

	pvva.setMirror(mirrorptr);//����

	pvva.setOut(&FieldBaseOut);//�����Ķ���
	
	if (AnalyticalMirror) pvva.CalZ0Theta();
	else pvva.CalMirrorCenter();

	pvva.ReflectCUDA();//������������1���䴫����������2�ϣ�����ֵ

	pvva.InterVal(); //��ֵ�õ�������2�Ͼ��ȷֲ��ĳ�

	pvva.Result();//������2�Ͼ��ȷֲ��ĳ�������Ŀ�����ϣ����

	FieldBaseOut = pvva.getFieldBaseout();
	//pvva.getField(&inputField);
}

FieldBase PVVADllAS::getFieldout() {
	return FieldBaseOut;
}

void PVVADllAS::SetReturnFloat(void(*returnFloat)(float, void *), void * _user)
{
	this->returnFloat = returnFloat;
	this->user = _user;
}
