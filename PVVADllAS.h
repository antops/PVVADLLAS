#pragma once
 //���ýӿ�
#ifndef PPVADLLAS_H
#define PPVADLLAS_H

#include <string>
#include "../Util/FieldBase.h"

class Mirror;
class _declspec(dllexport) PVVADllAS
{
public:
	PVVADllAS();
	~PVVADllAS();

	void setInField(FieldBase * _fin);
	void setOutField(FieldBase * _fout);
	void setAnalyticalModelFile(std::string _mirrorfile);
	void setSTLModelFile(std::string _stlfile);
	FieldBase getFieldout();
	void setSingleTilt(bool _in);				//���λ��Ƕ����б
	void setThreadNum(int _threadNum) { ThreadNum = _threadNum; }
	void calculate(double fre);	
	void setCATRMode(bool in);
	void SetReturnFloat(void(*returnFloat)(float, void*), void*_user);// ע��ص�����


private:
	FieldBase FieldBaseIn;
	FieldBase FieldBaseOut;
	Mirror *mirrorptr;

	std::string inputFieldFile;
	std::string stlMirrorFile;
	bool oneside;
	void(*returnFloat)(float, void*);
	void *user; // �ص���������ָ��'
	bool AnalyticalMirror;
	bool CATRMode;
	int ThreadNum;
};


#endif // PPVADLL_H