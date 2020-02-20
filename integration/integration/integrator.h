#pragma once
#include "rightPatrs.h"
#include <fstream>

class integrator
{
private:
	long double t, tk, dt;
	long double eps;
	std::vector<std::vector<long double>> K;//������� ������������� �
	std::vector<long double> c, b1, b2; // �� ������� �������
	std::vector<std::vector<long double>> a; //�� ������� �������
	long double machineZero; //�������� ����
	std::ofstream resFile;
	std::ofstream xtFile;
public:
	std::vector<long double> PowerSpector(long double omega0, long double omega1, long double omegaStep);
	std::vector<std::vector<long double>> resMatrix;
	integrator(long double t0, long double tk, long double dt);
	void integrate(rightPatrs& rp);
	void setPrecision(long double eps); //������ ��� �������
	long double max(long double a, long double b); //���������� ������������ �� � � b
	long double min(long double a, long double b);
	~integrator();
};

