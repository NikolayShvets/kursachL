#pragma once
#include "rightPatrs.h"
#include <fstream>

class integrator
{
private:
	long double t, tk, dt;
	long double eps;
	std::vector<std::vector<long double>> K;//матрица коэффициентов К
	std::vector<long double> c, b1, b2; // из таблицы Бутчера
	std::vector<std::vector<long double>> a; //из таблицы Бутчера
	long double machineZero; //машинный нуль
	std::ofstream resFile;
	std::ofstream xtFile;
public:
	std::vector<long double> PowerSpector(long double omega0, long double omega1, long double omegaStep);
	std::vector<std::vector<long double>> resMatrix;
	integrator(long double t0, long double tk, long double dt);
	void integrate(rightPatrs& rp);
	void setPrecision(long double eps); //сеттер для эпсилон
	long double max(long double a, long double b); //возвращает максимальное из а и b
	long double min(long double a, long double b);
	~integrator();
};

