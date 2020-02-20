#include "pch.h"
#include "integrator.h"

integrator::integrator(long double t0, long double tk, long double dt)
{
	this->resFile.open("result.txt");
	this->xtFile.open("xtresult.txt");
	//начало, конец, шаг интегрирования
	this->t = t0;
	this->tk = tk;
	this->dt = dt;
	//матрица промежуточных вычислений коэффициентов методов Рунге-Кутты
	this->K.resize(7);
	//коэффициенты из таблицы Бутчера
	this->b1.resize(7);
	this->b2.resize(7);
	this->c.resize(7);
	this->a.resize(7);

	for (auto& elem : a) {
		elem.resize(6);
	}
	//коэффициенты из таблицы Бутчера
	b1 = { 35. / 384, 0., 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0 };
	b2 = { 5179. / 57600, 0., 7571. / 16695, 393. / 640, -92097. / 339200, 187. / 2100, 1. / 40 };
	c = { 0, 3.0, 3. / 10, 4. / 5, 8. / 9, 1., 1. };
	//непосредственно матрица коэффициентов Бутчера
	a = {
	{ 0. },
	{ 1. / 5 },
	{ 3. / 40, 9. / 40 },
	{ 44. / 45, -56. / 15, 32. / 9 },
	{ 19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729 },
	{ 9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656 },
	{ 35. / 384, 0., 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84 }
	};
	//вычисление машинного нуля
	long double temp{ 0.1 };
	while (this->machineZero + 1.0 > 1) {
		this->machineZero = temp;
		temp /= 2.0;
	}
}

std::vector<long double> integrator::PowerSpector(long double omega0, long double omega1, long double omegaStep)
{
	//const long double M_PI = 3.14159265358979323846264;
	long double ts{ -2 * M_PI }, tst{ 0.001 }, tfin{ 2*M_PI };
	resMatrix.clear();
	int i{ 0 };
	while (ts <= tfin) {
		resMatrix.resize(i + 1);
		resMatrix[i].push_back(1*sinl(ts));
		resMatrix[i].push_back(ts);
		ts += tst;
		++i;
	}
	for (const auto & elem : resMatrix) {
		std::cout << elem[0] << " " << elem[1] << std::endl;
	}
	
	std::vector<long double> result;
	long double A{ 0.0 };
	long double B{ 0.0 };
	/*long double leftT{0.0}, rightT{0.0};
	
	leftT = this->resMatrix[0][1];
	rightT = this->resMatrix[1][1];

	std::vector<std::vector<long double>> resMatrixNorm;
	for (int i = 1; i < resMatrix.size(); ++i) {
		if ( (abs(rightT - leftT - 1.0)) < 0.1) {
			resMatrixNorm.push_back(resMatrix[i]);
			leftT = resMatrix[i-1][1];
		}
		rightT = resMatrix[i][1];
	}

	for (const auto& elem : resMatrixNorm) {
		std::cout << elem[0] << std::endl;
	}*/

	while (omega0 <= omega1) {
		for (int i = 1; i < resMatrix.size(); ++i) {
		//for (int i = resMatrix.size()-1; i > 0; --i){
			A += 0.5L * (resMatrix[i - 1][0] * cosl(omega0 * resMatrix[i - 1][1]) +
				resMatrix[i][0] * cosl(omega0 * resMatrix[i][1])) * (resMatrix[i][1] - resMatrix[i - 1][1]);

			B += 0.5L * (resMatrix[i - 1][0] * sinl(omega0 * resMatrix[i - 1][1]) +
				resMatrix[i][0] * sinl(omega0 * resMatrix[i][1])) * (resMatrix[i][1] - resMatrix[i - 1][1]);
		}
		result.push_back(powl(A, 2.L) + powl(B, 2.L));
		xtFile << result.back() << "|" << omega0 << std::endl;
		A = B = 0.0L;
		omega0 += omegaStep;
	}


	/*while (omega0 <= omega1) {
		for (int i = 1; i < this->resMatrix.size(); ++i) {
			if ((resMatrix[i][1] - resMatrix[i - 1][1] - 1.0) < 0.0000000001) {
				flag = true;
				A += 0.5 * (resMatrix[i - 1][0] * cos(omega0 * resMatrix[i - 1][1]) +
					resMatrix[i][0] * cos(omega0 * resMatrix[i][1])) * (resMatrix[i][1] - resMatrix[i - 1][1]);

				B += 0.5 * (resMatrix[i - 1][0] * sin(omega0 * resMatrix[i - 1][1]) +
					resMatrix[i][0] * sin(omega0 * resMatrix[i][1])) * (resMatrix[i][1] - resMatrix[i - 1][1]);
				result.push_back(pow(A, 2.) + pow(B, 2.));
			}
			else {
				flag = false;
				continue;
			}
		}
		if(flag) omega0 += omegaStep;
	}*/
	//long double temp = 0.0;
	/*for (const auto& r : result) {
		std::cout << r << "|" << std::endl;
		xtFile << r << "|" <<temp<< std::endl;
		temp += omegaStep;
	}*/
	return result;
}

/*void integrator::getMZ() {
	long double temp{ 0.1 };
	while (this->machineZero + 1.0 > 1) {
		this->machineZero = temp;
		temp /= 2.0;
	}
}*/

void integrator::setPrecision(long double eps) {
	this->eps = eps;
}

void integrator::integrate(rightPatrs& rp) {
	std::vector<long double>
		X = rp.X, //вектор начальных условий
		Y(X.size()), //вектор правых частей
		X4(X.size()),//решение 4 порядка
		X5(X.size()), //решение 5 порядка
		Xres(X.size()); //результат после плотной выдачи
	long double newDt = dt, e = 0.0L, step, newStep = dt, T = t;

	for (auto& k : K) { //заводим по 7 коэффициентов К для каждого уравнения в системе
		k.resize(X.size());
	}
	int i{0};
	while (this->t < this->tk) {//пока не вышло время интегрирования
		step = newStep;
	//пробегаем в цикле : 
		for (int j = 0; j < K.size(); ++j) { //матрицу 7хn К по строкам
			for (int k = 0; k < X.size(); ++k) {//каждое уравнение в системе
				Y[k] = X[k]; //забираем начальное условие
				for (int s = 0; s < j; ++s) { 
					Y[k] += K[s][k] * a[j][s] * step; //собираем правую часть уравнения К = F(X + h*a*K)
				}
				//на выходе - вектор правых частей уравнений для К
			}
			//закидываем это в RightPArts, результат вычисления правых частей по вектору X (Y = X 68стр) присваивается K[j]
			rp.getRP(Y, K[j]);//правые части вызываеются 7 раз, пеесчитывая коэффициенты К
		}
	
		e = 0.0;

		for (int i = 0; i < X.size(); ++i){
			X4[i] = X[i];
			X5[i] = X[i];
			for (int j = 0; j < b1.size(); ++j) {
				X4[i] += K[j][i] * b1[j] * step;
				X5[i] += K[j][i] * b2[j] * step;
			}
			e += powl(step * (X4[i] - X5[i]) / 
				max(fmaxl(abs(X[i]), fabsl(X4[i])), max(1e-5L, 2 * machineZero / eps)), 2.0);
		}
		e = sqrtl(e / X.size());
		newStep = step / max(0.1, min(5., powl(e / eps, 0.2) / 0.9));

		if (e > eps) continue;

		while ( (T < t + step) && (T <= tk) ){
			long double theta = (T - t) / step;
			std::vector<long double> b(6);
			b[0] = theta * (1 + theta * (-1337. / 480. + theta * (1039. / 360. + theta * (-1163. / 1152.))));
			b[1] = 0.0;
			b[2] = 100. * powl(theta, 2.) * (1054. / 9275. + theta * (-4682. / 27825. + theta * (379. / 5565.))) / 3.;
			b[3] = -5. * powl(theta, 2.) * (27. / 40. + theta * (-9. / 5. + theta * (83. / 96.))) / 2.;
			b[4] = 18225. * powl(theta, 2.) * (-3. / 250. + theta * (22. / 375. + theta * (-37. / 600.))) / 848.;
			b[5] = -22. * powl(theta, 2.) * (-3. / 10. + theta * (29. / 30. + theta * (-17. / 24.))) / 7.;

			for (int i = 0; i < X.size(); ++i)
			{
				long double sum = 0;
				for (int j = 0; j < b.size(); ++j) {
					sum += b[j] * K[j][i];
				}
				Xres[i] = X[i] + step * sum;
			}
			T += dt;
		}

		for (const auto& x :Xres) {
			std::cout << x << " ";
			resFile << x << "|";
		}
		resMatrix.resize(i + 1);
		resMatrix[i].push_back(Xres.front());
		resMatrix[i].push_back(t);
		std::cout << t <<std::endl;
		resFile << t << std::endl;
		X = X4;
		t += step;
		++i;
	}
}

long double integrator::max(long double a, long double b) {
	if (a > b) return a;
	else return b;
}
long double integrator::min(long double a, long double b) {
	if (a < b) return a;
	else return b;
}

integrator::~integrator()
{
}
