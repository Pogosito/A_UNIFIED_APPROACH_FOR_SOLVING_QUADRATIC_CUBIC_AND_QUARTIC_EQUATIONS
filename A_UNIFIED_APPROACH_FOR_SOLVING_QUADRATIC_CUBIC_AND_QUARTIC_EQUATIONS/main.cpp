//
//  main.cpp
//  A_UNIFIED_APPROACH_FOR_SOLVING_QUADRATIC_CUBIC_AND_QUARTIC_EQUATIONS
//
//  Created by Pogos Anesyan on 09.04.2023.
//

/*
	Реализация метода вычисления корней полинома четвертой степени с помощью метода "A unified approach for solving quadratic cubic and quartic equations".
	Автор: A.A.UNGAR
	Ссылка на статью: https://core.ac.uk/download/pdf/82351384.pdf
	Описание опечаток стати можно найти на: https://disk.yandex.ru/i/otdqTx2_2Ax_lw
 */

#define PRINT true
#include <iostream>
#include <vector>
#include "excerpt.h"
#include <limits>
#include <cmath>

using namespace  std;

// MARK: - fms/fmac/fmsc

template<typename fp_t>
inline fp_t fms(fp_t a, fp_t b, fp_t c, fp_t d)
{
	fp_t cd = -c * d;

	return fma(a, b, cd) - fma(c, d, cd);
}

template<typename fp_t>
complex<fp_t> fmac(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c){
	fp_t ar, ai, br, bi, cr, ci, p11, p1, p21, p2;
	ar = a.real(); ai = a.imag();
	br = b.real(); bi = b.imag();
	cr = c.real(); ci = c.imag();
	p11 = fma(ar, br, cr);
	p1 = fma(-ai, bi, p11);
	p21 = fma(ai, br, ci);
	p2 = fma(ar, bi, p21);
	complex<fp_t> num(p1, p2);
	return num;
}

template<typename fp_t>
complex<fp_t> fmsc(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d)
{
	complex<fp_t> cd = -c * d;

	return fmac(a, b, cd) - fmac(c, d, cd);
}

// MARK: - Helper methods

/// Выясняет является ли переданное число `истинно` комплексным
/// - Parameter x: Число для которого выясняем является ли оно комплексным.
template<typename fp_t>
inline bool isComplex(const complex<fp_t>& x)
{
	return abs(x) * numeric_limits<fp_t>::epsilon() <= abs(x.imag());
}

/// Кубический корень от комплексного числа
/// - Parameter number: Число, у которого выделяем кубический корень
template<typename fp_t>
vector<complex<fp_t>> cubeRoot(complex<fp_t> number)
{
	vector<complex<fp_t>> result;
	static const fp_t _2PI = static_cast<fp_t>(numbers::pi) * static_cast<fp_t>(2.0L);
	static const fp_t ONE_THIRDS = static_cast<fp_t>(1.0L / 3.0L);

	fp_t argValue = arg(number);
	fp_t modulus = abs(number);

	for (int k = 0; k < 3; ++k) {

		fp_t realOfCube = cos(fma(_2PI, k, argValue) * ONE_THIRDS);
		fp_t imagOfCube = sin(fma(_2PI, k, argValue) * ONE_THIRDS);

		if (abs(imagOfCube) < std::numeric_limits<fp_t>::epsilon()) { imagOfCube = 0; }

		complex<fp_t> root_k = pow(modulus, ONE_THIRDS) * complex<fp_t>(realOfCube, imagOfCube);

		result.push_back(root_k);
	}

	return result;
}

/// Метод для нахождения индексов у массивов кубических корней от `betta0, gamma0`, таких что выполнилось условие 4.10
/// - Parameters:
///   - betta0CubeRoots: Кубические корни betta0
///   - gamma0CubeRoots: Кубические корни gamma0
///   - Q: Расчетный коэффициент Q
///   - beta0_I: Индекс корня от `beta0`
///   - gamma0_J: Индекс корня от `gamma0`
template<typename fp_t>
void check4_10Condition(vector<complex<fp_t>> betta0CubeRoots,
						vector<complex<fp_t>> gamma0CubeRoots,
						fp_t Q, int &beta0_I, int &gamma0_J) {
	fp_t min_r = numeric_limits<fp_t>::infinity();
	fp_t min_rr = numeric_limits<fp_t>::infinity();
	
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			complex<fp_t> helper = betta0CubeRoots[i] * gamma0CubeRoots[j];
			if (abs(helper.imag()) <= min_rr)
			{
				min_rr = abs(helper.imag());
				if (abs(helper.real() - Q) <= min_r)
				{
					min_r = abs(helper.real() - Q);
					beta0_I = i;
					gamma0_J = j;
					break;
				}
			}
		}
	}
}

/// Метод для нахождения индексов у массивов корней от `betta, gamma delta`, таких что выполнилось условие 4.11
/// - Parameters:
///   - bettaSqrts: Корни от `betta`
///   - gammaSqrts: Корни от `gamma`
///   - deltaSqrts: Корни от `delta`
///   - P: Расчетный коэффициент P
///   - bettaI: Индекс корня от `betta`
///   - gammaJ: Индекс корня от `gamma`
///   - deltaK: Индекс корня от `delta`
template<typename fp_t>
void check4_11Condition(vector<complex<fp_t>> bettaSqrts,
						vector<complex<fp_t>> gammaSqrts,
						vector<complex<fp_t>> deltaSqrts,
						fp_t P, int& bettaI, int& gammaJ,
						int& deltaK)
{

	fp_t min_r1 = numeric_limits<fp_t>::infinity();
	fp_t min_r2 = numeric_limits<fp_t>::infinity();

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				complex<fp_t> bettaGammaDeltaSqrts = complex<fp_t>(bettaSqrts[i] * gammaSqrts[j] * deltaSqrts[k]);
				if (abs(bettaGammaDeltaSqrts.imag()) <= min_r2 && abs(bettaGammaDeltaSqrts.real() + P) <= min_r1) {
					min_r2 = abs(bettaGammaDeltaSqrts.imag());
					min_r1 = abs(bettaGammaDeltaSqrts.real() + P);
					bettaI = i;
					gammaJ = j;
					deltaK = k;
					break;
				}
			}
		}
	}
}

// MARK: - A UNIFIED APPROACH FOR SOLVING QUADRATIC CUBIC AND QUARTIC EQUATIONS

template<typename fp_t>
unsigned int unifiedApproachForQuartic(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{

	unsigned int numberOfRoots = 0;

	// MARK: - Объявление констант

	static const fp_t EIGHT_THIRDS = static_cast<fp_t>(8.0L / 3.0L);
	static const fp_t ONE_FOURTH = static_cast<fp_t>(0.25L);

	static const complex<fp_t> COMPLEX_ONE_FOURTH = complex<fp_t>(ONE_FOURTH, static_cast<fp_t>(0.0L));
	static const complex<fp_t> COMPLEX_FOUR_THIRDS = complex<fp_t>(static_cast<fp_t>(4.0L / 3.0L), static_cast<fp_t>(0.0L));
	static const complex<fp_t> COMPLEX_HALF = complex<fp_t>(static_cast<fp_t>(0.5L), static_cast<fp_t>(0.0L));

	static const complex<fp_t> q0 = complex<fp_t>(static_cast<fp_t>(1.0L), static_cast<fp_t>(0.0L));
	static const complex<fp_t> q1 = complex<fp_t>(static_cast<fp_t>(-0.5L), static_cast<fp_t>(0.5L) * sqrt(static_cast<fp_t>(3.0L)));
	static const complex<fp_t> q2 = complex<fp_t>(static_cast<fp_t>(-0.5L), static_cast<fp_t>(-0.5L) * sqrt(static_cast<fp_t>(3.0L)));

	// MARK: - Вычисляем расчетные коэффициенты P, Q, R, alpha0, betta0, gamma0

	// P
	fp_t P_helper = fms(static_cast<fp_t>(4.0L), -b, -a, a);
	fp_t P = fms(static_cast<fp_t>(8.0L), c, -a, P_helper);

	// Q
	fp_t Q_helper = fms(static_cast<fp_t>(4.0L), d, a, c);
	fp_t Q = fms(static_cast<fp_t>(3.0L), Q_helper, -b, b);

	// R
	fp_t R_helper0 = fms(static_cast<fp_t>(3.0L) * d, a, b, c);
	fp_t R_helper1 = fms(R_helper0, a, c * c, -static_cast<fp_t>(3.0L));
	fp_t R_helper2 = fms(b, b, static_cast<fp_t>(36.0L), d);
	fp_t R = fms(static_cast<fp_t>(2.0L) * R_helper2, b, -R_helper1, static_cast<fp_t>(9.0L));

	fp_t alpha0 = fms(a, a, EIGHT_THIRDS,b);

	fp_t betta_gamma_helper = fma(ONE_FOURTH, R * R, -pow(Q, static_cast<fp_t>(3.0L)));

	complex<fp_t> sqrt_betta_gamma_helper = (betta_gamma_helper > 0) ? sqrt(betta_gamma_helper) : complex<fp_t>(0, sqrt(abs(betta_gamma_helper)));
	complex<fp_t> betta_helper_1 = fmac(COMPLEX_HALF, complex<fp_t>(R, 0), sqrt_betta_gamma_helper);
	complex<fp_t> gamma_helper_1 = fmac(COMPLEX_HALF, complex<fp_t>(R, 0), -sqrt_betta_gamma_helper);

	// MARK: Находим все значения кубического корня (R + sqrt(R^2 - 4Q^3)) / 2. В статье ошибка (см. файл с опечатками)
	vector<complex<fp_t>> betta0CubeRoots = cubeRoot(betta_helper_1);
	vector<complex<fp_t>> gamma0CubeRoots = cubeRoot(gamma_helper_1);

	// MARK: - Согласно теореме 5.1 считаем количество комплексных и действительных корней.
	// 1) R*R - 4*pow(Q,3) > 0  --- > 2 действит. + 2 комплексных корня
	// 2) R*R - 4*pow(Q,3) = 0  --- > 2 действит. + (2 комплексных <=> T >= 0 для всех cuberoots(betta_helper_1))
	// 3) R*R - 4*pow(Q,3) < 0  --- > 3.1) 4 действит. <=> T >= 0 для всех cuberoots(betta_helper_1)
	int cnt_negative_T = 0;
	int cnt_pozitive_T = 0;

	std::vector<fp_t> T;

	for (int i = 0; i < 3; ++i)
	{
		fp_t _3a1_8b = fms(static_cast<fp_t>(3.0L) * a, a, static_cast<fp_t>(8.0L), b);
		fp_t t = fma(static_cast<fp_t>(8.0L), real(betta0CubeRoots[i]), _3a1_8b);
		// MARK: Используем round т.к никогда не получаем чистый 0. (Проверка с eps не проходит. При вычитании двух одинаковых чисел получали не 0, но при этом результат больше eps)
		if (round(abs(static_cast<fp_t>(8.0L) * real(betta0CubeRoots[i]) + _3a1_8b)) <= numeric_limits<fp_t>::epsilon()) { t = 0; }
		T.push_back(t);
		if (T[i] < 0) cnt_negative_T++;
		if (T[i] >= 0) cnt_pozitive_T++;
	}

	// MARK: Используем round т.к никогда не получаем чистый 0. (Проверка с eps не проходит)
	if (round(abs(betta_gamma_helper)) <= numeric_limits<fp_t>::epsilon() && cnt_pozitive_T == 3)
	{
		numberOfRoots = 4;
	} else if (betta_gamma_helper > 0) {
		numberOfRoots = 2;
	} else if (betta_gamma_helper < 0) {
		numberOfRoots = cnt_pozitive_T == 3 ? 4 : 0;
	}

	// Провекра условия 4.10
	int beta0_I = 0;
	int gamma0_J = 0;
	check4_10Condition(betta0CubeRoots, gamma0CubeRoots, Q, beta0_I, gamma0_J);

	complex<fp_t> alpha0_complex = complex<fp_t>(alpha0, static_cast<fp_t>(0.0L));
	vector<complex<fp_t>> bettaSqrts = {sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q0, betta0CubeRoots[beta0_I], -q0, gamma0CubeRoots[gamma0_J]), alpha0_complex) ),
										-sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q0, betta0CubeRoots[beta0_I], -q0, gamma0CubeRoots[gamma0_J]), alpha0_complex) )};

	vector<complex<fp_t>> gammaSqrts = {sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q1, betta0CubeRoots[beta0_I], -q2, gamma0CubeRoots[gamma0_J]), alpha0_complex)),
										-sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q1, betta0CubeRoots[beta0_I], -q2, gamma0CubeRoots[gamma0_J]), alpha0_complex))};

	vector<complex<fp_t>> deltaSqrts = {sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q2, betta0CubeRoots[beta0_I], -q1, gamma0CubeRoots[gamma0_J]), alpha0_complex)),
										-sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q2, betta0CubeRoots[beta0_I], -q1, gamma0CubeRoots[gamma0_J]), alpha0_complex))};

	// Провекра условия 4.11
	int bettaI = 0;
	int gammaJ = 0;
	int deltaK = 0;
	check4_11Condition(bettaSqrts, gammaSqrts, deltaSqrts, P, bettaI, gammaJ, deltaK);

	// MARK: - Находим корни

	fp_t alpha = -ONE_FOURTH * a;
	complex<fp_t> alpha_complex = complex<fp_t>(alpha, static_cast<fp_t>(0.0L));

	complex<fp_t> w0 = fmac(COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));
	complex<fp_t> w1 = fmac(-COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(-COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));
	complex<fp_t> w2 = fmac(-COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(-COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));
	complex<fp_t> w3 = fmac(COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(-COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(-COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));

	if (numberOfRoots == 4)
	{
		roots[0] = w0.real();
		roots[1] = w1.real();
		roots[2] = w2.real();
		roots[3] = w3.real();
	} else if (numberOfRoots == 2) {
		if (!isComplex(w0)) roots[0] = w0.real();
		if (!isComplex(w1)) roots[1] = w1.real();
		if (!isComplex(w2)) roots[2] = w2.real();
		if (!isComplex(w3)) roots[3] = w3.real();
	}
	return numberOfRoots;
}

// MARK: - Метод для теста

template<typename fp_t>
void testQuarticPolynomial(int testCount, long double maxDistance)
{
	unsigned P = 4; // Степень исходного полинома
	fp_t low = -1, high = 1; // Интервал на котором заданы корни полинома
	fp_t absMaxError, relMaxError; // Абсолютная и относительная погрешность по итогам пройденного теста
	fp_t absMaxErrorTotal = -1, relMaxErrorTotal = -1; // Итоговая максимальная абсолютная и относительная погрешность по итогам всех тестов
	long double absErrorAvg = 0, relErrorAvg = 0; // Средняя абсолютная и относительная погрешность по итогам всех тестов
	unsigned numberOfFoundRoots; // Количество найденных корней
	unsigned cantFind = 0; // Счетчик количества ситуаций, когда методу не удалось найти корни (numberOfFoundRoots == 0)
	vector<fp_t> coefficients(P + 1); // Вектор коэффициентов полинома
	unsigned count = 0; // Счетчик количества ситуаций, когда относительная погрешность больше определенного числа (relMaxError > n)
	int countExcessRoots = 0;
	int countLostRoots = 0;

	for (size_t i = 0; i < testCount; ++i)
	{
		vector<fp_t> foundRoots(P);
		vector<fp_t> trueRoots(P);
		int excessRoots = 0;
		int lostRoots = 0;

		generate_polynomial<fp_t>(P, 0, 0, 0, static_cast<fp_t>(maxDistance), low, high, trueRoots, coefficients);

		numberOfFoundRoots = unifiedApproachForQuartic<fp_t>(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

		if (numberOfFoundRoots > 0)
		{
			compare_roots<fp_t>(numberOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError, excessRoots, lostRoots);

			absMaxErrorTotal = absMaxError > absMaxErrorTotal ? absMaxError : absMaxErrorTotal;
			absErrorAvg += absMaxError;

			relMaxErrorTotal = relMaxError > relMaxErrorTotal ? relMaxError : relMaxErrorTotal;
			relErrorAvg += relMaxError;

			countExcessRoots += excessRoots;
			countLostRoots += lostRoots;

			count += relMaxError > 1 ? 1 : 0;
		} else {
			countLostRoots += 4;
			cantFind += 1;
		}
	}

	absErrorAvg /= (testCount - cantFind);
	relErrorAvg /= (testCount - cantFind);

	if (PRINT)
	{
		cout << "QUARTIC TEST RESULTS" << endl;
		cout << "========================================" << endl;
		cout << "Max distance: " << maxDistance << endl;
		cout << "Total count of tests: " << testCount << endl;
		cout << "Couldn't find roots: " << cantFind << " times " << endl;
		cout << "----------------------------------------" << endl;
		cout << "Average absolute error: " << absErrorAvg << endl;
		cout << "Total maximum absolute error: " << absMaxErrorTotal << endl;
		cout << "Average relative error: " << relErrorAvg << endl;
		cout << "Total maximum relative error: " << relMaxErrorTotal << endl;
		cout << "----------------------------------------" << endl;
		cout << "Total count of lost roots: " << countLostRoots << endl;
		cout << "Total count of excess roots: " << countExcessRoots << endl;
		cout << "----------------------------------------" << endl;
		cout << "relMaxError > 1: " << count << " times" << endl;
		cout << "========================================" << endl;
	}
}

int main() {
	testQuarticPolynomial<float>(1'000'000, 1e-5);
}
