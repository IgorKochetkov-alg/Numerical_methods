#define _CRT_CESURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define A -3.4
#define B 1.1
#define EPS 1e-13
#define FVALUE 69.93014887001885


double Func(double x)
{
	return (pow(x, 5) - 7.2 * pow(x, 3) + 8.5 * pow(x, 2) - 7 * x - 4.5) * cos(0.4 * x);
}

double IntFunc(double x)
{
	return pow(x, 6) / 6 - 7.2 * pow(x, 4) / 4 + 8.5 * pow(x, 3) / 3 - 7 * x * x / 2 - 4.5 * x;
}

double Simpson(double a, double b, int n)
{
	int newN = 2 * n;
	double h = (b - a) / newN;
	double sum1 = 0;
	double sum2 = 0;
	double result = 0;

	for (int k = 1; k <= n; k++)
	{
		sum1 += Func(a + h * (2 * k - 1));
	}
	
	for (int k = 1; k < n; k++)
	{
		sum2 += Func(a + h * (2 * k));
	}

	result = (h / 3) * (Func(a) + Func(b) + 4 * sum1 + 2 * sum2);

	return result;
}

double Simp(double a, double b, int n) 
{
	double h = (b - a) / n;
	double sum1 = 0;
	double sum2 = 0;

	for (int i = 1; i < n; i += 2) 
	{
		sum1 += Func(a + i * h);
		sum2 += Func(a + (i + 1) * h);
	}

	return h / 3 * (Func(a) + 4 * sum1 + 2 * sum2);
}

int main(void)
{
	double e = 0.1;

	while (e > EPS)
	{
		int n = 1;
		int i = 1;
		double S2N = 0, SN = 0;
		int func_count = 0;

		S2N = Simpson(A, B, n);
		func_count += 2 * n + 1;

		do
		{
			SN = S2N;
			i++;
			n = n * 2;
			S2N = Simpson(A, B, n);
			func_count += 2 * n + 1;
		} while (fabs(S2N - SN) / 15 > e);

		//printf("%d; ", i);
		//printf("%d; ", func_count);
		//double result = IntFunc(B) - IntFunc(A);
		//printf("%.16Lf; ", fabs(S2N - FVALUE));
		printf("%.16Lf; ", e);
		e = e * 0.1;
	}

	//int i = 1;
	//while (i <= 15)
	//{
	//	int n = (int)pow(2, i);
	//	double S2N = 0, SN = 0;

	//	S2N = Simpson(A, B, n);

	//	//printf("%d; ", i);
	//	double result = IntFunc(B) - IntFunc(A);
	//	//printf("%.16Lf; ", fabs(S2N - result));
	//	printf("%.16Lf; ", (B - A) / (2 * n));

	//	i++;
	//}

	return 0;
}