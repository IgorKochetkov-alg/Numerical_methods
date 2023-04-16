#define _CRT_SECURE_NOWARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LEFT 1.1
#define RIGHT 4.0
#define EPS 1e-12
#define Y 0.4582575694
#define TEST_EPS 0.00001

double Func(double x, double y)
{
	return ((x * y) / (pow(x, 2) - 1) + (x / y)) / 2;
}
//+
double TrueFunc(double x)
{
	return pow(pow(x, 2) - 1, 0.5);
}
//+
double RungeKuttaStep(double prevX, double prevY, double h)
{
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	double curY = 0;

	k1 = Func(prevX, prevY);
	k2 = Func(prevX + h / 3, prevY + k1 * h / 3);
	k3 = Func(prevX + 2 * h / 3, prevY - k1 * h / 3 + k2 * h);
	k4 = Func(prevX + h, prevY + k1 * h - k2 * h + k3 * h);
	curY = prevY + h * (k1 + 3 * k2 + 3 * k3 + k4) / 8.0;

	return curY;
}
//+
double RK(double prevX, double prevY, double h, double eps)
{
	double tmpH = h;
	double curY = 0, lastY = 0;
	int iter = 1;
	curY = RungeKuttaStep(prevX, prevY, h);
	double x = 0;

	do
	{
		int count = pow(2, iter);
		iter++;
		tmpH = tmpH / 2.0;
		lastY = curY;
		curY = prevY;
		x = prevX;

		for (int i = 0; i < count; i++)
		{
			curY = RungeKuttaStep(x, curY, tmpH);
			x += tmpH;
		}

	} while ((fabs(curY - lastY) / 15) > eps);

	//printf("%d; ", iter);

	return curY;
}



void EPS_RungeKutta(double a, double b, int n, double eps, double y0)
{
	double curY = y0;
	double x = LEFT;
	double infNorm = 0;
	int itr = 0;
	double h = (b - a) / (double)n;

	for (int i = 1; i <= n; i++)
	{
		curY = RK(x, curY, h, eps);

		x += h;

		if (fabs(curY - TrueFunc(x)) > infNorm)
		{
			infNorm = fabs(curY - TrueFunc(x));
		}
	}

	printf("%.15Lf; ", infNorm * 10);

	//printf("%.10Lf; ", maxIter);
}

void H_RungeKutta(double a, double b, int n)
{
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	double h = (b - a) / (double)n;
	double curY = Y;
	double prevY = 0;
	double x = 0;

	//printf("%.10Lf; ", curY);
	printf("%.10Lf; ", curY - TrueFunc(LEFT));


	for (int i = 1; i <= n; i++)
	{
		x = a + (i - 1) * h;
		//printf("%.4Lf; ", x);
		k1 = Func(x, curY);
		k2 = Func(x + h / 3, curY + k1 * h / 3);
		k3 = Func(x + 2 * h / 3, curY - k1 * h / 3 + k2 * h);
		k4 = Func(x + h, curY + k1 * h - k2 * h + k3 * h);
		prevY = curY;
		curY = prevY + h * (k1 + 3 * k2 + 3 * k3 + k4) / 8;
		//printf("%.10Lf; ", curY);

		printf("%.10Lf; ", curY - TrueFunc(x + h));
	}

	return;
}

int main()
{
	//H_RungeKutta(LEFT, RIGHT, 58);

	double eps = 0;
	double y0 = Y;

	for (int i = 1; i < 13; i++)
	{
		eps = pow(0.1, i);

		y0 = Y * (1 + eps);
		EPS_RungeKutta(LEFT, RIGHT, 29, TEST_EPS, y0);
		//printf("\n\n");
		//printf("%.13Lf; ", y0 - Y);
	}

	return 0;
}