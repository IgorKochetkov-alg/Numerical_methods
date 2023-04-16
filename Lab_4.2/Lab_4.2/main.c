#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>

#define N 128
#define M 8

#define PI 3.1415926535897932

#define LEFT -1
#define RIGHT 1

double Func(double x)
{
    return pow(x, 2) + 1 - tan(x);
}

double Cheb(double x, int n)
{
    double result = 0;

    if (n < 0)
    {
        return 0;
    }

    if (n == 0)
    {
        return 1;
    }

    if (n == 1)
    {
        return x;
    }
    
    result = 2 * x * Cheb(x, n - 1) - Cheb(x, n - 2);

    return result;
}

void FindLU(double A[M + 1][M + 1], double L[M + 1][M + 1], double U[M + 1][M + 1])
{
	double sumU = 0;
	double sumL = 0;

	for (int i = 0; i <= M; ++i)
	{
		L[i][i] = 1;

		for (int j = 0; j <= M; ++j)
		{
			if (i <= j)
			{
				for (int k = 0; k < i; ++k)
				{
					sumU += L[i][k] * U[k][j];
				}

				U[i][j] = A[i][j] - sumU;
				sumU = 0;
			}

			if (i > j)
			{
				for (int k = 0; k < j; ++k)
				{
					sumL += L[i][k] * U[k][j];
				}

				L[i][j] = (A[i][j] - sumL) / U[j][j];
				sumL = 0;
			}
		}
	}
}

void Lyb(double L[M + 1][M + 1], double b[M + 1], double y[M + 1])
{
	double rightPart = 0;

	for (int i = 0; i <= M; i++)
	{
		rightPart = b[i];

		for (int j = 0; j <= M; j++)
		{
if (j < i)
{
	rightPart -= L[i][j] * y[j];
}
		}

		y[i] = rightPart;

		rightPart = 0;
	}
}

void Uxy(double U[M + 1][M + 1], double y[M + 1], double x[M + 1])
{
	double rightPart = 0;
	double k = 0;

	for (int i = M; i >= 0; i--)
	{
		rightPart = y[i];

		for (int j = 0; j <= M; j++)
		{
			if (j == i)
			{
				k = U[i][j];
			}

			if (j > i)
			{
				rightPart -= U[i][j] * x[j];
			}
		}

		x[i] = rightPart / k;

		rightPart = 0;
		k = 0;
	}
}

double LSM(double a, double b, double x, int n, int m)
{
	double nodes[N + 1] = { 0 };
	double h = (b - a) / n;
	double result = 0;
	double alfa[M + 1] = { 0 };

	for (int i = 0; i <= N; i++)
	{
		nodes[i] = a + h * i;
	}

	double R[M + 1] = { 0 };

	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			R[i] += Func(nodes[j]) * Cheb(nodes[j], i);
		}
	}

	double G[M + 1][M + 1] = { 0 };

	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			for (int k = 0; k <= n; k++)
			{
				G[i][j] += Cheb(nodes[k], i) * Cheb(nodes[k], j);
			}
		}
	}

	double L[M + 1][M + 1] = { 0 };
	double U[M + 1][M + 1] = { 0 };
	double y[M + 1] = { 0 };

	//double T[M + 1][M + 1] = { 1, 1, 3, 1, 2, 3, 3, 3, 9 };

	FindLU(G, L, U);
	Lyb(L, R, y);
	Uxy(U, y, alfa);

	for (int i = 0; i <= M; i++)
	{
		result += alfa[i] * Cheb(x, i);
	}

	return result;
}

int main()
{
	double result = 0;
	double h = (4.0 - 2.0) / 256;
	double x = 0;
	double max = 0;

		for (int i = 0; i <= 256; i++)
		{
			x = 2.0 + h * i;
			result = LSM(LEFT, RIGHT, x, N, M);
			/*result = fabs(result - Func(x));
			if (result > max)
			{
				max = result;
			}*/
		}
		printf("%.10Lf; ", result);
		max = 0;

	

    return 0;
}