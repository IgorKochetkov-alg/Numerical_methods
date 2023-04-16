//#define _CRT_SECURE_NO_WARNINGS
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//
//#define LEFT 2 //границы
//#define RIGHT 4
//#define PI 3.1415926535897932
//
//#define N 128 //любые
//#define M 3

//double Cheb(double x, int n)
//{
//	double result = 0;
//
//	if (n < 0)
//	{
//		return 0;
//	}
//
//	if (n == 0)
//	{
//		return 1;
//	}
//
//	if (n == 1)
//	{
//		return x;
//	}
//
//	result = 2 * x * Cheb(x, n - 1) - Cheb(x, n - 2);
//
//	return result;
//}
//
//void FindLU(double** A, double** L, double** U, int m)
//{
//	double sumU = 0;
//	double sumL = 0;
//
//	for (int i = 0; i <= m; ++i)
//	{
//		L[i][i] = 1;
//
//		for (int j = 0; j <= m; ++j)
//		{
//			if (i <= j)
//			{
//				for (int k = 0; k < i; ++k)
//				{
//					sumU += L[i][k] * U[k][j];
//				}
//
//				U[i][j] = A[i][j] - sumU;
//				sumU = 0;
//			}
//
//			if (i > j)
//			{
//				for (int k = 0; k < j; ++k)
//				{
//					sumL += L[i][k] * U[k][j];
//				}
//
//				L[i][j] = (A[i][j] - sumL) / U[j][j];
//				sumL = 0;
//			}
//		}
//	}
//}
//
//void Lyb(double** L, double* b, double* y, int m)
//{
//	double rightPart = 0;
//
//	for (int i = 0; i <= m; i++)
//	{
//		rightPart = b[i];
//
//		for (int j = 0; j <= m; j++)
//		{
//			if (j < i)
//			{
//				rightPart -= L[i][j] * y[j];
//			}
//		}
//
//		y[i] = rightPart;
//
//		rightPart = 0;
//	}
//}
//
//void Uxy(double** U, double* y, double* x, int m)
//{
//	double rightPart = 0;
//	double k = 0;
//
//	for (int i = m; i >= 0; i--)
//	{
//		rightPart = y[i];
//
//		for (int j = 0; j <= m; j++)
//		{
//			if (j == i)
//			{
//				k = U[i][j];
//			}
//
//			if (j > i)
//			{
//				rightPart -= U[i][j] * x[j];
//			}
//		}
//
//		x[i] = rightPart / k;
//
//		rightPart = 0;
//		k = 0;
//	}
//}

//double LSM(double* nodes, double x, int n, int m)
//{
//	double result = 0;
//	double* alfa = (double*)calloc((m + 1), sizeof(double));
//	
//	//столбец свободных членов
//	double* R = (double*)calloc((m + 1), sizeof(double));
//
//	for (int i = 0; i <= m; i++)
//	{
//		for (int j = 0; j <= n; j++)
//		{
//			R[i] += Function(nodes[j]) * Cheb(nodes[j], i);
//		}
//	}
//	//матрица Грамма
//	double** G = (double**)calloc((m + 1), sizeof(double*));
//
//	for (int i = 0; i <= m; i++)
//	{
//		G[i] = (double*)calloc((m + 1), sizeof(double));
//	}
//
//	for (int i = 0; i <= m; i++)
//	{
//		for (int j = 0; j <= m; j++)
//		{
//			for (int k = 0; k <= n; k++)
//			{
//				G[i][j] += Cheb(nodes[k], i) * Cheb(nodes[k], j);
//			}
//
//		//printf("%.4Lf ", G[i][j]);
//		}
//		//printf("\n");
//	}
//
//	double** L = (double**)calloc((m + 1), sizeof(double*));
//	double** U = (double**)calloc((m + 1), sizeof(double*));
//
//	for (int i = 0; i <= m; i++)
//	{
//		L[i] = (double*)calloc((m + 1), sizeof(double));
//		U[i] = (double*)calloc((m + 1), sizeof(double));
//	}
//
//	double* y = (double*)calloc((m + 1), sizeof(double));
//
//	FindLU(G, L, U, m);
//	Lyb(L, R, y, m);
//	Uxy(U, y, alfa, m);
//
//	for (int i = 0; i <= m; i++)
//	{
//		//printf("%Lf\t", alfa[i]);
//		result += alfa[i] * Cheb(x, i);
//	}
//
//	return result;
//}

//int main(void)
//{
//	double* nodes = (double*)calloc(N + 1, sizeof(double));
//	double* allNodes = (double*)calloc(2 * N + 1, sizeof(double));
//
//	double result = 0;
//	double x = 0;
//	double max = 0;
//
//	for (int i = 0; i <= N; i++)
//	{
//		nodes[i] = (LEFT + RIGHT) / 2 + (RIGHT - LEFT) / 2 * cos((2 * i + 1) * PI / (2 * 256 + 2));
//	}
//
//	for (int i = 0; i <= 2 * N; i++)
//	{
//		nodes[i] = (LEFT + RIGHT) / 2 + (RIGHT - LEFT) / 2 * cos((2 * i + 1) * PI / (2 * 256 + 2));
//	}
//
//	for (int i = 0; i <= 256; i++)
//	{
//		result = LSM(nodes[i], N, M);
//		result = fabs(result - Function(x));
//
//		if (result > max)
//		{
//			max = result;
//		}
//	}
//	printf("%.10Lf; ", max);
//	max = 0;
//	return 0;
//}