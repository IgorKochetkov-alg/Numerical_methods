//#define _CRT_SECURE_NO_WARNINGS
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//
//#define LEFT 2
//#define RIGHT 4
//
//#define N 150
//#define M 30
//
//double Function(double x)
//{
//	return pow(x, 2) + 1 - tan(x);
//}
////+
//double Cheb(double x, int n)
//{
//	double prevRes = 1;
//	double result = x;
//	int count = 1;
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
//	while (count < n)
//	{
//		double save = result;
//		result = 2 * x * result - prevRes;
//		prevRes = save;
//		count++;
//	}
//
//	return result;
//}
////+
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
////+
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
////+
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
////+
//void FillNet(double* net, int n, double a, double b)
//{
//	double h = (b - a) / n;
//
//	for (int i = 0; i <= n; i++)
//	{
//		net[i] = a + h * i;
//		//printf("%.10Lf; ", net[i]);
//	}
//
//	return;
//}
////+
//void FillX(double* x, int n, double a, double b)
//{
//	double h = (b - a) / (2 * n);
//
//	for (int i = 0; i <= 2 * n; i++)
//	{
//		x[i] = a + i * h;
//		//printf("%.10Lf; ", x[i]);
//	}
//}
////+
//void LSM(double* net, double* x, int n, int m)
//{
//	double result = 0;
//	double max = 0;
//	double* alfa = (double*)calloc((m + 1), sizeof(double));
//
//	//столбец свободных членов
//	double* R = (double*)calloc((m + 1), sizeof(double));
//
//	for (int i = 0; i <= m; i++)
//	{
//		for (int j = 0; j <= n; j++)
//		{
//			R[i] += Function(net[j]) * Cheb(net[j], i);
//		}
//
//		printf("%.5Lf; ", R[i]);
//	}
//	printf("\n\n");
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
//				G[i][j] += Cheb(net[k], i) * Cheb(net[k], j);
//			}
//
//			printf("%.10Lf ", G[i][j]);
//		}
//		printf("\n");
//	}
//	printf("\n\n");
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
//	//for (int i = 0; i < (2 * n + 1); i++)
//	//{
//	//	result = 0;
//
//	//	for (int j = 0; j <= m; j++)
//	//	{
//	//		//printf("%Lf\t", alfa[i]);
//	//		result += alfa[j] * Cheb(x[i], j);
//	//	}
//
//	//	printf("%.16Lf; ", result);
//	//}
//
//	//	if (fabs(result - Function(x[i])) > max)
//	//	{
//	//		max = fabs(result - Function(x[i]));
//	//	}
//	//}
//
//	for (int j = 0; j <= m; j++)
//	{
//		printf("%.15Lf; ", alfa[j]);
//	}
//	printf("\n\n");
//
//	double sum = 0;
//
//	for (int i = 0; i <= n; i++)
//	{
//		result = 0;
//
//		for (int j = 0; j <= m; j++)
//		{
//			result += alfa[j] * Cheb(net[i], j);
//		}
//
//		sum += pow((result - Function(net[i])), 2);
//	}
//
//	printf("%.15Lf; ", (sum / (n + 1)));
//
//	free(alfa);
//	free(R);
//	free(y);
//
//	for (int i = 0; i <= m; i++)
//	{
//		free(G[i]);
//		free(L[i]);
//		free(U[i]);
//	}
//
//	free(G);
//	free(L);
//	free(U);
//
//	return;
//}
//
//int main()
//{
//	double n = 128;
//	double m = 4;
//
//	double* net = (double*)calloc(n + 1, sizeof(double));
//
//	FillNet(net, 128, -1, 1);
//
//	double* x = (double*)calloc(2 * n + 1, sizeof(double));
//
//	FillX(x, 128, 2.0, 4.0);
//
//	for (int i = 3; i < 5; i++)
//	{
//		LSM(net, x, 128, i);
//	}
//
//	/*double* chebNet = (double*)calloc(n + 1, sizeof(double));
//	FillNet(chebNet, n, -1, 1);
//
//	double* x = (double*)calloc(2 * n + 1, sizeof(double));
//	FillX(2, 4, x, n);
//
//	LSM(chebNet, x, n, m);
//	printf("\n\n");
//	free(chebNet);
//	free(x);*/
//	
//
//	/*while (m <= 4)
//	{
//		double* chebNet = (double*)calloc(n + 1, sizeof(double));
//		FillNet(chebNet, n, -1, 1);
//
//		double* x = (double*)calloc(2 * n + 1, sizeof(double));
//		FillX(chebNet, x, n);
//
//		LSM(chebNet, x, n, m);
//		printf("\n\n");
//		m = m + 1;
//		free(chebNet);
//		free(x);
//	}*/
//	//LSM(chebNet, x, n, m);
//}