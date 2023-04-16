#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LEFT 2
#define RIGHT 4

#define N 150
#define M 30

double Function(double x)
{
	return pow(x, 2) + 1 - tan(x);
}
//+
double Cheb(double x, int n)
{
	double prevRes = 1;
	double result = x;
	int count = 1;

	if (n == 0)
	{
		return 1;
	}

	if (n == 1)
	{
		return x;
	}

	while (count < n)
	{
		double save = result;
		result = 2 * x * result - prevRes;
		prevRes = save;
		count++;
	}

	return result;
}
//+
void FillNet(double* net, int n, double a, double b)
{
	double h = (b - a) / n;

	for (int i = 0; i <= n; i++)
	{
		net[i] = a + h * i;
		//printf("%.10Lf; ", net[i]);
	}

	return;
}
//+
void FillX(double* x, int n, double a, double b)
{
	double h = (b - a) / (2 * n);

	for (int i = 0; i <= 2 * n; i++)
	{
		x[i] = a + i * h;
		//printf("%.10Lf; ", x[i]);
	}
}
//+
void doCholesky(double** A, double** L, int m) 
{
	for (int i = 0; i <= m; i++) 
	{
		double sum;
		for (int j = 0; j < i; j++) 
		{
			sum = 0;
			for (int k = 0; k < j; k++)
				sum += L[i][k] * L[j][k];
			L[i][j] = (A[i][j] - sum) / L[j][j];
		}
		sum = A[i][i];
		for (int k = 0; k < i; k++)
			sum -= L[i][k] * L[i][k];
		L[i][i] = sqrt(sum);
	}
}

void getTrans(double** L, double** Lt, int m) {
	for (int i = 0; i <= m; i++)
		for (int j = 0; j <= m; j++)
			Lt[i][j] = L[j][i];
}

void getSolution(double** L, double* B, double* X, int m) 
{
	double** Lt = (double**)calloc((m + 1), sizeof(double*));
	double* v = (double*)calloc((m + 1), sizeof(double));

	for (int i = 0; i <= m; i++)
	{
		Lt[i] = (double*)calloc((m + 1), sizeof(double));
	}

	getTrans(L, Lt, m);

	// Lv=B
	for (int i = 0; i <= m; i++) {
		double alpha = B[i];
		for (int j = 0; j < i; j++)
			alpha -= L[i][j] * v[j];
		v[i] = alpha / L[i][i];
	}
	// LtX=y
	for (int i = m; i >= 0; i--) {
		double alpha = v[i];
		for (int j = m; j > i; j--)
			alpha -= Lt[i][j] * X[j];
		X[i] = alpha / Lt[i][i];
	}
}

void LSM(double* net, double* x, int n, int m)
{
	double result = 0;
	double max = 0;
	double* alfa = (double*)calloc((m + 1), sizeof(double));

	//столбец свободных членов
	double* R = (double*)calloc((m + 1), sizeof(double));

	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			R[i] += Function(net[j]) * Cheb(net[j], i);
		}

		printf("%.5Lf; ", R[i]);
	}
	printf("\n\n");
	//матрица Грамма
	double** G = (double**)calloc((m + 1), sizeof(double*));

	for (int i = 0; i <= m; i++)
	{
		G[i] = (double*)calloc((m + 1), sizeof(double));
	}

	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			for (int k = 0; k <= n; k++)
			{
				G[i][j] += Cheb(net[k], i) * Cheb(net[k], j);
			}

			printf("%.10Lf ", G[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");

	double** L = (double**)calloc((m + 1), sizeof(double*));

	for (int i = 0; i <= m; i++)
	{
		L[i] = (double*)calloc((m + 1), sizeof(double));
	}

	doCholesky(G, L, m);
	getSolution(L, R, alfa, m);

	for (int j = 0; j <= m; j++)
	{
		printf("%.15Lf; ", alfa[j]);
	}
	printf("\n\n");

	double sum = 0;

	for (int i = 0; i <= n; i++)
	{
		result = 0;

		for (int j = 0; j <= m; j++)
		{
			result += alfa[j] * Cheb(x[i], j);
		}

		sum += pow((result - Function(net[i])), 2);

		/*if (fabs(result - Function(x[i])) > max)
		{
			max = fabs(result - Function(x[i]));
		}*/
	}

	//printf("%.15Lf; ", max);
	printf("%.15Lf; ", (sum / (n + 1)));

	free(alfa);
	free(R);

	for (int i = 0; i <= m; i++)
	{
		free(G[i]);
		free(L[i]);
	}

	free(G);
	free(L);

	return;
}

int main()
{
	double n = 128;
	double m = 4;

	double* net = (double*)calloc(n + 1, sizeof(double));

	FillNet(net, 128, -1, 1);

	double* x = (double*)calloc(2 * n + 1, sizeof(double));

	FillX(x, 128, -1.0, 1.0);

	for (int i = 2; i < 5; i++)
	{
		LSM(net, x, 128, i);
		printf("\n\n");
	}

	/*double* chebNet = (double*)calloc(n + 1, sizeof(double));
	FillNet(chebNet, n, -1, 1);

	double* x = (double*)calloc(2 * n + 1, sizeof(double));
	FillX(2, 4, x, n);

	LSM(chebNet, x, n, m);
	printf("\n\n");
	free(chebNet);
	free(x);*/
	

	/*while (m <= 4)
	{
		double* chebNet = (double*)calloc(n + 1, sizeof(double));
		FillNet(chebNet, n, -1, 1);

		double* x = (double*)calloc(2 * n + 1, sizeof(double));
		FillX(chebNet, x, n);

		LSM(chebNet, x, n, m);
		printf("\n\n");
		m = m + 1;
		free(chebNet);
		free(x);
	}*/
	//LSM(chebNet, x, n, m);
}