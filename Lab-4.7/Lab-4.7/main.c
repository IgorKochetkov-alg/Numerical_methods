#define _CRT_SECURE_NOWARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LEFT 1.0
#define RIGHT 2.0 //y'' + y = 0 sinx
#define MAX 4096
#define TEST_N 8
#define MAX_ERROR 12

double yFunc(double x)
{
	return log(x);
	//return sin(x);
}

double dyFunc(double x, double y, double z)
{
	return z;
}

double U_d2yFunc(double x, double y, double z)
{
	return ((2 * log(x) - 2 * y - z) / x);
	//return -y;
}

double V_d2yFunc(double x, double y, double z)
{
	return (( - 2 * y - z) / x);
	//return -y;
}

void U_RungeKutta(double a, double b, double n, double y0, double z0, double* u, double* du) // n - number of split segments
{
	double h = (b - a) / n;
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0; // for y
	double q1 = 0, q2 = 0, q3 = 0, q4 = 0; // for z
	double x = a, y = y0, z = z0;
	u[0] = y;
	du[0] = z;
	//printf("%.15Lf; ", x);
	//printf("%.15Lf; ", y);
	//printf("%.15Lf; ", z);

	for (int i = 1; i <= n; i++)
	{
		k1 = dyFunc(x, y, z);
		q1 = U_d2yFunc(x, y, z);

		k2 = dyFunc(x + h / 3, y + k1 * h / 3, z + q1 * h / 3);
		q2 = U_d2yFunc(x + h / 3, y + k1 * h / 3, z + q1 * h / 3);

		k3 = dyFunc(x + 2 * h / 3, y - k1 * h / 3 + k2 * h, z - q1 * h / 3 + q2 * h);
		q3 = U_d2yFunc(x + 2 * h / 3, y - k1 * h / 3 + k2 * h, z - q1 * h / 3 + q2 * h);

		k4 = dyFunc(x + h, y + k1 * h - k2 * h + k3 * h, z + q1 * h - q2 * h + q3 * h);
		q4 = U_d2yFunc(x + h, y + k1 * h - k2 * h + k3 * h, z + q1 * h - q2 * h + q3 * h);

		y = y + h * (k1 + 3 * k2 + 3 * k3 + k4) / 8.0;
		z = z + h * (q1 + 3 * q2 + 3 * q3 + q4) / 8.0;

		x = x + h;
		u[i] = y;
		du[i] = z;
		//printf("%.15Lf; ", x);
		//printf("%.15Lf; ", u[i]);
		//printf("%.15Lf; ", z);
	}
	return;
}

void V_RungeKutta(double a, double b, double n, double y0, double z0, double* u, double* du) // n - number of split segments
{
	double h = (b - a) / n;
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0; // for y
	double q1 = 0, q2 = 0, q3 = 0, q4 = 0; // for z
	double x = a, y = y0, z = z0;
	u[0] = y;
	du[0] = z;
	//printf("%.15Lf; ", x);
	//printf("%.15Lf; ", y);
	//printf("%.15Lf; ", z);

	for (int i = 1; i <= n; i++)
	{
		k1 = dyFunc(x, y, z);
		q1 = V_d2yFunc(x, y, z);

		k2 = dyFunc(x + h / 3, y + k1 * h / 3, z + q1 * h / 3);
		q2 = V_d2yFunc(x + h / 3, y + k1 * h / 3, z + q1 * h / 3);

		k3 = dyFunc(x + 2 * h / 3, y - k1 * h / 3 + k2 * h, z - q1 * h / 3 + q2 * h);
		q3 = V_d2yFunc(x + 2 * h / 3, y - k1 * h / 3 + k2 * h, z - q1 * h / 3 + q2 * h);

		k4 = dyFunc(x + h, y + k1 * h - k2 * h + k3 * h, z + q1 * h - q2 * h + q3 * h);
		q4 = V_d2yFunc(x + h, y + k1 * h - k2 * h + k3 * h, z + q1 * h - q2 * h + q3 * h);

		y = y + h * (k1 + 3 * k2 + 3 * k3 + k4) / 8.0;
		z = z + h * (q1 + 3 * q2 + 3 * q3 + q4) / 8.0;

		x = x + h;
		u[i] = y;
		du[i] = z;
		//printf("%.15Lf; ", x);
		//printf("%.15Lf; ", u[i]);
		//printf("%.15Lf; ", z);
	}
	return;
}


int main(void)
{
	int n = 2;
	double result = 0;
	double h = 0;
	double c = 0;
	double maxError = 0;

	//while (n <= MAX)
	//{
	//	double* u = (double*)calloc(n + 1, sizeof(double));
	//	double* v = (double*)calloc(n + 1, sizeof(double));
	//	double* du = (double*)calloc(n + 1, sizeof(double));
	//	double* dv = (double*)calloc(n + 1, sizeof(double));

	//	maxError = 0;
	//	h = (RIGHT - LEFT) / n;

	//	U_RungeKutta(LEFT, RIGHT, n, 0, 1, u, du);
	//	V_RungeKutta(LEFT, RIGHT, n, 0, -1, v, dv);

	//	c = ((log(RIGHT) + (1 / RIGHT)) - (du[n] + u[n])) / (dv[n] + v[n]);

	//	for (int i = 0; i <= n; i++)
	//	{
	//		//printf("%.15Lf; ", LEFT + i * h);
	//		result = c * v[i] + u[i];

	//		if (fabs(result - yFunc(LEFT + i * h)) > maxError)
	//		{
	//			maxError = fabs(result - yFunc(LEFT + i * h));
	//		}

	//		
	//	}
	//	//printf("%.15Lf; ", maxError);
	//	printf("%.15Lf; ", h);
	//	n = n * 2;
	//	free(u);
	//	free(du);
	//	free(v);
	//	free(dv);
	//}

	//возмущение
	//n = TEST_N;
	double eps = 0;

	for (int i = 1; i <= 12; i++)
	{
		eps = pow(10, -i);
		double* u = (double*)calloc(n + 1, sizeof(double));
		double* v = (double*)calloc(n + 1, sizeof(double));
		double* du = (double*)calloc(n + 1, sizeof(double));
		double* dv = (double*)calloc(n + 1, sizeof(double));

		maxError = 0;
		h = (RIGHT - LEFT) / n;

		U_RungeKutta(LEFT, RIGHT, 2, 0, 1, u, du);
		V_RungeKutta(LEFT, RIGHT, 2, 0, -1, v, dv);

		c = ((log(RIGHT) + (1 / RIGHT)) - (du[n] + u[n])) / (dv[n] + v[n]);

		for (int i = 0; i <= n; i++)
		{
			//printf("%.15Lf; ", LEFT + i * h);
			result = c * v[i] + u[i];

			if (fabs(result - yFunc(LEFT + i * h)) > maxError)
			{
				maxError = fabs(result - yFunc(LEFT + i * h));
			}


		}
		printf("%.15Lf; ", maxError);
		//printf("%.15Lf; ", eps);
		free(u);
		free(du);
		free(v);
		free(dv);
	}

	return 0;
}