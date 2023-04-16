#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932
#define NUM 8
#define A 2.0
#define B 4.0


double Func(double x) 
{
    return pow(x, 2) + 1 - tan(x);
}
double DifFunc(double x) 
{
    return 2 * x - (1 / pow(cos(x), 2));
}

double SqrPol(double x)
{
    double result = fabs((x - 2.0) * (x - 3.0) * (x - 4.0) / 3);
    return result;
}

double Hermit(double a, double b, double x, double N) 
{
    double result = 0;
    double multiplication = 1, sum = 0;
    double h = (b - a) / N;

    for (int j = 0; j <= N; j++) 
    {
        multiplication = 1;
        sum = 0;

        for (int k = 0; k <= N; k++)
        {
            if (k != j)
            {
                sum += (x - (a + h * j)) / ((a + h * j) - (a + h * k));
            }
        }
            
        for (int i = 0; i <= N; i++)
        {
            if (i != j)
            {
                multiplication *= (x - (a + h * i)) / ((a + h * j) - (a + h * i)) * (x - (a + h * i)) / ((a + h * j) - (a + h * i));
            }
        }
            
        result += (((x - (a + h * j)) * DifFunc(a + h * j) + (1 - 2 * sum) * Func(a + h * j)) * multiplication);
    }

    return result;
}

int main(void) 
{
    double result = 0;
    double h = (B - A) / 32;
    double x = 0;
    double trueValue = 0;
    double max = 0;
    double maxX = 0;

    //for (int j = 1; j < 50; j++)
    //{
    //    max = -1000;
    //    maxX = 0;

    //    for (int i = 0; i < 33; i++)
    //    {
    //        x = A + i * h;
    //        result = Hermit(A, B, x, j);
    //        trueValue = Func(x);
    //        //fprintf(stdout, "%.10Lf; ", x);
    //        //fprintf(stdout, "%.10Lf; ", result);
    //        if (fabs(result - trueValue) > max)
    //        {
    //            max = fabs(result - trueValue);
    //            maxX = x;
    //        }
    //        //fprintf(stdout, "%.10Lf; ", fabs(result - trueValue));
    //    }

    //    fprintf(stdout, "%d; ", j);
    //}

    for (int i = 0; i < 33; i++)
    {
        x = A + i * h;
        result = Hermit(A, B, x, 2);
        trueValue = Func(x);
        result = SqrPol(x);
        //fprintf(stdout, "%.10Lf; ", x);
        fprintf(stdout, "%.10Lf; ", result);
        if (fabs(result - trueValue) > max)
        {
            max = fabs(result - trueValue);
            maxX = x;
        }
        //fprintf(stdout, "%.10Lf; ", fabs(result - trueValue));
    }
    
    return 0;
}