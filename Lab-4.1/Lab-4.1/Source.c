#include <stdio.h>
#include <math.h>

#define MAS_LENGTH 3
#define A 0
#define B 6
#define M_PI 3.14159265358979323846L
#define double long double

double x[MAS_LENGTH], y[MAS_LENGTH];

double func(double point) {
    return point * log(point + 1.0) - 0.5;
}

double polynom(double point, double* m, int size) {
    double p = 0, q = 1;

    for (int i = 0; i < size; i++) {
        q = 1;

        for (int j = 0; j < i; j++)
            q *= (point - x[size - 1 - j]);

        p += m[i] * q;
    }

    return p;
}

int main() {
    double prev[MAS_LENGTH] = { 0 }, cur[MAS_LENGTH] = { 0 }, ans[MAS_LENGTH] = { 0 };

    int n;
    scanf("%d", &n);
    if (n < 2)
        return 0;

    for (int i = 0; i < n; i++) {
        x[i] = ((double)(A + B)) / 2.0L + ((double)(B - A)) / 2.0L * cosl(M_PI * (2.0L * (n - 1 - i) + 1.0L) / (2.0L * n));
        //         x[i] = A + (double)(B - A) * (i) / (n-1);
        y[i] = func(x[i]);
        prev[i] = y[i];
    }

    //cur[MAS_LENGTH - 1] = prev[MAS_LENGTH - 1];
    cur[0] = y[n - 1];

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++) {
            cur[i] = (prev[j] - prev[j + 1]) / (x[j] - x[j + i]);
            prev[j] = cur[i];
        }
    }

    printf("X fun_val poly_val\n");
    for (int i = 0; i < n; i++) {
        printf("%0.25Lf  %0.25Lf  %0.25Lf\n", x[i], y[i], polynom(x[i], cur, n));
    }
    //     for (int i = 0; i < n; i++) {
    //         printf("%0.40Lf\n", cur[i]);
    //     }
    return 0;
}