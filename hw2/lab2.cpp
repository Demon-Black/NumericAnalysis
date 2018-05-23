#include <cstring>
#include <cstdio>
#include <cmath>
using namespace std;

#define N 7

void gauss(double* a, double* b, int n) {
    int i, j, k;
    double c[n], x[n];
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++)
            c[i] = (*(a + i * n + k)) / (*(a + k * n + k));
        for (i = k + 1; i < n; i++) {
            for (j = 0; j < n; j++)
                (*(a + i * n + j)) = (*(a + i * n + j)) - c[i] * (*(a + k * n + j));
            b[i] = b[i] - c[i] * b[k];
        }
    }

    x[n - 1] = b[n - 1] / *(a + n * n - 1);
    for (i = n - 2; i >= 0; i--) {
        double sum = 0;
        for (j = i + 1; j < n; j++)
            sum += (*(a + i * n + j)) * x[j];
        x[i] = (b[i] - sum) / (*(a + i * n + i));
    }

    for (int i = 0; i < n; i++)
        printf("%lf\n", x[i]);
}

int main(int argc, char const *argv[]) {
    double x0[N] = {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
    double y0[N] = {-4.467, -0.452, 0.551, 0.048, -0.447, 0.549, 4.552};
    int n;
    scanf("%d", &n);
    double a[n + 1][n + 1];
    double p[2 * n + 1];
    memset(p, 0.0, sizeof(double) * (2 * n + 1));
    double b[n + 1];
    memset(b, 0.0, sizeof(double) * (n + 1));
    for (int i = 1; i < 2 * n + 1; i++)
        for (int j = 0; j < N; j++)
            p[i] += pow(x0[j], i);
    p[0] = N;

    for (int i = 0; i < n + 1; i++)
        for (int j = 0;j < n + 1; j++)
            a[i][j] = p[i + j];
    for (int i = 0; i < n + 1; i++)
        for (int j = 0; j < N; j++)
            b[i] += pow(x0[j], i) * y0[j];

    gauss((double*)a, b, n + 1);

    return 0;
}