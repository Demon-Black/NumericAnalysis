#include <cstdio>
#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

void LU(int n, double* a, double* b, double* x) {
	double l[n][n] = {0};
	double u[n][n] = {0};
	for (int i = 0; i < n; i++)
		u[0][i] = a[i];
	for (int i = 0; i < n; i++)
		l[i][0] = a[i * n] / u[0][0];
	for (int i = 1; i < n; i++) {
		for (int j = i; j < n; j++) {
			double sum = 0.0;
			for (int k = 0; k < i; k++)
				sum += l[i][k] * u[k][j];
			u[i][j] = a[i * n + j] - sum;
		}

		for (int j = i + 1; j < n; j++) {
			double sum = 0.0;
			for (int k = 0; k < i; k++)
				sum += l[j][k] * u[k][i];
			l[j][i] = (a[j * n + i] - sum) / u[i][i];
		}
	}

	double y[n] = {0};
	y[0] = b[0];
	for (int i = 1; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += l[i][j] * y[j];
		y[i] = b[i] - sum;
	}

	x[n - 1] = y[n - 1] / u[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < n; j++)
			sum += u[i][j] * x[j];
		x[i] = (y[i] - sum) / u[i][i];
	}
}

void inverse(int n, double* a, double* a_i) {
	double l[n][n] = {0};
	double u[n][n] = {0};
	for (int i = 0; i < n; i++)
		u[0][i] = a[i];
	for (int i = 0; i < n; i++)
		l[i][0] = a[i * n] / u[0][0];
	for (int i = 1; i < n; i++) {
		for (int j = i; j < n; j++) {
			double sum = 0.0;
			for (int k = 0; k < i; k++)
				sum += l[i][k] * u[k][j];
			u[i][j] = a[i * n + j] - sum;
		}

		for (int j = i + 1; j < n; j++) {
			double sum = 0.0;
			for (int k = 0; k < i; k++)
				sum += l[j][k] * u[k][i];
			l[j][i] = (a[j * n + i] - sum) / u[i][i];
		}
	}

	for (int s = 0;s < n; s++) {
		double b[n] = {0};
		b[s] = 1.0;

		double y[n] = {0};
		y[0] = b[0];
		for (int i = 1; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < i; j++)
				sum += l[i][j] * y[j];
			y[i] = b[i] - sum;
		}

		a_i[s + n * (n - 1)] = y[n - 1] / u[n - 1][n - 1];
		for (int i = n - 2; i >= 0; i--) {
			double sum = 0;
			for (int j = i + 1; j < n; j++)
				sum += u[i][j] * a_i[s + n * j];
			a_i[i * n + s] = (y[i] - sum) / u[i][i];
		}
	}
} 

void hilbert(int n, double* h) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			h[i * n + j] = 1.0 / (i + j + 1);
} 

double condition(int n, double* h) {
	double h_i[n][n];
	inverse(n, h, (double*)h_i);
	double c = 0.0;
	for (int i = 0; i < n; i++) {
		double line_sum = 0.0;
		for (int j = 0; j < n; j++)
			line_sum += abs(h[i * n + j]);
		if (line_sum > c)
			c = line_sum;
	}

	double c_i = 0.0;
	for (int i = 0; i < n; i++) {
		double line_sum = 0.0;
		for (int j = 0; j < n; j++)
			line_sum += abs(h_i[i][j]);
		if (line_sum > c_i)
			c_i = line_sum;
	}

	return c * c_i;
}

void cholesky(int n, double* a, double* b, double* x) {
    double t[n][n];
    double d[n];
    double l[n][n];
    for (int i = 0; i < n; i++) {  
        double sum = 0.0;
        for (int j = 0; j <= i - 1; j++) {
            for (int k = 0; k <= j - 1; k++)
                sum += t[i][k] * l[j][k];
            t[i][j] = a[i * n + j] - sum;
            l[i][j] = t[i][j] / d[j];
        }
        double sum_ = 0.0;
        for (int k = 0; k <= i - 1; k++)
            sum_ += t[i][k] * l[i][k];
        d[i] = a[i * n + i] - sum_;
    }

    double y[n];
    y[0] = b[0];
    for (int i = 1; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k <= i - 1; k++)
            sum += l[i][k] * y[k];
        y[i] = b[i] - sum;
    }
    
    x[n - 1] = y[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        double sum = 0.0;
        for (int k = i + 1; k < n; k++)
            sum += l[k][i] * x[k];
        x[i] = y[i] / d[i] - sum;
    }
}

void solve(int n, double r) {
	double hn[n][n];
	hilbert(n, (double*)hn);
	double b[n];
	double x[n];
	for (int i = 0; i < n; i++) {
		x[i] = 1.0;
		b[i] = r;
	}
	double x_[n] = {0.0};
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			b[i] += x[j] * hn[i][j];

	LU(n ,(double*)hn, b, x_);
	// for (int i = 0; i < n; i++)
	// 	cout << x_[i] << endl;
	double delta1[n], delta2[n];

	for (int i = 0; i < n; i++) {
		double t = 0.0;
		for (int j = 0; j < n; j++)
			t += x_[j] * hn[i][j];
		delta1[i] = b[i] - t;
		delta2[i] = x_[i] - x[i];
	}

	double d1 = 0.0;
	double d2 = 0.0;

	for (int i = 0; i < n; i++) {
		if (abs(delta1[i]) > d1)
			d1 = abs(delta1[i]);
		if (abs(delta2[i]) > d2)
			d2 = abs(delta2[i]);
	}

	printf("%.08lf %.08lf\n", d1, d2);
}

void test_time(int n) {
	double hn[n][n];
	hilbert(n, (double*)hn);
	double b[n];
	double x[n];
	for (int i = 0; i < n; i++) {
		x[i] = 1.0;
		b[i] = 0.0;
	}
	double x_[n] = {0.0};
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			b[i] += x[j] * hn[i][j];

	clock_t start = clock();
	for (int i = 0; i < 1000000; i++)
		LU(n ,(double*)hn, b, x_);
    clock_t ends = clock();
    cout <<"LU Running Time : "<< (double)(ends - start) / CLOCKS_PER_SEC << endl;
    start = clock();
	for (int i = 0; i < 1000000; i++)
		cholesky(n ,(double*)hn, b, x_);
    ends = clock();
    cout <<"cholesky Running Time : "<< (double)(ends - start) / CLOCKS_PER_SEC<< endl;
}

int main() {
	// double t[3][3] = {
	// 				{1.0, 1 / 2.0, 1 / 3.0}, 
	// 				{1 / 2.0, 1 / 3.0, 1 / 4.0}, 
	// 				{1 / 3.0, 1 / 4.0, 1 / 5.0}
	// 			};
	// double t_i[3][3];
	// inverse(3, (double*)t, (double*)t_i);
	// cout << t_i[0][0] << " " << t_i[0][1] << " " << t_i[0][2] <<endl;
	// cout << t_i[1][0] << " " << t_i[1][1] << " " << t_i[1][2] <<endl;
	// cout << t_i[2][0] << " " << t_i[2][1] << " " << t_i[2][2] <<endl;
	
	double h3[3][3];
	double h4[4][4];
	hilbert(3, (double*)h3);
	hilbert(4, (double*)h4);
	printf("%.08lf\n", condition(3, (double*)h3));
	printf("%.08lf\n", condition(4, (double*)h4));

	solve(10, 0.0);

	test_time(10);

	solve(10, 1e-7);
	for (int i = 11; i < 20; i++) {
		printf("%d ", i);
		solve(i, 0.0);
	}
	return 0;
}