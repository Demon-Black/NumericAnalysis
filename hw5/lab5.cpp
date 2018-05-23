#include <iostream>
#include <cmath>
using namespace std;

// double* jacobi(double *a, double *b, int n, double* result, int t) {
// 	double *r = new double[n];
// 	for (int i = 0; i < n; i++) {
// 		r[i] = b[i];
// 		for (int j = 0; j < n; j++) {
// 			if (i == j)
// 				continue;
// 			r[i] = r[i] - a[i * n + j] * result[j];
// 		}
// 		r[i] = r[i] / a[i * n + i];
// 	}
// 	delete result;
// 	if (++t == 10000)
// 		return r;
// 	else
// 		return jacobi(a, b, n, r, t);
// }

double* jacobi(double *a, double *b, int n, int t) {
	double *r = new double[n];
	for (int i = 0; i < n; i++)
		r[i] = 0.0;
	for (int k = 0; k < t; k++) {
		for (int i = 0; i < n; i++) {
			double sum = 0.0;
			for (int j = 0; j < n; j++) 
				if (i != j)
					sum += a[i * n + j] * r[j];
			r[i] = (b[i] - sum) / a[i * n + i];
		}
	}
	return r;
}

// double* gauss_seidel(double *a, double *b, int n, double* result, int t) {
// 	double *r = new double[n];
// 	for (int i = 0; i < n; i++) {
// 		r[i] = b[i];
// 		for (int j = 0; j < i; j++)
// 			r[i] -= a[i * n + j] * r[j];
// 		for (int j = i + 1; j < n; j++) 
// 			r[i] -= a[i * n + j] * result[j];
// 		r[i] /= a[i * n + i];
// 	}
// 	delete result;
// 	if (++t == 10000)
// 		return r;
// 	else
// 		return gauss_seidel(a, b, n, r, t);
// }

double* gauss_seidel(double *a, double *b, int n, int t) {
	double *r = new double[n];
	for (int i = 0; i < n; i++)
		r[i] = 0.0;
	for (int k = 0; k < t; k++) {
		for (int i = 0; i < n; i++) {
			double s1 = 0.0, s2 =0.0;
			for (int j = 0; j < i; j++)
				s1 += a[i * n + j] * r[j];
			for (int j = i + 1; j < n; j++)
				s2 += a[i * n + j] * r[j];
			r[i] = (b[i] - s1 - s2) / a[i * n + i];
		}
	}
	return r;
}

double* SOR(double *a, double *b, int n, int t) {
	double w = 0.5;
	double *r = new double[n];
	for (int i = 0; i < n; i++)
		r[i] = 0.0;
	for (int k = 0; k < t; k++) {
		double r_[n];
		for (int i = 0; i < n; i++)
			r_[i] = r[i];
		for (int i = 0; i < n; i++) {
			double s1 = 0.0, s2 =0.0;
			for (int j = 0; j < i - 1; j++)
				s1 += a[i * n + j] * r[j];
			for (int j = i - 1; j < n; j++)
				s2 += a[i * n + j] * r_[j];
			r[i] = r_[i] + w * (b[i] - s1 - s2) / a[i * n + i];
		}
	}
	return r;
}

double f(double x) {
	return 0.5 * (1 - exp(-x)) / (1 - exp(-1)) + 0.5 * x;
}

int main(int argc, char const *argv[]) {
	const int n = 100;
	double h = 1.0 / n;
	const double epsilon = 1.0;
	double A[n][n] = {0.0};
	double b[n] = {0.0};
	
	for (int i = 0; i < n - 1; i++) {
		A[i][i] = - (2 * epsilon + h);
		A[i][i + 1] = epsilon + h;
		A[i + 1][i] = epsilon; 
	}
	A[n - 1][n - 1] = - (2 * epsilon + h);

	for (int i = 0; i < n; i++)
		b[i] = 0.5 * h * h;

	double *result = new double[n];
	int t = 10;

	for (int i = 0; i < n; i++)
		result[i] = 0.0;
	cout.precision(4);
	double *r =  jacobi((double *)A, b, n, t);
	for (int i = 0; i < n; i++)
		cout << r[i] << " ";
	cout << endl << endl;

	for (int i = 0; i < n; i++)
		result[i] = 0.0;
	r =  gauss_seidel((double *)A, b, n, t);
	for (int i = 0; i < n; i++)
		cout << r[i] << " ";
	cout << endl << endl;

	for (int i = 0; i < n; i++)
		result[i] = 0.0;
	r =  SOR((double *)A, b, n, t);
	for (int i = 0; i < n; i++)
		cout << r[i] << " ";
	cout << endl;


	for (int i = 0; i < n; i++)
		cout << f(i * 0.01) << " ";
	cout << endl;
	return 0;
}