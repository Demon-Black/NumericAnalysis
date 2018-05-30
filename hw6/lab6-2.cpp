#include <cmath>
#include <iostream>

using namespace std;

inline double phi1(double x) {
	return pow((x + 1) / 2, 1 / 3);
}

inline double phi2(double x) {
	return 2 * pow(x, 3) - 1;
}

double f_newton(double x){
	return x - (pow(x, 3) - x - 1) / (3 * pow(x, 2) - 1);
}

void newton(double x, double (*phi)(double)) {
	double x_k_1 = x;
	double x_k = double(x);
	for (int i = 0; i < 9; i++) {
		x_k_1 = x_k;
		x_k = phi(x_k_1);
	}
	cout << "after 10 iterations, the approximate solution of the equation is "
		 << x_k << endl;
}

void newton(double x, double epsilon) {
	double x_k_1 = x;
	double x_k = f_newton(x);
	int n = 1;
	while (abs(x_k - x_k_1) > epsilon) {
		x_k_1 = x_k;
		x_k = f_newton(x_k_1);
		n++;
	}
	cout << "after " << n << " iterations, the approximate solution of the equation is "
		 << x_k << endl;
}

int main() {
	newton(0, phi1);
	newton(0, phi2);
	newton(0, 1e-5);
	newton(1.5, 1e-5);
}