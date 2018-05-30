#include <cmath>
#include <iostream>

using namespace std;

inline double f(double x) {
	return 2 * pow(x, 3) - pow(x, 2) + 3 * x - 1;
}

inline double f_1(double x) {
	return 6 * pow(x, 2) - 2 * x + 3;
}

inline double f_newton(double x) {
	return x - f(x) / f_1(x);
}

void dichotomy(double a, double b, double epsilon) {
	cout << "dichotomy initial value: a = " << a << " b = " << b <<endl;
	int n = 1;
	while (abs(a - b) > epsilon) {
		if (f((a + b) / 2) == 0)
			break;
		if (f((a + b) / 2) * f(b) < 0)
			a = (a + b) / 2;
		else 
			b = (a + b) / 2;
		cout << "iteration " << n << " : a = " << a << " b = " << b << endl;
		n++;
	}
	cout << "after " << n << " iterations, the approximate solution of the equation is "
		 << (a + b) / 2 << endl;
}

void newton(double x, double epsilon) {
	double x_k_1 = x;
	double x_k = f_newton(x);
	int n = 1;
	cout << "newton initial value: x0 = " << x << endl;
	while (abs(x_k - x_k_1) > epsilon) {
		cout << "iteration " << n << " : x = " << x_k << endl;
		x_k_1 = x_k;
		x_k = f_newton(x_k_1);
		n++;
	}
	cout << "after " << n << " iterations, the approximate solution of the equation is "
		 << x_k << endl;
}

int main() {
	dichotomy(-3, 3, 1e-5);
	newton(0, 1e-5);
}