#include <cstdio>
#include <cmath>
#include <iostream>
using namespace std;

void T() {
	int n = 496;
	double h = 1.0 / n;
	double result = 0.0;
	for (int i = 1; i < n; i++)
		result += 2 * exp(i * h);
	result = 0.5 * h * (result + exp(0) + exp(1));
	printf("%.12lf\n", result);
	printf("%.12lf\n", result - exp(1) + 1);
}

void S() {
	int n = 6;
	double h = 1.0 / n;
	double result = 0.0;

	for (int i = 0; i < n; i++)
		result += 4 * exp((i + 0.5) * h);
	for (int i = 1; i < n; i++)
		result += 2 * exp(i * h);

	result = h * (result + exp(0) + exp(1)) / 6;
	printf("%.12lf\n", result);
	printf("%.12lf\n", result - exp(1) + 1);
}

void R() {
	double h[100] = {0};
	double t[100][100] = {0}; 
	h[0] = 1.0;
	h[1] = 0.5;
	t[0][0] = h[0] * (exp(0) + exp(1)) / 2;
	t[1][0] = t[0][0] / 2 + h[1] * (exp(h[1]));
	t[0][1] = 4 * t[1][0] / 3 - 1 * t[0][0] / 3;

	int n = 1;

	while(abs(t[0][n] - t[0][n - 1]) > 1e-6) {
		n++;
		double t1 = n;
		double t2 = pow(2, t1 - 1);
		double t3 = 0;
		h[n] = h[n - 1] / 2;
		t[n][0] = t[n - 1][0] / 2;
		
		for(int i = 0; i < t2; i++)
			t[n][0] += h[n] * exp((i + 0.5) * h[n - 1]);
		for(int i = 1; i <= n; i++){
			t3 = i;
			t[n - i][i] = pow(4, t3) * t[n - i + 1][i - 1] / (pow(4, t3) - 1) - t[n - i][i - 1] / (pow(4, t3) - 1);
		}
	}

	printf("%.12lf\n", t[0][n]);
	printf("%.12lf\n", t[0][n] - exp(1) + 1);
	printf("%d\n", n);
}
 
int main() {
	T();
	cout << endl;
	S();
	cout << endl;
	R();
}