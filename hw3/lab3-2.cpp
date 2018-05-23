#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

#define f(x) (4 / (1 + (x) * (x)))
#define pi 3.1415926535897932384626433

int main () {
	int n = 13;
	double h = 1.0 / n;
	double result = 0.0;
	for(int i = 0; i < n; i++) 
		result += f((i + 0.5) * h - h / (2 * sqrt(3))) + f((i + 0.5)* h + h / (2 * sqrt(3)));
	result = result * h / 2;
	printf("%.12lf\n", result);
	printf("%.12lf\n", result - pi);
	return 0;
}