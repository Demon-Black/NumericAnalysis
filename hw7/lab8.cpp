#include <iostream>
#include <cmath>
using namespace std;

const double A1[3][3] = {
	{ 5, -4,  1},
	{-4,  6, -4},
	{ 1, -4,  7}
};

const double A2[4][4] = {
	{ 25, -41,  10, -6},
	{-41,  68, -17, 10},
	{ 10, -17,   5, -3},
	{ -6,  10,  -3,  2}
};

void matrixx(double *A, double *x, double *v, int n) {
    for(int i = 0; i < n; i++) {
        v[i] = 0;
        for(int j = 0; j < n; j++)
            v[i] += A[i * n + j] * x[j];
    }
}

double slove(double *v, int n) {
    double result;
    for(int i = 0; i < n - 1; i++) 
    	result = v[i] > v[i+1] ? v[i] : v[i+1];
    return result;
}

void get_result(int n, double *A) {
	double x[n], v[n], u[n], p[n];
	for (int i = 0; i < n; i++) {
		x[i] = 1.0;
		v[i] = u[i] = p[i] = 0.0;
	}

	double epsilon = 1e-5, delta = 1;
	int k = 0;
	while (delta >= epsilon) {
		for (int q = 0;q < n; q++) 
			p[q] = v[q];
        matrixx(A, x, v, n);
        for(int i = 0; i < n; i++)
        	u[i] = v[i] / slove(v, n);
        delta = fabs(slove(v, n)-slove(p, n));
        k++;
        for(int l = 0; l < n; l++) 
        	x[l] = u[l];
	}

	cout << "iteration times:" << k << endl;
    cout << "result:" << slove(v, n) << endl;
    for(int i = 0; i < n;i++) 
    	cout << u[i] << " " ;
    cout << endl << endl;
    return;
}

int main() {
    get_result(3, (double*) A1);
    get_result(4, (double*) A2);
    return 0;
} 
