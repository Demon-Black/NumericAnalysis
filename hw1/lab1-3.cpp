#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

float matrix[100][100] = {};

void print_matrix(int n) {
	for(int i = 0;i <= n;i++) {
		for(int j = 0;j <= n + 1;j++)
			cout<<matrix[i][j] << ' ';
		cout << endl;
	}
	cout << endl;
}

int main() {
	int n;
	cin >> n;
	float h = 10.0 / n;
	float x_spot[n + 1];

	for(int i = 0;i <= n;i++)
		x_spot[i] = -5.0 + i*h;
	float y_spot[n + 1];
	for(int i = 0;i <= n;i++)
		y_spot[i] = 1.0/(1 + 16 * x_spot[i] * x_spot[i]);
	
	matrix[0][0] = 2.0;
	matrix[0][1] = 0.5;
	
	for(int i = 1;i < n;i++) {
		matrix[i][i - 1] = 0.5;
		matrix[i][i] = 2;
		matrix[i][i + 1] = 0.5;
	}
	
	matrix[n][n - 1] = 0.5;
	matrix[n][n] = 2;
	float f1, f2;
	f1 = -32.0 * (-5) / (1 + 32 * 25 + 256 * 625);
	f2 = -f1;
	
	for(int i = 1;i < n;i++) 
		matrix[i][n + 1] = 6.0 * (y_spot[i + 1] - y_spot[i] - y_spot[i] + y_spot[i - 1]) / (2 *h * h);
	
	matrix[0][n + 1] = 6.0 / h * ((y_spot[1] - y_spot[0]) / h - f1);
	matrix[n][n + 1] = 6.0 / h * (f2 - (y_spot[n] - y_spot[n - 1]) / h);
	
	print_matrix(n);
	
	for(int i = 0;i < n;i++) {
		float k1 = matrix[i][i];
		float k2 = matrix[i + 1][i];
		float mul = k2 / k1;
		matrix[i + 1][i] = 0;
		matrix[i + 1][i + 1] -= mul * matrix[i][i + 1];
		matrix[i + 1][n + 1] -= mul * matrix[i][n + 1];
	}
	
	print_matrix(n);
	
	for(int i = n;i > 0;i--) {
		float k1 = matrix[i][i];
		float k2 = matrix[i - 1][i];
		float mul = k2 / k1;
		matrix[i - 1][i] = 0;
		matrix[i - 1][i - 1] -= mul * matrix[i][i - 1];
		matrix[i - 1][n + 1] -= mul * matrix[i][n + 1];
	}
	
	print_matrix(n);
	
	for(int i = 0;i <= n;i++) {
		matrix[i][n + 1] /= matrix[i][i];
		matrix[i][i] = 1;
	}
	
	print_matrix(n);
	
	for(int i = 0;i <= n;i++)
		cout << matrix[i][n + 1] <<", " ;
	cout << endl;
	
	for(int i = 0;i <= n;i++)
		cout << y_spot[i] << endl;
	cout << endl;
}