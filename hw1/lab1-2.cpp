#include <cstdio>
#include <cstring>
#define function(x) (1 / (16 * (x) * (x) + 1))
using namespace std;

double x1[11], x2[21];

void Lagrange(double* point,int number) {
	double c[number];
	memset(c, 0.0, sizeof(double) * (number));
	for (int i = 0;i < number;i++) {
		double y = function(point[i]);
		double d = 1.0;
		for (int j = 0;j < number;j++)
			if (i != j)
				d *= (point[i] - point[j]);
		double p[number], q[number];
		memset(p, 0.0, sizeof(double) * (number));
		memset(q, 0.0, sizeof(double) * (number));
		q[number - 1] = 1.0;
		for (int j = 0;j < number;j++) {
			if (i != j) {
				for (int k = number - 1;k > 0;k--)
					p[k - 1] = q[k];
				p[number - 1] = 0.0;
				for (int k = number - 1;k >= 0;k--)
					p[k] -= point[j] * q[k];
				for (int k = number - 1;k >= 0;k--)
					q[k] = p[k];
			}
		}

		for (int k = 0;k < number;k++)
			c[k] += q[k] * y / d;
	}
	for (int i = 0;i < number;i++)
		printf("%.8lf\n", c[i]);
}

int main(int argc, char** argv) {
	for (int i = 0;i < 11;i++)
		x1[i] = -5.0 + i * 1.0;
	for (int i = 0;i < 21;i++)
		x2[i] = -5.0 + i * 0.5;
	Lagrange(x2, 21);
	return 0;
}