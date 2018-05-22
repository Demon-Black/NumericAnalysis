#include <iostream>
#include <cmath>
#include <cstdio>
using namespace std;

float ln2 = 0.693147190546;
float epsilon = 0.00005;

int main(int argc, char** argv) {
	float f = 0.0;
	for (float i = 1.0;;i += 1.0) {
		if (int(i) % 2)
			f += (1 / i);
		else
			f -= (1 / i);
		if (abs(f - ln2) < epsilon) {
			printf("%d\n", int(i));
			break;
		}
	}
	return 0;
}