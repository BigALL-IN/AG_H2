#include "functions.h"
#include <iostream>

double Rastrigin(std::vector<double> v) {
	double res = 10 * v.size();
	for (int i = 0; i < v.size(); i++) {
		res += v[i] * v[i] - 10 * cos(2 * M_PI * v[i]);
	}
	return res;
}

double Michalewicz(std::vector<double> v) {
	double res = 0;
	for (int i = 0; i < v.size(); i++) {
		res += sin(v[i]) * pow(sin(((i + 1) * v[i] * v[i]) / M_PI), 20);
	}
	return -res;
}

double Sphere(std::vector<double> v) {
	double res = 0;
	for (int i = 0; i < v.size(); i++) {
		res += v[i] * v[i];
	}
	return res;
}

double Sumsq(std::vector<double> v) {
	double res = 0;
	for (int i = 0; i < v.size(); i++) {
		res += (i + 1) * v[i] * v[i];
	}
	return res;
}