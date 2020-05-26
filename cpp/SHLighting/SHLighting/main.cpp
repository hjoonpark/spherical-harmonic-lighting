#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>

using std::cout;
using std::endl;
#define PI 3.14159265359
struct SHSample {
	Eigen::Vector3d sph;
	Eigen::Vector3d vec;
	double* coeff;

	SHSample(const int n_coeff) {
		coeff = new double[n_coeff];
	}
};
double P(int l, int m, double x)
{
	// evaluate an Associated Legendre Polynomial P(l,m,x) at x
	double pmm = 1.0;
	if (m > 0) {
		double somx2 = sqrt((1.0 - x) * (1.0 + x));
		double fact = 1.0;
		for (int i = 1; i <= m; i++) {
			pmm *= (-fact) * somx2;
			fact += 2.0;
		}
	}
	if (l == m) {
		return pmm;
	}

	double pmmp1 = x * (2.0 * m + 1.0) * pmm;
	if (l == m + 1) {
		return pmmp1;
	}

	double pll = 0.0;
	for (int ll = m + 2; ll <= l; ++ll) {
		pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
		pmm = pmmp1;
		pmmp1 = pll;
	}
	return pll;
}
int factorial(int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}
double K(int l, int m)
{
	// renormalisation constant for SH function
	int factorial1 = factorial(l - m);
	int factorial2 = factorial(l + m);
	double temp = ((2.0 * l + 1.0) * factorial1) / (4.0 * PI * factorial2);
	return sqrt(temp);
}
double SH(int l, int m, double theta, double phi)
{
	// return a point sample of a Spherical Harmonic basis function
	// l is the band, range [0..N]
	// m in the range [-l..l]
	// theta in the range [0..Pi]
	// phi in the range [0..2*Pi]
	const double sqrt2 = sqrt(2.0);
	if (m == 0) {
		double p = P(l, m, cos(theta));
		double k = K(l, 0);
		return k * p;
	}
	else if (m > 0) {
		double k = K(l, m);
		double p = P(l, m, cos(theta));
		return sqrt2 * k * cos(m * phi) * p;
	}
	else {
		double k = K(l, -m);
		double p = P(l, -m, cos(theta));
		return sqrt2 * k * sin(-m * phi) * p;
	}
}
void SH_setup_spherical_samples(std::vector<SHSample> &samples, int sqrt_n_samples, int n_bands, int n_coeff)
{
	// fill an N*N*2 array with uniformly distributed
	// samples across the sphere using jittered stratification
	int i = 0; // array index
	double oneoverN = 1.0 / sqrt_n_samples;
	for (int a = 0; a < sqrt_n_samples; a++) {
		for (int b = 0; b < sqrt_n_samples; b++) {
			// generate unbiased distribution of spherical coords
			double randx = ((double)rand() / (RAND_MAX));
			double randy = ((double)rand() / (RAND_MAX));

			double x = (a + randx) * oneoverN; // do not reuse results
			double y = (b + randy) * oneoverN; // each sample must be random
			if (randx > 1.0 || randx < 0.0 || randy > 1.0 || randy < 0.0) {
				cout << "break" << endl;
				std::getchar();
			}
			double theta = 2.0 * acos(sqrt(1.0 - x));
			double phi = 2.0 * PI * y;
			SHSample sample(n_coeff);
			sample.sph = Eigen::Vector3d(theta, phi, 1.0);
			// convert spherical coords to unit vector
			sample.vec = Eigen::Vector3d(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
			// precompute all SH coefficients for this sample
			for (int l = 0; l < n_bands; ++l) {
				for (int m = -l; m <= l; ++m) {
					int index = l * (l + 1) + m;
					sample.coeff[index] = SH(l, m, theta, phi);
				}
			}
			samples.push_back(sample);
			++i;
		}
	}
	printf("%d == %d\n", samples.size(), i);
}
typedef double (*SH_polar_fn)(double theta, double phi);
void SH_project_polar_function(SH_polar_fn fn, const std::vector<SHSample> &samples, const int n_samples, const int n_coeff, double result[])
{
	const double weight = 4.0 * PI;
	// for each sample
	for (int i = 0; i < n_samples; ++i) {
		double theta = samples[i].sph[0];
		double phi = samples[i].sph[1];
		for (int n = 0; n < n_coeff; ++n) {
			result[n] += fn(theta, phi) * samples[i].coeff[n];
		}
	}
	// divide the result by weight and number of samples
	double factor = weight / n_samples;
	for (int i = 0; i < n_coeff; ++i) {
		result[i] = result[i] * factor;
	}
}
double light_fn(double theta, double phi) {
	double L = std::max(0.0, 5.0 * cos(theta) - 4.0) + std::max(0.0, -4.0 * sin(theta - PI) * cos(phi - 2.5) - 3.0);
	return L;
}
int main() {
	int n_bands = 4;
	int n_samples = 200;
	int n_coeff = (n_bands - 1) * (n_bands - 1 + 1) + n_bands;
	printf("%d, %d, %d\n", n_bands, n_samples, n_coeff);
	double* results = new double[n_coeff];
	for (int i = 0; i < n_coeff; i++) {
		results[i] = 0.0;
	}
	cout << "1" << endl;
	std::vector<SHSample> samples;
	samples.reserve(n_samples);

	cout << "2" << endl;
	SH_setup_spherical_samples(samples, n_samples, n_bands, n_coeff);
	cout << "3" << endl;
	SH_project_polar_function(light_fn, samples, n_samples, n_coeff, results);
	cout << "4" << endl;

	for (int l = 0; l < n_bands; l++) {
		for (int m = -l; m <= l; m++) {
			int i = l * (l + 1) + m;
			printf("\t%f", results[i]);
		}
		printf("\n");
	}

	delete[] results;
	return 0;
}