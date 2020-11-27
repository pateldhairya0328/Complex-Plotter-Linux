#pragma once
#include <complex>
#include "complex_bessel.h"
#define PI 3.14159265358979323846
#define EPS 2.2204460493e-16

//size of array needed for lanczos approx
const int pSize = 9;
//array needed for gamma function by lanczos approx
const double p[9] = { 0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };

std::complex<double> gamma(std::complex<double> z);
std::complex<double> digamma(std::complex<double> z);
std::complex<double> zeta(std::complex<double> z);
std::complex<double> airy(std::complex<double> z);
std::complex<double> biry(std::complex<double> z);
std::complex<double> cylBesselJ(double nu, std::complex<double> z);
std::complex<double> cylBesselY(double nu, std::complex<double> z);
std::complex<double> sphBesselJ(double nu, std::complex<double> z);
std::complex<double> sphBesselY(double nu, std::complex<double> z);
