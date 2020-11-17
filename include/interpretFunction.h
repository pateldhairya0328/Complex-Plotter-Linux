#pragma once
#include <complex>
#include <string.h>
#include <vector>
#include <queue>
#include <stack>
#include <sstream>
#include <algorithm>

constexpr double PI = 3.14159265358979323846;
constexpr double E = 2.71828182845904523536;
constexpr std::complex<double> I = std::complex<double>(0, 1.0);

//needed for gamma function by lanczos approx
constexpr double p[9] = { 0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
constexpr int pSize = 9;

//struct used for creating postfix expression, since the expression
//needs to handle operators, functions and numbers. I could just use
//a string to represent all of them, and convert to correct format 
//when necessary, but that would be a significant tax on program speed
//this approach uses more memory but will be faster when evaluating 
//the function expression
struct Token {
	int type = 0;
	std::complex<double> num;//type = 0: constant, type1: variable (z)
	int op;//type = 2: function, type = 3: operator

	/*
	op codes:
	-3 - (
	-2 - )

	-1 - not an operator or function or ()

	0  - (+)
	1  - (-)
	2  - (*)
	3  - (/)
	4  - (^)

	5  - Re(z)
	6  - Im(z)
	7  - abs(z)
	8  - arg(z)
	9  - conj(z)

	10  - cos(z)
	11 - sin(z)
	12 - tan(z)
	13 - sec(z)
	14 - csc(z)
	15 - cot(z)

	16 - acos(z)
	17 - asin(z)
	18 - atan(z)
	19 - asec(z)
	20 - acsc(z)
	21 - acot(z)

	22 - cosh(z)
	23 - sinh(z)
	24 - tanh(z)
	25 - sech(z)
	26 - csch(z)
	27 - coth(z)

	28 - acosh(z)
	29 - asinh(z)
	30 - atanh(z)
	31 - asech(z)
	32 - acsch(z)
	33 - acoth(z)

	34 - exp(z)
	35 - ln(z)
	36 - log(z) (base 10)

	37 - step(|z|) (step greater than equal to 0)
	38 - delta(|z|) (1 around |z|, 0 elsewhere)

	39 - gamma(z)
	*/
};

void setStep(double argstep);
void initFunc(std::string infix);
int getOpCode(std::string& token);
int getOp(std::string& infix, int n); 
std::complex<double> f(std::complex<double> z);
std::complex<double> evalFunc(int opCode, std::complex<double> z);
std::complex<double> gamma(std::complex<double> z);
std::complex<double> zeta(std::complex<double> z);
//std::complex<double> bessel_J(int alpha, std::complex<double> z);
