#pragma once
#include <complex>
#include <string.h>
#include <vector>
#include <stack>
#include <sstream>
#include <algorithm>
#include <iostream>

const double PI = 3.14159265358979323846;
const double E = 2.71828182845904523536;
const std::complex<double> I = std::complex<double>(0, 1.0);

//needed for gamma function by lanczos approx
const double p[9] = { 0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
const int pSize = 9;

//struct used for creating postfix expression, since the expression
//needs to handle operators, functions and numbers. I could just use
//a string to represent all of them, and convert to correct format 
//when necessary, but that would be a significant tax on program speed
//this approach uses more memory but will be faster when evaluating 
//the function expression
struct Token {
	int type = 0;
	std::complex<double> num;//type = 0: constant, type = 1: variable (z)
	int op;//type = 2: function, type = 3: operator
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
