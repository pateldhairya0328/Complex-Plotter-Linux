#pragma once
#include <complex>
#include <string.h>
#include <vector>
#include <stack>
#include <sstream>
#include <iostream>
#include <algorithm>

#define PI 3.14159265358979323846
#define E 2.71828182845904523536

const std::complex<double> I = std::complex<double>(0, 1.0);

enum otherTokens { LBRACKET = -3, RBRACKET = -2, OTHER = -1};
enum operation { ADD, SUB, MUL, DIV, POW, RE, IM, ABS, ARG, CONJ, EXP, LOG, COS, SIN, TAN, SEC, CSC, COT, ACOS, ASIN, ATAN, COSH, SINH, TANH, ACOSH, ASINH, ATANH, STEP, DELTA, GAMMA, ZETA, DIGAMMA, AIRY, BIRY, BESSELJ, BESSELY, SPHBESSELJ, SPHBESSELY };
enum tokenType { FUNCTION = -3, BIN_OPERATOR = -2, IMMEDIATE = -1, VARIABLE = 0 };

//struct used for creating postfix expression, since the expression
//needs to handle operators, functions and numbers. I could just use
//a string to represent all of them, and convert to correct format 
//when necessary, but that would be a significant tax on program speed
//this approach uses more memory but will be faster when evaluating 
//the function expression
struct Token {
    int type = 0;
    std::complex<double> num;
    int op;
};

void setStep(double argstep);
void initFunc(std::string infix);
int getOpCode(std::string& token);
int getOp(std::string& infix, int n); 
std::complex<double> f(std::complex<double> z);
std::complex<double> evalFunc(int opCode, std::complex<double> z);
std::complex<double> evalFuncTwoArg(int opcode, std::complex<double> arg1, std::complex<double> arg2);
