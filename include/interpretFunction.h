#pragma once
#include <complex>
#include <string.h>
#include <vector>
#include <stack>
#include <sstream>
#include <iostream>
#include <algorithm>

#define PI 3.14159265358979323846 // pi
#define E 2.71828182845904523536 // Euler's number
#define EPS 2.2204460493e-16 // machine epsilon for double precision, to help the optimal step for basic central difference numerical derivative

const std::complex<double> I = std::complex<double>(0, 1.0);

enum otherTokens { LBRACKET = -3, RBRACKET = -2, OTHER = -1};
enum operation { ADD, SUB, MUL, DIV, POW, RE, IM, ABS, ARG, CONJ, EXP, LOG, COS, SIN, TAN, SEC, CSC, COT, ACOS, ASIN, ATAN, COSH, SINH, TANH, ACOSH, ASINH, ATANH, STEP, DELTA, GAMMA, ZETA, DIGAMMA, AIRY, BIRY, BESSELJ, BESSELY, SPHBESSELJ, SPHBESSELY };
enum tokenType { DIFF = -5, SUM = -4, FUNCTION = -3, BIN_OPERATOR = -2, IMMEDIATE = -1, VARIABLE = 0 };

//struct used for creating postfix expression, since the expression
//needs to handle operators, functions and numbers. I could just use
//a string to represent all of them, and convert to correct format 
//when necessary, but that would be a significant tax on program speed
//this approach uses more memory but will be faster when evaluating 
//the function expression
struct Token {
    int type = 0;
    std::complex<double> num; //if type >= -1
    int op; //operation if tokenType = FUNCTION or BIN_OPERATOR, variable name if tokenType = SUM, order of derivative if tokenType = DIFF
    std::vector<int> coeffs; //bounds if tokenType = SUM, finite difference derivative coeffs f tokenType = DIFF
    std::vector<Token> postfix; //if tokenType = SUM
};

void setStep(double argstep);
void initFunc(std::string &infix);
std::vector<Token> parseFunc(std::string &infix);
std::vector<std::string> getInfixVec(std::string &infix);
int getOpCode(std::string &token);
int getOp(std::string &infix, int n); 
std::complex<double> f(std::complex<double> z);
std::complex<double> f(std::vector<Token> &postfix);
std::complex<double> evalFunc(int opCode, std::complex<double> z);
std::complex<double> evalFuncTwoArg(int opcode, std::complex<double> arg1, std::complex<double> arg2);
