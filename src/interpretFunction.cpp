#include "interpretFunction.h"

std::vector<Token> expr;
std::vector<std::complex<double>> evalStack;
double step = 0;

void setStep(double argstep) {
	step = argstep;
}

void initFunc(std::string infix) {
	std::vector<std::string> infixVec;
	std::stack<int> opStack;
	
	std::replace(infix.begin(), infix.end(), '(', '{');
	std::replace(infix.begin(), infix.end(), ')', '}');
    
	//Isolate each operator, function, number and variable as an individual token string
	for (size_t i = 0; i < infix.size(); i++) {
		if (infix[i] == '('){
			infix[i] = '{';
		}
		else if (infix[i] == ')'){
			infix[i] = '}';
		}
	    
		if (infix[i] == '\\') {
			int j = getOp(infix, i);
			std::string tempStr = infix.substr(i + 1, j - i - 1);
			infixVec.push_back(tempStr);
			i = j - 1;
		}
		else if (i > 0 && infix[i] == '-' && infix[i - 1] == '{') {
			infixVec.push_back("0");
			infixVec.push_back("-");
		}
		else if (i == 0 && infix[i] == '-') {
			infixVec.push_back("0");
			infixVec.push_back("-");
		}
		else if (infix.substr(i, 2) == "pi") {
			infixVec.push_back("pi");
			i++;
		}
		else if (std::isdigit(infix[i]) || infix[i] == '.') {
			int j = i;
			std::string tempStr = "";
			while (j < infix.size() && (std::isdigit(infix[j]) || infix[j] == '.')) {
			    tempStr += infix[j];
				j++;
			}
			if (j < infix.size() && infix[j] == 'i') {
				infixVec.push_back("(0,"+tempStr+")");
				j++;
			}
			else {
				infixVec.push_back(tempStr);
			}
			i = j - 1;
		}
		else if (infix[i] == '['){
			int j = i + 1;
			std::string tempStr = "(";
			while (infix[j] != ']'){
				tempStr += infix[j];
				j++;
			}
			infixVec.push_back(tempStr+")");
			i = j;
		}
		else {
			infixVec.push_back(infix.substr(i, 1));
		}
	}

	//Put expression into reverse polish notation using the shunting-yard algorithm with functions
	for (std::vector<std::string>::iterator it = infixVec.begin(); it != infixVec.end(); it++) {
		std::string token = *it;
		std::complex<double> value;
		int code = getOpCode(token);

		if (code == OTHER) {
			if (token == "z") {
				Token t = { VARIABLE, 0, -1 };
				expr.push_back(t);
			}
			else if (token == "pi") {
				Token t = { IMMEDIATE, PI, -1 };
				expr.push_back(t);
			}
			else if (token == "e") {
				Token t = { IMMEDIATE, E, -1 };
				expr.push_back(t);
			}
			else if (token == "i") {
				Token t = { IMMEDIATE, I, -1 };
				expr.push_back(t);
			}
			else {
				std::istringstream iss(token);
				iss >> value;
				expr.push_back({ IMMEDIATE, value, -1 });
			}
		}
		else if (code >= 5) {
			opStack.push(code);
		}
		else if (code >= 0 && code <= 4) {
			while (!opStack.empty() &&
				((opStack.top() >= code)
					|| ((code == ADD || code == SUB) && (opStack.top() == 0))
					|| ((code == MUL || code == DIV) && (opStack.top() == 2)))) {
				expr.push_back({ opStack.top() <= 4 ? BIN_OPERATOR : FUNCTION, 0.0, opStack.top() });
				opStack.pop();
			}
			opStack.push(code);
		}
		else if (code == LBRACKET) {
			opStack.push(code);
		}
		else if (code == RBRACKET) {
			while (!opStack.empty() && opStack.top() != LBRACKET) {
				expr.push_back({ opStack.top() <= 4 ? BIN_OPERATOR : FUNCTION, 0.0, opStack.top() });
				opStack.pop();
			}
			if (opStack.top() == LBRACKET) {
				opStack.pop();
			}
		}
	}
	while (!opStack.empty()) {
		expr.push_back({ opStack.top() <= 4 ? BIN_OPERATOR : FUNCTION, 0.0, opStack.top() });
		opStack.pop();
	}
	
	//pre-evaluates all constant valued sub-expressions in the postfox expression
	//this means that time is saved not re-evaluating constant expressions thousands if not millions of times
	//it also allows users to input constant expressions into the program without risk of slower processing.
	//constant expressions are often easier to write than the long decimal numerical form (ex. \cos(2) is much 
	//easier to input than −0.416146836547.
	bool changed = true;
	while(changed){
		changed = false;
		
		//iterator that lags behind by 1, needed for the rotations
		std::vector<Token>::iterator it = expr.begin();
		
		for (size_t i = 1; i < expr.size(); i++) {
			if(expr[i].type == FUNCTION && expr[i - 1].type == IMMEDIATE){
				expr[i - 1].num = evalFunc(it->op, expr[i - 1].num);
				
				//need to shift expression left by 1, to fill empty slot left by evaluating function
				if (i != expr.size() - 1){
				    std::rotate(it + 1, it + 2, expr.end());
				}
				expr.pop_back();
				changed = true;
			}
			else if(expr[i].type == BIN_OPERATOR && expr[i - 1].type == IMMEDIATE && expr[i - 2].type == IMMEDIATE){
			    switch (it->op) {
			        case ADD:
				        expr[i - 2].num = expr[i - 2].num + expr[i - 1].num;
				        break;
			        case SUB:
				        expr[i - 2].num = expr[i - 2].num - expr[i - 1].num;
				        break;
			        case MUL:
				        expr[i - 2].num = expr[i - 2].num * expr[i - 1].num;
				        break;  
			        case DIV:
				        expr[i - 2].num = expr[i - 2].num / expr[i - 1].num;
				        break;
			        case POW:
				        expr[i - 2].num = std::pow(expr[i - 2].num , expr[i - 1].num);
				        break;
			        default:
				        expr[i - 2].num = 0;
				        break;
			    }
				
				//need to shift expression left by 2, to fill empty slot by evaluating binary operator
				if (i != expr.size() - 1){
				    std::rotate(it, it + 2, expr.end());
				}
				expr.pop_back();
				expr.pop_back();
				changed = true;
			}
            it++;
		}
	}
	
	//initializes vector to necessary size that will be used as a stack for function evaluations
	evalStack.insert(evalStack.begin(), expr.size(), 0);
}

//gets operation/function/number indices that can be used to isolate token 
int getOp(std::string& infix, int n) {
	for (size_t i = n + 1; i < infix.size(); i++) {
		if (infix[i] == '\\' || infix[i] == '-' || infix[i] == '+' || infix[i] == '*' || infix[i] == '/' || infix[i] == '^' || infix[i] == '{' || infix[i] == '}') {
			if (i != n + 1) {
				return i;
			}
		}
	}
	return infix.size();
}

//Gets operation code
int getOpCode(std::string& token) {
	if (token == "{")
		return LBRACKET;
	else if (token == "}")
		return RBRACKET;
	else if (token == "+")
		return ADD;
	else if (token == "-")
		return SUB;
	else if (token == "*")
		return MUL;
	else if (token == "/")
		return DIV;
	else if (token == "^")
		return POW;
	else if (token == "re")
		return RE;
	else if (token == "im")
		return IM;
	else if (token == "abs")
		return ABS;
	else if (token == "arg")
		return ARG;
	else if (token == "conj")
		return CONJ;
	else if (token == "exp")
		return EXP;
	else if (token == "ln" || token == "log")
		return LOG;
	else if (token == "cos")
		return COS;
	else if (token == "sin")
		return SIN;
	else if (token == "tan")
		return TAN;
	else if (token == "sec")
		return SEC;
	else if (token == "csc")
		return CSC;
	else if (token == "cot")
		return COT;
	else if (token == "acos")
		return ACOS;
	else if (token == "asin")
		return ASIN;
	else if (token == "atan")
		return ATAN;
	else if (token == "cosh")
		return COSH;
	else if (token == "sinh")
		return SINH;
	else if (token == "tanh")
		return TANH;
	else if (token == "acosh")
		return ACOSH;
	else if (token == "asinh")
		return ASINH;
	else if (token == "atanh")
		return ATANH;
	else if (token == "step")
		return STEP;
	else if (token == "delta")
		return DELTA;
	else if (token == "gamma")
		return GAMMA;
	else if (token == "zeta")
		return ZETA;
	else
		return OTHER;
}

//evaluates a function based on an operation code that identifies function and argument
std::complex<double> evalFunc(int opCode, std::complex<double> z) {
	switch (opCode) {
	case RE:
		return z.real();
	case IM:
		return z.imag();
	case ABS:
		return std::abs(z);
	case ARG:
		return std::arg(z);
	case CONJ:
		return std::conj(z);
	case EXP:
		return std::exp(z);
	case LOG:
		return std::log(z);
	case COS:
		return std::cos(z);
	case SIN:
		return std::sin(z);
	case TAN:
		return std::tan(z);
	case SEC:
		return 1.0 / std::cos(z);
	case CSC:
		return 1.0 / std::sin(z);
	case COT:
		return 1.0 / std::tan(z);
	case ACOS:
		return std::acos(z);
	case ASIN:
		return std::asin(z);
	case ATAN:
		return std::atan(z);
	case COSH:
		return std::cosh(z);
	case SINH:
		return std::sinh(z);
	case TANH:
		return std::tanh(z);
	case ACOSH:
		return std::acosh(z);
	case ASINH:
		return std::asinh(z);
	case ATANH:
		return std::atanh(z);
	case STEP:
		return z.real() >= 0 ? 1 : 0;
	case DELTA:
		return std::abs(z.real()) <= 5*step;
	case GAMMA:
		return gamma(z);
	case ZETA:
		return zeta(z);
	default:
		return 0.0;
	}
}

//evaluates the overall expression using a stack
//this is often the bottleneck on the program speed
std::complex<double> f(std::complex<double> z) {
	int stackCounter = 0;
	std::complex<double> temp1, temp2;

	for (std::vector<Token>::iterator it = expr.begin(); it != expr.end(); it++) {
		if (it->type == IMMEDIATE) {
			evalStack[stackCounter++] = it->num;
		}
		else if (it->type == VARIABLE) {
			evalStack[stackCounter++] = z;
		}
		else if (it->type == FUNCTION) {
			temp1 = evalStack[--stackCounter];
			evalStack[stackCounter++] = evalFunc(it->op, temp1);
		}
		else if (it->type == BIN_OPERATOR) {
			temp1 = evalStack[--stackCounter];
			temp2 = evalStack[--stackCounter];

			switch (it->op) {
			case ADD:
				evalStack[stackCounter++] = temp2 + temp1;
				break;
			case SUB:
				evalStack[stackCounter++] = temp2 - temp1;
				break;
			case MUL:
				evalStack[stackCounter++] = temp2 * temp1;
				break;
			case DIV:
				evalStack[stackCounter++] = temp2 / temp1;
				break;
			case POW:
				evalStack[stackCounter++] = std::pow(temp2, temp1);
				break;
			default:
				evalStack[stackCounter++] = 0;
				break;
			}
		}
	}
	return evalStack[--stackCounter];
}

//uses Lanczos approximation
std::complex<double> gamma(std::complex<double> z) {
	std::complex<double> y, x, t;
	if (z.real() < 0.5) {
		y = PI / (std::sin(PI * z) * gamma(1.0 - z)); //reflection formula
	}
	else {
		z -= 1.0;
		x = p[0];
		for (size_t i = 1; i < pSize; i++) {
			x += p[i] / (z + (double)(i - 1) + 1.0);
		}
		t = z + (double)pSize - 0.5 - 1.0;
		y = std::sqrt(2 * PI) * std::pow(t, z + 0.5) * std::exp(-t) * x;
	}

	return y;
}

//using formula 21 on https://mathworld.wolfram.com/RiemannZetaFunction.html
std::complex<double> zeta(std::complex<double> z) {
	if (z.real() < 0) {
		return 2.0 * std::pow(2.0 * PI, z - 1.0) * std::sin(0.5 * PI * z) * gamma(1.0 - z) * zeta(1.0 - z);
	} 
	else if (z.real() > 0.25) {
		std::complex<double> tot = 0, partialTot = 1.0;
		double div = -1, error = 0.00001, N = std::pow(error, -1 / z.real());
		if (N <= 0.5 * std::max(460.0, 4 * z.imag() * z.imag() + 6 * z.imag())) {
			for (double n = 1; n < 20; n++) {
				div *= -1;
				partialTot = div * std::pow(n, -z);
				tot += partialTot;
			}

			return tot / (1.0 - std::pow(2.0, 1.0 - z));
		}
		else {
			int binCoeff = 1, n = 0;
			std::complex<double> tot = 0, partialTot = 1.0;
			double div = 0.5;

			while (n <= std::max(20.0, 2 * z.imag())) {
				partialTot = 1.0;
				binCoeff = 1;
				for (int k = 1; k <= n; k++) {
					binCoeff = binCoeff * (k - n - 1) / k;
					partialTot += ((double)binCoeff) * std::pow(k + 1.0, -z);
				}
				partialTot *= div;
				tot += partialTot;
				div *= 0.5;
				n++;
			}
	
			return tot / (1.0 - std::pow(2.0, 1.0 - z));
		}
	}
	else {
		std::complex<double> a = 2.0 * std::pow(2.0 * PI, z - 1.0) * std::sin(0.5 * PI * z) * gamma(1.0 - z) * zeta(1.0 - z);
		std::complex<double> b = a;
		std::complex<double> tot = 0, partialTot = 1.0;
		double div = -1, error = 0.00001, N = std::pow(error, -1 / z.real());
		if (N <= 0.5 * std::max(460.0, 4 * z.imag() * z.imag() + 6 * z.imag())) {
			for (double n = 1; n < 20; n++) {
				div *= -1;
				partialTot = div * std::pow(n, -z);
				tot += partialTot;
			}

			b = tot / (1.0 - std::pow(2.0, 1.0 - z));
		}
		else {
			int binCoeff = 1, n = 0;
			std::complex<double> tot = 0, partialTot = 1.0;
			double div = 0.5;

			while (n <= std::max(20.0, 2 * z.imag())) {
				partialTot = 1.0;
				binCoeff = 1;
				for (int k = 1; k <= n; k++) {
					binCoeff = binCoeff * (k - n - 1) / k;
					partialTot += ((double)binCoeff) * std::pow(k + 1.0, -z);
				}
				partialTot *= div;
				tot += partialTot;
				div *= 0.5;
				n++;
			}

			b = tot / (1.0 - std::pow(2.0, 1.0 - z));
		}

		return 4.0 * (z.real() * a + (0.25 - z.real()) * b);
	}
}
