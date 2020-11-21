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
<<<<<<< HEAD
	
	std::replace(infix.begin(), infix.end(), '(', '{');
	std::replace(infix.begin(), infix.end(), ')', '}');
	
	//Isolate each operator, function, number and variable as an individual token string
	for (int i = 0; i < infix.size(); i++) {
=======

	//Isolate each operator, function and number as an individual token string
	for (int i = 0; i < infix.size(); i++) {
		if (infix[i] == '('){
			infix[i] = '{';
		}
		else if (infix[i] == ')'){
			infix[i] = '}';
		}
	    
>>>>>>> 143c7c45e707cb67173091374fab8be96c0284c1
		if (infix[i] == '\\') {
			int j = getOp(infix, i);
			std::string tempStr = infix.substr(i + 1, j - i - 1);
			infixVec.push_back(tempStr);
			i = j - 1;
		}
		else if (i > 0 && infix[i] == '-' && infix[i-1] == '{') {
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
<<<<<<< HEAD
				tempStr += infix[j];
=======
			    tempStr += infix[j];
>>>>>>> 143c7c45e707cb67173091374fab8be96c0284c1
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

		if (code == -1) {
			if (token == "z") {
				Token t = { 1, 0, -1 };
				expr.push_back(t);
			}
			else if (token == "pi") {
				Token t = { 0, PI, -1 };
				expr.push_back(t);
			}
			else if (token == "e") {
				Token t = { 0, E, -1 };
				expr.push_back(t);
			}
			else if (token == "i") {
				Token t = { 0, I, -1 };
				expr.push_back(t);
			}
			else {
				std::istringstream iss(token);
				iss >> value;
				expr.push_back({ 0, value, -1 });
			}
		}
		else if (code >= 5) {
			opStack.push(code);
		}
		else if (code >= 0 && code <= 4) {
			while (!opStack.empty() &&
				((opStack.top() >= code)
					|| ((code == 0 || code == 1) && (opStack.top() == 0))
					|| ((code == 2 || code == 3) && (opStack.top() == 2)))) {
				expr.push_back({ opStack.top() <= 4 ? 3 : 2, 0.0, opStack.top() });
				opStack.pop();
			}
			opStack.push(code);
		}
		else if (code == -3) {
			opStack.push(code);
		}
		else if (code == -2) {
			while (!opStack.empty() && opStack.top() != -3) {
				expr.push_back({ opStack.top() <= 4 ? 3 : 2, 0.0, opStack.top() });
				opStack.pop();
			}
			if (opStack.top() == -3) {
				opStack.pop();
			}
		}
	}
	while (!opStack.empty()) {
		expr.push_back({ opStack.top() <= 4 ? 3 : 2, 0.0, opStack.top() });
		opStack.pop();
	}
	
	//pre-evaluates all constant sub-expressions in the postfox expression
	//this means that time is saved not re-evaluating constant expressions thousands if not millions of times
	//it also allows users to input constant expressions into the program without any risk of slower processing
	//constant expressions are often easier to write than the long decimal numerical form (ex. \cos(2) is much 
	//easier to input than −0.416146836547.
	bool changed = true;
	while(changed){
		changed = false;
		std::vector<Token>::iterator it = expr.begin()+1;
		for (; it != expr.end(); it++) {
			if(it->type == 2 && (it - 1)->type == 0){
				(it - 1)->num = evalFunc(it->op, (it - 1)->num);
				std::rotate(it, std::next(it), expr.end());
				expr.pop_back();
				changed = true;
			}
			else if(it->type == 3 && (it - 1)->type == 0 && (it - 2)->type == 0){
				switch (it->op) {
				case 0:
					(it - 2)->num = (it - 2)->num + (it - 1)->num;
					break;
				case 1:
					(it - 2)->num = (it - 2)->num - (it - 1)->num;
					break;
				case 2:
					(it - 2)->num = (it - 2)->num * (it - 1)->num;
					break;
				case 3:
					(it - 2)->num = (it - 2)->num / (it - 1)->num;
					break;
				case 4:
					(it - 2)->num = std::pow((it - 2)->num , (it - 1)->num);
					break;
				default:
					(it - 2)->num = 0;
					break;
				}
				
				std::rotate(std::prev(it), std::next(it), expr.end());
				expr.pop_back();
				expr.pop_back();
				changed = true;
			}
		}
	}
	
	//initializes vector to necessary size that will be used as a stack for function evaluations
	evalStack.insert(evalStack.begin(), expr.size(), 0);
}

//gets operation/function/number indices that can be used to isolate token 
int getOp(std::string& infix, int n) {
	for (int i = n + 1; i < infix.size(); i++) {
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
		return -3;
	else if (token == "}")
		return -2;
	else if (token == "+")
		return 0;
	else if (token == "-")
		return 1;
	else if (token == "*")
		return 2;
	else if (token == "/")
		return 3;
	else if (token == "^")
		return 4;
	else if (token == "Re" || token == "re")
		return 5;
	else if (token == "Im" || token == "im")
		return 6;
	else if (token == "abs")
		return 7;
	else if (token == "arg")
		return 8;
	else if (token == "conj")
		return 9;
	else if (token == "cos")
		return 10;
	else if (token == "sin")
		return 11;
	else if (token == "tan")
		return 12;
	else if (token == "sec")
		return 13;
	else if (token == "csc")
		return 14;
	else if (token == "cot")
		return 15;
	else if (token == "acos")
		return 16;
	else if (token == "asin")
		return 17;
	else if (token == "atan")
		return 18;
	else if (token == "asec")
		return 19;
	else if (token == "acsc")
		return 20;
	else if (token == "acot")
		return 21;
	else if (token == "cosh")
		return 22;
	else if (token == "sinh")
		return 23;
	else if (token == "tanh")
		return 24;
	else if (token == "sech")
		return 25;
	else if (token == "csch")
		return 26;
	else if (token == "coth")
		return 27;
	else if (token == "acosh")
		return 28;
	else if (token == "asinh")
		return 29;
	else if (token == "atanh")
		return 30;
	else if (token == "asech")
		return 31;
	else if (token == "acsch")
		return 32;
	else if (token == "acoth")
		return 33;
	else if (token == "exp")
		return 34;
	else if (token == "ln")
		return 35;
	else if (token == "log")
		return 36;
	else if (token == "step")
		return 37;
	else if (token == "delta" || token == "Delta")
		return 38;
	else if (token == "gamma" || token == "Gamma")
		return 39;
	else if (token == "zeta")
		return 40;
	else
		return -1;
}

//evaluates a function based on an operation code that identifies function and argument
std::complex<double> evalFunc(int opCode, std::complex<double> z) {
	switch (opCode) {
	case 5:
		return z.real();
	case 6:
		return z.imag();
	case 7:
		return std::abs(z);
	case 8:
		return std::arg(z);
	case 9:
		return std::conj(z);
	case 10:
		return std::cos(z);
	case 11:
		return std::sin(z);
	case 12:
		return std::tan(z);
	case 13:
		return 1.0 / std::cos(z);
	case 14:
		return 1.0 / std::sin(z);
	case 15:
		return 1.0 / std::tan(z);
	case 16:
		return std::acos(z);
	case 17:
		return std::asin(z);
	case 18:
		return std::atan(z);
	case 19:
		return std::acos(1.0 / z);
	case 20:
		return std::asin(1.0 / z);
	case 21:
		return std::atan(1.0 / z);
	case 22:
		return std::cosh(z);
	case 23:
		return std::sinh(z);
	case 24:
		return std::tanh(z);
	case 25:
		return 1.0 / std::cosh(z);
	case 26:
		return 1.0 / std::sinh(z);
	case 27:
		return 1.0 / std::tanh(z);
	case 28:
		return std::acosh(z);
	case 29:
		return std::asinh(z);
	case 30:
		return std::atanh(z);
	case 31:
		return std::acosh(1.0 / z);
	case 32:
		return std::asinh(1.0 / z);
	case 33:
		return std::atanh(1.0 / z);
	case 34:
		return std::exp(z);
	case 35:
		return std::log(z);
	case 36:
		return std::log10(z);
	case 37:
		return z.real() >= 0 ? 1 : 0;
	case 38:
		return std::abs(z.real()) <= 5*step;
	case 39:
		return gamma(z);
	case 40:
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
		if (it->type == 0) {
			evalStack[stackCounter++] = it->num;
		}
		else if (it->type == 1) {
			evalStack[stackCounter++] = z;
		}
		else if (it->type == 2) {
			temp1 = evalStack[--stackCounter];
			evalStack[stackCounter++] = evalFunc(it->op, temp1);
		}
		else if (it->type == 3) {
			temp1 = evalStack[--stackCounter];
			temp2 = evalStack[--stackCounter];

			switch (it->op) {
			case 0:
				evalStack[stackCounter++] = temp2 + temp1;
				break;
			case 1:
				evalStack[stackCounter++] = temp2 - temp1;
				break;
			case 2:
				evalStack[stackCounter++] = temp2 * temp1;
				break;
			case 3:
				evalStack[stackCounter++] = temp2 / temp1;
				break;
			case 4:
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
		for (int i = 1; i < pSize; i++) {
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
