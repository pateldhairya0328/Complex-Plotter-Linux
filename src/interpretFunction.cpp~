#include "interpretFunction.h"
#include "specialFunctions.h"

std::vector<Token> expr;
std::stack<std::complex<double>> evalStack;

std::vector<std::complex<double>> variables;
std::vector<std::string> variableNames;

double step = 0;

void setStep(double argstep) {
    step = argstep;
}

// Initialize the user inputed function
// Initializes the z variable which is varied through the execution of the program
// Uses parseFunc to parse the entire input expression
void initFunc(std::string infix){
    variables.push_back(0.0);
    variableNames.push_back("z");
    expr = parseFunc(infix);
}

// Parse function by taking an infix string, and returning a vector of tokens that reprsents the
// postfix/ reverse polish notation equivalent expression of the infix
std::vector<Token> parseFunc(std::string infix) {
    //convert infix string to a vector of strings of meaningful string tokens
    std::vector<std::string> infixVec = getInfixVec(infix);
    
    std::stack<int> opStack;
    std::vector<Token> emptyVec = {};
    std::vector<Token> postfix;
    
    //Put expression into reverse polish notation using the shunting-yard algorithm with functions
    for (std::vector<std::string>::iterator it = infixVec.begin(); it != infixVec.end(); it++) {
        std::string token = *it;
        std::complex<double> value;
        int code = getOpCode(token);
        
        // The other code is meant to indicate a variable or immediate quantity (such as z, i, pi, e, or 
        // some user input number). Derivatives and sums are also treated as one, because during execution
        // once the sum or derivative is evaluated, it will simply be replaced by a number in the 
        // evaluation stack. 
        if (code == OTHER) {
            if (token.substr(0, 3) == "sum"){
                std::vector<std::string> args[4];
                std::string argString = token.substr(4, token.size() - 5);// isolates string of arguments
                int lowerBound, upperBound;
                std::string subInfix;
                std::string sumIndexVar;
                size_t k = 0, prev = 0;
                
                // extracts the four arguments. it is done the long way with a for loop as the last argument
                // itself can have a nested sum with its own semicolons, so simply string splitting on 
                // semicolons does not work
                for (size_t i = 0; i < argString.size() && k < 3; i++){
                    if (argString[i] == ';'){
                        switch(k){
                            case 0: // first argument is name of index variable
                                sumIndexVar = argString.substr(0, i);
                                prev = i + 1;
                                k = 1;
                                break;
                            case 1: // second argument is lower bound of sum
                                lowerBound = std::stoi(argString.substr(prev, i - prev));
                                prev = i + 1;
                                k = 2;
                                break;
                            case 2: // third argument is upper bound of sum
                                upperBound = std::stoi(argString.substr(prev, i - prev));
                                prev = i + 1;
                                // fourth argument is the expression inside the sum
                                subInfix = argString.substr(prev, argString.size() - prev);
                                k = 4;
                                break;
                            default:
                                break;
                        }
                    }
                }
                
                // initializes the variable in the vector of variables, if it does not already exist
                bool varExists = false;
                int varIndex = -1;
                for (int j = 0; j < variables.size(); j++){
                    if (variableNames[j] == sumIndexVar){
                        varExists = true;
                        varIndex = j; // set variable to the one that already exists
                        break;
                    }
                }
                if (!varExists){
                    variableNames.push_back(sumIndexVar);
                    variables.push_back(0.0);
                    varIndex = (int)(variables.size() - 1); // varIndex is the last one (newest entry)
                }
                
                // parse expression inside the sum, which will be evaluated when the program is running
                std::vector<Token> subExpr = parseFunc(subInfix);
                
                // make the token and push it into the postfix expression
                Token t = { SUM, 0, varIndex, {lowerBound, upperBound}, subExpr };
                postfix.push_back(t);
            }
            else if (token.substr(0, 4) == "diff") {
                std::string argString = token.substr(5, token.size() - 6);// isolates string of arguments
                size_t k = 0, j = 0;
                // only three arguments; finds the first semicolon
                while (argString[k] != ';' && k < argString.size()){
                    k++;
                }
                // first argument is the variable to do derivative wrt
                std::string varName = argString.substr(0, k);
                
                j = k + 1;
                while (argString[j] != ';' && k < argString.size()){
                    j++;
                }
                // second argument is the order of the derivative
                int order = std::stoi(argString.substr(k + 1, j - k - 1));
                
                std::vector<int> coeffs = getFiniteDiffCoeffs(order);                
                coeffs.push_back(0); //derivative along real axis
                if (varName == "y"){
                    coeffs[coeffs.size() - 1] = -1; //derivative along imaginary axis
                }
                
                // get the third and last argument, and parse it into a postfix expression
                std::string subInfix = argString.substr(j + 1, argString.size() - j - 1);
                std::vector<Token> subExpr = parseFunc(subInfix);
                // make the token and add it to the main postfix expression
                Token t = { DIFF, 0, order, coeffs, subExpr };
                postfix.push_back(t);
            }
            //check if token was one of the three specific immediates
            else if (token == "pi") {
                Token t = { IMMEDIATE, PI, -1, {}, emptyVec };
                postfix.push_back(t);
            }
            else if (token == "e") {
                Token t = { IMMEDIATE, E, -1, {}, emptyVec };
                postfix.push_back(t);
            }
            else if (token == "i") {
                Token t = { IMMEDIATE, I, -1, {}, emptyVec };
                postfix.push_back(t);
            }
            else {
                bool tokenAdded = false;
                //check if token is a variable, and add it
                for (int i = 0; i < variables.size(); i++){
                    if (token == variableNames[i]){
                        Token t = { VARIABLE + i, 0, -1, {}, emptyVec };
                        postfix.push_back(t);
                        tokenAdded = true;
                        break;
                    }
                }
                
                //if token was not a variable, it must have been an immediate, so add that instead
                if (!tokenAdded){
                    std::istringstream iss(token);
                    iss >> value;
                    postfix.push_back({ IMMEDIATE, value, -1, {}, emptyVec });
                }
            }
        }
        // it is a function, which we also can treat like an immediate, as it will evaluate to one value and result in
        // pushing a single constant onto the stack at the end, so we instantly push it onto the stack like an immediate
        else if (code >= 5) {
            opStack.push(code);
        }
        // it is a binary operator, so we deal with it based on operator precendence
        else if (code >= 0 && code <= 4) {
            while (!opStack.empty() &&
                ((opStack.top() >= code)
                    || ((code == ADD || code == SUB) && (opStack.top() == 0))
                    || ((code == MUL || code == DIV) && (opStack.top() == 2)))) {
                postfix.push_back({ opStack.top() <= 4 ? BIN_OPERATOR : FUNCTION, 0.0, opStack.top(), {}, emptyVec });
                opStack.pop();
            }
            opStack.push(code);
        }
        // push left bracket
        else if (code == LBRACKET) {
            opStack.push(code);
        }
        // if we see a right bracket, keep popping until we meet a left bracket, or the stack is emptied
        else if (code == RBRACKET) {
            while (!opStack.empty() && opStack.top() != LBRACKET) {
                postfix.push_back({ opStack.top() <= 4 ? BIN_OPERATOR : FUNCTION, 0.0, opStack.top(), {}, emptyVec });
                opStack.pop();
            }
            if (opStack.top() == LBRACKET) {
                opStack.pop();
            }
        }
    }
    // pop out all the values still left on the stack, to complete the postfix expression. It can only either be a binary
    // operator, or a function, which are the only two things not instantly pushed onto the postfix expression, but rather
    // onto the stack first
    while (!opStack.empty()) {
        postfix.push_back({ opStack.top() <= 4 ? BIN_OPERATOR : FUNCTION, 0.0, opStack.top(), {}, emptyVec });
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
        
        //iterator that lags behind by 1, needed for the rotations when expressions are simplified
        std::vector<Token>::iterator it = postfix.begin();
        
        for (size_t i = 1; i < postfix.size(); i++) {
            // we find a bessel function and two immediates, so it can be simplified to one immediate = bessel(imm1, imm2)
            if (postfix[i].type == FUNCTION && postfix[i].op >= BESSELJ && postfix[i - 1].type == IMMEDIATE && postfix[i - 2].type == IMMEDIATE){
                postfix[i-2].num = evalFuncTwoArg(postfix[i].op, postfix[i-2].num.real(), postfix[i-1].num);
                
                //need to shift expression left by 1, to fill empty slot left by evaluating function with two args
                if (i != postfix.size() - 1) {
                    std::rotate(it, it + 2, postfix.end());
                }
                postfix.pop_back();
                postfix.pop_back();
                changed = true;
            }
            // we find a function of one variable and one immediate, so it can be simplified to one immediate = func(imm1)
            if(postfix[i].type == FUNCTION && postfix[i - 1].type == IMMEDIATE){
                postfix[i - 1].num = evalFunc(postfix[i].op, postfix[i - 1].num);
                
                //need to shift expression left by 1, to fill empty slot left by evaluating function
                if (i != postfix.size() - 1){
                    std::rotate(it + 1, it + 2, postfix.end());
                }
                postfix.pop_back();
                changed = true;
            }
            // we find a binary operator and two immediates, so it can be simplified to one immediate = imm1 (operator) imm2
            else if(postfix[i].type == BIN_OPERATOR && postfix[i - 1].type == IMMEDIATE && postfix[i - 2].type == IMMEDIATE){
                switch (postfix[i].op) {
                    case ADD:
                        postfix[i - 2].num = postfix[i - 2].num + postfix[i - 1].num;
                        break;
                    case SUB:
                        postfix[i - 2].num = postfix[i - 2].num - postfix[i - 1].num;
                        break;
                    case MUL:
                        postfix[i - 2].num = postfix[i - 2].num * postfix[i - 1].num;
                        break;  
                    case DIV:
                        postfix[i - 2].num = postfix[i - 2].num / postfix[i - 1].num;
                        break;
                    case POW:
                        postfix[i - 2].num = std::pow(postfix[i - 2].num , postfix[i - 1].num);
                        break;
                    default:
                        postfix[i - 2].num = 0;
                        break;
                }
                
                //need to shift expression left by 2, to fill empty slot by evaluating binary operator
                if (i != postfix.size() - 1){
                    std::rotate(it, it + 2, postfix.end());
                }
                postfix.pop_back();
                postfix.pop_back();
                changed = true;
            }
            it++;
        }
    }
    
    return postfix;
}

// converts infix string to a vector of strings, to make it easier to parse in later steps
std::vector<std::string> getInfixVec(std::string infix){
    std::vector<std::string> infixVec;
    
    // ( and ) are used in strings that represent complex numbers - a + ib is written as (a, b) when its a string for C++
    // to correctly parse it into a complex number. Thus, normal ( and ) would conflict, and thus we replace them with 
    // { and } which cause no issues  
    std::replace(infix.begin(), infix.end(), '(', '{');
    std::replace(infix.begin(), infix.end(), ')', '}');
    
    //Isolate each operator, function, number and variable as an individual token string
    for (size_t i = 0; i < infix.size(); i++) {
        if (infix[i] == '\\') {
            // if it is a sum, we simply take the entire string related to the sum and let the parser deal with it later
            if (infix.substr(i + 1, 3) == "sum"){
                int j = 0;
                size_t k = i;
                for (; k < infix.size(); k++){
                    if (infix[k] == '{'){
                        j++;
                    }
                    else if (infix[k] == '}'){
                        j--;
                        if (j == 0){
                            std::string tempStr =  infix.substr(i + 1, k - i);
                            infixVec.push_back(tempStr);
                            i = k;
                            break;
                        }
                    }
                }
            }
            // if it is a diff, we simply take the entire string related to the derivative and let the parser deal with it later
            else if (infix.substr(i + 1, 4) == "diff"){
                int j = 0;
                size_t k = i;
                for (; k < infix.size(); k++){
                    if (infix[k] == '{'){
                        j++;
                    }
                    else if (infix[k] == '}'){
                        j--;
                        if (j == 0){
                            std::string tempStr =  infix.substr(i + 1, k - i);
                            infixVec.push_back(tempStr);
                            i = k;
                            break;
                        }
                    }
                }
            }
            // it is some other function, for which the end of the function is found using getOp
            else{
                int j = getOp(infix, i);
                std::string tempStr = infix.substr(i + 1, j - i - 1);
                infixVec.push_back(tempStr);
                i = j - 1;
            }
        }
        // we find a negative sign right after a bracket, so we push a 0 in front for unambiguity
        else if (i > 0 && infix[i] == '-' && infix[i - 1] == '{') {
            infixVec.push_back("0");
            infixVec.push_back("-");
        }
        // we find a negative sign at start of expression, so we push a 0 in front for unambiguity
        else if (i == 0 && infix[i] == '-') {
            infixVec.push_back("0");
            infixVec.push_back("-");
        }
        // we find pi, so we push pi
        else if (infix.substr(i, 2) == "pi") {
            infixVec.push_back("pi");
            i++;
        }
        // we find a number, so we look for the end of the number, which is when we would find the first
        // non-digit or non-period character
        else if (std::isdigit(infix[i]) || infix[i] == '.') {
            int j = i;
            std::string tempStr = "";
            while (j < infix.size() && (std::isdigit(infix[j]) || infix[j] == '.')) {
                tempStr += infix[j];
                j++;
            }
            // if there is an i at the end, then it's an imaginary number
            if (j < infix.size() && infix[j] == 'i') {
                infixVec.push_back("(0,"+tempStr+")");
                j++;
            }
            else {
                infixVec.push_back(tempStr);
            }
            i = j - 1;
        }
        // find the other way to input complex numbers, [a,b]
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
        // other cases - either we have a variable, or we have a single character quantity (i, e, +, -, *, /, '(', ')', ',')
        else {
            bool tokenAdded = false;
            //check if token is a variable, and add it
            for (int j = 0; j < variables.size(); j++){
                if (i < infix.size() - variableNames[j].size() - 1 && infix.substr(i, variableNames[j].size()) == variableNames[j]){
                    infixVec.push_back(variableNames[j]);
                    i += variableNames[j].size() - 1;
                    tokenAdded = true;
                    break;
                }
            }
            
            // we skip commas, which seperate arguments in bessel functions, and we do nothing if a token was already added as a variable
            if (infix[i] != ',' && !tokenAdded){
                infixVec.push_back(infix.substr(i, 1));
            }
        }
    }
    
    return infixVec;
}

// for the nth derivative (order = n), get central finite difference coeffs (times 2)
std::vector<int> getFiniteDiffCoeffs(int order){
    // initialize vector of finite difference coefficients of exact size needed to numerically 
    // compute derivative of certain order
    std::vector<int> coeffs((order + 1) / 2 * 2 + 1, 0);
    
    // the stored coefficients are double the actual coefficients this is so the coefficients 
    // are all integers, as some of them have decimal part 0.5. Thus, my storing the double of
    // the coefffs, we can store them using integers instead of doubles or floats
    // can find coeffs for even orders using Pascal's triangle 
    if (order % 2 == 0) {
        coeffs[0] = 2;
        coeffs[order] = 2;
        for (int i = 1; i <= order/2; i++) {
            coeffs[i] = - (order - i + 1) * coeffs[i - 1];
            coeffs[i] /= i;
            coeffs[order - i] = coeffs[i];
        }
    }
    // for odd numbers, need to find even coeffs for previous order (which is always even) then sum 
    // coeffs from previous order to go to next order as is done to find the next order central finite
    // difference typically. instead of storing the even order coeffs, we simply only keep the 
    // previously computed one, as that is the only one needed, and we do the sums to create the next
    // row (the one we need) along the way.
    else {
        int prevEvenCoeff = 1;
        for (int i = 0; i < coeffs.size()/2; i++) {
            coeffs[i] -= prevEvenCoeff;
            coeffs[i + 2] += prevEvenCoeff;
            // in the last iteration, coeffs[i] = coeffs[size - i - 3], and coeffs[i + 2] = 
            // coeffs[size - i - 3], and we would double up the last coefficient in the last
            // iteration, so to prevent that, we do not run these two lines if it is the last loop
            if (i != coeffs.size()/2 - 1){
                coeffs[coeffs.size() - i - 1] += prevEvenCoeff;
                coeffs[coeffs.size() - i - 3] -= prevEvenCoeff;
            }
            std::cout << std::endl;
            prevEvenCoeff = - (coeffs.size() - i - 3) * prevEvenCoeff;
            prevEvenCoeff /= i + 1;
        }
    }
    
    return coeffs;
}

//gets operation/function/number indices that can be used to isolate token 
int getOp(std::string& infix, int n) {
    for (size_t i = n + 1; i < infix.size(); i++) {
        if (infix[i] == '\\' || infix[i] == '-' || infix[i] == '+' || infix[i] == '*' || infix[i] == '/' || infix[i] == '^' || infix[i] == '{' || infix[i] == '}' || infix[i] == ',') {
            if (i != n + 1) {
                return i;
            }
        }
    }
    return infix.size();
}

//Gets operation code
int getOpCode(std::string token) {
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
    else if (token == "digamma")
        return DIGAMMA;
    else if (token == "airy")
        return AIRY;
    else if (token == "biry")
        return BIRY;
    else if (token == "besselj")
        return BESSELJ;
    else if (token == "bessely")
        return BESSELY;
    else if (token == "sphbesselj")
        return SPHBESSELJ;
    else if (token == "sphbessely")
        return SPHBESSELY;
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
    case DIGAMMA: 
        return digamma(z);
    case AIRY:
        return airy(z);
    case BIRY:
        return biry(z);
    default:
        return 0.0;
    }
}

//evaluates functions with two arguments (the bessel functions, for now)
std::complex<double> evalFuncTwoArg(int opCode, std::complex<double> arg1, std::complex<double> arg2) {
    switch(opCode){
        case BESSELJ:
            return cylBesselJ(arg1.real(), arg2);
        case BESSELY:
            return cylBesselY(arg1.real(), arg2);
        case SPHBESSELJ:
            return sphBesselJ(arg1.real(), arg2);
        case SPHBESSELY:
            return sphBesselY(arg1.real(), arg2);
        default:
            return 0.0;
    }
}

//evaluates the overall expression using a stack
//this is often the bottleneck on the program speed
std::complex<double> f(std::complex<double> z) {
    variables[0] = z;
    return f(expr);
}

//evaluate a given expression using the variable values stored in the variables vector
//deals with sums recursively
std::complex<double> f(std::vector<Token> &postfix){
    std::complex<double> temp1, temp2;
    
    for (std::vector<Token>::iterator it = postfix.begin(); it != postfix.end(); it++) {
        // if we encounter immediate, push it onto stack
        if (it->type == IMMEDIATE) {
            evalStack.push(it->num);
        }
        // if we encounter variable, find its value and push value onto stack
        else if (it->type >= VARIABLE) {
            evalStack.push(variables[it->type]);
        }
        // if we encounter function, pop required number of numbers, eval function, and push result
        else if (it->type == FUNCTION) {
            if (it->op >= BESSELJ){
                temp1 = evalStack.top();
                evalStack.pop();
                temp2 = evalStack.top();
                evalStack.pop();
                evalStack.push(evalFuncTwoArg(it->op, temp2, temp1));
            }
            else{
                temp1 = evalStack.top();
                evalStack.pop();
                evalStack.push(evalFunc(it->op, temp1));
            }
        }
        // if we encounter binary operator, pop two numbers, apply binary operator, and push result
        else if (it->type == BIN_OPERATOR) {
            temp1 = evalStack.top();
            evalStack.pop();
            temp2 = evalStack.top();
            evalStack.pop();
            
            switch (it->op) {
                case ADD:
                    evalStack.push(temp2 + temp1);
                    break;
                case SUB:
                    evalStack.push(temp2 - temp1);
                    break;
                case MUL:
                    evalStack.push(temp2 * temp1);
                    break;
                case DIV:
                    evalStack.push(temp2 / temp1);
                    break;
                case POW:
                    evalStack.push(std::pow(temp2, temp1));
                    break;
                default:
                    evalStack.push(0);
                    break;
            }
        }
        // if we encounter a sum, set the sum index variable, and use a for loop and recursive function calls
        // to evaluate the sum
        else if (it->type == SUM){
            temp1 = 0.0;
            for (int i = it->coeffs[0]; i <= it->coeffs[1]; i++){
                variables[it->op] = (std::complex<double>)i;
                temp1 += f(it->postfix);
            }
            evalStack.push(temp1);
        }
        // if we encounter a derivative, use the precalculated finite difference coefficients to evaluate the
        // derivative, by recursively evaluating the function at the required values 
        else if (it->type == DIFF){
            temp1 = 0.0;
            temp2 = variables[0]; // stores the value we will be manipulating
            
            int n = (it->coeffs.size() - 2)/2;
            
            // order size of optimal step is machine constant^(n+2) for nth derivative
            std::complex<double> h = std::pow(EPS, 1.0/(it->op + 2.0));
            if (it->coeffs[it->coeffs.size() - 1] == -1){
                h *= I; //if doing y derivative, delta is in imaginary direction
            }
            
            for (int i = 0; i < it->coeffs.size() - 1; i++){
                if (it->coeffs[i] != 0){
                    variables[0] = temp2 + (double)(i - n) * h;
                    temp1 += f(it->postfix) * (double)it->coeffs[i];
                }
            }
            
            // divide the total difference by step size to the appropriate power
            // divide by 2, as the coeffs we used were actually double the actual coefficients needed 
            temp1 /= 2.0*std::pow(h, (double)it->op);
            variables[0] = temp2; // restores the z value to original
            evalStack.push(temp1); // pushes result of derivative
        }
    }
    
    auto temp = evalStack.top();
    evalStack.pop();
    return temp;
}
