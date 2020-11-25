#include "complexPlotter.h"
#include "interpretFunction.h"

int main(int argc, char *argv[])
{
	double realMin = -1, imagMin = -1, realMax = 1, imagMax = 1;
	double step = 2/3840.0, xstep = -1, ystep = -1;
	bool ystepmod = false, xstepmod = false;
	bool equiAngleLines = false, axis = false;
	int colorScheme = 0;
	int lineType = -1; //only -1, 0-10
	std::string func = "z";
	std::string name = "saved";
	bool help = false;
	std::string temp = "";

	if (argc == 1) {
		std::cout << "Enter Function: ";
		std::cin >> func;
		std::cout << "\nDomain Parameters" << std::endl;
		std::cout << "Enter real min: ";
		std::cin >> realMin;
		std::cout << "Enter real max: ";
		std::cin >> realMax;
		std::cout << "Enter imag min: ";
		std::cin >> imagMin;
		std::cout << "Enter imag max: ";
		std::cin >> imagMax;
		std::cout << "\nPlot Parameters" << std::endl;
		std::cout << "x-step: ";
		std::cin >> xstep;
		xstepmod = true;
		std::cout << "y-step: ";
		std::cin >> ystep;
		ystepmod = true;
		std::cout << "Lines of constant real and imaginary components (f/0-10)?: ";
		std::cin >> temp;
		if (temp=="f") {
			lineType = -1;
		}
		else{
		    for (int j = 0; j < 11; j++){
		        if (temp == std::to_string(j)) {
			        lineType = j;
		        }
		    }
		}
		std::cout << "Equal phase lines (t/f)?: ";
		std::cin >> temp;
		equiAngleLines = temp == "t" ? true : false;
		std::wcout << "Grid Parallel to axes lines (t/f)?: ";
		std::cin >> temp;
		axis = temp == "t" ? true : false;
		std::cout << "Use the smooth (s) or regular HSV (h) color scheme: ";
		std::cin >> temp;
		if (temp == "s") {
			colorScheme = 0;
		}
		else if (temp == "h") {
			colorScheme = 1;
		}
		std::cout << "\nFile name (Do NOT include file type like .png): ";
		std::cin >> name;
	}

	//Parsing Command Line Arguments
	for (int i = 1; i < argc; i++) {
		if (!std::strcmp(argv[i], "-realMin") || !std::strcmp(argv[i], "-x1"))
			realMin = std::strtod(argv[i + 1], 0);
		else if (!std::strcmp(argv[i], "-realMax") || !std::strcmp(argv[i], "-x2"))
			realMax = std::strtod(argv[i + 1], 0);
		else if (!std::strcmp(argv[i], "-imagMin") || !std::strcmp(argv[i], "-y1"))
			imagMin = std::strtod(argv[i + 1], 0);
		else if (!std::strcmp(argv[i], "-imagMax") || !std::strcmp(argv[i], "-y2"))
			imagMax = std::strtod(argv[i + 1], 0);
		else if (!std::strcmp(argv[i], "-step") || !std::strcmp(argv[i], "-s"))
			step = std::strtod(argv[i + 1], 0);
		else if (!std::strcmp(argv[i], "-xstep") || !std::strcmp(argv[i], "-xs")) {
			xstepmod = true;
			xstep = std::strtod(argv[i + 1], 0);
		}
		else if (!std::strcmp(argv[i], "-ystep") || !std::strcmp(argv[i], "-ys")) {
			ystepmod = true;
			ystep = std::strtod(argv[i + 1], 0);
		}
		else if (!std::strcmp(argv[i], "-name") || !std::strcmp(argv[i], "-n"))
			name = argv[i + 1];
		else if (!std::strcmp(argv[i], "-axes") || !std::strcmp(argv[i], "-x")) {
			if (!std::strcmp(argv[i + 1], "t")) {
				axis = true;
			}
			if (!std::strcmp(argv[i + 1], "f")) {
				axis = false;
			}
		}
		else if (!std::strcmp(argv[i], "-grid") || !std::strcmp(argv[i], "-g")) {
			if (!std::strcmp(argv[i + 1], "f")) {
				lineType = -1;
			}
			else{
			    temp = argv[i+1];
			    for (int j = 0; j < 11; j++){
		            if (temp == std::to_string(j)) {
			            lineType = j;
		            }
		        }
			}
		}
		else if (!std::strcmp(argv[i], "-angleLines") || !std::strcmp(argv[i], "-l")) {
			if (!std::strcmp(argv[i + 1], "f")) {
				equiAngleLines = false;
			}
			if (!std::strcmp(argv[i + 1], "t")) {
				equiAngleLines = true;
			}
		}
		else if (!std::strcmp(argv[i], "-func") || !std::strcmp(argv[i], "-f")) {
			func = argv[i + 1];
		}
		else if (!std::strcmp(argv[i], "-cScheme") || !std::strcmp(argv[i], "-c")) {
			if (!std::strcmp(argv[i + 1], "s")) {
				colorScheme = 0;
			}
			else if (!std::strcmp(argv[i + 1], "h")) {
				colorScheme = 1;
			}
		}
		else if (!std::strcmp(argv[i], "-help") || !std::strcmp(argv[i], "-h")) {
			help = true;
			std::cout << "Command Line Arguments: " << std::endl;
			std::cout << "\t-help (-h)" << std::endl;
			std::cout << "\t-realMin (-x1) [decimal number]" << std::endl;
			std::cout << "\t-realMax (-x2) [decimal number]" << std::endl;
			std::cout << "\t-imagMin (-y1) [decimal number]" << std::endl;
			std::cout << "\t-imagMax (-y2) [decimal number]" << std::endl;
			std::cout << "\t-step (-s) [decimal number]" << std::endl;
			std::cout << "\t-xstep (-xs) [decimal number]" << std::endl;
			std::cout << "\t-ystep (-ys) [decimal number]" << std::endl;
			std::cout << "\t-grid (-g) [f/0-10]" << std::endl;
			std::cout << "\t-angleLines (-l) [t/f]" << std::endl;
			std::cout << "\t-axes (-x) [t/f]" << std::endl;
			std::cout << "\t-cScheme (-c) [s/h]" << std::endl;
			std::cout << "\t-name (-n) [string name]" << std::endl;
			std::cout << "\t-func (-f) [function expression]" << std::endl;
			std::cout << "\nFor more information about how to use the arguments and other features of the program, visit https://github.com/pateldhairya0328/Complex-Plotter" << std::endl;
		}
	}

	if (help) {
		return 0;
	}

	if (!xstepmod) {
		xstep = step;
	}
	if (!ystepmod) {
		ystep = step;
	}
    
    //convert function to all lower case
    std::transform(func.begin(), func.end(), func.begin(), [](unsigned char const &c){return std::tolower(c);});    
    
	const int realNum = std::round(((realMax - realMin) / xstep));
	const int imagNum = std::round(((imagMax - imagMin) / ystep));

	//Writing the processed arguments to console so user can verify
	std::cout << "Input Function: " << func << std::endl;
	std::cout << "Domain: [" << realMin << ", " << realMax << "]x[" << imagMin << ", " << imagMax << "]"<< std::endl;
	std::cout << "x-step: " << xstep << "; y-step: " << ystep << "; # Pixels real: "<< realNum << "; # Pixels imag: " << imagNum << std::endl;
	std::cout << "Grid: " << lineType << "; Equal Angle lines: " << (equiAngleLines ? "true" : "false") << "; Axes: " << (axis ? "true" : "false") << std::endl;

	//Initializing the Function Interpreter to set up postfix expression of function
	setStep(std::max(xstep, ystep));
	initFunc(func);
	std::cout << std::endl;

	//Making and saving the plot
	cv::Mat plot(imagNum, realNum, CV_8UC3);
	clock_t t = clock();
	plot_func(realNum, imagNum, realMin, imagMin, xstep, ystep, equiAngleLines, lineType, axis, colorScheme, plot);

	//Writing out time it took to complete the plot in milliseconds
	std::cout << "Time taken: " << 0.000001*(float)(clock() - t) << " s" << std::endl;
	std::cout << "Saving...";
	cv::imwrite(name+".png", plot);
	std::cout << "\rPlot saved as " << name << ".png" << std::endl;

	return 0;
}

//Calculates every pixel's RGB value and saves it into the matrix plot
int plot_func(int realNum, int imagNum, double realMin, double imagMin, double xstep, double ystep, bool equiAngleLines, int lineType, bool axis, int colorScheme, cv::Mat& plot) {
	cv::Mat grid(imagNum, realNum, CV_8UC1);
	cv::Mat A;

	double meanVal = 0, max = std::abs(f(std::complex<double>(realMin, imagMin))), min = max;
	int n = 0;

	clock_t time = clock();
	//Find an approximate average value for the function over the region in order to find an appropriate grid size
	for (int y = 0; y <= imagNum; y+=imagNum/10) {
		double im = imagMin + ystep * y;
		for (int x = 0; x <= realNum; x+=realNum/10) {
			double w = std::abs(f(std::complex<double>(realMin + xstep * x, im)));
			if (!std::isnan(w) && !std::isinf(w) && w > 0) {
				meanVal += w;
				if (w > max) {
					max = w;
				}
				if (w < min) {
					min = w;
				}
				n++;
			}
		}
	}

	meanVal = std::pow(meanVal/max, 1/(double)n);
	
	//declaring some variabless we will use
	std::complex<double> w;
	double arg = 0, mag = 1, im = 0, argn = 0, temp = 0, temp1 = 0, sat = 0;

	//higher value of linetype means larger value of modval, which means the end grid is more separated
	double modVal = meanVal * std::pow(2, lineType - 5), modVal2 = modVal / 2.0; 

	for (int y = 0; y < imagNum; y++) {
		cv::Vec3b* ptr = plot.ptr<cv::Vec3b>(imagNum - y - 1);
		im = imagMin + ystep * y;

		for (int x = 0; x < realNum; x++) {
			//evaluate function
			w = f(std::complex<double>(realMin + xstep * x, im));
			if (!std::isnan(std::abs(w))) {
				arg = std::arg(w);
				mag = std::abs(w);
			}

			//positive mod of argument
			argn = arg + 2 * PI;
			if (arg < 0)
				arg += 2 * PI;

			//decreases the saturation near regions where the argument is close to a multiple of pi/6
			temp = std::abs(fmod(argn - PI / 12, PI / 6) - PI / 12);
			sat = 0.9;
			if (temp < 0.02 && equiAngleLines)
				sat = temp * 45;

			//makes black and white checkerboard, which will be used to make the grid later on
			if (lineType != -1) {
				temp = fmod(w.real() - modVal2, modVal);
				if (temp < 0)
					temp += modVal;
				temp -= modVal2;
				temp1 = fmod(w.imag() - modVal2, modVal);
				if (temp1 < 0)
					temp1 += modVal;
				temp1 -= modVal2;

				if (temp * temp1 > 0) {
					grid.at<uchar>(imagNum - y - 1, x) = 255;
				}
				else {
					grid.at<uchar>(imagNum - y - 1, x) = 0;
				}
			}

			//work with mod of log of the magnitude, as it allows for a larger range of magnitudes
			//to be easily distinguished in the final plot
			mag = std::log(mag) / std::log(1.5);
			mag = fmod(mag, 1.5);
			if (mag < 0)
				mag += 1.5;

			mag = (mag + 1.5 * (1 - sat)) / (2 - sat);
			mag += 1.5;

			ptr[x] = getColor(colorScheme, sat, mag, arg);
		}
		if (y % (imagNum/100) == 0) {
			std::cout << "\rProgress: " << 100 * y / imagNum << "%\t\u2502";
			for (int k = 0; k < 20; k++) {
				std::cout << (k < 20 * y / imagNum ? "\u2588": " ");
			}
			std::cout << "\u2502     ; Time taken so far: " << 0.000001*(float)(clock() - time) << " s";
		}
	}

	if (lineType != -1) {
		//Edge detects checkerboard pattern, and then widens the edges using a blur
		int kernelSize = std::min(realNum, imagNum) / 500;
		if (kernelSize % 2 == 0) {
			kernelSize += 1;
		}

		cv::Mat A;
		//cv::GaussianBlur(grid, grid, cv::Size(kernelSize * 2 + 1, kernelSize * 2 + 1), 0, 0, cv::BORDER_DEFAULT);
		cv::Canny(grid, A, 0, 50, 3);
		cv::GaussianBlur(A, A, cv::Size(kernelSize, kernelSize), 0, 0, cv::BORDER_DEFAULT);

		//Makes regions where an edge is appear darker
		for (int y = 0; y < A.rows; y++) {
			cv::Vec3b* ptr = plot.ptr<cv::Vec3b>(y);
			const unsigned char* row = A.ptr<unsigned char>(y);
			for (int x = 0; x < A.cols; x++) {
				if (row[x] == 0) {
					ptr[x] /= 0.85;
				}
			}
		}
		
		cv::GaussianBlur(plot, plot, cv::Size(3, 3), 0, 0, cv::BORDER_DEFAULT);
	}

	if (axis) {
		int width = std::round(std::max(realNum, imagNum)/1000.0);
		
		double axisStep = xstep * realNum / 10.0;
		double scale = std::pow(10.0, std::ceil(std::log10(std::abs(axisStep))) + 1);
		axisStep = std::round(axisStep * 2 * scale) / (2 * scale * xstep);
		double x = std::round(realMin * 2 * scale) / (2 * scale);
		if (x == realMin) {
			x += axisStep * xstep;
		}

		for (x = (x - realMin) / xstep; x < realNum; x += axisStep) {
			cv::line(plot, cv::Point(x, 0), cv::Point(x, imagNum - 1), cv::Scalar(50, 50, 50), width);
		}
		
		axisStep = ystep * imagNum / 10.0;
		scale = std::pow(10.0, std::ceil(std::log10(std::abs(axisStep))) + 1);
		axisStep = std::round(axisStep * 2 * scale) / (2 * scale * ystep);
		x = std::round(imagMin * 2 * scale) / (2 * scale);
		if (x == imagMin) {
			x += axisStep * ystep;
		}

		for (x = (x - imagMin) / ystep; x < imagNum; x += axisStep) {
			cv::line(plot, cv::Point(0, x), cv::Point(realNum - 1, x), cv::Scalar(50, 50, 50), width);
		}
	}

	std::cout << "\rProgress: 100%\t\u2502\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2502" << std::endl;

	return 0;
}

cv::Vec3b getColor(int colorScheme, double sat, double mag, double arg) {
	if (colorScheme == 0) {
		//Saturation, Value, Hue
		float S = sat;
		float V = 1 - std::pow(0.5, mag);//0.5 is just a magic number, can choose numbers in (0,1) for varying gradients of values
		float H = 128 * arg / PI;
		int HLo = std::floor(H), HHi = std::ceil(H);

		//Explanation of how the color is calculated
		//(HHi - H) * colormap[HLo] + (H - HLo) * colormap[HHi] linearly interpolates the color associated with the hue 
		//From the look up table of colors
		//S * color + (1 - S) * 255 finds the weighted average of the color and white, with the color being whiter as 
		//the value of S goes down (so higher S results in a more saturated color, S = 1 is most saturated)
		//V * color at the end averages the color with black (the (1 - V) * 0 is always 0 so it is left out), so the color 
		//becomes blacker as the value decreases

		//The color is converted to HSV in this method instead of a more typical HSV to RGB conversion since the normal
		//method results in a very perceptually uneven output (there are visible streaks at yellow, cyan and magenta). 
		//This method instead uses a lookup table of values which are not exactly the correct RGB values, but are a 
		//smoothed out version, which results in a more perceptually smoother plot that will convey less incorrect 
		//information about the plotted data. The HSLuv was not used to create a perceptually uniform colormap, since it
		//looks kinda ugly
		if (HLo != HHi) {
			cv::Vec3b color;
			color[0] = V * (S * ((HHi - H) * colormap[HLo][2] + (H - HLo) * colormap[HHi][2]) + (1 - S) * 255);
			color[1] = V * (S * ((HHi - H) * colormap[HLo][1] + (H - HLo) * colormap[HHi][1]) + (1 - S) * 255);
			color[2] = V * (S * ((HHi - H) * colormap[HLo][0] + (H - HLo) * colormap[HHi][0]) + (1 - S) * 255);
			return color;
		}
		else {
			cv::Vec3b color;
			color[0] = V * (S * colormap[HLo][2] + (1 - S) * 255);
			color[1] = V * (S * colormap[HLo][1] + (1 - S) * 255);
			color[2] = V * (S * colormap[HLo][0] + (1 - S) * 255);
			return color;
		}
	}
	else if (colorScheme == 1) {
		//Saturation, Value, Hue
		double S = sat;
		double V = 1 - std::pow(0.5, mag);//0.5 is just a magic number, can choose numbers in (0,1) for varying gradients of values
		double H = 3 * arg / PI;
		cv::Vec3b color;
		color[0] = hsvToRGB(1, H, S, V);
		color[1] = hsvToRGB(3, H, S, V);
		color[2] = hsvToRGB(5, H, S, V);
		return color;
	}
	return 0;
}

//Converts HSV to one of R, G or B based on argument n
//Conversion done as per the process on https://en.wikipedia.org/wiki/HSL_and_HSV#HSV_to_RGB
int hsvToRGB(int n, double H, double S, double V) {
	double k = fmod(n + H, 6);
	return (int)(255.0 * (V - V * S * std::max(std::min(std::min(k, 4 - k), 1.0), 0.0)));
}
