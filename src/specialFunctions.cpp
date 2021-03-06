#include "specialFunctions.h"

//uses Lanczos approximation to calculate gamma function
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

std::complex<double> digamma(std::complex<double> z){
    std::complex<double> temp1, temp2;
    temp1 = gamma(z + EPS);
    temp2 = gamma(z - EPS);
    return (temp1 - temp2)/(EPS * (temp1 + temp2));
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

// following special functions implemented by calls to external functions, through a wrapper for a 
// FORTRAN library. Need to suppress console output which may show up during errors when calling 
// external functions, which are not important in the grand scheme of things, but can be disruptive
// in the  console output
std::complex<double> airy(std::complex<double> z){
    std::cout.setstate(std::ios_base::failbit); //suppress console output from external function calls
    return sp_bessel::airy(z);
    std::cout.clear(); //unsuppress console output
}

std::complex<double> biry(std::complex<double> z){
    std::cout.setstate(std::ios_base::failbit); //suppress console output from external function calls
    return sp_bessel::biry(z);
    std::cout.clear(); //unsuppress console output
}

std::complex<double> cylBesselJ(double nu, std::complex<double> z){
    std::cout.setstate(std::ios_base::failbit); //suppress console output from external function calls
    return sp_bessel::besselJ(nu, z);
    std::cout.clear(); //unsuppress console output
}

std::complex<double> cylBesselY(double nu, std::complex<double> z){
    std::cout.setstate(std::ios_base::failbit); //suppress console output from external function calls
    return sp_bessel::besselI(nu, z);
    std::cout.clear(); //unsuppress console output
}

std::complex<double> sphBesselJ(double nu, std::complex<double> z){
    std::cout.setstate(std::ios_base::failbit); //suppress console output from external function calls
    return sp_bessel::sph_besselJ(nu, z);
    std::cout.clear(); //unsuppress console output
}

std::complex<double> sphBesselY(double nu, std::complex<double> z){
    std::cout.setstate(std::ios_base::failbit); //suppress console output from external function calls
    return sp_bessel::sph_besselY(nu, z);
    std::cout.clear(); //unsuppress console output
}
