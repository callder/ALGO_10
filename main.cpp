#include <iostream>
#include <cmath>
#include <complex>

using namespace std;


const double pi = acos(-1);

void DFT(const double* input, complex<double>* result, int N) {
    for (int k = 0; k < N; ++k) {
        result[k] = 0.0;
        for (int n = 0; n < N; ++n) {
            double realPart = cos(-2.0 * pi * k * n / N);
            double imagPart = sin(-2.0 * pi * k * n / N);
            complex<double> w(realPart, imagPart);
            result[k] += input[n] * w;
        }
    }
}
void FFT(const double* input, complex<double>* result, int N) {
    if (N <= 1) {
        result[0] = input[0];
    } else {
        double* parzyste = new double[N/2];
        double* nieparzyste = new double[N / 2];
        complex<double>* resultParzyste = new complex<double>[N / 2];
        complex<double>* resultNieparzyste = new complex<double>[N/2];

        for (int i = 0; i < N/2; ++i) {
            parzyste[i] = input[2*i];
            nieparzyste[i] = input[2 * i + 1];
        }

        FFT(parzyste, resultParzyste, N / 2);
        FFT(nieparzyste, resultNieparzyste, N / 2);

        for (int k = 0; k < N/2; ++k) {
            double rzeczywista = cos(-2.0 * M_PI * k / N);
            double urojona = sin(-2.0 * M_PI * k / N);
            complex<double> liczbaZesp(rzeczywista, urojona);

            result[k] = resultParzyste[k] + liczbaZesp * resultNieparzyste[k];
            result[k + N/2] = resultParzyste[k] - liczbaZesp * resultNieparzyste[k];
        }

        delete[] parzyste;
        delete[] nieparzyste;
        delete[] resultParzyste;
        delete[] resultNieparzyste;
    }
}
int main() {
    const int MAX_ORDER = 13;
    const bool PRINT_COEFS = false;

    for (int o = 1; o <= MAX_ORDER; o++) {
        const int N = 1 << o;
        printf("N: %i \n", N);

        double* f = new double[N];
        for (int n = 0; n < N; n++)
            f[n] = n / (double)N;

        clock_t t1 = clock();
        complex<double>* cDFT = new complex<double>[N];
        DFT(f, cDFT, N);
        clock_t t2 = clock();
        double dft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
        printf("DFT time [ms]: %f \n", dft_time);

        t1 = clock();
        complex<double>* cFFT = new complex<double>[N];
        FFT(f, cFFT, N);
        t2 = clock();
        double fft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
        printf("FFT time [ms]: %f \n", fft_time);

        if (PRINT_COEFS) {
            for (int k = 0; k < N; k++)
                cout << "c_" << k << " DFT = " << cDFT[k] << " FFT = " << cFFT[k] << endl;
            cout << "------" << endl;
        }

        delete[] f;
        delete[] cDFT;
        delete[] cFFT;
    }

    return 0;
};