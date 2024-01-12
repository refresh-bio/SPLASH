#pragma once
#include <cmath>

#include "../../libs/cdflib/cdflib.hpp"

template <typename T, int MaxN>
class NChooseK {
	T tab[MaxN + 1][MaxN + 1]{ 0 };

public:
	NChooseK() {
		tab[0][0] = 1;
		for (int n = 1; n <= MaxN; ++n) {
			tab[n][0] = 1;
			for (int k = 1; k <= n; ++k) {
				tab[n][k] = tab[n - 1][k - 1] + tab[n - 1][k];
			}
		}
	}

	T operator()(int n, int k) const {
		return tab[n][k];
	}
};


template <int MaxN>
class Binomial {
	NChooseK<uint64_t, MaxN> n_choose_k;
	double tab[MaxN + 1][MaxN + 1]{ 0 };
	double p;

public:
	Binomial(double p_success) : p(p_success) {
		for (int n = 0; n <= MaxN; ++n) {
			for (int k = 0; k <= n; ++k) {
				tab[n][k] = (double)n_choose_k(n, k) * std::pow(p, k) * std::pow(1.0 - p, n - k);
			}
		}
	}

	double pdf(int n, int k) const { return tab[n][k]; }

};


template <int MaxK>
class Poisson {
	double* v_pdf;
	double* v_cdf;
	double lambda;

public:
	Poisson(double lambda) : v_pdf(new double[MaxK + 1]), v_cdf(new double[MaxK + 1]), lambda(lambda) {
		double total_p = 0;
		for (int k = 0; k <= MaxK; k++) {
			double p = std::exp((double)k * std::log(lambda) - std::lgamma((double)k + 1.0) - lambda);
			v_pdf[k] = p;
			
			total_p += p;
			v_cdf[k] = total_p;
		}
	}

	~Poisson() {
		delete[] v_pdf;
		delete[] v_cdf;
	}

	double pdf(int k) const { return v_pdf[k]; }

	double cdf(int k) const { return v_cdf[k]; }
};



inline double poissonCdf(double lambda, double s) {
	
	int which = 1; //  Calculate P and Q from S and XLAM
	double p;
	double q;
	int status;
	double bound;

	cdfpoi(&which, &p, &q, &s, &lambda, &status, &bound);

	return p;
}