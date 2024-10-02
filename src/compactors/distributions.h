#pragma once
#include <cmath>
#include <cinttypes>

#include "../../libs/cdflib/cdflib.hpp"

/*
template <typename T, int MaxN>
class NChooseK {
	const int SIZE = (MaxN + 1);
	T* tab;

public:
	NChooseK() {
		tab = new T[SIZE * SIZE];
		std::fill(tab, tab + SIZE * SIZE, 0);

		tab[0] = 1;
		for (int n = 1; n <= MaxN; ++n) {
			tab[n * SIZE + 0] = 1;
			for (int k = 1; k <= n; ++k) {
				tab[n * SIZE + k] = tab[(n - 1) * SIZE + k - 1] + tab[(n - 1) * SIZE + k];
			}
		}
	}

	T operator()(int n, int k) const {
		return tab[n * SIZE + k];
	}

	~NChooseK() {
		delete[] tab;
	}
};
*/

template <int MaxN>
class NChooseK_log {
	double ret[32 * 32];
	double* prefix_logs;

	double calc(int n, int k) const
	{
		double nominator = prefix_logs[n] - prefix_logs[n - k];
		double denominator = prefix_logs[k];
		return std::exp(nominator - denominator);
	}

public:
	NChooseK_log() : prefix_logs(new double [MaxN + 1]) {

		prefix_logs[0] = prefix_logs[1] = 0;
		double current = 0;
		for (int n = 2; n <= MaxN; ++n) {
			current += std::log((double)n);
			prefix_logs[n] = current;
		}

		for (int n = 0; n < 32; ++n)
			for (int k = 0; k <= n; ++k)
				ret[n * 32 + k] = calc(n, k);
	}

	~NChooseK_log() { delete[] prefix_logs; }

	double operator()(int n, int k) const {
		if (n < 32)
			return ret[n * 32 + k];

		return calc(n, k);
	}
};


template <int MaxN>
class Binomial {
	NChooseK_log<MaxN> n_choose_k;
	double* pows;
	double* inv_pows;
	double p;

public:
	Binomial(double p_success) : p(p_success), pows(new double[MaxN + 1]), inv_pows(new double[MaxN + 1]) {
		for (int i = 0; i <= MaxN; ++i) {
			pows[i] = std::pow(p, i);
			inv_pows[i] = std::pow(1.0 - p, i);
		}
	}

	~Binomial() {
		delete [] pows;
		delete [] inv_pows;
	}

	double pdf(int n, int k) const { return n_choose_k(n, k) * pows[k] * inv_pows[n - k]; }

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
