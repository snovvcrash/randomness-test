/**
 * randomness_test.hxx
 *
 * A Series of Randmoness Tests for Binary Sequence Validation
 * by snovvcrash
 * 04.2017
 */

/* A statistical test suite for random and pseudorandom number generators for cryptographic applications
	http://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
GSL - GNU Scientific Library - GNU Project - Free Software Foundation
	https://www.gnu.org/software/gsl */

#pragma once
#ifndef RANDOMNESS_TEST_HXX
#define RANDOMNESS_TEST_HXX

#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <gsl/gsl_sf_erf.h>   // gsl_sf_erfc
#include <gsl/gsl_sf_gamma.h> // gsl_sf_gamma_inc_Q
#include <gsl/gsl_fft_real.h> // gsl_fft_real_wavetable_alloc,
                              // gsl_fft_real_workspace_alloc,
                              // gsl_fft_real_transform,
                              // gsl_fft_real_wavetable_free,
                              // gsl_fft_real_workspace_free
#include <gsl/gsl_cdf.h>      // gsl_cdf_ugaussian_P

#define ALPHA 0.01

namespace statistics {

	class tests {
		std::vector<bool> m_seq;
		std::vector<bool> m_verdict;
		int m_numOnes;
		int m_numZeros;

	public:
		tests(std::vector<bool> seq_) : m_seq(seq_), m_numOnes(0), m_numZeros(0) {
			for (const auto& bit : m_seq)
				if (bit)
					++m_numOnes;
				else
					++m_numZeros;
		}

		void print_verdict() {
			if (m_verdict.empty()) return;

			// Verdict v1

			/*int count = 0;

			for (const auto& ans : m_verdict)
				count += 2*ans - 1;

			if (count >= 0)
				std::cout << "good" << std::endl;
			else
				std::cout << "bad" << std::endl;*/

			// Verdict v2

			for (const auto& ans : m_verdict)
				if (!ans) {
					std::cout << "bad" << std::endl;
					return;
				}

			std::cout << "good" << std::endl;
		}

		// ------------------------------------------------------
		// ----------------- MONOBIT FREQUENCY ------------------
		// ------------------------------------------------------

		void monobit_frequency() {
			int sum = m_numOnes - m_numZeros;

			double s_obs = std::abs(sum) / std::sqrt(m_seq.size());
			double P_value = gsl_sf_erfc(s_obs / std::sqrt(2));

			#ifdef DEBUG
			std::cout << "Monobit Frequency:   " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ------------------ BLOCK FREQUENCY -------------------
		// ------------------------------------------------------

		void block_frequency(int M) {
			int N = m_seq.size() / M;

			if (!N) {
				m_verdict.push_back(0);
				return;
			}

			int    onesCount = 0; // number of ones in block
			double sum = 0.0;
			double pi_i = 0.0;

			for (int i = 0; i < N; ++i) {
				onesCount = 0;
				pi_i = 0.0;
				for (int j = 0; j < M; ++j)
					if (m_seq[j + i*M])
						++onesCount;
				pi_i = (double)onesCount / (double)M;
				sum += std::pow(pi_i - 0.5, 2);
			}

			double chi_square = 4 * M * sum;
			double P_value = gsl_sf_gamma_inc_Q(N / 2.0, chi_square / 2.0);

			#ifdef DEBUG
			std::cout << "Block Frequency:     " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ------------------------ RUNS ------------------------
		// ------------------------------------------------------

		void runs() {
			double pi = (double)m_numOnes / (double)m_seq.size();

			if (std::abs(pi - 0.5) >= 2/std::sqrt(m_seq.size())) {
				m_verdict.push_back(0);
				return;
			}

			int V_n = 0;
			for (size_t i = 0; i < m_seq.size()-1; ++i)
				if (m_seq[i] != m_seq[i+1])
					++V_n;
			++V_n;

			double numerator = std::abs(V_n - 2 * m_seq.size() * pi * (1 - pi));
			double denominator = 2 * std::sqrt(2*m_seq.size()) * pi * (1 - pi);
			if (!denominator) {
				m_verdict.push_back(0);
				return;
			}
			double P_value = gsl_sf_erfc(numerator / denominator);

			#ifdef DEBUG
			std::cout << "Runs:                " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ---------------- LONGEST RUN OF ONES -----------------
		// ------------------------------------------------------

		void longest_run_of_ones(int M) {
			int N = m_seq.size() / M;

			if (!N) {
				m_verdict.push_back(0);
				return;
			}

			std::vector<double> pi(7, 0.0);
			int K = 0;

			switch (M) {
				case 8 :
					pi[0] = 0.2148;
					pi[1] = 0.3672;
					pi[2] = 0.2305;
					pi[3] = 0.1875;
					K = 3;
					break;
				case 128 :
					pi[0] = 0.1174;
					pi[1] = 0.2430;
					pi[2] = 0.2493;
					pi[3] = 0.1752;
					pi[4] = 0.1027;
					pi[5] = 0.1124;
					K = 5;
					break;
				case 512 :
					pi[0] = 0.1170;
					pi[1] = 0.2460;
					pi[2] = 0.2523;
					pi[3] = 0.1755;
					pi[4] = 0.1027;
					pi[5] = 0.1124;
					K = 5;
					break;
				case 10000 :
					pi[0] = 0.0882;
					pi[1] = 0.2092;
					pi[2] = 0.2483;
					pi[3] = 0.1933;
					pi[4] = 0.1208;
					pi[5] = 0.0675;
					pi[6] = 0.0727;
					K = 6;
					break;
				default :
					// std::cout << "statistical::longest_run_of_ones: Invalid block size" << std::endl;
					m_verdict.push_back(0);
					return;
			}

			std::vector<int> v(6, 0);

			for (int i = 0; i < N; ++i) {
				int runCurr = 0;
				int runMax = 0;

				for (int j = 0; j < M; ++j)
					if (m_seq[j + i*M]) {
						++runCurr;
						runMax = std::max(runCurr, runMax);
					}
					else runCurr = 0;

				switch (M) {
					case 8 :
						if (runMax <= 1)
							++v[0];
						else if (runMax == 2)
							++v[1];
						else if (runMax == 3)
							++v[2];
						else if (runMax >= 4)
							++v[3];
						break;
					case 128 :
						if (runMax <= 4)
							++v[0];
						else if (runMax == 5)
							++v[1];
						else if (runMax == 6)
							++v[2];
						else if (runMax == 7)
							++v[3];
						else if (runMax == 8)
							++v[4];
						else if (runMax >= 9)
							++v[5];
						break;
					case 512 :
						if (runMax <= 6)
							++v[0];
						else if (runMax == 7)
							++v[1];
						else if (runMax == 8)
							++v[2];
						else if (runMax == 9)
							++v[3];
						else if (runMax == 10)
							++v[4];
						else if (runMax >= 11)
							++v[5];
						break;
					case 10000 :
						if (runMax <= 10)
							++v[0];
						else if (runMax == 11)
							++v[1];
						else if (runMax == 12)
							++v[2];
						else if (runMax == 13)
							++v[3];
						else if (runMax == 14)
							++v[4];
						else if (runMax == 15)
							++v[5];
						else if (runMax >= 16)
							++v[6];
						break;
				}
			}

			double chi_square = 0.0;
			for (int i = 0; i <= K; ++i)
				chi_square += std::pow((double)v[i] - (double)N*pi[i], 2) / ((double)N*pi[i]);

			double P_value = gsl_sf_gamma_inc_Q(K / 2.0, chi_square / 2.0);

			#ifdef DEBUG
			std::cout << "Longest Run of Ones: " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ---------------------- SPECTRAL ----------------------
		// ------------------------------------------------------

		void spectral() {
			double* data = new double[m_seq.size()];
			for (size_t i = 0; i < m_seq.size(); ++i)
				data[i] = 2*m_seq[i] - 1;

			gsl_fft_real_wavetable *wavetable;
			gsl_fft_real_workspace *work;

			wavetable = gsl_fft_real_wavetable_alloc((double)m_seq.size()); // prepare trigonometric lookup tables for an FFT
			work      = gsl_fft_real_workspace_alloc((double)m_seq.size()); // prepare a workspace for a real transform

			gsl_fft_real_transform(data, 1, m_seq.size(), wavetable, work); // 1 = step (stride)

			gsl_fft_real_wavetable_free(wavetable);
			gsl_fft_real_workspace_free(work);

			int halfLen = m_seq.size()/2 + 1;
			double* M = new double[halfLen];
			M[0] = std::abs(data[0]);

			if (m_seq.size() % 2 == 0) {
				for (int i = 0; i < halfLen-2; ++i)
					M[i+1] = std::sqrt(std::pow(data[2*i+1], 2) + std::pow(data[2*i+2], 2));
				M[m_seq.size()/2] = data[m_seq.size()-1];
			}
			else
				for (int i = 0; i < halfLen-1; ++i)
					M[i+1] = std::sqrt(std::pow(data[2*i+1], 2) + std::pow(data[2*i+2], 2));

			delete [] data;

			double T = std::sqrt(std::log(1/0.05) * m_seq.size());
			double N_0 = 0.95 * m_seq.size() / 2;

			double N_1 = 0.0;
			for (size_t i = 1; i < m_seq.size()/2; ++i)
				if (M[i] < T)
					++N_1;

			delete [] M;

			double d = (N_1 - N_0) / std::sqrt(m_seq.size() * 0.95 * 0.05 / 4);
			double P_value = gsl_sf_erfc(std::abs(d) / std::sqrt(2));

			#ifdef DEBUG
			std::cout << "Spectral (FFT):      " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ----------------------- SERIAL -----------------------
		// ------------------------------------------------------

		double compute_psi(int m) {
			if ((m == 0) || (m == -1)) return 0.0;

			int s = (1 << (m+1)) - 1; // 2^(m+1) - 1
			std::vector<int> P(s, 0);

			for (size_t i = 0; i < m_seq.size(); ++i) {
				int k = 1;

				for (int j = 0; j < m; ++j) {
					if (m_seq[(i+j) % m_seq.size()] == 0)
						k = 2*k;
					else // if (m_seq[(i+j) % m_seq.size()] == 1)
						k = 2*k + 1;
				}

				++P[k-1];
			}

			double sum = 0.0;
			for (auto it = P.begin(); it != P.end(); ++it)
				if (*it) sum += std::pow(*it, 2);

			double psi_square = ((double)(1 << m)/(double)m_seq.size())*sum - m_seq.size();
			return psi_square;
		}

		void serial(int m) {
			double psi2_m  = compute_psi(m);
			double psi2_m1 = compute_psi(m-1);
			double psi2_m2 = compute_psi(m-2);

			double delta_psi2_m  = psi2_m - psi2_m1;
			double delta2_psi2_m = psi2_m - 2*psi2_m1 + psi2_m2;

			double arg1 = (double)(1 << (m-2));
			double P_value1 = gsl_sf_gamma_inc_Q(arg1, delta_psi2_m / 2.0);
			double P_value2 = gsl_sf_gamma_inc_Q(arg1 / 2.0, delta2_psi2_m / 2.0);

			#ifdef DEBUG
			std::cout << "Serial:              " << P_value1 << " " << P_value2 << std::endl;
			#endif // DEBUG

			if (P_value1 >= ALPHA && P_value2 >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ---------------- APPROXIMATE ENTROPY -----------------
		// ------------------------------------------------------

		double compute_phi(int m) {
			int s = (1 << (m+1)) - 1; // 2^(m+1) - 1
			std::vector<int> P(s, 0);
			std::vector<double> C(s, 0.0);

			for (size_t i = 0; i < m_seq.size(); ++i) {
				int k = 1;

				for (int j = 0; j < m; ++j) {
					if (m_seq[(i+j) % m_seq.size()] == 0)
						k = 2*k;
					else // if (m_seq[(i+j) % m_seq.size()] == 1)
						k = 2*k + 1;
				}

				++P[k-1];
			}

			for (size_t i = 0; i < P.size(); ++i)
				if (P[i]) C[i] = (double)P[i] / (double)m_seq.size();

			double phi = 0.0;
			for (auto it = C.begin(); it != C.end(); ++it)
				if (*it) phi += (*it) * std::log(*it);

			return phi;
		}

		void approximate_entropy(int m) {
			double phi_m  = compute_phi(m);
			double phi_m1  = compute_phi(m+1);

			double chi_square = 2 * m_seq.size() * (std::log(2) - (phi_m-phi_m1));
			double P_value = gsl_sf_gamma_inc_Q((double)(1 << (m-1)), chi_square / 2.0);
			
			#ifdef DEBUG
			std::cout << "Approximate Entropy: " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}

		// ------------------------------------------------------
		// ------------------ CUMULATIVE SUMS -------------------
		// ------------------------------------------------------

		void cumulative_sums() {
			int cusumCurr = 0;
			int cusumMax  = 0;

			for (const auto& bit : m_seq) {
				cusumCurr += 2*bit - 1;
				cusumMax = std::max(std::abs(cusumCurr), cusumMax);
			}

			int z = cusumMax;

			if (!z) {
				m_verdict.push_back(0);
				return;
			}

			int size = m_seq.size();

			int start = (-size/z + 1) / 4;
			int end   = (size/z - 1) / 4;
			
			double sum_1 = 0.0;
			for (int k = start; k <= end; ++k)
				sum_1 += gsl_cdf_ugaussian_P((4*k+1)*z/std::sqrt(size)) - gsl_cdf_ugaussian_P((4*k-1)*z/std::sqrt(size));

			start = (-size/z - 3) / 4;
			end   = (size/z - 1) / 4;

			double sum_2 = 0.0;
			for (int k = start; k <= end; ++k)
				sum_2 += gsl_cdf_ugaussian_P((4*k+3)*z/std::sqrt(size)) - gsl_cdf_ugaussian_P((4*k+1)*z/std::sqrt(size));

			double P_value = 1 - sum_1 + sum_2;

			#ifdef DEBUG
			std::cout << "Cumulative Sums:     " << P_value << std::endl;
			#endif // DEBUG

			if (P_value >= ALPHA)
				m_verdict.push_back(1);
			else
				m_verdict.push_back(0);
		}
	};

}

#endif // RANDOMNESS_TEST_HXX
