#define ARMA_64BIT_WORD
#include <iostream>
#include <fstream>
#include <sstream>
#include <armadillo>
#include <cmath>
#include <random>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include "argh.hpp"

using namespace std;
using namespace arma;

typedef unsigned long long ullong;

double wf_sampling_coefficient(ullong i, ullong Ne, double s = 0, double h = 0.5, double u = 0, double v = 0) {
    double w_11 = 1 + s;
    double w_12 = 1 + (s * h);
    double w_22 = 1;
    double a = w_11 * i * i;
    double b = w_12 * i * (Ne - i);
    double c = w_22 * (Ne - i) * (Ne - i);
    double w_bar = a + (2 * b) + c;
    return (((a + b) * (1 - u)) + ((b + c) * v)) / w_bar;
}

mat wf_transition_matrix(ullong N, double s = 0, double h = 0.5, double u = 0, double v = 0) {
    ullong Ne = 2 * N;
    ullong size = Ne + 1;
    ullong st = size - 2;
    mat Q(st, st);
    
    for (ullong i = 0; i < st; i++) {
        double psi = wf_sampling_coefficient(i + 1, Ne, s, h, u, v);
        for (ullong j = 0; j < st; j++) {
            Q(i, j) = gsl_ran_binomial_pdf(j + 1, psi, Ne);
        }
    }
    return Q;
}

mat cumulative_row_sums(mat P) {
    mat T(P.n_rows, P.n_cols);
    for (ullong i = 0; i < P.n_rows; i++) {
        double s = 0;
        for (ullong j = 0; j < P.n_cols; j++) {
            s += P(i, j);
            T(i, j) = s;
        }
    }
    return T;
}

mat reversed_transient_matrix(mat Q, colvec N1) {
    mat Qp(Q.n_rows, Q.n_cols);
    for (ullong i = 0; i < Q.n_rows; i++) {
        for (ullong j = 0; j < Q.n_cols; j++) {
            Qp(i, j) = Q(j, i) * N1(j) / N1(i);
        }
    }
    return Qp;
}

// Thread-safe rng
// each rng is thread-local
// each thread gets a different seed based on thread number
mt19937_64* get_thread_rng(ullong seed) {
    static thread_local mt19937_64* rng = nullptr;
    if (!rng) {
        ullong thread_seed = seed ^ omp_get_thread_num();
        rng = new mt19937_64(thread_seed);
    }
    return rng;
}

uvec simulate_allele_age_parallel(const mat cQ, const ullong observed, ullong replicates, ullong seed) // {{{
{
    
    //uvec ages(replicates);
    vector<ullong> ages(replicates, 0);
    uniform_real_distribution<double> unif(0, 1);
    
    #pragma omp parallel for
    for(ullong i = 0; i < replicates; i++) {
        ullong state = observed;
        ullong j = 0;
        double u = 0;
        ullong age = 0;

        while (state != 0) {
            // Draw U and linear search
            mt19937_64* rng = get_thread_rng(seed);
            u = unif(*rng);
            try {
                for(j = 0; u > cQ(state, j); j++);
                state = j;
                age++;
            } catch (logic_error& e) {
                cerr << "Out of bound at state: " << state << "; j: " << j << "; u: " << u << endl;
                cerr << cQ.row(state) << endl;
                exit(3);
            }
        }
        ages[i] = age - 1;
    }
    return uvec(ages);
}
// }}}

uvec simulate_allele_freq_trajectory(const mat cQ, const ullong observed) // {{{
{
    
    vector<ullong> freq;
    
    ullong state = observed;
    ullong j = 0;
    double u = 0;

    while (state != 0) {
        // Draw U and linear search
        u = (double) rand() / (double)(RAND_MAX);
        try {
            for(j = 0; u > cQ(state, j); j++);
            state = j;
            freq.push_back(state);
        } catch (logic_error& e) {
            cerr << "Out of bound at state: " << state << "; j: " << j << "; u: " << u << endl;
            cerr << cQ.row(state) << endl;
            cerr << sum(cQ.row(state)) << endl;
            exit(3);
        }
    }
    return uvec(freq);
} // }}}

int main(int argc, char *argv[])
{
    // {{{ Parse args
    ullong Ne, x, r, seed;
    double t, h;
    string out_suffix;

    argh::parser opts(argc, argv);

    if(!(opts("population_size") >> Ne)) {
        cerr << "Must provide population_size" << endl;

        cerr << "USAGE:" << endl <<
            "./allele_age_simulator --population_size=N [--theta=0.0] [--dominance=0.5] [--obsereved=1] [--replicates=1] [--stdout] [--suffix=out]" << endl <<
            "population_size    number of individuals in a population" << endl <<
            "theta              population-scaled mutation rate. Assumes bidirectinally equal mutation rates" << endl <<
            "dominance          Ewens' dominance coefficient h" << endl <<
            "observed           number of allele copies observed" << endl <<
            "replicates         number of simulation replicates" << endl << 
            "seed               random number seed" << endl <<
            "stdout             output results to STDOUT. If not specified, `.csv` with input parameters is used as output" << endl <<
            "suffix             suffix to append to out file name. Ignored if `stdout` is present.`" << endl;

        exit(2);
    }

    opts("theta", 0.0) >> t;
    opts("dominance", 0.5) >> h;
    opts("observed", 1) >> x;
    opts("replicates", 1) >> r;
    opts("suffix", "out") >> out_suffix;

    if (!(opts("seed") >> seed)) {
        seed = time(NULL);
    }
    // }}}

    cerr << "Using random seed " << seed << endl;
    
    double mu = t / (4 * Ne);
    
    mat Q = wf_transition_matrix(Ne, 0, h, mu, mu);
    ullong size = Q.n_rows;
    mat I = eye(size, size);

    cerr << "Solving" << endl;
    vec N = solve(trans(I - Q), I.col(0));

    double R = 1 / N(0);

    cerr << "Building matrix Qp" << endl;
    mat Qp(size + 1, size + 1, fill::zeros);
    Qp(0, 0) = 1;
    Qp(1, 0) = R;
    Qp.submat(1, 1, size, size) = reversed_transient_matrix(Q, N);
    
    cerr << "Building matrix cQ" << endl;
    mat cQ = cumulative_row_sums(Qp);
    
    cerr << "Simulating" << endl;
    //uvec freq = simulate_allele_freq_trajectory(cQ, x);
    uvec ages = simulate_allele_age_parallel(cQ, x, r, seed);

    if (opts["stdout"]) {
        ages.print();
    } else {
        stringstream out_name;
        out_name << "ages_N_" << Ne << "_theta_" << t << "_h_" << h << "_x_" << x << "_" << out_suffix << ".csv";
        ages.save(out_name.str(), raw_ascii);
    }
    
    return 0;
}
