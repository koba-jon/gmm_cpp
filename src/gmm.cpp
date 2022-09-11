#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <random>
#include <cmath>
#include "gmm.hpp"
#include "parameter.hpp"


// ---------------------------------------
// namespace{math} -> function{solve_LU}
// ---------------------------------------
void math::solve_LU(Parameter &mat, Parameter &L, Parameter &U, Parameter &P){

    long long D = mat.shape[0];
    /********************************/
    L.create({D, D}, 0.0);
    U = mat;
    P.create({D, D}, 0.0);
    for (long long i = 0; i < D; i++){
        P.at(i, i) = 1.0;
    }

    double hold_val;
    double cat;
    long long hold_index;
    /********************************/
    for (long long k = 0; k < D - 1; k++){

        hold_val = 0.0;
        hold_index = 0;
        for (long long j = k; j < D; j++){
            if (hold_val < std::abs(U.at(j, k))){
                hold_val = std::abs(U.at(j, k));
                hold_index = j;
            }
        }

        if (hold_index != k){
            for (long long i = 0; i < D; i++){
                std::swap(L.at(hold_index, i), L.at(k, i));
                std::swap(U.at(hold_index, i), U.at(k, i));
                std::swap(P.at(hold_index, i), P.at(k, i));
            }
        }

        for (long long j = 0; j < k; j++){
            L.at(j, k) = 0.0;
        }
        L.at(k, k) = 1.0;

        for (long long j = k + 1; j < D; j++){
            cat = U.at(j, k) / U.at(k, k);
            L.at(j, k) = cat;
            for (long long i = 0; i < D; i++){
                U.at(j, i) -= U.at(k, i) * cat;
            }
        }

    }

    L.at(D - 1, D - 1) = 1.0;

    return;

}


// -----------------------------------------
// namespace{math} -> function{calc_L_inv}
// -----------------------------------------
Parameter math::calc_L_inv(Parameter &L){

    long long D;
    double tmp;
    Parameter col, L_inv;
    /********************************/
    D = L.shape[0];
    L_inv.create({D, D}, 0.0);
    /********************************/
    for (long long k = D - 1; k >= 0; k--){

        col.create({D}, 0.0);
        col.at(k) = 1.0;
        for (long long j = k + 1; j < D; j++){
            
            tmp = 0.0;
            for (long long i = 0; i < j; i++){
                tmp += L.at(j, i) * col.at(i);
            }
            col.at(j) = -tmp;

        }

        for (long long j = 0; j < D; j++){
            L_inv.at(j, k) = col.at(j);
        }

    }

    return L_inv;

}


// -----------------------------------------
// namespace{math} -> function{calc_U_inv}
// -----------------------------------------
Parameter math::calc_U_inv(Parameter &U){

    long long D;
    double tmp;
    Parameter col, U_inv;
    /********************************/
    D = U.shape[0];
    U_inv.create({D, D}, 0.0);
    /********************************/
    for (long long k = D - 1; k >= 0; k--){

        col.create({D}, 0.0);
        col.at(k) = 1.0 / U.at(k, k);
        for (long long j = k - 1; j >= 0; j--){
            
            tmp = 0.0;
            for (long long i = D - 1; i > j; i--){
                tmp += U.at(j, i) * col.at(i);
            }
            col.at(j) = -tmp / U.at(j, j);

        }

        for (long long j = 0; j < D; j++){
            U_inv.at(j, k) = col.at(j);
        }

    }

    return U_inv;

}


// ------------------------------------------
// namespace{math} -> function{calc_LU_inv}
// ------------------------------------------
Parameter math::calc_LU_inv(Parameter &L_inv, Parameter &U_inv, Parameter &P){

    Parameter dot, LU_inv;
    /********************************/
    dot = math::calc_dot(U_inv, L_inv);
    LU_inv = math::calc_dot(dot, P);

    return LU_inv;

}


// ------------------------------------------
// namespace{math} -> function{calc_LU_det}
// ------------------------------------------
double math::calc_LU_det(Parameter &L, Parameter &U, Parameter &P){

    long long D;
    double det;
    /********************************/
    D = L.shape[0];
    det = 1.0;
    /********************************/
    for (long long i = 0; i < D; i++){
        det *= U.at(i, i);
    }

    constexpr double eps = 0.0001;
    size_t count;
    /********************************/
    count = 0;
    /********************************/
    for (long long i = 0; i < D; i++){
        if (P.at(i, i) < eps){
            count++;
        }
    }
    /********************************/
    if ((count % 2 == 0) && (count % 4 != 0)){
        det *= -1.0;
    }

    return det;

}


// ---------------------------------------------
// namespace{math} -> function{calc_vec_trans}
// ---------------------------------------------
Parameter math::calc_vec_trans(Parameter &vec){

    Parameter vec_T;
    /********************************/
    if (vec.shape[0] == 1){
        vec_T.create({vec.shape[1], 1});
    }
    else{
        vec_T.create({1, vec.shape[0]});
    }
    /********************************/
    for (size_t i = 0; i < vec.numel; i++){
        vec_T.data[i] = vec.data[i];
    }

    return vec_T;

}


// ---------------------------------------
// namespace{math} -> function{calc_dot}
// ---------------------------------------
Parameter math::calc_dot(Parameter &mat1, Parameter &mat2){

    long long rows, cols, length;
    Parameter dot;
    /********************************/
    rows = mat1.shape[0];
    cols = mat2.shape[1];
    length = mat1.shape[1];
    /********************************/
    dot.create({rows, cols}, 0.0);
    /********************************/
    for (long long j = 0; j < rows; j++){
        for (long long i = 0; i < cols; i++){
            for (long long k = 0; k < length; k++){
                dot.at(j, i) += mat1.at(j, k) * mat2.at(k, i);
            }
        }
    }

    return dot;

}


// --------------------------------------------
// namespace{math} -> function{calc_norm_pdf}
// --------------------------------------------
double math::calc_norm_pdf(const Parameter &x_, const Parameter &mu_, const Parameter &sigma_){

    long long D;
    Parameter x, mu, sigma;
    /********************************/
    x = x_;
    mu = mu_;
    sigma = sigma_;
    D = x.shape[0];

    Parameter L, U, P;
    Parameter L_inv, U_inv;
    Parameter sigma_inv;
    double sigma_det;
    /********************************/
    math::solve_LU(sigma, L, U, P);
    L_inv = math::calc_L_inv(L);
    U_inv = math::calc_U_inv(U);
    sigma_inv = math::calc_LU_inv(L_inv, U_inv, P);
    sigma_det = math::calc_LU_det(L, U, P);

    double numerator;
    Parameter dev, dev_T, dot;
    /********************************/
    dev = x - mu;
    dev_T = math::calc_vec_trans(dev);
    dot = math::calc_dot(dev_T, sigma_inv);
    numerator = std::exp(-0.5 * math::calc_dot(dot, dev).at(0, 0));

    double denominator;
    /********************************/
    denominator = std::sqrt(std::pow(2.0 * PI, D) * sigma_det);

    double pdf;
    /********************************/
    pdf = numerator / denominator;

    return pdf;

}


// ------------------------------------------
// namespace{math} -> function{show_matrix}
// ------------------------------------------
void math::show_matrix(const Parameter &mat_){

    Parameter mat = mat_;
    long long rows = mat.shape[0];
    long long cols = mat.shape[1];

    for (long long j = 0; j < rows; j++){
        for (long long i = 0; i < cols; i++){
            std::cout << mat.at(j, i) << " ";
        }
        std::cout << std::endl;
    }

    return;

}


// ---------------------------
// class{GMM} -> constructor
// ---------------------------
GMM::GMM(const bool verbose_, const double eps_){
    this->verbose = verbose_;
    this->eps = eps_;
}


// -----------------------------
// class{GMM} -> function{log}
// -----------------------------
void GMM::log(const std::string &str){
    if (this->verbose){
        std::cout << str << std::flush;
    }
    return;
}


// -----------------------------------------
// class{GMM} -> function{init_parameters}
// -----------------------------------------
void GMM::init_parameters(const std::vector<std::vector<double>> &data){

    // Set 'x'
    for (long long n = 0; n < this->N; n++){
        for (long long d = 0; d < this->D; d++){
            this->x.at(n, d, 0) = data[n][d];
        }
    }

    // Calculate 'mu_gauss' and 'sigma_gauss' for initialization
    Parameter x_sum, x2_sum;
    Parameter mu_gauss, sigma_gauss;
    /********************************/
    x_sum.create({this->D, 1}, 0.0);
    x2_sum.create({this->D, 1}, 0.0);
    for (long long n = 0; n < this->N; n++){
        x_sum += this->x(n);
        x2_sum += this->x(n) * this->x(n);
    }
    mu_gauss = x_sum / double(this->N);
    sigma_gauss = x2_sum / double(this->N) - mu_gauss * mu_gauss;

    // Set 'pi', 'mu' and 'sigma'
    std::mt19937 mt(0);
    /********************************/
    this->pi = 1.0 / double(this->K);
    /********************************/
    for (long long k = 0; k < this->K; k++){
        for (long long d = 0; d < this->D; d++){
            std::normal_distribution<double> dist(mu_gauss.at(d, 0), std::sqrt(sigma_gauss.at(d, 0)));
            this->mu.at(k, d, 0) = dist(mt);
        }
    }
    /********************************/
    for (long long k = 0; k < this->K; k++){
        for (long long d2 = 0; d2 < this->D; d2++){
            for (long long d1 = 0; d1 < this->D; d1++){
                if (d1 == d2){
                    this->sigma.at(k, d2, d1) = sigma_gauss.at(d1, 0);
                }
                else{
                    this->sigma.at(k, d2, d1) = 0.0;
                }
            }
        }
    }

    // Set 'ppdf' and 'ppdf_sum'
    for (long long n = 0; n < this->N; n++){
        this->ppdf_sum.at(n) = 0.0;
        for (long long k = 0; k < this->K; k++){
            this->ppdf.at(n, k) = this->pi.at(k) * math::calc_norm_pdf(this->x(n), this->mu(k), this->sigma(k));
            this->ppdf_sum.at(n) += this->ppdf.at(n, k);
        }
    }

    // Set 'L'
    this->L = 0.0;
    for (long long n = 0; n < this->N; n++){
        this->L += std::log(this->ppdf_sum.at(n));
    }

    return;

}


// ----------------------------------
// class{GMM} -> function{apply_EM}
// ----------------------------------
void GMM::apply_EM(){

    Parameter mu_;
    Parameter sigma_, dev, dev_T, dot;
    double L_old;

    do{

        //////////////////////////////
        // (1) E-step (Expectation) //
        //////////////////////////////

        // (1.1) Calculate probability (gamma)
        for (long long n = 0; n < this->N; n++){
            for (long long k = 0; k < this->K; k++){
                this->gamma.at(n, k) = this->ppdf.at(n, k) / this->ppdf_sum.at(n);
            }
        }

        // (1.2) Calculate sum of probability (Nk)
        for (long long k = 0; k < this->K; k++){
            this->Nk.at(k) = 0.0;
            for (long long n = 0; n < this->N; n++){
                this->Nk.at(k) += this->gamma.at(n, k);
            }
        }


        ///////////////////////////////
        // (2) M-step (Maximization) //
        ///////////////////////////////

        // (2.1) Calculate mixing coefficient (pi)
        this->pi = this->Nk / double(this->N);

        // (2.2) Calculate mean (mu)
        for (long long k = 0; k < this->K; k++){
            mu_.create({this->D, 1}, 0.0);
            for (long long n = 0; n < this->N; n++){
                mu_ += this->x(n) * this->gamma.at(n, k);
            }
            mu_ /= this->Nk.at(k);
            this->mu.inplace({k}, mu_);
        }

        // (2.3) Calculate covariance matrix (sigma)
        for (long long k = 0; k < this->K; k++){
            sigma_.create({this->D, this->D}, 0.0);
            for (long long n = 0; n < this->N; n++){
                dev = this->x(n) - this->mu(k);
                dev_T = math::calc_vec_trans(dev);
                dot = math::calc_dot(dev, dev_T);
                sigma_ += dot * this->gamma.at(n, k);
            }
            sigma_ /= this->Nk.at(k);
            this->sigma.inplace({k}, sigma_);
        }

        // (2.4) Calculate joint probability of mixing coefficient and normal distribution-PDF (ppdf, ppdf_sum)
        for (long long n = 0; n < this->N; n++){
            this->ppdf_sum.at(n) = 0.0;
            for (long long k = 0; k < this->K; k++){
                this->ppdf.at(n, k) = this->pi.at(k) * math::calc_norm_pdf(this->x(n), this->mu(k), this->sigma(k));
                this->ppdf_sum.at(n) += this->ppdf.at(n, k);
            }
        }

        // (2.5) Update likelihood (L)
        L_old = this->L;
        this->L = 0.0;
        for (long long n = 0; n < this->N; n++){
            this->L += std::log(this->ppdf_sum.at(n));
        }
        this->log("likelihood = " + std::to_string(this->L) + "\n");


    }while(this->L - L_old > this->eps);


    return;

}


// -------------------------------
// class{GMM} -> function{train}
// -------------------------------
void GMM::train(const std::vector<std::vector<double>> &data, const long long D_, const long long K_){

    // (1) Create parameters
    this->N = data.size();
    this->K = K_;
    this->D = D_;
    this->x.create({this->N, this->D, 1});
    this->pi.create({this->K});
    this->mu.create({this->K, this->D, 1});
    this->sigma.create({this->K, this->D, this->D});
    this->ppdf.create({this->N, this->K});
    this->ppdf_sum.create({this->N});
    this->gamma.create({this->N, this->K});
    this->Nk.create({this->K});

    // (2) Initialize parameters
    this->init_parameters(data);

    // (3) Training based on EM algorithm for GMM
    this->apply_EM();

    return;

}

