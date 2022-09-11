#ifndef GMM_HPP
#define GMM_HPP

#include <string>
#include <vector>
#include "parameter.hpp"
#define PI 3.14159265358979


// -----------------
// namespace{math}
// -----------------
namespace math{
    void solve_LU(Parameter &mat, Parameter &L, Parameter &U, Parameter &P);
    Parameter calc_L_inv(Parameter &L);
    Parameter calc_U_inv(Parameter &U);
    Parameter calc_LU_inv(Parameter &L_inv, Parameter &U_inv, Parameter &P);
    double calc_LU_det(Parameter &L, Parameter &U, Parameter &P);
    Parameter calc_vec_trans(Parameter &vec);
    Parameter calc_dot(Parameter &mat1, Parameter &mat2);
    double calc_norm_pdf(const Parameter &x_, const Parameter &mu_, const Parameter &sigma_);
    void show_matrix(const Parameter &mat_);
}


// ------------
// class{GMM}
// ------------
class GMM{

private:

    // Member variable
    bool verbose;
    double eps;
    long long N, K, D;
    Parameter x;  // x: {N, D, 1}
    Parameter pi, mu, sigma;  // pi: {K}, mu: {K, D, 1}, sigma: {K, D, D}
    Parameter ppdf, ppdf_sum;  // ppdf: {N, K}, ppdf_sum: {N}
    Parameter gamma, Nk;  // gamma: {N, K}, Nk: {K}
    double L;
    
    // Function
    void log(const std::string &str);
    void init_parameters(const std::vector<std::vector<double>> &data);
    void apply_EM();

public:
    
    // Constructor
    GMM() = delete;
    GMM(const bool verbose_ = true, const double eps_ = 0.000001);

    // Function
    void train(const std::vector<std::vector<double>> &data, const long long D_, const long long K_ = 4);

};


#endif