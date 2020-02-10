library(parallelDist)
library(RcppArmadillo)
library(RcppXPtrUtils)

# Forstner distance
forstnerFuncPtr <- cppXPtr(
"double customDist(const arma::mat &A, const arma::mat &B) {
double tol = 1e-3;
arma::vec eigval = real(eig_pair(A, B));
eigval = eigval.elem( find_finite(eigval) );
arma::vec leigval = log(eigval.elem(find(eigval > tol)));
return sqrt(accu(square(leigval)));
}", depends = c("RcppArmadillo"))
forstner.dist <- function(x) parDist(x, method="custom", func = forstnerFuncPtr)


# Riemannian metric
riemannFuncPtr <- cppXPtr(
"double customDist(const arma::mat &A, const arma::mat &B) {
#include <cmath>
double tol = 1e-3;
double k = 0.5;
arma::mat AU, BU, V;
arma::vec Ad, Bd;
//svd_econ(AU, Ad, V, A, \"left\");
svd(AU, Ad, V, A);
//svd_econ(BU, Bd, V, B, \"left\");
svd(BU, Bd, V, B);
Ad = Ad.elem(find(Ad > tol));
Bd = Bd.elem(find(Bd > tol));
int n = std::min(Ad.n_elem, Bd.n_elem);
Ad = Ad.subvec(0,n-1);
Bd = Bd.subvec(0,n-1);
arma::mat AD = diagmat(Ad);
arma::mat BD = diagmat(Bd);
AU = AU.cols(0,n-1);
BU = BU.cols(0,n-1);
arma::mat BD_isqrt = diagmat(1/sqrt(Bd));
arma::mat M = BU - AU * (AU.t() * BU);
arma::vec Md = svd(M);
Md.elem(find(Md > 1)).ones();
arma::vec angles = asin(Md)/1.570796;
double p2 = norm(logmat(BD_isqrt*AD*BD_isqrt), \"fro\");
return k*accu(square(angles)) + (1-k)*p2*p2;
}", depends = c("RcppArmadillo"))
riemannian.dist <- function(x) parDist(x, method="custom", func = riemannFuncPtr)

# Frobenius distance
frobeniusFuncPtr <- cppXPtr(
"double customDist(const arma::mat &A, const arma::mat &B) {
return norm(A-B, \"fro\");
}", depends = c("RcppArmadillo"))
frobenius.dist <- function(x) parDist(x, method="custom", func = frobeniusFuncPtr)

## Spectral distance
spectralFuncPtr <- cppXPtr(
"double customDist(const arma::mat &A, const arma::mat &B) {
arma::mat x = symmatu(A-B);
arma::mat C = symmatu(x.t() * x);
arma::vec eigval = eig_sym(C);
return sqrt(eigval[eigval.n_elem-1]);
}", depends = c("RcppArmadillo"))
spectral.dist <- function(x) parDist(x, method="custom", func = spectralFuncPtr)

cmdFuncPtr <- cppXPtr(
"double customDist(const arma::mat &A, const arma::mat &B) {
return 1 - trace( A * B )/ (norm(A, \"fro\")*norm(B, \"fro\"));
}", depends = c("RcppArmadillo"))
corrmat.dist <- function(x) parDist(x, method="custom", func = cmdFuncPtr)
