// We can now use the BH package
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/factorials.hpp>

using namespace Rcpp;

double w_k3(double beta1, double beta2, double k, double lambda){
  if(k <= lambda){
    return exp(-1 * beta1 * (lambda - k));
  } else {
    return exp(-1 * beta2 * (k - lambda));
  }
}//EWP_3 weights


double W_inner3(double beta1, double beta2, double k, double lambda){
    return exp(-1 * lambda) * pow(lambda, k) * w_k3(beta1, beta2, k, lambda)/boost::math::factorial<double>(k);
  }



double W3(double beta1, double beta2, double lambda, int sum_limit){
 double out = 0;
 for (int i = 0; i < sum_limit+1; i++){
    out += W_inner3(beta1, beta2, i,lambda);
 }
 return out;
}


//' Probability mass function of the three-parameter EWP
//'
//' @param x vector of (positive integer) quantiles.
//' @param lambda centrality parameter
//' @param beta1 lower-tail dispersion parameter
//' @param beta2 upper tail dispersion parameter
//' @param sum_limit summation limit for the normalizing factor
//'
//' @return a probability mass
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dewp3_cpp(Rcpp::IntegerVector x, double lambda, double beta1, double beta2, int sum_limit){
  int n = x.size();
  //Rcout << "The receivec size of x : " << n << "\n";
  Rcpp::NumericVector res(n);
  for (int i=0; i<n; i++){
    res[i] = exp(-1 * lambda) * pow(lambda, x[i]) * w_k3(beta1, beta2, x[i], lambda)/(W3(beta1, beta2, lambda, sum_limit)*boost::math::factorial<double>(x[i]));
  }
  return res;
}

// non-vectorized form because I failed to understand how to recast Rcpp::Vector subsets
// [[Rcpp::export]]
double dewp3_cpp_nv(int x, double lambda, double beta1, double beta2, int sum_limit){
  //int n = x.size();
  //Rcout << "The receivec size of x : " << n << "\n";
  //Rcpp::NumericVector res(n);
  //for (int i=0; i<n; i++){
    //res[i] =
    return exp(-1 * lambda) * pow(lambda, x) * w_k3(beta1, beta2, x, lambda)/(W3(beta1, beta2, lambda, sum_limit)*boost::math::factorial<double>(x));
  //}
  //return res;
}


// [[Rcpp::export]]
double pllik3_part_cpp(Rcpp::IntegerVector X, Rcpp::NumericVector lambda, double beta1, double beta2, int sum_limit){

  Rcpp::NumericVector ll(X.size());
  for (int i = 0; i < X.size(); i++){
    // Rcout << "The value of i : " << i << "\n";
    // Rcout << "The value of ll : " << ll << "\n";
    // Rcout << "The value of X : " << X << "\n";
    // Rcout << "The value of X[i] : " << X[i] << "\n";
    //Rcout << "The size of X[i] : " << X[i].size() << "\n";
    // Rcout << "Type of X" << typeid(X).name() << '\n';
    // Rcout << "Type of X[i]" << typeid(X[i]).name() << '\n';
    // Rcout << "The value of lambda : " << lambda << "\n";
    // Rcout << "The value of lambda[i] : " << lambda[i] << "\n";
    // Rcout << "Type of lambda[i]" << typeid(lambda[i]).name() << '\n';
    //Rcout << "raw calc:" << dewp3_cpp(as<IntegerVector>(X[i]),lambda(i),beta1,beta2) << "\n";
    ll[i] = log(dewp3_cpp_nv(X(i),lambda(i),beta1,beta2,sum_limit));//ugly? hack to force type conversion from Rcpp::NumericVector to std::double
    //Rcout << "The value of ll : " << ll << "\n";
  }
  return(-1*sum(ll));
}
