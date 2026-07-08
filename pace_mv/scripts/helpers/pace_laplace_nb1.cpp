// pace_laplace_nb1.cpp — TMB Laplace inner solve for one gene under PACE's
// joint random-effect structure.
//
// Model:
//   y_i ~ NB1(mu_i, alpha)     where var(y_i) = mu_i * (1 + alpha)
//   log(mu_i) = X_i beta + Z_i u + offset_i
//   u_k ~ N(0, 1 / tau_inv_k)  independent across k
//
// We profile u out via Laplace and optimise beta (and optionally log_alpha).
// PACE's outer loop updates tau_inv between calls; alpha is updated separately
// after each outer iteration via NB1 MLE on residuals (so log_alpha is fixed
// during this inner solve).

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // -------- data --------
  DATA_VECTOR(y);            // (n) observed counts
  DATA_MATRIX(X);            // (n x p) fixed-effect design (dense, small p)
  DATA_SPARSE_MATRIX(Z);     // (n x q) random-effect design (sparse, block)
  DATA_VECTOR(offset);       // (n) log offset (library + E^tech)
  DATA_VECTOR(tau_inv);      // (q) per-column ridge precision = 1 / tau_g
  DATA_SCALAR(alpha);        // NB1 dispersion, FIXED (PACE updates outside)

  // -------- parameters --------
  PARAMETER_VECTOR(beta);    // (p)
  PARAMETER_VECTOR(u);       // (q)  RANDOM

  Type nll = Type(0);

  // -------- linear predictor --------
  vector<Type> eta = X * beta + Z * u + offset;
  vector<Type> mu  = exp(eta);

  // -------- NB1 log-likelihood --------
  // var = mu * (1 + alpha)  -> classic NB(size = mu/alpha, mu = mu)
  // TMB's dnbinom2(x, mu, var) uses the (mean, variance) parameterisation
  Type one_plus_alpha = Type(1) + alpha;
  for (int i = 0; i < y.size(); i++) {
    Type var_i = mu(i) * one_plus_alpha;
    nll -= dnbinom2(y(i), mu(i), var_i, true);
  }

  // -------- Ridge prior on u --------
  // u_k ~ N(0, tau_k)  =>  -log p(u) = 0.5 * sum(tau_inv * u^2) + const
  // Drop the normalising constant since it does not affect beta/u optimisation
  for (int k = 0; k < u.size(); k++) {
    nll += Type(0.5) * tau_inv(k) * u(k) * u(k);
  }

  return nll;
}
