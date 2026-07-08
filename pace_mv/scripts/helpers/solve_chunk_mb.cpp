// solve_chunk_mb.cpp -- C++ Eigen + OpenMP per-gene Cholesky kernel for the
// chunked multi-block PACE-MV solve. Drop-in replacement for the inner
// per-gene loop in .solve_genes_chunk_multiblock().
//
// Design:
//   * R precomputes the chunk-level tensors (Stage 1 + Stage 2 in
//     .solve_genes_chunk_multiblock) and reshapes them to 2D matrices.
//   * R calls solve_chunk_mb_cpp() once per chunk; C++ extracts all data
//     to native std::vector<MatrixXd> structures BEFORE the parallel loop
//     (R API is not thread-safe inside OpenMP).
//   * Per-gene loop runs in parallel via OpenMP; each thread builds A, b
//     from precomputed tensors, runs Eigen LLT, computes Ainv diagonal.
//   * Output: 3 dense matrices (B, U, Ainv_diag), one column per gene.

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp,cpp14)]]

using namespace Eigen;
using namespace Rcpp;

// Per-(cross-block) cached layout
struct CrossBlock {
  int K_t_1, K_t_2, K_g_1, K_g_2;
  int g1_zero, g2_zero;        // 0-indexed group indices
  int co1, co2;
  MatrixXd arr;                // (K_t_1 * K_t_2) × G_chunk
};

// [[Rcpp::export]]
List solve_chunk_mb_cpp(
    const Eigen::Map<MatrixXd>& X_fixed,           // n × p
    const Eigen::Map<MatrixXd>& w_chunk,            // n × G_chunk
    const Eigen::Map<MatrixXd>& z_chunk,            // n × G_chunk
    const Eigen::Map<MatrixXd>& lam_diag_chunk,    // q × G_chunk
    int q_total,
    const List& blocks,                              // per-block: list with col_offset, K_terms, K_groups
    const List& ZtWZ_within,                         // per-block: list of K_t² × G_chunk matrices (one per group)
    const List& XtWZ_within,                         // per-block: list of (p*K_t) × G_chunk matrices (one per group)
    const List& ZtWz_within,                         // per-block: list of K_t × G_chunk matrices (one per group)
    const List& cross_blocks,                        // list of cross-block entries
    int n_threads
) {
  int n = X_fixed.rows();
  int p = X_fixed.cols();
  int G_chunk = w_chunk.cols();
  int B_n = blocks.size();

  // ---- Pre-extract block metadata (thread-safe before OMP) ----
  std::vector<int> col_offsets(B_n), K_t_v(B_n), K_g_v(B_n);
  for (int b = 0; b < B_n; b++) {
    List blk = blocks[b];
    col_offsets[b] = as<int>(blk["col_offset"]);
    K_t_v[b]       = as<int>(blk["K_terms"]);
    K_g_v[b]       = as<int>(blk["K_groups"]);
  }

  // ---- Pre-extract per-(block, group) tensors ----
  std::vector<std::vector<MatrixXd>> ZtWZ_data(B_n);
  std::vector<std::vector<MatrixXd>> XtWZ_data(B_n);
  std::vector<std::vector<MatrixXd>> ZtWz_data(B_n);
  for (int b = 0; b < B_n; b++) {
    List ZtWZ_b = ZtWZ_within[b];
    List XtWZ_b = XtWZ_within[b];
    List ZtWz_b = ZtWz_within[b];
    int K_g_b = K_g_v[b];
    ZtWZ_data[b].reserve(K_g_b);
    XtWZ_data[b].reserve(K_g_b);
    ZtWz_data[b].reserve(K_g_b);
    for (int g = 0; g < K_g_b; g++) {
      ZtWZ_data[b].push_back(as<MatrixXd>(ZtWZ_b[g]));
      XtWZ_data[b].push_back(as<MatrixXd>(XtWZ_b[g]));
      ZtWz_data[b].push_back(as<MatrixXd>(ZtWz_b[g]));
    }
  }

  // ---- Pre-extract cross blocks ----
  std::vector<CrossBlock> cross_data;
  cross_data.reserve(cross_blocks.size());
  for (int cb_idx = 0; cb_idx < cross_blocks.size(); cb_idx++) {
    List cb = cross_blocks[cb_idx];
    CrossBlock cd;
    cd.K_t_1   = as<int>(cb["K_t_1"]);
    cd.K_t_2   = as<int>(cb["K_t_2"]);
    cd.K_g_1   = as<int>(cb["K_g_1"]);
    cd.K_g_2   = as<int>(cb["K_g_2"]);
    cd.g1_zero = as<int>(cb["g1"]) - 1;
    cd.g2_zero = as<int>(cb["g2"]) - 1;
    cd.co1     = as<int>(cb["col_offset_1"]);
    cd.co2     = as<int>(cb["col_offset_2"]);
    cd.arr     = as<MatrixXd>(cb["array"]);
    cross_data.push_back(std::move(cd));
  }

  // ---- Output matrices ----
  MatrixXd B_out         = MatrixXd::Constant(p, G_chunk, NA_REAL);
  MatrixXd U_out         = MatrixXd::Constant(q_total, G_chunk, NA_REAL);
  MatrixXd Ainv_diag_out = MatrixXd::Constant(p + q_total, G_chunk, NA_REAL);

  // ---- Parallel per-gene loop ----
#ifdef _OPENMP
  #pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
  for (int gi = 0; gi < G_chunk; gi++) {
    // Per-gene weights and working response
    VectorXd w_g = w_chunk.col(gi);
    VectorXd z_g = z_chunk.col(gi);

    // XtWX (p × p), XtWz (p)
    MatrixXd Xw      = X_fixed.array().colwise() * w_g.array();
    MatrixXd XtWX   = X_fixed.transpose() * Xw;
    VectorXd XtWz   = X_fixed.transpose() * (w_g.cwiseProduct(z_g));

    // Allocate ZtWZ, XtWZ, ZtWz
    MatrixXd ZtWZ = MatrixXd::Zero(q_total, q_total);
    MatrixXd XtWZ = MatrixXd::Zero(p, q_total);
    VectorXd ZtWz = VectorXd::Zero(q_total);

    // Within-block contributions
    for (int b = 0; b < B_n; b++) {
      int K_t = K_t_v[b];
      int K_g = K_g_v[b];
      int co  = col_offsets[b];

      const std::vector<MatrixXd>& ZtWZ_b = ZtWZ_data[b];
      const std::vector<MatrixXd>& XtWZ_b = XtWZ_data[b];
      const std::vector<MatrixXd>& ZtWz_b = ZtWz_data[b];

      for (int g_idx = 0; g_idx < K_g; g_idx++) {
        const MatrixXd& Zb_g = ZtWZ_b[g_idx];   // K_t² × G_chunk
        const MatrixXd& Xb_g = XtWZ_b[g_idx];   // (p*K_t) × G_chunk
        const MatrixXd& zb_g = ZtWz_b[g_idx];   // K_t × G_chunk

        // Compute global column indices for this (block, group)
        // R's column order: col_offset + (t-1)*K_g + g (1-indexed)
        // 0-indexed: co + t * K_g + g_idx for t in 0..K_t-1
        for (int t1 = 0; t1 < K_t; t1++) {
          int col1 = co + t1 * K_g + g_idx;
          // ZtWZ block diagonal entries (within group g_idx, term pair t1, t2)
          // In Zb_g, column gi holds the K_t² entries for this gene in column-major order
          // i.e., Zb_g(t1 + t2 * K_t, gi) corresponds to entry [t1, t2] of the K_t × K_t block
          for (int t2 = 0; t2 < K_t; t2++) {
            int col2 = co + t2 * K_g + g_idx;
            ZtWZ(col1, col2) = Zb_g(t1 + t2 * K_t, gi);
          }
          // XtWZ: Xb_g column gi holds (p × K_t) entries column-major
          // entry [pi, t1] at position pi + t1 * p
          for (int pi = 0; pi < p; pi++) {
            XtWZ(pi, col1) = Xb_g(pi + t1 * p, gi);
          }
          ZtWz(col1) = zb_g(t1, gi);
        }
      }
    }

    // Cross-block contributions (b1 < b2)
    for (const CrossBlock& cd : cross_data) {
      int K_t_1 = cd.K_t_1, K_t_2 = cd.K_t_2;
      // Column indices in the global q-vector
      // cols1[t1] = co1 + t1 * K_g_1 + g1_zero
      // cols2[t2] = co2 + t2 * K_g_2 + g2_zero
      for (int t1 = 0; t1 < K_t_1; t1++) {
        int c1 = cd.co1 + t1 * cd.K_g_1 + cd.g1_zero;
        for (int t2 = 0; t2 < K_t_2; t2++) {
          int c2 = cd.co2 + t2 * cd.K_g_2 + cd.g2_zero;
          // arr column gi: K_t_1 × K_t_2 flat (column-major)
          // entry [t1, t2] at position t1 + t2 * K_t_1
          double val = cd.arr(t1 + t2 * K_t_1, gi);
          ZtWZ(c1, c2) = val;
          ZtWZ(c2, c1) = val;  // symmetric
        }
      }
    }

    // Add 1/tau penalty to ZtWZ diagonal (lam_diag_chunk is 1/tau)
    for (int i = 0; i < q_total; i++) {
      ZtWZ(i, i) += lam_diag_chunk(i, gi);
    }

    // Build full A = [XtWX, XtWZ; XtWZ', ZtWZ]
    MatrixXd A(p + q_total, p + q_total);
    A.topLeftCorner(p, p)             = XtWX;
    A.topRightCorner(p, q_total)      = XtWZ;
    A.bottomLeftCorner(q_total, p)    = XtWZ.transpose();
    A.bottomRightCorner(q_total, q_total) = ZtWZ;

    VectorXd b(p + q_total);
    b.head(p)         = XtWz;
    b.tail(q_total)   = ZtWz;

    // Cholesky-based solve
    LLT<MatrixXd> chol(A);
    if (chol.info() != Success) continue;  // leave NA in output

    VectorXd sol = chol.solve(b);
    B_out.col(gi) = sol.head(p);
    U_out.col(gi) = sol.tail(q_total);

    // Diagonal of A^-1 via solve against identity
    MatrixXd I = MatrixXd::Identity(p + q_total, p + q_total);
    Ainv_diag_out.col(gi) = chol.solve(I).diagonal();
  }

  return List::create(
    Named("B") = B_out,
    Named("U") = U_out,
    Named("Ainv_diag") = Ainv_diag_out
  );
}
