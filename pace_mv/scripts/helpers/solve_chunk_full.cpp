// solve_chunk_full.cpp -- Full C++ port of .solve_genes_chunk_multiblock.
//
// Subsumes Stage 1 (per-(block, group) within-block tensors), Stage 2
// (cross-block tensors over (g1, g2) intersections), and Stage 3 (per-gene
// Cholesky / WLS / diag-inverse). Single C++ entry point per chunk; OpenMP
// parallelism on the per-gene loop. Numerically equivalent to both the
// R-only path and the Stage-3-only solve_chunk_mb_cpp.
//
// Memory: tensors materialised per chunk (same as R path) but built via
// Eigen dgemm with no R-side allocations. At DKD scale (G_chunk=256,
// q=623, K_t=36, K_g=17) total ~85 MB per chunk.
//
// Stage 3 has two paths, selected by `stage3_mode`:
//   0 (auto, default): Schur-partitioned solve when B_n == 2; dense fallback
//                       otherwise.
//   1 (force dense)  : original (p+q)×(p+q) Cholesky path.
//   2 (force Schur)  : Schur-partitioned path; errors if B_n != 2.
// Schur exploits the fact that Z'WZ + Λ is block-diagonal *within* each RE
// block (every cell belongs to exactly one group per block) and that the
// only dense coupling is the cross-block (block 0 × block 1) corner.
// Cost drops from O((p+q)^3) per gene to O(K_g·K_t^3 + q1·q2^2 + p·q^2),
// and per-gene scratch memory drops by ~3x (no full A or full A^{-1}).
//
// Interior numeric precision is selected by `interior_precision`:
//   0 (default): all of Stage 1+2+3 in double.
//   1          : Stage 1+2+3 in float (single precision). The R caller
//                should run interior IRLS iterations with this on, then
//                flip to 0 for the final iter so se_B / se_U inherit
//                full double precision. ~1.5-1.7x speedup on Apple Silicon
//                NEON; numerical perturbation ~1e-6 per element which
//                damps out across IRLS updates.

#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <unordered_map>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp,cpp14)]]

using namespace Eigen;
using namespace Rcpp;

template <typename T>
struct CrossBlockTplt {
  int K_t_1, K_t_2, K_g_1, K_g_2;
  int g1, g2;
  int co1, co2;
  Eigen::Matrix<T, -1, -1> arr;  // (K_t_1 * K_t_2) × G_chunk
};

// Templated implementation. Inputs come in already cast to T.
// Outputs are returned in T precision; the wrapper casts back to double.
template <typename T>
List solve_chunk_full_impl(
    const Eigen::Matrix<T, -1, -1>& X_fixed,           // n × p
    const Eigen::Matrix<T, -1, -1>& w_chunk,            // n × G_chunk
    const Eigen::Matrix<T, -1, -1>& z_chunk,            // n × G_chunk
    const Eigen::Matrix<T, -1, -1>& lam_diag_chunk,    // q × G_chunk
    int q_total,
    const std::vector<int>& col_offsets,
    const std::vector<int>& K_t_v,
    const std::vector<int>& K_g_v,
    const std::vector<Eigen::Matrix<T, -1, -1>>& X_terms_vec,
    const std::vector<VectorXi>& cell_grp_vec,
    const std::vector<std::vector<std::vector<int>>>& cells_by_grp_vec,
    int n_threads,
    int stage3_mode,
    int B_n
) {
  using MatT = Eigen::Matrix<T, -1, -1>;
  using VecT = Eigen::Matrix<T, -1, 1>;

  int n = X_fixed.rows();
  int p = X_fixed.cols();
  int G_chunk = w_chunk.cols();

  // ============================================================
  // Stage 1: per-(block, group) within-block tensors
  // ============================================================
  std::vector<std::vector<MatT>> ZtWZ_w(B_n);
  std::vector<std::vector<MatT>> XtWZ_w(B_n);
  std::vector<std::vector<MatT>> ZtWz_w(B_n);

  for (int b = 0; b < B_n; b++) {
    int K_t = K_t_v[b];
    int K_g = K_g_v[b];
    const MatT& Xt_b = X_terms_vec[b];

    ZtWZ_w[b].resize(K_g);
    XtWZ_w[b].resize(K_g);
    ZtWz_w[b].resize(K_g);

#ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
    for (int g = 0; g < K_g; g++) {
      const std::vector<int>& idx = cells_by_grp_vec[b][g];
      int n_g = idx.size();
      if (n_g == 0) {
        ZtWZ_w[b][g] = MatT::Zero(K_t * K_t, G_chunk);
        XtWZ_w[b][g] = MatT::Zero(p * K_t, G_chunk);
        ZtWz_w[b][g] = MatT::Zero(K_t, G_chunk);
        continue;
      }

      MatT Xt_g(n_g, K_t);
      MatT X_g(n_g, p);
      for (int i = 0; i < n_g; i++) {
        int row = idx[i];
        Xt_g.row(i) = Xt_b.row(row);
        X_g.row(i)  = X_fixed.row(row);
      }
      // COLUMN-wise gather of w_g / z_g (sequential within each source column):
      // cache-friendly relative to the per-row gather above.
      MatT w_g(n_g, G_chunk), z_g(n_g, G_chunk);
      for (int j = 0; j < G_chunk; j++) {
        const T* ws = w_chunk.col(j).data(); const T* zs = z_chunk.col(j).data();
        T* wd = w_g.col(j).data(); T* zd = z_g.col(j).data();
        for (int i = 0; i < n_g; i++) { int row = idx[i]; wd[i]=ws[row]; zd[i]=zs[row]; }
      }

      MatT M_pair(n_g, K_t * K_t);
      for (int t2 = 0; t2 < K_t; t2++) {
        for (int t1 = 0; t1 < K_t; t1++) {
          M_pair.col(t2 * K_t + t1) =
            Xt_g.col(t1).cwiseProduct(Xt_g.col(t2));
        }
      }
      ZtWZ_w[b][g].noalias() = M_pair.transpose() * w_g;

      MatT M_xt(n_g, p * K_t);
      for (int t = 0; t < K_t; t++) {
        for (int pi = 0; pi < p; pi++) {
          M_xt.col(t * p + pi) =
            X_g.col(pi).cwiseProduct(Xt_g.col(t));
        }
      }
      XtWZ_w[b][g].noalias() = M_xt.transpose() * w_g;

      MatT WZ = w_g.cwiseProduct(z_g);
      ZtWz_w[b][g].noalias() = Xt_g.transpose() * WZ;
    }
  }

  // ============================================================
  // Stage 2: cross-block tensors (b1 < b2)
  // ============================================================
  std::vector<CrossBlockTplt<T>> cross_data;
  if (B_n >= 2) {
    for (int b1 = 0; b1 < B_n - 1; b1++) {
      int K_t_1 = K_t_v[b1];
      int K_g_1 = K_g_v[b1];
      const MatT& Xt_1 = X_terms_vec[b1];

      for (int b2 = b1 + 1; b2 < B_n; b2++) {
        int K_t_2 = K_t_v[b2];
        int K_g_2 = K_g_v[b2];
        const MatT& Xt_2 = X_terms_vec[b2];
        const VectorXi& grp2 = cell_grp_vec[b2];

        for (int g1 = 0; g1 < K_g_1; g1++) {
          const std::vector<int>& idx1 = cells_by_grp_vec[b1][g1];
          if (idx1.empty()) continue;

          std::unordered_map<int, std::vector<int>> buckets;
          for (int cell : idx1) {
            int g2_idx = grp2[cell] - 1;
            buckets[g2_idx].push_back(cell);
          }

          for (const auto& kv : buckets) {
            int g2 = kv.first;
            const std::vector<int>& idx = kv.second;
            int n_idx = idx.size();
            if (n_idx == 0) continue;

            MatT Xt_1_sub(n_idx, K_t_1);
            MatT Xt_2_sub(n_idx, K_t_2);
            MatT w_sub(n_idx, G_chunk);
            for (int i = 0; i < n_idx; i++) {
              Xt_1_sub.row(i) = Xt_1.row(idx[i]);
              Xt_2_sub.row(i) = Xt_2.row(idx[i]);
              w_sub.row(i)    = w_chunk.row(idx[i]);
            }

            MatT M_cross(n_idx, K_t_1 * K_t_2);
            for (int t2 = 0; t2 < K_t_2; t2++) {
              for (int t1 = 0; t1 < K_t_1; t1++) {
                M_cross.col(t2 * K_t_1 + t1) =
                  Xt_1_sub.col(t1).cwiseProduct(Xt_2_sub.col(t2));
              }
            }

            CrossBlockTplt<T> cb;
            cb.K_t_1 = K_t_1; cb.K_t_2 = K_t_2;
            cb.K_g_1 = K_g_1; cb.K_g_2 = K_g_2;
            cb.g1 = g1; cb.g2 = g2;
            cb.co1 = col_offsets[b1]; cb.co2 = col_offsets[b2];
            cb.arr.noalias() = M_cross.transpose() * w_sub;
            cross_data.push_back(std::move(cb));
          }
        }
      }
    }
  }

  // ============================================================
  // Stage 3: per-gene parallel solve (Schur or dense, dispatched per gene)
  // ============================================================
  // Output matrices in T (cast back to double in wrapper). NaN sentinel via
  // a quiet NaN (cast from double NA_REAL is fine because NaN bit pattern
  // is preserved across float<->double cast for the float canonical NaN).
  T NaN_T = std::numeric_limits<T>::quiet_NaN();
  MatT B_out         = MatT::Constant(p, G_chunk, NaN_T);
  MatT U_out         = MatT::Constant(q_total, G_chunk, NaN_T);
  MatT Ainv_diag_out = MatT::Constant(p + q_total, G_chunk, NaN_T);

  bool use_schur;
  if (stage3_mode == 1) {
    use_schur = false;
  } else if (stage3_mode == 2) {
    if (B_n != 2) stop("stage3_mode = 2 (force Schur) requires B_n == 2 (got %d)", B_n);
    use_schur = true;
  } else {
    use_schur = (B_n == 2);
  }

  int K_t1_g = 0, K_g1_g = 0, q1_g = 0, co0_g = 0;
  int K_t2_g = 0, K_g2_g = 0, q2_g = 0, co1_g = 0;
  if (use_schur) {
    K_t1_g = K_t_v[0]; K_g1_g = K_g_v[0]; q1_g = K_t1_g * K_g1_g; co0_g = col_offsets[0];
    K_t2_g = K_t_v[1]; K_g2_g = K_g_v[1]; q2_g = K_t2_g * K_g2_g; co1_g = col_offsets[1];
  }

#ifdef _OPENMP
  #pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
  for (int gi = 0; gi < G_chunk; gi++) {
    VecT w_g_v = w_chunk.col(gi);
    VecT z_g_v = z_chunk.col(gi);

    MatT Xw       = X_fixed.array().colwise() * w_g_v.array();
    MatT XtWX    = X_fixed.transpose() * Xw;
    VecT XtWz    = X_fixed.transpose() * (w_g_v.cwiseProduct(z_g_v));

    if (use_schur) {
      // ---------- Schur-partitioned Stage 3 (B_n == 2) ----------
      const int K_t1 = K_t1_g, K_g1 = K_g1_g, q1 = q1_g, co0 = co0_g;
      const int K_t2 = K_t2_g, K_g2 = K_g2_g, q2 = q2_g, co1 = co1_g;

      std::vector<LLT<MatT>> G11_chol;
      G11_chol.reserve(K_g1);
      bool ok = true;
      for (int g = 0; g < K_g1; g++) {
        MatT Bg(K_t1, K_t1);
        const MatT& Zg = ZtWZ_w[0][g];
        for (int t2 = 0; t2 < K_t1; t2++) {
          for (int t1 = 0; t1 < K_t1; t1++) {
            Bg(t1, t2) = Zg(t1 + t2 * K_t1, gi);
          }
        }
        for (int t = 0; t < K_t1; t++) {
          Bg(t, t) += lam_diag_chunk(co0 + t * K_g1 + g, gi);
        }
        G11_chol.emplace_back(Bg);
        if (G11_chol.back().info() != Success) { ok = false; break; }
      }
      if (!ok) continue;

      std::vector<MatT> G22_blocks(K_g2);
      for (int g = 0; g < K_g2; g++) {
        MatT Bg(K_t2, K_t2);
        const MatT& Zg = ZtWZ_w[1][g];
        for (int t2 = 0; t2 < K_t2; t2++) {
          for (int t1 = 0; t1 < K_t2; t1++) {
            Bg(t1, t2) = Zg(t1 + t2 * K_t2, gi);
          }
        }
        for (int t = 0; t < K_t2; t++) {
          Bg(t, t) += lam_diag_chunk(co1 + t * K_g2 + g, gi);
        }
        G22_blocks[g] = std::move(Bg);
      }

      MatT G12 = MatT::Zero(q1, q2);
      for (const CrossBlockTplt<T>& cd : cross_data) {
        int row_off = cd.g1 * K_t1;
        int col_off = cd.g2 * K_t2;
        for (int t2 = 0; t2 < K_t2; t2++) {
          for (int t1 = 0; t1 < K_t1; t1++) {
            G12(row_off + t1, col_off + t2) = cd.arr(t1 + t2 * K_t1, gi);
          }
        }
      }

      MatT Tmat(q1, q2);
      for (int g = 0; g < K_g1; g++) {
        int off = g * K_t1;
        Tmat.middleRows(off, K_t1) = G11_chol[g].solve(G12.middleRows(off, K_t1));
      }

      MatT SGC = -(G12.transpose() * Tmat);
      for (int g = 0; g < K_g2; g++) {
        SGC.block(g * K_t2, g * K_t2, K_t2, K_t2) += G22_blocks[g];
      }
      LLT<MatT> SGC_chol(SGC);
      if (SGC_chol.info() != Success) continue;

      MatT XtWZ_perm(p, q_total);
      for (int g = 0; g < K_g1; g++) {
        const MatT& Xg = XtWZ_w[0][g];
        int col_off = g * K_t1;
        for (int t = 0; t < K_t1; t++) {
          for (int pi = 0; pi < p; pi++) {
            XtWZ_perm(pi, col_off + t) = Xg(pi + t * p, gi);
          }
        }
      }
      for (int g = 0; g < K_g2; g++) {
        const MatT& Xg = XtWZ_w[1][g];
        int col_off = q1 + g * K_t2;
        for (int t = 0; t < K_t2; t++) {
          for (int pi = 0; pi < p; pi++) {
            XtWZ_perm(pi, col_off + t) = Xg(pi + t * p, gi);
          }
        }
      }

      VecT ZtWz_perm(q_total);
      for (int g = 0; g < K_g1; g++) {
        const MatT& zg = ZtWz_w[0][g];
        for (int t = 0; t < K_t1; t++) ZtWz_perm(g * K_t1 + t) = zg(t, gi);
      }
      for (int g = 0; g < K_g2; g++) {
        const MatT& zg = ZtWz_w[1][g];
        for (int t = 0; t < K_t2; t++) ZtWz_perm(q1 + g * K_t2 + t) = zg(t, gi);
      }

      VecT ws(q1);
      for (int g = 0; g < K_g1; g++) {
        int off = g * K_t1;
        ws.segment(off, K_t1) = G11_chol[g].solve(ZtWz_perm.segment(off, K_t1));
      }
      VecT r2 = ZtWz_perm.tail(q2) - G12.transpose() * ws;
      VecT q_eff_2 = SGC_chol.solve(r2);
      VecT q_eff_1 = ws - Tmat * q_eff_2;

      MatT ws_p(q1, p);
      for (int g = 0; g < K_g1; g++) {
        int off = g * K_t1;
        ws_p.middleRows(off, K_t1) =
          G11_chol[g].solve(XtWZ_perm.block(0, off, p, K_t1).transpose());
      }
      MatT r2_p = XtWZ_perm.rightCols(q2).transpose() - G12.transpose() * ws_p;
      MatT Q2 = SGC_chol.solve(r2_p);
      MatT Q1 = ws_p - Tmat * Q2;

      MatT S_beta = XtWX;
      S_beta.noalias() -= XtWZ_perm.leftCols(q1)  * Q1;
      S_beta.noalias() -= XtWZ_perm.rightCols(q2) * Q2;

      VecT b_eff = XtWz;
      b_eff.noalias() -= XtWZ_perm.leftCols(q1)  * q_eff_1;
      b_eff.noalias() -= XtWZ_perm.rightCols(q2) * q_eff_2;

      LLT<MatT> S_beta_chol(S_beta);
      if (S_beta_chol.info() != Success) continue;
      VecT beta = S_beta_chol.solve(b_eff);

      VecT u_perm_1 = q_eff_1 - Q1 * beta;
      VecT u_perm_2 = q_eff_2 - Q2 * beta;

      B_out.col(gi) = beta;
      for (int g = 0; g < K_g1; g++) {
        int perm_off = g * K_t1;
        for (int t = 0; t < K_t1; t++) {
          U_out(co0 + t * K_g1 + g, gi) = u_perm_1(perm_off + t);
        }
      }
      for (int g = 0; g < K_g2; g++) {
        int perm_off = g * K_t2;
        for (int t = 0; t < K_t2; t++) {
          U_out(co1 + t * K_g2 + g, gi) = u_perm_2(perm_off + t);
        }
      }

      MatT S_beta_inv = S_beta_chol.solve(MatT::Identity(p, p));
      for (int i = 0; i < p; i++) Ainv_diag_out(i, gi) = S_beta_inv(i, i);

      MatT K = SGC_chol.solve(Tmat.transpose());

      VecT diag_G11inv(q1);
      for (int g = 0; g < K_g1; g++) {
        MatT Linv = G11_chol[g].matrixL().solve(MatT::Identity(K_t1, K_t1));
        int perm_off = g * K_t1;
        for (int t = 0; t < K_t1; t++)
          diag_G11inv(perm_off + t) = Linv.col(t).squaredNorm();
      }

      MatT L_SGC_inv = SGC_chol.matrixL().solve(MatT::Identity(q2, q2));

      for (int g = 0; g < K_g1; g++) {
        int perm_off = g * K_t1;
        for (int t = 0; t < K_t1; t++) {
          int i = perm_off + t;
          T d_Ginv  = diag_G11inv(i) + Tmat.row(i).dot(K.col(i));
          T d_QSinv = Q1.row(i) * S_beta_inv * Q1.row(i).transpose();
          int orig_col = co0 + t * K_g1 + g;
          Ainv_diag_out(p + orig_col, gi) = d_Ginv + d_QSinv;
        }
      }
      for (int g = 0; g < K_g2; g++) {
        int perm_off = g * K_t2;
        for (int t = 0; t < K_t2; t++) {
          int i = perm_off + t;
          T d_Ginv  = L_SGC_inv.col(i).squaredNorm();
          T d_QSinv = Q2.row(i) * S_beta_inv * Q2.row(i).transpose();
          int orig_col = co1 + t * K_g2 + g;
          Ainv_diag_out(p + orig_col, gi) = d_Ginv + d_QSinv;
        }
      }

      continue;
    }

    // ---------- Dense (legacy) Stage 3 path ----------
    MatT ZtWZ = MatT::Zero(q_total, q_total);
    MatT XtWZ = MatT::Zero(p, q_total);
    VecT ZtWz = VecT::Zero(q_total);

    for (int b = 0; b < B_n; b++) {
      int K_t = K_t_v[b];
      int K_g = K_g_v[b];
      int co  = col_offsets[b];
      const std::vector<MatT>& Zb = ZtWZ_w[b];
      const std::vector<MatT>& Xb = XtWZ_w[b];
      const std::vector<MatT>& zb = ZtWz_w[b];

      for (int g_idx = 0; g_idx < K_g; g_idx++) {
        const MatT& Zb_g = Zb[g_idx];
        const MatT& Xb_g = Xb[g_idx];
        const MatT& zb_g = zb[g_idx];

        for (int t1 = 0; t1 < K_t; t1++) {
          int col1 = co + t1 * K_g + g_idx;
          for (int t2 = 0; t2 < K_t; t2++) {
            int col2 = co + t2 * K_g + g_idx;
            ZtWZ(col1, col2) = Zb_g(t1 + t2 * K_t, gi);
          }
          for (int pi = 0; pi < p; pi++) {
            XtWZ(pi, col1) = Xb_g(pi + t1 * p, gi);
          }
          ZtWz(col1) = zb_g(t1, gi);
        }
      }
    }

    for (const CrossBlockTplt<T>& cd : cross_data) {
      for (int t1 = 0; t1 < cd.K_t_1; t1++) {
        int c1 = cd.co1 + t1 * cd.K_g_1 + cd.g1;
        for (int t2 = 0; t2 < cd.K_t_2; t2++) {
          int c2 = cd.co2 + t2 * cd.K_g_2 + cd.g2;
          T val = cd.arr(t1 + t2 * cd.K_t_1, gi);
          ZtWZ(c1, c2) = val;
          ZtWZ(c2, c1) = val;
        }
      }
    }

    for (int i = 0; i < q_total; i++) ZtWZ(i, i) += lam_diag_chunk(i, gi);

    MatT A(p + q_total, p + q_total);
    A.topLeftCorner(p, p)             = XtWX;
    A.topRightCorner(p, q_total)      = XtWZ;
    A.bottomLeftCorner(q_total, p)    = XtWZ.transpose();
    A.bottomRightCorner(q_total, q_total) = ZtWZ;

    VecT b(p + q_total);
    b.head(p)         = XtWz;
    b.tail(q_total)   = ZtWz;

    LLT<MatT> chol(A);
    if (chol.info() != Success) continue;

    VecT sol = chol.solve(b);
    B_out.col(gi) = sol.head(p);
    U_out.col(gi) = sol.tail(q_total);

    int N = p + q_total;
    MatT Linv = chol.matrixL().solve(MatT::Identity(N, N));
    Ainv_diag_out.col(gi) = Linv.colwise().squaredNorm();
  }

  return List::create(
    Named("B")          = B_out.template cast<double>(),
    Named("U")          = U_out.template cast<double>(),
    Named("Ainv_diag")  = Ainv_diag_out.template cast<double>()
  );
}


// [[Rcpp::export]]
List solve_chunk_full_cpp(
    const Eigen::Map<MatrixXd>& X_fixed,           // n × p
    const Eigen::Map<MatrixXd>& w_chunk,            // n × G_chunk
    const Eigen::Map<MatrixXd>& z_chunk,            // n × G_chunk
    const Eigen::Map<MatrixXd>& lam_diag_chunk,    // q × G_chunk
    int q_total,
    const List& blocks,
    const List& X_terms_list,
    const List& cells_by_grp_list,
    const List& cell_grp_list,
    int n_threads,
    int stage3_mode = 0,
    int interior_precision = 0
) {
  int B_n = blocks.size();

  // ---- Pre-extract block metadata + design data (shared across precisions) ----
  std::vector<int> col_offsets(B_n), K_t_v(B_n), K_g_v(B_n);
  for (int b = 0; b < B_n; b++) {
    List blk = blocks[b];
    col_offsets[b] = as<int>(blk["col_offset"]);
    K_t_v[b]       = as<int>(blk["K_terms"]);
    K_g_v[b]       = as<int>(blk["K_groups"]);
  }

  std::vector<VectorXi> cell_grp_vec(B_n);
  std::vector<std::vector<std::vector<int>>> cells_by_grp_vec(B_n);
  for (int b = 0; b < B_n; b++) {
    cell_grp_vec[b] = as<VectorXi>(cell_grp_list[b]);
    List cbg = cells_by_grp_list[b];
    cells_by_grp_vec[b].resize(K_g_v[b]);
    for (int g = 0; g < K_g_v[b]; g++) {
      IntegerVector iv = cbg[g];
      cells_by_grp_vec[b][g].resize(iv.size());
      for (int i = 0; i < iv.size(); i++) {
        cells_by_grp_vec[b][g][i] = iv[i] - 1;
      }
    }
  }

  if (interior_precision == 0) {
    // ---- double path: views directly over caller-owned MatrixXd buffers ----
    std::vector<MatrixXd> X_terms_vec(B_n);
    for (int b = 0; b < B_n; b++) {
      X_terms_vec[b] = as<MatrixXd>(X_terms_list[b]);
    }
    return solve_chunk_full_impl<double>(
      X_fixed, w_chunk, z_chunk, lam_diag_chunk,
      q_total, col_offsets, K_t_v, K_g_v,
      X_terms_vec, cell_grp_vec, cells_by_grp_vec,
      n_threads, stage3_mode, B_n);
  } else if (interior_precision == 1) {
    // ---- float path: cast inputs once at the boundary ----
    MatrixXf X_fixed_f        = X_fixed.cast<float>();
    MatrixXf w_chunk_f        = w_chunk.cast<float>();
    MatrixXf z_chunk_f        = z_chunk.cast<float>();
    MatrixXf lam_diag_chunk_f = lam_diag_chunk.cast<float>();
    std::vector<MatrixXf> X_terms_vec_f(B_n);
    for (int b = 0; b < B_n; b++) {
      X_terms_vec_f[b] = as<MatrixXd>(X_terms_list[b]).cast<float>();
    }
    return solve_chunk_full_impl<float>(
      X_fixed_f, w_chunk_f, z_chunk_f, lam_diag_chunk_f,
      q_total, col_offsets, K_t_v, K_g_v,
      X_terms_vec_f, cell_grp_vec, cells_by_grp_vec,
      n_threads, stage3_mode, B_n);
  } else {
    stop("interior_precision must be 0 (double) or 1 (float), got %d", interior_precision);
  }
}
