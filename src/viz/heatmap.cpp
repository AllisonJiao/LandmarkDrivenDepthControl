#include "heatmap.h"

#include <cmath>
#include <limits>

using namespace std;
using namespace Eigen;

// matplotlib RdBu_r sampled at 9 evenly spaced stops
const double kHeatmapTable[9][3] = {
    {0.020, 0.188, 0.380}, {0.129, 0.400, 0.675}, {0.263, 0.576, 0.765},
    {0.573, 0.773, 0.871}, {0.969, 0.969, 0.969}, {0.957, 0.647, 0.510},
    {0.839, 0.376, 0.302}, {0.698, 0.094, 0.169}, {0.404, 0.000, 0.122},
};

VectorXd displacementField(const MatrixXd &VRest, const MatrixXd &VDefault) {
  VectorXd dh = VectorXd::Zero(VRest.rows());
  if (VRest.rows() != VDefault.rows()) return dh;

  // dh_i = |h_i| - |h0_i|
  dh = VRest.col(2).array().abs() - VDefault.col(2).array().abs();

  return dh;
}

VectorXd pressureField(const VectorXd &dpRaw, const MatrixXd &VDefault,
                       double c, const vector<bool> &fixed) {
  const int n = VDefault.rows();
  VectorXd dp = VectorXd::Zero(n);
  if (dpRaw.size() != n || c == 0) return dp;

  // dp_i = sgn(h0_i) * dpRaw_i / c
  dp = VDefault.col(2).array().sign() * dpRaw.array() / c;

  // Outlines
  for (int i = 0; i < n; i++)
    if (fixed[i]) dp[i] = numeric_limits<double>::quiet_NaN();

  return dp;
}

void scalarsToColors(const VectorXd &f, double vmax, bool symLog,
                     double symLogF0, MatrixXd &C, double &vmaxUsed) {
  const int n = f.size();
  C.resize(n, 3);
  vmaxUsed = vmax;

  VectorXd g = f;
  const double f0 = symLogF0 > 0.0 ? symLogF0 : 1.0;

  // log scaling
  if (symLog) {
    g = f.array().sign() * (f.array().abs() / f0 + 1).log();
  }

  double range = 0.0;
  if (vmax > 0.0) {
    range = symLog ? log1p(vmax / f0) : vmax;
  } else {
    for (int i = 0; i < n; ++i) {
      if (isfinite(g[i])) range = max(range, abs(g[i]));
    }
    if (range == 0.0) range = 1.0;
  }
  vmaxUsed = symLog ? f0 * expm1(range) : range;

  const double invRange = 0.5 / range;
  for (int i = 0; i < n; i++) {
    if (!isfinite(g[i])) {
      C.row(i).setConstant(0.6);
      continue;
    }
    const double t = min(1.0, max(0.0, 0.5 + g[i] * invRange));
    const double tablePos = t * 8.0;
    const int lo = min(7, static_cast<int>(tablePos));
    const double alpha = tablePos - lo;
    for (int k = 0; k < 3; k++) {
      C(i, k) = (1.0 - alpha) * kHeatmapTable[lo][k] +
                alpha * kHeatmapTable[lo + 1][k];
    }
  }
}
