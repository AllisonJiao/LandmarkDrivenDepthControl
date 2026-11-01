#include "landmarksolver.h"

#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>

#include "json.hpp"

using namespace std;
using namespace Eigen;

SparseMatrix<double> buildEnergy(const SparseMatrix<double> &L,
                                 const SparseMatrix<double> &M,
                                 const SparseMatrix<double> &Minv) {
  // Alternative triharmonic. Smoother but less local
  // return L.transpose() * Minv * (-1) * L * Minv * L.transpose();
  return L.transpose() * Minv * L;
}

const char *handleMethodName(HandleMethod method) {
  switch (method) {
    case HandleMethod::G1: return "G1";
    case HandleMethod::G2: return "G2";
    case HandleMethod::P1: return "P1";
    case HandleMethod::P2: return "P2";
    case HandleMethod::P3: return "P3";
  }
  return "?";
}

void LandmarkSolver::prepare(const MatrixXd &VDefault, const MatrixXi &F,
                             const vector<vector<int>> &bnds) {
  this->VDefault = VDefault;
  this->F = F;
  ready = factored = false;
  numConstraints = 0;
  du.resize(0); isFixed.clear();
  dp.resize(0);
  if (VDefault.rows() == 0 || F.rows() == 0) return;

  // operators of the flat mesh, as in the inflation step of reconstruction()
  MatrixXd VFlat = VDefault;
  VFlat.col(2).setZero();
  SparseMatrix<double> Minv;
  igl::cotmatrix(VFlat, F, L);
  igl::massmatrix(VFlat, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, Minv);
  Q = buildEnergy(L, M, Minv);

  const VectorXd h = VDefault.col(2);
  u0 = h.array().sign() * h.array().square();

  set<int> fixedSet;
  for (const auto &bnd : bnds)
    for (int i : bnd) fixedSet.insert(i);
  fixedIdx.resize(fixedSet.size());
  int c = 0;
  for (int i : fixedSet) fixedIdx(c++) = i;
  isFixed.assign(VDefault.rows(), false);
  for (int i: fixedSet) isFixed[i] = true;

  // TODO(human) [prepare]: build what forward() and implied() need, once per
  // mesh: the interior vertices (all vertices not in fixedIdx) with a map
  // from vertex id to interior index, L_II = L restricted to the interior,
  // and a factorization of -L_II. Return without setting ready if the
  // factorization fails.

  ready = true;
}

void LandmarkSolver::setMethod(HandleMethod method, double radius) {
  this->method = method;
  this->radius = radius;
}

bool LandmarkSolver::setConstraints(const vector<Landmark> &lms,
                                    double weight) {
  factored = false;
  this->weight = weight;
  numConstraints = lms.size();
  if (!ready || lms.empty()) return false;

  // One row per landmark with a single 1 at its vertex, so (B u)_k reads u
  // there and B * u0 is the default u at the landmarks.
  const int n = VDefault.rows();
  vector<Triplet<double>> triplets;
  for (int k = 0; k < lms.size(); k++)
    triplets.emplace_back(k, lms[k].vertexId, 1.0);
  B = SparseMatrix<double>(lms.size(), n);
  B.setFromTriplets(triplets.begin(), triplets.end());

  // TODO(human) [D6 weight]: the baselines hit their targets exactly by
  // construction, so `weight` (soft vs hard) only affects P3 here. The UI side
  // of this decision is in src/ui/main.js (#landmarkMethod handler).
  if (method == HandleMethod::P3)
    factored = setConstraintsP3();
  else
    factored = setConstraintsBaseline(lms);
  if (!factored)
    cerr << "landmark solve (" << handleMethodName(method)
         << "): setup failed" << endl;
  return factored;
}

bool LandmarkSolver::setConstraintsP3() {
  if (isinf(weight))
    return igl::min_quad_with_fixed_precompute(Q, fixedIdx, B, false, data);
  SparseMatrix<double> A = Q + weight * SparseMatrix<double>(B.transpose()) * B;
  return igl::min_quad_with_fixed_precompute(A, fixedIdx,
                                             SparseMatrix<double>(), true, data);
}

bool LandmarkSolver::setConstraintsBaseline(const vector<Landmark> &lms) {
  // Everything that depends on the landmark set, the method or the radius,
  // but not on the factors, belongs here: solveBaseline() runs on every drag
  // frame and should only rescale what is computed here.
  //
  // Decisions that apply to several methods:
  // TODO(human) [D2 sheet]: the front and back sheets share x,y, so a disk
  //   around the handle (G2, P2) also covers the other sheet. Restrict the
  //   falloff to the handle's sheet, e.g. sign(u0(j)) == lms[k].side for a
  //   one-part shape, or connected components without the fixed vertices for
  //   several parts.
  // TODO(human) [D3 radius]: `radius` is in canvas pixels, like
  //   VDefault.col(0) and col(1). Its range and default are set by the
  //   #landmarkRadius slider in src/ui/myshell.html.
  // TODO(human) [D5 anchors]: an anchor is a landmark with a fixed factor
  //   (lms[k].anchor). With several landmarks, G2/P1/P2 need one basis column
  //   per landmark and a k x k solve so that each one hits its own target.
  // [D7 outline] pickLandmark and rebindLandmarks never put a handle on a
  //   fixed vertex: its target r would be 0 and forward() is 0 there, which
  //   makes the k x k system singular. A check of isFixed here would only
  //   guard landmarks that come from somewhere else.
  switch (method) {
    case HandleMethod::G1:
      // TODO(human) [G1]: du = r at each handle vertex, 0 elsewhere. Store
      // the handle vertex ids.
      break;
    case HandleMethod::G2:
      // TODO(human) [G2]: falloff w_j = (1 - (d_j / radius)^2)^3 for
      // d_j < radius around each handle (2D distance on VDefault), with
      // du_j = u0_j * sum_k c_k w_kj. Precompute w and the k x k matrix that
      // maps the coefficients c to du at the handles.
      break;
    case HandleMethod::P1:
      // TODO(human) [P1]: basis dp = e_i per handle; precompute
      // g_k = forward(e_i) and the k x k matrix (B g).
      break;
    case HandleMethod::P2:
      // TODO(human) [P2]: basis dp = falloff w around each handle; precompute
      // g_k = forward(w_k) and the k x k matrix (B g).
      break;
    case HandleMethod::P3:
      return false;  // handled by setConstraintsP3()
  }
  // placeholder until the cases above are implemented
  cerr << handleMethodName(method)
       << ": not implemented yet, showing the default inflation" << endl;
  return true;
}

bool LandmarkSolver::solve(const vector<Landmark> &lms, MatrixXd &VRestOut) {
  VRestOut = VDefault;
  const int n = VDefault.rows();
  du = dp = VectorXd::Zero(n);
  if (!ready || lms.empty()) return true;  // default inflation
  if (!factored || lms.size() != numConstraints) return false;

  const auto tStart = chrono::high_resolution_clock::now();

  // Targets in u-space. Scaling a height by f scales u by f^2, and taking the
  // default u at the landmark keeps factor 1 an exact no-op.
  const VectorXd uAtLandmarks = B * u0;
  VectorXd r(lms.size());
  for (int k = 0; k < lms.size(); k++)
    r(k) = (lms[k].factor * lms[k].factor - 1) * uAtLandmarks(k);

  VectorXd duNew = VectorXd::Zero(n);
  const bool ok = method == HandleMethod::P3 ? solveP3(r, duNew)
                                             : solveBaseline(lms, r, duNew);
  if (!ok) return false;
  du = duNew;
  dp = implied(du);

  const VectorXd u = u0 + du;
  VRestOut.col(2) = u.array().sign() * u.array().abs().sqrt();

  const double ms = chrono::duration<double, milli>(
                        chrono::high_resolution_clock::now() - tStart)
                        .count();
  cout << "landmark solve " << handleMethodName(method) << ": " << ms
       << " ms (" << n << " vertices, " << lms.size() << " landmarks)"
       << endl;
  return true;
}

bool LandmarkSolver::solveP3(const VectorXd &r, VectorXd &duOut) {
  const int n = VDefault.rows();
  const VectorXd Y = VectorXd::Zero(fixedIdx.size());
  if (isinf(weight))
    return igl::min_quad_with_fixed_solve(data, VectorXd::Zero(n), Y, r, duOut);
  // 0.5 du'(Q + w B'B) du - w du' B' r
  const VectorXd lin = -weight * (B.transpose() * r);
  return igl::min_quad_with_fixed_solve(data, lin, Y, VectorXd(), duOut);
}

bool LandmarkSolver::solveBaseline(const vector<Landmark> &lms,
                                   const VectorXd &r, VectorXd &duOut) {
  // r(k) is the required change of u at landmark k. duOut arrives as zeros
  // with one entry per vertex and must stay 0 on the fixed vertices.
  // TODO(human) [D4 signs]: u0 and r are negative on the back sheet
  //   (lms[k].side == -1). Check that each method moves a back handle away
  //   from the drawing plane for factors > 1.
  switch (method) {
    case HandleMethod::G1:
      // TODO(human) [G1]: duOut(vertex of k) = r(k).
      break;
    case HandleMethod::G2:
      // TODO(human) [G2]: solve the k x k system for c, then
      // duOut = u0 .* sum_k c_k w_k.
      break;
    case HandleMethod::P1:
    case HandleMethod::P2:
      // TODO(human) [P1, P2]: alpha = (B g)^-1 r, duOut = sum_k alpha_k g_k.
      break;
    case HandleMethod::P3:
      return false;  // handled by solveP3()
  }
  return true;
}

VectorXd LandmarkSolver::forward(const VectorXd &dp) const {
  // TODO(human) [forward]: du_I = L_II^-1 M_II dp_I with the factorization
  // from prepare() (mind the sign: it factors -L_II), du = 0 on the fixed
  // vertices.
  return VectorXd::Zero(dp.size());
}

VectorXd LandmarkSolver::implied(const VectorXd &du) const {
  // TODO(human) [implied]: dp_I = (L du)_I / M_II, dp = 0 on the fixed
  // vertices. Used for every method to report the pressure change.
  return VectorXd::Zero(du.size());
}

bool LandmarkSolver::exportScene(const string &path,
                                 const vector<Landmark> &lms) const {
  if (!ready) return false;
  using nlohmann::json;

  json V = json::array();
  for (int i = 0; i < VDefault.rows(); i++)
    V.push_back({VDefault(i, 0), VDefault(i, 1), VDefault(i, 2)});
  json Fj = json::array();
  for (int i = 0; i < F.rows(); i++) Fj.push_back({F(i, 0), F(i, 1), F(i, 2)});

  json Li = json::array(), Lj = json::array(), Lv = json::array();
  for (int k = 0; k < L.outerSize(); k++)
    for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      Li.push_back(it.row());
      Lj.push_back(it.col());
      Lv.push_back(it.value());
    }
  const VectorXd Mdiag = M.diagonal();

  // read by sds/scene.py (Scene.landmark_matrix builds the same B)
  json landmarks = json::array();
  for (const Landmark &lm : lms)
    landmarks.push_back({{"vertex", lm.vertexId},
                         {"factor", lm.factor},
                         {"side", lm.side},
                         {"anchor", lm.anchor},
                         {"baseHeight", lm.baseHeight}});

  json scene = {{"version", 1},
                {"V", V},
                {"F", Fj},
                {"fixed", vector<int>(fixedIdx.data(),
                                      fixedIdx.data() + fixedIdx.size())},
                {"L", {{"i", Li}, {"j", Lj}, {"v", Lv}}},
                {"M", vector<double>(Mdiag.data(),
                                     Mdiag.data() + Mdiag.size())},
                {"landmarks", landmarks},
                {"method", handleMethodName(method)},
                {"radius", radius},
                {"du", vector<double>(du.data(), du.data() + du.size())},
                {"dp", vector<double>(dp.data(), dp.data() + dp.size())}};

  ofstream out(path);
  if (!out) return false;
  out << scene.dump();
  return out.good();
}
