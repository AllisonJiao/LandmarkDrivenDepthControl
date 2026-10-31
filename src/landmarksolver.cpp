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

void LandmarkSolver::prepare(const MatrixXd &VDefault, const MatrixXi &F,
                             const vector<vector<int>> &bnds) {
  this->VDefault = VDefault;
  this->F = F;
  ready = factored = false;
  numConstraints = 0;
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

  ready = true;
}

bool LandmarkSolver::setConstraints(const vector<Landmark> &lms,
                                    double weight) {
  factored = false;
  this->weight = weight;
  numConstraints = lms.size();
  if (!ready || lms.empty()) return false;

  const int n = VDefault.rows();
  vector<Triplet<double>> triplets;
  for (int k = 0; k < lms.size(); k++)
    for (int j = 0; j < 3; j++)
      triplets.emplace_back(k, F(lms[k].faceId, j), lms[k].bary(j));
  B = SparseMatrix<double>(lms.size(), n);
  B.setFromTriplets(triplets.begin(), triplets.end());

  if (isinf(weight)) {
    factored = igl::min_quad_with_fixed_precompute(Q, fixedIdx, B, false, data);
  } else {
    SparseMatrix<double> A = Q + weight * SparseMatrix<double>(B.transpose()) * B;
    factored = igl::min_quad_with_fixed_precompute(
        A, fixedIdx, SparseMatrix<double>(), true, data);
  }
  if (!factored) cerr << "landmark solve: factorization failed" << endl;
  return factored;
}

bool LandmarkSolver::solve(const vector<Landmark> &lms, MatrixXd &VRestOut) {
  VRestOut = VDefault;
  if (!ready || lms.empty()) return true;  // default inflation
  if (!factored || lms.size() != numConstraints) return false;

  const auto tStart = chrono::high_resolution_clock::now();

  // Targets in u-space. Scaling a height by f scales u by f^2, and taking the
  // default u at the landmark keeps factor 1 an exact no-op.
  const VectorXd uAtLandmarks = B * u0;
  VectorXd r(lms.size());
  for (int k = 0; k < lms.size(); k++)
    r(k) = (lms[k].factor * lms[k].factor - 1) * uAtLandmarks(k);

  const int n = VDefault.rows();
  const VectorXd Y = VectorXd::Zero(fixedIdx.size());
  VectorXd du;
  bool ok;
  if (isinf(weight)) {
    ok = igl::min_quad_with_fixed_solve(data, VectorXd::Zero(n), Y, r, du);
  } else {
    // 0.5 du'(Q + w B'B) du - w du' B' r
    const VectorXd lin = -weight * (B.transpose() * r);
    ok = igl::min_quad_with_fixed_solve(data, lin, Y, VectorXd(), du);
  }
  if (!ok) return false;

  const VectorXd u = u0 + du;
  VRestOut.col(2) = u.array().sign() * u.array().abs().sqrt();

  const double ms = chrono::duration<double, milli>(
                        chrono::high_resolution_clock::now() - tStart)
                        .count();
  cout << "landmark solve: " << ms << " ms (" << n << " vertices, "
       << lms.size() << " landmarks)" << endl;
  return true;
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

  json landmarks = json::array();
  for (const Landmark &lm : lms)
    landmarks.push_back({{"face", lm.faceId},
                         {"bary", {lm.bary(0), lm.bary(1), lm.bary(2)}},
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
                {"landmarks", landmarks}};

  ofstream out(path);
  if (!out) return false;
  out << scene.dump();
  return out.good();
}
