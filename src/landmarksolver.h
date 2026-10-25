// Full-dimensional landmark solve: finds the per-vertex change of inflation
// pressure that moves the surface to the landmark targets while staying
// closest to the default inflation (biharmonic constrained interpolation,
// cf. Zhang et al. 2001, Botsch & Kobbelt 2004, Jacobson et al. 2010).
//
// Works on the pre-sqrt field u = sgn(h) * h^2 of the rest mesh, where the
// inflation is linear: L u = M p on the flat mesh with u = 0 on the outline.

#ifndef LANDMARKSOLVER_H
#define LANDMARKSOLVER_H

// clang-format off
#include <Eigen/Core>
#include <Eigen/src/SparseCore/SparseSolverBase-mod.h>
#include <Eigen/Sparse>
#include <Eigen/src/OrderingMethods/Amd-mod.h>
#include <Eigen/src/OrderingMethods/Ordering-mod.h>
#include <Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include <Eigen/src/SparseCholesky/SimplicialCholesky_impl-mod.h>
#include <igl/min_quad_with_fixed.h>
// clang-format on

#include <Eigen/Dense>
#include <vector>

#include "landmarks.h"

// Quadratic form Q of the correction du (energy 0.5 * du' Q du), built from
// the flat cotangent Laplacian L (negative semi-definite) and the diagonal
// mass matrix M (Minv is its inverse).
Eigen::SparseMatrix<double> buildEnergy(const Eigen::SparseMatrix<double> &L,
                                        const Eigen::SparseMatrix<double> &M,
                                        const Eigen::SparseMatrix<double> &Minv);

class LandmarkSolver {
 public:
  // Called after every reconstruction with the default (unedited) rest mesh.
  void prepare(const Eigen::MatrixXd &VDefault, const Eigen::MatrixXi &F,
               const std::vector<std::vector<int>> &bnds);
  // Rebuilds the constraint matrix and factors the system. Call when
  // landmarks are added, removed or toggled, or when the weight changes.
  // weight = infinity means hard constraints.
  bool setConstraints(const std::vector<Landmark> &lms, double weight);
  // Solves for the current landmark factors and writes the new rest heights
  // into VRestOut (x,y are copied from the default mesh).
  bool solve(const std::vector<Landmark> &lms, Eigen::MatrixXd &VRestOut);

  bool isReady() const { return ready; }
  const Eigen::MatrixXd &defaultMesh() const { return VDefault; }

 private:
  Eigen::MatrixXd VDefault;
  Eigen::MatrixXi F;
  Eigen::SparseMatrix<double> Q, B;
  Eigen::VectorXd u0;
  Eigen::VectorXi fixedIdx;
  igl::min_quad_with_fixed_data<double> data;
  double weight = 0;
  int numConstraints = 0;
  bool ready = false;
  bool factored = false;
};

#endif  // LANDMARKSOLVER_H
