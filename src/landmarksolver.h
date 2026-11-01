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
#include <string>
#include <vector>

#include "landmarks.h"

// Quadratic form Q of the correction du (energy 0.5 * du' Q du), built from
// the flat cotangent Laplacian L (negative semi-definite) and the diagonal
// mass matrix M (Minv is its inverse).
Eigen::SparseMatrix<double> buildEnergy(const Eigen::SparseMatrix<double> &L,
                                        const Eigen::SparseMatrix<double> &M,
                                        const Eigen::SparseMatrix<double> &Minv);

// How landmark targets become a new rest shape. Every method hits the target
// u at each handle exactly (P3 in hard mode); they differ in which field
// they edit and where:
//   G1  du at the handle vertex only
//   G2  du with a radial falloff of radius R
//   P1  dp at the handle vertex only, du = forward(dp)
//   P2  dp with a radial falloff of radius R, du = forward(dp)
//   P3  dp everywhere, minimizing the energy of buildEnergy (KKT solve)
// The values match the method dropdown in the web UI.
enum class HandleMethod { G1 = 0, G2 = 1, P1 = 2, P2 = 3, P3 = 4 };

const char *handleMethodName(HandleMethod method);

class LandmarkSolver {
 public:
  // Called after every reconstruction with the default (unedited) rest mesh.
  void prepare(const Eigen::MatrixXd &VDefault, const Eigen::MatrixXi &F,
               const std::vector<std::vector<int>> &bnds);
  // Rebuilds the constraint matrix and factors the system. Call when
  // landmarks are added, removed or toggled, or when the weight, the method
  // or the radius changes. weight = infinity means hard constraints.
  bool setConstraints(const std::vector<Landmark> &lms, double weight);
  // Solves for the current landmark factors and writes the new rest heights
  // into VRestOut (x,y are copied from the default mesh).
  bool solve(const std::vector<Landmark> &lms, Eigen::MatrixXd &VRestOut);

  // Takes effect at the next setConstraints. radius is in the units of the
  // rest mesh's x,y (canvas pixels) and is used by G2 and P2.
  void setMethod(HandleMethod method, double radius);
  HandleMethod getMethod() const { return method; }
  double getRadius() const { return radius; }

  // Writes the default mesh, the flat operators L and M, the outline
  // indices, the landmarks and the fields of the last solve as JSON for the
  // external tools in sds/.
  bool exportScene(const std::string &path,
                   const std::vector<Landmark> &lms) const;

  bool isReady() const { return ready; }
  const Eigen::MatrixXd &defaultMesh() const { return VDefault; }
  // vertices held at u = 0: the outlines of all parts
  const Eigen::VectorXi &fixedVertices() const { return fixedIdx; }
  // the same set as one flag per vertex
  const std::vector<bool> &fixedMask() const { return isFixed; }
  // Fields of the last solve, per vertex: du = change of u = sgn(h) h^2,
  // dp = implied change of pressure (0 on the outline). Empty before a solve.
  const Eigen::VectorXd &lastDu() const { return du; }
  const Eigen::VectorXd &lastDp() const { return dp; }

 private:
  // P3: factors the KKT system of buildEnergy with the landmark rows B
  bool setConstraintsP3();
  bool solveP3(const Eigen::VectorXd &r, Eigen::VectorXd &duOut);
  // G1, G2, P1, P2: per-landmark-set precomputation and per-frame solve
  bool setConstraintsBaseline(const std::vector<Landmark> &lms);
  bool solveBaseline(const std::vector<Landmark> &lms, const Eigen::VectorXd &r,
                     Eigen::VectorXd &duOut);
  // factors A = B * basis; false if it is singular
  bool factorBasis();

  // Pressure -> displacement: du = L_II^-1 M_II dp on the interior vertices,
  // 0 on the fixed ones. dp and the result have one entry per vertex.
  Eigen::VectorXd forward(const Eigen::VectorXd &dp) const;
  // Displacement -> pressure: dp = M_II^-1 L_I du on the interior vertices,
  // 0 on the fixed ones.
  Eigen::VectorXd implied(const Eigen::VectorXd &du) const;

  Eigen::MatrixXd VDefault;
  Eigen::MatrixXi F;
  Eigen::SparseMatrix<double> L, M;  // flat cotangent Laplacian, Voronoi mass
  Eigen::SparseMatrix<double> Q, B;
  Eigen::VectorXd u0;
  Eigen::VectorXi fixedIdx;
  igl::min_quad_with_fixed_data<double> data;
  double weight = 0;
  int numConstraints = 0;
  bool ready = false;
  bool factored = false;

  HandleMethod method = HandleMethod::P3;
  double radius = 50;
  Eigen::VectorXd du, dp;  // fields of the last solve
  std::vector<bool> isFixed;  // per vertex: in fixedIdx

  // once per mesh (prepare): the interior vertices, their inverse map (-1 on
  // the fixed vertices), the interior mass, and -L_II factored (it is SPD)
  std::vector<int> interiorIdx;
  Eigen::VectorXi toInterior;
  Eigen::VectorXd massI;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;

  // once per landmark set (setConstraintsBaseline), read on every drag frame:
  // one column of du per handle (n x k; P1: forward(e_v)), and the LU of
  // A = B * basis (k x k), so that du = basis * A^-1 r hits every target
  Eigen::MatrixXd basis;
  Eigen::FullPivLU<Eigen::MatrixXd> basisLU;
};

#endif  // LANDMARKSOLVER_H
