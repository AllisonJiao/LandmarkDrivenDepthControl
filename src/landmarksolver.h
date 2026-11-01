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

  // Take effect at the next setConstraints. The falloff radius R of G2 and P2
  // is in mean edge lengths of the rest mesh, so one value gives the same
  // locality on a coarse or a fine mesh.
  void setMethod(HandleMethod method);
  void setFalloffRadius(double edges);
  HandleMethod getMethod() const { return method; }
  double getFalloffRadius() const { return radiusEdges; }
  // R in canvas pixels on the current mesh
  double falloffRadiusPixels() const { return radiusEdges * meanEdge; }

  // Writes the default mesh, the flat operators L and M, the outline
  // indices, the landmarks and the fields of the last solve as JSON for the
  // external tools in sds/ (the M1 contract's export; sds/scene.py reads it).
  // inflationAmount is c, the default pressure that dp is reported against.
  bool exportScene(const std::string &path, const std::vector<Landmark> &lms,
                   const std::string &meshName, double inflationAmount) const;

  bool isReady() const { return ready; }
  const Eigen::MatrixXd &defaultMesh() const { return VDefault; }
  // vertices held at u = 0: the outlines of all parts
  const Eigen::VectorXi &fixedVertices() const { return fixedIdx; }
  // the same set as one flag per vertex
  const std::vector<bool> &fixedMask() const { return isFixed; }
  // Fields of the last solve, per vertex, in the solver's signs: du = change
  // of u = sgn(h) h^2, dp = implied change of pressure (0 on the outline).
  // Zero until the first solve after prepare().
  const Eigen::VectorXd &lastDu() const { return du; }
  const Eigen::VectorXd &lastDp() const { return dp; }
  // The same as the M1 contract reports them, "outward positive": thickening
  // is positive on both sheets. du * sgn(u0), dp * -sgn(u0), and
  // dh = sgn(h0) (h - h0), so a sheet crossing shows as a large negative dh.
  const Eigen::VectorXd &lastDuOutward() const { return duOut; }
  const Eigen::VectorXd &lastDpOutward() const { return dpOut; }
  const Eigen::VectorXd &lastDhOutward() const { return dhOut; }
  // false until the first solve after prepare() and after a failed solve
  // (rejected landmark, singular setup); the fields are 0 then
  bool lastSolveOk() const { return solveOk; }

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
  // (1 - (d / R)^2)^3 on the vertices of v's sheet within R of v (x,y
  // distance d, R in pixels), 0 elsewhere; 1 at v. The bump of G2 and P2.
  Eigen::VectorXd falloff(int v) const;

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
  double radiusEdges = 4;  // falloff radius R in mean edge lengths
  double meanEdge = 0;     // mean x,y edge length of the rest mesh (prepare)
  Eigen::VectorXd du, dp;  // fields of the last solve
  Eigen::VectorXd duOut, dpOut, dhOut;  // the same, outward positive
  Eigen::VectorXd sheetSign;            // sgn(u0) = sgn(h0), 0 on the outline
  bool solveOk = false;
  std::vector<bool> isFixed;  // per vertex: in fixedIdx

  // once per mesh (prepare): the interior vertices, their inverse map (-1 on
  // the fixed vertices), the interior mass, and -L_II factored (it is SPD)
  std::vector<int> interiorIdx;
  Eigen::VectorXi toInterior;
  Eigen::VectorXd massI;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  // per vertex: its sheet, a connected component of the interior (each part's
  // front and back are separate sheets), -1 on the fixed vertices
  Eigen::VectorXi sheetLabels;

  // once per landmark set (setConstraintsBaseline), read on every drag frame:
  // one column of du per handle (n x k; G2: u0 .* falloff(v), P1: forward(e_v),
  // P2: forward(falloff(v))), and the LU of A = B * basis (k x k), so that
  // du = basis * A^-1 r hits every target
  Eigen::MatrixXd basis;
  Eigen::FullPivLU<Eigen::MatrixXd> basisLU;
};

#endif  // LANDMARKSOLVER_H
