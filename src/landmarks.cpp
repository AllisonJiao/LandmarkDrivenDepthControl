#include "landmarks.h"

#include <miscutils/mesh3d.h>

#include <cmath>

#include "def3dsdl.h"

using namespace std;
using namespace Eigen;

namespace {

Vector2d projectPoint(const Matrix4d &proj, const Vector3d &p) {
  return (proj * p.homogeneous()).hnormalized().head(2);
}

// Face under (x,y) in V together with its barycentric weights. The weights
// that getMeshFace can output belong to the last face it tested rather than
// the returned one, so they are recomputed here from the face's xy triangle.
int faceWithBary(const MatrixXd &V, const MatrixXi &F, double x, double y,
                 bool highestZ, bool reverse, Vector3d &bary) {
  const int face = Mesh3D::getMeshFace(V, F, x, y, highestZ, reverse);
  if (face == -1) return -1;
  const Vector2d a = V.row(F(face, 0)).head(2);
  const Vector2d b = V.row(F(face, 1)).head(2);
  const Vector2d c = V.row(F(face, 2)).head(2);
  const Vector2d p(x, y);
  auto cross = [](const Vector2d &u, const Vector2d &v) {
    return u(0) * v(1) - u(1) * v(0);
  };
  const double area = cross(b - a, c - a);
  if (area == 0) return -1;
  bary(1) = cross(p - a, c - a) / area;
  bary(2) = cross(b - a, p - a) / area;
  bary(0) = 1 - bary(1) - bary(2);
  return face;
}

Vector3d blend(const MatrixXd &V, const MatrixXi &F, int face,
               const Vector3d &bary) {
  Vector3d p = Vector3d::Zero();
  for (int k = 0; k < 3; k++) p += bary(k) * V.row(F(face, k)).transpose();
  return p;
}

}  // namespace

bool pickLandmark(const MatrixXd &VCurr, const MatrixXd &VRest,
                  const MatrixXi &F, const Matrix4d &proj,
                  const Vector2d &mouse, Landmark &out) {
  if (VCurr.rows() == 0 || VRest.rows() != VCurr.rows()) return false;
  const MatrixXd VProj =
      (VCurr.rowwise().homogeneous() * proj.transpose()).rowwise().hnormalized();
  Vector3d bary;
  // same face selection as a left click on a control point in Animate mode
  const int face =
      faceWithBary(VProj, F, mouse(0), mouse(1), false, false, bary);
  if (face == -1) return false;

  const Vector3d rest = blend(VRest, F, face, bary);
  out = Landmark();
  out.faceId = face;
  out.bary = bary;
  out.sketchPos = rest.head(2);
  out.side = rest(2) >= 0 ? 1 : -1;
  out.baseHeight = abs(rest(2));
  return true;
}

Vector3d surfacePoint(const MatrixXd &V, const MatrixXi &F,
                      const Landmark &lm) {
  return blend(V, F, lm.faceId, lm.bary);
}

Vector3d targetPoint(const MatrixXd &VDefault, const MatrixXi &F,
                     const Landmark &lm) {
  Vector3d p = surfacePoint(VDefault, F, lm);
  p(2) += lm.side * (lm.factor - 1) * lm.baseHeight;
  return p;
}

int hitTestLandmarks(const LandmarkData &data, const MatrixXd &V,
                     const MatrixXd &VDefault, const MatrixXi &F,
                     const Matrix4d &proj, const Vector2d &mouse,
                     double radius) {
  int best = -1;
  double bestDist = radius;
  for (int i = 0; i < data.items.size(); i++) {
    const Landmark &lm = data.items[i];
    for (const Vector3d &p :
         {targetPoint(VDefault, F, lm), surfacePoint(V, F, lm)}) {
      const double d = (projectPoint(proj, p) - mouse).norm();
      if (d <= bestDist) {
        bestDist = d;
        best = i;
      }
    }
  }
  return best;
}

double factorFromDrag(const Vector2d &screenAt1, const Vector2d &screenAt2,
                      const Vector2d &mouse, double currentFactor) {
  // project the mouse onto the screen line from screenAt1 (factor 1) towards
  // screenAt2 (factor 2)
  Vector2d d = screenAt2 - screenAt1;
  if (d.squaredNorm() == 0) return currentFactor;
  currentFactor = 1.0 + (mouse - screenAt1).dot(d) / d.squaredNorm();
  return currentFactor;
}

void dragLandmark(Landmark &lm, const MatrixXd &VDefault, const MatrixXi &F,
                  const Matrix4d &proj, const Vector2d &mouse) {
  const Vector3d base = surfacePoint(VDefault, F, lm);
  const Vector3d doubled = base + Vector3d(0, 0, lm.side * lm.baseHeight);
  lm.factor = factorFromDrag(projectPoint(proj, base),
                             projectPoint(proj, doubled), mouse, lm.factor);
}

void rebindLandmarks(LandmarkData &data, const MatrixXd &VRest,
                     const MatrixXi &F) {
  vector<Landmark> kept;
  for (Landmark lm : data.items) {
    // The rest mesh is viewed orthographically from the sketch plane; front
    // and back sheets wind in opposite directions, so try both and keep the
    // outermost face on the landmark's side.
    int bestFace = -1;
    Vector3d bestBary;
    double bestHeight = -1;
    for (bool reverse : {false, true}) {
      Vector3d bary;
      const int face = faceWithBary(VRest, F, lm.sketchPos(0), lm.sketchPos(1),
                                    lm.side > 0, reverse, bary);
      if (face == -1) continue;
      const double z = blend(VRest, F, face, bary)(2);
      if (z * lm.side < 0) continue;
      if (abs(z) > bestHeight) {
        bestHeight = abs(z);
        bestFace = face;
        bestBary = bary;
      }
    }
    if (bestFace == -1) continue;
    lm.faceId = bestFace;
    lm.bary = bestBary;
    lm.baseHeight = bestHeight;
    kept.push_back(lm);
  }
  data.items = kept;
  data.selected = -1;
  data.dragging = false;
}

void drawLandmarks(const LandmarkData &data, const MatrixXd &V,
                   const MatrixXd &VDefault, const MatrixXi &F,
                   const Matrix4d &proj,
                   MyPainter &painter) {
  const Cu colorBlack{0, 0, 0, 255};
  const Cu colorRed{255, 0, 0, 255};
  const Cu colorHandle{30, 110, 230, 255};
  const Cu colorAnchor{90, 90, 90, 255};
  for (int i = 0; i < data.items.size(); i++) {
    const Landmark &lm = data.items[i];
    const Vector3d base = surfacePoint(V, F, lm);
    const Vector3d target = targetPoint(VDefault, F, lm);
    const Vector2d b = projectPoint(proj, base);
    const Vector2d t = projectPoint(proj, target);
    const Cu &color = lm.anchor ? colorAnchor : colorHandle;

    painter.setColor(color);
    painter.drawLine(b(0), b(1), t(0), t(1), 2);
    painter.filledEllipse(b(0), b(1), 3, 3);
    if (i == data.selected) {
      drawControlPoint(VectorXd(target), painter, 10, proj, 3, colorRed,
                       color);
    } else {
      drawControlPoint(VectorXd(target), painter, 7, proj, 1, colorBlack,
                       color);
    }
  }
}
