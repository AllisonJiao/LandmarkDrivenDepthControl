// Landmark handles for depth control: a landmark sits on the inflated surface
// and stores a target height relative to the default inflation
// (factor 1.5 means 50% more protrusion at that point).

#ifndef LANDMARKS_H
#define LANDMARKS_H

#include <Eigen/Dense>
#include <vector>

#include "mypainter.h"

struct Landmark {
  Eigen::Vector2d sketchPos;  // x,y on the rest mesh, used to re-attach
  int side = 1;               // sign of rest z, keeps it on the same sheet
  int faceId = -1;
  Eigen::Vector3d bary = Eigen::Vector3d::Zero();  // weights of F.row(faceId)
  double baseHeight = 0;  // |z| of the rest surface at the landmark
  double factor = 1;      // target height = factor * baseHeight
  bool anchor = false;    // anchors hold the current depth and can't be dragged
};

struct LandmarkData {
  std::vector<Landmark> items;
  int selected = -1;
  bool dragging = false;
};

// Picks the front-most face under the screen position `mouse` (the model is
// projected with `proj`) and fills the landmark's attachment data. Heights
// are taken from the rest mesh `VRest`, which shares faces with `VCurr`.
bool pickLandmark(const Eigen::MatrixXd &VCurr, const Eigen::MatrixXd &VRest,
                  const Eigen::MatrixXi &F, const Eigen::Matrix4d &proj,
                  const Eigen::Vector2d &mouse, Landmark &out);

// Point on the (possibly deformed) surface V at the landmark.
Eigen::Vector3d surfacePoint(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                             const Landmark &lm);

// Point on the default (unedited) surface VDefault moved along z so that its
// height becomes factor * baseHeight.
Eigen::Vector3d targetPoint(const Eigen::MatrixXd &VDefault,
                            const Eigen::MatrixXi &F, const Landmark &lm);

// Index of the landmark whose target or surface point (on the current
// surface V) is within `radius` pixels of `mouse`, or -1.
int hitTestLandmarks(const LandmarkData &data, const Eigen::MatrixXd &V,
                     const Eigen::MatrixXd &VDefault, const Eigen::MatrixXi &F,
                     const Eigen::Matrix4d &proj, const Eigen::Vector2d &mouse,
                     double radius);

// Converts a mouse position into a new factor. `screenAt1` and `screenAt2` are
// the screen positions of the landmark's target at factor 1 and factor 2.
double factorFromDrag(const Eigen::Vector2d &screenAt1,
                      const Eigen::Vector2d &screenAt2,
                      const Eigen::Vector2d &mouse, double currentFactor);

// Updates the factor of landmark `lm` from a mouse position. Measured on the
// default surface so that edits of the current surface don't feed back.
void dragLandmark(Landmark &lm, const Eigen::MatrixXd &VDefault,
                  const Eigen::MatrixXi &F, const Eigen::Matrix4d &proj,
                  const Eigen::Vector2d &mouse);

// Re-attaches landmarks to a newly reconstructed rest mesh by their sketch
// position; landmarks that fall outside the new shape are removed.
void rebindLandmarks(LandmarkData &data, const Eigen::MatrixXd &VRest,
                     const Eigen::MatrixXi &F);

void drawLandmarks(const LandmarkData &data, const Eigen::MatrixXd &V,
                   const Eigen::MatrixXd &VDefault, const Eigen::MatrixXi &F, const Eigen::Matrix4d &proj,
                   MyPainter &painter);

#endif  // LANDMARKS_H
