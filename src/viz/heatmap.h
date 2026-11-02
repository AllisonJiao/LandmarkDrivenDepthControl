// Heat maps of the landmark edit (Part B3): per-vertex displacement and
// pressure change, mapped to colours with a diverging colormap.
//
// Sign convention (shared with the offline plots): "outward positive", i.e.
// thickening and more inflation are positive on both the front and the back
// sheet. Colours use the same table as scripts/ plotting (RdBu_r, 9 stops).

#ifndef HEATMAP_H
#define HEATMAP_H

#include <Eigen/Dense>
#include <vector>

enum class HeatField { Off = 0, Displacement = 1, Pressure = 2 };

struct HeatmapOptions {
  HeatField field = HeatField::Off;
  bool symLog = false;    // symmetric-log scaling for fields spanning decades
  double symLogF0 = 1.0;  // linear range of the symmetric log, in field units
  double lockedVmax = 0;  // > 0: fixed colour range, otherwise automatic
};

// Diverging colour table, blue (negative) -> white (0) -> red (positive),
// RGB in [0,1]. Keep in sync with the Python plotting script.
extern const double kHeatmapTable[9][3];

// Displacement of the rest surface relative to the default inflation,
// dh_i = |h_i| - |h0_i|, with h = VRest.col(2), h0 = VDefault.col(2).
Eigen::VectorXd displacementField(const Eigen::MatrixXd &VRest,
                                  const Eigen::MatrixXd &VDefault);

// Pressure change relative to the base inflation, dp_i = side_i * dpRaw_i / c
// with side_i = sgn(h0_i). `dpRaw` is the solver's dp (empty if not
// available yet); entries of `fixed` vertices are undefined and set to NaN.
Eigen::VectorXd pressureField(const Eigen::VectorXd &dpRaw,
                              const Eigen::MatrixXd &VDefault, double c,
                              const std::vector<bool> &fixed);

// Maps values to RGB rows of C (n x 3, in [0,1]). vmax <= 0 means automatic
// (max |f| over finite entries); NaN entries are drawn grey. The range that
// was actually used is returned in vmaxUsed.
void scalarsToColors(const Eigen::VectorXd &f, double vmax, bool symLog,
                     double symLogF0, Eigen::MatrixXd &C, double &vmaxUsed);

#endif  // HEATMAP_H
