// Headless landmark solve for the offline scripts (scripts/cpp_bridge.py):
// runs the app's LandmarkSolver on a scene file and writes one result per
// method in the format of LandmarkSolver::exportScene (scene version 2), so
// Python plots exactly what the app computes.
//
//   node landmark_cli.js --scene in.json --out-dir DIR
//        [--methods G1,G2,P1,P2,P3] [--radius 4] [--weight W] [--prefix NAME]
//
// The input needs V, F, fixed and landmarks [{vertex, factor, anchor}] (a
// base scene plus resolved handles; an app export works as is), and
// optionally c and mesh_name. Output: DIR/NAME_<method>.json for each method.
// The exit code is 0 if every method solved, 1 if any failed (its file is
// still written, with "ok": false), 2 on bad input.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "json.hpp"
#include "landmarksolver.h"

using namespace std;
using namespace Eigen;
using nlohmann::json;

namespace {

struct Args {
  string scene, outDir = ".", prefix = "result";
  vector<string> methods{"G1", "G2", "P1", "P2", "P3"};
  double radius = 4;
  double weight = numeric_limits<double>::infinity();  // hard
};

vector<string> splitList(const string &s) {
  vector<string> out;
  stringstream ss(s);
  string item;
  while (getline(ss, item, ','))
    if (!item.empty()) out.push_back(item);
  return out;
}

bool parseArgs(int argc, char *argv[], Args &args) {
  for (int i = 1; i < argc; i++) {
    const string key = argv[i];
    if (i + 1 >= argc) {
      cerr << "missing value for " << key << endl;
      return false;
    }
    const string value = argv[++i];
    if (key == "--scene") args.scene = value;
    else if (key == "--out-dir") args.outDir = value;
    else if (key == "--prefix") args.prefix = value;
    else if (key == "--methods") args.methods = splitList(value);
    else if (key == "--radius") args.radius = stod(value);
    else if (key == "--weight") args.weight = stod(value);
    else {
      cerr << "unknown option " << key << endl;
      return false;
    }
  }
  if (args.scene.empty()) cerr << "--scene is required" << endl;
  return !args.scene.empty();
}

bool methodFromName(const string &name, HandleMethod &method) {
  for (HandleMethod m : {HandleMethod::G1, HandleMethod::G2, HandleMethod::P1,
                         HandleMethod::P2, HandleMethod::P3})
    if (name == handleMethodName(m)) {
      method = m;
      return true;
    }
  return false;
}

// the bundled nlohmann json predates contains()
bool has(const json &j, const char *key) { return j.find(key) != j.end(); }

// parse without exceptions: the Emscripten build doesn't catch them
bool loadScene(const string &path, MatrixXd &V, MatrixXi &F, vector<int> &fixed,
               json &landmarks, double &c, string &meshName) {
  ifstream in(path);
  if (!in) {
    cerr << "cannot open " << path << endl;
    return false;
  }
  const json d = json::parse(in, nullptr, false);
  if (d.is_discarded() || !d.is_object() || !has(d, "V") || !has(d, "F") ||
      !has(d, "fixed")) {
    cerr << path << ": not a scene (needs V, F, fixed)" << endl;
    return false;
  }
  const json &Vj = d["V"], &Fj = d["F"];
  V.resize(Vj.size(), 3);
  for (int i = 0; i < V.rows(); i++)
    for (int k = 0; k < 3; k++) V(i, k) = Vj[i][k].get<double>();
  F.resize(Fj.size(), 3);
  for (int i = 0; i < F.rows(); i++)
    for (int k = 0; k < 3; k++) F(i, k) = Fj[i][k].get<int>();
  fixed = d["fixed"].get<vector<int>>();
  landmarks = d.value("landmarks", json::array());
  c = d.value("c", 2.0);
  meshName = d.value("mesh_name", string());
  return true;
}

// the same attachment data pickLandmark fills, from a vertex id
bool makeLandmarks(const json &lj, const MatrixXd &V, vector<Landmark> &lms) {
  for (const json &l : lj) {
    if (!has(l, "vertex") || !has(l, "factor")) {
      cerr << "landmark needs vertex and factor" << endl;
      return false;
    }
    Landmark lm;
    lm.vertexId = l["vertex"].get<int>();
    if (lm.vertexId < 0 || lm.vertexId >= V.rows()) {
      cerr << "landmark vertex " << lm.vertexId << " out of range" << endl;
      return false;
    }
    lm.factor = l["factor"].get<double>();
    lm.anchor = l.value("anchor", false);
    const double z = V(lm.vertexId, 2);
    lm.sketchPos = V.row(lm.vertexId).head<2>().transpose();
    lm.side = z >= 0 ? 1 : -1;
    lm.baseHeight = abs(z);
    lms.push_back(lm);
  }
  return true;
}

}  // namespace

int main(int argc, char *argv[]) {
  Args args;
  if (!parseArgs(argc, argv, args)) return 2;

  MatrixXd V;
  MatrixXi F;
  vector<int> fixed;
  json landmarksJson;
  double c;
  string meshName;
  if (!loadScene(args.scene, V, F, fixed, landmarksJson, c, meshName)) return 2;
  vector<Landmark> lms;
  if (!makeLandmarks(landmarksJson, V, lms)) return 2;

  // prepare() recomputes L and a Voronoi M from V and ignores the scene's own
  // L and M. scripts/make_grid_scenes.py solves the synthetic ellipse's
  // heights with that same Voronoi M, so there p0 = M^-1 L u0 is exactly +-c.
  std::error_code ec;
  filesystem::create_directories(args.outDir, ec);  // NODERAWFS: the host dir

  LandmarkSolver solver;
  solver.prepare(V, F, {fixed});
  if (!solver.isReady()) {
    cerr << "solver setup failed (factorization of -L_II)" << endl;
    return 2;
  }

  bool allOk = true;
  for (const string &name : args.methods) {
    HandleMethod method;
    if (!methodFromName(name, method)) {
      cerr << "unknown method " << name << endl;
      return 2;
    }
    solver.setMethod(method);
    solver.setFalloffRadius(args.radius);
    MatrixXd VRest;
    // solve() runs even after a failed setup: it resets the fields and
    // lastSolveOk(), so a failed method never exports the previous method's
    // fields. No landmarks is a valid run (setConstraints refuses it, solve
    // returns the default inflation).
    const bool constrained = solver.setConstraints(lms, args.weight);
    const bool solved = solver.solve(lms, VRest);
    const bool ok = (constrained || lms.empty()) && solved && solver.lastSolveOk();
    const string out = args.outDir + "/" + args.prefix + "_" + name + ".json";
    if (!solver.exportScene(out, lms, meshName, c)) {
      cerr << "cannot write " << out << endl;
      return 2;
    }
    cout << name << ": " << (ok ? "ok" : "FAILED") << " -> " << out << endl;
    allOk &= ok;
  }
  return allOk ? 0 : 1;
}
