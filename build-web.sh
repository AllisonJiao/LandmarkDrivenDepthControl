#!/usr/bin/env bash
# Build the Monster Mash web app with Emscripten (Windows Git Bash / Linux / macOS).
# Expects emsdk and a cmake+ninja venv in ../tools (see LOCAL_SETUP.md); override with EMSDK / BUILD_TOOLS.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
EMSDK="${EMSDK:-$ROOT/../tools/emsdk}"
BUILD_TOOLS="${BUILD_TOOLS:-$ROOT/../tools/buildenv}"
NODE_DIR="$(ls -d "$EMSDK"/node/*/ | head -1)"
export EMSDK
export PATH="$BUILD_TOOLS/Scripts:$BUILD_TOOLS/bin:$EMSDK:$EMSDK/upstream/emscripten:$NODE_DIR:$NODE_DIR/bin:$PATH"

BUILD_TYPE="${1:-Release}"
mkdir -p "$ROOT/build/web"
cmake -S "$ROOT/src" -B "$ROOT/build/web" -G Ninja \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DCMAKE_TOOLCHAIN_FILE="$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake"
cmake --build "$ROOT/build/web"
echo "Built. Serve with: python -m http.server 8080 --directory \"$ROOT/build/web\""
