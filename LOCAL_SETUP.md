# Local setup (Windows)

This checkout builds and runs the full Monster Mash web app (the same UI as monstermash.zone) locally.

## Layout

```
workspace/
  monster-mash/          this repo
  tools/emsdk/           Emscripten SDK (6.0.13)
  tools/buildenv/        Python venv with cmake + ninja (pip install cmake ninja)
```

## Rebuild

From Git Bash:

```bash
./build-web.sh
```

Output goes to `build/web/` (`index.html`, `monstermash.{js,wasm,data}`, `includes/`, `imgs/`).

## Run

```bash
python -m http.server 8080 --directory build/web
```

Then open http://localhost:8080.

## What was added or changed compared to upstream

- `src/CMakeLists.txt`
  - Detect Emscripten with `if (EMSCRIPTEN)`. Upstream matched the compiler name `em++`, which fails on Windows (`em++.bat`).
  - `-sSTACK_SIZE=5MB`. Emscripten >= 3.1.27 shrank the default stack to 64KB, which crashes the app with "memory access out of bounds" as soon as a project is inflated.
  - Use `string(TIMESTAMP)` for the version instead of the Unix `date` command.
  - Preload `data/examples` into `/tmp/examples` when the folder exists, so the Help-page examples work.
  - Post-build step that copies the UI files (`includes/`, `imgs/`) next to the output and writes `index.html`.
- `third_party/triangle/`: Triangle library from netlib (gitignored because of its license).
- `data/examples/*.zip`: the 8 example projects, taken from the live monstermash.zone data bundle.
- `web-assets/imgs/`: UI icons, help videos and example thumbnails, taken from monstermash.zone (they aren't in the repo).
