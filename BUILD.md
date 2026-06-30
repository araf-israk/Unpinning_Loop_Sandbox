# The Unpinning Game - build & run

## Just build it
```
make
```
That's the whole story. The first `make` compiles **raylib from the bundled
source** (`raylib/`) with your own compiler, then builds and links the game into
a single standalone `game.exe` (on Windows) or `game` (Linux/macOS). There is
**nothing else to install** - no raylib package, no DLLs, no vendored binary that
can mismatch your toolchain.

```
make clean    # remove build/, lib/, and the binary
make          # rebuild (raylib is cached in lib/, so this is fast after the first run)
```

## Why it's built this way
Earlier attempts linked a system or prebuilt raylib, which broke when moved to
another machine (`cannot find -lraylib`, or `undefined reference to __imp_glfw*`
from a raylib `.a` that wanted a separate GLFW DLL). Building raylib from the
bundled source with the same compiler that builds the game removes all of that:
the ABI always matches, GLFW is compiled in, and on Windows the result is fully
static (`-static`), so stray `libstdc++ / libwinpthread` DLLs on PATH (e.g. from
an STM32 toolchain) can't break it and `game.exe` runs on any Windows machine.

## Requirements
- A C/C++ toolchain with **C++20** (MSYS2 UCRT64 `g++`/`gcc`, or any MinGW-w64).
- `make`, `ar` (both come with MSYS2/MinGW).
- Linux/macOS only: the usual OpenGL + windowing dev headers
  (e.g. `libgl1-mesa-dev libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev
  libxi-dev libasound2-dev` on Debian/Ubuntu). Windows needs none of these - the
  system already ships `opengl32`/`gdi32`/`winmm`.

## Options
```
make RAYLIB_SRC=0   # skip building bundled raylib; link a system/prebuilt -lraylib instead
```

## Layout
- `*.cpp` / `*.h`  - the game (entry point: `Sandbox_UPG.cpp`).
- `raylib/`        - bundled raylib 5.5 source (built into `lib/libraylib.a`).
- `build/`, `lib/` - generated; safe to delete (`make clean`).

## Verified
The full flow (build raylib from `raylib/` source -> link the game) was tested
end-to-end and produces a working executable; the menu accepts a vertex sigma
permutation and starts the game.
