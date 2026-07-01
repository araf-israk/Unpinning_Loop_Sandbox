# ============================================================================
# The Unpinning Game - self-contained Makefile
#
# Download the project, run `make`, get a working game(.exe). raylib is BUILT
# FROM THE BUNDLED SOURCE in raylib/ with your own compiler, then linked - no
# raylib install, no DLL, no ABI mismatch.
#
#   make            build raylib (once) + the game
#   make clean      remove build/, lib/, and the binary
#   make RAYLIB_SRC=0   link a system/prebuilt -lraylib instead of building it
#
# Portable across cmd.exe, PowerShell, MSYS2, and w64devkit. Directory creation
# uses a flagless `mkdir` (valid for both Windows and Unix mkdir) run once via
# real directory targets, so it never depends on which shell make happens to use.
# ============================================================================

CXX = g++
CC  = gcc

TARGET   = game
BUILDDIR = build

# ---- our sources -----------------------------------------------------------
SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,$(BUILDDIR)/%.o,$(SRC))
HDR = $(wildcard *.h)

# raylib.h lives in the bundled raylib/ folder; -Iraylib lets <raylib.h> resolve.
CXXFLAGS = -Wall -Wextra -std=c++20 -O2 -I. -Iraylib

# ---- bundled raylib (built from source) ------------------------------------
RAYLIB_DIR     = raylib
RAYLIB_A       = lib/libraylib.a
RAYLIB_MODULES = rcore rshapes rtextures rtext rmodels raudio utils rglfw
# Flat object names (build/ray_rcore.o) so only single-level dirs are needed.
RAYLIB_OBJ     = $(patsubst %,$(BUILDDIR)/ray_%.o,$(RAYLIB_MODULES))
# Desktop GLFW backend, OpenGL 3.3. -w silences raylib/GLFW's own warnings.
RAYLIB_CFLAGS  = -O2 -std=gnu99 -DPLATFORM_DESKTOP -DGRAPHICS_API_OPENGL_33 -w \
                 -I$(RAYLIB_DIR) -I$(RAYLIB_DIR)/external/glfw/include

RAYLIB_SRC ?= 1          # 1 = build raylib from source (default); 0 = system -lraylib

LDFLAGS = -Llib
LDLIBS  = -lraylib

# ---- per-OS settings (NOT shell-dependent) ---------------------------------
ifeq ($(OS),Windows_NT)
	TARGET := game.exe
	SYSLIBS = -lopengl32 -lgdi32 -lwinmm
	RAYLIB_PLATDEF =
	# Fully static on Windows -> standalone game.exe, no DLLs, immune to stray
	# libstdc++/libwinpthread on PATH.
	LDFLAGS += -static -static-libgcc -static-libstdc++
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		SYSLIBS = -lGL -lm -lpthread -ldl -lrt -lX11 -lXrandr -lXinerama -lXcursor -lXi
		RAYLIB_PLATDEF = -D_GLFW_X11
	endif
	ifeq ($(UNAME_S),Darwin)
		SYSLIBS = -framework CoreVideo -framework IOKit -framework Cocoa -framework OpenGL
		RAYLIB_PLATDEF = -D_GLFW_COCOA
	endif
endif

ifeq ($(RAYLIB_SRC),0)
	RAYLIB_A =
endif

# ---- rules -----------------------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ) $(RAYLIB_A)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS) $(SYSLIBS)

# our C++ -> object (rebuild a TU if any of our headers change)
$(BUILDDIR)/%.o: %.cpp $(HDR) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# raylib C module -> object (the -I flags are for the desktop backend's GLFW header)
$(BUILDDIR)/ray_%.o: $(RAYLIB_DIR)/%.c | $(BUILDDIR)
	$(CC) $(RAYLIB_CFLAGS) $(RAYLIB_PLATDEF) -c $< -o $@

# bundle raylib objects into a static library
$(RAYLIB_A): $(RAYLIB_OBJ) | lib
	ar rcs $@ $(RAYLIB_OBJ)

# Directory targets: the recipe runs only when the directory does not yet exist,
# so `mkdir` is invoked exactly once (no repeated "already exists"). Flagless
# `mkdir DIR` is valid for both Windows mkdir and Unix mkdir; the leading '-'
# ignores a harmless failure if it already exists. NOTE: no `mkdir -p` (Windows
# mkdir would treat -p as a folder name) and no cmd `if not exist` (breaks in sh).
$(BUILDDIR):
	-mkdir $(BUILDDIR)
lib:
	-mkdir lib

# Clean tries both the Unix (rm) and Windows (rmdir/del) forms and ignores
# whichever your shell doesn't provide; a few "not found" lines here are harmless.
clean:
	-rm -rf $(BUILDDIR) lib game game.exe
	-rmdir /s /q $(BUILDDIR)
	-rmdir /s /q lib
	-del /q game.exe

.PHONY: all clean
