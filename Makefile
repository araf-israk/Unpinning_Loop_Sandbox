# ============================================================================
# The Unpinning Game - self-contained Makefile
#
# Goal: download the project, run `make`, get a working game(.exe) with NOTHING
# else to install. raylib is BUILT FROM THE BUNDLED SOURCE in raylib/ using your
# own compiler, then linked statically - so there is no raylib install step, no
# DLL to ship, no ABI mismatch, and the resulting binary is standalone.
#
#   make            build raylib (once) + the game -> standalone game(.exe)
#   make clean      remove build/ , lib/ , and the binary
#   make RAYLIB_SRC=0   link a system/prebuilt -lraylib instead of building it
#
# Requirements: a C/C++ toolchain (MSYS2 UCRT64 g++/gcc, or any MinGW-w64; on
# Linux/macOS the system gcc plus the usual OpenGL dev headers). C++20.
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
RAYLIB_OBJ     = $(patsubst %,$(BUILDDIR)/raylib/%.o,$(RAYLIB_MODULES))
# Desktop GLFW backend, OpenGL 3.3. -w silences raylib/GLFW's own warnings.
RAYLIB_CFLAGS  = -O2 -std=gnu99 -DPLATFORM_DESKTOP -DGRAPHICS_API_OPENGL_33 -w \
                 -I$(RAYLIB_DIR) -I$(RAYLIB_DIR)/external/glfw/include

RAYLIB_SRC ?= 1          # 1 = build raylib from source (default); 0 = system -lraylib

LDFLAGS = -Llib
LDLIBS  = -lraylib

# ---- per-OS settings -------------------------------------------------------
ifeq ($(OS),Windows_NT)
	SYSLIBS = -lopengl32 -lgdi32 -lwinmm
	RAYLIB_PLATDEF =
	# Fully static on Windows: bundles libstdc++/libgcc/libwinpthread + our raylib
	# into game.exe, so no DLLs are needed and stray ones on PATH can't break it.
	LDFLAGS += -static -static-libgcc -static-libstdc++
	TARGET := $(TARGET).exe
	RM_ALL  = if exist $(BUILDDIR) rmdir /s /q $(BUILDDIR) & if exist lib rmdir /s /q lib & if exist $(TARGET) del /q $(TARGET)
	MKBUILD = if not exist $(BUILDDIR)\raylib mkdir $(BUILDDIR)\raylib
	MKLIB   = if not exist lib mkdir lib
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
	RM_ALL  = rm -rf $(BUILDDIR) lib game game.exe
	MKBUILD = mkdir -p $(BUILDDIR)/raylib
	MKLIB   = mkdir -p lib
endif

# If RAYLIB_SRC=0, don't build the bundled lib; rely on system/prebuilt -lraylib.
ifeq ($(RAYLIB_SRC),0)
	RAYLIB_A =
endif

# ---- rules -----------------------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ) $(RAYLIB_A)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS) $(SYSLIBS)

# our C++ -> object (rebuild a TU if any of our headers change)
$(BUILDDIR)/%.o: %.cpp $(HDR) | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

# raylib C module -> object (the -I flags are for the desktop backend's GLFW header)
$(BUILDDIR)/raylib/%.o: $(RAYLIB_DIR)/%.c | dirs
	$(CC) $(RAYLIB_CFLAGS) $(RAYLIB_PLATDEF) -c $< -o $@

# bundle raylib objects into a static library
$(RAYLIB_A): $(RAYLIB_OBJ) | dirs
	ar rcs $@ $(RAYLIB_OBJ)

dirs:
	$(MKLIB)
	$(MKBUILD)

clean:
	$(RM_ALL)

.PHONY: all clean dirs
