# ============================================================================
# The Unpinning Game - Makefile
#
# Flat layout: all .cpp / .h live in this directory. Builds every .cpp, links
# against raylib. Vendored raylib is expected in include/ (raylib.h) and lib/
# (libraylib.a); a system raylib also works (the -I/-L just no-op then).
#
# Targets:  make            build game(.exe)
#           make clean       remove build/ and the binary (Windows + Unix)
#           make             builds a self-contained (static) game.exe by default
#           make STATIC=0    smaller dynamic build (needs MSYS2 DLLs on PATH)
# ============================================================================

CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++20 -O2 -I. -Iinclude
TARGET   = game
BUILDDIR = build

# Every source in this folder (Sandbox_UPG.cpp holds main()).
SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,$(BUILDDIR)/%.o,$(SRC))
HDR = $(wildcard *.h)

LDFLAGS = -Llib
LDLIBS  = -lraylib

# Self-contained linking is ON by default (a statically linked game.exe is
# immune to stray libstdc++/libwinpthread DLLs on PATH and runs on machines
# without MSYS2). Opt out with `make STATIC=0` for a smaller dynamic build.
STATIC ?= 1
ifeq ($(STATIC),1)
	LDFLAGS += -static -static-libgcc -static-libstdc++
endif

# ---- per-OS settings -------------------------------------------------------
ifeq ($(OS),Windows_NT)
	SYSLIBS = -lopengl32 -lgdi32 -lwinmm
	TARGET := $(TARGET).exe
	# cmd.exe-friendly clean (no 'rm' on Windows)
	RM_DIR  = if exist $(BUILDDIR) rmdir /s /q $(BUILDDIR)
	RM_BIN  = if exist $(TARGET) del /q $(TARGET)
	MKDIR   = if not exist $(BUILDDIR) mkdir $(BUILDDIR)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		SYSLIBS = -lGL -lm -lpthread -ldl -lrt -lX11
	endif
	ifeq ($(UNAME_S),Darwin)
		SYSLIBS = -framework CoreVideo -framework IOKit -framework Cocoa -framework OpenGL
	endif
	RM_DIR  = rm -rf $(BUILDDIR)
	RM_BIN  = rm -f game game.exe
	MKDIR   = mkdir -p $(BUILDDIR)
endif

# ---- rules -----------------------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS) $(SYSLIBS)

# Recompile a TU if any header changes (coarse but safe).
$(BUILDDIR)/%.o: %.cpp $(HDR) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR):
	$(MKDIR)

clean:
	$(RM_DIR)
	$(RM_BIN)

.PHONY: all clean
