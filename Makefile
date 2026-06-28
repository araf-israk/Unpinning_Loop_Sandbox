# Compiler and application settings
CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++20 -O2 -Iinclude
TARGET   = game
SRC      = src/main.cpp

# Base Linker Flags (look inside the local lib folder)
LDFLAGS  = -Llib -lraylib

# OS Detection
ifeq ($(OS),Windows_NT)
	# Windows dependencies
	LDFLAGS += -lopengl32 -lgdi32 -lwinmm
	TARGET  := $(TARGET).exe
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S), Linux)
		# Linux dependencies
		LDFLAGS += -lGL -lm -lpthread -ldl -lrt -lX11
	endif
	ifeq ($(UNAME_S), Darwin)
		# macOS dependencies & Frameworks
		LDFLAGS += -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL
	endif
endif

# Build Rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)