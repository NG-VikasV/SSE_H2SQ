# ===============================
# Compiler
# ===============================
CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -pedantic-errors -g \
            -msse2 -march=native -funroll-loops -ffast-math \
            -fomit-frame-pointer -fstrict-aliasing

# ===============================
# Source directory (unchanged)
# ===============================
SDIR := ../src
SRCS := $(wildcard $(SDIR)/*.cpp)

# ===============================
# Object files in PWD
# ===============================
OBJS := $(notdir $(SRCS:.cpp=.o))

TARGET := main

# ===============================
# Default target
# ===============================
all: $(TARGET)
	@echo Compilation Successful

# ===============================
# Link
# ===============================
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# ===============================
# Compile EACH source separately
# ===============================
%.o: $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ===============================
# Clean (PWD only)
# ===============================
clean:
ifeq ($(OS),Windows_NT)
	@cmd /c "del /Q *.o 2>NUL"
	@cmd /c "del /Q $(TARGET).exe 2>NUL"
else
	@rm -f *.o $(TARGET)
endif

.PHONY: all clean
