# ===============================
# Compiler
# ===============================
CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -pedantic-errors -g \
            -msse2 -march=native -funroll-loops -ffast-math \
            -fomit-frame-pointer -fstrict-aliasing
LDFLAGS := -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive,-Bdynamic

# ===============================
# Source directory (unchanged)
# ===============================
SDIR := ../src
SRCS := $(wildcard $(SDIR)/*.cpp)

# ===============================
# Object files in PWD
# ===============================
OBJS := $(notdir $(SRCS:.cpp=.o))

TARGET := sse_sim_test

# ===============================
# Default target
# ===============================
all: $(TARGET)
	@echo Compilation Successful

# ===============================
# Link
# ===============================
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

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

# ===============================
# Debug target — enables SSE_DEBUG
# Recompiles everything from scratch with -DSSE_DEBUG so the diagnostic
# output (init check, op-string dump, linked list, loop BFS trace) is active.
# Usage:  make debug     (from backend/build/)
# ===============================
.PHONY: debug
debug:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DSSE_DEBUG"

# ===============================
# Verbose target — enables VERBOSE print output
# Recompiles with -DVERBOSE so parameter dumps and [SSE-CHECK] lines are shown.
# Usage:  make verbose   (from backend/build/)
# ===============================
.PHONY: verbose
verbose:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DVERBOSE"

.PHONY: all clean
