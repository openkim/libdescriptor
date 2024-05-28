# Compiler
CXX := clang++
CXXFLAGS := -std=c++17 -O3 -shared -fpic -fuse-ld=lld -flto  -Iinclude

#Enzyme lib location
ENZYME_LIB := /opt/enzyme/enzyme/build/Enzyme/LLDEnzyme-13.so

# Directories
SRCDIR := src
INCDIR := include
OBJDIR := obj
BINDIR := bin

# Output library
TARGET := $(BINDIR)/libdescriptor.so

# Source files
SRCFILES := $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/*/*.cpp)
# Object files
OBJFILES := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES))

# Create directories if they don't exist
$(shell mkdir -p $(BINDIR) $(OBJDIR) $(OBJDIR)/maths)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJFILES)
	$(CXX)  -flto $(CXXFLAGS) -o $@ $^  -Wl,--lto-legacy-pass-manager -Wl,-mllvm=-load=$(ENZYME_LIB) -Wl,-mllvm=-enzyme-loose-types

# Compile source files into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean
#
#/opt/homebrew/Cellar/llvm@14/14.0.6/bin/clang++ -O3 -fuse-ld=lld -std=c++14 -dynamiclib  `python3.12 -m pybind11 --includes` python_bindings.cpp `python3.12-config --ldflags` -L/opt/homebrew/Cellar/python@3.12/3.12.3/lib  -lpython3.12  -I./include -L/Users/ec2-user/libdescriptor/bin -ldescriptor -o libdescriptor`python3.12-config --extension-suffix`
#
#/opt/homebrew/Cellar/llvm@14/14.0.6/bin/clang++ -O3 -std=c++14 -dynamiclib  `python3.12 -m pybind11 --includes` python_bindings.cpp `python3.12-config --ldflags` -L/opt/homebrew/Cellar/python@3.12/3.12.3/lib  -lpython3.12  -I./include -L/Users/ec2-user/libdescriptor/bin -ldescriptor -o libdescriptor`python3.12-config --extension-suffix`