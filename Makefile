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