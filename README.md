Libdescriptor
======================================

Libdescriptor is a high performance descriptor library for providing access to fully differentiable descriptor functions.
While `libdescriptor` is a general purpose descriptor library, it's API compatible with KIM models and associated projects.
This will also provide uniform access to various selected descriptors for KLIFF using Pybind11 ports.
For gradient calculations, Libdescriptor relies on [Enzyme AD](https://github.com/EnzymeAD/Enzyme), which provides it with capability to trivially generate near analytical performance gradient functions.
Use of Enzyme AD enables Libdescriptor to not only provide gradients against coordinates, but against hyperparameters as well, thus opening way for better optimized descriptors.
This should enable rapid development, extension and deployment of various descriptors.

## Compiling 
At present, it needs functioning Enzyme compiler environment, in future it will be provided as binary package or a conda environment.
For Installing Enzyme, simply follow the instructions given on enzyme page. At the time of writing these instructions, we are using
- LLVM 12.0.1
- Enzyme v0.0.41

Although Enzyme recommends compiling LLVM/Clang from scratch, in our experience precompiled stock binaries also work fine.
Though depending on your platform, your mileage may wary.

Steps to compile **Enzyme**:
```shell
# Get the Clang/LLVM binaries
wget https://github.com/llvm/llvm-project/releases/download/llvmorg-12.0.1/clang+llvm-12.0.1-x86_64-linux-gnu-ubuntu-16.04.tar.xz
# Untar and export environment variables
tar -xvf clang+llvm-12.0.1-x86_64-linux-gnu-ubuntu-16.04.tar.xz
export PATH="/path/to/clang/extract/bin:"$PATH
export INCLUDE="/path/to/clang/extract/include:"$INCLUDE
export LD_LIBRARY_PATH="/path/to/clang/extract/lib:"$LD_LIBRARY_PATH

# Get lit (Enzyme needs it for running tests)
wget https://files.pythonhosted.org/packages/7c/0c/2d58790cb0fa24812382289a584e05dd1df4b30ccf5e2218ee5a556a0529/lit-12.0.1.tar.gz
tar -xvf lit-12.0.1.tar.gz
cd lit-12.0.1
pip install . --user

# Clone and compile Enzyme
git clone https://github.com/EnzymeAD/Enzyme
cd Enzyme/enzyme
git checkout v0.0.41
mkdir build; cd build
CC=clang CXX=clang++ cmake .. -DLLVM_DIR=/path/to/clang/lib/cmake/llvm -DLLVM_EXTERNAL_LIT=/path/to/lit-12.0.1/lit.py
make 
# optional make install if you want
```
Now your enzyme is ready to use. For compiling your code, you would need to know the location of your compiled Enzyme libraries. 
In the build folder you should see 3 shared objects: `ClangEnzyme-12.so`, `LLDEnzyme-12.so`, and `LLVMEnzyme-12.so`.
As a rule of thumb, you need `ClangEnzyme` for compiling single file programs, whereas more complicated build schemes (such as used in this program),
needs creating derivatives at link time, using `LLDEnzyme` file. 

You need to provide the location of `LLDEnzyme` to the `Cmake` file for successful build.
This can be done simply by defining Cmake variable. To build `libdescritpor`
```shell
git clone https://github.com/ipcamit/colabfit-descriptor-library
cd colabfit-descriptor-library
mkdir build; cd build
cmake .. -DENZYME_LIB=/path/to/*Enzyme.so/files
make
```
Your build folder should now contain `libdescriptor.so` file, which you can link against your own projects.


## Python bindings
**WIP**

Libdescriptor also provides build target for making `descriptor.cpython-cp3.xx.so` Python module using Pybind11.
To use it, you need to install Pybind11 on your system
> Note: pip installation of Pybind11 can cause issues with compiling (it is a known Pybind11 limitation) 
> so it is recommended to use either system installer like `apt` or `rpm`  or use `conda`.

## Descriptors supported (or planned)
- [x] Behler Symmetry Functions
- [x] Bispectrum 
- [ ] SOAP (WIP)
- [ ] ACE

## Extending Libdescriptor
New descriptors can be added by extending the `DescriptorKind` class, and implementing its `compute` function.
`compute` function will compute the descriptor for single atom. It is also advised to write a `clone_empty` function if you want
to take derivatives against hyperparamters. Without an empty clone, Enzyme can segfault at times when default constructor does not initialize
all fields of class properly.

> When in doubt, make your code more "C" like for higher success rate.

TODO: More documentation for extending.