from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "libdescriptor",
        sorted(glob("*.cpp") + glob("**/*.cpp")),  # Sort source files for reproducibility
        cxx_std=17,
        include_dirs=[".","./autodiff"],
    ),
]

setup(
    name="libdescriptor",
    version="0.5",
    author="Amit Gupta",
    author_email="gupta839@umn.edu",
    url="https://github.com/ipcamit/libdescriptor",
    description="auto differentiated descriptor library",
    long_description="",
    ext_modules=ext_modules,
    # extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    # cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
