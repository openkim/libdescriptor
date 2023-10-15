from distutils.sysconfig import get_config_vars
from setuptools import Extension, find_packages, setup

# remove `-Wstrict-prototypes' that is for C not C++
cfg_vars = get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str and "-Wstrict-prototypes" in value:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


class get_pybind11_includes:
    """
    Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11 until it is actually
    installed, so that the ``get_include()`` method can be invoked.

    see:
    https://github.com/pybind/python_example/blob/master/setup.py
    https://github.com/pybind/python_example/issues/32
    """

    def __str__(self):
        import pybind11

        return pybind11.get_include()


def get_includes():
    return [get_pybind11_includes()]


def get_extra_compile_args():
    return ["-std=c++17"]

neighlist = Extension(
    "libdescriptor.neighbor.neighlist",
    sources=[
        "libdescriptor/neighbor/neighbor_list.cpp",
        "libdescriptor/neighbor/neighbor_list_bind.cpp",
    ],
    include_dirs=get_includes(),
    extra_compile_args=get_extra_compile_args(),
    language="c++",
)



setup(
    name="libdescriptor",
    version="0.0.6",
    packages=find_packages(),
    ext_modules=[neighlist],
    install_requires=[
        "numpy",
        "ase",
        "pybind11",
    ],
    package_data={
        'libdescriptor': [
                        "__init__.py",
                        "libc.so.6",
                        "libdescriptor.cpython-39-x86_64-linux-gnu.so",
                        "libdescriptor.so",
                        "libgcc_s.so.1",
                        "libm.so.6",
                        "libstdc++.so.6"
                        ],
    },
    author="Amit Gupta",
    include_package_data=True,
)
