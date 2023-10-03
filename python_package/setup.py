from setuptools import setup, find_packages

setup(
    name="libdescriptor_dev",
    version="0.0.3",
    packages=find_packages(),
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
    include_package_data=True,
)
