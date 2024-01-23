Installation
============

Installation from PyPI
----------------------

For ease of use, Libdescriptor is available on PyPI and can be installed with pip:

.. code-block:: bash

    pip install libdescriptor

This method relies on binary wheels, which will have limited availability for some platforms.
If you encounter problems, please file an issue. At present it only support Linux on x86_64.

.. note::

    Please note that **this should be the preferred way of installing** libdescriptor if possible
    as compiling from source requires a valid installation of LLVM and clang compilers, and is not trivial.

Installation from source
------------------------
For unsupported platforms, or if you want to hack on libdescriptor, you can install from source.
This requires a valid installation of LLVM >= 13.0.0 and clang >= 13.0.0 compilers.

First step requires that you install the `Enzyme AD library <https://github.com/enzymeAD/enzyme>`__ as per
the instructions provided by the library. This will install the `enzyme` python package.
If Enzyme is installed in some non standard location, you should specify the environment variable
``ENZYME_LIB`` to point to the installation directory.

1. Then, clone the repository:

.. code-block:: bash

    git clone https://github.com/ipcamit/libdescriptor
    cd libdescriptor

2. Build the library:

.. code-block:: bash

    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DENZYME_LIB=/path/to/enzyme/installation
    make

3. Once you have ``libdescriptor.so`` and ``libdescriptor.cpython-*.so`` in the build directory,
   you can copy these file in the ``python_package`` directory of the repository and install the
   python package:

.. code-block:: bash

    cp libdescriptor.so libdescriptor.cpython-*.so ../python_package/libdescriptor
    cd ../python_package
    pip install .

