.. Libdescriptor documentation master file, created by
   sphinx-quickstart on Tue Oct 17 11:07:51 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Libdescriptor: Fully differentiable atomic descriptors
======================================================
Libdescriptor is a library for computing atomic descriptors. It is written in C++ and has
extensive Python bindings. It is designed to be fast and flexible, and to be easily extended
with new descriptors. Its major design goal include using automatic differentiation to compute
derivatives of the descriptors, and to be able to compute derivatives of the descriptors with
respect to the descriptor parameters. It uses state of the art `Enzyme AD <https://github.com/EnzymeAD/Enzyme>`__  library for automatic
differentiation.

Libdescriptor comes bundled with minimal support for neighbor lists, and ability to compute descriptors for
ASE atoms objects. Libdescriptor provides a unified API to compute and train the descriptors in Python
and transfer the descriptors to C++ for production MD simulations.

.. image:: /_static/libdescriptor.svg
    :align: center
    :alt: libdescriptor

Contents:
--------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   installation
   design
   enzyme
   libdescriptor
   neighbor


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
