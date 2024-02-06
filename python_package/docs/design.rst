Design
======

The python package is just a pybind11 wrapper around the C++ code.
The C++ code is designed in the following hierarchy:

.. figure:: /_static/libdescriptorDesign.svg
   :align: center
   :width: 100%

   UML diagram of the C++ code

Here namespace ``Descriptor`` contains the main ``DescriptorKind`` class, enumerated types for
supported descriptor types (``AvailableDescriptors`` and all the utility functions for descriptor
generation and manipulation.

The ``DescriptorKind`` class is the main class that is used to generate descriptors, and it is
interfaced by the individual descriptor classes. It is done by inheriting from the ``DescriptorKind``
and implementing the pure virtual function, ``DescriptorKind::compute``. Hence all of the computational
work is done in the ``DescriptorKind::compute`` function. ``compute`` implements calculation for a single atomic environment,
and gradients are generated automatically in the ``Descriptor::gradient`` and ``Descriptor::gradient_single_atom`` functions.

To ensure that your implemented descriptor class does get differentiated, you need to hook it in the ``Descriptor::compute...`` and
``Descriptor::gradient...`` functions. This is done by adding a new case in the switch statement in the ``Descriptor::compute...`` and
``Descriptor::gradient...`` functions. You would also need to implement appropriate ``DescriptorKind::initDescriptor`` overloaded function.
This is done so that at runtime any descriptor can be initialized with the ``DescriptorKind`` class, and the individual
descriptors implementations are abstracted away. This ensures a more portable implementation for end user cases.

For hooking the gradient function you also need to implement a ``clone_empty`` function in your descriptor class. This is used to
generate a ``__enzyme_virtual`` class for the gradient computation. This class contains gradients against all of your class members
but currently this class is deleted post computation. In future releases I plan to provide a way to access these gradients as well.

Please let me know if you have any questions or need help with implementing a new descriptor.
