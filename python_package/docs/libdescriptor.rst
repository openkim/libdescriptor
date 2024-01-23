Libdescriptor API
===================================

The core libdescriptor pybind11 module can be imported as

.. code-block:: python

    from libdescriptor import libdescriptor as c_lds

Usually you would not need it, as the libdescriptor module provides
a more convenient interface to the library, but it is documented here
for lower level access.
Below is a list of the routines and classes that are available in the
libdescriptor library python module.

.. class:: AvailableDescriptors

    An enum class that represents the available descriptors.

    .. attribute:: SymmetryFunctions

        The enum for symmetry functions descriptor, mapped to ``Descriptor::AvailableDescriptor::KindSymmetryFunctions`` in C++. In python Enum, this value is mapped to integer ``0``.

    .. attribute:: Bispectrum

        The bispectrum descriptor, mapped to ``Descriptor::AvailableDescriptor::KindBispectrum`` in C++. In python Enum, this value is mapped to integer ``1``.

    .. attribute:: SOAP

        The SOAP descriptor, mapped to ``Descriptor::AvailableDescriptor::KindSOAP`` in C++. In python Enum, this value is mapped to integer ``2``.

    .. attribute:: Xi

        The Xi descriptor, mapped to ``Descriptor::AvailableDescriptor::KindXi`` in C++. In python Enum, this value is mapped to integer ``3``.



.. class:: DescriptorKind

    A class that represents a descriptor. It is directly exposing the
    C++ ``DescriptorKind`` class.

    :ivar kind:

        The ``enum`` kind of the descriptor.

    :ivar width:

        The width of the descriptor.

    :ivar param_file:

        The parameter file of the descriptor, this file is used in initialization of the class in C++ and archiving.


    .. method:: compute(self, index: int, species: np.ndarray[int32], neighbors: np.ndarray[int32], coordinates: np.ndarray)

            Compute the descriptor for the single atom with given atom index, species, neighbors and coordinates.

            :param index: The index of the atom for which the descriptor is computed.
            :param species: The species of the atoms.
            :param neighbors: The neighbors index of the atoms.
            :param coordinates: The coordinates of the atoms.

            :return: The descriptor as a numpy array.


    .. method:: init_descriptor(kind: AvailableDescriptors) -> DescriptorKind

            Initialize the empty descriptor class with the given kind.

            :param kind: The kind of the descriptor.

            :return: The initialized descriptor.

    .. method:: init_descriptor(file_name: str, kind: AvailableDescriptors) -> DescriptorKind

            Initialize the descriptor class with the given kind and parameter file.

            :param file_name: The parameter file name.
            :param kind: The kind of the descriptor.

            :return: The initialized descriptor.

    .. method:: init_descriptor(kind: AvailableDescriptor, species: List[str], cutoff_fun: str, cutoff_mat: np.ndarray, sym_fun_list: List[str],sym_fun_sizes: List[int],sym_fun_param: List[float]) -> DescriptorKind

            Initialize the descriptor class with the symmetry functions and parameters.

            :param kind: The kind of the descriptor, should be :attr:`Available.SymmetryFunctions` kind.
            :param species: The species of the atoms.
            :param cutoff_fun: The cutoff function.
            :param cutoff_mat: The cutoff matrix.
            :param sym_fun_list: The symmetry functions ``g1`` ... ```g5`.
            :param sym_fun_sizes: The symmetry function sizes.
            :param sym_fun_param: The symmetry function parameters.

            :return: The initialized descriptor.

    .. method:: init_descriptor(kind: AvailableDescriptors, rfac0: float, twojmax: int, diagonalstyle: int, shared_array: int, rmin0: float, switch_flag: int, bzero_flag: int, cutoff_array: np.ndarray, species: List[str], weights: List[float]) -> DescriptorKind

            Initialize the descriptor class with the bispectrum parameters. They follow exact same meaning as their LAMMPS
            counterparts.

            :param kind: The kind of the descriptor, should be :attr:`Available.Bispectrum` kind.
            :param rfac0: The rfac0 parameter.
            :param twojmax: The twojmax parameter.
            :param diagonalstyle: The diagonalstyle parameter.
            :param shared_array: The shared_array parameter.
            :param rmin0: The rmin0 parameter.
            :param switch_flag: The switch_flag parameter.
            :param bzero_flag: The bzero_flag parameter.
            :param cutoff_array: The cutoff_array parameter.
            :param species: The species of the atoms.
            :param weights: The weights of the atoms.

            :return: The initialized descriptor.

    .. method:: init_descriptor(kind: AvailableDescriptors, n_max: int, l_max: int, cutoff: float, species: List[str], radial_basis: str, eta: float) -> DescriptorKind

            Initialize the descriptor class with the SOAP parameters. They follow exact same meaning as their LAMMPS
            counterparts.

            :param kind: The kind of the descriptor, should be :attr:`Available.SOAP` kind.
            :param n_max: Number of radial basis functions to use.
            :param l_max: Maximum degree of spherical harmonics.
            :param cutoff: The cutoff parameter.
            :param species: List of species to consider.
            :param radial_basis: Radial basis function to use. Currently supported are "polynomial".
            :param eta: The gaussian width parameter of the radial basis function.

            :return: The initialized descriptor.

    .. method:: init_descriptor(kind: AvailableDescriptors, l_max: int, cutoff: float, species: List[str], radial_basis: str) -> DescriptorKind

            Initialize the descriptor class with the Xi parameters. They follow exact same meaning as their LAMMPS
            counterparts.

            :param kind: The kind of the descriptor, should be :attr:`Available.Xi` kind.
            :param l_max: Maximum degree of spherical harmonics.
            :param cutoff: The cutoff parameter.
            :param species: List of species to consider.
            :param radial_basis: Radial basis function to use. Currently supported are "bessel".

            :return: The initialized descriptor.

.. method:: compute_single_atom(descriptor_class: DescriptorKind, index: int, species: np.ndarray[int32], neighbor_idx: np.ndarray[int32], coordinates: np.ndarray) -> np.ndarray

    Compute the descriptor for the single atom with given atom index, species, neighbors and coordinates.

    :param descriptor_class: The initialized descriptor class.
    :param index: The index of the atom for which the descriptor is computed.
    :param species: The array of species indexes of the atoms.
    :param neighbor_idx: The array of neighbors indexes of the atoms.
    :param coordinates: The array of coordinates of the atoms.

    :return: The descriptor as a numpy array.


.. method:: gradient_single_atom(descriptor_class: DescriptorKind, index: int, species: np.ndarray[int32], neighbor_idx: np.ndarray[int32], coordinates: np.ndarray, computed_desc: np.ndarray, dE_dzeta: np.ndarray) -> np.ndarray

    Compute the gradient of the descriptor for the single atom with given atom index, species, neighbors and coordinates.
    This method computes the vector-Jacobian product of the descriptor function with respect to incoming :math:`\frac{dE}{d\zeta}` vector.

    :param descriptor_class: The initialized descriptor class.
    :param index: The index of the atom for which the descriptor is computed.
    :param species: The array of species indexes of the atoms.
    :param neighbor_idx: The array of neighbors indexes of the atoms.
    :param coordinates: The array of coordinates of the atoms.
    :param computed_desc: The computed descriptor of the environment.
    :param dE_dzeta: The gradient of the energy with respect to the descriptor.

    :return: The gradient of the descriptor as a numpy array.


.. method:: compute(descriptor_class: DescriptorKind, n_atoms: int, species: np.ndarray[int32], neighbor_idx: np.ndarray[int32], num_neighbors: np.ndarray[int32], coordinates: np.ndarray) -> np.ndarray

    Compute the descriptor for the atoms with given species, neighbors and coordinates. This method is more efficient than
    calling :meth:`compute_single_atom` for each atom.

    :param descriptor_class: The initialized descriptor class.
    :param n_atoms: The number of atoms.
    :param species: The array of species indexes of the atoms.
    :param neighbor_idx: The array of neighbors indexes of the atoms.
    :param num_neighbors: The array of number of neighbors of the atoms.
    :param coordinates: The array of coordinates of the atoms.

    :return: The descriptor as a numpy array.

.. method:: gradient(descriptor_class: DescriptorKind, n_atoms: int, species: np.ndarray[int32], neighbor_idx: np.ndarray[int32], num_neighbors: np.ndarray[int32], coordinates: np.ndarray, computed_desc: np.ndarray, dE_dzeta: np.ndarray) -> np.ndarray

    Compute the gradient of the descriptor for the atoms with given species, neighbors and coordinates. This method is more efficient than
    calling :meth:`gradient_single_atom` for each atom.

    :param descriptor_class: The initialized descriptor class.
    :param n_atoms: The number of atoms.
    :param species: The array of species indexes of the atoms.
    :param neighbor_idx: The array of neighbors indexes of the atoms.
    :param num_neighbors: The array of number of neighbors of the atoms.
    :param coordinates: The array of coordinates of the atoms.
    :param computed_desc: The computed descriptor of the environment.
    :param dE_dzeta: The gradient of the energy with respect to the descriptor.

    :return: The gradient of the descriptor as a numpy array.

.. method:: jacobian(descriptor_class: DescriptorKind, n_atoms: int, species: np.ndarray[int32], neighbor_idx: np.ndarray[int32], num_neighbors: np.ndarray[int32], coordinates: np.ndarray) -> np.ndarray

    Compute the jacobian of the descriptor for the atoms with given species, neighbors and coordinates. Usually you should not
    calculate the jacobian of the descriptor, but rather use the :meth:`gradient` method, which is more efficient. This is only
    useful if you want to calculate the jacobian of the descriptor for some other reason.

    :param descriptor_class: The initialized descriptor class.
    :param n_atoms: The number of atoms.
    :param species: The array of species indexes of the atoms.
    :param neighbor_idx: The array of neighbors indexes of the atoms.
    :param num_neighbors: The array of number of neighbors of the atoms.
    :param coordinates: The array of coordinates of the atoms.

    :return: The jacobian of the descriptor as a numpy array.