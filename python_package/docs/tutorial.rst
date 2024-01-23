Tutorial
========

This tutorial is a step-by-step guide for computing supported descriptors in libdescriptor.
Currently, the following descriptors are supported:

1. Behler-Parrinello symmetry functions
2. Bispectrum coefficients
3. SOAP descriptors
4. Xi [Experimental]


List available descriptors
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can get the supported descriptors by using the enumerated lists in the available.

.. code-block:: python

    import libdescriptor as lds

    for i in range(5):
        print(lds.AvailableDescriptors(i))

The output will be:

.. code-block:: python

    AvailableDescriptors.SymmetryFunctions
    AvailableDescriptors.Bispectrum
    AvailableDescriptors.SOAP
    AvailableDescriptors.Xi
    AvailableDescriptors.???

``???`` indicates that only first 4 descriptors are supported, and there is no 5th descriptor.
In future there will be better utility functions to get the supported descriptors.


Example evaluation of descriptors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us create a simple system to compute the descriptors.

.. code-block:: python

    import numpy as np
    from ase import Atoms

    # Create a simple system
    n_atoms = 10
    positions = np.random.rand(n_atoms, 3) * 10
    symbols = ['Si'] * n_atoms
    species = ['Si']
    cell = np.eye(3) * n_atoms

    # Create an ASE atoms object
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)


Now, we can compute the descriptors. To do so, we need to get all the neighbors for each contributing atoms.
Libdescriptor provides a neighbor list class to do so, :class:`libdescriptor.NeighborList`.
Pass a ASE atoms object to the constructor of the neighbor list class. You can use the function
:func:`libdescriptor.NeighborList.get_numneigh_and_neighlist_1D` to get the number of neighbors and the neighbor list.

.. code-block:: python

    # Get neighbor list
    from libdescriptor.neighbor import NeighborList

    neighbor_list = NeighborList(atoms, 4.0)
    element_map = {"Si": 0}
    num_neighs, neigh_list = neighbor_list.get_numneigh_and_neighlist_1D()
    coordinates = neighbor_list.get_coords() # get padded coordinates
    species_indices = neighbor_list.get_species_code(element_map) # get species indices


Now, we can create a descriptor object. From the enumerated list of available descriptors, first pick the desired descriptor.
Then, use the function :func:`libdescriptor.DescriptorKind.init_descriptor` along with either with
**filename of the descriptor file containing parameters** or **a dictionary of hyperparameter arguments**.
To see list of supported hyperparameters, see the documentation of the descriptor class :class:`libdescriptor.DescriptorKind`.
After creating the descriptor object, we can compute the descriptor by using the function :func:`libdescriptor.compute`.

.. code-block:: python

    # Create a descriptor object [SOAP]
    import libdescriptor as lds

    n_max = 4
    l_max = 3
    cutoff = 4.0
    species = ['Si']
    basis = 'polynomial'
    eta = 0.1
    kind = lds.AvailableDescriptors(2) # SOAP

    descriptor = lds.DescriptorKind.init_descriptor(kind, n_max, l_max, cutoff, species, basis, eta)

You can check the width of the descriptor by checking the property :attr:`libdescriptor.DescriptorKind.width`.

.. code-block:: python

    print(f"SOAP width for l={l_max} and n_max={n_max} = {descriptor.width}")
    # should print 40

Now, we can compute the descriptor by using the function :func:`libdescriptor.compute`.

.. code-block:: python

    # Compute the descriptor
    desc = lds.compute(descriptor, 10, # number of atoms
                species_indices.astype(np.intc), # species indices
                neigh_list.astype(np.intc), # neighbor list
                num_neighs.astype(np.intc), # number of neighbors
                coordinates)

For computing reverse-mode autodiﬀ, or vector-Jacobian product, we need to call the gradient function with
the vector as input (for which we need to compute the gradient, in normal ML applications this would be
derivative of descriptor against the energy). You would also need to pass the previously computed descriptor
to the gradient function.

.. code-block:: python

    # Compute the gradient
    dE_dzeta = np.ones(descriptor.width)
    gradient = lds.gradient(descriptor, 10, # number of atoms
                 species_indices.astype(np.intc), # species indices
                 neigh_list.astype(np.intc), # neighbor list
                 num_neighs.astype(np.intc), # number of neighbors
                 coordinates, desc , dE_dzeta)

    print(gradient)
    # should print the gradient of the descriptor