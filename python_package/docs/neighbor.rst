Neighbor API
============

This is a bare bones neighbor lists support so that for basic utilization there is no external library needed.
This module was taken from the KLIFF package and modified to support ASE Atoms object, instead of KLIFF Configuration object.
This makes ASE a hard dependency.

NeighborList can be imported from the libdescriptor module as


.. code-block:: python

    from libdescriptor.neighbor import NeighborList


Classes
-------

.. class:: NeigborList

    This class accepts an ASE atoms configuration and cutoff.

    .. method:: get_neigh(index)

        Get neighbors for index

    .. method:: get_numneigh_and_neighlist_1D()

        Return list number of neighbors and linear neighbor list for the configuration

    .. method:: get_coords()

        Get the complete padded coordinates.

    .. method:: get_species()

    .. method:: get_species_code(mapping)

    .. method:: get_image()

        Get padded images

    .. method:: get_padding_coords()

    .. method:: get_padding_species()

    .. method:: get_padding_species_code(mapping)

    .. method:: get_padding_image()

