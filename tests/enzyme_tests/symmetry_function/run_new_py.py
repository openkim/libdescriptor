import libdescriptor
ds = libdescriptor.DescriptorKind.init_descriptor("descriptor.dat", libdescriptor.AvailableDescriptors.SymmetryFunctions)
from kliff.neighbor import NeighborList
from kliff.dataset import Dataset
data = Dataset("test_00000.xyz")
config = data.get_configs()
nl = NeighborList(config[0], 3.77)
neig, coords, species = nl.get_neigh(0)
import numpy as np
desc = libdescriptor.compute_single_atom(ds, 0, np.zeros(len(nl.coords), dtype=int), neig, nl.coords)
print(desc)

grad = libdescriptor.gradient_single_atom(ds, 0, np.zeros(len(nl.coords), dtype=int), neig, nl.coords, desc, np.ones_like(desc))
print("AD: ", grad[0])

num_grad = libdescriptor.num_gradient_single_atom(ds, 0, np.zeros(len(nl.coords), dtype=int), neig, nl.coords, np.ones_like(desc))
print("NUM: ", num_grad)
print("MAD: ", np.mean(np.abs(num_grad - grad[0])))
