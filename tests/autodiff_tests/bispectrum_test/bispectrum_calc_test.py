import libdescriptor as lds
import numpy as np
from collections import OrderedDict

def test(tol=1e-6):
    n_atoms = 8
    coords = np.loadtxt('data/coords.dat')
    n_neigh = np.loadtxt('data/n_neigh.dat', dtype=int)
    neigh = np.loadtxt('data/neigh.dat', dtype=int)
    fwd = np.loadtxt('data/fwd.dat')
    rev = np.loadtxt('data/rev.dat')
    image = np.loadtxt('data/image.dat', dtype=int)

    # Create descriptor object
    cutoff = 3.77
    hyp_param = OrderedDict({"jmax": 4, "rfac0": 0.99363, "diagonalstyle": 3, "rmin0": 0,
            "switch_flag": 1, "bzero_flag": 0, "use_shared_array": False, "weights": None})
    species = ["Si"]
    weights = np.ones(len(species))
    cutoff_array = np.ones((len(species), len(species))) * cutoff
    desc = lds.DescriptorKind.init_descriptor(lds.AvailableDescriptors(1), hyp_param["rfac0"], 
            2 * hyp_param["jmax"], hyp_param["diagonalstyle"], 
            1 if hyp_param["use_shared_array"] else 0, hyp_param["rmin0"],
            hyp_param["switch_flag"], hyp_param["bzero_flag"], cutoff_array, species, weights)
    species_arr = np.zeros(coords.shape[0], dtype=np.intc)

    err = 0
    for i in range(n_atoms):
        fwd_new = lds.compute_single_atom(desc ,i, species_arr, np.array(neigh[i], dtype=np.intc), coords)
        err += np.sum(np.abs(fwd_new - fwd[i]))

    rev = np.loadtxt("data/rev.dat")
    
    err_grad = 0
    rev_new = np.zeros((np.max(image)+1 ,3))
    
    for i in range(n_atoms):
        d_coords = lds.gradient_single_atom(desc ,i, species_arr, np.array(neigh[i], dtype=np.intc), coords, fwd_new, np.ones_like(fwd_new))
        for j,k in enumerate(image):
            rev_new[k,:] += d_coords[j,:]
    
    err_grad = np.sum(np.abs((rev_new - rev)))

    print("Error in descriptor: ", err)
    print("Error in gradient: ", err_grad)

    assert err < tol, "Error in forward pass is too large"
    assert err_grad < tol, "Error in backward pass is too large"

    # TODO: numerical gradient check
    # err_grad_num = 0
    # rev_new_num = np.zeros((np.max(image)+1 ,3))
    # for i in range(n_atoms):
    #     d_coords = lds.num_gradient_single_atom(desc, i, n_atoms,  species_arr, np.array(neigh[i], dtype=np.intc), coords, np.ones_like(fwd_new))
    #     rev_new_num[i] += d_coords

    # err_grad_num = np.sum(np.abs((rev_new_num - rev)))
    # print(err_grad_num)

if __name__ == "__main__":
    test()
