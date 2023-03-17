import libdescriptor as lds
import numpy as np
from collections import OrderedDict

def test(tol = 1e-10):
    n_atoms = 8
    coords = np.loadtxt('data/coords.dat')
    n_neigh = np.loadtxt('data/n_neigh.dat', dtype=int)
    neigh = np.loadtxt('data/neigh.dat', dtype=int)
    fwd = np.loadtxt('data/fwd.dat')
    rev = np.loadtxt('data/rev.dat')
    image = np.loadtxt('data/image.dat', dtype=int)
    
    # Create the descriptor object
    hyp_param =  OrderedDict([('g2',
                  [{'eta': 0.0035710676725828126, 'Rs': 0.0},
                   {'eta': 0.03571067672582813, 'Rs': 0.0},
                   {'eta': 0.07142135345165626, 'Rs': 0.0},
                   {'eta': 0.12498736854039845, 'Rs': 0.0},
                   {'eta': 0.21426406035496876, 'Rs': 0.0},
                   {'eta': 0.3571067672582813, 'Rs': 0.0},
                   {'eta': 0.7142135345165626, 'Rs': 0.0},
                   {'eta': 1.428427069033125, 'Rs': 0.0}]),
                 ('g4',
                  [{'zeta': 1, 'lambda': -1, 'eta': 0.00035710676725828126},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.00035710676725828126},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.00035710676725828126},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.00035710676725828126},
                   {'zeta': 1, 'lambda': -1, 'eta': 0.010713203017748437},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.010713203017748437},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.010713203017748437},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.010713203017748437},
                   {'zeta': 1, 'lambda': -1, 'eta': 0.0285685413806625},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.0285685413806625},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.0285685413806625},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.0285685413806625},
                   {'zeta': 1, 'lambda': -1, 'eta': 0.05356601508874219},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.05356601508874219},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.05356601508874219},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.05356601508874219},
                   {'zeta': 4, 'lambda': -1, 'eta': 0.05356601508874219},
                   {'zeta': 4, 'lambda': 1, 'eta': 0.05356601508874219},
                   {'zeta': 16, 'lambda': -1, 'eta': 0.05356601508874219},
                   {'zeta': 16, 'lambda': 1, 'eta': 0.05356601508874219},
                   {'zeta': 1, 'lambda': -1, 'eta': 0.08927669181457032},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.08927669181457032},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.08927669181457032},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.08927669181457032},
                   {'zeta': 4, 'lambda': -1, 'eta': 0.08927669181457032},
                   {'zeta': 4, 'lambda': 1, 'eta': 0.08927669181457032},
                   {'zeta': 16, 'lambda': -1, 'eta': 0.08927669181457032},
                   {'zeta': 16, 'lambda': 1, 'eta': 0.08927669181457032},
                   {'zeta': 1, 'lambda': -1, 'eta': 0.16069804526622655},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.16069804526622655},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.16069804526622655},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.16069804526622655},
                   {'zeta': 4, 'lambda': -1, 'eta': 0.16069804526622655},
                   {'zeta': 4, 'lambda': 1, 'eta': 0.16069804526622655},
                   {'zeta': 16, 'lambda': -1, 'eta': 0.16069804526622655},
                   {'zeta': 16, 'lambda': 1, 'eta': 0.16069804526622655},
                   {'zeta': 1, 'lambda': -1, 'eta': 0.28568541380662504},
                   {'zeta': 1, 'lambda': 1, 'eta': 0.28568541380662504},
                   {'zeta': 2, 'lambda': -1, 'eta': 0.28568541380662504},
                   {'zeta': 2, 'lambda': 1, 'eta': 0.28568541380662504},
                   {'zeta': 4, 'lambda': -1, 'eta': 0.28568541380662504},
                   {'zeta': 4, 'lambda': 1, 'eta': 0.28568541380662504},
                   {'zeta': 16, 'lambda': 1, 'eta': 0.28568541380662504}])])
    
    species = ["Si"]
    species_arr = np.zeros(coords.shape[0], dtype=np.intc)
    cutoff = 3.77
    cutoff_array = np.ones((len(species), len(species))) * cutoff
    symmetry_function_types = list(hyp_param.keys())
    symmetry_function_sizes = []
    
    cutoff_function = "cos"
    
    symmetry_function_param_matrices = []
    param_num_elem = 0
    width = 0
    for function in symmetry_function_types:
        if function.lower() not in ["g1", "g2", "g3", "g4", "g5"]:
            ValueError("Symmetry Function provided, not supported")
    
        if function.lower() == "g1":
            rows = 1
            cols = 1
            params_mat = np.zeros((1, 1), dtype=np.double)
        else:
            params = hyp_param[function]
            rows = len(params)
            cols = len(list(params[0].keys()))
            params_mat = np.zeros((rows, cols), dtype=np.double)
    
            for i in range(rows):
                if function.lower() == "g2":
                    params_mat[i, 0] = params[i]["eta"]
                    params_mat[i, 1] = params[i]["Rs"]
                elif function.lower() == "g3":
                    params_mat[i, 0] = params[i]["kappa"]
                elif function.lower() == "g4":
                    params_mat[i, 0] = params[i]["zeta"]
                    params_mat[i, 1] = params[i]["lambda"]
                    params_mat[i, 2] = params[i]["eta"]
                elif function.lower() == "g5":
                    params_mat[i, 0] = params[i]["zeta"]
                    params_mat[i, 1] = params[i]["lambda"]
                    params_mat[i, 2] = params[i]["eta"]
        symmetry_function_sizes.extend([rows, cols])
        symmetry_function_param_matrices.append(params_mat)
        param_num_elem += rows * cols
        width += rows
    
    symmetry_function_param = np.zeros((param_num_elem,), dtype=np.double)
    k = 0
    for i in range(len(symmetry_function_types)):
        symmetry_function_param[k: k + symmetry_function_sizes[2 * i] * symmetry_function_sizes[2 * i + 1]] = \
            symmetry_function_param_matrices[i].reshape(1, -1)
        k += symmetry_function_sizes[2 * i] * symmetry_function_sizes[2 * i + 1]
    
    descriptor_kind = lds.AvailableDescriptors(0)
    desc = lds.DescriptorKind.init_descriptor(descriptor_kind, species, cutoff_function, cutoff_array,
                        symmetry_function_types, symmetry_function_sizes, symmetry_function_param)
    
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

    # TODO: numerical derivative test
    # err_grad_num = 0
    # rev_new_num = np.zeros((np.max(image)+1 ,3))
    # for i in range(n_atoms):
    #     d_coords = lds.num_gradient_single_atom(desc, i, n_atoms,  species_arr, np.array(neigh[i], dtype=np.intc), coords, np.ones_like(fwd_new))
    #     rev_new_num[i] += d_coords

    # err_grad_num = np.sum(np.abs((rev_new_num - rev)))
    # print(err_grad_num)

    return 0


if __name__ == "__main__":
    test()
