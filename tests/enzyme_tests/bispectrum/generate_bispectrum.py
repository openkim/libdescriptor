import numpy as np
from kliff.descriptors import Bispectrum
from kliff.dataset import Dataset

cut_name = "cos"
cut_dists = {"Si-Si": 3.77}
hyper_params = {'jmax': 4, 'weight':{"Si": 1.0}}

descriptor = Bispectrum(cut_dists, cut_name, hyper_params)

data = Dataset("test_00000.xyz")
config = data.get_configs()[0]

desc = descriptor.transform(config)

np.savetxt("kliff_desc.dat", desc[0])

params = descriptor.get_hyperparams()

with open("hyper_param.dat", 'w') as f:
    for param in params:
        f.write(param + "\n")
        f.write(str(params[param]) + "\n")
        f.write("\n")

