import numpy as np
from kliff.descriptors import SymmetryFunction
from kliff.dataset import Dataset

descriptor = SymmetryFunction(
    cut_name="cos", cut_dists={"Si-Si": 3.77}, hyperparams="set51", normalize=False
)

data = Dataset("test_00000.xyz")
config = data.get_configs()[0]

desc = descriptor.transform(config)

print(desc[0])
np.savetxt("kliff_desc.dat", desc[0])

params = descriptor.get_hyperparams()

with open("hyper_param.dat", 'w') as f:
    for param in params:
        f.write(param + "\n")
        f.write(str(params[param]) + "\n")
        f.write("\n")

