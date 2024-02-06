import numpy
from typing import ClassVar, List, overload

class AvailableDescriptors:
    __members__: ClassVar[dict] = ...  # read-only
    Bispectrum: ClassVar[AvailableDescriptors] = ...
    SOAP: ClassVar[AvailableDescriptors] = ...
    SymmetryFunctions: ClassVar[AvailableDescriptors] = ...
    Xi: ClassVar[AvailableDescriptors] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class DescriptorKind:
    kind: AvailableDescriptors
    param_file: str
    width: int
    def __init__(self) -> None: ...
    def compute(self, arg0: int, arg1: numpy.ndarray[numpy.int32], arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
    @overload
    def init_descriptor(self) -> DescriptorKind: ...
    @overload
    def init_descriptor(self, arg0: AvailableDescriptors) -> DescriptorKind: ...
    @overload
    def init_descriptor(self, arg0: List[str], arg1: str, arg2: float, arg3: List[str], arg4: List[int], arg5: List[float]) -> DescriptorKind: ...
    @overload
    def init_descriptor(self, arg0: float, arg1: int, arg2: int, arg3: int, arg4: float, arg5: int, arg6: int, arg7: float, arg8: List[str], arg9: List[float]) -> DescriptorKind: ...
    @overload
    def init_descriptor(self, arg0: int, arg1: float, arg2: List[str], arg3: str) -> DescriptorKind: ...
    @overload
    def init_descriptor(self, arg0: int, arg1: int, arg2: float, arg3: List[str], arg4: str, arg5: float) -> DescriptorKind: ...

class NeighList:
    def __init__(self) -> None: ...
    def build(self, coords: numpy.ndarray[numpy.float64], influence_distance: float, cutoffs: numpy.ndarray[numpy.float64], need_neigh: numpy.ndarray[numpy.int32]) -> None: ...
    def get_neigh(self, cutoffs: numpy.ndarray[numpy.float64], neighbor_list_index: int, particle_number: int) -> tuple: ...

def compute(arg0: DescriptorKind, arg1: int, arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.int32], arg4: numpy.ndarray[numpy.int32], arg5: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
def compute_single_atom(arg0: DescriptorKind, arg1: int, arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.int32], arg4: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
def create() -> NeighList: ...
def create_paddings(influence_distance: float, cell: numpy.ndarray[numpy.float64], pbc: numpy.ndarray[numpy.int32], coords: numpy.ndarray[numpy.float64], species: numpy.ndarray[numpy.int32]) -> tuple: ...
def get_neigh_kim() -> capsule: ...
def gradient(arg0: DescriptorKind, arg1: int, arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.int32], arg4: numpy.ndarray[numpy.int32], arg5: numpy.ndarray[numpy.float64], arg6: numpy.ndarray[numpy.float64], arg7: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
def gradient_single_atom(arg0: DescriptorKind, arg1: int, arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.int32], arg4: numpy.ndarray[numpy.float64], arg5: numpy.ndarray[numpy.float64], arg6: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
def jacobian(arg0: DescriptorKind, arg1: int, arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.int32], arg4: numpy.ndarray[numpy.int32], arg5: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
def num_gradient_single_atom(arg0: DescriptorKind, arg1: int, arg2: numpy.ndarray[numpy.int32], arg3: numpy.ndarray[numpy.int32], arg4: numpy.ndarray[numpy.float64], arg5: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]: ...
