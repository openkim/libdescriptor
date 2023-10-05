import os
import sys

current_folder = os.path.dirname(os.path.realpath(__file__))
if current_folder not in sys.path:
    sys.path.append(current_folder)

from . import libdescriptor
