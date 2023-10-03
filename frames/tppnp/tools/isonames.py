import numpy as np
import os as _os

encd = "ascii"

# create an element dictionary to map to/from z/symbol
element_file = _os.path.join(_os.path.dirname(__file__), "elements.txt")
element = dict()
with open(element_file,"r") as f:
    n = int( f.readline() )
    for i in range(n):
        spl = f.readline().split()
        z, el = np.int32(spl[0]), spl[1].encode(encd)
        element[z] = el
        element[el] = z


@np.vectorize
def ppn_name_to_a_z(iso):
    iso = iso.encode("ascii") if not isinstance(iso,bytes) else iso
    if iso == b"NEUT ":
        a = 1; z = 0
    elif iso == b"PROT ":
        a = 1; z = 1
    else:
        z = element[ iso[:2].capitalize().strip() ]
        a = np.int32(iso[2:])

    return a, z

@np.vectorize
def a_z_to_ppn_name(a,z):
    if a == 1 and z == 0:
        iso = b"NEUT "
    elif a == 1 and z == 1:
        iso = b"PROT "
    else:
        iso = element[z].upper().ljust(2) + str(a).encode(encd).rjust(3)

    return iso
