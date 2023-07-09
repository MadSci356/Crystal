# running in Python 3
# ----------script for creating crystals----------#
#
# Uses Lattice and Atom classes
# runs in docker with:
# exec(open("crystal.py").read())
#
# Author: Sayam Patel (sdpate11@ncsu.edu)
# Date created: 8/3/2018
#
# ------------------------------------------------#
from math import sqrt
from crystal.classes.lattice import Lattice
from crystal.classes.latticeio import LatticeIO
from crystal.classes.atom import Atom

# making a lattice basis vectors
angles = (0, 0, 120)  # angles for the lattice vector
sqrt3 = sqrt(3)
basis_vectors = [(0, 0, 0), (.5, sqrt3/6, 0)]

dims = (1, 1, 0)  # dimenstions of the lattice (x, y, z) steps
a = 8  # unit cell size


# getting a lattice object from angles
# by default, basis vector of origin is included when lattice object is made
graphene = Lattice.from_angles("graphene", a, angles, dims)
graphene.basis_vectors = basis_vectors

# creating atom to add to lattice
radius = 1
carbon_a = Atom("CarbonA", radius)
carbon_b = Atom("CarbonB", radius + .5)

# adding CarbonB bond to carbonA (note: don't need to add a to b, adding b to a
# is good enough
carbon_a.bonds_to = [carbon_b]

atom_list = [carbon_a, carbon_b]  # can be in a list (for multiple atoms)
graphene.atoms = atom_list

# print(graphene.info())
# making io object
latticeio = LatticeIO(graphene)
# filename of the output obj file
filename = "testing/object_oriented/graphene11b.obj"

# making obj file
latticeio.make_obj_file(filename)
