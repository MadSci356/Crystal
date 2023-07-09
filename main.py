# running in Python 3
# ----------script for creating crystals----------#
#
# Uses Lattice and Atom classes
# runs in docker with:
# exec(open("main.py").read())
#
# run PyMesh docker with sudo docker run -v path/to/crystal:/root -it qnzhou/pymesh
# Author: Sayam Patel (sdpate11@ncsu.edu)
# Date created: 8/3/2018
#
# ------------------------------------------------#
from math import sqrt
from lib.crystal.lattice import Lattice
from lib.crystal.latticeio import LatticeIO
from lib.crystal.atom import Atom


def make_graphene():

    # making a lattice basis vectors
    angles = (0, 0, 120)  # angles for the lattice vector
    dims = (4, 4, 0)  # dimenstions of the lattice (x, y, z) steps
    a = 8  # unit cell size

    # creating atom to add to lattice
    radius = 1
    carbon_a = Atom("CarbonA", radius, 6)
    carbon_b = Atom("CarbonC", radius, 6.1)
    bond_pattern = "CarbonA-CarbonC"
    atom_list = [carbon_a, carbon_b]  # can be in a list (for multiple atoms)
    basis_vectors = [(0, 0, 0), (.5, sqrt(3)/6, 0)]

    # getting a lattice object from angles
    graphene = Lattice.from_angles("graphene", a, angles, dims)
    graphene.basis_vectors = basis_vectors
    graphene.atoms = atom_list # adding atoms to graphene

    # making IO object
    latticeio = LatticeIO(graphene, bond_pattern)
    # filename of the output obj file
    filename = "graphene.obj"
    # making obj file
    latticeio.make_obj_file(filename)

    # CarbonB bond to carbonA
    # (note: don't need to add a to b, adding b to a
    # is good enough
    # carbon_a.bonds_to = [carbon_b]
    #  by default, basis vector of origin is included in lattice obj



def make_SrTiO3():
    # lattice params
    angles = (90, 90, 90)  # angles for the lattice vector
    basis_vectors = [(0, 0, 0), (.5, .5, .5), (.5, 0, .5),
                     (0, .5, 0), (.5, .5, 0)]
    dims = (3, 1, 3)
    a = 10
    # atom params
    r = 1.5 #radius
    strontium = Atom("Sr", r+.4, 38)
    titanium = Atom("Ti", r-.6, 22)
    oxygen1 = Atom("O1", r-.4, 8)
    oxygen2 = Atom("O2", r-.4, 8)
    oxygen3 = Atom("O3", r-.4, 8)
    atom_list = [strontium, titanium, oxygen1, oxygen2, oxygen3]
    bond_pattern = "Sr-Sr, Ti-O1, Ti-O2, Ti-O3"
    # make lattice
    perovskite = Lattice.from_angles("SrTiO3", a, angles, dims)
    perovskite.basis_vectors = basis_vectors # adding crystal basis
    perovskite.atoms = atom_list # adding atoms
    # file IO
    latticeio = LatticeIO(perovskite, bond_pattern)
    filename = "perovskite.obj"
    latticeio.make_obj_file(filename, mtl=True)


def processfile():
    latticeio = LatticeIO.from_file("input/lno_5uc.XYZ", radius=.85)
    filename = "xyzlatticemtlscale.obj"
    latticeio.make_obj_file(filename, mtl=True)

# def clean_mtl_files():
    # removes the default mtl files in the


make_SrTiO3()
#processfile()
#make_graphene()
