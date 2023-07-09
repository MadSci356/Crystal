# running in Python 3
# ----------Atom Class----------#
#
# Class definitions for better flexibility
# in creating bonds between atom. Bonds
# are 'matched' based on an atom's bonds_to list.
#
# Author: Sayam Patel (sdpate11@ncsu.edu)
# Date created: 7/30/2018
#
# ----------------------------------#
import pymesh as pm
from numpy import linalg, asarray
from numpy import round as npround

"""
Atom object that defines properties of the atom in the crystal and calculates
bond_pairs  with other atoms in its bonds_to atom list.

"""


class Atom:
    # default refinement/resolution for atom meshes
    _REFINE = 3
    # number of decimal points to round to (for vecs and dists)
    _ROUND = 3
    # max atomic number
    _MAX = 118
    # min atomic number
    _MIN = 1

    def __init__(self, name, radius, z, is_xyz=False):
        """ Attributes
        name: name of atom
        radius: radius of the atom
        z: atomic number of the atom has to be between 1 and 118
        centers: list of positions where this atom will be in the lattice
            numpy.array([xpos, ypos, zpos])
        bonds_to = list of Atoms that it's bonded to within a unit cell
        cell_position: position of the atom within the cell (3D cartesian)
        is_xyz: True if atom is made from an xyz file. Is not relevant
        currently but maybe in the future."""
        self.name = name
        self.radius = radius
        self.is_xyz = is_xyz and True  # this needs to be before other attrib
        self.number = z
        self.cell_position = None
        self.bonds_to = []
        self.centers = []

    def info(self):
        """Returns a information string for this atom"""
        info_str = "Name: {0}\nNumber:{3}\nRadius: {1}\nCell Position: {2}\n"
        info_str = info_str.format(
            self.name, self.number, self.radius, self.cell_position)
        bonds_str = "Bonds to: "
        for atom in self.bonds_to:
            bonds_str += atom.name + " "
        return info_str + bonds_str

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        """sets the name if it's a string. Make it unique if there are multiple
        atom of the same element in single script."""
        if type(new_name) is not str:
            raise ValueError(new_name + " is not a string")
        self._name = new_name

    @property
    def number(self):
        return self._z

    @number.setter
    def number(self, z):
        if z > 118 or z < 1:
            raise ValueError("{0} is not a valid atomic number.".format(z))
        self._z = z

    @property
    def centers(self):
        return self._centers

    @centers.setter
    def centers(self, pos_list):
        """sets the position list for this atom"""
        if not isinstance(pos_list, list):
            p_str = "{0} ".format(pos_list)
            raise ValueError(p_str + " in center list is not a valid list.")
        self._centers = pos_list

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, r):
        self._radius = r

    @property
    def cell_position(self):
        return self._cell_position

    @cell_position.setter
    def cell_position(self, vector):
        if (not self.is_xyz and vector is not None):
            assert(len(vector) == 3)
        self._cell_position = asarray(vector)

    @property
    def bonds_to(self):
        return self._bonds_to

    def will_bond_to(self, atom):
        """ Adds an atom to the bonds_to list of this atom"""
        assert(isinstance(atom, Atom))
        self._bonds_to.append(atom)

    # def add_centers(self, center):
    #    self._centers.append(center)

    @bonds_to.setter
    def bonds_to(self, bond_list):
        """sets the bonds_to list for this atom"""
        for atom in bond_list:
            if not isinstance(atom, Atom):
                raise ValueError(atom + " in bond list is not an Atom.")
        self._bonds_to = bond_list

    def get_atom_mesh(self, refinement=_REFINE, i=0):
        """ Return a mesh of atom at the ith center in centers
            refinement: resolution of the meshes, each int increment is 4X data
            created. ie object file created from refinement 4 will be 4X as big
            as that generated with refinement of 3.
            Return: mesh of this atom at a position given"""
        c = self._centers[i]
        return pm.generate_icosphere(self._radius, c, refinement)

    @classmethod
    def _get_pair(self, pointa, pointb):
        """ Private method returns a tuple of two points.
            The point with smaller norm first."""
        norm_a = linalg.norm(pointa)
        norm_b = linalg.norm(pointb)
        if norm_a < norm_b:
            return (pointa, pointb)
        return (pointb, pointa)

    def get_bond_pairs(self, unit_size, bonds_range=None):
        """ Gets all the bond_pairs for each pos for this atom using the
        bonds_to list. If the bonded atom is less than 1 unit cell away,
        pos of self and the bonded atom is added to the list.
        unit_size: length of a unit within the crystal
        returns: returns a list of tuples containing bond endpoints.
        Uses a brute force method to get the bonds. # of centers shouldn't be a
        large number.
        Note: If A is bonded to B, B does not need to have A in its bonds_to
        list. If it does, will have an uncessarily larger mesh output data."""
        bond_pairs = []
        for atom in self.bonds_to:  # for each bonded atom
            bonded_centers = atom.centers  # list of centers of bonded
            atom_cp = atom.cell_position
            #  distance to draw bonds at
            bond_dist = linalg.norm(self.cell_position - atom_cp)
            bond_dist = npround(bond_dist, self._ROUND)
            for center in self.centers:
                for bonded_center in bonded_centers:
                    dist = linalg.norm(center - bonded_center)
                    dist = npround(dist, self._ROUND)
                    if dist == bond_dist:
                        pair = self._get_pair(center, bonded_center)
                        bond_pairs.append(pair)
        return bond_pairs

    def get_simple_bond_pairs(self, unit_size, bonds_range=1):
        """ Gets simple bond_pairs that bond each atom to atom of its own kind
            unit_size: length of a unit within the crystal
            bond_range: can be 0, 1, or 2. See Lattice Params for details.
                      not used currently
            returns: bond_pairs to self atoms that are 1 unit cell away"""
        centers_list = list(self.centers)
        list_length = len(centers_list)
        bond_pairs = []
        for i in range(list_length):
            center = self.centers[i]
            for j in range(i+1, list_length):
                test_point = self.centers[j]
                dist = linalg.norm(center - test_point)
                dist = npround(dist, self._ROUND)
                if dist == unit_size:  # dist is a unit cell away
                    pair = self._get_pair(center, test_point)
                    bond_pairs.append(pair)
        return bond_pairs
