# running in Python 3
# ----------Lattice Class----------#
#
# Class for managing Atoms, their bonds, and meshes
#
# Author: Sayam Patel (sdpate11@ncsu.edu)
# Date created: 8/1/2018
#
# ----------------------------------#

from lib.crystal.atom import Atom
import pymesh as pm
import numpy as np
from math import sin, cos, radians, pi


class Lattice:
    """
    Crystal object
    """
    # default refinement/resolution for atom meshes
    _REFINE = Atom._REFINE
    # number of decimal points to round to (for vecs and dists)
    _ROUND = 3
    # origin numpy vector
    _ORIGIN = np.array([0, 0, 0])
    # factor by which to reduce the unit cell length to get bond radius
    _BOND_RADIUS = 30
    # list of unique atom in this crystal (used to count faces for color)
    _atom_models = []
    # list of bond meshes in this crystal
    _bond_models = []

    def __init__(self, name, unit_size, vectors, dims, is_xyz=False):
        """ Attributes
            name: name of crystal
            vectors: (r1, r2, r3) 3 primitive vectors define a bravais lattice
            dims: dimensions of the lattice (x, y, z)
            unit_size: size of a unit/basis cell in the crystal
            atoms: list of Atom objects in this crystal
            basis_vectors: vectors define the postion of atoms within unit cell
            is_xyz: True if lattice is made from an xyz file.
            #bond_pairs: list of endpoints of defining bonds
                        [(pt1, [pt2]), (p3, p4), (p5, p6), ... ]"""
        self.name = name
        self._updated = True
        self.is_xyz = is_xyz and True
        self._atoms = []
        self.dims = dims
        if not is_xyz:
            self.unit_size = unit_size
            self.vectors = vectors
        self.basis_vectors = [Lattice._ORIGIN]
        self._bond_pairs = []

    def info(self):
        """returns a string of information for this lattice"""
        info_str = "Name: {0}\nUnit Size: {1}\nDimensions: {2}\n"
        info_str = info_str.format(self.name, self.unit_size, self.dims)
        vector_str = "Vectors:\n"
        for v in self.vectors:
            vector_str = vector_str + "{0}\n".format(v)
        atoms_str = "Atoms:\n"
        for atom in self.atoms:
            atoms_str = atoms_str + "\n" + atom.info() + "\n"
        return info_str + vector_str + atoms_str

    @classmethod
    def from_angles(cls, name, unit_size, angles, dims):
        """Calculates 3 bravais lattice from 3 given angles between them"""
        msg = "Angles need to be in a tuple of 3 angles (In degrees)."
        cls._check_tuple(angles, msg)
        alpha = radians(angles[0])  # angle between r2 and r3
        beta = radians(angles[1])  # angle between r1 and r3
        gamma = radians(angles[2])  # angle between r1 and r2
        theta = pi/2 - alpha  # angle to rotate pre_r3
        # making vectors
        r1 = np.array([1, 0, 0])  # 'starting' vector in x hat
        # rotate xhat by gamma to get r2
        r2 = np.array([cos(gamma), sin(gamma), 0])
        # rotating xhat by beta to partial r3
        pre_r3 = np.array([cos(beta), 0, sin(beta)])
        # getting the two components of r3 after rotating pre_r3
        r3_component1 = sin(theta) * r2
        r3_component2 = cos(theta) * pre_r3
        r3 = r3_component1 + r3_component2
        # normalizing r3
        r3_norm = np.linalg.norm(r3)
        if r3_norm != 0:
            r3 = r3 / r3_norm

        vectors = (r1, r2, r3)
        vectors = [np.round(v, cls._ROUND) for v in vectors]
        return cls(name, unit_size, vectors, dims)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        """sets the name if it's a string"""
        if type(new_name) is not str:
            raise ValueError(new_name + " is not a string")
        self._name = new_name

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atom_list):
        """Sets the atoms list for this lattice. Checks correct input list
        element type (Atom). Asserts that number of atoms being set is equal
        to the number of basis vectors for this lattice. Then assigns those
        basis_vectors (as cell_positions) to the atoms. Assigment is done in
        order that the basis vectors/atoms are added ie. first atom in the
        input list gets the first (origin) vector."""
        if not self.is_xyz:
            try:
                assert(len(self.basis_vectors) > 0)
            except AssertionError:
                print("No basis vectors found. Cannot add atoms.")
                raise
            try:
                assert(len(self.basis_vectors) == len(atom_list))
            except AssertionError:
                print("# of basis vectors not equal to # of atoms being added.")
                raise
            for i in range(len(atom_list)):
                atom = atom_list[i]
                if not isinstance(atom, Atom):
                    raise ValueError(atom + " in atom list is not an Atom.")
                    # assign cell positions to atoms
                atom.cell_position = self.basis_vectors[i] * self.unit_size
        else:  # when atoms are from file
            for i in range(len(atom_list)):
                atom = atom_list[i]
                if not isinstance(atom, Atom):
                    raise ValueError(atom + " in atom list is not an Atom.")
        # checking for unique names
        atom_names = []
        for atom in atom_list:
            atom_names.append(atom.name)
        for name in atom_names:
            if atom_names.count(name) > 1:
                msg = name + """ is repeated more than once. Make sure atoms in
                 a basis have unique names. Don't use ","s or "-"s in names."""
                raise ValueError(msg)
        self._atoms = atom_list
        self._updated = True

    @property
    def unit_size(self):
        return self._unit_size

    @unit_size.setter
    def unit_size(self, unit_size):
        """sets unit size"""
        self._unit_size = unit_size
        self.updated = True

    @property
    def dims(self):
        return self._dims

    @dims.setter
    def dims(self, dims):
        """sets the dims if given a tuple of 3 numbers"""
        self._check_tuple(dims, "Needs to be: (length, width, height)")
        self._dims = dims
        self._updated = True

    @property
    def vectors(self):
        return (self._r1, self._r2, self._r3)

    @vectors.setter
    def vectors(self, vectors):
        """sets the lattice vectors if a tuple of 3 tuples of 3 numbers"""
        for v in vectors:
            self._check_tuple(v, "Tuple needs to be a vector in cartesian.")
        self._r1 = np.asarray(vectors[0])
        self._r2 = np.asarray(vectors[1])
        self._r3 = np.asarray(vectors[2])
        self._updated = True

    @property
    def basis_vectors(self):
        return self._basis_vectors

    @basis_vectors.setter
    def basis_vectors(self, basis_vectors):
        """sets the basis vector for atoms in a unit cell"""
        self._basis_vectors = []
        for bv in basis_vectors:
            msg = "Basis vector tuple needs to be a vector in cartesian."
            self._check_tuple(bv, msg)
            self._basis_vectors.append(np.asarray(bv))
            self._updated = True

    @property
    def bond_pairs(self):
        return self._bond_pairs

    @bond_pairs.setter
    def bond_pairs(self, bond_pairs):
        self._bond_pairs = bond_pairs

    def add_atom(self, atom):
        if not isinstance(atom, Atom):
            raise ValueError(atom + " in atom list is not an Atom.")
        self._atoms.append(atom)
        self.updated = True

    @classmethod
    def _check_tuple(self, alleged_tuple, msg=""):
        """Checks if the input is a tuple of size 3. Raises appropriate
            exceptions if not valid tuple along with err msg"""
        is_not_tuple = type(alleged_tuple) is not tuple
        is_not_np_array = type(alleged_tuple) is not np.ndarray
        if is_not_tuple and is_not_np_array:
            t_str = "{0} ".format(alleged_tuple)
            raise ValueError(t_str + "is not a tuple or numpy array " + msg)
        elif len(alleged_tuple) != 3:
            t_str = "{0} ".format(alleged_tuple)
            raise ValueError(t_str + "doesn't have 3 entries.")
        return None

    def get_lattice_mesh(self, refinement=_REFINE):
        """ Return a complete mesh of this lattice. Uses assign_points_to_atoms
            refinement: resolution of the meshes, each int increment is 4X data
            created. ie object file created from refinement 4 will be 4X as big
            as that generated with refinement of 3. Looks for a name if given,
            if not given, uses the index value, default index value is 0
            Return: mesh of this crystal's lattice with bonds
        """
        if self._updated:  # change was made to the atoms/dims/vectors
            # assigning positions to atoms
            self.assign_points_to_atoms()
            # getting all the new bonds
            self.bond_pairs = []
            if (len(self.atoms) == 1):
                # when there is only 1 atom in lattice gets simple bonds
                atom = self.atoms[0]
                self.bond_pairs += atom.get_simple_bond_pairs(self.unit_size)
            else:
                for atom in self.atoms:
                    self.bond_pairs += atom.get_bond_pairs(self.unit_size)
            self._updated = False
        # getting the spheres
        spheres = []  # [atom 1 meshes, atom 2 meshes, ...]
        cylinders = []  # bond meshes
        cr = self.unit_size / Lattice._BOND_RADIUS
        for i in range(len(self.atoms)):
            spheres.append(self.get_simple_mesh(refinement, i))
        for bp in self._bond_pairs:
            cylinders.append(pm.generate_cylinder(bp[0], bp[1], cr, cr))

        return pm.merge_meshes(spheres+cylinders)

    def get_lattice_points(self):
        """ Generates centers for spheres (atoms) from this lattice's 3
            primitive bravais vectors in 3D cartesian coordinates.
            r1,r2,r3: tuple of 1st, 2nd, 3rd primitive vector
            r1,r2,r3steps: steps to take in r1, r2, r3
            Returns a list of points in the lattice
            """
        r1steps = self.dims[0]
        r2steps = self.dims[1]
        r3steps = self.dims[2]
        # getting the complete vectors
        r1 = self._r1 * self.unit_size
        r2 = self._r2 * self.unit_size
        r3 = self._r3 * self.unit_size
        # center lists
        r1centers = []
        r2centers = []
        r3centers = []

        # adding points in dir of 1st primitive vec (making a line)
        for i in range(r1steps+1):
            new_center = np.array([])
            if i == 0:
                new_center = np.array([0, 0, 0])  # initial center at origin
            else:
                new_center = r1centers[-1] + r1
            r1centers.append(new_center)

            # add points in r3 direction from line of points crossing origin
            for k in range(r3steps):
                if k == 0:
                    new_center = r1centers[-1] + r3
                else:
                    new_center = r3centers[-1] + r3
                r3centers.append(new_center)
            # adding points in dir of 2nd primitive vec (making a 2d grid)
            for j in range(r2steps):
                if j == 0:
                    new_center = r1centers[-1] + r2
                else:
                    new_center = r2centers[-1] + r2
                r2centers.append(new_center)

                # adding points in r3 direction to the points in the grid
                # (making a 3d lattice)
                for k in range(r3steps):
                    if k == 0:
                        new_center = r2centers[-1] + r3
                    else:
                        new_center = r3centers[-1] + r3
                    r3centers.append(new_center)
        lattice_points = r1centers + r2centers + r3centers
        return lattice_points

    def assign_points_to_atoms(self):
        """
        Uses get_lattice_points to get the basic points, then assigns each
        atom in the atoms list its positions calculated using BL and basic
        vectors.
        Each i'th atom's postion is calulated as:
        n*r1 + m*r2 + l*r3 + basis_vec[i]
        Where n, m, and l are integer steps <= r1,r2,r3steps but >= 0.
        """
        # same number of atoms as positioned by basis_vectors
        assert(len(self.atoms) == len(self.basis_vectors))
        assert(len(self.atoms) > 0)

        lattice_points = self.get_lattice_points()
        for i in range(len(self.atoms)):
            # calculating positions: lattice pt vec + basic vec
            points_to_add = []
            for pt in lattice_points:
                points_to_add.append(pt + self.basis_vectors[i]*self.unit_size)
            # assigning the points to atom
            self.atoms[i].centers = points_to_add

    def get_simple_mesh(self, refinement=_REFINE, index=0, name_of_atom=''):
        """ Return a mesh of an atom in this crystal at all centers of that
            atom.
            refinement: resolution of the meshes, each int increment is 4X
            data created. ie object file created from refinement 4 will be 4X
            as big as that generated with refinement of 3. Looks for a name if
            given, if not given, uses the index value, default index value is 0
            i: ith atom in the atoms list for this crystal
            name_of_atom: name of atom to get mesh of
            Return: mesh of atom at positions defined by atom centers, no bonds
            """
        atom = None
        if len(name_of_atom) > 0:
            for atom_object in self.atoms:
                if name_of_atom == atom_object.name():
                    atom = atom_object
                    break
            if atom is None:
                msg = " not found in lattice."
                raise ValueError("Name: " + name_of_atom + msg)
        else:
            atom = self.atoms[index]
        # getting all the positions to draw the atom at
        centers = atom.centers
        spheres = []
        print("Making atom mesh for:", atom.name)
        for i in range(len(centers)):
            spheres.append(atom.get_atom_mesh(refinement, i))
        return pm.merge_meshes(spheres)
