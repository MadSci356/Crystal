# running in Python 3
# ----------LatticeIO Class----------#
#
# Object for lattice IO
#
# Author: Sayam Patel (sdpate11@ncsu.edu)
# Date created: 8/6/2018
#
# ----------------------------------#
import pymesh as pm
import numpy as np
from os import path, listdir, stat
from shutil import copyfile
from lib.crystal.lattice import Lattice
from lib.crystal.atom import Atom


class LatticeIO:

    # default refinement/resolution for atom meshes
    _REFINE = Atom._REFINE
    # number of decimal points to round to (for dims)
    _ROUND = Lattice._ROUND
    # default radius of atoms
    _RADIUS = 1
    # factor that scales the distance between atoms from xyz file
    _DIST_FACTOR = 1.2
    # list of colors for lattice, last color is for bond always
    _COLORS = ("red", "green", "blue", "yellow", "orange",
               "purple", "magenta")
    _BOND_COLOR = ["bondgray"]
    # material file name
    _MTL_FILE = "lattice-colors.mtl"

    # material file source
    _MTL_DIR = "./data/" # path from crystal directory to dir of mtl file

    def __init__(self, lattice, bond_pattern):
        """Attributes:
            lattice: Lattice object for creating mesh
            bond_pattern: Bond Pattern to make bonds in atoms in this lattice
            See assign_bonds method below for more details."""
        self.lattice = lattice
        if (not lattice.is_xyz):
            self._assign_bonds(bond_pattern)

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, new_lattice):
        """Sets lattice for LatticeIO if correct instance."""
        if (isinstance(new_lattice, Lattice)):
            self._lattice = new_lattice
        else:
            msg = "{0} is not a Lattice object.".format(new_lattice)
            raise ValueError(msg)

    @classmethod
    def _check_name(self, atoms, name, radius):
        """Checks if name is in the list of atoms given. If so, returns that
            atom. If not, makes a new atom object and gives it a number
            that is len(atom) before adding the new atom. So a new atom will
            have number 0. Second new atom will have 2...
            atoms: list of atoms to check the name against
            name: name of atom to check in atoms list
            return: Atom object with same name as name input if found. If not
            found, will make a new object and will add it to atoms list."""
        for atom in atoms:
            if atom.name == name:
                return atom
        new_atom = Atom(name, radius, len(atoms) + 1, is_xyz=True)
        atoms.append(new_atom)
        return new_atom

    @classmethod
    def from_file(cls, xyzfile, radius=_RADIUS):
        if not (xyzfile.endswith(".XYZ")):
            msg = "{0} is not a valid file. Make sure it's a .xyz file."
            raise ValueError(msg.format(xyzfile))
        f = open(xyzfile, "r")
        lines = f.readlines()
        atoms = []  # atom in this file
        for line in lines[2:]:
            line_data = line.split()
            name = line_data[0]
            # make sure element symbol is capitalized and 2nd letter is lower
            name = name.capitalize()
            if len(name) == 2 and name.isalpha():
                name = name[0] + name[1].lower() # in case 2nd letter is upper

            position = []  # extracting position data
            for x in line_data[1:4]:
                position.append(float(x) * cls._DIST_FACTOR)
            center = np.array(position)  # making it np array
            atom = cls._check_name(atoms, name, radius)
            atom.centers.append(center)
        xyz_name = lines[1].split()[0]  # name of file or crystal
        num_atoms = int(lines[0])  # num of atoms
        dim = round(num_atoms**(1/3), cls._ROUND)  # assuming it's N^3 atoms
        dims = (dim, dim, dim)
        xyz_lattice = Lattice(xyz_name, None, None, dims, is_xyz=True)
        xyz_lattice.atoms = atoms
        return LatticeIO(xyz_lattice, "")

    def make_obj_file(self, filename, mtl=True):
        """Uses _get_lattice_mesh to create mesh from the atoms in the lattice
        and save in a folder called output. Prints some lattice and .obj
        file stats"""

        # docker image is mounted on the folder with script and output dir
        relative_path = "output/"  # relative path of main.py to output dir
        current_dir = path.dirname(path.realpath('__file__'))
        output_file = path.join(current_dir, relative_path, filename)

        # culmination of my work in the next line
        lattice_mesh = None
        if not self.lattice.is_xyz:
            lattice_mesh = self.lattice.get_lattice_mesh()
        else:  # if from xyz file
            atom_meshes = []
            for i in range(len(self.lattice.atoms)):
                atom_meshes.append(self.lattice.get_simple_mesh(index=i))
            lattice_mesh = pm.merge_meshes(atom_meshes)
            #lattice_mesh = pm.remove_duplicated_faces(lattice_mesh)

        # saving (colorless) mesh to a file
        pm.save_mesh(output_file, lattice_mesh, ascii=True)

        dims = self.lattice.dims
        num_points = len(self.lattice.atoms[0].centers)

        num_atoms = len(self.lattice.atoms)
        if not self.lattice.is_xyz:
            print("Dimensions of the lattice:",
              dims[0], "X", dims[1], "X", dims[2])
        else:
            print("(approx) Dimensions of the lattice:",
              dims[0], "X", dims[1], "X", dims[2])
        print("Number of lattice points:", num_points)
        print("Number of unique atoms:", num_atoms)
        print("Number of verticies:", lattice_mesh.num_vertices)
        print("Size of file:", round(path.getsize(
                output_file) / 10.0**6, 2), "MB")
        print("Saved to:", output_file)

        # adding mtl
        if mtl:
            self._addmtls(output_file, filename, lattice_mesh)

    def _assign_bonds(self, bond_pattern):
        """Method to assign bonds to atoms in this lattice.
        bond_pattern: string pattern defining atom bond pairs within the basis
        of this lattice. Pattern: "A-B, C-D, E-F, ..., X-Y" where the capital
        letters are names of the atoms in the lattice. The names MUST match
        those of the atoms. Don't use hypens (-) or commas (,) names. Anything
        you can't name a regular python variable, go ahead and apply those
        no-nos to naming your atoms too. A space should follow every comma.
        No comma after the last bond-pair."""
        try:
            assert(len(self.lattice.atoms) > 1)
        except AssertionError:
            msg = """{0} has less than or equal to 1 atom. Bond pattern
             will not be used. Simple bonds will form. Add atoms to
             Lattice object if bond pattern is to be used."""
            msg = msg.format(self.lattice.name)
            print(msg)
            return
        # basic checking for valid input, further Errors can occur later
        num_comma = bond_pattern.count(", ")
        num_hypens = bond_pattern.count("-")
        try:
            assert(num_comma + 1 == num_hypens)
        except AssertionError:
            print("',' and '-' pattern inconsistent, check bond pattern.")
            exit(1)
        # bond pattern processing
        bond_strs = bond_pattern.split(", ")  # ['A-B', 'C-D', ...']
        for bond_str in bond_strs:
            atom_names = bond_str.split("-")  # ['A', 'B']
            atom1_name = atom_names[0]
            atom2_name = atom_names[1]
            atom1, atom2 = None, None
            for atom in self.lattice.atoms:
                if atom.name == atom1_name:
                    atom1 = atom
                if atom.name == atom2_name:
                    atom2 = atom
            if atom1 is None or atom2 is None:
                msg = "Bond pattern had a name that doesn't match that of" + \
                    " any atom in {0}. Please make sure bond pattern has" + \
                    " name matching those of the atoms added to the lattice."
                msg = msg.format(self.lattice.name)
                raise ValueError(msg)
            atom1.will_bond_to(atom2)  # bonding atoms

    @classmethod
    def _get_mtl_list(cls, color_strs):
        """Given a list of colors, generates a list of grouping and usemtl
        lines to be inserted into the obj so it can use the material file
        See colors section for material file details.
        Example Input is ["red", "green"]
        Output: [\ng red\nusemtl red\n\n, \ng green\nusemtl green\n\n]
        When a mtl string into the obj file, it will look like:
        >...obj file data
        >
        >g red
        >usemtl red
        >
        >...obj file data
        Returns: List of mtl group and usemtl string for each color"""
        mtl_list = []
        for color in color_strs:
            mtlstr = "\ng " + color + "\n" + "usemtl " + color + "\n\n"
            mtl_list.append(mtlstr)
        return mtl_list

    def _get_num_faces(self):
        """Returns a list of numbers indicating the number of faces in the mesh
        of in each of the atoms that are in the lattice."""
        num_faces = []
        for atom in self.lattice.atoms:
            atom_mesh = atom.get_atom_mesh()
            atom_mesh.add_attribute("face_index")
            faces = len(atom_mesh.get_attribute("face_index"))
            num_faces.append(faces)
        return num_faces

    def _addmtls(self, objfile_path, filename, lattice_mesh):
        """Adds mtls to obj_file from the mtl_list. Mtl names are in the
        _COLORS class variable"""

        print("Adding Materials to atoms")

        # getting all the things ready for writing
        # which mtl file to use at top
        mtl_header = "\n# mtls added by crystal"
        mtl_header += "\nmtllib {0} \n\n".format(filename[:-4] + ".mtl" )
        mtl_list = self._get_mtl_list(LatticeIO._COLORS)

        # same elements will have same color if <= # of colors available
        atomic_nums = [atom.number for atom in self.lattice.atoms]
        unique_nums = []
        unique_nums = [num for num in atomic_nums if num not in unique_nums]
        if len(unique_nums) > len(LatticeIO._COLORS):
            print("""Warning: There are more unique atoms in this lattice than
                   there are unique colors available. Unique atoms will have
                   same colors.""")

        # faces to be colored in objfile
        # 1 bc 1st line comment in objfile
        faces_start_index = lattice_mesh.num_vertices + 1
        faces_per_atom = self._get_num_faces()
        dims = self.lattice.dims
        num_atoms = (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1)

        # editing the file
        f = open(objfile_path, "r+")
        lines = f.readlines()
        lines[0] = lines[0] + mtl_header

        # adding mtl for faces of atoms
        # faces are ordered by atoms. So faces for atom1 then atom 2...atom3...
        face_to_edit = faces_start_index
        num_colors = len(LatticeIO._COLORS)
        for i in range(len(self.lattice.atoms)):
            atom = self.lattice.atoms[i]
            # modulus incase more atom types than colors available
            color_index = unique_nums.index(atom.number) % num_colors
            newstr = mtl_list[color_index]
            lines[face_to_edit] = newstr + lines[face_to_edit]
            if (self.lattice.is_xyz):
                face_to_edit += faces_per_atom[i]*(len(atom.centers))
                print(len(atom.centers), faces_per_atom[i])
            else:
                face_to_edit += faces_per_atom[i]*num_atoms

        if (not self.lattice.is_xyz):
            print("Adding bond mtl...")
            # adding a separate bond mtl after stoms are done
            bondmtl = self._get_mtl_list(LatticeIO._BOND_COLOR)[0]  # bond mtl

            # face_to_edit should be at the first face of the bonds
            lines[face_to_edit] = bondmtl + lines[face_to_edit]

        lines = "".join(lines)
        f.seek(0)
        f.write(lines)
        f.close()

        # copying mtl file
        output_dir = path.dirname(objfile_path)
        # output_dir += "/" + self._MTL_FILE
        output_dir += "/" + filename[:-4] + ".mtl"

        # mtl_file_exists = False
        # for f in listdir(output_dir):
        #     if f.endswith(".mtl"):
        #         mtl_file_exists = True
        # if not mtl_file_exists:
        copyfile(self._MTL_DIR + self._MTL_FILE, output_dir)
        #mode = stat(output_dir).filemode('-rw-rw-rw-')



        print("Done.")
