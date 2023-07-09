#  running in Python 3
#  ----------Bravais Lattices----------#
#
#  Algorithms for creating lattice structure from
#  bravais lattice values
#
#  Author: Sayam Patel (sdpate11@ncsu.edu)
#  Date: 5/16/2018
#
#  ----------------------------------#

#  in docker, this file is excecuted as:
#  exec(open("name-of-file.py").read())

import pymesh as pm
import numpy as np
import os
#  --------------------Lattice functions--------------------#


def getpair(pointa, pointb):
    """Returns a tuple of two points.
    The point with smaller norm first"""
    norm_a = np.linalg.norm(pointa)
    norm_b = np.linalg.norm(pointb)
    if norm_a < norm_b:
        return (pointa, pointb)
    return (pointb, pointa)


def getbondpairs(allcenters, a, basis_vectors):
    """ Given a list of centers and unit distance between
    returns a list of tuples containing bond endpoints.
    Uses a brute force method to get the bonds. #  of centers
    shouldn't be a large number.
    allcenters: list of all center points in the lattice
    a: unit of distance between lattice points. See Lattice Params.
    bondtype: can be 0, 1, or 2. See Lattice Params for details.
    basis_vectors: list of vectors defining atoms in a unit cell
    Returns: List of tuples signifying bonds:
    [(pointa, pointb), (pointc, pointd),...]"""
    basis_lengths = []
    for vector in basis_vectors:
        vec_len = np.linalg.norm(vector)
        basis_lengths.append(vec_len*a)
    bondpairs = []
    list_length = len(allcenters)
    for i in range(list_length):
        center = allcenters[i]
        for j in range(i+1, list_length):
            test_point = allcenters[j]
            dist = np.linalg.norm(center - test_point)
            if dist in basis_lengths:  # dist is within a unit cell
                pair = getpair(center, test_point)
                bondpairs.append(pair)
    return bondpairs


def bravais_centers(r1, r2, r3, r1steps, r2steps, r3steps, a):
    """Generates the centers for the spheres (atoms) from the given
        3 primitive bravais lattice vectors 3D cartesian coordinates.
        r1: numpy array of 1st primitive vector
        r2: numpy array of 2nd primitive vector
        r3: numpy array of 3rd primitive vector
        r1steps: steps to take in r1
        r2steps: steps to take in r2
        r3steps: steps to take in r3
        a: length/space between atoms
        Returns a list of points in the lattice"""
    #  getting the complete vectors
    r1 = r1*a
    r2 = r2*a
    r3 = r3*a
    #  center lists
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
            #  bondpairs.append((r1centers[-1], new_center))
        r1centers.append(new_center)

        #  adding points in r3 direction from line of points crossing origin
        for k in range(r3steps):
            if k == 0:
                new_center = r1centers[-1] + r3
                #  bond1 = r1centers[-1]
            else:
                new_center = r3centers[-1] + r3
                # bond1 = r3centers[-1]
            # bondpairs.append((bond1, new_center))
            r3centers.append(new_center)
        # adding points in dir of 2nd primitive vec (making a 2d grid)
        for j in range(r2steps):
            if j == 0:
                new_center = r1centers[-1] + r2
                # bond1 = r1centers[-1]
            else:
                new_center = r2centers[-1] + r2
                # bond1 = r2centers[-1]
            # bondpairs.append((bond1, new_center))
            r2centers.append(new_center)

            # adding points in r3 direction to the points in the grid
            # (making a 3d lattice)
            for k in range(r3steps):
                if k == 0:
                    new_center = r2centers[-1] + r3
                    # bond1 = r2centers[-1]
                else:
                    new_center = r3centers[-1] + r3
                    # bond1 = r3centers[-1]
                # bondpairs.append((bond1, new_center))
                r3centers.append(new_center)
    lattice_points = r1centers + r2centers + r3centers
    basis_centers = []
    # adding the additional vector points defined by basisvectors
    for point in lattice_points:
        for vector in basis_vectors:
            basis_centers.append(point + a*vector)
    allcenters = lattice_points + basis_centers
    return allcenters

# --------------------Mesh functions--------------------#


"""Generates lattice atoms an bonds at the given locations
    and merges them. Adds one of each type of spheres
    to the global mesh_atoms and mesh_bonds list (used to add
    materials)
    centers: list of centers to create spheres at
    bondpairs: a list of 2 element tuples containing
    two points for creating bond cylinders
    r: radius of each sphere
    refinement: resolution of sphere, 5 is fine for now
    Returns a merged Mesh object of lattice
"""


def mesh_lattice(centers, bondpairs, r, refinement):
    # this is fine for now but will have to see if
    # it's better to copy and translate a sphere
    # from origin and then merge
    cr = r/8  # cylinder radius

    spheres = []
    for c in centers:
        spheres.append(pm.generate_icosphere(r, c, refinement))
    cylinders = []
    for bp in bondpairs:
        cylinders.append(pm.generate_cylinder(bp[0], bp[1], cr, cr))

    # adding samples to model lists
    atom_models.append(spheres[0])
    if len(cylinders) > 0:
        bond_models.append(cylinders[0])

    return pm.merge_meshes(spheres + cylinders)


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
    Returns: List of mtl group and usemtl string for each color
"""


def get_mtl_list(colorlist):
    mtl_list = []
    for color in colorlist:
        mtlstr = "\ng " + color + "\n" + "usemtl " + color + "\n\n"
        mtl_list.append(mtlstr)
    return mtl_list


def get_num_faces(models):
    """Returns a list of numbers indicating the number of faces
        in each of the meshes in the input list
        models: list of meshes
        Returns: a list of number of faces in each mesh in models
        """
    num_faces = []
    for atom in atom_models:
        atom.add_attribute("face_index")
        faces = len(atom.get_attribute("face_index"))
        num_faces.append(faces)
    return num_faces


"""Adds mtls to obj_file from the mtl_list in a specific pattern"""


def addmtls(obj_file, colorlist, lattice_mesh, pattern):
    print("Adding Materials...")

    # faces to be colored in objfile
    # 1 bc 1st line comment in objfile
    faces_start_index = lattice_mesh.num_vertices + 1
    num_faces_per_atom = get_num_faces(atom_models)[0]  # 0 bc only 1 atom
    # num_faces_per_bond = get_num_faces(bond_models)[0]
    # 0 bc only 1 type of bond
    num_atoms = (r1steps + 1) * (r2steps + 1) * (r3steps + 1)

    # editing file
    f = open(obj_file, "r+")
    lines = f.readlines()
    # writing which mtl file to use at top
    mtl_header = "\n# mtls added by crystal\nmtllib " + mtlfilename + "\n\n"
    mtl_list = get_mtl_list(colorlist)

    print("Materials:", len(mtl_list))

    lines[0] = lines[0] + mtl_header
    face_to_edit = -1
    for i in range(num_atoms):
        # pattern input can be used to color lattice, basic pattern for now
        face_to_edit = faces_start_index + i*num_faces_per_atom
        is_even = int(i % 2 == 0)
        newstr = mtl_list[is_even]
        lines[face_to_edit] = newstr + lines[face_to_edit]

    # adding a separate bond mtl after stoms are done
    bondmtl = mtl_list[-1]  # bond material
    bonds_index_start = face_to_edit + num_faces_per_atom
    lines[bonds_index_start] = bondmtl + lines[bonds_index_start]

    lines = "".join(lines)
    f.seek(0)
    f.write(lines)
    f.close()
# --------------Atom params----------------#


radius = 1  # radius of atoms
refinement = 3  # resolution for generating the models

# ------------Lattice params---------------#

"""Notes on steps:
- Each direction will have (step+1) number of atoms
- For eg a 1x1x1 lattice will be a cube with side size a and 8 atoms
- They have to be positive whole numbers 0 to inf
- Direction below is that of the respective primitive vectors"""
r1steps = 1  # steps in first dir
r2steps = 1  # steps in second dir
r3steps = 0  # steps in third dir

'''
can be 0, 1, or 2.
- 0: only finds bonds less than 1 unit distance
- 1: only finds bonds at 1 unit distance
- 2: finds bonds <= 1 unit
'''
bondtype = 0

# primitive vector below in [x_hat, y_hat, z_hat]
# primitive vectors are unit vectors
'''
r1 = np.array([0, 0.5, 0.5])
r2 = np.array([0.5, 0, 0.5])
r3 = np.array([0.5, 0.5, 0])
'''
sqrt2 = 2**.5
r1 = np.array([1, 0, 0])
r2 = np.array([sqrt2/2, sqrt2/2, 0])
r3 = np.array([0, 0, 0])

'''List of vectors that define atom positions relative to a lattice point
calculated from the Bravais vectors. Basis vectors will be < 1 unit away.
note: there is one 'standard' atom at each of the lattice point. Basis vectors
are additonal vectors in a unit cells.
'''

basis_vectors = [np.array([1/2, 1/6, 0])]  # ~20 deg from x axis

a = 8  # length/space between atoms

# ----------------Colors--------------------#
"""
in the list below will be colors that are written in
the lattice-colors.mtl file. make sure the names match
only put colors here that are in the file
"""

colorlist = ["red", "green", "bondcolor"]  # , "blue"]
# material file name. generally it's in the same
# directory as the obj file
mtlfilename = "lattice-colors.mtl"
# list of different atom meshes in the lattice
# a list because in the future there might be different
# types of atoms/bonds
atom_models = []
bond_models = []  # hehe
# ----------End Parameters------------------#

"""
Brings everything together and outputs .obj file
"""


def main():
    # get centers and bondpairs from constants
    centers = bravais_centers(r1, r2, r3, r1steps, r2steps, r3steps, a)
    bondpairs = getbondpairs(centers, a, basis_vectors)
    # create mesh
    lattice_mesh = mesh_lattice(centers, bondpairs, radius, refinement)

    filename = "testing/basistest/graphitetest.obj"

    # docker image is mounted on the folder with this script and output dir
    relative_path = "output/"  # relative path from current dir to output dir
    current_dir = os.path.dirname(os.path.realpath('__file__'))
    output_file = os.path.join(current_dir, relative_path, filename)

    # saving (colorless) mesh to a file
    pm.save_mesh(output_file, lattice_mesh, ascii=True)

    print("Dimensions of the lattice:", r1steps, "X", r2steps, "X", r3steps)
    print("Number of atoms in lattice:", len(centers))
    print("Number of verticies:", lattice_mesh.num_vertices)
    size = round(os.path.getsize(output_file) / 10.0**6, 2)
    print("Size of file:", size, "MB")
    print("Saved to:", output_file)

    # addmtls(output_file, colorlist, lattice_mesh, pattern=0)


main()
