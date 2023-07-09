#runnning in Python 3
#----------Cubic Lattices----------#
#
#Algorithms for cubic crystal structures
#
#Author: Sayam Patel (sdpate11@ncsu.edu)
#Date: 5/16/2018
#
#----------------------------------#

import pymesh as pm
import numpy as np
import os

#in docker, this file is excecute as:
#exec(open("cubic.py").read())

"""Cubic Lattices:

Notes:
#i, j, k are 3 unit orthogonal vectors in 3D cartesian coordinates.
	Primitive vectors for lattices will be defined by a combination of these.
	As a sum, each crystal forms a unique Bravais Lattice made of of 3 primitive vectors
#each 3D lattice will have 3 primitive vectors that it will be defined by
#x, y, z will be the amount of discreet steps in each of the 3 primitive vectors of a
	lattice

1) Simple Cubic: 

	Atoms are at the 8 corners of a cube.
	The primitive lattice  R = a(i + j + k)
    where a is the distance between centers of two atoms on a primitive vector

2) Body Centered
	
	...

3) Face Centered
	
	...
  
"""

#--------Atom params--------#
radius = 1
refinement = 3

#--------Lattice params-----#

#note: each direction will have (step+1) number of atoms
#for eg a 1x1x1 lattice will be a cube with side size a and 8 atoms 
#they have to be positive whole numbers 0 to inf
xlen = 10  #steps in 1st primitive vector
ylen = 3 #steps in 2nd primitive vector
zlen = 4 #steps in 3rd primitive vector

dimensions = (xlen, ylen, zlen)

a = 5 #size of the cube
#---------------------------#

"""Creating simple cubic centers for spheres
	
dims: tuple 3 numbers defining the steps 
in each of the 3 primitive vectors for the Simple 
Cubic crystal. Going assume that it fills up from 
left to right. So for example if y = 0 and x=z=4, 
dims = (4, 4, 0)
	
a: size of the sides of this cube 
    
Returns a list of centers for atoms in simple cubic lattice
""" 
def simple_cubic_centers(dims, a):
	centers = []
	
	#create atoms in i direction with xsteps
	#xsteps + 1 so at least there is 1 atom at origin if xstep is 0
	for x in range(dims[0] + 1):
		centers.append(np.array([x*a, 0, 0]))
	
	
	#adding centers in the j direction on a line in i direction
	for x in range(dims[0] + 1):
		xcenter = centers[x]
		#start from 1 b/c 0 is already set by previous loop
		for y in range(1, dims[1] + 1) :
			xycenter = np.array(xcenter)
			xycenter[1] = y*a
			centers.append(xycenter)
	
	size_xygrid = (dims[0] + 1) * (dims[1] + 1) #size of grid
	
	#now taking each point on the ij grid and iterating in k direction	
	for xy in range(size_xygrid):
		xycenter = centers[xy]
		#start from 1 b/c 0 is already set by previous loops
		for z in range(1, dims[2] + 1):
			xyzcenter = np.array(xycenter)
			xyzcenter[2] = z*a
			centers.append(xyzcenter)
			
	return centers

"""Generates spheres and cylinders at the given locations and 
merges them together.
centers: list of centers to create crystals
r: radius of each sphere
refinement: resolution of sphere, 5 is fine for now
Returns a merged Mesh object of all the spheres
"""
def mesh_objects(centers, r, refinement):
	#this is fine for now but will have to see if 
	#it's better to copy and translate a sphere 
    #from origin and then merge	
	
	spheres = []
	
	for c in centers:
		spheres.append(pm.generate_icosphere(r, c, refinement))
	
	#making bonds from one atom to the next
	cylinders = []
	i = 0
	while i < len(centers):
		c = pm.generate_cylinder(centers[i], centers[i+1], cylinder_radius, cylinder_radius)
		i += 2
		cylinders.append(c)	
	
	return pm.merge_meshes(spheres + cylinders )
		

"""Brings everything together and outputs .obj file
"""
def main():
	centers = simple_cubic_centers(dimensions, a)
		
	mesh = mesh_objects(centers, radius, refinement)
	
	filename = "simple-cubic-lattice.obj"

	#docker image is mounted on the folder with this script and output dir
	relative_path = "output" #relative path from this current dir to output dir 
	current_dir = os.path.dirname(os.path.realpath('__file__'))
	output_file = os.path.join(current_dir, relative_path, filename)
	
	#saving mesh to a file
	pm.save_mesh(output_file, mesh, ascii=True)	
	
	print("Number of atoms:", len(centers))
	print("Number of verticies:", mesh.num_vertices)
	print("Size of file:", round(os.path.getsize(output_file) / 10.0**6, 2), "MB") 		
	print("Saved to:", output_file)

main()
