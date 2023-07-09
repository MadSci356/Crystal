# import pymesh as pm
# import numpy as np
# import os
# from math import cos, sin, pi, radians
#
# #in docker, this file is excecute as:
# #exec(open("crystal.py").read())
#
# #setting the radius and resolution
# radius = 1
# cylinder_radius = radius/5
# refinement = 3	 #refinement order (graphics resolution)
#
# #centers
# c0 = np.zeros(3)
# c1 = c0 + 3
# c2 = c1 + 3
# c3 = c2 + 3
# #c2 = c0 #np.array([3,0,0])
# centers = [c0, c1, c2, c3] #list of centers
#
# #whatever folder you may want to put it in, must exist before writing to it
# filename = "testing/basistest/graphitetest.obj"
#
# #docker image is mounted on the folder with this script and output dir
# relative_path = "output/" #relative path from this script dir to output dir
# current_dir = os.path.dirname(os.path.realpath('__file__'))
# print(current_dir)
# output_file = os.path.join(current_dir, relative_path, filename)
# print(output_file)
#
# #sphere = pm.generate_icosphere(radius, c0, refinement) #sphere at origin
# #sphere1 = pm.generate_icosphere(radius, c1, refinement) #sphere at 3,3,3
#
# spheres = [] #list of spheres
# cylinders = []
# for center in centers:
# 	#creating a mesh of sphere at a given center and adding to a list
# 	s = pm.generate_icosphere(radius, center, refinement)
# 	#s.add_attribute("color");
# 	#s.add_attribute("face_index")
# 	#print(s.get_attribute("face_index"))
# 	#print('')
# 	spheres.append(s)
#
#
#
# #making bonds from one atom to the next
# '''i = 0
# while i < len(centers):
# 	c = pm.generate_cylinder(centers[i], centers[i+1], cylinder_radius, cylinder_radius)
# 	i += 2
# 	cylinders.append(c)'''
# '''for i in range(len(centers)):
# 	print ("index:", i)
# 	c = pm.generate_cylinder(centers[i], centers[i+1], cylinder_radius, cylinder_radius)
# 	i += 2
# 	cylinders.append(cylinders)'''
#
# #combining all the sphere and cylinders meshes into a single mesh
# mesh = pm.merge_meshes(spheres + cylinders)
#
# #saving mesh to a file
# pm.save_mesh(output_file, mesh, ascii=True)
#
# #adding material
# color1 = "green"
# color2 = "red"
# color3 = "blue"
# mtlfilename = "lattice-colors.mtl"
# mtlstr1 = "\ng " + color1 + "\n" + "usemtl " + color1 + "\n\n"
# mtlstr2 = "\ng " + color2 + "\n" + "usemtl " + color2 + "\n\n"
# mtl_list = [mtlstr1 , mtlstr2]
#
# # f = open(output_file, "r") #opening file
# # contents = f.readlines()
# # f.close()
# # #writing which mtl file to use
# # contents[0] = contents[0] + "\n# mtls added by Sayam Patel\nmtllib " + mtlfilename + "\n\n"
# # spheres[0].add_attribute("face_index")
# # num_faces = len(spheres[0].get_attribute("face_index"))
# #
# # faces_start_index = mesh.num_vertices + 1 # +1 bc first line is a comment added by PyMesh
# #
# # #making groups and adding mtl lines
# # for i in range(len(spheres)):
# # 	index = faces_start_index + i*num_faces
# # 	if i % 2 == 0:
# # 		contents[index] = mtl_list[0] + contents[index]
# # 	else:
# # 		contents[index] = mtl_list[1] + contents[index]
# #
# # contents = "".join(contents)
# # f = open(output_file, "w")
# # f.write(contents)
# # f.close()
#
# # print(mesh.get_attribute_names())
# # print("Saved to:", output_file)
# # print("Verticies: ", mesh.num_vertices)
# # print("Materials: ", len(mtl_list))
# # print("Size of file:", round(os.path.getsize(output_file) / 10.0**6, 2), "MB")
#
# def calculate_vectors(angles):
#     """Calculates 3 bravais lattice from 3 given angles between them"""
#     msg = "Angles need to be in a tuple of 3 angles (In degrees)."
#     #_check_for_tuple_and_size(angles, msg)
#     alpha = radians(angles[0]) #angle between r2 and r3
#     beta = radians(angles[1]) #angle between r1 and r3
#     gamma = radians(angles[2]) #angle between r1 and r2
#     theta = pi/2 - alpha #angle to rotate pre_r3
#     #making vectors
#     r1 = np.array([1, 0, 0]) #'starting' vector in x hat
#     #rotate xhat by gamma to get r2
#     r2 = np.array([cos(gamma), sin(gamma), 0])
#     #rotating xhat by beta to partial r3
#     pre_r3 = np.array([cos(beta), 0, sin(beta)])
#     #getting the two components of r3 after rotating pre_r3
#     r3_component1 = sin(theta) * r2
#     r3_component2 = cos(theta) * pre_r3
#     r3 = _normalize(r3_component1 + r3_component2)
#     return (r1, r2, r3)
#
# def _normalize(vector):
#     """returns a normalized vector. Input, output is an numpy array."""
#     norm = np.linalg.norm(vector)
#     if norm == 0:
#         return vector
#     return vector / norm
#
#
# angles = (0, 0, 120)
#
# # vectors = calculate_vectors(angles)
# # for v in vectors:
# #     print(v)
# #     print('')
#
radii = [53, 31, 167, 112, 87, 67, 56, 48, 42, 38, 190, 154, 118, 111, 98,
         88, 79, 71, 243, 194, 184, 176, 171, 166, 161, 156, 152, 149, 145,
         142, 136, 125, 114, 103, 94, 88, 265, 219, 212, 206, 198, 190, 183,
         178, 173, 169, 165, 161, 156, 145, 133, 123, 115, 108, 298, 253,
         0, 0, 247, 206, 205, 238, 231, 233, 225, 228, 226, 226, 222,
         217, 208, 200, 193, 188, 185, 180, 177, 174, 171, 156, 154, 143,
         135, 127, 120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

relative = list(radii)
sum = 0
count = 0
for i in range(len(relative)):
    relative[i] /= 298 #or use smallest 31

avg = 158.723
print(relative)
#
#
#
#
# '''
# if bondtype == 0 or bondtype == 2:
#     for i in range(0, list_length, 1): # not skipping every other
#         center = allcenters[i]
#         for j in range(i+1, list_length):
#             test_point = allcenters[j]
#             dist = np.linalg.norm(center - test_point)
#             if dist < a: # dist less than a unit
#                 pair = getpair(center, test_point)
#                 bondpairs.append(pair)
# if bondtype == 1 or bondtype == 2:
#     for i in range(list_length):
#         center = allcenters[i]
#         for j in range(i+1, list_length):
#             test_point = allcenters[j]
#             dist = np.linalg.norm(center - test_point)
#             if dist == a: # dist equal to a unit
#                 pair = getpair(center, test_point)
#                 bondpairs.append(pair)
# '''
