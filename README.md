# crystal
Created while working as an undergrad researcher for Kumah Lab at NC State University (2017-2018)
Python Algorithms for creating 3D crystal models and lattice structures for use with Captive.

## This repository is archived. Newer version under development.

How to install on a Unix system:
I will be making it into a library for now this is how to install it.

1) Clone the crystal repo
2) Install docker
2) Start PyMesh docker with `docker run -v script/folder:/root -it qnzhou/pymesh` might have to sudo this
3) Then you can run your script with `exec(open("script.py").read())` in the interpreter. I suggest making `yourscript.py` in the   crystal repo folder (where main.py is) and then running that in that docker image.
4) Additional note, where ever you put the output obj file, make sure you put the `lattice-color.mtls` file there too. 
##### New note: You won't have to add the .mtl manually, I have automated this feature.


## Demo: 

https://github.com/MadSci356/Crystal/assets/5410205/2962ca48-854f-4a17-a575-d3bcf5c91bc4



![FCC Lattice](https://github.com/MadSci356/Crystal/assets/5410205/a64eacf1-d675-43f7-9b18-1ab962a3f5ea)

![Cubic Bonds](https://github.com/MadSci356/Crystal/assets/5410205/693cab2d-ec65-4cca-b0f2-350ac2c8cb90)

![Graphene](https://github.com/MadSci356/Crystal/assets/5410205/71a99104-d66d-408f-865e-1919eeef9ae6)
