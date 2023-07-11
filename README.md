# crystal
Python Algorithms for creating 3D crystal models and lattice structures for use with Captive.

How to install on a Unix system:
I will be making it into a library for now this is how to install it.

1) Clone the crystal repo
2) Install docker
2) Start PyMesh docker with `docker run -v script/folder:/root -it qnzhou/pymesh` might have to sudo this
3) Then you can run your script with `exec(open("script.py").read())` in the interpreter. I suggest making `yourscript.py` in the   crystal repo folder (where main.py is) and then running that in that docker image.
4) Additional note, where ever you put the output obj file, make sure you put the `lattice-color.mtls` file there too. 
##### New note: You won't have to add the .mtl manually, I have automated this feature.
