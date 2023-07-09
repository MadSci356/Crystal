import numpy as np

'''def addvectors(a, b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]'''

#direction is that of the respective primitive vectors
r1len = 2 #length in first dir
r2len = 2 #length in second dir
r3len = 2 #length in third dir

a = 1 #length/space between atoms

#primitive vector below in [x_hat, y_hat, z_hat]
#primitive vectors are unit vectors
#also will be np arrays 
r1 = np.array(1.0, 0.0, 0.0) #[3.0, 0.0, 0.0]
r2 = np.array(0.0, 1.0, 0.0)
r3 = np.array(0.0, 0.0, 1.0)

r1 = r1*a
r2 = r2*a
r3 = r3*a

r1centers = []
r2centers = []
r3centers = []

for i in range(r1len+1):  #adding points in dir of 1st primitive vec (making a line)
	new_center = np.array()
	if i == 0:
		new_center = np.array(0,0,0) #initial center at origin
	else:
		#new_center = addvectors(r1centers[-1], r1)
		new_cener = r1centers[-1] + r1
	r1centers.append(new_center)
    
	for k in ange(r2len): #adding points r3 dir from line of points crossing origin
		if k == 0:
            #new_center = addvectors(r1centers[-1], r3)
	    	new_center = r1centers[-1] + r3
        else:
            #new_center = addvectors(r3centers[-1], r3)
	    	new_center = r3centers[-1] + r3
        r3centers.append(new_center)
        
    for j in range(r2len): #adding points in dir of 2nd primitive vec (making a 2d grid)
        if j == 0:
            #new_center = addvectors(r1centers[-1], r2)
	    	new_center = r1centers[-1] + r2
        else:
            #new_center = addvectors(r2centers[-1], r2)
			new_center = r2centers[-1] + r2
        r2centers.append(new_center)
        
		#adding points in r3 direction 
		#to the points in the grid (making a 3d lattice)
        for k in range(r2len): 
            if k == 0:
                #new_center = addvectors(r2centers[-1], r3)
				new_center = r2centers[-1] + r3
            else:
                #new_center = addvectors(r3centers[-1], r3)
				new_center = r3centers[-1] + r3
            r3centers.append(new_center)
    
print (r1centers)
print ('')
print (r2centers)
print ('')
print (r3centers)
'''
i = 0
for center in r2centers:
    print (center)
    i+=1
    if i > 1:
        print ('\n')
        i=0
print ('end r2')'''
'''i = 0
for center in r3centers:
    print center
    i+=1
    if i > 1:
        print '\n'
        i=0'''
        
print ("Number of points in lattice", len(r1centers) + len(r2centers) + len(r3centers))
