def addvectors(a, b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]

#direction is that of the respective primitive vectors
r1len = 5 #length in first dir
r2len = 5 #length in second dir
r3len = 5 #length in third dir

a = 3 #length/space between atoms

#primitive vector below in [x_hat, y_hat, z_hat]
#primitive vectors are unit vectors
#also will be np arrays 
r1 = [3, 0, 0]
r2 = [0, 3, 0]
r3 = [0, 0, 3]

#r1 = r1*a
#r2 = r2*a
#r3 = r3*a

r1centers = [[0,0,0]]
r2centers = []
r3centers = []

for i in range(r1len):
    new_center = addvectors(r1centers[-1], r1)
    r1centers.append(new_center)

print r1centers

for r1center in r1centers:
    for j in range(r2len):
        new_center = []
        if j == 0:
            new_center = addvectors(r1center, r2)
        else:
            new_center = addvectors(r2centers[-1], r2)
        r2centers.append(new_center)

i = 0
for center in r2centers:
    print center
    i+=1
    if i > 4:
        print '\n'
        i=0
        
        

