import math
from math import sqrt

def distance(a1, a2):
    file = open('1ubq.pdb','r')
    lines = file.readlines()
    file.close()
    
    x = []
    y = []
    z = []
    
    for i in lines:
        x.append(i.split()[6]) 
        y.append(i.split()[7])
        z.append(i.split()[8])
        
    x1 = float(x[a1-1])
    y1 = float(x[a1-1])
    z1 = float(x[a1-1])
    x2 = float(x[a2-1])
    y2 = float(x[a2-1])
    z2 = float(x[a2-1])
    
    distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    return distance

print(distance(1, 2))
