#Projectiles program

#Import
import math
import numpy as np
import matplotlib.pyplot as plt

#Defining parameters
theta = int(input("Enter theta"))
velocity = int(input("Enter velocity"))
height = int(input("Enter height"))

#Calculations

horizontalvel = velocity*np.cos(math.radians(theta))
verticalvel = velocity*np.sin(math.radians(theta))
if height == 0 and theta>=0:
    time = 2*((0-verticalvel)/(-9.81))
    vertical = time*horizontalvel
    maxheight = 0.5*(verticalvel)*time
    if theta == 0:
        vertical /= 2
    print("Vertical: ", vertical)
    print("Max height: ", maxheight)
    print("Time: ", time)
elif height>0 and theta >=0:
    halftime1 = (0-verticalvel)/(-9.81)
    maxheight = height+verticalvel*halftime1+0.5*(-9.81)*(halftime1**2)
    halftime2 = math.sqrt(4*0.5*9.81*maxheight)/9.81
    time = halftime1+halftime2
    vertical=horizontalvel*time
    print("Vertical: ", vertical)
    print("Max height: ", maxheight)
    print("Half time 1: ", halftime1, "Half time 2: ", halftime2, "Total time: ", time)
elif height>0 and theta<0:
    maxheight=height
    verticalvel = abs(verticalvel)
    finalvelocity = math.sqrt((verticalvel**2)+(2*9.81*height))
    time = (finalvelocity-verticalvel)/9.81
    vertical = time*horizontalvel
    print("Vertical: ", vertical)
    print("Max height: ", maxheight)
    print("Time: ", time)
    
#Plot projectile

#Define quadratic
def f(x,b,k):
    return k*(x)*(x-b)
xvertex = vertical/2
k = maxheight/((xvertex)*(xvertex-vertical))
ylist=f(xlist,(vertical/2), k)

#Plot
plt.plot(xlist, ylist)
plt.xlabel("Vertical (meters)")
plt.ylabel("Height (meters)")
plt.show()
