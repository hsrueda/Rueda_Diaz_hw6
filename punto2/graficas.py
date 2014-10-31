from pylab import *
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re

arc=sys.argv[1]
datos=np.loadtxt(arc)
kalpha=re.findall("\d+",arc)

t=datos[:,0]
x=datos[:,1]
y=datos[:,2]
z=datos[:,3]
Rt=6378137 #Radio de la tierra

#Grafica en plano XY
plt.plot(x,y)
plt.title('Posiciones de Proton X vs Y')
plt.xlabel('Trayectoria en eje X')
plt.ylabel('Trayectoria en eje Y')
plt.grid()
plt.savefig("trayectoria_xy_Ec="+str(kalpha[0])+"_alpha="+str(kalpha[1]))+".png", format='png', bbox_inches='tight', transparent=True)

#Grafica en 3D de trayectoria
fig=figure()
ax=Axes3D(fig)
ax.set_aspect('equal')
ax.set_xlabel("$x$",fontsize=20)
ax.set_ylabel("$y$",fontsize=20)
ax.set_zlabel("$z$",fontsize=20)
ax.set_title("$Trayectoria$", fontsize=25)
ax.plot(x/Rt,y/Rt,z/Rt)
plt.savefig('trayectoria3D_Ec='+str(kalpha[0])+"_alpha="+str(kalpha[1]))+'.png',format ='png', transparent=True)
