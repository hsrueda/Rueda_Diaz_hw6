# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 19:16:41 2014

@author: Sony
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# Se selecciona el archivo
archivo = sys.argv[1]
documento = archivo.split(".")

# Se cargan los datos
datos = np.loadtxt(archivo)

# Grafica de cada archivo 
plt.scatter(datos[:,1],datos[:,2])
plt.title(documento[0] )
plt.savefig(documento[0]+".pdf")