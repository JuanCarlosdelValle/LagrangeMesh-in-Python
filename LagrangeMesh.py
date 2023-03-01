#    Lagrange Mesh in Python v1.0
#            26 Feb 2023
#     Juan Carlos del Valle Rosales
#          jcdevaller@gmail.com

import numpy as np # Linear Algebra Package

import mpmath as mp  # Multiple Precision Arithmetics
from mpmath import *


# Mesh Data 
nroots =2  # Size
WP = 10      # Working Precision
mp.dps=WP    # Assingning WP
mp.pretty = True # Show Numbers in Stardard Notation

#Scaling Pameter
h=1

# Coulomb Potential
def VC(x):
    return -1/x



#File's Name with Roots
FileRoots = "Laguerre_"+str(nroots)+"_WP_"+str(WP)+".dat"  


# Importing and Assigning Roots 
file = open(FileRoots, "r")
lines = file.read().splitlines() 
r = list(map(mpmathify, lines)) # Strings to MP Numbers



# Kinetic Matrix

T = matrix(nroots) # Empty Kinetic Matrix T

# Insering Elements to T
for j in range(0,nroots):
    for i in range(0,j+1):
        if i<j:
            T[i,j] = (-1)**(i-j)*(r[i]+r[j])/((r[i]*r[j])**(1/2)*(r[i]-r[j])**2) # Non-Diagonal Elements
            T[j,i] = T[i,j] # 
        else:
            T[i,j]= (4+(4*nroots+2)*r[i]-r[i]**2)/(12*r[i]**2) # Diagonal Elements
            
# Potential Matrix

V = matrix(nroots) # Empty Potential Matrix V

# Insering Elements to T
for i in range(0,nroots):
    V[i,i] = VC(h*r[i])  # Diagonal Elements


# Hamiltionan Matrix 

H = (1/(2*h**2))*T + V

# Diagonalization of H
Energies, S = mp.eigsy(H) # S is such that S^T*H*S = diag(energies), where T stands for transposition 


#Exporting Matrix S
CS=[]
for i in range(0,nroots):
    Crow=[]
    for j  in range(0,nroots):
        Crow.append(S[i,j])
    CS.append(Crow)    
    
np.savetxt("MatrixS.csv",CS,delimiter =", ",fmt ='% s')    
    
    
# Exporting Energies   
CEnergies=[]
for i in range(0,nroots):
    CEnergies.append(Energies[i])      
np.savetxt("Energies.csv",CEnergies,delimiter =", ",fmt ='% s')            
     