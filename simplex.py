# -*- coding: utf-8 -*-
"""
Created on Fri May 13 14:17:26 2016

@author: maxflugelman, lucasheinrich, fernandokabas
"""

import numpy as np


############################################################################

def pfourelement(t, r, l, c, k, v):

    # devuelve presion arterial (Pout) de windkessel de 4 elementos, k es rl

    u = (l*k+r*l)/(k*l*c*r*2)
    if t==0.0:
        # la funcion original incluye un delta que lo solucione con este if por que con metodo daba error
    
        return k*v*(1+1/k)
    else:
        return k*v*(np.exp(-u*t))*(np.cos((np.sqrt(u))*t)+np.sin((np.sqrt(u))*t))


############################################################################
def sort(m):
    #ordena el simplex de acuerdo a los valores de la funcion
    for i in range(0,len(m[0])-1):
        for j in range(len(m[0])-1,0,-1):
            if(func(m[:,j])<func(m[:,j-1])):
                temp = np.copy(m[:, j])
                m[:, j] = m[:, j-1]
                m[:, j-1] = temp
    return m;

############################################################################

def func (x):
    #funcion a minimizar
    result = -np.cos(x[0])*np.cos(x[1])*np.exp(-((x[0]-np.pi)**2+(x[1]-np.pi)**2))
    return result;
############################################################################
def getM(s):
    M = np.zeros(len(s))
    for i in range(0,len(s)):
        for j in range(0,len(s)):
            M[i]= M[i]+np.copy(s[i,j])
    
    M = M/len(s)
    return M;

############################################################################
    #reflection
def getR(m,p):
    R = 2*np.copy(m)-np.copy(p)
    
    return R;
    
############################################################################
    #expansion, m es el centro de simplex
def getE(m, r):
    E = 2*np.copy(r)-np.copy(m)
    return E;
############################################################################
    #outside contraction
def getCO(m, r):
    CO = 0.5*(np.copy(r)+np.copy(m))  
    return CO;
    
############################################################################
    #inside contraction
def getCI(m, r):
    CI = 0.5*(np.copy(m)-np.copy(r))
    return CI;
    
############################################################################
    #shrink simplex
def getS(s):
    for i in range(1,len(s)+1):
        s[:,i] = 0.5*(s[:,0]+s[:,i])
    return s;
    

############################################################################
def Nelder (x0, iter):
    size = len(x0)+1
    t = len(x0)
    simplex = np.zeros((len(x0),size),float) #genera una matriz de ceros
    for i in range (0,len(x0)):
        for j in range (0,size):
            if (j==0):      #asigna el vector inicial a la primer columna
                simplex[i,j]=x0[i]
            elif(j-i==1):   #asigna el corrimiento a cada direccion
                simplex[i,j]=x0[i]+0.05
            else:   #completa la matriz con los valores que no se modificaban
                simplex[i,j]=x0[i]
    #ordeno el simplex de acuerdo al minimo valor de la funcion objetivo
    
    simplex = sort(simplex)


    for i in range (0,iter):
         
        P = simplex[:,t]   # peor punto
        O = simplex[:,0]      # mejor punto
        M = getM(simplex)
        #punto de reflexion
            
        R = getR(M,P)
        
        if (func(R)<func(simplex[:,1])):

            if (func(getE(M,R))<func(O)):
                simplex[:,t] = np.copy(getE(M,R))
            else:
                simplex[:,t] = np.copy(R)
        else:
            if(func(R)<func(P)):  
                
                if(func(getCO(M,R))<func(R)):
                    simplex[:,t] = np.copy(getCO(M,R))
                    
                    
                elif(func(getCI(M,R))<func(R)):
                    simplex[:,t] = np.copy(getCI(M,R))
            else:
                simplex = np.copy(getS(simplex))
            simplex = sort(simplex)
    
                
                    
        
        
        
        

    
    print (simplex[:,0])
        

    return simplex[:,0]
    
 ############################################################################   
    




 

#####################M     MAIN   RUN    ####################################

x0 = np.array ([3., 3.,0.])

s = Nelder(x0,500)

print (func(s))
