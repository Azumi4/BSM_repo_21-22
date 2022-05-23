# -*- coding: utf-8 -*-
"""
Created on Mon May 23 09:54:54 2022

@author: azumi
"""

import numpy as np
import matplotlib.pyplot as plt
import json

def myplot(x,y,xlabel,ylabel,title):
    plt.figure()
    plt.plot(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
fname = open ("HHParamaters.json")
params  = json.load(fname)

gna_max = float(params['gna_max'])
gk_max = float(params['gk_max'])
geq_max = float(params['geq_max'])
C = float(params['C'])
Ena = float(params['Ena'])
Ek = float(params['Ek'])
Eeq = float(params['Eeq'])
V0 = float(params['V0'])
T0 = float(params['T0'])
Tend = float(params['Tend'])
dt = float(params['dt'])
I = float(params['i'])

t=np.arange (T0, Tend+dt, dt)
L= len(t)
n = np.zeros(L)
m = np.zeros(L)
h = np.zeros(L)

ah = np.zeros(L)
an = np.zeros(L)
am = np.zeros(L)
bh = np.zeros(L)
bn = np.zeros(L)
bm = np.zeros(L)
V = np.zeros(L)

V[0]= V0
am[0]= 0.1*(-40-V[0])/(np.exp((-40-V0)/10)-1)
bm[0]= 4*np.exp((-65-V0)/18)
ah[0]= 0.07*np.exp((-65-V0)/20)
bh[0]=1/(np.exp((-35-V0)/10) + 1)
an[0]=0.01*(-55-V0)/(np.exp((-55-V0)/10)-1)
bn[0]=0.125*np.exp((-65-V0)/80)

m[0]=am[0]/(am[0]+bm[0]) 
n[0]=an[0]/(an[0]+bn[0])
h[0]=ah[0]/(ah[0]+bh[0])

for i in range(L-1):
    dV = (I- gna_max*m[i]**3 *h[i]*(V[i]-Ena) - gk_max*n[i]**4 *(V[i]-Ek) - geq_max*(V[i]-Eeq))/C
    dm = am[i]*(1-m[i]) - bm[i]*m[i]
    dn = an[i]*(1-n[i]) - bn[i]*n[i]
    dh = ah[i]*(1-h[i]) - bh[i]*h[i]
    
    
    m[i+1] = m[i]+dt*dm
    V[i+1] = V[i]+dt*dV
    n[i+1] = n[i]+dt*dn
    h[i+1] = h[i]+dt*dh
    
    am[i+1]= 0.1*(-40-V[i+1])/(np.exp((-40-V[i+1])/10)-1)
    bm[i+1]= 4*np.exp((-65-V[i+1])/18)
    ah[i+1]= 0.07*np.exp((-65-V[i+1])/20)
    bh[i+1]=1/(np.exp((-35-V[i+1])/10) + 1)
    an[i+1]=0.01*(-55-V[i+1])/(np.exp((-55-V[i+1])/10)-1)
    bn[i+1]=0.125*np.exp((-65-V[i+1])/80)

myplot(t,V,"t","V","VT")
myplot(t,h,"t","h","hT")
myplot(t,m,"t","m","mT")
myplot(t,n,"t","n","nT")
    
        
    


