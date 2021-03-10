# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 09:27:15 2021

@author: Eric
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

#Parameters of interest
communitySize = 50 #number of populations in community
K = 100 #Max number of cells in droplet
N0 = 1 #Initial starting condition // may change later to be different per pop
maxTime = 168
minDetectCell = 25

#Functions - Scale data x between values a,b
def scaleData(x,a,b):
    return((b-a)*((x-np.min(x))/ (np.max(x) - np.min(x)) + a ))

"""
We'll start with the parameters shape = 0.6 and scale = 0.9 for defining our
growth rates in the community - see exploreDistro.py to explore the parameter
space of the gamma distribution

"""

growthRates = np.random.gamma(.6,.9,communitySize)

growthRates = scaleData(growthRates,0,1)

#time points
t = np.linspace(0,maxTime)

#define logistic growth model with carrying capacity max droplet capacity
def model(N,t,r,K):
    dydt = r*((K-N)/K)*N
    return dydt

#Solve for all members of community
NSolved = np.array([])
for i in growthRates:
    xs = np.asarray(odeint(model,N0,t,args=(i,K,))).transpose()
    NSolved = np.vstack([NSolved,xs]) if NSolved.size else xs

'''
Plotting
'''

#multipanel plot
fig, axs = plt.subplots(2)
axs[0].set_title("Growth rate distribution")
axs[0].hist(growthRates,bins=75,edgecolor='black')
axs[0].set(xlabel="cells h-1",ylabel="count")

#Visualize number of cells over time
for pop in NSolved:
    axs[1].plot(t,pop,alpha=0.25)

plt.xlabel("Time (hr)")
plt.ylabel("Cells")

#Calculate the % of species that are above minDetectCell cells at a given timepoint
p = []
for i in NSolved[:,].transpose():
    p.append(((i > minDetectCell).sum())/K)

#change to percentage to fit on existing plot axis
p= [i * 100 for i in p]

plt.plot(t,p,linewidth=2,c= 'black')
plt.hlines(minDetectCell,t.min(),t.max(),colors='black',linestyle='dashed')

#finish plot and output to file
plt.savefig('fig1.png')
plt.show()

'''
Writing a table to file
'''
#Write a table with number of cells on ech day after 7 days
hrs = [24,48,72,96,120,144,168]
inds = np.in1d(t,hrs).nonzero() #get indices
inds = np.asarray(inds).squeeze() #get rid of extra stuff

NSolved = np.around(NSolved,decimals=0).astype(int)

df = pd.DataFrame(data = NSolved[:,inds], columns=hrs, index=range(1,communitySize+1))

df.to_csv("table1.txt",sep='\t')
