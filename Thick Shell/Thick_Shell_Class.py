# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 02:48:13 2021

@author: samuel
"""
import dill
import numpy as np
import matplotlib.pyplot as plt

with open('sigma2_function', 'rb') as file:
    s2app_func = dill.load(file)
    
with open('V_function', 'rb') as file:
    V_combine_func = dill.load(file)
    
with open('Vprime_function', 'rb') as file:
    Vprime_combine_func = dill.load(file)
    
with open('b1_function', 'rb') as file:
    b1_func = dill.load(file)
    
with open('b2_function', 'rb') as file:
    b2_func = dill.load(file)
    
with open('b3_function', 'rb') as file:
    b3_func = dill.load(file)

with open('b4_function', 'rb') as file:
    b4_func = dill.load(file)
    
with open('s1_function', 'rb') as file:
    s1_func = dill.load(file)

class Dark_Photon_Thick_Shell():
    def __init__(self,epsilon,m,R1,R2,s2,fast=False):
        #free parameters
        self.epsilon = epsilon
        self.m = m
        self.R1 = R1
        self.R2 = R2
        self.s2 = s2
        #derived parameters
        self.m_reduced = m/np.sqrt(1+epsilon**2)
        self.mb = self.m_reduced #alias
        if not fast: #don't compute these constants for faster code
            self.s1 = s1_func(epsilon,m,R1,R2,s2)
            self.b1 = b1_func(epsilon,m,R1,R2,s2)
            self.b2 = b2_func(epsilon,m,R1,R2,s2)
            self.b3 = b3_func(epsilon,m,R1,R2,s2)
            self.b4 = b4_func(epsilon,m,R1,R2,s2)
            self.Q = self.get_charge()
        #vectorized functions
        self.V = np.vectorize(self.V_unvectorize)
        self.Vprime = np.vectorize(self.Vprime_unvectorize)
        self.rho_bulk = np.vectorize(self.rho_bulk_unvectorize)
        self.delta = np.vectorize(self.delta_unvectorize,excluded=['r0','dr'])
        
    def get_charge(self):
        b3_coeff = np.exp(self.mb*self.R2) * (self.mb*self.R2 - 1) \
            - np.exp(self.mb*self.R1) * (self.mb*self.R1 - 1)
        b4_coeff = np.exp(-self.mb*self.R2) * (self.mb*self.R2 + 1) \
            - np.exp(-self.mb*self.R1) * (self.mb*self.R1 + 1)
        bracket = (self.b3 / (self.mb**3)) * b3_coeff - (self.b4 / (self.mb**2)) * b4_coeff
        Q_nocoeff = self.epsilon * (self.mb**2) * bracket \
            + self.s1 * (self.R1**2) + self.s2 * (self.R2**2)
        return 4*np.pi*Q_nocoeff
    
    def V_unvectorize(self,r):
        return V_combine_func(r,self.epsilon,self.m,self.R1,self.R2,self.s2)
        
    def Vprime_unvectorize(self,r):
        return Vprime_combine_func(r,self.epsilon,self.m,self.R1,self.R2,self.s2)
    
    def rho_bulk_unvectorize(self,r):
        if self.R1 < r < self.R2:
            return self.epsilon * (self.mb**2) * self.Vprime_unvectorize(r)
        else:
            return 0
        
    def delta_unvectorize(self,r,r0,dr):
        if  r0 - dr/2 <= r <= r0 + dr/2:
            return 1/dr
        else:
            return 0
    
    def rho_EM(self,r,dr):
        return self.rho_bulk(r) + self.s1*self.delta(r,self.R1,dr) \
            + self.s2*self.delta(r,self.R2,dr)
    
    def Vpp(self,r):
        return self.V(r) + self.epsilon*self.Vprime(r)
    

# epsilon=10
# m=20
# R1=1
# R2=2
# s2=10
# DPTS = Dark_Photon_Thick_Shell(epsilon,m,R1,R2,s2)
# r_linspace=np.linspace(0,3,10000)
# dr = r_linspace[1] - r_linspace[0]

#plt.plot(r_linspace,DPTS.V(r_linspace))


