# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 21:13:33 2021

@author: samuel
"""
from Thick_Shell_Class import Dark_Photon_Thick_Shell
from Thick_Shell_Class import s2app_func
import numpy as np
import matplotlib.pyplot as plt

class Dark_Photon_Thick_Shell_Experiment():
    def __init__(self,v_sense,Vapp,R0_ft,R1_ft,thickness_in_cm,
                 number_of_epsilon,number_of_m,epsilon_lowlim=-10,
                 epsilon_highlim=0,m_lowlim=-11,m_highlim=-4.5):
        ft_to_ev_inverse = 1544855.50049
        self.v_sense = v_sense
        self.Vapp = Vapp
        self.R0_ft = R0_ft
        self.R1_ft = R1_ft
        self.thickness_in_cm = thickness_in_cm
        self.R0 = R0_ft * ft_to_ev_inverse
        self.R1 = R1_ft * ft_to_ev_inverse
        self.R2 = (R1_ft + thickness_in_cm * 0.0328084) * ft_to_ev_inverse
        self.number_of_epsilon = number_of_epsilon
        self.number_of_m = number_of_m
        self.epsilon_list = np.logspace(epsilon_lowlim,epsilon_highlim,
                                        number_of_epsilon)
        #this is the base 10 exponent range
        self.m_list = np.logspace(m_lowlim,m_highlim,number_of_m)

    def Delta_V(self,epsilon,m):
        s2app = s2app_func(epsilon,m,self.R1,self.R2,self.Vapp)
        DPTS = Dark_Photon_Thick_Shell(epsilon,m,self.R1,self.R2,s2app,
                                       fast=True)
        return self.Vapp - DPTS.Vpp(self.R0)

    def get_Delta_V_evaluated(self):
        Delta_V_evaluated = np.zeros(shape=(self.number_of_epsilon,
                                            self.number_of_m))
        for (i,m) in enumerate(self.m_list):
            for (j,epsilon) in enumerate(self.epsilon_list):
                Delta_V_evaluated[i][j] = self.Delta_V(epsilon,m)
        return Delta_V_evaluated
        
    def plot_phase_space(self):
        Delta_V_evaluated = self.get_Delta_V_evaluated()
        #check whether the potential difference is larger than v_sensitivity.
        #If so, already eliminated by experiment
        #if still smaller, then it's still possible
        Delta_V_possible = Delta_V_evaluated < self.v_sense
        plt.pcolormesh(self.m_list,self.epsilon_list,Delta_V_possible.T,
                       cmap="gray") #always plot the transpose with colormesh
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel('$\epsilon$')
        plt.xlabel('m (eV)')
        title = "Plimpton & Lawton \n (vsense={} V, Vapp={} V, R0={} ft, R1={} ft, thickness={} cm)".format(str(self.v_sense),
    str(self.Vapp),str(self.R0_ft),str(self.R1_ft),str(self.thickness_in_cm))
        plt.title(title)
        title = "Plimpton & Lawton (vsense={} V, Vapp={} V, R0={} ft, R1={} ft, thickness={} cm)".format(str(self.v_sense),
    str(self.Vapp),str(self.R0_ft),str(self.R1_ft),str(self.thickness_in_cm))
        plt.title(title)
        plt.savefig(title+".png",dpi=300)
        
if __name__=="__main__":
    v_sense = 10**(-6)
    Vapp = 3000
    R0_ft = 2
    R1_ft = 2.5
    thickness_in_cm = 0.1
    number_of_epsilon = 300
    number_of_m = 300
    m_lowlim=-25
    m_highlim=-5
    DPTSE = Dark_Photon_Thick_Shell_Experiment(v_sense,Vapp,R0_ft,R1_ft,
                                       thickness_in_cm,number_of_epsilon,
                                       number_of_m,m_lowlim=m_lowlim,
                                       m_highlim=m_highlim)
    DPTSE.plot_phase_space()