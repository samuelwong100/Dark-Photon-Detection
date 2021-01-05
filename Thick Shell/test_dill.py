# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 02:22:32 2021

@author: samue
"""
import dill

with open('sigma2_function', 'rb') as file:
    s2app_func = dill.load(file)

