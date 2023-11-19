import pandas as pd
import streamlit as st
from scipy import optimize
import math
import matplotlib.pyplot as plt
import numpy as np

#hydrolysis constants for zirconium
k1=0.3
k2=-1.7
k3=-5.1
k4=-9.7
k5=-16
k6=-0.6
k7=-3.7
k8=6
k9=-1.9
k10=-2.57

#kw value
kw=-14

def calc_zr_conc(x, ck, ccl):
# x is the molar concentration of variable x12 - [OH-]
  #print(ck)
  zr = np.zeros(14)
  zr[12] = math.log10(x) #lgoh
  zr[13] = math.log10(ck) #lgk
  zr[10] = kw - zr[12] #lgh
  zr[1] = k9 + 4*zr[10]
  zr[2] = k1 + k9 + 3*zr[10]
  zr[3] = k2 + k9 + 2*zr[10]
  zr[4] = k3 + k9 + zr[10]
  zr[5] = k4 + k9
  zr[6] = k5 + k9 - zr[10]
  zr[7] = k6 + 3*k9 + 8*zr[10]
  zr[8] = k7 + 3*k9 + 7*zr[10]
  zr[9] = k8 + 4*k9 + 8*zr[10]
  zr[11] = math.log10(2*ccl)
  return(zr)

def calc_zr_solubility(zr):
  szr = 10**zr[1] + 10**zr[2] + 10**zr[3] + 10**zr[4] + 10**zr[5] + 10**zr[6] + 10**zr[7] + 10**zr[8] + 10**zr[9]
  return(szr)

