import pandas as pd
import streamlit as st
from scipy import optimize
import math
import matplotlib.pyplot as plt
import numpy as np


#hydrolysis constants for yttrium
m1=-7.7
m2=-36.5
m3=-14.23
m4=-31.6
m5=-17.5
m6=-8.5

#kw value
kw=-14

def calc_yt_conc(y, ck, cn):
    # y is the variable y10 for log[OH-]
    # print(ck)
    yt = np.zeros(11)
    yt[8] = math.log10(ck)  # lgk
    yt[9] = math.log10(3 * cn)  # lgn
    yt[10] = math.log10(y)  # lgoh, y is not the same as
    yt[7] = kw - yt[10]  # lgh
    yt[1] = 3 * yt[7] - m5
    yt[2] = m1 + yt[1] - yt[7]
    yt[3] = m2 + yt[1] - 4 * yt[7]
    yt[4] = m3 + 2 * yt[1] - 2 * yt[7]
    yt[5] = m4 + 3 * yt[1] - 5 * yt[7]
    yt[6] = m6
    return (yt)

def calc_yt_netcharge(y):
  # y is the variable y10 for log[OH-]
  yt = np.zeros(11)
  yt = calc_yt_conc(y)
  pos = 3*10**yt[1] + 2*10**yt[2] + 4*10**yt[4] + 4*10**yt[5] + 10**yt[7] + 10**yt[8]
  neg = 10**yt[3] + 10**yt[9] + 10**yt[10]
  net = pos - neg
  return(net)

def calc_yt_solubility(yt):
  syt = 10**yt[1] + 10**yt[2] + 10**yt[3] + 10**yt[4] + 10**yt[5] + 10**yt[6]
  return(syt)

def calc_dict(cn):
    yt = np.zeros(10)
    pH_values = []
    ck_values = []
    cn_values = []
    syt_values = []
    dss_yt_values = []

    for x in range(1, 5000):
        ck = x / 1000  # between 0.001M and 5M increasing in 0.0001M each iteration, KOH added to solution
        a = 1e-20  # lower bound concentration
        b = ck  # upper bound concentration
        try:
            # for each given KOH concentration, what is the OH concentration at which net charge is 0 ("correct" amount of dissociation of KOH)
            # lets us find the [OH] and pH of the solution for each [KOH] iterated through
            # root found by algo is [OH]

            root, results = optimize.toms748(calc_yt_netcharge, a, b, full_output=True)

            yt = calc_yt_conc(root)  # finds concs of all yt species for the concentration of KOH added
            pH = -(kw - yt[10])  # finds pH for the concentration of KOH added
            syt = calc_yt_solubility(
                yt)  # finds syt for the conc of KOH added -> when plotted in contour, will be diff value for each [KOH] (directly affects by pH) and cn (affects by net charge)
            dss_yt = cn / syt  # finds dss_yt

            cn_values.append(cn)  # constant at 0.5
            ck_values.append(ck)  # depends on iteration in for loop
            pH_values.append(pH)
            syt_values.append(syt)
            dss_yt_values.append(dss_yt)
            # print('ck =', ck ,' ,cn=' , cn , ' pH =', pH , 'syt=',syt,' net=', calc_netcharge(root))

        except Exception as e:
            # print('ck =',ck,' ,cn=',cn,', y10 =', yt[10], ' pH =', pH, 'net=', calc_netcharge(root), 'Skipping...')
            #print(results)  # print to see number of iterations, confirm convergence
            st.write(e)
            continue

    data = {
        'pH': pH_values,
        'ck': ck_values,
        'cn': cn_values,
        'syt': syt_values,
        'dss_yt': dss_yt_values
    }
    return data



