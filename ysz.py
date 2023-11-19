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

#defining preliminary function and calculations
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

def calc_yt_solubility(yt):
  syt = 10**yt[1] + 10**yt[2] + 10**yt[3] + 10**yt[4] + 10**yt[5] + 10**yt[6]
  return(syt)

def calc_zr_solubility(zr):
  szr = 10**zr[1] + 10**zr[2] + 10**zr[3] + 10**zr[4] + 10**zr[5] + 10**zr[6] + 10**zr[7] + 10**zr[8] + 10**zr[9]
  return(szr)

def calc_yt_zr_netcharge(x, ck, ccl, cn):
  # x is the molar concentration of variable x12 - [OH-]
  zr = np.zeros(14)
  zr = calc_zr_conc(x, ck, ccl)
  pos_zr = 4*10**zr[1] + 3*10**zr[2] + 2*10**zr[3] + 10**zr[4] + 8*10**zr[7] + 7*10**zr[8] + 8*10**zr[9] + 10**zr[13] # Remove + 10**zr[10] - avoid double count [H+] in Yt
  neg_zr = 10**zr[6] + 10**zr[11] # Remove + 10**zr[12] - avoid double count [OH-] in Yt

  # x is the variable y10 for log[OH-]
  yt = np.zeros(11)
  yt = calc_yt_conc(x, ck, cn)
  pos_yt = 3*10**yt[1] + 2*10**yt[2] + 4*10**yt[4] + 4*10**yt[5] + 10**yt[7] + 10**yt[8]
  neg_yt = 10**yt[3] + 10**yt[9] + 10**yt[10]

  net = pos_zr - neg_zr + pos_yt - neg_yt
  return(net)

def calc_dict(cn, ccl):
    # define lists to be used in plotting functions
    yt = np.zeros(10)
    zr = np.zeros(14)
    pH_values = []  # list of pH values
    ck_values = []  # list of KOH concentrations
    cn_values = []  # list of Y(NO3)3 concentrations
    ccl_values = []  # list of ZrCl3 concentrations
    szr_values = []  # solubility of zirconium
    syt_values = []  # solubility of yttrium
    dss_zr_values = []  # degree of supersaturation - Zirconium
    dss_yt_values = []  # degree of supersaturation - Yttrium

    for x in range(1, 250):
        ck = x / 500  # 0.002M to 1M
        a = 1e-20  # lower bound concentration
        b = ck  # upper bound concentration
        try:
            root, results = optimize.toms748(calc_yt_zr_netcharge, a, b, args=(ck, ccl, cn), full_output=True)
            # finding correct [OH] for given KOH in this iteration
            # solubility lists for zr and yt
            zr = calc_zr_conc(root, ck, ccl)
            yt = calc_yt_conc(root, ck, cn)
            pH = -(kw - zr[12])
            szr = calc_zr_solubility(zr)
            dss_zr = ccl / szr
            syt = calc_yt_solubility(yt)
            dss_yt = cn / syt

            ccl_values.append(ccl)  # constant
            ck_values.append(ck)  # value found by iteration
            cn_values.append(cn)  # constant
            pH_values.append(pH)  # correct dissociated OH for this KOH conc

            szr_values.append(szr)
            dss_zr_values.append(dss_zr)

            syt_values.append(syt)
            dss_yt_values.append(dss_yt)
            #print('cn=', cn, 'ck =', ck, ' ,ccl=', ccl, ' [OH-]=', zr[12], ' pH =', pH, 'szr=', szr, ' net=', calc_yt_zr_netcharge(root))
        except Exception as e:
            #print('cn=', cn, 'ck =', ck, ' ,ccl=', ccl, ', [OH-] =', zr[12], ' pH =', pH, 'net=', calc_yt_zr_netcharge(root), 'Skipping...')
            # print(results) #print to see number of iterations, confirm convergence
            print(e)
            continue

    data = {
        'pH': pH_values,
        'ck': ck_values,
        'cn': cn_values,
        'ccl': ccl_values,
        'szr':szr_values,
        'syt':syt_values,
        'dss_zr':dss_zr_values,
        'dss_yt':dss_yt_values
    }
    return data

def contour_yt(ccl):
    #net charge graph only has root when ZrOCl2 = 3.50M, otherwise no crossing-over point.

    #defining lists
    pH_values = []  # list of pH values
    ck_values = []  # list of KOH concentrations
    cn_values = []  # list of Y(NO3)3 concentrations
    syt_values = []  # solubility of yttrium
    dss_yt_values = []  # degree of supersaturation - Yttrium

    for x in range(1,11):
        for y in range(1,11):
            ck = x/100 #improve acuracy of plot by dividing this more
            cn = y/100
            a = 1e-20
            b = ck
            try:
                root, results = optimize.toms748(calc_yt_zr_netcharge, a, b, args=(ck, ccl, cn), full_output=True)
                yt = calc_yt_conc(root, ck, cn)
                zr = calc_zr_conc(root, ck, cn)

                pH = -(kw - zr[12])
                syt = calc_yt_solubility(yt)
                dss_yt = cn / syt

                ck_values.append(ck)
                cn_values.append(cn)
                pH_values.append(pH)

                syt_values.append(syt)
                dss_yt_values.append(dss_yt)

                #print('cn=', cn, ' ck =', ck, ' ,ccl=', ccl, ' pH =', pH, 'szr=', szr, ' net=', calc_yt_zr_netcharge(root))
            except Exception as e:
                #print('ck =', ck, ' ,ccl=', ccl, ', [OH-] =', zr[12], ' pH =', pH, 'net=', calc_yt_zr_netcharge(root), 'Skipping...')
                # print(results) #print to see number of iterations, confirm convergence
                syt_values.append(None)
                pH_values.append(None)
                st.write(e)
                break



    X1 = np.reshape(ck_values, (10,10))
    X2 = np.reshape(pH_values, (10,10))
    Y = np.reshape(cn_values, (10,10))
    Z = np.reshape(syt_values, (10,10))

    data = {
        'KOH': X1,
        'pH': X2,
        'Yt': Y,
        'syt': Z
    }

    return data

def contour_zr(cn):
    #defining lists
    pH_values = []  # list of pH values
    ck_values = []  # list of KOH concentrations
    ccl_values = []  # list of Y(NO3)3 concentrations
    szr_values = []  # solubility of yttrium
    dss_zr_values = []  # degree of supersaturation - Yttrium

    for x in range(1,11):
        for y in range(1,11):
            ck = x/100 #improve acuracy of plot by dividing this more
            ccl = y/10
            a = 1e-20
            b = ck
            try:
                root, results = optimize.toms748(calc_yt_zr_netcharge, a, b, args=(ck, ccl, cn), full_output=True)
                zr = calc_zr_conc(root, ck, ccl)

                pH = -(kw - zr[12])
                szr = calc_yt_solubility(zr)
                dss_zr = ccl / szr

                ck_values.append(ck)
                ccl_values.append(ccl)
                pH_values.append(pH)

                szr_values.append(szr)
                dss_zr_values.append(dss_zr)

                #print('cn=', cn, ' ck =', ck, ' ,ccl=', ccl, ' pH =', pH, 'szr=', szr, ' net=', calc_yt_zr_netcharge(root))
            except Exception as e:
                #print('ck =', ck, ' ,ccl=', ccl, ', [OH-] =', zr[12], ' pH =', pH, 'net=', calc_yt_zr_netcharge(root), 'Skipping...')
                # print(results) #print to see number of iterations, confirm convergence
                szr_values.append(None)
                pH_values.append(None)
                st.write(e)
                break

    X1 = np.reshape(ck_values, (10,10))
    X2 = np.reshape(pH_values, (10,10))
    Y = np.reshape(ccl_values, (10,10))
    Z = np.reshape(szr_values, (10,10))

    data = {
        'KOH': X1,
        'pH': X2,
        'Zr': Y,
        'szr': Z
    }
    return data

def contour_ysz(ck):
    pH_values = []  # list of pH values
    ck_values = []  # list of KOH concentrations
    cn_values = []  # list of Y(NO3)3 concentrations
    ccl_values = []  # list of ZrCl3 concentrations
    szr_values = []  # solubility of zirconium
    syt_values = []  # solubility of yttrium
    dss_zr_values = []  # degree of supersaturation - Zirconium
    dss_yt_values = []  # degree of supersaturation - Yttrium

    for x in range(1, 11):
        for y in range(1, 11):
            ccl = x / 10
            cn = y / 10
            a = 1e-20  # lower bound concentration
            b = ck  # upper bound concentration
            try:
                root, results = optimize.toms748(calc_yt_zr_netcharge, a, b, args=(ck, ccl, cn), full_output=True)
                zr = calc_zr_conc(root, ck, ccl)
                yt = calc_yt_conc(root, ck, cn)
                pH = -(kw - zr[12])
                szr = calc_zr_solubility(zr)
                dss_zr = ccl / szr
                syt = calc_yt_solubility(yt)
                dss_yt = cn / syt

                if ccl not in ccl_values:
                    ccl_values.append(ccl)

                if cn not in cn_values:
                    cn_values.append(cn)

                ck_values.append(ck)
                pH_values.append(pH)

                szr_values.append(szr)
                dss_zr_values.append(dss_zr)

                syt_values.append(syt)
                dss_yt_values.append(dss_yt)

                #print('cn=', cn, ' ck =', ck, ' ,ccl=', ccl, ' pH =', pH, 'szr=', szr, ' net=', calc_yt_zr_netcharge(root))

            except Exception as e:
                pH_values.append(0.0)
                szr_values.append(0.0)
                syt_values.append(0.0)
                #st.write(e)
                continue

    try:
        X, Y = np.meshgrid(ccl_values, cn_values)
        print(len(szr_values), len(syt_values), len(pH_values))
        Z1 = np.reshape(szr_values, (len(X), len(Y)))
        Z2 = np.reshape(syt_values, (len(X), len(Y)))
        Z3 = np.reshape(pH_values, (len(X), len(Y)))

        """
        X = np.reshape(ccl_values, (10,10))
        Y = np.reshape(cn_values, (10,10))
        Z1 = np.reshape(szr_values, (len(X), len(Y)))
        Z2 = np.reshape(syt_values, (len(X), len(Y)))
        Z3 = np.reshape(pH_values, (len(X), len(Y)))
        """

        data = {
            'X': X,
            'Y': Y,
            'szr': Z1,
            'syt': Z2,
            'pH': Z3
        }
        return data

    except Exception as e:
        st.write("reshaping error")

