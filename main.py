import pandas as pd
import streamlit as st
import ysz
from scipy import optimize
import math
import matplotlib.pyplot as plt
import numpy as np

#st.set_page_config(layout='wide')

st.title("Solubility Maps for YSZ")
st.subheader("This model will generate solubility maps for the hydrothermal synthesis of YSZ based on the input concentrations entered below.")
st.subheader('\n This model is currently assumed to be operated at 25°C.')

col1, col2 = st.columns(2)
with col1:
    cn = st.number_input('Input [Y(NO3)3]/M to be added:', value=1.00, min_value=0.0, max_value=5.0, step=0.0001)
    ccl = st.number_input('Input [ZrOCl2]/M to be added:', value = 1.00, min_value=0.0, max_value=5.0, step=0.0001)
    ck = st.number_input('Input [KOH]/M to be added:', min_value=0.0, max_value=1., value=0.5, step=0.0001)

# plotting the solubility&dss against pH&KOH curves
st.header("Graphs of Solubilities of Yt and Zr")
st.write('Plotted against pH and [KOH]/M')
dict = ysz.calc_dict(cn, ccl)
solubility_data = pd.DataFrame(dict)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
plt.rcParams["figure.figsize"] = (20, 20)
# titlestr = 'Vary [KOH], cn=' + str(cn) + ', ccl=' + str(ccl)

ax1.semilogy(dict['pH'], dict['syt'], label='Yt')
ax1.semilogy(dict['pH'], dict['szr'], label='Zr')
ax1.set(xlabel='pH', ylabel='Solubility')
ax1.grid()
ax1.legend(['Y', 'Zr'])

ax2.semilogy(dict['pH'], dict['dss_zr'], dict['pH'], dict['dss_yt'])
ax2.set(xlabel='pH', ylabel='Deg of Supersaturation')
ax2.grid()
ax2.legend(['Y', 'Zr'])

ax3.semilogy(dict['ck'], dict['syt'], dict['ck'], dict['szr'])
ax3.set(xlabel='[KOH]/M', ylabel='Solubility')
ax3.grid()
ax3.legend(['Y', 'Zr'])

ax4.semilogy(dict['ck'], dict['dss_yt'], dict['ck'], dict['dss_zr'])
ax4.set(xlabel='[KOH]/M', ylabel='Degree of Supersaturation')
ax4.grid()
ax4.legend(['Y', 'Zr'])
st.pyplot(fig)

def plot_contour_yt(ccl):
    # plotting contour plot for Y(NO3)3 against KOH, with syt lines
    st.header("Contour plots for Y(NO3)3 against KOH and pH")
    st.write('Blue: Syt, Red: pH')
    data = ysz.contour_yt(ccl)
    fig2, ax1 = plt.subplots(1, 1)

    CS1 = ax1.contour(data['KOH'], data['Yt'], data['syt'], colors='blue')
    ax1.clabel(CS1, CS1.levels, inline=True, fontsize=20)

    CS2 = ax1.contour(data['KOH'], data['Yt'], data['pH'], colors='red')
    ax1.clabel(CS2, CS2.levels, inline=True, fontsize=20)

    ax1.set_ylabel("[Y(NO3)3]/M")
    ax1.set_xlabel("[KOH]/M")
    st.pyplot(fig2, figsize=(10, 5))

def plot_contour_zr(cn):
    # plotting contour plot for ZrOCl2 afainst KOH, with szr lines

    st.header("Contour plots for ZrOCl2 against KOH and pH.")
    st.write('Green: S_Zr, Red: pH')
    data = ysz.contour_zr(cn)
    fig3, ax2 = plt.subplots(1, 1)

    ax2.set_ylabel("[ZrOCl2]/M", fontsize=30)
    ax2.set_xlabel("[KOH]/M", fontsize=30)
    plt.rcParams.update({'font.size': 20})

    CS1 = ax2.contour(data['KOH'], data['Zr'], data['szr'], colors='green')
    ax2.clabel(CS1, CS1.levels, inline=True, fontsize=20)

    CS2 = ax2.contour(data['KOH'], data['Zr'], data['pH'], colors='red')
    ax2.clabel(CS2, CS2.levels, inline=True, fontsize=20)

    ax2.legend(['S_Zr', 'pH'])
    st.pyplot(fig3, figsize=(40, 40))

def plot_contour_ysz(ck):
# contour plot of yt/zr with syt and szr and pH
    st.header("Final Solubility Contour Surface Plot")
    #st.write("Blue: S_yt")
    #st.write("Green: S_zr")
    #st.write("Red: pH")

    data = ysz.contour_ysz(ck)

    length = len(data['X'])

    #checking to see all lengths are the same size so that reshaping will work
    # for i in range(length):
        # print('length of array number ' + str(i) + ' for X:' + str(len(data['X'][i])))
        # print('length of array number ' + str(i) + ' for Y:' + str(len(data['Y'][i])))
        # print('length of array number ' + str(i) + ' for szr:' + str(len(data['szr'][i])))
        # print('length of array number ' + str(i) + ' for syt:' + str(len(data['syt'][i])))
        # print('length of array number ' + str(i) + ' for pH:' + str(len(data['pH'][i])))

    #print(data['syt'])
    #print(data['szr'])


    #ignore the plot stuff for now 
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel("[ZrOCl2]/M", fontsize=30)
    ax.set_ylabel("[Y(NO3)3]/M", fontsize=30)
    plt.rcParams.update({'font.size': 15})

    cs1 = ax.contour(data['X'], data['Y'], data['szr'], colors='green')
    ax.clabel(cs1, cs1.levels, inline=True, fontsize=20)

    cs2 = ax.contour(data['X'], data['Y'], data['syt'], colors='blue')
    ax.clabel(cs2, cs2.levels, inline=True, fontsize=20)

    cs3 = ax.contour(data['X'], data['Y'], data['pH'], colors='red')
    ax.clabel(cs3, cs3.levels, inline=True, fontsize=20)

    st.pyplot(fig)


    #st.write(data)
    #print(data)

# problematic section -> trying to split system into yt and zr if either one is 0
if cn>0 and ccl>0:
    plot_contour_yt(ccl)
    plot_contour_zr(cn)
    plot_contour_ysz(ck)
# elif cn==0 and ccl>0:
#     #zirconia system
#     plot_contour_zr(cn)
#     plot_contour_ysz(ck)
# elif ccl==0 and cn>0:
#     #yttria system
#     plot_contour_yt(ccl)
#     plot_contour_ysz(ck)