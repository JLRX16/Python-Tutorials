# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 11:33:08 2021

@author: Jacques
"""
import scipy.integrate as it
import numpy as np
import scipy.optimize as opt
tb=298      #   K
te=330      #   K               #Random values
n1=1        #   kmol
n2=2        #   kmol
Cv= 30      #   kJ/kmol
R=8.314     #   kJ/kmol.K
Cp= Cv+R    #   kj/kmol         #Cp - Cv = R
gamma= Cp/Cv 
vb=14       #   m3
ve=25       #   m3
dHrxn=-20   #   kJ/kmol
pb=100      #   kPa
pe= 220     #   kPa
# =============================================================================
# Delta U integration (nCvdT)
# =============================================================================

def delU1(x):
    return n1*Cv
del_u1= it.quad(delU1,tb,te)[0]
print(f'delU1= {np.round(del_u1,4)} kJ/mol')

# =============================================================================
# Delta H integration   (nCpdT)
# =============================================================================
def delH1(x):
    return n1*Cp
del_h1 = it.quad(delH1, tb, te)[0]
print(f'delH1= {np.round(del_h1,4)} kJ/mol')

# =============================================================================
# Work integration
# =============================================================================

def W2(v):
    return n1*R*tb/v # =Pbound for a specific case this can change Pb

W2 = -(it.quad(W2,vb,ve )[0])
print(f'W2= {np.round(W2,4)} kJ/mol')

# =============================================================================
# Find values at adiabatic and Reversible
# =============================================================================

print(f'Gamma ={gamma}')
   #V for T known
def v2func(v2):
    return tb*((vb/v2)**(gamma-1))-te              

v2 = opt.fsolve(v2func, 16 )[0]
print(f'V2 = {np.round(v2,4)} m3/mol')


def v1func(v1):
    return tb*((v1/ve)**(gamma-1))-te  

v1 = opt.fsolve(v1func,16)[0]
print(f'V1 = {np.round(v1,4)} m3/mol')

#P for T known
def p2func(p2):
    return tb*((p2/pb)**((gamma-1)/gamma))-te  

p2 = opt.fsolve(p2func,16)[0]
print(f'P2 = {np.round(p2,4)} kPa')


def p1func(p1):
    return tb*((pe/p1)**((gamma-1))/gamma)-te  

p1 = opt.fsolve(p1func,16)[0]
print(f'P1 = {np.round(p1,4)} kPa')

#P for V known 

def p_2func(p_2):
    return (pb*(vb**gamma))/(ve**gamma)-p_2
p_2 = opt.fsolve(p_2func, 16)[0]
print(f'P2 = {np.round(p_2,4)} kPa')


def p_1func(p_1):
    return (p_1*(vb**gamma))/(ve**gamma)-pe
p_1 = opt.fsolve(p_1func, 16)[0]

#V for P known
def v_1func(v_1):
    return (pb*(v_1**gamma))/(ve**gamma)-pe
v_1 = opt.fsolve(v_1func, 16)[0]

print(f'V1 = {np.round(v_1,4)} m3')


# Issue at the moment
# def v_2func(v_2):
#     return (pb*(vb**gamma))/(v_2**gamma)-pe
# v_2 = opt.fsolve(v_2func, 16)[0]

# print(f'V2 = {np.round(v_2,4)} m3')


#T for P known
def T2func(t2):
    return tb*((pe/pb)**((gamma-1)/gamma))-t2  

t2 = opt.fsolve(T2func,300)[0]
print(f'T2 = {np.round(t2,4)} K')

def T1func(t1):
    return t1*((pe/pb)**((gamma-1)/gamma))-te 

t1 = opt.fsolve(T1func,300)[0]
print(f'T1 = {np.round(t1,4)} K')

#T for V known 
def T_2func(t_2):
    return tb*((vb/ve)**(gamma-1))-t_2              

t_2 = opt.fsolve(T_2func, 16 )[0]
print(f'T_2 = {np.round(t_2,4)} K')

def T_1func(t_1):
    return t_1*((vb/ve)**(gamma-1))-te              

t_1 = opt.fsolve(T_1func, 16 )[0]
print(f'T_1 = {np.round(t_1,4)} K')

# =============================================================================
# Urxn from Hrxn  (dHrxn =dUrxn + d(PV)rxn) PV =nRT tf d(PV) = RT(dn) 
# 
# NOTE: dU = dUrxn + nCvdT
#       dH = dHrxn + nCpdT
# =============================================================================

Urxn = dHrxn-(R*tb*(n2-n1))          #T1 specified as condition of dHrxn and is CONSTANT THROUGHOUT REACTION
print(f'delUrxn = {np.round(Urxn,4)} kJ')
