# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 09:03:04 2021

@author: Jacques
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

# =============================================================================
# IMPORTANT!!
# P in Pascal
# T in Kelvin 
# U.H in J/kg 
# p(rho) in kg/m^3
# =============================================================================

#Specify parameters that we know:

Vtot = 0.1 #m^#
Vliq_tot = 0.1*Vtot     #m3 10@of the total volume
Vvap_tot = Vtot - Vliq_tot
P_beg = 2*100*1000 #2 bar in Pa

# =============================================================================
# Use CoolProp to find values of saturated liquid and vapour 
# =============================================================================

#Densities of V and L
rho_1 = PropsSI('D','P',P_beg, 'Q', 0, 'water')  #('What we want', 'What we have', 'What we have value'
 #For liquid Quality =0                                                   #,'What we have', 'What we have value','Fluid name')
print(f'Desntiy of water liquid at 2 bar = {np.round(rho_1,0)} lg/m3')

rho_2 = PropsSI('D','P',P_beg, 'Q', 1, 'water')
#For liquid Quality =1
print(f'Desntiy of water vapour at 2 bar = {np.round(rho_2,4)} lg/m3')

#Now the mass of water and vapoutr can be found 

m_1 = rho_1*Vliq_tot
print(f'The initial liquid mass is m_1= {np.round(m_1,4)}kg') 

m_v=rho_2*Vvap_tot
print(f'The initial vapour mass is m_v= {np.round(m_v,4)}kg')
m_tot = m_1+m_v 
print(f'The total mass is m_tot= {np.round(m_tot,4)}kg')  

#Finding the U of liq and vap at beginning

u_1 = PropsSI('U','P',P_beg,'Q',0,'water')
print(f'The internal energy of the sat_liq(2bar) = {np.round(u_1/1000,2)} kJ/kg')

u_v = PropsSI('U','P',P_beg,'Q',1,'water')
print(f'The internal energy of the sat_vap(2bar) = {np.round(u_v/1000,2)} kJ/kg')

u_tot_beg= m_1*u_1 + m_v*u_v
print(f'The total internal energy initially = {np.round(u_tot_beg/1000,2)} kJ')

#Calculate density after superheating to fix system and find Uend
rho_end = m_tot/Vtot
u_end = PropsSI('U','D',rho_end,'Q',1,'water') #U at only vapour after heating 
print(f'The internal energy of the sat_vap after heating = {np.round(u_end/1000,2)} kJ/kg')
u_tot_end = m_tot*u_end 
print(f'The total internal energy of steam after heating = {np.round(u_tot_end/1000,2)} kJ')

#Now we can finally calculate Q (Q=delU)
Q= u_tot_end-u_tot_beg
print(f'THE AMOUNT OF HEAT FLOW INTO THE SYSTEM TO HAVE STATE CHANGE= {np.round(Q/1000,2)} kJ')

P_end = PropsSI('P','D', rho_end, 'Q', 1, 'water')
print(f'The pressure at the end of state change Pend= {np.round(P_end/1000,2)} kPa')

T_end = PropsSI('T','D',rho_end,'Q', 1, 'water')
print(f'The temperature at the end of phase change= {np.round(T_end,2)} K')
