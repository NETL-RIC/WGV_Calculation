import numpy as np
from scipy.optimize import newton
import sys
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox

def preos(Tc, T, Pc, P, omega):
    """
    Peng-Robinson equation of state (PREOS)
    :param compound: string compound of interest - hydrogen or methane
    :param Tc: float critical temperature in Kelvin
    :param T: float temperature in Kelvin
    :param Pc: float critical pressure in bar
    :param P: float pressure in bar
    :param omega: float accentric factor
    Returns the density of a gas at the given T and P.
    """

    # PREOS
    Tr = T / Tc
    a = 0.457235 * R ** 2 * Tc ** 2 / Pc
    b = 0.0777961 * R * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
    alpha = (1 + kappa * (1 - np.sqrt(Tr))) ** 2

    A = a * alpha * P / R ** 2 / T ** 2
    B = b * P / R / T

    # build cubic polynomial
    def g(z):
        """
        Cubic polynomial in z from EOS. This should be zero.
        :param z: float compressibility factor
        """
        return z ** 3 - (1 - B) * z ** 2 + (A - 2 * B - 3 * B ** 2) * z - (
                A * B - B ** 2 - B ** 3)

    # Solve cubic polynomial for the compressibility factor
    z = newton(g, 1.0)  # compressibility factor
    rho = P / (R * T * z)  # density (mol/m3)

    return [rho, z]

def convert_output(wgm_H2, wgm_CH4, wgm, wge_H2, wge_CH4, wge, wgv_H2, wgv_CH4, wgv, units):
    """
    Output conversion function
    :param wgm_H2: calculated working gas mass for H2 in kg
    :param wgm_CH4: calculated working gas mass for CH4 in kg
    :param wgm: calculated total working gas mass in kg
    :param wge_H2: calculated working gas energy for H2 in TWh
    :param wge_CH4: calculated working gas energy for CH4 in TWh
    :param wge: calculated total working gas energy in TWh
    :param wgv_H2: calculated working gas volume for H2 in m3
    :param wgv_CH4: calculated working gas volume for CH4 in m3
    :param wgv: calculated total working gas volume in m3
    Returns a dictionary of outputs with labels and units
    """

    # Conversion factors for output
    m3_to_mmcf = 0.00003531
    kg_to_tons = 0.001102
    twh_to_tj = 3600  # tera joules in a tera-watt hour
    kg_to_MT = 1E-9  # kg to metric megatons
    m3_to_bcm = 1E-9  # m3 to billion cubic meters

    # Create outputs
    if units == 'Imperial':
        # Convert outputs to SI
        wgm_H2 = wgm_H2 * kg_to_tons
        wgm_CH4 = wgm_CH4 * kg_to_tons
        wgm = wgm * kg_to_tons

        wgv_H2 = wgv_H2 * m3_to_mmcf
        wgv_CH4 = wgv_CH4 * m3_to_mmcf
        wgv = wgv * m3_to_mmcf

        output = {'Working Gas Mass': wgm, 'Working Gas Mass of Hydrogen': wgm_H2,
                  'Working Gas Mass of Methane': wgm_CH4,
                  'Working Gas Energy': wge, 'Working Gas Energy of Hydrogen': wge_H2,
                  'Working Gas Energy of Methane': wge_CH4,
                  'Working Gas Volume': wgv, 'Working Gas Volume of Hydrogen': wgv_H2,
                  'Working Gas Volume of Methane': wgv_CH4}

    elif units == 'SI':
        wgm_H2 = wgm_H2 * kg_to_MT
        wgm_CH4 = wgm_CH4 * kg_to_MT
        wgm = wgm * kg_to_MT

        wge_H2 = wge_H2 * twh_to_tj
        wge_CH4 = wge_CH4 * twh_to_tj
        wge = wge * twh_to_tj

        wgv_H2 = wgv_H2 * m3_to_bcm
        wgv_CH4 = wgv_CH4 * m3_to_bcm
        wgv = wgv * m3_to_bcm

        output = {'Working Gas Mass': wgm, 'Working Gas Mass of Hydrogen': wgm_H2,
                  'Working Gas Mass of Methane': wgm_CH4,
                  'Working Gas Energy': wge, 'Working Gas Energy of Hydrogen': wge_H2,
                  'Working Gas Energy of Methane': wge_CH4,
                  'Working Gas Volume': wgv, 'Working Gas Volume of Hydrogen': wgv_H2,
                  'Working Gas Volume of Methane': wgv_CH4}
    return output

def wge_calc_ngsf(T_res, P_res, wgv_in, H2_frac_stp, units):
    """
    Volumetric calculation of the working gas energy for a natural gas storage site
    :param T_res: float temperature in the reservoir (C or F)
    :param P_res: float pressure in the reservoir at storage conditions (kPa or psi)
    :param wgv_in: float working gas in (m3 or MCF)
    :param hyd_frac: float fraction of hydrogen in gas mixture at STP
    :param units: string options: SI or Imperial - sets units for input and output
    Returns the total working gas energy (TJ or  or , hydrogen working gas energy, methane working gas energy, total working gas volume, hydrogen working gas volume, and methane working gas volume of the storage site
    """
    # Conversion factors:
    mcf_to_m3 = 28.32
    F_to_K = 255.9
    psi_to_bar = 0.06895
    kPa_to_bar = 0.01

    if units == 'Imperial':
        # Convert inputs to SI
        T_res = T_res * F_to_K
        P_res = P_res * psi_to_bar
        wgv_in = wgv_in * mcf_to_m3

    elif units == 'SI':
        P_res = P_res * kPa_to_bar
        T_res = T_res + 273.15  # C to K

    # Hydrogen properties
    Tc_H2 = 33.18  # K
    Pc_H2 = 13  # bar
    omega_H2 = -0.220
    mw_H2 = 0.002016  # kg/mol
    lhv_H2 = 3.332e-8  # TWh/kg

    # Methane properties
    Tc_CH4 = 190.6  # K
    Pc_CH4 = 46.1  # bar
    omega_CH4 = 0.011
    mw_CH4 = 0.01604  # kg/mol
    lhv_CH4 = 1.389e-8  # TWh kg

    # STP conditions
    P_stp = 1  # bar
    T_stp = 273.15  # K

    # Fraction of CH4 in gas
    CH4_frac_stp = 1 - H2_frac_stp

    # Density of H2 and CH4 at STP in kg/m^3
    rho_H2_stp = preos(Tc_H2, T_stp, Pc_H2, P_stp, omega_H2)[0] * mw_H2
    rho_CH4_stp = preos(Tc_CH4, T_stp, Pc_CH4, P_stp, omega_CH4)[0] * mw_CH4

    # Density of H2 and CH4 at reservoir conditions
    rho_H2_res = preos(Tc_H2, T_res, Pc_H2, P_res, omega_H2)[0] * mw_H2
    rho_CH4_res = preos(Tc_CH4, T_res, Pc_CH4, P_res, omega_CH4)[0] * mw_CH4

    # Hydrogen fraction in gas at reservoir conditions
    if H2_frac_stp == 1:
        H2_frac_res = H2_frac_stp
        CH4_frac_res = CH4_frac_stp
    else:
        H2_frac_res = ((rho_H2_stp / rho_H2_res) * H2_frac_stp) / (
                ((rho_H2_stp / rho_H2_res) * H2_frac_stp) + ((rho_CH4_stp / rho_CH4_res) * CH4_frac_stp))
        CH4_frac_res = 1 - H2_frac_res

    # Calculate working gas mass - results in kg
    wgm_H2 = (rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in)
    wgm_CH4 = (rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in)
    wgm = wgm_H2 + wgm_CH4

    # Calculate working gas energy - results in TWh
    wge_H2 = lhv_H2 * rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in
    wge_CH4 = lhv_CH4 * rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in
    wge = wge_H2 + wge_CH4

    # Calculate working gas volume - results in m3
    wgv_H2 = ((rho_H2_res / rho_H2_stp) * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in)
    wgv_CH4 = ((rho_CH4_res / rho_CH4_stp) * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in)
    wgv = wgv_H2 + wgv_CH4

    
    output = convert_output(wgm_H2, wgm_CH4, wgm, wge_H2, wge_CH4, wge, wgv_H2, wgv_CH4, wgv, units)

    return output

def wge_calc_dgf(T_res, P_res, V_res, cg_frac, H2_frac_stp, units):
    """
    Volumetric calculation of the working gas energy for a depleted gas reservoir
    :param T_res: float temperature in the reservoir (C or F)
    :param P_res: float pressure in the reservoir at storage conditions (kPa or psi)
    :param V_res: float volume of gas produced from field in (m3 or MCF)
    :param cg_frac: float fraction of reservoir gas used for cushion gas
    :param hyd_frac: float fraction of hydrogen in gas mixture at STP
    :param units: string options: SI or Imperial - sets units for input and output  
    Returns the total working gas energy, hydrogen working gas energy, methane working gas energy, total working gas volume, hydrogen working gas volume, and methane working gas volume of the storage site
    """

    #Conversion factors:
    mcf_to_m3 = 28.32
    F_to_K = 255.9 
    psi_to_bar = 0.06895
    kPa_to_bar = 0.01

    if units == 'Imperial':
        #Convert inputs to SI
        T_res = T_res*F_to_K
        P_res = P_res*psi_to_bar
        wgv_in = V_res * (1 - cg_frac) * mcf_to_m3

    elif units == 'SI':
        P_res = P_res * kPa_to_bar
        T_res = T_res + 273.15 #C to K
        wgv_in = V_res * (1 - cg_frac) 

    # Hydrogen properties
    Tc_H2 = 33.18  # K
    Pc_H2 = 13  # bar
    omega_H2 = -0.220
    mw_H2 = 0.002016  # kg/mol
    lhv_H2 = 3.332e-8  # TWh/kg

    # Methane properties
    Tc_CH4 = 190.6  # K
    Pc_CH4 = 46.1  # bar
    omega_CH4 = 0.011
    mw_CH4 = 0.01604  # kg/mol
    lhv_CH4 = 1.389e-8  # TWh kg

    # STP conditions
    P_stp = 1  # bar
    T_stp = 273.15  # K

    # Fraction of CH4 in gas
    CH4_frac_stp = 1 - H2_frac_stp

    # Density of H2 and CH4 at STP in kg/m^3
    rho_H2_stp = preos(Tc_H2, T_stp, Pc_H2, P_stp, omega_H2)[0] * mw_H2
    rho_CH4_stp = preos(Tc_CH4, T_stp, Pc_CH4, P_stp, omega_CH4)[0] * mw_CH4

    # Density of H2 and CH4 at reservoir conditions
    rho_H2_res = preos(Tc_H2, T_res, Pc_H2, P_res, omega_H2)[0] * mw_H2
    rho_CH4_res = preos(Tc_CH4, T_res, Pc_CH4, P_res, omega_CH4)[0] * mw_CH4

    # Hydrogen fraction in gas at reservoir conditions
    if H2_frac_stp == 1:
        H2_frac_res = H2_frac_stp
        CH4_frac_res = CH4_frac_stp
    else:
        H2_frac_res = ((rho_H2_stp / rho_H2_res) * H2_frac_stp) / (
                    ((rho_H2_stp / rho_H2_res) * H2_frac_stp) + ((rho_CH4_stp / rho_CH4_res) * CH4_frac_stp))
        CH4_frac_res = 1 - H2_frac_res

    #Calculate working gas mass - results in kg
    wgm_H2 = (rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in) 
    wgm_CH4 = (rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in) 
    wgm = wgm_H2 + wgm_CH4

    #Calculate working gas energy - results in TWh
    wge_H2 = lhv_H2 * rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in
    wge_CH4 = lhv_CH4 * rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in
    wge = wge_H2 + wge_CH4

    #Calculate working gas volume - results in m3
    wgv_H2 = ((rho_H2_res / rho_H2_stp) * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in) 
    wgv_CH4 = ((rho_CH4_res / rho_CH4_stp) * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in)
    wgv = wgv_H2 + wgv_CH4


    output = convert_output(wgm_H2, wgm_CH4, wgm, wge_H2, wge_CH4, wge, wgv_H2, wgv_CH4, wgv, units)


    return output

def wge_calc_saq(T_res, P_res, A_res, H_res, phi_res, ef, cg_frac, H2_frac_stp, units):
    """
    Volumetric calculation of the working gas energy for a saline aquifer
    :param T_res: float temperature in the reservoir (C or F)
    :param P_res: float pressure in the reservoir at storage conditions (kPa or psi)
    :param A_res: float reservoir area (m2 or ft^2)
    :param H_res: float reservoir thickness (ft or m)
    :param phi_res: float reservoir porosity 
    :param ef: float efficiency factor for storage in reservoir
    :param cg_frac: float fraction of reservoir gas used for cushion gas
    :param hyd_frac_stp: float fraction of hydrogen in gas mixture at STP
    :param units: string options: SI or Imperial - sets units for input and output  
    Returns the total working gas energy, hydrogen working gas energy, methane working gas energy, total working gas volume, hydrogen working gas volume, and methane working gas volume of the storage site
    """

    #Conversion factors:
    mcf_to_m3 = 28.32
    F_to_K = 255.9 
    psi_to_bar = 0.06895
    kPa_to_bar = 0.01
    ft2_to_m2 = 0.0929
    ft_to_m = 0.3048 

    if units == 'Imperial':
        #Convert inputs to SI
        T_res = T_res * F_to_K
        P_res = P_res * psi_to_bar
        A_res = A_res * ft2_to_m2
        H_res = H_res * ft_to_m
        V_res = A_res * H_res * phi_res * ef
        wgv_in = V_res * (1 - cg_frac) 

    elif units == 'SI':
        P_res = P_res * kPa_to_bar
        T_res = T_res + 273.15 #C to K   
        V_res = A_res * H_res * phi_res * ef
        wgv_in = V_res * (1 - cg_frac) 

    # Hydrogen properties
    Tc_H2 = 33.18  # K
    Pc_H2 = 13  # bar
    omega_H2 = -0.220
    mw_H2 = 0.002016  # kg/mol
    lhv_H2 = 3.332e-8  # TWh/kg

    # Methane properties
    Tc_CH4 = 190.6  # K
    Pc_CH4 = 46.1  # bar
    omega_CH4 = 0.011
    mw_CH4 = 0.01604  # kg/mol
    lhv_CH4 = 1.389e-8  # TWh kg

    # STP conditions
    P_stp = 1  # bar
    T_stp = 273.15  # K

    # Fraction of CH4 in gas
    CH4_frac_stp = 1 - H2_frac_stp

    # Density of H2 and CH4 at STP in kg/m^3
    rho_H2_stp = preos(Tc_H2, T_stp, Pc_H2, P_stp, omega_H2)[0] * mw_H2
    rho_CH4_stp = preos(Tc_CH4, T_stp, Pc_CH4, P_stp, omega_CH4)[0] * mw_CH4

    # Density of H2 and CH4 at reservoir conditions
    rho_H2_res = preos(Tc_H2, T_res, Pc_H2, P_res, omega_H2)[0] * mw_H2
    rho_CH4_res = preos(Tc_CH4, T_res, Pc_CH4, P_res, omega_CH4)[0] * mw_CH4

    # Hydrogen fraction in gas at reservoir conditions
    if H2_frac_stp == 1:
        H2_frac_res = H2_frac_stp
        CH4_frac_res = CH4_frac_stp
    else:
        H2_frac_res = ((rho_H2_stp / rho_H2_res) * H2_frac_stp) / (
                    ((rho_H2_stp / rho_H2_res) * H2_frac_stp) + ((rho_CH4_stp / rho_CH4_res) * CH4_frac_stp))
        CH4_frac_res = 1 - H2_frac_res

    #Calculate working gas mass - results in kg
    wgm_H2 = (rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in) 
    wgm_CH4 = (rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in) 
    wgm = wgm_H2 + wgm_CH4

    #Calculate working gas energy - results in TWh
    wge_H2 = lhv_H2 * rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in
    wge_CH4 = lhv_CH4 * rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in
    wge = wge_H2 + wge_CH4

    #Calculate working gas volume - results in m3
    wgv_H2 = ((rho_H2_res / rho_H2_stp) * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in) 
    wgv_CH4 = ((rho_CH4_res / rho_CH4_stp) * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv_in)
    wgv = wgv_H2 + wgv_CH4

    output = convert_output(wgm_H2, wgm_CH4, wgm, wge_H2, wge_CH4, wge, wgv_H2, wgv_CH4, wgv, units)

    return output

def wge_calc_salt_spec(T_res, H_cav, D_cav, h_n, h_o, g_min, cg_frac, units):
    """
    Volumetric calculation of the working gas energy for a salt cavern with specified dimensions
    :param T_res: float temperature in the cavern (C or F)
    :param H_cav: float height of the cavern (m or ft)
    :param D_cav: float diameter of the cavern (m or ft)
    :param h_n: float depth to top of cavern neck (m or ft)
    :param h_c: float depth to center of cavern (m or ft)
    :param h_o: float depth of cavern center which can be emptied to zero pressure (m or ft)
    :param g_min: float gradient of minimum storage pressure (based on rock salt strength and stress distribution in cavern) (MPa/m or psi/ft)
    :param cg_frac: float cushion gas fraction in cavern
    :param units: string options: SI or Imperial - sets units for input and output  
    Returns the total working gas energy, total working gas volume, and working gas mass for a specific salt cavern storage site

    """

    #Conversion factors:
    F_to_K = 255.9 
    ft_to_m = 0.3048 
    psift_to_mpam = 0.02262

    if units == 'Imperial':
        T_res = T_res * F_to_K
        H_cav = H_cav * ft_to_m
        D_cav = D_cav * ft_to_m
        h_n = h_n * ft_to_m

        h_o = h_o * ft_to_m
        g_min = g_min * psift_to_mpam

    elif units == 'SI':
        T_res = T_res + 273.15 #C to K

    # Conversion factors
    lhv_H2 = 3.332e-8  # TWh/kg

    # Constants
    # Fracture gradient
    g_f = 0.016  # MPa/m

    # Hydrogen properties
    Tc_H2 = 33.18  # K
    Pc_H2 = 13  # bar
    omega_H2 = -0.220
    mw_H2 = 0.002016  # kg/mol
    lhv_H2 = 3.332e-8  # TWh/kg
    
    # Individual gas constant for hydrogen (J/kg K)
    R_h2 = 4121.73  # K/kg K
    rho_h2_n = 0.089  # kg/m3 Density of H2 at normal conditions
    
    # calculate depth to center of cavern
    h_c = h_n + H_cav/2

    # Volume of cavern
    V_cav = (np.pi/ 12) * (D_cav ** 2)*(3 * H_cav - D_cav)  # m^3

    # Maximum pressure
    P_max = g_f * h_n

    # Minimum pressure
    P_min = g_min * (h_c - h_o)  # MPa

    # Calculate compressibilities
    Z_max = preos(Tc_H2, T_res, Pc_H2, P_max, omega_H2)[1]
    Z_min = preos(Tc_H2, T_res, Pc_H2, P_min, omega_H2)[1]

    # Maximum and minimum mass stored
    M_max = (1 - cg_frac) * ((P_max * V_cav) / R_h2 * T_res * Z_max)
    M_min = (1 - cg_frac) * ((P_min * V_cav) / R_h2 * T_res * Z_min)

    # Working gas mass
    wgm = M_max - M_min  # kg

    # Working gas energy (TWh)
    wge = wgm * lhv_H2

    # Working gas volume (m3)
    wgv = (wgm / rho_h2_n) 

    output = convert_output(wgm, 0 , wgm, wge, 0, wge, wgv, 0, wgv, units)

    return output

def wge_calc_salt_max(T_res, F_thick, F_area, h_n, h_o, g_min, cg_frac, units):
    """
    Volumetric calculation of the working gas energy for maximum number of salt caverns fitted in an area
    :param T_res: float temperature in the cavern (C or F)
    :param F_thick: float formation thickness in (m or ft)
    :param F_area: float total area of the formation (hectare or acres)  
    :param h_n: float depth to top of cavern neck (m or ft)
    :param h_c: float depth to center of cavern (m or ft)
    :param h_o: float depth of cavern center which can be emptied to zero pressure (m or ft)
    :param g_min: float gradient of minimum storage pressure (based on rock salt strength and stress distribution in cavern) (MPa/m or psi/ft)
    :param cg_frac: float cushion gas fraction in cavern
    :param units: string options: SI or Imperial - sets units for input and output  
    Returns the total working gas energy, total working gas volume, and working gas mass for a multiple salt cavern storage sites over an entire formation
    """

    #Conversion factors:
    F_to_K = 255.9
    kg_to_tons = 0.001102
    m3_to_mmcf = 0.00003531
    ft_to_m = 0.3048 
    psift_to_mpam = 0.02262
    acre_to_m2 = 4047
    ft2_to_acre = 43560
    ha_to_m2 = 10000
    lhv_H2 = 3.332e-8  # TWh/kg

    if units == 'Imperial':
        T_res = T_res * F_to_K
        F_thick = F_thick * ft_to_m
        F_area = F_area * acre_to_m2
        h_n = h_n * ft_to_m

        h_o = h_o * ft_to_m
        g_min = g_min * psift_to_mpam

    elif units == 'SI':
        T_res = T_res + 273.15 #C to K
        F_area = F_area * ha_to_m2

    # Cavern density in space
    cav_den = 0.0000185 #caverns/m**2 (1 cavern/5.4 hectare) - Lankof and Tarkowski

    # Constants
    # Fracture gradient
    g_f = 0.016  # MPa/m

    # Hydrogen properties
    Tc_H2 = 33.18  # K
    Pc_H2 = 13  # bar
    omega_H2 = -0.220
    mw_H2 = 0.002016  # kg/mol
    lhv_H2 = 3.332e-8  # TWh/kg
    # Individual gas constant for hydrogen (J/kg K)
    R_h2 = 4121.73  # K/kg K
    rho_h2_n = 0.089  # kg/m3 Density of H2 at normal conditions
    


    # First check to see if the formation is thick enough
    if unit=='SI' and F_thick <= 65:
        tkinter.messagebox.showinfo(title='Warning',
                                    message='Formation thickness must be greater than 65 m to contain a storage cavern')

    if unit=='Imperial' and F_thick <= 213:
        tkinter.messagebox.showinfo(title='Warning',
                                    message='Formation thickness must be greater than 213 ft to contain a storage cavern')        
        

    H_cav = F_thick - 65
    D_cav = H_cav * (2/3)
    
    # calculate depth to center of cavern
    h_c = h_n + H_cav/2
    
    
    if unit=='SI' and F_area < 16 * D_cav**2:
            tkinter.messagebox.showinfo(title='Warning',
                            message="Based on the formation thickness, area must be greater than {0:.2f} ha to contain multiple storage caverns".format(16 * D_cav**2/ha_to_m2))

    if unit=='Imperial' and F_area < 16 * D_cav**2:
            tkinter.messagebox.showinfo(title='Warning',
                            message="Based on the formation thickness, area must be greater than {0:.2f} acre to contain multiple storage caverns".format(16 * D_cav**2/acre_to_m2))
        

    # Volume of cavern
    V_cav = (np.pi / 12) * (D_cav ** 2)*(3 * H_cav - D_cav)  # m^3

    # Maximum pressure
    P_max = g_f * h_n

    # Minimum pressure
    P_min = g_min * (h_c - h_o)  # MPa

    # Calculate compressibilities
    Z_max = preos(Tc_H2, T_res, Pc_H2, P_max, omega_H2)[1]
    Z_min = preos(Tc_H2, T_res, Pc_H2, P_min, omega_H2)[1]

    # Maximum and minimum mass stored
    M_max = (1 - cg_frac) * ((P_max * V_cav) / R_h2 * T_res * Z_max)
    M_min = (1 - cg_frac) * ((P_min * V_cav) / R_h2 * T_res * Z_min)

    # Working gas mass per cavern
    wgm_pc = M_max - M_min  # kg

    wgm_tons_pc = wgm_pc * kg_to_tons  # tons (short tons)

    # Working gas energy (TWh)
    wge_pc = wgm_pc * lhv_H2

    # Working gas volume (MMCF)
    wgv_pc = (wgm_pc / rho_h2_n) * m3_to_mmcf

    #Number of caverns in formation
    n_cav = F_area*cav_den

    #Calculate total energy, volume, and mass for formation
    wgm = wgm_pc * n_cav
    wge = wge_pc * n_cav
    wgv = wgv_pc * n_cav

    output = convert_output(wgm, 0 , wgm, wge, 0, wge, wgv, 0, wgv, units)

    return output

def calc_run():
    
    """
    Function to calculate working gas volume, mass, and energy and populate the final results   
    """
    click = clicked.get()
    
    if click == "Natural Gas Storage Site":
        results = wge_calc_ngsf(float(t_entry.get()), float(p_entry.get()),
                                    float(wvg_entry.get()), float(hf_entry.get()), unit)
    elif click == "Depleted Gas Field":    
        results = wge_calc_dgf(float(t_entry.get()), float(p_entry.get()),
                                    float(cum_entry.get()), float(cu_entry.get()), float(hf_entry.get()), unit)
    elif click == "Saline Aquifer":
        results = wge_calc_saq(float(t_entry.get()), float(p_entry.get()),
                               float(ra_entry.get()), float(rt_entry.get()),
                               float(po_entry.get()), float(ef_entry.get()),
                               float(cu_entry.get()), float(hf_entry.get()), unit)
    
    elif ((click == "Salt Cavern")&(salt==0)):       
        results = wge_calc_salt_spec(float(t_entry.get()), float(h_entry.get()),
                                     float(di_entry.get()), float(dt_entry.get()),
                                     float(dz_entry.get()),
                                     float(gm_entry.get()), float(cu_entry.get()), unit)
    elif ((click == "Salt Cavern")&(salt==1)):       
        results = wge_calc_salt_max(float(t_entry.get()), float(ft_entry.get()),
                                     float(fa_entry.get()), float(dt_entry.get()),
                                     float(dz_entry.get()),
                                     float(gm_entry.get()), float(cu_entry.get()), unit)
        
 
    # results
    results_label = ttk.Label(botframe,text="Results:", font=('Helvetica 10 underline'))
    results_label.grid(row=1, column=0, sticky=tk.W, pady=5)   

    # working gas mass
    WGMass = ttk.Label(botframe, text="Working gas mass")
    WGMass.grid(row=2, column=2, sticky=tk.E, padx=20)

    hyd = ttk.Label(botframe,text="Hydrogen:")
    hyd.grid(row=3, column=1, sticky=tk.W, pady=2)
    h1_entry = ttk.Entry(botframe, width=10)
    h1_entry.grid(row=3, column=2, pady=2)
    
    met = ttk.Label(botframe,text="Methane:")
    met.grid(row=4, column=1, sticky=tk.W, pady=2)
    m1_entry = ttk.Entry(botframe, width=10)
    m1_entry.grid(row=4, column=2, pady=2)
    
    tot = ttk.Label(botframe,text="Total:")
    tot.grid(row=5, column=1, sticky=tk.W, pady=2)
    t1_entry = ttk.Entry(botframe, width=10)
    t1_entry.grid(row=5, column=2, pady=2)
    
    # working gas volume
    WGVol = ttk.Label(botframe, text="Working gas volume")
    WGVol.grid(row=2, column=3, sticky=tk.E,padx=20)

    h2_entry = ttk.Entry(botframe, width=10)
    h2_entry.grid(row=3, column=3, pady=2)

    m2_entry = ttk.Entry(botframe, width=10)
    m2_entry.grid(row=4, column=3, pady=2)

    t2_entry = ttk.Entry(botframe, width=10)
    t2_entry.grid(row=5, column=3, pady=2)
    
    # working gas energy
    WGVol = ttk.Label(botframe, text="Working gas energy")
    WGVol.grid(row=2, column=4, sticky=tk.E,padx=20)

    h3_entry = ttk.Entry(botframe, width=10)
    h3_entry.grid(row=3, column=4, pady=2)

    m3_entry = ttk.Entry(botframe, width=10)
    m3_entry.grid(row=4, column=4, pady=2)

    t3_entry = ttk.Entry(botframe, width=10)
    t3_entry.grid(row=5, column=4, pady=2)
    
    # set output units
    if unit == 'SI':
        WGMass = ttk.Label(botframe, text="Working gas mass (MT)")
        WGMass.grid(row=2, column=2, sticky=tk.E, padx=20)
        
        WGVol = ttk.Label(botframe, text="Working gas volume (BCM)")
        WGVol.grid(row=2, column=3, sticky=tk.E,padx=20)
        
        WGVol = ttk.Label(botframe, text="Working gas energy (TJ)")
        WGVol.grid(row=2, column=4, sticky=tk.E,padx=20)
        
    if unit == 'Imperial':
        WGMass = ttk.Label(botframe, text="Working gas mass (tones)")
        WGMass.grid(row=2, column=2, sticky=tk.E, padx=20)
        
        WGVol = ttk.Label(botframe, text="Working gas volume (MMCF)")
        WGVol.grid(row=2, column=3, sticky=tk.E,padx=20)
        
        WGVol = ttk.Label(botframe, text="Working gas energy (Twh)")
        WGVol.grid(row=2, column=4, sticky=tk.E,padx=20)   
    


    run_text.set("Running...")
    


    #results_dic = wge_calc_ngsf(float(t_entry.get()), float(p_entry.get()), float(wgv_entry.get()),float(h2_entry.get()))

    h1_entry.delete(0, tk.END)
    h2_entry.delete(0, tk.END)
    h3_entry.delete(0, tk.END)
    m1_entry.delete(0, tk.END)
    m2_entry.delete(0, tk.END)
    m3_entry.delete(0, tk.END)
    t1_entry.delete(0, tk.END)
    t2_entry.delete(0, tk.END)
    t3_entry.delete(0, tk.END)
    


    t1_entry.insert(0, round(results.get('Working Gas Mass'), 6))
    h1_entry.insert(0, round(results.get('Working Gas Mass of Hydrogen'), 6))
    m1_entry.insert(0, round(results.get('Working Gas Mass of Methane'), 6))
    
    t2_entry.insert(0, round(results.get('Working Gas Volume'), 6))
    h2_entry.insert(0, round(results.get('Working Gas Volume of Hydrogen'), 6))
    m2_entry.insert(0, round(results.get('Working Gas Volume of Methane'), 6))
    
    t3_entry.insert(0, round(results.get('Working Gas Energy'), 6))
    h3_entry.insert(0, round(results.get('Working Gas Energy of Hydrogen'), 6))
    m3_entry.insert(0, round(results.get('Working Gas Energy of Methane'), 6))
    
    

    run_text.set("Run")
    
    
    
    

    
def cleartopframe():
    for widget in topframe.winfo_children():
        widget.destroy()

def clearFrame():
    # destroy all widgets from frame
    for widget in frame.winfo_children():
        widget.destroy()
    for widget in botframe.winfo_children():
        widget.destroy()
def clearradio():
    for widget in radioframe.winfo_children():
        widget.destroy()  
        
def show(value):
    
    """
    Function to show the input values
    """
    click = clicked.get()
    global results
    results={}
    
    global t_entry, p_entry, wvg_entry, hf_entry
    global cum_entry, cu_entry
    global ra_entry, rt_entry, po_entry, ef_entry
    global h_entry, di_entry, dt_entry, dc_entry, dz_entry, gm_entry
 
    if click == "Natural Gas Storage Site":
        clearFrame()
        clearradio()
        
     
        inputParam = ttk.Label(frame,text="Input Parameters:", font=('Helvetica 10 underline'))
        inputParam.grid(row=2, column=0, sticky=tk.W, pady=5)   
        
        inputPres = ttk.Label(frame,text="Pressure:")
        inputPres.grid(row=3, column=0, sticky=tk.W)
        p_entry = ttk.Entry(frame, width=10)
        p_entry.grid(row=3, column=1, sticky=tk.E, pady=2, padx=5)
        

        inputTemp = ttk.Label(frame,text="Temperature:")
        inputTemp.grid(row=4, column=0, sticky=tk.W, pady=2)
        t_entry = ttk.Entry(frame, width=10)
        t_entry.grid(row=4, column=1, sticky=tk.E, pady=2, padx=5)
        

        inputWVG = ttk.Label(frame,text="Reported working gas volume:")
        inputWVG.grid(row=5, column=0, sticky=tk.W, pady=2)
        wvg_entry = ttk.Entry(frame, width=10)
        wvg_entry.grid(row=5, column=1, sticky=tk.E, pady=2, padx=5)
        

        inputH2Frac = ttk.Label(frame,text="Hydrogen fraction:")
        inputH2Frac.grid(row=6, column=0, sticky=tk.W, pady=2)
        hf_entry = ttk.Entry(frame, width=10)
        hf_entry.grid(row=6, column=1, sticky=tk.E, pady=2, padx=5)
        
        
        # set units
        if unit == 'SI':
            p_unit = ttk.Label(frame,text="KPa")
            p_unit.grid(row=3, column=2, sticky=tk.W, pady=2)
            p_entry.insert(0,76.92)
            
            t_unit = ttk.Label(frame,text="K")
            t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
            t_entry.insert(0,305.11)
            
            wvg_unit = ttk.Label(frame,text="m3")
            wvg_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
            wvg_entry.insert(0,2200000)
            
            hf_entry.insert(0,0.1)
            
            
        elif unit == 'Imperial':
            p_unit = ttk.Label(frame,text="psi")
            p_unit.grid(row=3, column=2, sticky=tk.W, pady=2)
            p_entry.insert(0,11.16)
            
            t_unit = ttk.Label(frame,text="F")
            t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
            t_entry.insert(0,89.53)
            
            wvg_unit = ttk.Label(frame,text="MCF")
            wvg_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
            wvg_entry.insert(0,77692.3)
            
            hf_entry.insert(0,0.1)


    if click == "Depleted Gas Field":
        clearFrame()
        clearradio()
                 
            
        
        inputParam = ttk.Label(frame,text="Input Parameters:", font=('Helvetica 10 underline'))
        inputParam.grid(row=2, column=0, sticky=tk.W, pady=5)   
        
        inputPres = ttk.Label(frame,text="Pressure:")
        inputPres.grid(row=3, column=0, sticky=tk.W)
        p_entry = ttk.Entry(frame, width=10)
        p_entry.grid(row=3, column=1, sticky=tk.E, pady=2, padx=5)

        inputTemp = ttk.Label(frame,text="Temperature:")
        inputTemp.grid(row=4, column=0, sticky=tk.W, pady=2)
        t_entry = ttk.Entry(frame, width=10)
        t_entry.grid(row=4, column=1, sticky=tk.E, pady=2, padx=5)

        inputCum = ttk.Label(frame,text="Cumulative produced gas:")
        inputCum.grid(row=5, column=0, sticky=tk.W, pady=2)
        cum_entry = ttk.Entry(frame, width=10)
        cum_entry.grid(row=5, column=1, sticky=tk.E, pady=2, padx=5)

        inputH2Frac = ttk.Label(frame,text="Hydrogen fraction:")
        inputH2Frac.grid(row=6, column=0, sticky=tk.W, pady=2)
        hf_entry = ttk.Entry(frame, width=10)
        hf_entry.grid(row=6, column=1, sticky=tk.E, pady=2, padx=5)
        
        inputCush = ttk.Label(frame,text="Cushion gas fraction:")
        inputCush.grid(row=7, column=0, sticky=tk.W, pady=2)
        cu_entry = ttk.Entry(frame, width=10)
        cu_entry.grid(row=7, column=1, sticky=tk.E, pady=2, padx=5)
        
        # set units
        if unit == 'SI':
            p_unit = ttk.Label(frame,text="KPa")
            p_unit.grid(row=3, column=2, sticky=tk.W, pady=2)
            p_entry.insert(0,76.92)
            
            t_unit = ttk.Label(frame,text="K")
            t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
            t_entry.insert(0,305.11)
            
            cum_unit = ttk.Label(frame,text="m3")
            cum_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
            cum_entry.insert(0,2200000)
            
            hf_entry.insert(0,0.1)
            cu_entry.insert(0,0.5)
            
            
            
        elif unit == 'Imperial':
            p_unit = ttk.Label(frame,text="psi")
            p_unit.grid(row=3, column=2, sticky=tk.W, pady=2)
            p_entry.insert(0,11.16)
            
            t_unit = ttk.Label(frame,text="F")
            t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
            t_entry.insert(0,89.53)
            
            cum_unit = ttk.Label(frame,text="MCF")
            cum_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
            cum_entry.insert(0,77692.3)
            
            hf_entry.insert(0,0.1)
            cu_entry.insert(0,0.5)

        
    
    if click == "Saline Aquifer":
        clearFrame()
        clearradio()
              
        
        inputParam = ttk.Label(frame,text="Input Parameters:", font=('Helvetica 10 underline'))
        inputParam.grid(row=2, column=0, sticky=tk.W, pady=5)   
        
        inputPres = ttk.Label(frame,text="Pressure:")
        inputPres.grid(row=3, column=0, sticky=tk.W)
        p_entry = ttk.Entry(frame, width=10)
        p_entry.grid(row=3, column=1, sticky=tk.E, pady=2, padx=5)

        inputTemp = ttk.Label(frame,text="Temperature:")
        inputTemp.grid(row=4, column=0, sticky=tk.W, pady=2)
        t_entry = ttk.Entry(frame, width=10)
        t_entry.grid(row=4, column=1, sticky=tk.E, pady=2, padx=5)

        inputResA = ttk.Label(frame,text="Reservoir area:")
        inputResA.grid(row=5, column=0, sticky=tk.W, pady=2)
        ra_entry = ttk.Entry(frame, width=10)
        ra_entry.grid(row=5, column=1, sticky=tk.E, pady=2, padx=5)

        inputResT = ttk.Label(frame,text="Reservoir thickness:")
        inputResT.grid(row=6, column=0, sticky=tk.W, pady=2)
        rt_entry = ttk.Entry(frame, width=10)
        rt_entry.grid(row=6, column=1, sticky=tk.E, pady=2, padx=5)
        
        inputPor = ttk.Label(frame,text="Porosity:")
        inputPor.grid(row=7, column=0, sticky=tk.W, pady=2)
        po_entry = ttk.Entry(frame, width=10)
        po_entry.grid(row=7, column=1, sticky=tk.E, pady=2, padx=5)
        
        inputEff = ttk.Label(frame,text="Efficiency factor:")
        inputEff.grid(row=3, column=2, sticky=tk.W, pady=2, padx=50)
        ef_entry = ttk.Entry(frame, width=10)
        ef_entry.grid(row=3, column=3, sticky=tk.E, pady=2, padx=5)
        
        inputH2Frac = ttk.Label(frame,text="Hydrogen fraction:")
        inputH2Frac.grid(row=4, column=2, sticky=tk.W, pady=2, padx=50)
        hf_entry = ttk.Entry(frame, width=10)
        hf_entry.grid(row=4, column=3, sticky=tk.E, pady=2, padx=5)
        
        inputCush = ttk.Label(frame,text="Cushion gas fraction:")
        inputCush.grid(row=5, column=2, sticky=tk.W, pady=2, padx=50)
        cu_entry = ttk.Entry(frame, width=10)
        cu_entry.grid(row=5, column=3, sticky=tk.E, pady=2, padx=5)
        
        
        # set units
        if unit == 'SI':
            p_unit = ttk.Label(frame,text="KPa")
            p_unit.grid(row=3, column=2, sticky=tk.W, pady=2)
            p_entry.insert(0,76.92)
            
            t_unit = ttk.Label(frame,text="K")
            t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
            t_entry.insert(0,305.11)
            
            ra_unit = ttk.Label(frame,text="m2")
            ra_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
            ra_entry.insert(0,100000)
            
            rt_unit = ttk.Label(frame,text="m")
            rt_unit.grid(row=6, column=2, sticky=tk.W, pady=2)
            rt_entry.insert(0,50)
            
            po_entry.insert(0,0.3)
            ef_entry.insert(0,0.1)
            hf_entry.insert(0,0.2)
            cu_entry.insert(0,0.8)
            
           
            
            
        elif unit == 'Imperial':
            p_unit = ttk.Label(frame,text="psi")
            p_unit.grid(row=3, column=2, sticky=tk.W, pady=2)
            p_entry.insert(0,11.16)
            
            t_unit = ttk.Label(frame,text="F")
            t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
            t_entry.insert(0,89.53)
            
            ra_unit = ttk.Label(frame,text="ft2")
            ra_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
            ra_entry.insert(0,1076391)
            
            rt_unit = ttk.Label(frame,text="ft")
            rt_unit.grid(row=6, column=2, sticky=tk.W, pady=2)
            rt_entry.insert(0,164)
            
            po_entry.insert(0,0.3)
            ef_entry.insert(0,0.1)
            hf_entry.insert(0,0.2)
            cu_entry.insert(0,0.8)
            
            
            
        
    if click == "Salt Cavern":
        clearFrame()
        global radioVar, unitVar
        
        
        def show_single():
            clearFrame()
            global salt
            salt = 0
            
            global t_entry, h_entry, di_entry, dt_entry, dz_entry, gm_entry,cu_entry
          
            
            
            
            inputParam = ttk.Label(frame,text="Input Parameters:", font=('Helvetica 10 underline'))
            inputParam.grid(row=3, column=0, sticky=tk.W, pady=5)   

            inputTemp = ttk.Label(frame,text="Temperature:")
            inputTemp.grid(row=4, column=0, sticky=tk.W)
            t_entry = ttk.Entry(frame, width=10)
            t_entry.grid(row=4, column=1, sticky=tk.W, pady=2, padx=5)

            inputhgt = ttk.Label(frame,text="Height of cavern:")
            inputhgt.grid(row=5, column=0, sticky=tk.W, pady=2)
            h_entry = ttk.Entry(frame, width=10)
            h_entry.grid(row=5, column=1, sticky=tk.W, pady=2, padx=5)

            inputDia = ttk.Label(frame,text="Diameter of cavern:")
            inputDia.grid(row=6, column=0, sticky=tk.W, pady=2)
            di_entry = ttk.Entry(frame, width=10)
            di_entry.grid(row=6, column=1, sticky=tk.W, pady=2, padx=5)

            inputdet = ttk.Label(frame,text="Depth to top of cavern:")
            inputdet.grid(row=7, column=0, sticky=tk.W, pady=2)
            dt_entry = ttk.Entry(frame, width=10)
            dt_entry.grid(row=7, column=1, sticky=tk.W, pady=2, padx=5)

            inputdcz = ttk.Label(frame,text="Depth of cavern center to zero pressure:")
            inputdcz.grid(row=4, column=2, sticky=tk.W, pady=2, padx=50)
            dz_entry = ttk.Entry(frame, width=10)
            dz_entry.grid(row=4, column=3, sticky=tk.W, pady=2, padx=5)

            inputGms = ttk.Label(frame,text="Gradient of minimum storage pressure:")
            inputGms.grid(row=5, column=2, sticky=tk.W, pady=2, padx=50)
            gm_entry = ttk.Entry(frame, width=10)
            gm_entry.grid(row=5, column=3, sticky=tk.W, pady=2, padx=5)
            

            inputCush = ttk.Label(frame,text="Cushion gas fraction:")
            inputCush.grid(row=6, column=2, sticky=tk.W, pady=2, padx=50)
            cu_entry = ttk.Entry(frame, width=10)
            cu_entry.grid(row=6, column=3, sticky=tk.W, pady=2, padx=5)
            
            
            
            
            # set units
            if unit == 'SI':
                t_unit = ttk.Label(frame,text="K")
                t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
                t_entry.insert(0,305.11)

                h_unit = ttk.Label(frame,text="m")
                h_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
                h_entry.insert(0,300)

                di_unit = ttk.Label(frame,text="m")
                di_unit.grid(row=6, column=2, sticky=tk.W, pady=2)
                di_entry.insert(0,50)

                dt_unit = ttk.Label(frame,text="m")
                dt_unit.grid(row=7, column=2, sticky=tk.W, pady=2)
                dt_entry.insert(0,1000)
                
#                 dc_unit = ttk.Label(frame,text="m")
#                 dc_unit.grid(row=8, column=2, sticky=tk.W, pady=2)
#                 dc_entry.insert(0,1150)
                
                dz_unit = ttk.Label(frame,text="m")
                dz_unit.grid(row=4, column=4, sticky=tk.W, pady=2)
                dz_entry.insert(0,1200)
                
                gm_unit = ttk.Label(frame,text="MPa/m")
                gm_unit.grid(row=5, column=4, sticky=tk.W, pady=2)
                gm_entry.insert(0,0.00835)
                
                cu_entry.insert(0,0.2)



            elif unit == 'Imperial':
                t_unit = ttk.Label(frame,text="F")
                t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
                t_entry.insert(0,89.53)

                h_unit = ttk.Label(frame,text="ft")
                h_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
                h_entry.insert(0,984)

                di_unit = ttk.Label(frame,text="ft")
                di_unit.grid(row=6, column=2, sticky=tk.W, pady=2)
                di_entry.insert(0,164)

                dt_unit = ttk.Label(frame,text="ft")
                dt_unit.grid(row=7, column=2, sticky=tk.W, pady=2)
                dt_entry.insert(0,3281)
                
                dz_unit = ttk.Label(frame,text="ft")
                dz_unit.grid(row=4, column=4, sticky=tk.W, pady=2)
                dz_entry.insert(0,3937)
                
                gm_unit = ttk.Label(frame,text="psi/ft")
                gm_unit.grid(row=5, column=4, sticky=tk.W, pady=2)
                gm_entry.insert(0,0.00835)
                
                cu_entry.insert(0,0.2)
            
            
            
        def show_multiple():
            clearFrame()
            global salt
            salt = 1
            
            global t_entry, ft_entry, fa_entry, dt_entry, dz_entry, gm_entry,cu_entry
            
 
            inputParam = ttk.Label(frame,text="Input Parameters:", font=('Helvetica 10 underline'))
            inputParam.grid(row=3, column=0, sticky=tk.W, pady=5)   

            inputTemp = ttk.Label(frame,text="Temperature:")
            inputTemp.grid(row=4, column=0, sticky=tk.W)
            t_entry = ttk.Entry(frame, width=10)
            t_entry.grid(row=4, column=1, sticky=tk.W, pady=2, padx=5)

            inputft = ttk.Label(frame,text="Formation thickness:")
            inputft.grid(row=5, column=0, sticky=tk.W, pady=2)
            ft_entry = ttk.Entry(frame, width=10)
            ft_entry.grid(row=5, column=1, sticky=tk.W, pady=2, padx=5)

            inputfa = ttk.Label(frame,text="Formation area:")
            inputfa.grid(row=6, column=0, sticky=tk.W, pady=2)
            fa_entry = ttk.Entry(frame, width=10)
            fa_entry.grid(row=6, column=1, sticky=tk.W, pady=2, padx=5)

            inputdet = ttk.Label(frame,text="Depth to top of cavern:")
            inputdet.grid(row=7, column=0, sticky=tk.W, pady=2)
            dt_entry = ttk.Entry(frame, width=10)
            dt_entry.grid(row=7, column=1, sticky=tk.W, pady=2, padx=5)

            inputdcz = ttk.Label(frame,text="Depth of cavern center to zero pressure:")
            inputdcz.grid(row=4, column=2, sticky=tk.W, pady=2, padx=50)
            dz_entry = ttk.Entry(frame, width=10)
            dz_entry.grid(row=4, column=3, sticky=tk.W, pady=2, padx=5)

            inputGms = ttk.Label(frame,text="Gradient of minimum storage pressure:")
            inputGms.grid(row=5, column=2, sticky=tk.W, pady=2, padx=50)
            gm_entry = ttk.Entry(frame, width=10)
            gm_entry.grid(row=5, column=3, sticky=tk.W, pady=2, padx=5)


            inputCush = ttk.Label(frame,text="Cushion gas fraction:")
            inputCush.grid(row=6, column=2, sticky=tk.W, pady=2, padx=50)
            cu_entry = ttk.Entry(frame, width=10)
            cu_entry.grid(row=6, column=3, sticky=tk.W, pady=2, padx=5)

            
            
            # set units
            if unit == 'SI':
                t_unit = ttk.Label(frame,text="K")
                t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
                t_entry.insert(0,305.11)

                ft_unit = ttk.Label(frame,text="m")
                ft_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
                ft_entry.insert(0,500)

                fa_unit = ttk.Label(frame,text="ha")
                fa_unit.grid(row=6, column=2, sticky=tk.W, pady=2)
                fa_entry.insert(0,180)

                dt_unit = ttk.Label(frame,text="m")
                dt_unit.grid(row=7, column=2, sticky=tk.W, pady=2)
                dt_entry.insert(0,1000)
                
#                 dc_unit = ttk.Label(frame,text="m")
#                 dc_unit.grid(row=8, column=2, sticky=tk.W, pady=2)
#                 dc_entry.insert(0,1150)
                
                dz_unit = ttk.Label(frame,text="m")
                dz_unit.grid(row=4, column=4, sticky=tk.W, pady=2)
                dz_entry.insert(0,1200)
                
                gm_unit = ttk.Label(frame,text="MPa/m")
                gm_unit.grid(row=5, column=4, sticky=tk.W, pady=2)
                gm_entry.insert(0,0.00835)
                
                cu_entry.insert(0,0.2)


            elif unit == 'Imperial':
                t_unit = ttk.Label(frame,text="F")
                t_unit.grid(row=4, column=2, sticky=tk.W, pady=2)
                t_entry.insert(0,89.53)
                

                ft_unit = ttk.Label(frame,text="ft")
                ft_unit.grid(row=5, column=2, sticky=tk.W, pady=2)
                ft_entry.insert(0,1640)

                fa_unit = ttk.Label(frame,text="acre")
                fa_unit.grid(row=6, column=2, sticky=tk.W, pady=2)
                fa_entry.insert(0,445)

                dt_unit = ttk.Label(frame,text="ft")
                dt_unit.grid(row=7, column=2, sticky=tk.W, pady=2)
                dt_entry.insert(0,3281)
                
                dz_unit = ttk.Label(frame,text="ft")
                dz_unit.grid(row=4, column=4, sticky=tk.W, pady=2)
                dz_entry.insert(0,3937)
                
                gm_unit = ttk.Label(frame,text="psi/ft")
                gm_unit.grid(row=5, column=4, sticky=tk.W, pady=2)
                gm_entry.insert(0,0.00835)
                
                cu_entry.insert(0,0.2)
            

        show_single()
        
        # salt cavern radiobar
        radioVar = tk.StringVar(value=0)
        R1 = ttk.Radiobutton(radioframe, text="Single cavern with specified information", variable=radioVar, value=0,command=show_single)
        R1.grid(row=1, column=1,sticky=tk.W, columnspan=2)      
        R2 = ttk.Radiobutton(radioframe, text="Multiple equally-sized caverns in a formation", variable=radioVar, value=1,command=show_multiple)
        R2.grid(row=2, column=1,sticky=tk.W, columnspan=2)

        
        
        
        
        
    return results

#####################################################
# main

R = 8.314e-5  # universal gas constant, m3-bar/K-mol
unit = 'SI'


root = tk.Tk()

# Adjust size
root.geometry( "700x500" )
root.resizable(0, 0)
root.title('Hydrogen Storage Capacity Calculation')

# set the frames
topframe = tk.Frame(root)
topframe.grid(row=0,column=0, sticky=tk.W, padx=15, pady=10)

radioframe = tk.Frame(root)
radioframe.grid(row=1,column=0, sticky=tk.W, padx=15)

frame = tk.Frame(root)
frame.grid(row=2,column=0, sticky=tk.W, padx=15)


botframe = tk.Frame(root)
botframe.grid(row=3, column=0, sticky=tk.W, padx=15, pady=30)
#botframe.config(bg="white")  



options = [
"Select a formation",
"Natural Gas Storage Site",
"Depleted Gas Field",
"Saline Aquifer",
"Salt Cavern"]

# run button  
run_text = tk.StringVar()
run_btn = tk.Button(topframe, textvariable=run_text, command=lambda: calc_run(),
                    font=('Raleway', 10, 'bold'), bg='#20bebe',
                    fg='white', height=1, width=7)
run_text.set("Run")
run_btn.grid(row=0, column=2, sticky=tk.W, padx=0, columnspan=2)



clicked = tk.StringVar()
clicked.set(options[0])  # default value

ttk.Label(topframe, text="Formation Type: ").grid(row=0, column=0)
someStyle=ttk.Style()
someStyle.configure('my.TMenubutton',font=('Futura',15))
drop = ttk.OptionMenu(topframe, clicked, *options, command=show, style='my.TMenubutton')
drop.grid(row=0, column=1)
drop.config(width=25) 



# set the unit
def switch_to_Imperial():
    global unit
    tkinter.messagebox.showinfo(title='Unit Change', message='Start Over!')
    clearFrame()
    clearradio()
    clicked.set(options[0])  # default value
    unit = 'Imperial'
    show
def switch_to_SI():
    global unit
    tkinter.messagebox.showinfo(title='Unit Change', message='Start Over!')
    clearFrame()
    clearradio()
    clicked.set(options[0])  # default value
    unit = 'SI'
    show

# unit radiobar
unitVar = tk.StringVar(value=0)
R1 = ttk.Radiobutton(topframe, text="SI units", variable=unitVar, value=0,command=switch_to_SI)
R1.grid(row=1, column=2,sticky=tk.W, columnspan=2)      
R2 = ttk.Radiobutton(topframe, text="Imperial units", variable=unitVar, value=1,command=switch_to_Imperial)
R2.grid(row=2, column=2,sticky=tk.W, columnspan=2)


root.mainloop()    
    
