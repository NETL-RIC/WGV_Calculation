import tkinter as tk
from PIL import Image, ImageTk
import numpy as np
from scipy.optimize import newton

R = 8.314e-5  # universal gas constant, m3-bar/K-mol

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
    a = 0.457235 * R**2 * Tc**2 / Pc
    b = 0.0777961 * R * Tc / Pc
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2

    A = a * alpha * P / R**2 / T**2
    B = b * P / R / T

    # build cubic polynomial
    def g(z):
        """
        Cubic polynomial in z from EOS. This should be zero.
        :param z: float compressibility factor
        """
        return z**3 - (1 - B) * z**2 + (A - 2*B - 3*B**2) * z - (
                A * B - B**2 - B**3)

    # Solve cubic polynomial for the compressibility factor
    z = newton(g, 1.0)  # compressibility factor
    rho = P / (R * T * z)  # density (mol/m3)

    return rho


def wge_calc(T_res, P_res, wgv_in, H2_frac_stp):

    """
    Volumetric calculation of the working gas energy for a natural gas storage site
    :param T_res: float temperature in the reservoir
    :param P_res: float pressure in the reservoir at storage conditions
    :param wgv: float working gas in MMCF
    :param hyd_frac: float fraction of hydrogen in gas mixture at STP
    Returns the total working gas energy, hydrogen working gas energy, methane working gas energy, total working gas volume, hydrogen working gas volume, and methane working gas volume of the storage site
    """
    #Conversion factors
    m3_to_mmcf = 0.00003531 
    mcf_to_m3 = 28.32

    #Working gas volume conversion
    wgv = wgv_in * mcf_to_m3

    #Hydrogen properties
    Tc_H2 = 33.18 #K
    Pc_H2 = 13 #bar
    omega_H2 = -0.220
    mw_H2 = 0.002016 #kg/mol
    lhv_H2 = 3.332e-8 #TWh/kg

    #Methane properties
    Tc_CH4 = 190.6 #K
    Pc_CH4 =  46.1 #bar
    omega_CH4 = 0.011
    mw_CH4 = 0.01604 #kg/mol
    lhv_CH4 = 1.389e-8 #TWh kg

    #STP conditions
    P_stp = 1 #bar
    T_stp = 273.15 #K

    #Fraction of CH4 in gas
    CH4_frac_stp = 1-H2_frac_stp

    #Density of H2 and CH4 at STP in kg/m^3
    rho_H2_stp = preos(Tc_H2, T_stp, Pc_H2, P_stp, omega_H2)*mw_H2
    rho_CH4_stp = preos(Tc_CH4, T_stp, Pc_CH4, P_stp, omega_CH4)*mw_CH4

    #Density of H2 and CH4 at reservoir conditions
    rho_H2_res = preos(Tc_H2, T_res, Pc_H2, P_res, omega_H2)*mw_H2
    rho_CH4_res = preos(Tc_CH4, T_res, Pc_CH4, P_res, omega_CH4)*mw_CH4

    #Hydrogen fraction in gas at reservoir conditions
    if H2_frac_stp == 1:
        H2_frac_res = H2_frac_stp
        CH4_frac_res = CH4_frac_stp
    else:
        H2_frac_res = ((rho_H2_stp / rho_H2_res) * H2_frac_stp) / (((rho_H2_stp / rho_H2_res) * H2_frac_stp) + ((rho_CH4_stp / rho_CH4_res) * CH4_frac_stp))
        CH4_frac_res = 1-H2_frac_res

    # print(rho_CH4_stp / rho_CH4_res)
    # print(wgv)
    # print(rho_H2_res)
        
    wge_H2 = lhv_H2 * rho_H2_res * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv
    wge_CH4 = lhv_CH4 * rho_CH4_res * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv  
    wge = wge_H2 + wge_CH4

    wgv_H2 =  ((rho_H2_res / rho_H2_stp) * H2_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv) * m3_to_mmcf
    wgv_CH4 = ((rho_CH4_res / rho_CH4_stp) * CH4_frac_res * (rho_CH4_stp / rho_CH4_res) * wgv) * m3_to_mmcf
    wgv_out = wgv_H2 + wgv_CH4

    return {'Working Gas Energy (TWh)' : wge, 'Working Gas Energy of Hydrogen (TWh)' : wge_H2, 'Working Gas Energy of Methane (TWh)' : wge_CH4,           
            'Working Gas Volume (MMCF)': wgv_out, 'Working Gas Volume of Hydrogen (MMCF)': wgv_H2, 'Working Gas Volume of Methane (MMCF)': wgv_CH4}


# --------------- Graphical User Interface ----------------
# Running and generating results in the blanks
def calc_run():
    run_text.set("Running...")
    

    results_dic = wge_calc(float(t_entry.get()), float(p_entry.get()), float(wgv_entry.get()), float(h2_entry.get()))
    
    wge_output.delete(0, tk.END)
    wge_H2_output.delete(0, tk.END)
    wge_CH4_output.delete(0, tk.END)
    wgv_out_output.delete(0, tk.END)
    wgv_H2_output.delete(0, tk.END)
    wgv_CH4_output.delete(0, tk.END)
    
    wge_output.insert(0,round(results_dic.get('Working Gas Energy (TWh)'),4))
    wge_H2_output.insert(0,round(results_dic.get('Working Gas Energy of Hydrogen (TWh)'),4))
    wge_CH4_output.insert(0,round(results_dic.get('Working Gas Energy of Methane (TWh)'),4))
    wgv_out_output.insert(0,round(results_dic.get('Working Gas Volume (MMCF)'),4))
    wgv_H2_output.insert(0,round(results_dic.get('Working Gas Volume of Hydrogen (MMCF)'),4))
    wgv_CH4_output.insert(0,round(results_dic.get('Working Gas Volume of Methane (MMCF)'),4))


    run_text.set("Run")

# create a GUI for input and output
    
root = tk.Tk()

root.title("Calculation of Working Gas Volume")

root.geometry('600x500')

# logo
basewidth = 250
logo = Image.open('logo.png')
wpercent = (basewidth/float(logo.size[0]))
hsize = int((float(logo.size[1])*float(wpercent)))
logo = logo.resize((basewidth,hsize), Image.ANTIALIAS)
logo = ImageTk.PhotoImage(logo)
logo_label = tk.Label(image=logo)
logo_label.image = logo
logo_label.pack()

                      
# Input label
title_msg = "Volumetric Calculation of \nThe Working Gas Energy for A Natural Gas Storage Site"
label = tk.Label(root, text=title_msg, font = ('Helvetica', 12, 'bold'))
label.place(relx=0.5,rely=0.20, anchor='center')

label_input = tk.Label(root, text='Enter Values:', font=('Helvetica 10 underline bold'))
label_input.place(relx=0.05,rely=0.30)

# inputs
t_label = tk.Label(root, text='Reservoir Temperature (K)')
t_label.place(relx=0.05,rely=0.35)
t_entry = tk.Entry(root, width=10)
t_entry.insert(0,'305.11')
t_entry.place(relx=0.35,rely=0.35)

p_label = tk.Label(root, text='Reservoir Pressure (bar)')
p_label.place(relx=0.05,rely=0.40)
p_entry = tk.Entry(root, width=10)
p_entry.insert(0,'76.92')
p_entry.place(relx=0.35,rely=0.40)

wgv_label = tk.Label(root, text='Workng Gas Volume (MCF)')
wgv_label.place(relx=0.05,rely=0.45)
wgv_entry = tk.Entry(root, width=10)
wgv_entry.insert(0,2200000)
wgv_entry.place(relx=0.35,rely=0.45)

h2_label = tk.Label(root, text='Hydrogen Fraction')
h2_label.place(relx=0.05,rely=0.50)
h2_entry = tk.Entry(root, width=10)
h2_entry.insert(0,0.1)
h2_entry.place(relx=0.35,rely=0.50)


# run button
run_text = tk.StringVar()
run_btn = tk.Button(root, textvariable = run_text, command=lambda:calc_run(), font=('Raleway',10, 'bold'), bg='#20bebe', fg='white', height=1, width=8)
run_text.set("Run")
run_btn.place(relx=0.05,rely=0.55)

# output label
label_output = tk.Label(root, text='Results:', font=('Helvetica 10 underline bold'))
label_output.place(relx=0.35,rely=0.60)


# output values
wge_label = tk.Label(root, text='Working Gas Energy (TWh)')
wge_label.place(relx=0.35,rely=0.65)
wge_output = tk.Entry(root, width=10)
wge_output.place(relx=0.80,rely=0.65)

wge_H2_label = tk.Label(root, text='Working Gas Energy of Hydrogen (TWh)')
wge_H2_label.place(relx=0.35,rely=0.70)
wge_H2_output = tk.Entry(root, width=10)
wge_H2_output.place(relx=0.80,rely=0.70)

wge_CH4_label = tk.Label(root, text='Working Gas Energy of Methane (TWh)')
wge_CH4_label.place(relx=0.35,rely=0.75)
wge_CH4_output = tk.Entry(root, width=10)
wge_CH4_output.place(relx=0.80,rely=0.75)

wgv_out_label = tk.Label(root, text='Working Gas Volume (MMCF)')
wgv_out_label.place(relx=0.35,rely=0.80)
wgv_out_output = tk.Entry(root, width=10)
wgv_out_output.place(relx=0.80,rely=0.80)

wgv_H2_label = tk.Label(root, text='Working Gas Volume of Hydrogen (MMCF)')
wgv_H2_label.place(relx=0.35,rely=0.85)
wgv_H2_output = tk.Entry(root, width=10)
wgv_H2_output.place(relx=0.80,rely=0.85)

wgv_CH4_label = tk.Label(root, text='Working Gas Volume of Methane (MMCF)')
wgv_CH4_label.place(relx=0.35,rely=0.90)
wgv_CH4_output = tk.Entry(root, width=10)
wgv_CH4_output.place(relx=0.80,rely=0.90)


root.resizable(False, False)

root.mainloop()

