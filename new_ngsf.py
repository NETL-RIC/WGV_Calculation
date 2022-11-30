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
