
# Imports 
import numpy as np
from scipy.optimize import fsolve
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import plotly.graph_objects as go

# units 
mm = 0.001
m = 1
N = .001
kN = 1
kPa = 1
MPa= 1000
# Inputs 


# Wall Properties

faim = 0.6
fais = 0.85
emu = 0.003
k= 1







def cross_section(t, s, db,fblock):
    
    if db == 10 * mm:
        As = 100 * mm**2
    elif db == 15 * mm:
        As = 200 * mm**2
    elif db == 20 * mm:
        As = 300 * mm**2
    elif db == 25 * mm:
        As = 500 * mm**2


    
    # Material properties

    if fblock == 10:
        fm_g = 6
        fm_ug = 7.5
    elif fblock == 15:
        fm_g = 9
        fm_ug = 12
    elif fblock == 20:
        fm_g = 12
        fm_ug = 15.5
    elif fblock == 25:
        fm_g = (12+16)/2
        fm_ug = (15.5+20.5)/2
    elif fblock == 30:
        fm_g = 16
        fm_ug = 20.5

    if t == 140 *mm: 
        tf = 26 *mm 
        rho_g = 3.03 *kN/m**2
        rho_ug = 1.67 *kN/m**2
    elif t == 190 *mm:  
        tf = 36.2 *mm
        rho_g = 4.12 *kN/m**2
        rho_ug = 2.19 *kN/m**2
    elif t == 240 *mm:  
        tf = 38.6 *mm
        rho_g = 5.22 *kN/m**2
        rho_ug = 2.62 *kN/m**2
    elif t == 290 *mm:  
        tf = 41.4 *mm
        rho_g = 6.32 *kN/m**2
        rho_ug = 3.059 *kN/m**2
    
    beff = np.minimum(6*t, s)
    beff_m_1= 1000 *mm
    beff_m_2= beff * 1/s
    Aseff_m = As *1/s
    bg_m = 200*mm * 1/s
    bug_m_1 = 1000*mm - bg_m  # Bars in compression 
    bug_m_2 = beff_m_2- bg_m  #Bars in tension

    # Area effective 
    A_gr= bg_m * t
    A_ug_1 = bug_m_1 * 2*tf    # Bars in compression 
    A_ug_2 = bug_m_2 * 2*tf   # Bars in tension

    Ae_1 = A_gr + A_ug_1 # Bars in compression 
    Ae_2 = A_gr + A_ug_2 # Bars in tension

    # fm effective

    fm_e_1 = (fm_g * A_gr + fm_ug * A_ug_1) / (Ae_1)
    fm_e_2 = (fm_g * A_gr + fm_ug * A_ug_2) / (Ae_2)
    E_m = 850 * fm_e_2
    d= t/2 
    # Inertia 
    I_gross_gr=  beff_m_1 * t**3 / 12
    I_gross_ug_1 = (beff_m_1 * tf**3 / 12 + beff_m_1 * tf * (t/2 - tf/2)**2) *2
    I_gross_eff = (I_gross_gr * A_gr + + I_gross_ug_1 * A_ug_1)/ (Ae_1)  
    S_eff = I_gross_eff / (t/2)
    n = 200000/E_m
    kd = sp.Symbol('kd', real=True)
    eq = n * Aseff_m * (d - kd) - (beff_m_2 * kd**2) / 2
    solutions  = sp.solve(eq, kd)
    kd= [sol.evalf() for sol in solutions if sol.is_real and sol > 0][0] 
    I_cr_eff = n * Aseff_m * (d- kd)**2 + (beff_m_2 * kd**3) / 3    
    ek = S_eff / Ae_1

   # Self Weight 
    rho_SW = (rho_g * A_gr + rho_ug * A_ug_1) / (Ae_1) 

    return  t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf


def solve_betaC(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t):
    # Compute Prmax
    Prmax = 0.8 * faim * 0.85 * fm_e_1 * Ae_1 * 1000  # kN
    
    # Define the equation to solve for betaC
    def equation(betaC):
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
        F_g_bottom = faim * 0.85 * fm_ug * (betaC - (t - tf)) * bug_m_1 * 1000 if betaC >= (t - tf) else 0
        
        if betaC >= (t - tf):
            Pr = Fg + F_ug_top + F_g_bottom
        else:
            Pr = Fg + F_ug_top
        return Pr - Prmax
    
    # Solve for betaC
    betaC_solution = fsolve(equation, t)  # Initial guess is t
    betaC = betaC_solution[0]
    
    # Recalculate forces with the obtained betaC

    Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
    F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
    F_g_bottom = faim * 0.85 * fm_ug * (betaC - (t - tf)) * bug_m_1 * 1000 if betaC >= (t - tf) else 0
    
    if betaC >= (t - tf):
        Mr = Fg *(t/2 - betaC/2) + F_ug_top * (t/2-tf/2) - F_g_bottom * (t/2 - (tf-(t-betaC))/2) # moment at center 
        Mr2 = Fg *betaC/2 + F_ug_top * tf/2 + F_g_bottom * (t+(tf-(t-betaC)/2)-tf) - Prmax * t/2 # moment at top 
    else : 
        Mr = Fg * (t/2 - betaC/2) + F_ug_top * (t/2 - tf/2) # moment at center
        Mr2 = Fg *betaC/2 + F_ug_top * tf/2  - Prmax * t/2 # moment at top 
    return betaC, Fg, F_ug_top, F_g_bottom, Prmax, Mr


def calculate_point2(betaC1, faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t, d, num_points=3):
    betaC_values = np.linspace(0.8*(t - d), betaC1, num=num_points)[::-1]
    results = []
    
    for betaC in betaC_values:
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_1 * 1000
        Pr = Fg + F_ug_top
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2)
        e = Mr / Pr 
        if e > t/3:
            Mr = Pr * t/3
        results.append((betaC, Fg, F_ug_top, Pr, Mr))
    
    return results

def calculate_point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t, d, Aseff_m, num_points):
    betaC_values = np.linspace(tf / 2, 0.8 * (t - d), num=num_points)[::-1]
    results = []
    
    Mr_y, Pr_y = None, None  # Initialize to None
    i = 0  # Initialize index
    
    while i < len(betaC_values):
        betaC = betaC_values[i]
        c = betaC / 0.8
        es = emu * (d - c) / c     
        fs = np.minimum(200000 * es, 400)
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        
        if betaC < tf:
            F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
        else: 
            F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
        
        Ts = fais * Aseff_m * fs * 1000
        Pr = Fg + F_ug_top - Ts 
        
        if betaC < tf:
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
        else: 
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
        
        e = Mr / Pr 

        # Store the first values where es >= 0.002
        if es >= 0.002 and Mr_y is None and Pr_y is None:
            ey, Mr_y, Pr_y = es, Mr, Pr
        
        if Pr <= 0:
            ey=0
            break  # Exit loop if Pr is non-positive

        results.append((betaC, Fg, F_ug_top, Pr, Mr))
        i += 1  # Increment index
    
    return results, Mr_y, Pr_y, ey

# Point 4 (pure tension moment)
def calculate_pure_moment(faim, Aseff_m, d, fm_g, bg_m, fm_ug, tf, bug_m_2, t):
    # Define the equation to solve for betaC
    def equation(betaC):
        Ptarget = 0.001 
        c = betaC / 0.8
        es = emu * (d - c) / c     
        fs = np.minimum(200000 * es, 400)
        Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
        
        if betaC < tf:
            F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
        else: 
            F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
        
        Ts = fais * Aseff_m * fs * 1000
        Pr = Fg + F_ug_top - Ts 
        if betaC < tf:
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
        else:
            Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
        return Pr - Ptarget
    
    # Solve for betaC
    betaC_solution = fsolve(equation, tf/2)  # Initial guess is t
    betaC = betaC_solution[0]
    
    c = betaC / 0.8
    es = emu * (d - c) / c     
    fs = np.minimum(200000 * es, 400)
    Fg = faim * 0.85 * fm_g * betaC * bg_m * 1000
    
    if betaC < tf:
        F_ug_top = faim * 0.85 * fm_ug * betaC * bug_m_2 * 1000
    else: 
        F_ug_top = faim * 0.85 * fm_ug * tf * bug_m_2 * 1000
    
    Ts = fais * Aseff_m * fs * 1000
    Pr = Fg + F_ug_top - Ts 
    if betaC < tf:
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - betaC / 2) - Ts * (t / 2 - d)
    else: 
        Mr = Fg * (t / 2 - betaC / 2) + F_ug_top * (t / 2 - tf / 2) - Ts * (t / 2 - d)
    return betaC, Fg, F_ug_top, Pr, Mr
# def PureMoment(faim, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t):
    


def Icr_function(Aseff_m, d, beff_m_2, beff_g, E_m,  PF, tf):
    fy = 400 *MPa 
    E_m = E_m * MPa
    n = 200000*MPa / E_m
    kd = sp.Symbol('kd', real=True, positive=True)
    fm = 0.002 * kd *E_m / (d-kd)

    # Try both equations
    eq1 = 0.5* kd * fm *beff_m_2 - ( PF + fy*Aseff_m)

    sol1 = [sol.evalf() for sol in sp.solve(eq1, kd) if sol.is_real and sol > 0]

    if not sol1:
        raise ValueError("No valid solution found for kd in eq1.")

    kd_val = sol1[0]

    if kd_val <= tf:
        I_cr_eff = n * Aseff_m * (d - kd_val)**2 + (beff_m_2 * kd_val**3) / 3    
    else:

        fm = 0.002 * kd *E_m / (d-kd)
        fm_g = fm * (kd - tf) / kd
        eq2 = 0.5 * (fm + fm_g) * tf * beff_m_2 + 0.5 * fm_g* (kd - tf) * beff_g - ( PF + fy*Aseff_m)
        sol2 = [sol.evalf() for sol in sp.solve(eq2, kd) if sol.is_real and sol > tf]


        if not sol2:
            raise ValueError("No valid solution found for kd in eq2.")
        kd_val = sol2[0]
        I_cr_eff = n * Aseff_m * (d- kd_val)**2 +beff_m_2 * tf**3/12 + (beff_m_2 * tf * (kd_val-tf/2))  + beff_g * (kd_val-tf)**3 / 3 



    return I_cr_eff, kd_val


def Moment_Calculation (t,e ,H, rho_SW, W, P_DL,P_LL, P_S,I_gross_eff, E_m, ek, I_cr_eff ):
    k=1
    E_m = E_m * MPa
    # Loads Calculation 
    e_P = np.maximum(0.1 * t, e*mm)
    e = e *mm 

    # Self Weight
    P_SW_mid = rho_SW  * H/2  # at mid span 
    P_SW_base = rho_SW  * H # at base
    lambda_h = H / t
    # moment due to wind 
    M_lateral = W * H**2 / 8
    M_unf = max((P_DL + P_LL + P_S) * e + M_lateral, (P_DL + P_LL+ P_S) * 0.1 * t)

    M_lateral_F1 = 0
    M_lateral_F2 = 0
    M_lateral_F3 = 0
    M_lateral_F4 = 1.4 *  M_lateral
    M_lateral_F5 = 1.4 *  M_lateral
    M_lateral_F6 = 0  # Snow + DL 
    M_lateral_F7 = 0  # Snow + DL 

    M_lateral_F = [M_lateral_F1, M_lateral_F2, M_lateral_F3, M_lateral_F4, M_lateral_F5, M_lateral_F6, M_lateral_F7]    

    # Factored Primary Moment at midspan 
    M_F1 = np.maximum(1.4 * P_DL * e /2, 1.4 * P_DL * 0.1 *t)
    M_F2 = np.maximum(1.25 * P_DL * e /2+ 1.5 * P_LL * e /2, 1.25 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t)
    M_F3 = np.maximum(0.9 * P_DL * e/2 + 1.5 * P_LL * e  /2, 0.9 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t)
    M_F4 = np.maximum(1.25 * P_DL * e /2+ 1.4 * M_lateral, 1.25 * P_DL * 0.1 *t)
    M_F5 = np.maximum(0.9 * P_DL * e /2+ 1.4 * M_lateral, 0.9 * P_DL * 0.1 *t)
    M_F6 = np.maximum(1.25 * P_DL * e /2 + 1.5 * P_S * e/2, 1.25 * P_DL * 0.1 *t + 1.5 * P_S * 0.1 *t)
    M_F7 = np.maximum(0.9 * P_DL * e /2 + 1.5 * P_S * e/2, 0.9 * P_DL * 0.1 *t + 1.5 * P_S * 0.1 *t)

    M_F = [M_F1, M_F2, M_F3, M_F4, M_F5, M_F6, M_F7]

    # Factored Primary Axial Load @ midspan
    P_F1 = 1.4 * (P_DL  + P_SW_mid)
    P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F4 = 1.25 * (P_DL  + P_SW_mid) 
    P_F5 = 0.9 * (P_DL  + P_SW_mid) 
    P_F6 = 1.25 * (P_DL  + P_SW_mid) + 1.5 *P_S # Snow load + DL
    P_F7 = 0.9 * (P_DL  + P_SW_mid) + 1.5* P_S   # Snow load + DL 

    P_F = [P_F1, P_F2, P_F3, P_F4, P_F5, P_F6, P_F7]
    # top Moment 
    M1_F1 =  np.maximum(1.4 * P_DL * e, 1.4 * P_DL *0.1*t)
    M1_F2 =  np.maximum(1.25 * P_DL * e+ 1.5 * P_LL * e, 1.25 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t) 
    M1_F3 =  np.maximum(0.9 * P_DL * e + 1.5 * P_LL * e, 0.9 * P_DL * 0.1 *t + 1.5 * P_LL * 0.1 *t)  
    M1_F4 =  np.maximum(1.25 * P_DL * e , 1.25 * P_DL * 0.1 *t)
    M1_F5 =  np.maximum(0.9 * P_DL * e  , 0.9 * P_DL * 0.1 *t)
    M1_F6 =  np.maximum(1.25 * P_DL * e+ 1.5 * P_S * e, 1.25 * P_DL * 0.1 *t + 1.5 * P_S * 0.1 *t) 
    M1_F7 =  np.maximum(0.9 * P_DL * e + 1.5 * P_S * e, 0.9 * P_DL * 0.1 *t + 1.5 * P_S * 0.1 *t)  

   

    Mtop_F= [M1_F1, M1_F2, M1_F3, M1_F4, M1_F5, M1_F6, M1_F7]

    # Bottom Moment 
    M2_F1 = np.maximum(0,  1.4 * (P_DL  + P_SW_base) *0.1*t)
    M2_F2 = np.maximum(0,  1.25 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_LL * 0.1 *t) 
    M2_F3 = np.maximum(0,  0.9 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_LL * 0.1 *t)  
    M2_F4 = np.maximum(0, 1.25 * (P_DL  + P_SW_base) * 0.1 *t)
    M2_F5 = np.maximum(0, 0.9 * (P_DL  + P_SW_base) * 0.1 *t)
    M2_F6 = np.maximum(0,  1.25 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_S * 0.1 *t) 
    M2_F7 = np.maximum(0,  0.9 * (P_DL  + P_SW_base) * 0.1 *t + 1.5 * P_S * 0.1 *t) 

    Mbot_F= [M2_F1, M2_F2, M2_F3, M2_F4, M2_F5, M2_F6, M2_F7]

    # Moment Ratio
    M_ratio =  np.minimum(np.array(Mtop_F), np.array(Mbot_F)) / np.maximum(np.array(Mtop_F), np.array(Mbot_F))
    I_cr_array = np.array(I_cr_eff)  # assuming I_cr_eff is a list
    ev = np.array([M_F1, M_F2, M_F3, M_F4, M_F5, M_F6,M_F7]) / np.array([P_F1, P_F2, P_F3, P_F4, P_F5, P_F6 ,P_F7])
    # Rigidty Coefficient 
    betad = P_DL * e_P / M_unf
    faiE = 0.75
    EI_eff_raw = E_m * (0.25 * I_gross_eff - (0.25 * I_gross_eff - I_cr_array)* ((ev - ek) / (2 * ek)))
    Min_EIeff= E_m * I_cr_array
    Max_EIeff=  0.25 * E_m * I_gross_eff
    EI_eff = np.clip(EI_eff_raw, Min_EIeff, Max_EIeff)

    Rigidty_c = faiE * EI_eff / (1+0.5*betad)

    Pcr = np.pi**2 * Rigidty_c/ (k*H)**2
    # Lateral Force Coefficient
    Cm = []
    for i in range(len(M_F)):
        if lambda_h < 30:
            Cm.append(np.where(M_lateral_F[i] / M_F[i] > 0.5, 1, np.maximum(0.6 + 0.4 * M_ratio[i], 0.4)))
        else:
            Cm.append(1)
    Cm = np.array(Cm)


    MagFactor=  Cm / (1-(np.array(P_F)/np.array(Pcr)))

    # Design Moment

    Mt_F = np.array(M_F) * np.array(MagFactor)

    # 
    Mt_F_list = [float(val) for val in Mt_F]
    M_F_list = [float(val) for val in M_F]
    P_F_list = [float(val) for val in P_F]
    Mag_list = [float(val) for val in MagFactor]
    Cm_list = [float(val) for val in Cm]
    EI_eff_list = [float(val) for val in EI_eff]
    Pcr_list = [float(val) for val in Pcr]
    Max_EIeff_list = [Max_EIeff] * len(I_cr_array)
    betad_list = [betad] * len(I_cr_array)

    return ev, Mt_F_list, M_F_list, P_F_list,Pcr_list, Mag_list, Cm_list,EI_eff_raw, EI_eff_list, betad, Min_EIeff , Max_EIeff_list, I_cr_array, betad_list


def draw_blocks_plotly(t_mm_actual, s, bar_diameter_mm=10, num_blocks=5, gap=10, width=350, fixed_height=160, fixed_faceshell=25):
    fig = go.Figure()

    # Determine the grouting pattern based on spacing s (assuming s is in mm)
    grout_pattern = s // 200  # Number of cells in a grouting cycle

    for i in range(num_blocks):
        x_offset = i * (width + gap)

        # Outer rectangle (block boundary)
        fig.add_shape(type="rect",
                      x0=x_offset, y0=0,
                      x1=x_offset + width, y1=fixed_height,
                      line=dict(width=2, color='black'),
                      fillcolor='lightgray',
                      layer="below") # Ensure block is below traces

        # Openings (inner rectangles)
        opening1_color = 'white'
        opening2_color = 'white'
        grout_center_x = []
        grout_center_y = []

        # Apply grouting pattern
        if (i * 2) % grout_pattern == 0:
            opening1_color = 'darkgrey'  # Grout the first opening
            gx1 = x_offset + width * 0.1 + width * 0.35 / 2
            gy1 = fixed_faceshell + (fixed_height - 2 * fixed_faceshell) / 2
            grout_center_x.append(gx1)
            grout_center_y.append(gy1)
        if (i * 2 + 1) % grout_pattern == 0:
            opening2_color = 'darkgrey'  # Grout the second opening
            gx2 = x_offset + width * 0.55 + width * 0.35 / 2
            gy2 = fixed_faceshell + (fixed_height - 2 * fixed_faceshell) / 2
            grout_center_x.append(gx2)
            grout_center_y.append(gy2)

        # Draw openings
        fig.add_shape(type="rect",
                      x0=x_offset + width * 0.1, y0=fixed_faceshell,
                      x1=x_offset + width * 0.1 + width * 0.35, y1=fixed_faceshell + (fixed_height - 2 * fixed_faceshell),
                      line=dict(width=2, color='black'),
                      fillcolor=opening1_color,
                      layer="below") # Ensure opening is below traces
        fig.add_shape(type="rect",
                      x0=x_offset + width * 0.55, y0=fixed_faceshell,
                      x1=x_offset + width * 0.55 + width * 0.35, y1=fixed_faceshell + (fixed_height - 2 * fixed_faceshell),
                      line=dict(width=2, color='black'),
                      fillcolor=opening2_color,
                      layer="below") # Ensure opening is below traces

        # Add larger dots to represent bars in the grouted cells - ADDED AFTER SHAPES
        if grout_center_x: # Only add if there was grouting in this block
            bar_size = bar_diameter_mm * 0.45 # Using a moderate size
            fig.add_trace(go.Scatter(x=grout_center_x, y=grout_center_y,
                                     mode='markers', marker=dict(color='red', size=bar_size), # Red for bars, larger size
                                     showlegend=False,
                                     hoverinfo='none'))


    # Add text under the figure to mention the actual block thickness
    fig.add_annotation(x=num_blocks * (width + gap) / 2,  # Center horizontally
                         y=-30,  # Position below the visual thickness text
                         text=f'Selected Block Thickness: {t_mm_actual} mm',
                         showarrow=False,
                         font=dict(size=12),
                         align="center")

    # Formatting
    fig.update_xaxes(range=[-1, num_blocks * (width + gap)], visible=False, fixedrange=True)
    fig.update_yaxes(range=[-35, fixed_height + 1], visible=False, scaleanchor="x", scaleratio=1, fixedrange=True) # Adjusted y-axis range to accommodate the new text
    fig.update_layout(width=3 * 4 * 50, height=2 * 50 * 1.2,  # Further adjusted height
                      margin=dict(l=0, r=0, b=0, t=0), # Increased bottom margin
                      plot_bgcolor='white',
                      showlegend=False)    

    return fig


def generate_side_view(H, P_DL, P_LL, P_S, e, W):
    fig = go.Figure()

    # --- Wall Rectangle ---
    fig.add_shape(
        type="rect",
        x0=0.4, x1=1,
        y0=0, y1=8,
        line=dict(color="black"),
        fillcolor="lightgrey"
    )


    # --- Dimension lines ---
    fig.add_shape(
        type="line",
        x0=1.5, x1=1.5,
        y0=0, y1=8,
        line=dict(color="grey")
    )

    fig.add_shape(
        type="line",
        x0=1.1, x1=1.6,
        y0=8, y1=8,
        line=dict(color="grey")
    )

    fig.add_shape(
        type="line",
        x0=1.2, x1=1.7,
        y0=0, y1=0,
        line=dict(color="grey")
    )

    # --- Center dashed line ---
    fig.add_shape(
        type="line",
        x0=0.7, x1=0.7,
        y0=0, y1=8+.5,
        line=dict(color="black", dash="dash")
    )

    # --- Dead Load (PDL) arrow ---
    fig.add_annotation(
        x=0.7+2*e/1000, y=8,
        text=f"<i>P</i><sub>DL</sub> = {P_DL:.2f} kN",
        showarrow=True,
        arrowhead=3, arrowsize=1, arrowwidth=2,
        ax=0.7-2+2*e/1000, ay=8-50,    # Base is above tip → points downward
        font=dict(color="blue", size=16, family="Times New Roman")
    )


    # --- Live Load (PLL) arrow ---
    fig.add_annotation(
        x=0.7+2*e/1000, y=8+1.3,
        text=f"<i>P</i><sub>LL</sub> = {P_LL:.2f} kN",
        font=dict(color="red",size=16, family="Times New Roman"),showarrow=False,
    )

        # --- Snow Load (P_S) arrow ---
    fig.add_annotation(
        x=0.7+2*e/1000, y=8+1.6,
        text=f"<i>P</i><sub>S</sub> = {P_S:.2f} kN",
        font=dict(color="black",size=16, family="Times New Roman"),showarrow=False,
    )

    # --- Wind Load arrows (left to right) ---
    num_arrows = 15 # Number of wind arrows

    for i in range(num_arrows):
        y_pos = (i + 1) * 8 / (num_arrows + 1)  # Even spacing along wall height

        # Horizontal shaft
        fig.add_shape(
            type="line",
            x0=0, y0=y_pos,
            x1=0.3, y1=y_pos,
            line=dict(color="blue", width=1)
        )

        # Arrowhead (triangle at tip x=0.4)
        fig.add_shape(
            type="path",
            path=f"M 0.25 {y_pos+0.05} L 0.25 {y_pos-0.05} L 0.38 {y_pos} Z",
            fillcolor="blue",
            line=dict(color="blue")
        )



    
    # --- Wind Load label ---
    fig.add_annotation(
        x=-0.2, y=8 / 2,
        text=f"<i>W</i> = {W:.2f} kpa",
        showarrow=False,
        font=dict(color="blue", size=14),
        xanchor="center", textangle=-90
    )

    # --- Support triangle (black) ---
    fig.add_shape(
        type="path",
        path="M 0.5 -0.4 L 0.9 -0.4 L 0.7 0 Z",
        fillcolor="black",
        line=dict(color="black")
    )
    # --- Roller Support at Top Right ---

    # Define the size of the roller triangle
    roller_size = 0.3  # You can adjust this

    # Add the roller (triangle pointing LEFT at top-right of wall)
    fig.add_shape(
        type="path",
        path=f"M {1.4} {8+0.15} L {1.4} {8+0.15 - roller_size} L {1.4 - roller_size} {8+0.15 - roller_size/2} Z",
        fillcolor="black",
        line=dict(color="black")
    )



    # --- Eccentricity label ---
    eccentricity_arrow_x = 0.7 + 2 * e / 1000  # Same as PDL arrow
    fig.add_annotation(
        x=eccentricity_arrow_x + 0.05,  # Shift slightly to the right
        y=8 + 0.35,
        text=f"e={e:.0f}",
        showarrow=False,
        font=dict(color="purple", size=14),xanchor="left"

    )


    # --- Wall height label ---
    fig.add_annotation(
        x=1.25, y=8/2,
        text=f"<i>h</i> = {H:.2f} m",
        showarrow=False,
        font=dict(size=14), textangle=-90
    )

    # --- Layout settings ---
    fig.update_layout(
        width=400,
        height=500,
         margin=dict(t=0, b=0, l=0, r=0),
        xaxis=dict(visible=False),
        yaxis=dict(range=[-0.5, 8 +3], visible=False),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent plot area 
        paper_bgcolor='rgba(0,0,0,0)'  # Transparent background
    )

    return fig

