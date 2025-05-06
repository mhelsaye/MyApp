import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from PIL import Image
import sympy as sp
import io
import base64
import numpy as np
from scipy.optimize import fsolve
import sympy as sp
from App_Functions import  solve_betaC,  calculate_point2, calculate_point3, calculate_pure_moment, generate_side_view, cross_section, Icr_function, Moment_Calculation, draw_blocks_plotly, generate_side_view

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


# units 
mm = 0.001
m = 1
N = .001
kN = 1
kPa = 1
MPa = 1000

# Inputs 
H = 8 *m 
t = 240 *mm  
s = 1200 *mm 
db = 25 *mm
d = t/2 
fblock= 25 
faim = 0.6
fais = 0.85
emu = 0.003
k= 1



# Default values
H_default = 8000 *mm   # mm
t_options = np.array([140, 190, 240, 290]) *mm  # mm
t_default = 240 *mm   # mm
fblock_options = np.array([10, 15, 20, 25, 30])   # MPa
fblock_default = 25  # MPa
S_options = np.array([200,400,600,800,1000,1200, 1400, 1600]) *mm  # mm
S_default = 1200 *mm  # mm
bar_options = np.array([10, 15, 20, 25]) *mm  # mm
bar_default = 25 *mm  # mm
P_DL_default = 10 *kN # kN
P_LL_default = 0  *kN  # kN
P_S_default = 0  *kN  # 
e_default = 0  *mm # mm
W_default = 1 *kPa # kPa



# Layout

def generate_layout():

    return dbc.Container([

        html.H1("Masonry Walls Design Tool", className="text-center my-4", style={"font-size": "32px"}),
        html.H1("(Out-of-Plane)", className="text-center my-4", style={"font-size": "28px"}),

        dbc.Row([

            dbc.Col([

                html.H5("Inputs", className="text-primary", style={"font-size": "24px"}),

                dbc.Label("Wall Height[m]", style={"font-size": "18px"}),

                dbc.Input(id="input-H", type="number", value=H_default, step=1, min=0, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Wall Thickness [mm]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-t", options=[{"label": f"{int(val*1000)} ", "value": val} for val in t_options], value=t_default, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Block Strength (fblock) [MPa]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-fblock", options=[{"label": f"{val} ", "value": val} for val in fblock_options], value=fblock_default, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Bar Spacing [mm]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-S", options=[{"label": f"{int(val*1000)}", "value": val} for val in S_options], value=S_default, style={"margin-bottom": "30px", "font-size": "16px"}),

                dbc.Label("Rebar Diameter [mm]", style={"font-size": "18px"}),

                dcc.Dropdown(id="dropdown-bar", options=[{"label": f"{int(val*1000)} ", "value": val} for val in bar_options], value=bar_default, style={"margin-bottom": "30px", "font-size": "16px"}),
            ], width=2),

            dbc.Col([

                html.H5("Loads", className="text-primary", style={"font-size": "24px"}),

                dbc.Label("Dead Load  [kN/m]", style={"font-size": "18px"}),

                dbc.Input(id="input-P_DL", type="number", value=P_DL_default, step=0.1, min=0, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Live Load  [kN/m]", style={"font-size": "18px"}),

                dbc.Input(id="input-P_LL", type="number", value=P_LL_default, step=0.1, min=0, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Snow Load  [kN/m]", style={"font-size": "18px"}),

                dbc.Input(id="input-P_S", type="number", value=P_S_default, step=0.1, min=0, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Eccentricity  [mm]", style={"font-size": "18px"}),

                dbc.Input(id="input-e", type="number", value=e_default, step=0.1, style={"margin-bottom": "30px", "font-size": "18px"}),

                dbc.Label("Wind Load [kPa]", style={"font-size": "18px"}),

                dbc.Input(id="input-W", type="number", value=W_default, step=0.1, style={"margin-bottom": "30px", "font-size": "18px"}),

               

            ], width=2),
            

            dbc.Col([
                html.Div([
                    html.Div(id="side-view", style={"maxWidth": "400px", "margin": "0 auto"}),
                    html.Button("Check", id="btn-check", className="btn btn-primary mt-4", style={
                        "font-size": "18px",
                        "display": "block",
                        "margin": "0 auto"
                    }),
                ])
            ], width=2),

            dbc.Col([

                html.H5("Interaction Diagram", className="text-primary", style={"textAlign": "center"}),

                html.Div(id="interaction-diagram",  style={"height": "600px"}),

            ], width=6),


        ], className="my-4"),


        dbc.Row([
            dbc.Col(
                html.Div(id="wall-image", style={
                    "overflowX": "auto",
                    "marginTop": "-72px",  # reduces vertical gap
                    "textAlign": "center"
                }),
                width=6,
                className="offset-0"  
            )
        ], justify="start", className="mt-0"),
        

        dbc.Row([
            dbc.Col(width=1),
            dbc.Col(
                dcc.Markdown(id="done-message", mathjax=True, style={
                                                                        "whiteSpace": "pre-wrap",
                                                                        "wordWrap": "break-word",
                                                                        "overflowX": "auto",
                                                                        "maxWidth": "100%",
                                                                        "fontSize": "20px",
                                                                        # "fontFamily": "monospace",
                                                                        # "textAlign":"center"
                                                                    }),
                width=11)
                ]),

            dbc.Row([
                dbc.Col(    
                html.Div(
                    id="effective_Section_image",
                    style={
                        "display": "flex",
                        "justifyContent": "center",   # horizontal centering
                        "alignItems": "center",       # optional, vertical centering
                        "overflowX": "auto",
                        "marginTop": "1px",
                    }
                ),
            width=12)]
        ),

            dbc.Row([
                dbc.Col(width=1),
                dbc.Col(
                    dcc.Markdown(id="done-message2", mathjax=True, style={
                                                                            "whiteSpace": "pre-wrap",
                                                                            "wordWrap": "break-word",
                                                                            "overflowX": "auto",
                                                                            "maxWidth": "100%",
                                                                            "fontSize": "20px",
                                                                            # "fontFamily": "monospace",
                                                                            # "textAlign":"center"
                                                                        }),
                    width=11)
                    ]),
                    dbc.Row([
                            dbc.Col(width=3),
                            dbc.Col(
                                html.Div(id="icr-table",style={
                                         "justifyContent": "center", "textAlign":"center"}),  # Placeholder for table content
                                width=5
                            )
                        ]),
                    dbc.Row([
                        dbc.Col(width=1),
                        dbc.Col(
                            dcc.Markdown(id="done-message3", mathjax=True, style={
                                                                            "whiteSpace": "pre-wrap",
                                                                            "wordWrap": "break-word",
                                                                            "overflowX": "auto",
                                                                            "maxWidth": "100%",
                                                                            "fontSize": "20px",
                                                                            "marginTop": "0px",  # reduces vertical gap
                                                                            # "fontFamily": "monospace",
                                                                            # "textAlign":"center"
                                                                        }),
                    width=11)
                    ]),                        
                    dbc.Row([
                            dbc.Col(width=1),
                            dbc.Col(
                                html.Div(id="moment-table",style={
                                         "justifyContent": "center", "textAlign":"center"}),  # Placeholder for table content
                                width=9
                            )
                        ]),    
                    dbc.Row([
                            dbc.Col(    
                            html.Div(
                                id="Icr_effective_Section_image1",
                                style={
                                    "display": "flex",
                                    "justifyContent": "center",   # horizontal centering
                                    "alignItems": "center",       # optional, vertical centering
                                    "overflowX": "auto",
                                    "marginTop": "1px",
                                    "marginTop": "-1820px",
                                }
                            ),
                        width=12)]),
                       
                    dbc.Row([dbc.Col(width=1),
                            dbc.Col(    
                            html.Div(
                                id="EquilbruimSection_image",
                                style={
                                    "display": "flex",
                                    "justifyContent": "center",   # horizontal centering
                                    "alignItems": "center",       # optional, vertical centering
                                    "overflowX": "auto",
                                    "marginTop": "1px",
                                    "marginTop": "0px",
                                }
                            ),
                        width=11)]),

                    dbc.Row([
                        dbc.Col(width=1),
                        dbc.Col(
                            dcc.Markdown(id="done-message4", mathjax=True, style={
                                                                            "whiteSpace": "pre-wrap",
                                                                            "wordWrap": "break-word",
                                                                            "overflowX": "auto",
                                                                            "maxWidth": "100%",
                                                                            "fontSize": "20px",
                                                                            "marginTop": "0px",  # reduces vertical gap
                                                                            # "fontFamily": "monospace",
                                                                            # "textAlign":"center"
                                                                        }),
                    width=11)
                    ]), 



                    dbc.Row([
                            dbc.Col(width=1),
                            dbc.Col(
                                html.Div(id="Mr-table",style={
                                         "justifyContent": "center", "textAlign":"center"}),  # Placeholder for table content
                                width=9
                            )
                        ]), 
                        
                    dcc.Store(id="moment-data"),
                    dcc.Store(id="axial-data"),
                    dcc.Store(id="Mt_F"),


    ])

@app.callback(
    Output("wall-image", "children"),
    Input("dropdown-t", "value"),
    Input("dropdown-S", "value"),
    Input("dropdown-bar", "value")
)
def update_wall_image(t, s, bar):
    if t is None:
        t = t_default
    if s is None:
        s = S_default
    if bar is None:
        bar = bar_default

    fig = draw_blocks_plotly(t_mm_actual=t*1000, s=s*1000, bar_diameter_mm=bar*1000)

    return dcc.Graph(figure=fig, config={"displayModeBar": False}, style={"height": "100%", "width": "100%", "padding": "0", "margin": "0", "backgroundColor": "transparent"})


@app.callback(
    Output("side-view", "children"),
    [Input("input-H", "value"),
     Input("input-P_DL", "value"),
     Input("input-P_LL", "value"),
     Input("input-P_S", "value"),
     Input("input-e", "value"),
     Input("input-W", "value")])
def update_side_view(H, P_DL, P_LL, P_S, e, W):
    if H is None:
        H = H_default
    if P_DL is None:
        P_DL = P_DL_default
    if P_LL is None:
        P_LL = P_LL_default
    if e is None:
        e = e_default
    if W is None:
        W = W_default
    if P_S is None:
        P_S = P_S_default

    fig = generate_side_view(H, P_DL, P_LL, P_S, e, W)
    return dcc.Graph(figure=fig, config={"displayModeBar": False}, style={"height": "100%", "width": "100%", "padding": "0", "margin": "0", "backgroundColor": "transparent"})



@app.callback(
    [Output("interaction-diagram", "children"),
     Output("moment-data", "data"),
     Output("axial-data", "data")],
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value")
)



def update_interaction_diagram(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):
    if not n_clicks:
        return "", [], []
        
    # Get cross section properties
    t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf=cross_section(t, S,bar,fblock)
    
    # Calculate maximum point
    PMax = solve_betaC(0.6,   fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t)
    betaC1 = float(PMax[0])
    
    # Calculate points for interaction diagram
    point2_results = calculate_point2(betaC1,faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t, d, num_points=15)
    point3_results, Mr_y, Pr_y, ey = calculate_point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t, d, Aseff_m, num_points=40)
    pure_moment = calculate_pure_moment(faim, Aseff_m, d, fm_g, bg_m, fm_ug, tf, bug_m_2, t)
    
    
    # Create arrays for plotting
    M = [0] + [PMax[-1]] + [pt[4] for pt in point2_results] + [pt[4] for pt in point3_results] + [pure_moment[4]]
    P = [PMax[-2]] + [PMax[-2]] + [pt[3] for pt in point2_results] + [pt[3] for pt in point3_results] + [pure_moment[3]]

    lambda_h = k *H/t 
     
    # Factored Loads
    P_SW_mid= rho_SW * H / 2  # at mid span
    P_F1 = 1.4 * (P_DL  + P_SW_mid)
    P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F4 = 1.25 * (P_DL  + P_SW_mid) 
    P_F5 = 0.9 * (P_DL  + P_SW_mid) 
    P_F6 = 1.25 * (P_DL  + P_SW_mid) + 1.5 *P_S # Snow load + DL
    P_F7 = 0.9 * (P_DL  + P_SW_mid) + 1.5* P_S   # Snow load + DL 

    P_F = [P_F1, P_F2, P_F3, P_F4, P_F5, P_F6, P_F7]

    Icr_results = []
    kd_vals = []
    for pf in P_F:
        I_cr_eff, kd_val = Icr_function(Aseff_m, d, beff_m_2, bg_m, E_m,  pf, tf)
        Icr_results.append(I_cr_eff)
        kd_vals.append(kd_val)

    ev, Mt_F_list, M_F_list, P_F_list,Pcr_list, Mag_list, Cm_list,EI_eff_raw, EI_eff_list, betad, Min_EIeff , Max_EIeff, I_cr_array, betad_list= Moment_Calculation (t,e ,H, rho_SW, W, P_DL,P_LL, P_S,I_gross_eff, E_m, ek, Icr_results )

    
    # Create the figure
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=M, y=P, mode="lines", name="Envelope"))

    fig.add_trace(go.Scatter(x=Mt_F_list, y=P_F_list, mode="markers", name="Total Moments"))
    fig.add_trace(go.Scatter(x=M_F_list, y=P_F_list, mode="markers",marker=dict(color='blue',symbol='x'), name="Primary Moment"))
    fig.update_layout(
        xaxis_title="Moment (kN⋅m)",
        yaxis_title="Axial Force (kN)",plot_bgcolor="white",
        width=800,       
        height=650,         
        xaxis=dict(gridcolor="lightgrey", range=[0, 1.1*max(max(M),max(Mt_F_list))],   
                   zeroline=True, zerolinecolor="black",zerolinewidth=1),
        yaxis=dict(gridcolor="lightgrey", range=[0, 1.1*max(max(P),max(P_F_list))],
                   zeroline=True, zerolinecolor="black",zerolinewidth=1),
        showlegend=True,
        margin=dict(l=0, r=0, b=0, t=0.1),  # Remove margins
                      paper_bgcolor='rgba(0,0,0,0)',
                          font=dict(
        size=20,  # Global font size for all text (titles, labels, legend)
        family="Times New Roman",  
    ),
    )
    fig.update_layout(
    legend=dict(
        x=0.98,
        y=0.98,
        xanchor='right',
        yanchor='top',
        bgcolor='rgba(255,255,255,0.7)',  # semi-transparent background

        font=dict(size=20)
    )
    )
    
    return dcc.Graph(figure=fig), M, P

@app.callback(
    Output("done-message", "children"),
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)
def show_done_message(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):
    if not n_clicks:
        return """
        return """
# Get cross section properties
    t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf=cross_section(t, S,bar,fblock)

    t *= 1000
    beff_m_1 *= 1000
    beff_m_2 *= 1000
    bg_m *= 1000
    bug_m_1 *= 1000
    bug_m_2 *= 1000
    kd *= 1000
    ek *= 1000
    tf *= 1000
    S *= 1000
    bar *= 1000
    A_gr *= 1e6
    A_ug_1 *= 1e6
    A_ug_2 *= 1e6
    Ae_1 *= 1e6
    Ae_2 *= 1e6
    As *= 1e6
    Aseff_m *= 1e6

    I_gross_gr *= 1e12
    I_gross_ug_1 *= 1e12
    I_gross_eff *= 1e12
    I_cr_eff *= 1e12

    h_over_t = H / t *1000

    # Conditional check
    if t < 140:
        t_check = f"**Thickness of the block ($t$):** {t:.0f} mm  ❌"
    else:
        t_check = f"**Thickness of the block ($t$):** {t:.0f} mm  ✅"

    report = f"""
## 🧱 ** Masonry Wall Summary Report **
This report summarizes the design of a masonry wall based on the provided inputs.  The design is based on the **CSA S304-24** standard.  The wall is designed to resist out-of-plane loads, including external dead loads, self weight, live loads, and wind loads. The design also considers the slenderness effects and the interaction between axial loads and moments. The design procdure follows Masonry Structures Behavior and Design 2nd Edition by Bennett & Drysdale. 
The inputs loads are unfactored. The Load combinations are obtained from **NBC 2020** as follows:
- **Load Combination 1**: $1.4\\, P_\\text{{DL}}$ 
- **Load Combination 2**: $1.25\\, P_\\text{{DL}} + 1.5\\, P_\\text{{LL}}$
- **Load Combination 3**: $0.9\\, P_\\text{{DL}} +  1.5\\, P_\\text{{LL}}$
- **Load Combination 4**: $1.25\\, P_\\text{{DL}}  + 1.4\\, W$ 
- **Load Combination 5**: $0.9\\, P_\\text{{DL}}  + 1.4\\, W$
- **Load Combination 6**: $1.25\\, P_\\text{{DL}} + 1.5\\, P_\\text{{S}}$
- **Load Combination 7**: $0.9\\, P_\\text{{DL}} + 1.5\\, P_\\text{{S}}$
This tool is intended for educational purposes only and should not be used for actual design without consulting a qualified engineer.  We are not responsible for any errors or omissions in the design.  Please use this tool at your own risk.
### **Wall Properties**
- **Height**: ${H:.0f}\\,\\text{{m}}$
- **Thickness**: ${t:.0f}\\,\\text{{mm}}$
- **Block Strength**: ${fblock:.0f}\\,\\text{{MPa}}$
- **Bar Spacing**: ${S:.0f}\\,\\text{{mm}}$
- **Bar Diameter**: ${bar:.0f}\\,\\text{{M}}$

### **Material Reduction Factors**
- $\\phi_\\text{{masonry}} = 0.6$
- $\\phi_\\text{{reinforcement}} = 0.85$

### **Applied Loads**
- **Dead Load ($P_\\text{{DL}}$)**: ${P_DL:.2f}\\,\\text{{kN}}$
- **Live Load ($P_\\text{{LL}}$)**: ${P_LL:.2f}\\,\\text{{kN}}$
- **Snow Load ($P_\\text{{S}}$)**: ${P_S:.2f}\\,\\text{{kN}}$
- **Eccentricity ($e$)**: ${e:.1f}\\,\\text{{mm}}$
- **Wind Load ($W$)**: ${W:.2f}\\,\\text{{kPa}}$
"""

    # Slenderness Check Section
    slenderness_check_section = f"""
###  **Slenderness Check**

- $\\dfrac{{h}}{{t}} = {h_over_t:.2f}$

**According to CSA S304-24:**  
$h/t > 30 \\Rightarrow \\text{{Slenderness effects must be considered}}$

In accordance with relevant provisions:

- The minimum thickness of the block is 140 mm (C11.7.4.6.2).
- Boundary conditions are assumed to be pinned unless a more detailed analysis is provided (C11.7.4.6.3).
- The area of the reinforcement must be less than the balanced area to ensure yielding of steel before masonry crushing (C11.7.4.6.4).
- The moment magnification method may be used, with a moment diagram factor $C_m = 1$ and effective height factor $k = 1$ (C11.7.4.6.3).

**Based on these provisions:**

- {t_check}
- **Boundary conditions:** pinned at the base and roller-supported at the top.✅
- **Moment diagram factor ($C_m$):** 1  
- **Effective height factor ($k$):** 1
"""
    effective_section = f"""
### **Effective Section Properties**
The width of the effective compression zone ($b$) is determined based on CSA S304-24 (C 11.6.1), using the spacing between mild reinforcement ($s$) and the thickness of the wall ($t$):

$$
b = \\min\\left\\{{
\\begin{{array}}{{ll}}
6t &= 6 \\times {t:.0f} = {6*t:.0f}\\,\\text{{mm}} \\\\
s &= {S:.0f}\\,\\text{{mm}}
\\end{{array}}
\\right\\}} = {beff_m_2:.0f}\\,\\text{{mm}}
$$

The effective width per unit length, $$b_{{\\text{{eff}}}}$$ is:

$$
b_{{\\text{{eff}}}} = \\frac{{b}}{{s}} = \\frac{{{beff_m_2:.0f}}}{{{S:.0f}}}\\times 1000  = {beff_m_2:.0f}\\text{{mm}}
$$

The effective compression zone consists of grouted ($b_{{gr_{{eff}}}}$) and ungrouted ($b_{{ug_{{eff}}}}$) portions, which can be determined as:

$$
b_{{gr_{{eff}}}} = b_g \\times \\frac{{1000\\,mm/m}}{{s}} = 200 \\times \\frac{{1000\\,mm/m}}{{{S:.0f}}} = {bg_m:.0f}\\,mm/m
$$

$$
b_{{ug_{{eff}}}} = b_{{eff}} - b_{{gr_{{eff}}}} = {beff_m_2:.0f} - {bg_m:.0f} = {bug_m_2:.0f}\\,mm
$$

With the effective compression zone length known, the effective area, $A_e$, can be determined. It includes grouted ($A_{{gr_{{eff}}}}$) and ungrouted ($A_{{ug_{{eff}}}}$) areas:

$$
A_{{gr_{{eff}}}} = b_{{gr_{{eff}}}} \\times t = {bg_m:.0f} \\times {t:.0f} = {A_gr:.0f}\\,mm^2
$$

$$
A_{{ug_{{eff}}}} = b_{{ug_{{eff}}}} \\times 2t_f = {bug_m_2:.0f} \\times (2 \\times {tf:.0f}) = {A_ug_2:.0f}\,\\text{{mm}}^2
$$

$$
A_e = A_{{gr_{{eff}}}} + A_{{ug_{{eff}}}} = {A_gr:.0f} + {A_ug_2:.0f} = {Ae_2:.0f}\\,mm^2
$$


The  effective compressive strength of the masonry assembly ($$f_{{m_{{eff}}}}$$) can be determined as: 
$$
f_{{m_{{eff}}}} = \\frac{{A_{{gr_{{eff}}}}}}{{A_e}} \\times f_{{m_{{gr}}}} + \\frac{{A_{{ug_{{eff}}}}}}{{A_e}} \\times f_{{m_{{ug}}}}
$$

$$
f_{{m_{{eff}}}} = \\frac{{{A_gr:.0f}}}{{{Ae_2:.0f}}} \\times {fm_g:.0f} + \\frac{{{A_ug_2:.2f}}}{{{Ae_2:.2f}}} \\times {fm_ug:.2f} = {fm_e_2:.2f}\\,\\text{{MPa}}
$$

The effective reinforcement area:

$$
A_{{s_{{eff}}}} = A_s \\times \\frac{{1000\\,mm/m}}{{s}} = {As:.0f} \\times \\frac{{1000}}{{{S:.0f}}} = {Aseff_m:.2f}\\,mm^2/m
$$



The effective cross-section per one-meter width is illustrated in the Figure below.
"""
    report += slenderness_check_section+effective_section
    return report





@app.callback(
    Output("effective_Section_image", "children"),
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)
def effective_Section_image (n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):
    if S == 200 *mm:
        img = Image.open(r"H:\My Drive\Python Codes to help you\MasonryWalls_ID App\assets\Effective Section _ FG.png")
    else:
        img = Image.open(r"H:\My Drive\Python Codes to help you\MasonryWalls_ID App\assets\Effective Section _ PG.png")

    width, height = img.size
    # Get cross section properties
    t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf=cross_section(t, S,bar,fblock)

    t *= 1000
    beff_m_1 *= 1000
    beff_m_2 *= 1000
    bg_m *= 1000
    bug_m_1 *= 1000
    bug_m_2 *= 1000
    kd *= 1000
    ek *= 1000
    tf *= 1000
    S *= 1000

    A_gr *= 1e6
    A_ug_1 *= 1e6
    A_ug_2 *= 1e6
    Ae_1 *= 1e6
    Ae_2 *= 1e6
    As *= 1e6
    Aseff_m *= 1e6

    I_gross_gr *= 1e12
    I_gross_ug_1 *= 1e12
    I_gross_eff *= 1e12
    I_cr_eff *= 1e12

    tf_text = f"t𝒻 = {tf:.1f} mm"
    Aseff_text = f"Aₛₑ𝒻𝒻 = {Aseff_m:.0f} mm²"
    beff_text = f"bₑ𝒻𝒻 = {beff_m_1:.0f} mm"
    t_text = f"t = {t:.0f} mm"
    bg_text = f"bg = {bg_m:.0f} mm"

    fig = go.Figure()

    # Add image as background
    fig.update_layout(
        images=[dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=height,  # Top-left corner
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below"
        )],
        xaxis=dict(visible=False, range=[0, width]),
        yaxis=dict(visible=False, range=[0, height], scaleanchor="x"),
        width=width,
        height=height,
        margin=dict(l=0, r=0, t=0, b=0)
    )

    # Add text annotation
    fig.add_annotation(
        x=780, y=80,  # Adjust based on image dimensions
        text=tf_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black"))


    fig.add_annotation(
        x=765, y=390,  # Adjust based on image dimensions
        text=Aseff_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=15, color="black"))

    fig.add_annotation(
        x=500, y=550,  # Adjust based on image dimensions
        text=beff_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black")
        
    )

    fig.add_annotation(
        x=50, y=300,  # Adjust based on image dimensions
        text=t_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black"), textangle=270
        
    )
    fig.add_annotation(
        x=500, y=120,  # Adjust based on image dimensions
        text=bg_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black")
        
    )

    scale_factor = 0.5  # 50% of original size

    fig.update_layout(
        images=[dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=height,
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below"
        )],
        xaxis=dict(visible=False, range=[0, width]),
        yaxis=dict(visible=False, range=[0, height], scaleanchor="x"),
        width=int(width * scale_factor),
        height=int(height * scale_factor),
        margin=dict(l=0, r=0, t=0, b=0)
    )
 

    return dcc.Graph(figure=fig)


@app.callback(
    [Output("done-message2", "children"),
     Output("icr-table", "children")],
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)
def show_done_message2(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):
    if not n_clicks:
        return """
        return """
# Get cross section properties
    t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf=cross_section(t, S,bar,fblock)
    # print("t",t)
    # print("beff_m_1",beff_m_1)
    # print("beff_m_2",beff_m_2)
    # print("As",As)
    # print("Aseff_m",Aseff_m)
    # print("bg_m",bg_m)
    # print("bug_m_1",bug_m_1)
    # print("bug_m_2",bug_m_2)
    # print("A_gr",A_gr)
    # print("A_ug_1",A_ug_1)
    # print("A_ug_2",A_ug_2)
    # print("Ae_1",Ae_1)
    # print("Ae_2",Ae_2)
    # print("fm_e_1",fm_e_1)
    # print("fm_e_2",fm_e_2)
    # print("I_gross_gr",I_gross_gr)
    # print("I_gross_ug_1",I_gross_ug_1)
    # print("I_gross_eff",I_gross_eff)
    # print("I_cr_eff",I_cr_eff)
    # print("kd",kd)
    # print("n",n)
    # print("E_m",E_m)

    d= t/2
    S_eff = I_gross_eff / (t/2)
    P_SW_mid= rho_SW * H / 2  # at mid span

    P_F1 = 1.4 * (P_DL  + P_SW_mid)
    P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F4 = 1.25 * (P_DL  + P_SW_mid) 
    P_F5 = 0.9 * (P_DL  + P_SW_mid) 
    P_F6 = 1.25 * (P_DL  + P_SW_mid) + 1.5 *P_S # Snow load + DL
    P_F7 = 0.9 * (P_DL  + P_SW_mid) + 1.5* P_S   # Snow load + DL 

    P_F = [P_F1, P_F2, P_F3, P_F4, P_F5, P_F6, P_F7]

    Icr_results = []
    kd_vals = []
    for pf in P_F:
        I_cr_eff, kd_val = Icr_function(Aseff_m, d, beff_m_2, bg_m, E_m,  pf, tf)
        Icr_results.append(I_cr_eff)
        kd_vals.append(kd_val)

    ev, Mt_F_list, M_F_list, P_F_list,Pcr_list, Mag_list, Cm_list,EI_eff_raw, EI_eff_list, betad, Min_EIeff , Max_EIeff, I_cr_array, betad_list= Moment_Calculation (t,e ,H, rho_SW, W, P_DL,P_LL, P_S, I_gross_eff, E_m, ek, Icr_results )

    h_over_t = H / t
    report = f"""
        The gross moment of inertia ($$I_o$$) consists of contributions from grouted cells ($$I_{{o_{{gr}}}}$$) and hollow cells ($$I_{{o_{{ug}}}}$$). 
        $$
        I_{{o_{{gr}}}} = b_{{eff}} \\times t ^3 / 12 
        $$

        $$
        I_{{o_{{gr}}}} =  {beff_m_1*1000:.0f} \\times {t*1000:.0f}^3 / 12 = {I_gross_gr*10**12:.0f}\\,mm^4
        $$

        $$ 
        I_{{o_{{ug}}}} = \\left( b_{{eff}} \\times t_f^3 / 12 + b_{{eff}} \\times t_f \\times (t/2 - t_f/2)^2 \\right) \\times 2
        $$

        $$ 
        I_{{o_{{ug}}}} = \\left( {beff_m_1*1000:.0f} \\times {tf*1000:.0f}^3 / 12 + {beff_m_1*1000:.0f} \\times {tf*1000:.0f} \\times ({t*1000:.0f}/2 - {tf*1000:.0f}/2)^2 \\right) \\times 2 = {I_gross_ug_1*10**12:.0f}\\,mm^4
        $$

        The effective gross moment of inertia is then calculated as:
                 
        $$
        I_{{\\text{{o}}}} = \\frac{{A_{{gr_{{eff}}}}}}{{A_e}} \\times I_{{o_{{gr}}}} + \\frac{{A_{{ug_{{eff}}}}}}{{A_e}} \\times I_{{o_{{ug}}}} 
        $$
        $$
        I_{{\\text{{o}}}} = \\frac{{{A_gr*10**6:.0f}}}{{{Ae_2*10**6:.0f}}} \\times {I_gross_gr*10**12:.0f}  + \\frac{{{A_ug_2*10**6:.0f}}}{{{Ae_2*10**6:.0f}}} \\times {I_gross_ug_1*10**12:.0f} = {I_gross_eff*10**12:.0f}\\,mm^4
        $$

       The section modulus ($$S_e$$) is:
       $$
        S_e = \\frac{{I_{{\\text{{o}}}}}}{{0.5t}} 
        $$
        $$
        S_e = \\frac{{{I_gross_eff*10**12:.0f}}}{{0.5\\times {t*10**3:.0f}}} = {S_eff*10**9:.0f}\\,mm^3
        $$

        Finally, the kern eccentricity ($$e_k$$) is calculated as:
        $$
        e_k = \\frac{{S_e}}{{A_e}}
        $$ 
        $$
        e_k = \\frac{{{S_eff*10**9:.0f}}}{{{Ae_2*10**6:.0f}}} = {ek*1000:.0f}\\,mm
        $$ 

        ### **Calculation of Self-weight at Midspan**
        For slender walls with For slender walls with $$h/t>30$$, CSA S304-24 (C 11.7.4.6.5) states that the self-weight must be included in the 
        total moment calculation, including both the primary moment due to wind and eccentric vertical loads, and the secondary moment 
        induced by vertical loads such as the self-weight and external vertical loads.  Therefore, it is necessary to calculate the self-weight.  

        Table B.1 in Drysdale has the weights per unit surface area of fully grouted blocks ($${{{rho_g:.2f}}}kN/m^2$$) and hollow blocks ($${{{rho_ug:.2f}}}kN/m^2$$).  

        Therefore, the  self  weight, $$P_{{sw}}$$, at midspan can be calculated as:
        $$ 
        P_{{sw}} = (\\frac{{A_{{gr_{{eff}}}}}}{{A_e}} \\times {rho_g:.2f} + \\frac{{A_{{ug_{{eff}}}}}}{{A_e}} \\times {rho_ug:.2f}) \\times h/2
        $$

        $$
        P_{{sw}} = (\\frac{{{A_gr*10**6:.0f}}}{{{Ae_2*10**6:.0f}}} \\times {rho_g:.2f} + \\frac{{{A_ug_2*10**6:.0f}}}{{{Ae_2*10**6:.0f}}} \\times {rho_ug:.2f}) \\times {H:.0f}/2= {rho_SW*H*0.5:.2f} kN/m
        $$

        """
    
    Cracked_Inertia_section = f"""
     ### ** Calculation of Cracked Moment of Inertia**
        CSA S304-24 (C11.7.4.4) defines cracked moment of inertia ($$I_{{cr}}$$) as:
        the transformed moment of inertia of the cracked section including the effects of the factored axial load, assuming a triangle stress 
        distribution in the masonry and strain in the extreme layer of tension reinforcement equal to the strain in the extreme layer of 
        tension reinforcement $$\\varepsilon_s = \\frac{{f_y}}{{E_s}}$$

        From strain comptability, the compressive strain in the masonry, \( $$\\varepsilon_m$$), is given by:

        $$
        \\varepsilon_m = \\frac{{\\varepsilon_y}}{{d - kd}} \\cdot kd
        $$

        The corresponding masonry stress, \\( $$f_m$$ \\), is:
        $$
        f_m = \\frac{{\\varepsilon_y}}{{d - kd}} \\cdot kd \\cdot E_m
        $$

        Equilibrium of internal forces requires:
        $$
        C_m = P + f_y A_s
        $$

        $$
        \\frac{1}{2} \\times kd \\times f_m \\times b_{{eff}} = P + f_y A_{{s_{{eff}}}}
        $$

        Solving the above equation allows for the determination of the neutral axis depth, \\( $$kd$$ \\).

        Once \\( $$kd$$ \\) is known, the cracked moment of inertia is calculated as:

        $$
        I_{{cr}} = \\frac{{b \\times kd^3}}{{3}} + n \\times A_{{s_{{eff}}}} \\times (d - kd)^2
        $$

        Since there are different load cases with different axial loads, the cracked moment of inertia is calculated for each load case.  The 
        following table summarizes the results:

    """
    input_data = [
    ("Aseff", f"{Aseff_m/mm**2:.2f}", "mm"),
    ("d", f"{d/mm:.3f}", "mm"),
    ("beff", f"{beff_m_2/mm:.3f}", "mm"),
    ("beffg", f"{bg_m/mm:.3f}", "mm"),
    ("Em", f"{E_m:.0f}", "MPa"),
    ("tf", f"{tf/mm:.2f}", "mm")
]

    input_table = dbc.Table([
        html.Thead(html.Tr([
            html.Th("Parameter"),
            html.Th("Value"),
            html.Th("Unit")
        ])),
        html.Tbody([
            html.Tr([
                html.Td(name),
                html.Td(value),
                html.Td(unit)
            ]) for name, value, unit in input_data
        ])
    ], bordered=True, size="sm", style={
        "fontSize": "14px",
        "marginTop": "15px",
        "maxWidth": "800px"
    })

    # Build table
    # Create new header with PF
    table_header = [
            html.Thead(html.Tr([
                html.Th("Load Case"),
                html.Th("kd [mm]"),
                html.Th("Icr [mm⁴]"),
                html.Th("P_F [kN]")  # 🔹 New column
            ]))
        ]

    # Create table body
    table_body = [
        html.Tbody([
            html.Tr([
                html.Td(str(i + 1)),
                html.Td(f"{kd_vals[i] * 1000:.2f}"),
                html.Td(f"{Icr_results[i]:.2e}"),
                html.Td(f"{P_F[i]:.2f}")
            ]) for i in range(len(Icr_results))
        ])
                    ]

    icr_table = dbc.Table(
    # Table content
    table_header + table_body,
    bordered=True,
    striped=True,
    hover=True,
    responsive=True,
    size="sm",  # 🔹 This makes the table more compact (smaller padding)
    style={
        "fontSize": "20px",      # 🔹 Smaller font
        "marginTop": "10px",
        "maxWidth": "700px",     # 🔹 Optional: restrict width
        "textAlign": "center"    # 🔹 Center content in cells
    }
    )
    
    return report + Cracked_Inertia_section, icr_table


@app.callback(
    [Output("done-message3", "children"),
     Output("moment-table", "children"),Output("Mt_F", "data")],
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)

def show_done_message3(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):
    if not n_clicks:
        return """""",[]
# Get cross section properties
    t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf=cross_section(t, S,bar,fblock)

    d= t/2
    S_eff = I_gross_eff / (t/2)
    P_SW_mid= rho_SW * H / 2  # at mid span

    P_F1 = 1.4 * (P_DL  + P_SW_mid)
    P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F4 = 1.25 * (P_DL  + P_SW_mid) 
    P_F5 = 0.9 * (P_DL  + P_SW_mid)
    P_F6 = 1.25 * (P_DL  + P_SW_mid) + 1.5 *P_S # Snow load + DL
    P_F7 = 0.9 * (P_DL  + P_SW_mid) + 1.5* P_S   # Snow load + DL  

    P_F = [P_F1, P_F2, P_F3, P_F4, P_F5, P_F6, P_F7]

    Icr_results = []
    kd_vals = []
    for pf in P_F:
        I_cr_eff, kd_val = Icr_function(Aseff_m, d, beff_m_2, bg_m, E_m,  pf, tf)
        Icr_results.append(I_cr_eff)
        kd_vals.append(kd_val)

    ev, Mt_F_list, M_F_list, P_F_list,Pcr_list, Mag_list, Cm_list,EI_eff_raw, EI_eff_list, betad, Min_EIeff , Max_EIeff, I_cr_array, betad_list= Moment_Calculation (t,e ,H, rho_SW, W, P_DL,P_LL, P_S, I_gross_eff, E_m, ek, Icr_results )

    h_over_t = H / t
    Icr_Image= """The figure below illustrate the internal forces for cracked moment of inertia calculation. 

       .
        
       
       
       
       
       
       
       
              
      . 


    """
        # Add new section for moment magnification
    moment_magnification_section = f"""
### **Determining the Total Applied Moment Magnification (C11.7.4.3)**

To determine the total applied moment, the moment magnification factor ($$\\Psi$$) is multiplied by the
primary factored moment ($$M_1$$), calculated as:


$$
M_1 = \\min\\left\\{{
\\begin{{array}}{{ll}}
w_F \\frac{{h^2}}{{8}} + P_F\\, e/2  \\\\
0.1 \\,P_F\\, t\\, (C11.7.2)
\\end{{array}}
\\right\\}}
$$

To obtain $$\\Psi$$, the critical buckling load ($$P_{{cr}}$$) shall be calculated using the effective stiffness ($$EI_{{eff}}$$).
The effective stiffness is calculated as:

$$
EI_{{eff}} = E_m \\left[0.25 I_{{cr}} - \\left(0.25 I_{{cr}} - I_{{cr}}\\right) \\left(\\frac{{e-e_k}}{{2e_k}}\\right)\\right] \\geq E_m I_{{cr}}
$$

The virtual eccentricity $$e_v$$ is calculated as:

$$
e_v = \\frac{{M_1}}{{P_F}}
$$

As per CSA S304-24 (C11.7.4.4), the effective stiffness must satisfy:

* Minimum limit:
    $$
    EI_{{eff}} \\geq E_{{m}} I_{{cr}}
    $$

* Maximum limit:
    $$
    EI_{{eff}} \\leq 0.25 E_{{m}} I_{{o}}
    $$

The critical buckling load can now be obtained as:

$$
P_{{cr}} = \\frac{{\\pi^2 \\cdot EI_{{eff}}}}{{(1+0.5\\beta_d)(kh)^2}}
$$

Finally, the moment magnification factor ($$\\Psi_{{sw}}$$) can be calculated using:

$$
\\Psi_{{sw}} = \\frac{{C_m}}{{1-\\frac{{P_F}}{{P_{{cr}}}}}}
$$

The total factored moment can now be determined as:

$$
M_{{total}} = \\Psi M_F
$$

The following table summarizes the results for each load case:
    """
    # Table summarizing full Moment Calculation results
    moment_calc_header = [
        html.Thead(html.Tr([
            html.Th("Load Case"),
            html.Th("Mₜ [kNm]"),
            html.Th("M₁ [kNm]"),
            html.Th("Pf [kN]"),
            html.Th("P꜀ᵣ [kN]"),
            html.Th("Ψ (Mag. Factor)"),
            html.Th("Cₘ"),
            html.Th("EIₑff [kNm²]"),
            html.Th("EIₑff_Eq [kNm²]"),
            html.Th("Min EIₑff [kNm²]"),
            html.Th("Max EIₑff [kNm²]"),
            html.Th("βd")
        ]))
    ]

    moment_calc_body = html.Tbody([
        html.Tr([
            html.Td(str(i + 1)),
            html.Td(f"{Mt_F_list[i]:.2f}"),
            html.Td(f"{M_F_list[i]:.2f}"),
            html.Td(f"{P_F_list[i]:.2f}"),
            html.Td(f"{Pcr_list[i]:.2f}"),
            html.Td(f"{Mag_list[i]:.2f}"),
            html.Td(f"{Cm_list[i]:.2f}"),
            html.Td(f"{EI_eff_list[i] :.2f}"),
            html.Td(f"{EI_eff_raw[i] :.2f}"),
            html.Td(f"{Min_EIeff[i] :.2f}"),
            html.Td(f"{Max_EIeff[i] :.2f}"),
            html.Td(f"{betad_list[i] :.2f}")    
        ]) for i in range(len(P_F_list))
    ])

    moment_table = dbc.Table(
        moment_calc_header + [moment_calc_body],
        bordered=True,
        striped=True,
        hover=True,
        responsive=True,
        size="sm",
        style={
            "fontSize": "20px",
            "marginTop": "1px",
            "maxWidth": "1000px",
            "textAlign": "center"
        }
    )
    
    return [
    Icr_Image + moment_magnification_section,
    html.Div([
        # html.Br(), html.Br(),  # space before table
        moment_table
    ]), Mt_F_list
]






@app.callback(
    Output("Icr_effective_Section_image1", "children"),
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)

def Icr_effective_Section_image1 (n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):

    # Get cross section properties
    t, beff_m_1, beff_m_2, As,Aseff_m, bg_m, bug_m_1, bug_m_2,A_gr,A_ug_1,A_ug_2 , Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, I_cr_eff, kd, n , E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf=cross_section(t, S,bar,fblock)
    
    if S == 200*mm:
        img = Image.open(r"H:\My Drive\Python Codes to help you\MasonryWalls_ID App\assets\Icr_insideFaceshell_Fully Grouted.png")
    else:
        img = Image.open(r"H:\My Drive\Python Codes to help you\MasonryWalls_ID App\assets\Icr_insideFaceshell.png")

    width, height = img.size
    t *= 1000
    beff_m_1 *= 1000
    beff_m_2 *= 1000
    bg_m *= 1000
    bug_m_1 *= 1000
    bug_m_2 *= 1000
    kd *= 1000
    ek *= 1000
    tf *= 1000
    S *= 1000

    A_gr *= 1e6
    A_ug_1 *= 1e6
    A_ug_2 *= 1e6
    Ae_1 *= 1e6
    Ae_2 *= 1e6
    As *= 1e6
    Aseff_m *= 1e6

    I_gross_gr *= 1e12
    I_gross_ug_1 *= 1e12
    I_gross_eff *= 1e12
    I_cr_eff *= 1e12

    tf_text = f"t𝒻 = {tf:.1f} mm"
    Aseff_text = f"Aₛₑ𝒻𝒻 = {Aseff_m:.0f} mm²"
    beff_text = f"bₑ𝒻𝒻 = {beff_m_1:.0f} mm"
    t_text = f"t/2 = {t/2:.0f} mm"
    bg_text = f"bg = {bg_m:.0f} mm"
    fig = go.Figure()

    # Add image as background
    fig.update_layout(
        images=[dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=height,  # Top-left corner
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below"
        )],
        xaxis=dict(visible=False, range=[0, width]),
        yaxis=dict(visible=False, range=[0, height], scaleanchor="x"),
        width=width,
        height=height,
        margin=dict(l=0, r=0, t=0, b=0)
    )

# Add text annotation
    fig.add_annotation(
        x=850, y=150,  # Adjust based on image dimensions
        text=Aseff_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black"))

    # Add text annotation
    fig.add_annotation(
        x=220, y=600,  # Adjust based on image dimensions
        text=tf_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black"),textangle=270)

    # Add text annotation
    fig.add_annotation(
        x=120, y=480,  # Adjust based on image dimensions
        text=t_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black"),textangle=270)

    fig.add_annotation(
        x=840, y=810,  # Adjust based on image dimensions
        text=beff_text,
        showarrow=False,
        arrowhead=2,
        font=dict(size=20, color="black"))

    scale_factor = 0.35 # 50% of original size

    fig.update_layout(
        images=[dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=height,
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below"
        )],
        xaxis=dict(visible=False, range=[0, width]),
        yaxis=dict(visible=False, range=[0, height], scaleanchor="x"),
        width=int(width * scale_factor),
        height=int(height * scale_factor),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    return dcc.Graph(figure=fig)


@app.callback(
    Output("EquilbruimSection_image", "children"),
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)

def EquilbruimSection_image (n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):

    if S == 200*mm:
         img = Image.open(r"H:\My Drive\Python Codes to help you\MasonryWalls_ID App\assets\EquilibruimSection_FG.png")  
    else:         
        img = Image.open(r"H:\My Drive\Python Codes to help you\MasonryWalls_ID App\assets\EquilibruimSection_PG.png")

    fig = go.Figure()

    width, height = img.size
    # Add image as background
    fig.update_layout(
        images=[dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=height,  # Top-left corner
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below"
        )],
        xaxis=dict(visible=False, range=[0, width]),
        yaxis=dict(visible=False, range=[0, height], scaleanchor="x"),
        width=width,
        height=height,
        margin=dict(l=0, r=0, t=0, b=0)
    )
    scale_factor = 0.35 # 50% of original size

    fig.update_layout(
        images=[dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=height,
            sizex=width,
            sizey=height,
            sizing="stretch",
            layer="below"
        )],
        xaxis=dict(visible=False, range=[0, width]),
        yaxis=dict(visible=False, range=[0, height], scaleanchor="x"),
        width=int(width * scale_factor),
        height=int(height * scale_factor),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    Moment_resistance_section_header = f"""
    ### **Determining Moment Resistance ($$M_r$$)**
    For illustration purpose, the following figure illustrates a case where the neutral axis exits in the faceshell. (Note that the NA axis locastion will change the calculation process):"""
    return html.Div([
                    dcc.Markdown(Moment_resistance_section_header, mathjax=True,style={
            "fontSize": "20px"
            
        }),
                    dcc.Graph(figure=fig)
                        ])


@app.callback(
    [Output("done-message4", "children"),
     Output("Mr-table", "children")],
    Input("btn-check", "n_clicks"),
    State("input-H", "value"),
    State("dropdown-t", "value"),
    State("dropdown-fblock", "value"),
    State("dropdown-S", "value"),
    State("dropdown-bar", "value"),
    State("input-P_DL", "value"),
    State("input-P_LL", "value"),
    State("input-P_S", "value"),
    State("input-e", "value"),
    State("input-W", "value"),
    prevent_initial_call=True
)

def Moment_resistanceText(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W):

    Moment_resistance_section = f"""
    

    The internal force equilibrium equation is expressed as:

    $$
    C_m = T_s + P_F
    $$

    Assuming that the mild steel has yielded and that the neutral axis depth $c \leq t_f$, the compression force in the masonry ($C_m$) is:

    $$
    C_m = \phi_m \\alpha_1 f_m^{{\\prime}} t_f b_{{eff}}
    $$

    The tensile force in the mild reinforcement ($T_s$) is:

    $$
    T_s = \phi_s A_s f_y
    $$

    Substituting all forces into the equilibrium equation:

    $$
    C_m = T_s + P_F
    $$

    $$
    \phi_m \\alpha_1 f_m^{{\\prime}} t_f b_{{eff}} = \phi_s A_s f_y + P_F
    $$

    Solving gives:

    $$
    c = \\frac{{\phi_s A_s f_y + P_F}}{{\phi_m \\alpha_1 f_m^{{\\prime}} b_{{eff}}}}
    $$

    Finally, the moment resistance is determined by summing moments about the extreme compression fiber (point o):

    $$
    M_r = -C_m \\frac{{\\beta_1 c}}{{2}} + T_s d + T_{{PT}} d_{{PT}} + (P_F) \\frac{{t}}{{2}}
    $$
    
    The following table $$Mr$$ the results for each load case:
    """
        # Factored Loads

    rho_SW=cross_section(t, S,bar,fblock)[23]

    P_SW_mid= rho_SW * H / 2  # at mid span
    P_F1 = 1.4 * (P_DL  + P_SW_mid)
    P_F2 = 1.25 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F3 = 0.9 * (P_DL  + P_SW_mid) + 1.5 * P_LL 
    P_F4 = 1.25 * (P_DL  + P_SW_mid) 
    P_F5 = 0.9 * (P_DL  + P_SW_mid) 
    P_F6 = 1.25 * (P_DL  + P_SW_mid) + 1.5 *P_S # Snow load + DL
    P_F7 = 0.9 * (P_DL  + P_SW_mid) + 1.5* P_S   # Snow load + DL 

    P_F = [P_F1, P_F2, P_F3, P_F4, P_F5, P_F6, P_F7]

    M=update_interaction_diagram(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W)[1]
    P=update_interaction_diagram(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W)[2]
    P_sorted, M_sorted = zip(*sorted(zip(P, M)))  # sort both lists based on P
    # Now interpolate moments for the given factored axial loads
    M_interpolated = np.interp(P_F, P_sorted, M_sorted)
    M_F_list = show_done_message3(n_clicks, H, t, fblock, S, bar, P_DL, P_LL, P_S, e, W)[2]
    # Table summarizing full Moment Resistance results
    moment_calc_header = [
        html.Thead(html.Tr([
            html.Th("Load Case"),
            html.Th("Pf [kNm]"),
            html.Th("Mₜ [kNm]"),
            html.Th("Mr [kNm]"),
        html.Th("Result")
    ]))
]


    moment_calc_body = html.Tbody([
        html.Tr([
            html.Td(str(i + 1)),
            html.Td(f"{P_F[i]:.2f}"),
            html.Td(f"{M_F_list[i]:.2f}"),
            html.Td(f"{M_interpolated[i]:.2f}"),
            html.Td(
            "✅ Pass" if M_interpolated[i] >= M_F_list[i] else "❌Fail",
            style={
                "color": "green" if M_interpolated[i] >= M_F_list[i] else "red",
                "fontWeight": "bold"
            })
        ]) for i in range(len(P_F))
    ])

    Mr_table = dbc.Table(
        moment_calc_header + [moment_calc_body],
        bordered=True,
        striped=True,
        hover=True,
        responsive=True,
        size="sm",
        style={
            "fontSize": "20px",
            "marginTop": "1px",
            "maxWidth": "1000px",
            "textAlign": "center"
        }
    )

    return Moment_resistance_section,html.Div([Mr_table])



app.layout = generate_layout()
if __name__ == "__main__":
    app.run(debug=True)
