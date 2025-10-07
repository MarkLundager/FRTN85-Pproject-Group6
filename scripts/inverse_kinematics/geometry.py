# geometry.py
import numpy as np

# --- Core dimensions (cm) ---
h = 1.64          # servo horn length
d = 17.0          # rod length (fixed)
r_top = 7.8       # top platform radius
r_bot = 9.3       # base radius
z0 = 16.5         # nominal height

# --- Pair separations (radians) ---
pDelta = np.deg2rad(15.25)   # platform pair separation
bDelta = np.deg2rad(23.58)   # base pair separation

# --- Build RAW layout angles ---
k = np.arange(6)
p_phi = (2*np.pi/3)*np.floor(k/2) - ((-1)**k)*(pDelta/2) + np.pi/3
b_phi = (2*np.pi/3)*np.floor((k+1)/2) + ((-1)**k)*(bDelta/2)
beta  = b_phi + (np.pi/2)*((-1)**k)   # servo tangent orientation in XY

# --- Anchors in local frames ---
B = np.column_stack([r_bot*np.cos(b_phi), r_bot*np.sin(b_phi), np.zeros(6)])
p = np.column_stack([r_top*np.cos(p_phi), r_top*np.sin(p_phi), np.zeros(6)])

# --- Limits and tolerances ---
SERVO_MIN = -90.0   # deg
SERVO_MAX =  90.0   # deg
IK_TOL  = 1e-9
ROD_TOL = 1e-3      # cm tolerance on |HP| == d
