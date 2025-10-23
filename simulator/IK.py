# ik.py
import numpy as np
from geometry import h, d, IK_TOL

def pose(roll, pitch, z0, p_local, R_y, R_x):
    R = R_y(pitch) @ R_x(roll)
    T = np.array([0.0, 0.0, z0])
    P = (R @ p_local.T).T + T
    return R, T, P

def compute_alpha_rad(L, beta):
    """
    Solve e*sin(a) + f*cos(a) = g for each leg.
    Returns alpha (rad) and per-leg flags ('' if OK).
    """
    alpha = np.empty(6)
    flags = ['']*6
    for i in range(6):
        lx, ly, lz = L[i]
        e = 2*h*lz
        f = 2*h*(np.cos(beta[i])*lx + np.sin(beta[i])*ly)
        g = np.dot(L[i], L[i]) - (d*d - h*h)
        A = np.hypot(e, f)

        if A < IK_TOL:
            alpha[i] = np.nan; flags[i] = 'degenerate (A≈0)'; continue

        val = g / A
        if abs(val) > 1.0 + 1e-6:
            alpha[i] = np.nan; flags[i] = 'IK infeasible (|g|>A)'; continue

        val = np.clip(val, -1.0, 1.0)
        alpha[i] = np.arcsin(val) - np.arctan2(f, e)
    return alpha, flags

def horn_endpoints(B, beta, alpha, h_len):
    """
    Horn tip H = B + [h*cos(a)*cos(b), h*cos(a)*sin(b), h*sin(a)]
    """
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    H = np.column_stack([
        B[:,0] + h_len * (ca * cb),
        B[:,1] + h_len * (ca * sb),
        B[:,2] + h_len * (sa)
    ])
    return H