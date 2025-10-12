# test_fk_vs_ik.py
import numpy as np

from geometry import (h, d, z0, B, p, beta)
from utils import R_x, R_y
from IK import pose, compute_alpha_rad, horn_endpoints   # your IK module

# --- small helpers ---
def skew(v):
    x,y,z = v
    return np.array([[0,-z, y],[z,0,-x],[-y,x,0]])

def so3_exp(omega):
    th = np.linalg.norm(omega)
    if th < 1e-12:
        return np.eye(3) + skew(omega)
    k = omega / th
    K = skew(k)
    return np.eye(3) + np.sin(th)*K + (1-np.cos(th))*(K@K)

def fk_iterative_6dof(H, p_top, d_len, R_init, T_init, iters=20, lam=1e-3, tol_r=1e-6, tol_dx=1e-8):
    R = R_init.copy()
    T = T_init.copy()
    for _ in range(iters):
        P = (R @ p_top.T).T + T        # (6,3)
        v = P - H                      # (6,3)
        L = np.linalg.norm(v, axis=1)  # (6,)
        r = L - d_len                  # (6,)

        if np.linalg.norm(r) < tol_r:
            break

        # Jacobian (6x6): rows = [ u^T(-R[p]_x)   u^T ]
        J = np.zeros((6,6))
        for i in range(6):
            ui = v[i] / (L[i] + 1e-12)
            J[i, 0:3] = ui @ (-R @ skew(p_top[i]))
            J[i, 3:6] = ui

        A = J.T @ J + (lam**2)*np.eye(6)
        dx = -np.linalg.solve(A, J.T @ r)

        dR = so3_exp(dx[0:3])
        R = R @ dR
        T = T + dx[3:6]

        if np.linalg.norm(dx) < tol_dx:
            break
    return R, T

def extract_roll_pitch_from_R(R):
    # R = R_y(pitch) @ R_x(roll)  (your convention)
    pitch = np.arcsin(-R[2,0])
    roll  = np.arctan2(-R[1,2], R[1,1])
    return roll, pitch

def main():
    # --- Desired pose for the IK (roll=30°, pitch=30°, z=z0)
    roll_d  = np.deg2rad(5.0)
    pitch_d = np.deg2rad(5.0)

    R_d, T_d, P_d = pose(roll_d, pitch_d, z0, p, R_y, R_x)  # your helper: returns R, T, top pts
    L = P_d - B                                             # leg vectors base->top anchors
    alpha_rad, flags = compute_alpha_rad(L, beta)           # IK → servo angles

    if any(flags):
        print("IK flags:", flags)

    # Horn tips from IK angles (world/base frame)
    H = horn_endpoints(B, beta, alpha_rad, h)

    # --- Run FK from a simple seed (I, [0,0,z0])
    R0 = np.eye(3)
    T0 = np.array([0.0, 0.0, z0])
    R_fk, T_fk = fk_iterative_6dof(H, p, d, R0, T0)

    roll_fk, pitch_fk = extract_roll_pitch_from_R(R_fk)

    # --- Report
    print("IK servo angles (deg):", np.round(np.degrees(alpha_rad), 3))
    print("FK recovered roll/pitch (deg):",
          np.round(np.degrees(roll_fk), 3),
          np.round(np.degrees(pitch_fk), 3))

if __name__ == "__main__":
    main()