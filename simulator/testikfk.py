import numpy as np

from geometry import (h, d, z0, B, p, beta)
from utils import R_x, R_y
from IK import pose, compute_alpha_rad, horn_endpoints


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
        P = (R @ p_top.T).T + T        
        v = P - H                      
        L = np.linalg.norm(v, axis=1)  
        r = L - d_len                 

        if np.linalg.norm(r) < tol_r:
            break

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
    pitch = np.arcsin(-R[2,0])
    roll  = np.arctan2(-R[1,2], R[1,1])
    return roll, pitch

def main():
    roll_d  = np.deg2rad(5.0)
    pitch_d = np.deg2rad(5.0)

    R_d, T_d, P_d = pose(roll_d, pitch_d, z0, p, R_y, R_x)  
    L = P_d - B                                             
    alpha_rad, flags = compute_alpha_rad(L, beta)           

    if any(flags):
        print("IK flags:", flags)

    H = horn_endpoints(B, beta, alpha_rad, h)

    R0 = np.eye(3)
    T0 = np.array([0.0, 0.0, z0])
    R_fk, T_fk = fk_iterative_6dof(H, p, d, R0, T0)

    roll_fk, pitch_fk = extract_roll_pitch_from_R(R_fk)

    print("IK servo angles (deg):", np.round(np.degrees(alpha_rad), 3))
    print("FK recovered roll/pitch (deg):",
          np.round(np.degrees(roll_fk), 3),
          np.round(np.degrees(pitch_fk), 3))

if __name__ == "__main__":
    main()