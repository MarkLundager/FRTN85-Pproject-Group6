import numpy as np

# ---------- FK helpers ----------
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

def fk_iterative_6dof(H, p_top, d_len, R_init, T_init,
                      iters=20, lam=1e-3, tol_r=1e-6, tol_dx=1e-8):
    """
    Solve for pose (R,T) s.t. ||R p_i + T - H_i|| = d_len, i=1..6.
      H: (6,3) horn tips in world
      p_top: (6,3) top anchors in top frame
      d_len: scalar rod length
      R_init: (3,3) initial rotation
      T_init: (3,)  initial translation
    """
    R = R_init.copy()
    T = T_init.copy()
    for _ in range(iters):
        P = (R @ p_top.T).T + T          # (6,3)
        v = P - H                        # (6,3)
        L = np.linalg.norm(v, axis=1)    # (6,)
        r = L - d_len                    # (6,)

        if np.linalg.norm(r) < tol_r:
            break

        # Jacobian (6x6): row i = [ u_i^T (-R [p_i]_x) , u_i^T ]
        J = np.zeros((6,6))
        for i in range(6):
            ui = v[i] / (L[i] + 1e-12)
            J[i, 0:3] = ui @ (-R @ skew(p_top[i]))
            J[i, 3:6] = ui

        A = J.T @ J + (lam**2) * np.eye(6)
        dx = -np.linalg.solve(A, J.T @ r)

        dR = so3_exp(dx[0:3])
        R = R @ dR
        T = T + dx[3:6]

        if np.linalg.norm(dx) < tol_dx:
            break
    return R, T

def extract_roll_pitch_from_R(R):
    # For convention R = R_y(pitch) @ R_x(roll)
    pitch = np.arcsin(-R[2,0])
    roll  = np.arctan2(-R[1,2], R[1,1])
    return roll, pitch