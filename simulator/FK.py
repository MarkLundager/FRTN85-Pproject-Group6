import numpy as np

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
            ui = v[i] / (L[i] + 1e-12)  # compute direction vector (normalized) which shows where we're going translation wise.
            J[i, 0:3] = ui @ (-R @ skew(p_top[i])) #Multiplying by R brings the point into the base frame and expresses movement based on pithc, roll and jaw. multiply by ui to see the movemnet in the leg length.
            J[i, 3:6] = ui 

        A = J.T @ J + (lam**2) * np.eye(6)  # we wish to solve r + J * delta x = 0, i.e. J*delta x = -r However, we use least squares problem to handle ill condiiton or non square matrices. which then becomes A = (JT*J+Y^2I)dX = -JT*r
        dx = -np.linalg.solve(A, J.T @ r)  # we solve it 

        dR = so3_exp(dx[0:3]) #we have rotation as [wx, wy, wz] -> transform into rotation matrix 3x3
        R = R @ dR  #correct
        T = T + dx[3:6] # correct

        if np.linalg.norm(dx) < tol_dx:
            break
    return R, T

def extract_roll_pitch_from_R(R):
    pitch = np.arcsin(-R[2,0])
    roll  = np.arctan2(-R[1,2], R[1,1])
    return roll, pitch