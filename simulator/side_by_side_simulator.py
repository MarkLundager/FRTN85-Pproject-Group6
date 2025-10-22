
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry import (h, d, r_top, r_bot, z0, B, p, beta, SERVO_MIN, SERVO_MAX)
from utils import R_x, R_y, circle_points
from IK import compute_alpha_rad
from FK import fk_iterative_6dof, extract_roll_pitch_from_R





def horn_endpoints_base(B_base, beta, alpha, h_len):
    
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    v = np.column_stack([h_len*(ca*cb), h_len*(ca*sb), h_len*sa])
    return B_base + v

def draw_frame(ax, R, T, scale=3.0, linewidth=1.8, alpha=1.0):
    origin = np.array(T)
    axes = np.eye(3) * scale
    cols = ['r','g','b']
    lines = []
    for i in range(3):
        vec = origin + R @ axes[:, i]
        line, = ax.plot([origin[0], vec[0]],
                        [origin[1], vec[1]],
                        [origin[2], vec[2]],
                        color=cols[i], lw=linewidth, alpha=alpha)
        lines.append(line)
    return lines

def angles_feasible(alpha):
    if not np.all(np.isfinite(alpha)):
        return False
    deg = np.degrees(alpha)
    return np.all((deg >= SERVO_MIN) & (deg <= SERVO_MAX))

def transform_points_only_R(R, P):
    return (R @ P.T).T

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


def jacobian_alpha_wrt_rpy(R_hat, eps_rot=1e-3):
    """
    Returns J (6x3) where columns are d(alpha)/d(rx), d(alpha)/d(ry), d(alpha)/d(rz)
    computed numerically via IK at pose (R_hat, z0).
    Also returns alpha_nom at R_hat.
    """
    P_nom = (R_hat @ p.T).T + np.array([0,0,z0])
    alpha_nom, _ = compute_alpha_rad(P_nom - B, beta)

    J = np.zeros((6,3))
    e = np.eye(3)
    for j in range(3):
        dR = so3_exp(e[:,j] * eps_rot)
        R_pert = R_hat @ dR
        P_pert = (R_pert @ p.T).T + np.array([0,0,z0])
        alpha_pert, _ = compute_alpha_rad(P_pert - B, beta)
        J[:, j] = (alpha_pert - alpha_nom) / eps_rot
    return J, alpha_nom





def run():
    
    DT = 0.02                       
    TAU_SERVO = 0.04                
    
    ALPHA_RATE_MAX_DEG_FK  = 360.0  
    
    ALPHA_RATE_MAX_DEG_INT = 600.0  

    
    ANGLE_NOISE_STD_DEG = 0.20      
    BIAS_RW_STD_DEG     = 0.01      
    JAC_EPS_ROT         = 2e-3      
    DAMP_JAC            = 5e-5      
    GEOM_SCALE_H        = 1.02      

    rng = np.random.default_rng()

    
    fig = plt.figure(figsize=(16, 8))
    ax_fk  = fig.add_subplot(121, projection='3d')
    ax_int = fig.add_subplot(122, projection='3d')
    plt.subplots_adjust(left=0.08, right=0.95, bottom=0.23, wspace=0.1)

    
    for ax in (ax_fk, ax_int):
        ax.set_axis_on(); ax.grid(True)
        ax.set_xlim(-12, 12); ax.set_ylim(-12, 12); ax.set_zlim(0, 25)
        ax.set_box_aspect([1,1,1]); ax.view_init(elev=20, azim=-60)
        ax.set_xlabel('X (cm)'); ax.set_ylabel('Y (cm)'); ax.set_zlabel('Z (cm)')

    ax_fk.set_title('FK-based (closure each tick)')
    ax_int.set_title('Jacobian integration (no FK) — medium drift')

    
    for ax in (ax_fk, ax_int):
        ax.plot([-12, 12], [0, 0], [0, 0], lw=1.2, ls='--', color='r', alpha=0.5)
        ax.plot([0, 0], [-12, 12], [0, 0], lw=1.2, ls='--', color='g', alpha=0.5)
        ax.plot([0, 0], [0, 0], [0, 25],     lw=1.2, ls='--', color='b', alpha=0.5)

    
    R_tb0 = np.eye(3)
    T_tb0 = np.array([0.0, 0.0, z0])

    
    P0 = (R_tb0 @ p.T).T + T_tb0
    alpha0, _ = compute_alpha_rad(P0 - B, beta)

    
    R_tb_fk = R_tb0.copy()
    T_tb_fk = T_tb0.copy()
    alpha_cmd_fk  = alpha0.copy()
    alpha_meas_fk = alpha0.copy()
    R_seed_fk = R_tb_fk.copy()
    T_seed_fk = T_tb_fk.copy()
    last_ok_fk = {
        "R_tb": R_tb_fk.copy(),
        "T_tb": T_tb_fk.copy(),
        "alpha_meas": alpha_meas_fk.copy(),
        "alpha_cmd": alpha_cmd_fk.copy()
    }

    
    R_hat_int = R_tb0.copy()        
    T_hat_int = T_tb0.copy()
    alpha_cmd_int  = alpha0.copy()
    alpha_meas_int = alpha0.copy()
    alpha_prev_int = alpha_meas_int.copy()
    bias_int = np.zeros(6)          

    
    base_circle_b = circle_points(r_bot)

    
    top_poly_fk  = Poly3DCollection([np.zeros((50,3))], alpha=0.25, facecolor=(0,0,1,0.15), edgecolor='none')
    base_poly_fk = Poly3DCollection([np.zeros((50,3))], alpha=0.18, facecolor=(0,0,0,0.10), edgecolor='none')
    ax_fk.add_collection3d(top_poly_fk); ax_fk.add_collection3d(base_poly_fk)
    world_frame_fk = draw_frame(ax_fk, np.eye(3), [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    base_frame_fk  = draw_frame(ax_fk, np.eye(3), [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    top_frame_fk   = draw_frame(ax_fk, R_tb_fk,  T_tb_fk,  scale=3.0, linewidth=2.0, alpha=1.0)
    top_pts_fk  = ax_fk.scatter([], [], [], s=8, c='b')
    base_pts_fk = ax_fk.scatter([], [], [], s=8, c='k')
    horns_fk = [ax_fk.plot([], [], [], lw=2)[0] for _ in range(6)]
    rods_fk  = [ax_fk.plot([], [], [], lw=2, color='k')[0] for _ in range(6)]
    imu_text_fk = fig.text(0.10, 0.96, "", color='black')

    
    top_poly_int  = Poly3DCollection([np.zeros((50,3))], alpha=0.25, facecolor=(0,0,1,0.15), edgecolor='none')
    base_poly_int = Poly3DCollection([np.zeros((50,3))], alpha=0.18, facecolor=(0,0,0,0.10), edgecolor='none')
    ax_int.add_collection3d(top_poly_int); ax_int.add_collection3d(base_poly_int)
    world_frame_int = draw_frame(ax_int, np.eye(3), [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    base_frame_int  = draw_frame(ax_int, np.eye(3), [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    top_frame_int   = draw_frame(ax_int, R_hat_int, T_hat_int, scale=3.0, linewidth=2.0, alpha=1.0)
    top_pts_int  = ax_int.scatter([], [], [], s=8, c='b')
    base_pts_int = ax_int.scatter([], [], [], s=8, c='k')
    horns_int = [ax_int.plot([], [], [], lw=2)[0] for _ in range(6)]
    rods_int  = [ax_int.plot([], [], [], lw=2, color='k')[0] for _ in range(6)]
    imu_text_int = fig.text(0.56, 0.96, "", color='black')

    
    status_text = fig.text(0.08, 0.02, "", color='red')

    
    ax_wroll  = plt.axes([0.18, 0.10, 0.60, 0.03])
    ax_wpitch = plt.axes([0.18, 0.06, 0.60, 0.03])
    s_wroll   = Slider(ax_wroll,  'World roll [°]',  -45, 45, valinit=0)
    s_wpitch  = Slider(ax_wpitch, 'World pitch [°]', -45, 45, valinit=0)

    
    
    
    def update(_):
        nonlocal R_tb_fk, T_tb_fk, R_seed_fk, T_seed_fk, alpha_cmd_fk, alpha_meas_fk, last_ok_fk
        nonlocal R_hat_int, T_hat_int, alpha_cmd_int, alpha_meas_int, alpha_prev_int, bias_int

        
        wroll  = np.radians(s_wroll.val)
        wpitch = np.radians(s_wpitch.val)
        R_bw = R_y(wpitch) @ R_x(wroll)   

        
        R_tb_des = R_bw.T
        T_tb_des = np.array([0.0, 0.0, z0])

        
        P_des_b = (R_tb_des @ p.T).T + T_tb_des
        alpha_des, flags = compute_alpha_rad(P_des_b - B, beta)
        ik_ok = angles_feasible(alpha_des) and (not any(flags))

        if not ik_ok:
            status_text.set_text("IK infeasible — holding commands")
        else:
            status_text.set_text("")

        
        rate_cap_fk = np.radians(ALPHA_RATE_MAX_DEG_FK) * DT
        if ik_ok:
            delta_fk = alpha_des - alpha_cmd_fk
            delta_fk = (delta_fk + np.pi) % (2*np.pi) - np.pi
            delta_fk = np.clip(delta_fk, -rate_cap_fk, +rate_cap_fk)
            alpha_cmd_fk = alpha_cmd_fk + delta_fk
        alpha_meas_fk = alpha_meas_fk + (alpha_cmd_fk - alpha_meas_fk) * (1 - np.exp(-DT/TAU_SERVO))

        
        H_b_fk = horn_endpoints_base(B, beta, alpha_meas_fk, h)
        R_try_b_fk, T_try_b_fk = fk_iterative_6dof(H_b_fk, p, d, R_seed_fk, T_seed_fk)

        
        P_try_b_fk = (R_try_b_fk @ p.T).T + T_try_b_fk
        L_fk = np.linalg.norm(P_try_b_fk - H_b_fk, axis=1)
        fk_ok = np.all(np.isfinite(L_fk)) and np.all(np.abs(L_fk - d) <= 1e-3)

        if fk_ok and angles_feasible(alpha_meas_fk):
            R_tb_fk, T_tb_fk = R_try_b_fk, T_try_b_fk
            R_seed_fk, T_seed_fk = R_tb_fk.copy(), T_tb_fk.copy()
            last_ok_fk = {
                "R_tb": R_tb_fk.copy(),
                "T_tb": T_tb_fk.copy(),
                "alpha_meas": alpha_meas_fk.copy(),
                "alpha_cmd": alpha_cmd_fk.copy()
            }
        else:
            R_tb_fk, T_tb_fk = last_ok_fk["R_tb"], last_ok_fk["T_tb"]
            alpha_meas_fk = last_ok_fk["alpha_meas"]
            alpha_cmd_fk  = last_ok_fk["alpha_cmd"]
            R_seed_fk, T_seed_fk = R_tb_fk.copy(), T_tb_fk.copy()
            status_text.set_text("FK infeasible — clamped (left)")

        
        rate_cap_int = np.radians(ALPHA_RATE_MAX_DEG_INT) * DT
        if ik_ok:
            delta_int = alpha_des - alpha_cmd_int
            delta_int = (delta_int + np.pi) % (2*np.pi) - np.pi
            delta_int = np.clip(delta_int, -rate_cap_int, +rate_cap_int)
            alpha_cmd_int = alpha_cmd_int + delta_int

        
        alpha_meas_int = alpha_meas_int + (alpha_cmd_int - alpha_meas_int) * (1 - np.exp(-DT/TAU_SERVO))

        
        bias_int += np.deg2rad(BIAS_RW_STD_DEG) * rng.standard_normal(6)
        alpha_meas_int += bias_int
        alpha_meas_int += np.deg2rad(ANGLE_NOISE_STD_DEG) * rng.standard_normal(6)

        
        delta_alpha = (alpha_meas_int - alpha_prev_int)
        alpha_prev_int = alpha_meas_int.copy()

        
        delta_alpha *= GEOM_SCALE_H

        J, _ = jacobian_alpha_wrt_rpy(R_hat_int, eps_rot=JAC_EPS_ROT)
        JTJ = J.T @ J + (DAMP_JAC**2) * np.eye(3)
        dtheta = np.linalg.solve(JTJ, J.T @ delta_alpha)
        R_hat_int = R_hat_int @ so3_exp(dtheta)
        T_hat_int = np.array([0.0, 0.0, z0])

        
        B_w_fk = transform_points_only_R(R_bw, B)
        H_w_fk = transform_points_only_R(R_bw, H_b_fk)
        P_w_fk = transform_points_only_R(R_bw, (R_tb_fk @ p.T).T + T_tb_fk)
        top_circle_w_fk  = transform_points_only_R(R_bw, (R_tb_fk @ circle_points(r_top).T).T + T_tb_fk)
        base_circle_w_fk = transform_points_only_R(R_bw, base_circle_b)

        top_poly_fk.set_verts([top_circle_w_fk])
        base_poly_fk.set_verts([base_circle_w_fk])

        for ln in base_frame_fk: ln.remove()
        for ln in top_frame_fk:  ln.remove()
        base_frame_fk[:] = draw_frame(ax_fk, R_bw, [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
        top_frame_fk[:]  = draw_frame(ax_fk, R_bw @ R_tb_fk, R_bw @ T_tb_fk, scale=3.0, linewidth=2.0, alpha=1.0)

        top_pts_fk._offsets3d  = (P_w_fk[:,0], P_w_fk[:,1], P_w_fk[:,2])
        base_pts_fk._offsets3d = (B_w_fk[:,0], B_w_fk[:,1], B_w_fk[:,2])
        for i in range(6):
            horns_fk[i].set_data_3d([B_w_fk[i,0], H_w_fk[i,0]],
                                    [B_w_fk[i,1], H_w_fk[i,1]],
                                    [B_w_fk[i,2], H_w_fk[i,2]])
            rods_fk[i].set_data_3d([H_w_fk[i,0], P_w_fk[i,0]],
                                   [H_w_fk[i,1], P_w_fk[i,1]],
                                   [H_w_fk[i,2], P_w_fk[i,2]])

        R_tw_fk = R_bw @ R_tb_fk
        roll_fk, pitch_fk = extract_roll_pitch_from_R(R_tw_fk)
        imu_text_fk.set_text(f"FK IMU: roll={np.degrees(roll_fk):.1f}°, pitch={np.degrees(pitch_fk):.1f}°")

        
        H_b_int = horn_endpoints_base(B, beta, alpha_meas_int, h)
        B_w_int = transform_points_only_R(R_bw, B)
        H_w_int = transform_points_only_R(R_bw, H_b_int)
        P_w_int = transform_points_only_R(R_bw, (R_hat_int @ p.T).T + T_hat_int)
        top_circle_w_int  = transform_points_only_R(R_bw, (R_hat_int @ circle_points(r_top).T).T + T_hat_int)
        base_circle_w_int = transform_points_only_R(R_bw, base_circle_b)

        top_poly_int.set_verts([top_circle_w_int])
        base_poly_int.set_verts([base_circle_w_int])

        for ln in base_frame_int: ln.remove()
        for ln in top_frame_int:  ln.remove()
        base_frame_int[:] = draw_frame(ax_int, R_bw, [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
        top_frame_int[:]  = draw_frame(ax_int, R_bw @ R_hat_int, R_bw @ T_hat_int, scale=3.0, linewidth=2.0, alpha=1.0)

        top_pts_int._offsets3d  = (P_w_int[:,0], P_w_int[:,1], P_w_int[:,2])
        base_pts_int._offsets3d = (B_w_int[:,0], B_w_int[:,1], B_w_int[:,2])
        for i in range(6):
            horns_int[i].set_data_3d([B_w_int[i,0], H_w_int[i,0]],
                                     [B_w_int[i,1], H_w_int[i,1]],
                                     [B_w_int[i,2], H_w_int[i,2]])
            rods_int[i].set_data_3d([H_w_int[i,0], P_w_int[i,0]],
                                    [H_w_int[i,1], P_w_int[i,1]],
                                    [H_w_int[i,2], P_w_int[i,2]])

        R_tw_int = R_bw @ R_hat_int
        pitch_now_int = np.arcsin(-R_tw_int[2,0])
        roll_now_int  = np.arctan2(-R_tw_int[1,2], R_tw_int[1,1])
        imu_text_int.set_text(f"INT IMU: roll={np.degrees(roll_now_int):.1f}°, pitch={np.degrees(pitch_now_int):.1f}°")

        fig.canvas.draw_idle()

    
    s_wroll.on_changed(update)
    s_wpitch.on_changed(update)

    plt.show()

if __name__ == "__main__":
    run()