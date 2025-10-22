
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry import (h, d, r_top, r_bot, z0, B, p, beta, SERVO_MIN, SERVO_MAX)
from utils import R_x, R_y, circle_points

from IK import compute_alpha_rad  





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
    ALPHA_RATE_MAX_DEG = 360.0      
    TAU_SERVO = 0.04                
    DAMP_JAC = 1e-4                 

    
    fig = plt.figure(figsize=(11, 8))
    ax3d = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.14, right=0.78, bottom=0.26)

    
    ax3d.set_axis_on()
    ax3d.grid(True)

    
    XMIN, XMAX = -12, 12
    YMIN, YMAX = -12, 12
    ZMIN, ZMAX = 0, 25

    ax3d.set_xlim(XMIN, XMAX)
    ax3d.set_ylim(YMIN, YMAX)
    ax3d.set_zlim(ZMIN, ZMAX)
    ax3d.set_box_aspect([1,1,1])
    ax3d.view_init(elev=20, azim=-60)
    ax3d.set_xlabel('X (cm)')
    ax3d.set_ylabel('Y (cm)')
    ax3d.set_zlabel('Z (cm)')

    
    R_world = np.eye(3)

    
    world_axis_lines = []
    world_axis_lines += ax3d.plot([XMIN, XMAX], [0, 0], [0, 0], lw=1.2, linestyle='--', color='r', alpha=0.5)
    world_axis_lines += ax3d.plot([0, 0], [YMIN, YMAX], [0, 0], lw=1.2, linestyle='--', color='g', alpha=0.5)
    world_axis_lines += ax3d.plot([0, 0], [0, 0], [ZMIN, ZMAX], lw=1.2, linestyle='--', color='b', alpha=0.5)

    
    R_hat = np.eye(3)               
    T_hat = np.array([0.0, 0.0, z0])

    
    P0 = (R_hat @ p.T).T + T_hat
    alpha0, _ = compute_alpha_rad(P0 - B, beta)

    
    alpha_cmd  = alpha0.copy()
    alpha_meas = alpha0.copy()
    alpha_prev = alpha_meas.copy()

    
    base_circle_b = circle_points(r_bot)

    
    top_poly  = Poly3DCollection([np.zeros((50,3))], alpha=0.25, facecolor=(0,0,1,0.15), edgecolor='none')
    base_poly = Poly3DCollection([np.zeros((50,3))], alpha=0.18, facecolor=(0,0,0,0.10), edgecolor='none')
    ax3d.add_collection3d(top_poly); ax3d.add_collection3d(base_poly)

    world_frame = draw_frame(ax3d, R_world, [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    base_frame  = draw_frame(ax3d, np.eye(3), [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    top_frame   = draw_frame(ax3d, R_hat, T_hat, scale=3.0, linewidth=2.0, alpha=1.0)

    top_pts  = ax3d.scatter([], [], [], s=8, c='b')
    base_pts = ax3d.scatter([], [], [], s=8, c='k')
    horns = [ax3d.plot([], [], [], lw=2)[0] for _ in range(6)]
    rods  = [ax3d.plot([], [], [], lw=2, color='k')[0] for _ in range(6)]

    
    ax_bar = fig.add_axes([0.82, 0.25, 0.16, 0.65])
    yidx = np.arange(6)
    bars = ax_bar.barh(yidx, np.degrees(alpha_meas), align='center')
    ax_bar.set_yticks(yidx); ax_bar.set_yticklabels([f'k={i}' for i in range(6)])
    ax_bar.axvline(0, lw=0.8); ax_bar.set_xlim(-90, 90)
    ax_bar.set_xlabel('alpha [deg]'); ax_bar.set_title("Servo angles (measured)")

    
    status_text = fig.text(0.15, 0.015, "", color='red')
    imu_text = fig.text(0.15, 0.97, "", color='black')

    
    ax_wroll  = plt.axes([0.20, 0.12, 0.45, 0.03])
    ax_wpitch = plt.axes([0.20, 0.07, 0.45, 0.03])
    s_wroll   = Slider(ax_wroll,  'World roll [°]',  -45, 45, valinit=0)
    s_wpitch  = Slider(ax_wpitch, 'World pitch [°]', -45, 45, valinit=0)

    
    
    
    def update(_):
        nonlocal alpha_meas, alpha_cmd, alpha_prev, R_hat, T_hat

        
        wroll  = np.radians(s_wroll.val)
        wpitch = np.radians(s_wpitch.val)
        R_bw = R_y(wpitch) @ R_x(wroll)   
        R_wb = R_bw.T

        
        
        R_tb_des = R_bw.T
        T_tb_des = np.array([0.0, 0.0, z0])

        
        P_des_b = (R_tb_des @ p.T).T + T_tb_des
        alpha_des, flags = compute_alpha_rad(P_des_b - B, beta)
        ik_ok = angles_feasible(alpha_des) and (not any(flags))

        if not ik_ok:
            status_text.set_text("IK infeasible — holding last command")
            
        else:
            
            delta_cmd = alpha_des - alpha_cmd
            delta_cmd = (delta_cmd + np.pi) % (2*np.pi) - np.pi  
            rate_cap = np.radians(ALPHA_RATE_MAX_DEG) * DT
            delta_cmd = np.clip(delta_cmd, -rate_cap, +rate_cap)
            alpha_cmd = alpha_cmd + delta_cmd
            status_text.set_text("")

        
        alpha_meas = alpha_meas + (alpha_cmd - alpha_meas) * (1 - np.exp(-DT/TAU_SERVO))

        
        
        delta_alpha = (alpha_meas - alpha_prev)
        alpha_prev = alpha_meas.copy()

        
        J, alpha_nom = jacobian_alpha_wrt_rpy(R_hat, eps_rot=1e-3)

        
        
        JTJ = J.T @ J + (DAMP_JAC**2) * np.eye(3)
        dtheta = np.linalg.solve(JTJ, J.T @ delta_alpha)

        
        R_hat = R_hat @ so3_exp(dtheta)
        T_hat = np.array([0.0, 0.0, z0])  

        
        
        H_b = horn_endpoints_base(B, beta, alpha_meas, h)
        P_b = (R_hat @ p.T).T + T_hat

        
        B_w = transform_points_only_R(R_bw, B)
        H_w = transform_points_only_R(R_bw, H_b)
        P_w = transform_points_only_R(R_bw, P_b)
        top_circle_w  = transform_points_only_R(R_bw, (R_hat @ circle_points(r_top).T).T + T_hat)
        base_circle_w = transform_points_only_R(R_bw, base_circle_b)

        
        top_poly.set_verts([top_circle_w])
        base_poly.set_verts([base_circle_w])

        
        for ln in base_frame: ln.remove()
        for ln in top_frame:  ln.remove()
        base_frame[:] = draw_frame(ax3d, R_bw, [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)            
        top_frame[:]  = draw_frame(ax3d, R_bw @ R_hat, R_bw @ T_hat, scale=3.0, linewidth=2.0, alpha=1.0)  

        
        top_pts._offsets3d  = (P_w[:,0], P_w[:,1], P_w[:,2])
        base_pts._offsets3d = (B_w[:,0], B_w[:,1], B_w[:,2])
        for i in range(6):
            horns[i].set_data_3d([B_w[i,0], H_w[i,0]],
                                 [B_w[i,1], H_w[i,1]],
                                 [B_w[i,2], H_w[i,2]])
            rods[i].set_data_3d([H_w[i,0], P_w[i,0]],
                                [H_w[i,1], P_w[i,1]],
                                [H_w[i,2], P_w[i,2]])

        
        a_deg = np.degrees(alpha_meas)
        for i, bar in enumerate(bars):
            bar.set_width(a_deg[i])
            bad = (a_deg[i] < SERVO_MIN) or (a_deg[i] > SERVO_MAX) or (not np.isfinite(a_deg[i]))
            bar.set_color('red' if bad else 'C0')

        
        R_tw_est = R_bw @ R_hat
        
        pitch_now = np.arcsin(-R_tw_est[2,0])
        roll_now  = np.arctan2(-R_tw_est[1,2], R_tw_est[1,1])
        imu_text.set_text(f"IMU (world est): roll={np.degrees(roll_now):.1f}°, pitch={np.degrees(pitch_now):.1f}°  |  Gravity → -Z_world")

        fig.canvas.draw_idle()

    
    s_wroll.on_changed(update)
    s_wpitch.on_changed(update)

    plt.show()

if __name__ == "__main__":
    run()