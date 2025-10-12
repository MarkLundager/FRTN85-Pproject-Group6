import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry import (h, d, r_top, r_bot, z0, B, p, beta, SERVO_MIN, SERVO_MAX)
from utils import R_x, R_y, circle_points
from IK import compute_alpha_rad
from FK import fk_iterative_6dof

def horn_endpoints_base(B_base, beta, alpha, h_len):
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    v = np.column_stack([h_len*(ca*cb), h_len*(ca*sb), h_len*sa])
    return B_base + v

def draw_frame(ax, R, T, scale=3.0, linewidth=1.5):
    origin = np.array(T); axes = np.eye(3)*scale; cols = ['r','g','b']; lines=[]
    for i in range(3):
        vec = origin + R @ axes[:, i]
        line, = ax.plot([origin[0], vec[0]],[origin[1], vec[1]],[origin[2], vec[2]],
                        color=cols[i], lw=linewidth)
        lines.append(line)
    return lines

def transform_points_only_R(R, P_dev):
    return (R @ P_dev.T).T

def angles_feasible(alpha):
    if not np.all(np.isfinite(alpha)): return False
    deg = np.degrees(alpha)
    return np.all((deg >= SERVO_MIN) & (deg <= SERVO_MAX))

def rods_feasible(R_tb, T_tb, alpha, B_dev, p_top, d_len, tol=1e-3):
    """Check |HP|≈d for all legs using current pose and servo angles."""
    H = horn_endpoints_base(B_dev, beta, alpha, h)
    P = (R_tb @ p_top.T).T + T_tb
    L = np.linalg.norm(P - H, axis=1)
    return np.all(np.isfinite(L)) and np.all(np.abs(L - d_len) <= tol)

def run():
    fig = plt.figure(figsize=(11, 8))
    ax3d = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.14, right=0.78, bottom=0.22)
    fig.patch.set_facecolor('white'); ax3d.set_facecolor('white'); ax3d.set_axis_off()

    B_dev = B.copy()

    R_tb = np.eye(3)
    T_tb = np.array([0.0, 0.0, z0])
    P_dev0 = (R_tb @ p.T).T + T_tb
    alpha_act, _ = compute_alpha_rad(P_dev0 - B_dev, beta)
    alpha_cmd = alpha_act.copy()

    last_ok = {
        "R_tb": R_tb.copy(),
        "T_tb": T_tb.copy(),
        "alpha_act": alpha_act.copy(),
        "alpha_cmd": alpha_cmd.copy()
    }

    H_dev = horn_endpoints_base(B_dev, beta, alpha_act, h)
    top_circle_dev  = (R_tb @ circle_points(r_top).T).T + T_tb
    base_circle_dev = circle_points(r_bot)

    top_poly  = Poly3DCollection([top_circle_dev],  alpha=0.25, facecolor=(0,0,1,0.15), edgecolor='none')
    base_poly = Poly3DCollection([base_circle_dev], alpha=0.18, facecolor=(0,0,0,0.10), edgecolor='none')
    ax3d.add_collection3d(top_poly); ax3d.add_collection3d(base_poly)

    base_frame = draw_frame(ax3d, np.eye(3), [0,0,0], scale=3.0)
    top_frame  = draw_frame(ax3d, R_tb, T_tb, scale=3.0)

    P_dev = (R_tb @ p.T).T + T_tb
    top_pts  = ax3d.scatter(P_dev[:,0], P_dev[:,1], P_dev[:,2], s=8, c='b')
    base_pts = ax3d.scatter(B_dev[:,0], B_dev[:,1], B_dev[:,2], s=8, c='k')
    horns = [ax3d.plot([B_dev[i,0], H_dev[i,0]],
                       [B_dev[i,1], H_dev[i,1]],
                       [B_dev[i,2], H_dev[i,2]], lw=2)[0] for i in range(6)]
    rods  = [ax3d.plot([H_dev[i,0], P_dev[i,0]],
                       [H_dev[i,1], P_dev[i,1]],
                       [H_dev[i,2], P_dev[i,2]], lw=2, color='k')[0] for i in range(6)]


    ax3d.set_xlim(-12, 12); ax3d.set_ylim(-12, 12); ax3d.set_zlim(0, 25)
    ax3d.set_box_aspect([1,1,1]); ax3d.view_init(elev=20, azim=-60)


    ax_bar = fig.add_axes([0.82, 0.25, 0.16, 0.65])
    yidx = np.arange(6)
    bars = ax_bar.barh(yidx, np.degrees(alpha_act), align='center')
    ax_bar.set_yticks(yidx); ax_bar.set_yticklabels([f'k={i}' for i in range(6)])
    ax_bar.axvline(0, lw=0.8); ax_bar.set_xlim(-90, 90)
    ax_bar.set_xlabel('alpha [deg]'); ax_bar.set_title("Servo angles")


    status_text = fig.text(0.15, 0.015, "", color='red')


    ax_broll  = plt.axes([0.20, 0.09, 0.45, 0.03])
    ax_bpitch = plt.axes([0.20, 0.04, 0.45, 0.03])
    s_broll   = Slider(ax_broll,  'Global roll [°]',  -45, 45, valinit=0)
    s_bpitch  = Slider(ax_bpitch, 'Global pitch [°]', -45, 45, valinit=0)


    DT  = 0.02
    TAU = 0.06

    R_seed = R_tb.copy()
    T_seed = T_tb.copy()

    def update(_):
        nonlocal alpha_act, alpha_cmd, R_tb, T_tb, R_seed, T_seed, last_ok


        broll  = np.radians(s_broll.val)
        bpitch = np.radians(s_bpitch.val)
        R_bw   = R_y(bpitch) @ R_x(broll)

        R_tb_des = R_bw.T
        T_tb_des = np.array([0.0, 0.0, z0])

        P_dev_des = (R_tb_des @ p.T).T + T_tb_des
        alpha_des, flags = compute_alpha_rad(P_dev_des - B_dev, beta)

        ik_ok = angles_feasible(alpha_des) and (not any(flags))
        if ik_ok:
            alpha_cmd = alpha_des
        else:
            status_text.set_text("IK infeasible — clamping to last feasible command")
            alpha_cmd = last_ok["alpha_cmd"]

        alpha_next = alpha_act + (alpha_cmd - alpha_act) * (1 - np.exp(-DT/TAU))

        H_try = horn_endpoints_base(B_dev, beta, alpha_next, h)
        R_try, T_try = fk_iterative_6dof(H_try, p, d, R_seed, T_seed)

        fk_ok = rods_feasible(R_try, T_try, alpha_next, B_dev, p, d)
        if fk_ok and angles_feasible(alpha_next):
            R_tb, T_tb = R_try, T_try
            alpha_act  = alpha_next
            R_seed, T_seed = R_tb.copy(), T_tb.copy()
            last_ok = {
                "R_tb": R_tb.copy(),
                "T_tb": T_tb.copy(),
                "alpha_act": alpha_act.copy(),
                "alpha_cmd": alpha_cmd.copy()
            }
            status_text.set_text("")
        else:
            R_tb, T_tb   = last_ok["R_tb"], last_ok["T_tb"]
            alpha_act    = last_ok["alpha_act"]
            alpha_cmd    = last_ok["alpha_cmd"]
            R_seed, T_seed = R_tb.copy(), T_tb.copy()
            if not ik_ok:
                status_text.set_text("IK infeasible — clamped")
            else:
                status_text.set_text("FK/rod constraint infeasible — clamped")

        H_dev = horn_endpoints_base(B_dev, beta, alpha_act, h)
        P_dev = (R_tb @ p.T).T + T_tb

        B_w = transform_points_only_R(R_bw, B_dev)
        H_w = transform_points_only_R(R_bw, H_dev)
        P_w = transform_points_only_R(R_bw, P_dev)

        top_circle_w  = transform_points_only_R(R_bw, (R_tb @ circle_points(r_top).T).T + T_tb)
        base_circle_w = transform_points_only_R(R_bw, base_circle_dev)

        top_poly.set_verts([top_circle_w])
        base_poly.set_verts([base_circle_w])

        for ln in base_frame: ln.remove()
        for ln in top_frame:  ln.remove()
        base_frame[:] = draw_frame(ax3d, R_bw, [0,0,0], scale=3.0, linewidth=1.5)
        top_frame[:]  = draw_frame(ax3d, R_bw @ R_tb, R_bw @ T_tb, scale=3.0, linewidth=1.5)

        top_pts._offsets3d  = (P_w[:,0], P_w[:,1], P_w[:,2])
        base_pts._offsets3d = (B_w[:,0], B_w[:,1], B_w[:,2])
        for i in range(6):
            horns[i].set_data_3d([B_w[i,0], H_w[i,0]],[B_w[i,1], H_w[i,1]],[B_w[i,2], H_w[i,2]])
            rods[i].set_data_3d([H_w[i,0], P_w[i,0]],[H_w[i,1], P_w[i,1]],[H_w[i,2], P_w[i,2]])

        a_deg = np.degrees(alpha_act)
        for i, bar in enumerate(bars):
            bar.set_width(a_deg[i])
            bad = (a_deg[i] < SERVO_MIN) or (a_deg[i] > SERVO_MAX) or (not np.isfinite(a_deg[i]))
            bar.set_color('red' if bad else 'C0')

        fig.canvas.draw_idle()

    s_broll.on_changed(update)
    s_bpitch.on_changed(update)
    plt.show()

if __name__ == "__main__":
    run()