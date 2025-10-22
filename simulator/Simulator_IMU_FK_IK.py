
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry import (h, d, r_top, r_bot, z0, B, p, beta, SERVO_MIN, SERVO_MAX)
from utils import R_x, R_y, circle_points
from IK import compute_alpha_rad
from FK import fk_iterative_6dof, extract_roll_pitch_from_R

# ---------------------------
# Helpers
# ---------------------------

def horn_endpoints_base(B_base, beta, alpha, h_len):
    """Horn tip H = B + [h*cos(a)*cos(b), h*cos(a)*sin(b), h*sin(a)]"""
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

# ---------------------------
# Main
# ---------------------------

def run():
    # --- Sim/servo params ---
    DT = 0.02                       # s, UI tick
    ALPHA_RATE_MAX_DEG = 360.0      # deg/s, command slew cap (~servo max speed)
    TAU_SERVO = 0.04                # s, internal servo following lag (sim)

    # --- Figure / axes ---
    fig = plt.figure(figsize=(11, 8))
    ax3d = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.14, right=0.78, bottom=0.26)

    # Show full coordinate axes with ticks & grid
    ax3d.set_axis_on()
    ax3d.grid(True)

    # Scene bounds
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

    # Global/world frame (fixed on screen). Gravity points along -Z_world.
    R_world = np.eye(3)

    # Draw long world axes lines through the whole window
    # X axis line (y=0, z=0), Y axis line (x=0, z=0), Z axis line (x=0, y=0)
    world_axis_lines = []
    world_axis_lines += ax3d.plot([XMIN, XMAX], [0, 0], [0, 0], lw=1.2, linestyle='--', color='r', alpha=0.5)
    world_axis_lines += ax3d.plot([0, 0], [YMIN, YMAX], [0, 0], lw=1.2, linestyle='--', color='g', alpha=0.5)
    world_axis_lines += ax3d.plot([0, 0], [0, 0], [ZMIN, ZMAX], lw=1.2, linestyle='--', color='b', alpha=0.5)

    # --- Initial top pose (relative to base) ---
    R_tb = np.eye(3)
    T_tb = np.array([0.0, 0.0, z0])

    # --- Initial IK to get starting angles (base frame) ---
    P0 = (R_tb @ p.T).T + T_tb
    alpha0, _ = compute_alpha_rad(P0 - B, beta)

    # Commanded and measured (encoder) angles
    alpha_cmd  = alpha0.copy()
    alpha_meas = alpha0.copy()

    # Seeds for FK (base frame)
    R_seed = R_tb.copy()
    T_seed = T_tb.copy()

    # --- Last feasible state for clamping ---
    last_ok = {
        "R_tb": R_tb.copy(),
        "T_tb": T_tb.copy(),
        "alpha_meas": alpha_meas.copy(),
        "alpha_cmd": alpha_cmd.copy()
    }

    # --- Static base geometry in base frame ---
    base_circle_b = circle_points(r_bot)

    # --- Initial drawing (we’ll transform to world each tick) ---
    # Create placeholders; we’ll set their data each update
    top_poly  = Poly3DCollection([np.zeros((50,3))], alpha=0.25, facecolor=(0,0,1,0.15), edgecolor='none')
    base_poly = Poly3DCollection([np.zeros((50,3))], alpha=0.18, facecolor=(0,0,0,0.10), edgecolor='none')
    ax3d.add_collection3d(top_poly); ax3d.add_collection3d(base_poly)

    # Frames
    world_frame = draw_frame(ax3d, R_world, [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    base_frame  = draw_frame(ax3d, np.eye(3), [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)
    top_frame   = draw_frame(ax3d, R_tb, T_tb, scale=3.0, linewidth=2.0, alpha=1.0)

    # Points & segments (placeholders)
    top_pts  = ax3d.scatter([], [], [], s=8, c='b')
    base_pts = ax3d.scatter([], [], [], s=8, c='k')
    horns = [ax3d.plot([], [], [], lw=2)[0] for _ in range(6)]
    rods  = [ax3d.plot([], [], [], lw=2, color='k')[0] for _ in range(6)]

    # --- Servo bars (measured) ---
    ax_bar = fig.add_axes([0.82, 0.25, 0.16, 0.65])
    yidx = np.arange(6)
    bars = ax_bar.barh(yidx, np.degrees(alpha_meas), align='center')
    ax_bar.set_yticks(yidx); ax_bar.set_yticklabels([f'k={i}' for i in range(6)])
    ax_bar.axvline(0, lw=0.8); ax_bar.set_xlim(-90, 90)
    ax_bar.set_xlabel('alpha [deg]'); ax_bar.set_title("Servo angles (measured)")

    # --- Status + IMU text ---
    status_text = fig.text(0.15, 0.015, "", color='red')
    imu_text = fig.text(0.15, 0.97, "", color='black')

    # --- Sliders: WORLD roll/pitch (tilt of base/ground relative to screen/world) ---
    ax_wroll  = plt.axes([0.20, 0.12, 0.45, 0.03])
    ax_wpitch = plt.axes([0.20, 0.07, 0.45, 0.03])
    s_wroll   = Slider(ax_wroll,  'World roll [°]',  -45, 45, valinit=0)
    s_wpitch  = Slider(ax_wpitch, 'World pitch [°]', -45, 45, valinit=0)

    # ---------------------------
    # Update loop
    # ---------------------------
    def update(_):
        nonlocal alpha_meas, alpha_cmd, R_tb, T_tb, R_seed, T_seed, last_ok

        # 1) Read world tilt (base->world rotation)
        wroll  = np.radians(s_wroll.val)
        wpitch = np.radians(s_wpitch.val)
        R_bw = R_y(wpitch) @ R_x(wroll)   # base -> world
        R_wb = R_bw.T                     # world -> base

        # 2) Controller goal: keep TOP LEVEL IN WORLD (parallel to XY_world).
        #    That means desired TOP ORIENTATION IN WORLD is Identity.
        #    Convert to a desired orientation in BASE frame for IK:
        #       R_bw * R_tb_des = I  => R_tb_des = R_bw^T
        R_tb_des = R_bw.T
        T_tb_des = np.array([0.0, 0.0, z0])

        # 3) IK in BASE frame to get desired joint angles
        P_des_b = (R_tb_des @ p.T).T + T_tb_des
        alpha_des, flags = compute_alpha_rad(P_des_b - B, beta)
        ik_ok = angles_feasible(alpha_des) and (not any(flags))

        if not ik_ok:
            status_text.set_text("IK infeasible — holding last command")
            # keep alpha_cmd unchanged
        else:
            # 4) Slew-limit command step (servo max speed)
            delta = alpha_des - alpha_cmd
            delta = (delta + np.pi) % (2*np.pi) - np.pi  # shortest turn
            rate_cap = np.radians(ALPHA_RATE_MAX_DEG) * DT
            delta = np.clip(delta, -rate_cap, +rate_cap)
            alpha_cmd = alpha_cmd + delta
            status_text.set_text("")

        # 5) Servo internal following (measured/encoder)
        alpha_meas = alpha_meas + (alpha_cmd - alpha_meas) * (1 - np.exp(-DT/TAU_SERVO))

        # 6) FK with measured angles (BASE frame) → achieved pose in BASE
        H_b = horn_endpoints_base(B, beta, alpha_meas, h)
        R_try_b, T_try_b = fk_iterative_6dof(H_b, p, d, R_seed, T_seed)

        # 7) Feasibility check (rods) in BASE frame
        P_try_b = (R_try_b @ p.T).T + T_try_b
        L = np.linalg.norm(P_try_b - H_b, axis=1)
        fk_ok = np.all(np.isfinite(L)) and np.all(np.abs(L - d) <= 1e-3)

        if fk_ok and angles_feasible(alpha_meas):
            R_tb, T_tb = R_try_b, T_try_b
            R_seed, T_seed = R_tb.copy(), T_tb.copy()
            last_ok = {
                "R_tb": R_tb.copy(),
                "T_tb": T_tb.copy(),
                "alpha_meas": alpha_meas.copy(),
                "alpha_cmd": alpha_cmd.copy()
            }
        else:
            # revert to last feasible
            R_tb, T_tb = last_ok["R_tb"], last_ok["T_tb"]
            alpha_meas = last_ok["alpha_meas"]
            alpha_cmd  = last_ok["alpha_cmd"]
            R_seed, T_seed = R_tb.copy(), T_tb.copy()
            status_text.set_text("FK/rod infeasible — clamped")

        # 8) Transform everything to WORLD for drawing:
        # Base anchors, horn tips, top anchors, and disks
        B_w = transform_points_only_R(R_bw, B)
        H_w = transform_points_only_R(R_bw, H_b)
        P_w = transform_points_only_R(R_bw, (R_tb @ p.T).T + T_tb)
        top_circle_w  = transform_points_only_R(R_bw, (R_tb @ circle_points(r_top).T).T + T_tb)
        base_circle_w = transform_points_only_R(R_bw, base_circle_b)

        # 9) Update polygons
        top_poly.set_verts([top_circle_w])
        base_poly.set_verts([base_circle_w])

        # 10) Frames in WORLD
        for ln in base_frame: ln.remove()
        for ln in top_frame:  ln.remove()
        base_frame[:] = draw_frame(ax3d, R_bw, [0,0,0], scale=3.0, linewidth=1.6, alpha=0.8)            # base frame (tilted)
        top_frame[:]  = draw_frame(ax3d, R_bw @ R_tb, R_bw @ T_tb, scale=3.0, linewidth=2.0, alpha=1.0)  # top frame

        # 11) Points & links in WORLD
        top_pts._offsets3d  = (P_w[:,0], P_w[:,1], P_w[:,2])
        base_pts._offsets3d = (B_w[:,0], B_w[:,1], B_w[:,2])
        for i in range(6):
            horns[i].set_data_3d([B_w[i,0], H_w[i,0]],
                                 [B_w[i,1], H_w[i,1]],
                                 [B_w[i,2], H_w[i,2]])
            rods[i].set_data_3d([H_w[i,0], P_w[i,0]],
                                [H_w[i,1], P_w[i,1]],
                                [H_w[i,2], P_w[i,2]])

        # 12) Servo angle bars (measured)
        a_deg = np.degrees(alpha_meas)
        for i, bar in enumerate(bars):
            bar.set_width(a_deg[i])
            bad = (a_deg[i] < SERVO_MIN) or (a_deg[i] > SERVO_MAX) or (not np.isfinite(a_deg[i]))
            bar.set_color('red' if bad else 'C0')

        # 13) IMU readout: achieved TOP orientation in WORLD should be ~ level
        R_tw = R_bw @ R_tb  # top relative to world
        roll_now, pitch_now = extract_roll_pitch_from_R(R_tw)
        imu_text.set_text(f"IMU (world): roll={np.degrees(roll_now):.1f}°, pitch={np.degrees(pitch_now):.1f}°  |  Gravity → -Z_world")

        fig.canvas.draw_idle()

    # Hook sliders
    s_wroll.on_changed(update)
    s_wpitch.on_changed(update)

    plt.show()

if __name__ == "__main__":
    run()