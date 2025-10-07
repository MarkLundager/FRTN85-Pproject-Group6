# visualize.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry import (h, d, r_top, r_bot, z0, B, p, beta,
                      SERVO_MIN, SERVO_MAX, ROD_TOL)
from utils import R_x, R_y, circle_points, transform_points
from inverse_kinematics import pose, compute_alpha_rad, horn_endpoints


def check_constraints(P, H, alpha_rad, flags_in):
    msgs = []
    ok = True
    alpha_deg = np.degrees(alpha_rad)
    for i in range(6):
        leg_msgs = []
        if flags_in[i]:
            ok = False
            leg_msgs.append(flags_in[i])
        if np.isfinite(alpha_rad[i]):
            rod_len = np.linalg.norm(P[i] - H[i])
            if abs(rod_len - d) > ROD_TOL:
                ok = False
                leg_msgs.append(f'|HP|={rod_len:.3f} (d={d})')
            if not (SERVO_MIN - 1e-6 <= alpha_deg[i] <= SERVO_MAX + 1e-6):
                ok = False
                leg_msgs.append(f'servo {alpha_deg[i]:.1f}° out [{SERVO_MIN:.0f},{SERVO_MAX:.0f}]')
        else:
            ok = False
            leg_msgs.append('alpha NaN')
        msgs.append('; '.join(leg_msgs))
    if np.any(P[:, 2] < 0.0):
        ok = False
        msgs.append('top anchor below base plane (z<0)')
    return ok, msgs


def draw_frame(ax, R, T, scale=3.0, linewidth=1.5):
    """Draws a small coordinate frame (X=red, Y=green, Z=blue)."""
    origin = np.array(T)
    axes = np.eye(3) * scale
    colors = ['r', 'g', 'b']
    lines = []
    for i in range(3):
        vec = origin + R @ axes[:, i]
        line, = ax.plot([origin[0], vec[0]],
                        [origin[1], vec[1]],
                        [origin[2], vec[2]],
                        color=colors[i], linewidth=linewidth)
        lines.append(line)
    return lines


def run():
    # --- Figure setup ---
    fig = plt.figure(figsize=(11, 8))
    ax3d = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.08, right=0.78, bottom=0.22)
    fig.patch.set_facecolor('white')
    ax3d.set_facecolor('white')
    ax3d.set_axis_off()

    # --- Initial pose ---
    roll0 = np.deg2rad(0.0)
    pitch0 = np.deg2rad(0.0)
    R, T, P = pose(roll0, pitch0, z0, p, R_y, R_x)
    L = P - B
    alpha_rad, flags0 = compute_alpha_rad(L, beta)
    H = horn_endpoints(B, beta, alpha_rad, h)

    # --- Base and top disks ---
    base_circle = circle_points(r_bot)
    base_poly = Poly3DCollection([base_circle], alpha=0.18,
                                 facecolor=(0, 0, 0, 0.10), edgecolor='none')
    ax3d.add_collection3d(base_poly)

    top_circle_local = circle_points(r_top)
    top_circle_world = transform_points(R, T, top_circle_local)
    top_poly = Poly3DCollection([top_circle_world], alpha=0.25,
                                facecolor=(0, 0, 1, 0.15), edgecolor='none')
    ax3d.add_collection3d(top_poly)

    # --- Coordinate frames (base fixed, top moves) ---
    draw_frame(ax3d, np.eye(3), [0, 0, 0], scale=3.0)  # base frame
    top_frame = draw_frame(ax3d, R, T, scale=3.0)      # top frame (dynamic)

    # --- Points and rods ---
    base_pts = ax3d.scatter(B[:, 0], B[:, 1], B[:, 2], s=8, c='k')
    top_pts = ax3d.scatter(P[:, 0], P[:, 1], P[:, 2], s=8, c='b')

    horns = [ax3d.plot([B[i, 0], H[i, 0]],
                       [B[i, 1], H[i, 1]],
                       [B[i, 2], H[i, 2]], linewidth=2)[0] for i in range(6)]
    rods = [ax3d.plot([H[i, 0], P[i, 0]],
                      [H[i, 1], P[i, 1]],
                      [H[i, 2], P[i, 2]], linewidth=2, color='k')[0] for i in range(6)]

    # --- Camera view ---
    ax3d.set_xlim(-12, 12)
    ax3d.set_ylim(-12, 12)
    ax3d.set_zlim(0, 25)
    ax3d.set_box_aspect([1, 1, 1])
    ax3d.view_init(elev=20, azim=-60)

    # --- Servo angle bar chart ---
    ax_bar = fig.add_axes([0.82, 0.25, 0.16, 0.65])
    yidx = np.arange(6)
    angles_deg0 = np.degrees(alpha_rad)
    bars = ax_bar.barh(yidx, angles_deg0, align='center')
    ax_bar.set_yticks(yidx)
    ax_bar.set_yticklabels([f'k={i}' for i in range(6)])
    ax_bar.axvline(0, linewidth=0.8)
    ax_bar.set_xlim(-90, 90)
    ax_bar.set_xlabel('alpha [deg]')
    ax_bar.set_title("Servo angles")

    status_text = fig.text(0.08, 0.02, "", color='red')

    # --- Sliders ---
    ax_roll = plt.axes([0.15, 0.14, 0.5, 0.03])
    ax_pitch = plt.axes([0.15, 0.09, 0.5, 0.03])
    slider_roll = Slider(ax_roll, 'Roll [°]', -15, 15, valinit=0)
    slider_pitch = Slider(ax_pitch, 'Pitch [°]', -15, 15, valinit=0)

    def update(_):
        roll = np.radians(slider_roll.val)
        pitch = np.radians(slider_pitch.val)

        R, T, P = pose(roll, pitch, z0, p, R_y, R_x)
        L = P - B
        alpha_rad, flags = compute_alpha_rad(L, beta)
        H = horn_endpoints(B, beta, alpha_rad, h)

        # Update top disk
        top_world = transform_points(R, T, top_circle_local)
        top_poly.set_verts([top_world])

        # Update moving top frame
        for line in top_frame:
            line.remove()
        top_frame[:] = draw_frame(ax3d, R, T, scale=3.0, linewidth=1.5)

        # Update points and rods
        top_pts._offsets3d = (P[:, 0], P[:, 1], P[:, 2])
        for i in range(6):
            if np.isfinite(alpha_rad[i]):
                horns[i].set_data_3d([B[i, 0], H[i, 0]],
                                     [B[i, 1], H[i, 1]],
                                     [B[i, 2], H[i, 2]])
                rods[i].set_data_3d([H[i, 0], P[i, 0]],
                                    [H[i, 1], P[i, 1]],
                                    [H[i, 2], P[i, 2]])
            else:
                horns[i].set_data_3d([B[i, 0], B[i, 0]],
                                     [B[i, 1], B[i, 1]],
                                     [B[i, 2], B[i, 2]])
                rods[i].set_data_3d([B[i, 0], B[i, 0]],
                                    [B[i, 1], B[i, 1]],
                                    [B[i, 2], B[i, 2]])

        # Update bars
        alpha_deg = np.degrees(alpha_rad)
        for i, bar in enumerate(bars):
            width = 0.0 if not np.isfinite(alpha_deg[i]) else alpha_deg[i]
            bar.set_width(width)
            if not np.isfinite(alpha_deg[i]):
                bar.set_color('gray')
            elif (alpha_deg[i] < SERVO_MIN) or (alpha_deg[i] > SERVO_MAX):
                bar.set_color('red')
            else:
                bar.set_color('C0')

        # Constraint checks
        ok, msgs = check_constraints(P, H, alpha_rad, flags)
        status_text.set_text('' if ok else '\n'.join([f'Leg {i}: {msgs[i]}' for i in range(6) if msgs[i]][:3]))

        fig.canvas.draw_idle()

    slider_roll.on_changed(update)
    slider_pitch.on_changed(update)
    plt.show()
