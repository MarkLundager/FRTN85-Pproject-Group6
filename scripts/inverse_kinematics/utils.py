# utils.py
import numpy as np

def R_x(a):
    return np.array([[1,0,0],
                     [0,np.cos(a),-np.sin(a)],
                     [0,np.sin(a), np.cos(a)]])

def R_y(a):
    return np.array([[ np.cos(a),0, np.sin(a)],
                     [0,          1, 0         ],
                     [-np.sin(a),0, np.cos(a) ]])

def rpy_to_R(roll, pitch):
    # yaw fixed to 0
    return R_y(pitch) @ R_x(roll)

def circle_points(radius, n=256):
    t = np.linspace(0, 2*np.pi, n, endpoint=True)
    return np.column_stack([radius*np.cos(t), radius*np.sin(t), np.zeros(n)])

def transform_points(R, T, Pts):
    # P' = R*P + T
    return (R @ Pts.T).T + T
