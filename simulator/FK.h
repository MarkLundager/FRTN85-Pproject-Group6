#pragma once
#include <Arduino.h>

// Forward kinematics via Gauss-Newton / Levenberg-Marquardt
bool fk_iterative_6dof(const float H[6][3], const float p_top[6][3],
                       const float d_len[6],
                       float R[9], float T[3],
                       int iters=25, float lam=1e-3f,
                       float tol_r=1e-3f, float tol_dx=1e-8f);