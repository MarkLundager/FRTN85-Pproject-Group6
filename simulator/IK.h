#pragma once
#include <Arduino.h>

// Closed-form IK for servo angles alpha
// v = P - B (leg vector in base frame), beta = servo tangent angle
// h_len = horn length, d_fix = rod length
// Outputs alpha[6] (radians) and feasible[6] (inside limits and solvable)
void compute_alpha_rad(const float v[6][3], const float beta[6],
                       float h_len, float d_fix,
                       float alpha[6], bool feasible[6]);