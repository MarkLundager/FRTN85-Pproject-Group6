#pragma once
#include <Arduino.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Geometry (cm / rad) like your geometry.py ---
extern const float h_servo;
extern const float d_rod;
extern const float r_top;
extern const float r_bot;
extern const float z0;

extern const float pDelta;     // radians
extern const float bDelta;     // radians
extern const float SERVO_MIN_DEG;
extern const float SERVO_MAX_DEG;

inline float deg2rad(float d){ return d * (float)M_PI / 180.0f; }
inline float rad2deg(float r){ return r * 180.0f / (float)M_PI; }

// Build base anchors B, top anchors p_top (in top frame), and servo tangent orientation beta
void build_B_p_beta(float B[6][3], float p_top[6][3], float beta[6]);

// Horn endpoints in base frame H = B + h * [cosβ cosα, sinβ cosα, sinα]
void horn_endpoints_base(const float B[6][3], const float beta[6],
                         const float alpha_rad[6], float h_len,
                         float H[6][3]);