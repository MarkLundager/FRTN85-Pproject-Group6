#pragma once
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Geometry constants (cm)
#define H_SERVO 1.64f   // servo horn length
#define D_ROD   17.0f   // rod length (fixed)
#define R_TOP   7.8f    // top platform radius
#define R_BOT   9.3f    // base radius
#define Z0      16.3f   // nominal height

// Pair separations (radians)
#define P_DELTA (15.25f * M_PI / 180.0f)  // platform pair separation
#define B_DELTA (23.58f * M_PI / 180.0f)  // base pair separation

// Limits and tolerances
#define SERVO_MIN -90.0f
#define SERVO_MAX  90.0f
#define IK_TOL     1e-9f
#define ROD_TOL    1e-3f   // cm tolerance

#define N_LEGS 6

// Global arrays
static float B[N_LEGS][3];     // base anchors
static float p[N_LEGS][3];     // top anchors
static float p_phi[N_LEGS];    // top platform angles
static float b_phi[N_LEGS];    // base platform angles
static float beta[N_LEGS];     // servo tangent orientations

// -----------------------------------------------------------------------------
// Geometry initialization
// -----------------------------------------------------------------------------
static inline void geometry_init() {
  for (int k = 0; k < N_LEGS; ++k) {
    float pfloor = floor(k / 2.0f);
    float bfloor = floor((k + 1) / 2.0f);
    float sign = (k % 2 == 0) ? 1.0f : -1.0f; // replaces pow(-1,k)

    p_phi[k] = (2.0f * M_PI / 3.0f) * pfloor - sign * (P_DELTA / 2.0f) + M_PI / 3.0f;
    b_phi[k] = (2.0f * M_PI / 3.0f) * bfloor + sign * (B_DELTA / 2.0f);
    beta[k]  = b_phi[k] + (M_PI / 2.0f) * sign;

    // Base anchors
    B[k][0] = R_BOT * cos(b_phi[k]);
    B[k][1] = R_BOT * sin(b_phi[k]);
    B[k][2] = 0.0f;

    // Top anchors
    p[k][0] = R_TOP * cos(p_phi[k]);
    p[k][1] = R_TOP * sin(p_phi[k]);
    p[k][2] = 0.0f;
  }
}