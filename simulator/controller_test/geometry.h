#pragma once
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#define H_SERVO 1.64f   
#define D_ROD   17.0f   
#define R_TOP   7.8f    
#define R_BOT   9.3f    
#define Z0      16.3f   


#define P_DELTA (15.25f * M_PI / 180.0f)  
#define B_DELTA (23.58f * M_PI / 180.0f)  


#define SERVO_MIN -90.0f
#define SERVO_MAX  90.0f
#define IK_TOL     1e-9f
#define ROD_TOL    1e-3f   

#define N_LEGS 6


static float B[N_LEGS][3];     
static float p[N_LEGS][3];     
static float p_phi[N_LEGS];    
static float b_phi[N_LEGS];    
static float beta[N_LEGS];     




static inline void geometry_init() {
  for (int k = 0; k < N_LEGS; ++k) {
    float pfloor = floor(k / 2.0f);
    float bfloor = floor((k + 1) / 2.0f);
    float sign = (k % 2 == 0) ? 1.0f : -1.0f; 

    p_phi[k] = (2.0f * M_PI / 3.0f) * pfloor - sign * (P_DELTA / 2.0f) + M_PI / 3.0f;
    b_phi[k] = (2.0f * M_PI / 3.0f) * bfloor + sign * (B_DELTA / 2.0f);
    beta[k]  = b_phi[k] + (M_PI / 2.0f) * sign;

    
    B[k][0] = R_BOT * cos(b_phi[k]);
    B[k][1] = R_BOT * sin(b_phi[k]);
    B[k][2] = 0.0f;

    
    p[k][0] = R_TOP * cos(p_phi[k]);
    p[k][1] = R_TOP * sin(p_phi[k]);
    p[k][2] = 0.0f;
  }
}