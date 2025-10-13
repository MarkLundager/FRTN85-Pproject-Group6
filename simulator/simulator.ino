#include <Arduino.h>
#include "geometry.h"
#include "IK.h"
#include "FK.h"

// simple rotations
static void R_x(float rx, float R[9]){
  float c=cosf(rx), s=sinf(rx);
  R[0]=1; R[1]=0; R[2]=0;
  R[3]=0; R[4]=c; R[5]=-s;
  R[6]=0; R[7]=s; R[8]=c;
}
static void R_y(float ry, float R[9]){
  float c=cosf(ry), s=sinf(ry);
  R[0]=c; R[1]=0; R[2]=s;
  R[3]=0; R[4]=1; R[5]=0;
  R[6]=-s;R[7]=0; R[8]=c;
}
static void mat3_mul(const float A[9], const float B[9], float C[9]){
  C[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
  C[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
  C[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
  C[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
  C[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
  C[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
  C[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
  C[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
  C[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}
static void mat3_vec_mul(const float A[9], const float v[3], float out[3]){
  out[0]=A[0]*v[0]+A[1]*v[1]+A[2]*v[2];
  out[1]=A[3]*v[0]+A[4]*v[1]+A[5]*v[2];
  out[2]=A[6]*v[0]+A[7]*v[1]+A[8]*v[2];
}
static void vsub3(const float a[3], const float b[3], float o[3]){ o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2]; }
static float vnorm3(const float v[3]){ return sqrtf(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

void setup(){
  Serial.begin(115200);
  while(!Serial){}  // some boards need this
  delay(200);

  // Build geometry
  float B[6][3], Ptop[6][3], Beta[6];
  build_B_p_beta(B, Ptop, Beta);
  float dlen[6]; for (int i=0;i<6;i++) dlen[i]=d_rod;

  // Choose a test pose
  float roll  = deg2rad(10.0f);
  float pitch = deg2rad(-7.0f);
  float Rx[9], Ry[9], R_true[9];
  R_x(roll, Rx); R_y(pitch, Ry);
  mat3_mul(Ry, Rx, R_true);
  float T_true[3] = { 0.0f, 0.0f, z0 };

  // P = R_true * p_top + T_true
  float P[6][3];
  for (int i=0;i<6;i++){
    float pr[3]; mat3_vec_mul(R_true, Ptop[i], pr);
    P[i][0]=pr[0]+T_true[0];
    P[i][1]=pr[1]+T_true[1];
    P[i][2]=pr[2]+T_true[2];
  }

  // IK: alpha from v = P - B
  float v[6][3]; for (int i=0;i<6;i++) vsub3(P[i], B[i], v[i]);
  float alpha[6]; bool feas[6];
  compute_alpha_rad(v, Beta, h_servo, d_rod, alpha, feas);

  // Horn endpoints from IK
  float H[6][3];
  horn_endpoints_base(B, Beta, alpha, h_servo, H);

  // FK: recover pose
  float R_fk[9] = {1,0,0, 0,1,0, 0,0,1};
  float T_fk[3] = {0.0f, 0.0f, z0};
  bool ok = fk_iterative_6dof(H, Ptop, dlen, R_fk, T_fk, 30, 1e-3f, 1e-3f, 1e-8f);

  // Compare
  float dT[3] = { T_fk[0]-T_true[0], T_fk[1]-T_true[1], T_fk[2]-T_true[2] };
  float pos_err = vnorm3(dT);

  // orientation error: angle of R_true^T * R_fk
  float Rt[9] = { R_true[0], R_true[3], R_true[6],
                  R_true[1], R_true[4], R_true[7],
                  R_true[2], R_true[5], R_true[8] };
  float Rrel[9]; mat3_mul(Rt, R_fk, Rrel);
  float tr = Rrel[0] + Rrel[4] + Rrel[8];
  float ang_rad = acosf(fmaxf(-1.0f, fminf(1.0f, (tr - 1.0f) * 0.5f)));

  // Print
  Serial.println(F("=== IK & FK comparison (modular) ==="));
  Serial.print(F("FK ok: ")); Serial.println(ok ? F("true") : F("false"));
  Serial.println(F("Alpha [deg] and feasibility:"));
  for (int i=0;i<6;i++){
    Serial.print(F("  k=")); Serial.print(i);
    Serial.print(F(": ")); Serial.print(rad2deg(alpha[i]), 3);
    Serial.print(F("  feasible=")); Serial.println(feas[i] ? F("true") : F("false"));
  }
  Serial.print(F("Pos error (cm): ")); Serial.println(pos_err, 6);
  Serial.print(F("Ori error (deg): ")); Serial.println(rad2deg(ang_rad), 6);

  Serial.print(F("T_true: ")); Serial.print(T_true[0],3); Serial.print(' ');
  Serial.print(T_true[1],3); Serial.print(' '); Serial.println(T_true[2],3);

  Serial.print(F("T_fk:   ")); Serial.print(T_fk[0],3); Serial.print(' ');
  Serial.print(T_fk[1],3); Serial.print(' '); Serial.println(T_fk[2],3);
}

void loop(){}