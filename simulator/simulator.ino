#include <Arduino.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------- Geometry (cm/rad) ----------
static const float h_servo = 1.64f;
static const float d_rod   = 17.0f;
static const float r_top   = 7.8f;
static const float r_bot   = 9.3f;
static const float z0      = 16.3f;

static const float pDelta  = 15.25f * (float)M_PI / 180.0f;
static const float bDelta  = 23.58f * (float)M_PI / 180.0f;
static const float SERVO_MIN_DEG = -90.0f, SERVO_MAX_DEG = 90.0f;

// ---------- Small LA helpers ----------
inline void mat3_identity(float R[9]) {
  R[0]=1;R[1]=0;R[2]=0; R[3]=0;R[4]=1;R[5]=0; R[6]=0;R[7]=0;R[8]=1;
}
inline void mat3_mul(const float A[9], const float B[9], float C[9]){
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
inline void mat3_vec_mul(const float A[9], const float v[3], float out[3]){
  out[0]=A[0]*v[0]+A[1]*v[1]+A[2]*v[2];
  out[1]=A[3]*v[0]+A[4]*v[1]+A[5]*v[2];
  out[2]=A[6]*v[0]+A[7]*v[1]+A[8]*v[2];
}
inline float vdot3(const float a[3], const float b[3]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
inline float vnorm3(const float v[3]){
  return sqrtf(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
inline void vsub3(const float a[3], const float b[3], float o[3]){
  o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
inline void skew(const float v[3], float S[9]){
  const float x=v[0], y=v[1], z=v[2];
  S[0]=0; S[1]=-z; S[2]=y; S[3]=z; S[4]=0; S[5]=-x; S[6]=-y; S[7]=x; S[8]=0;
}
inline void so3_exp(const float w[3], float R[9]){
  float th = vnorm3(w);
  if (th < 1e-12f) { mat3_identity(R); float K[9]; skew(w,K); for(int i=0;i<9;i++) R[i]+=K[i]; return; }
  float k[3] = { w[0]/th, w[1]/th, w[2]/th };
  float K[9], K2[9]; skew(k,K); mat3_mul(K,K,K2);
  mat3_identity(R);
  float s=sinf(th), c=cosf(th), a=(1.0f-c);
  for (int i=0;i<9;i++) R[i]+= s*K[i] + a*K2[i];
}

// ---------- 6x6 solver ----------
bool solve6x6(float A[36], float b[6], float x[6]) {
  float M[6][7];
  for (int r=0;r<6;r++){ for(int c=0;c<6;c++) M[r][c]=A[r*6+c]; M[r][6]=b[r]; }
  for (int k=0;k<6;k++){
    int piv=k; float maxv=fabsf(M[k][k]);
    for (int r=k+1;r<6;r++){ float v=fabsf(M[r][k]); if(v>maxv){maxv=v;piv=r;} }
    if (maxv < 1e-12f) return false;
    if (piv!=k){ for(int c=k;c<7;c++){ float t=M[k][c]; M[k][c]=M[piv][c]; M[piv][c]=t; } }
    float d=M[k][k];
    for (int c=k;c<7;c++) M[k][c]/=d;
    for (int r=k+1;r<6;r++){ float f=M[r][k]; for(int c=k;c<7;c++) M[r][c]-=f*M[k][c]; }
  }
  for (int r=5;r>=0;r--){ float s=M[r][6]; for(int c=r+1;c<6;c++) s-=M[r][c]*x[c]; x[r]=s; }
  return true;
}

// ---------- FK (Gauss-Newton) ----------
bool fk_iterative_6dof(
  const float H[6][3], const float p_top[6][3], const float d_len[6],
  float R[9], float T[3], int iters=25, float lam=1e-3f, float tol_r=1e-3f, float tol_dx=1e-8f
){
  for (int it=0; it<iters; ++it){
    float v[6][3], L[6], r[6], P[3];
    for (int i=0;i<6;i++){
      float pt[3] = { p_top[i][0], p_top[i][1], p_top[i][2] };
      mat3_vec_mul(R, pt, P);
      P[0]+=T[0]; P[1]+=T[1]; P[2]+=T[2];
      vsub3(P, H[i], v[i]);
      L[i] = fmaxf(vnorm3(v[i]), 1e-12f);
      r[i] = L[i] - d_len[i];
    }
    float maxr=0; for (int i=0;i<6;i++) maxr=fmaxf(maxr, fabsf(r[i]));
    if (maxr < tol_r) return true;

    float J[6][6];
    for (int i=0;i<6;i++){
      float ui[3] = { v[i][0]/L[i], v[i][1]/L[i], v[i][2]/L[i] };
      float S[9], RS[9], M[9];
      float pti[3] = { p_top[i][0], p_top[i][1], p_top[i][2] };
      skew(pti, S);
      // M = -R*S
      // RS = R*S
      mat3_mul(R, S, RS);
      for (int k=0;k<9;k++) M[k] = -RS[k];
      float c0[3]={M[0],M[3],M[6]}, c1[3]={M[1],M[4],M[7]}, c2[3]={M[2],M[5],M[8]};
      J[i][0]=vdot3(ui,c0); J[i][1]=vdot3(ui,c1); J[i][2]=vdot3(ui,c2);
      J[i][3]=ui[0]; J[i][4]=ui[1]; J[i][5]=ui[2];
    }
    float A[36]={0}, bvec[6]={0};
    for (int m=0;m<6;m++){
      for (int n=0;n<6;n++){
        float s=0; for (int i=0;i<6;i++) s+=J[i][m]*J[i][n];
        A[m*6+n]=s + ((m==n)?(lam*lam):0.0f);
      }
      float sb=0; for (int i=0;i<6;i++) sb+=J[i][m]*r[i];
      bvec[m] = -sb;
    }
    float dx[6]={0}; if(!solve6x6(A,bvec,dx)) return false;

    float dR[9], Rn[9], w[3]={dx[0],dx[1],dx[2]};
    so3_exp(w,dR); mat3_mul(R,dR,Rn); for(int i=0;i<9;i++) R[i]=Rn[i];
    T[0]+=dx[3]; T[1]+=dx[4]; T[2]+=dx[5];

    float dn = sqrtf(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+dx[3]*dx[3]+dx[4]*dx[4]+dx[5]*dx[5]);
    if (dn < tol_dx) return true;
  }
  return true;
}

// ---------- Geometry builders ----------
void build_B_p_beta(float B[6][3], float p_top[6][3], float beta[6]) {
  for (int k=0;k<6;k++){
    float fk2 = floorf(k/2.0f);
    float fk1 = floorf((k+1)/2.0f);
    float sgn = (k%2==0)? 1.0f : -1.0f;
    float p_phi = (2.0f*(float)M_PI/3.0f)*fk2 - sgn*(pDelta*0.5f) + (float)M_PI/3.0f;
    float b_phi = (2.0f*(float)M_PI/3.0f)*fk1 + sgn*(bDelta*0.5f);
    beta[k]     = b_phi + (float)M_PI*0.5f*sgn;
    B[k][0]=r_bot*cosf(b_phi); B[k][1]=r_bot*sinf(b_phi); B[k][2]=0.0f;
    p_top[k][0]=r_top*cosf(p_phi); p_top[k][1]=r_top*sinf(p_phi); p_top[k][2]=0.0f;
  }
}

void horn_endpoints_base(const float B[6][3], const float beta[6],
                         const float alpha[6], float h_len, float H[6][3]){
  for (int i=0;i<6;i++){
    float ca=cosf(alpha[i]), sa=sinf(alpha[i]);
    float cb=cosf(beta[i]), sb=sinf(beta[i]);
    H[i][0] = B[i][0] + h_len*(ca*cb);
    H[i][1] = B[i][1] + h_len*(ca*sb);
    H[i][2] = B[i][2] + h_len*(sa);
  }
}

// ---------- Closed-form IK for alpha ----------
// Given v = P - B (leg direction from base anchor) and beta,
// solve A cosα + B sinα = S, where
// A = cb*v_x + sb*v_y, B = v_z, S = (h^2 + ||v||^2 - d^2)/(2h).
// Returns alpha (rad) and a feasibility flag per leg.
void compute_alpha_rad(
  const float v[6][3], const float beta[6],
  float h_len, float d_fix,
  float alpha[6], bool feasible[6]
){
  for (int i=0;i<6;i++){
    float cb=cosf(beta[i]), sb=sinf(beta[i]);
    float vx=v[i][0], vy=v[i][1], vz=v[i][2];
    float A = cb*vx + sb*vy;
    float Bz = vz;
    float vv = vx*vx + vy*vy + vz*vz;
    float S = (h_len*h_len + vv - d_fix*d_fix) / (2.0f*h_len);
    float R = sqrtf(A*A + Bz*Bz);

    if (R < 1e-12f) { feasible[i]=false; alpha[i]=0; continue; }
    float cval = S / R;
    if (cval < -1.0f || cval > 1.0f) { feasible[i]=false; alpha[i]=0; continue; }

    float phi = atan2f(Bz, A);         // cos phi = A/R, sin phi = B/R
    float acosv = acosf(cval);

    // Two solutions: alpha = phi ± acos(S/R). Pick one within servo limits (prefer +).
    float a1 = phi + acosv;
    float a2 = phi - acosv;
    auto inLim = [](float a)->bool{
      float deg = a * 180.0f / (float)M_PI;
      return deg >= SERVO_MIN_DEG && deg <= SERVO_MAX_DEG;
    };
    if (inLim(a1)) { alpha[i]=a1; feasible[i]=true; }
    else if (inLim(a2)) { alpha[i]=a2; feasible[i]=true; }
    else { alpha[i]=a1; feasible[i]=false; }
  }
}

// ---------- Convenience rotations ----------
void R_x(float rx, float R[9]){
  float c=cosf(rx), s=sinf(rx);
  R[0]=1; R[1]=0; R[2]=0;
  R[3]=0; R[4]=c; R[5]=-s;
  R[6]=0; R[7]=s; R[8]=c;
}
void R_y(float ry, float R[9]){
  float c=cosf(ry), s=sinf(ry);
  R[0]=c; R[1]=0; R[2]=s;
  R[3]=0; R[4]=1; R[5]=0;
  R[6]=-s;R[7]=0; R[8]=c;
}

// ---------- Compare IK→FK ----------
float R_true[9], T_true[3];

void setup(){
  Serial.begin(115200);
  while(!Serial){} // for boards that need it
  delay(200);

  float B[6][3], Ptop[6][3], Beta[6];
  build_B_p_beta(B, Ptop, Beta);
  float dlen[6]; for (int i=0;i<6;i++) dlen[i]=d_rod;

  // Choose a test pose (adjust as desired)
  float roll = 10.0f * (float)M_PI/180.0f;
  float pitch= -7.0f * (float)M_PI/180.0f;

  float Rx[9], Ry[9], Rtmp[9];
  R_x(roll, Rx); R_y(pitch, Ry);                // R = Ry * Rx
  mat3_mul(Ry, Rx, R_true);
  T_true[0]=0.0f; T_true[1]=0.0f; T_true[2]=z0; // at nominal height

  // P = R*p_top + T
  float P[6][3];
  for (int i=0;i<6;i++){
    float pr[3]; mat3_vec_mul(R_true, Ptop[i], pr);
    P[i][0]=pr[0]+T_true[0]; P[i][1]=pr[1]+T_true[1]; P[i][2]=pr[2]+T_true[2];
  }

  // v = P - B
  float v[6][3]; for (int i=0;i<6;i++) vsub3(P[i], B[i], v[i]);

  // ---- IK: compute alpha ----
  float alpha[6]; bool feas[6];
  compute_alpha_rad(v, Beta, h_servo, d_rod, alpha, feas);

  // Horn tips H from IK alpha
  float H[6][3];
  horn_endpoints_base(B, Beta, alpha, h_servo, H);

  // ---- FK: recover pose from H ----
  float R_fk[9]; mat3_identity(R_fk);
  float T_fk[3] = {0.0f, 0.0f, z0};   // seed near truth
  bool ok = fk_iterative_6dof(H, Ptop, dlen, R_fk, T_fk, 30, 1e-3f, 1e-3f, 1e-8f);

  // ---- Compare poses ----
  // Position error
  float dT[3] = { T_fk[0]-T_true[0], T_fk[1]-T_true[1], T_fk[2]-T_true[2] };
  float pos_err = vnorm3(dT);

  // Orientation error (angle of R_true^T * R_fk)
  // Compute R_rel = R_true^T * R_fk
  float Rt[9] = { R_true[0], R_true[3], R_true[6],
                  R_true[1], R_true[4], R_true[7],
                  R_true[2], R_true[5], R_true[8] };
  float Rrel[9]; mat3_mul(Rt, R_fk, Rrel);
  float tr = Rrel[0] + Rrel[4] + Rrel[8];
  float angle_err = acosf(fmaxf(-1.0f, fminf(1.0f, (tr - 1.0f) * 0.5f))); // radians

  // ---- Print results ----
  Serial.println(F("=== IK → FK comparison ==="));
  Serial.print(F("FK ok: ")); Serial.println(ok ? F("true") : F("false"));

  Serial.println(F("Alpha [deg] and feasibility:"));
  for (int i=0;i<6;i++){
    float adeg = alpha[i]*180.0f/(float)M_PI;
    Serial.print(F("  k=")); Serial.print(i);
    Serial.print(F(": ")); Serial.print(adeg, 3);
    Serial.print(F(" deg  feasible="));
    Serial.println(feas[i] ? F("true") : F("false"));
  }

  Serial.print(F("Position error |T_fk - T_true| = "));
  Serial.print(pos_err, 6); Serial.println(F(" cm"));

  Serial.print(F("Orientation error angle = "));
  Serial.print(angle_err*180.0f/(float)M_PI, 6); Serial.println(F(" deg"));

  Serial.println(F("\nTrue T (cm):"));
  Serial.print(T_true[0],6); Serial.print(' '); Serial.print(T_true[1],6); Serial.print(' '); Serial.println(T_true[2],6);

  Serial.println(F("FK T (cm):"));
  Serial.print(T_fk[0],6); Serial.print(' '); Serial.print(T_fk[1],6); Serial.print(' '); Serial.println(T_fk[2],6);
}

void loop(){}