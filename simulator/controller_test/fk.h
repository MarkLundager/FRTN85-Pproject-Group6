#pragma once
#include <math.h>


static inline float dot3(const float a[3], const float b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline float norm3(const float a[3]) { return sqrtf(dot3(a,a)); }

static inline void add3(const float a[3], const float b[3], float out[3]) {
  out[0]=a[0]+b[0]; out[1]=a[1]+b[1]; out[2]=a[2]+b[2];
}
static inline void sub3(const float a[3], const float b[3], float out[3]) {
  out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2];
}
static inline void scale3(const float a[3], float s, float out[3]) {
  out[0]=a[0]*s; out[1]=a[1]*s; out[2]=a[2]*s;
}


static inline void mat3_mul_vec(const float A[3][3], const float v[3], float out[3]) {
  out[0] = A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
  out[1] = A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
  out[2] = A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2];
}

static inline void mat3_mul(const float A[3][3], const float B[3][3], float out[3][3]) {
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      out[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
    }
  }
}

static inline void mat3_eye(float I[3][3]) {
  I[0][0]=1; I[0][1]=0; I[0][2]=0;
  I[1][0]=0; I[1][1]=1; I[1][2]=0;
  I[2][0]=0; I[2][1]=0; I[2][2]=1;
}

static inline void skew3(const float v[3], float S[3][3]) {
  const float x=v[0], y=v[1], z=v[2];
  S[0][0]=0;   S[0][1]=-z;  S[0][2]= y;
  S[1][0]= z;  S[1][1]= 0;  S[1][2]=-x;
  S[2][0]=-y;  S[2][1]= x;  S[2][2]= 0;
}


static inline void so3_exp(const float omega[3], float dR[3][3]) {
  float th = norm3(omega);
  float K[3][3];
  if (th < 1e-12f) {
    
    float S[3][3]; skew3(omega,S);
    mat3_eye(dR);
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) dR[i][j] += S[i][j];
    return;
  }
  float k[3] = { omega[0]/th, omega[1]/th, omega[2]/th };
  skew3(k, K);
  float KK[3][3]; mat3_mul(K,K,KK);
  float I[3][3]; mat3_eye(I);
  float s = sinf(th), c = cosf(th);
  
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      dR[i][j] = I[i][j] + s*K[i][j] + (1.0f - c)*KK[i][j];
    }
  }
}


static inline void mat3_right_mul(float R[3][3], const float dR[3][3]) {
  float tmp[3][3]; mat3_mul(R,dR,tmp);
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) R[i][j]=tmp[i][j];
}


static inline bool solve6(float A[6][6], float b[6], float x[6]) {
  
  float M[6][7];
  for (int i=0;i<6;i++){
    for (int j=0;j<6;j++) M[i][j]=A[i][j];
    M[i][6]=b[i];
  }
  
  for (int col=0; col<6; col++){
    
    int piv = col;
    float best = fabsf(M[col][col]);
    for (int r=col+1; r<6; r++){
      float v = fabsf(M[r][col]);
      if (v > best){ best=v; piv=r; }
    }
    if (best < 1e-12f) return false; 
    if (piv != col){
      for (int j=col; j<7; j++){
        float t = M[piv][j]; M[piv][j]=M[col][j]; M[col][j]=t;
      }
    }
    
    float diag = M[col][col];
    for (int j=col; j<7; j++) M[col][j] /= diag;
    
    for (int r=col+1; r<6; r++){
      float f = M[r][col];
      if (f==0) continue;
      for (int j=col; j<7; j++) M[r][j] -= f * M[col][j];
    }
  }
  
  for (int i=5; i>=0; i--){
    float s = M[i][6];
    for (int j=i+1; j<6; j++) s -= M[i][j]*x[j];
    x[i] = s; 
  }
  return true;
}


/**
 * Inputs:
 *  H[6][3]      : horn tip positions (base frame), from current servo angles
 *  p_top[6][3]  : top anchor positions (top frame)
 *  d_len        : rod length (same for all legs)
 *  R_init[3][3] : initial rotation estimate
 *  T_init[3]    : initial translation estimate
 * Params:
 *  iters, lam, tol_r, tol_dx
 * Outputs:
 *  R_out[3][3], T_out[3]
 */
static inline void fk_iterative_6dof(
  const float H[6][3],
  const float p_top[6][3],
  float d_len,
  const float R_init[3][3],
  const float T_init[3],
  int iters,
  float lam,
  float tol_r,
  float tol_dx,
  float R_out[3][3],
  float T_out[3]
){
  
  float R[3][3]; for (int i=0;i<3;i++) for (int j=0;j<3;j++) R[i][j]=R_init[i][j];
  float T[3] = { T_init[0], T_init[1], T_init[2] };

  const float eps = 1e-12f;

  for (int it=0; it<iters; ++it){
    
    float P[6][3];
    for (int i=0;i<6;i++){
      float Rp[3]; mat3_mul_vec(R, p_top[i], Rp);
      add3(Rp, T, P[i]);
    }

    
    float v[6][3], L[6], r[6];
    float rnorm2 = 0.0f;
    for (int i=0;i<6;i++){
      sub3(P[i], H[i], v[i]);
      L[i] = norm3(v[i]);
      r[i] = L[i] - d_len;
      rnorm2 += r[i]*r[i];
    }
    if (sqrtf(rnorm2) < tol_r) break;

    
    float J[6][6];
    for (int i=0;i<6;i++){
      
      float ui[3] = {0,0,0};
      float denom = (L[i] > eps) ? L[i] : eps;
      ui[0] = v[i][0]/denom; ui[1] = v[i][1]/denom; ui[2] = v[i][2]/denom;

      float S[3][3]; skew3(p_top[i], S);
      float RS[3][3]; mat3_mul(R, S, RS);
      float M[3][3];
      for (int a=0;a<3;a++) for (int b=0;b<3;b++) M[a][b] = -RS[a][b];

      for (int j=0;j<3;j++){
        float col[3] = { M[0][j], M[1][j], M[2][j] };
        J[i][j] = dot3(ui, col);
      }
      J[i][3] = ui[0];
      J[i][4] = ui[1];
      J[i][5] = ui[2];
    }

    
    float A[6][6]; for (int a=0;a<6;a++) for (int b=0;b<6;b++) A[a][b]=0.0f;
    float g[6] = {0,0,0,0,0,0};
    for (int i=0;i<6;i++){
      for (int a=0;a<6;a++){
        g[a] += J[i][a] * r[i];
        for (int b=0;b<6;b++){
          A[a][b] += J[i][a] * J[i][b];
        }
      }
    }
    for (int a=0;a<6;a++) A[a][a] += lam*lam;

    
    float rhs[6]; for (int a=0;a<6;a++) rhs[a] = -g[a];
    float dx[6] = {0,0,0,0,0,0};
    bool ok = solve6(A, rhs, dx);
    if (!ok) break;

    
    float dR[3][3]; so3_exp(dx, dR);
    mat3_right_mul(R, dR);
    T[0] += dx[3]; T[1] += dx[4]; T[2] += dx[5];

    
    float dxnorm = 0.0f; for (int a=0;a<6;a++) dxnorm += dx[a]*dx[a];
    if (sqrtf(dxnorm) < tol_dx) break;
  }

  
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) R_out[i][j]=R[i][j];
  T_out[0]=T[0]; T_out[1]=T[1]; T_out[2]=T[2];
}


static inline void extract_roll_pitch_from_R(const float R[3][3], float &roll, float &pitch) {
  pitch = asinf(-R[2][0]);
  roll  = atan2f(-R[1][2], R[1][1]);
}