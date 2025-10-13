#include "FK.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- minimal LA helpers (private to FK.cpp) ---
static inline void mat3_identity(float R[9]){ R[0]=1;R[1]=0;R[2]=0; R[3]=0;R[4]=1;R[5]=0; R[6]=0;R[7]=0;R[8]=1; }
static inline void mat3_mul(const float A[9], const float B[9], float C[9]){
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
static inline void mat3_vec_mul(const float A[9], const float v[3], float out[3]){
  out[0]=A[0]*v[0]+A[1]*v[1]+A[2]*v[2];
  out[1]=A[3]*v[0]+A[4]*v[1]+A[5]*v[2];
  out[2]=A[6]*v[0]+A[7]*v[1]+A[8]*v[2];
}
static inline float vnorm3(const float v[3]){ return sqrtf(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
static inline float vdot3(const float a[3], const float b[3]){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static inline void vsub3(const float a[3], const float b[3], float o[3]){ o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2]; }
static inline void skew(const float v[3], float S[9]){
  const float x=v[0], y=v[1], z=v[2];
  S[0]=0;S[1]=-z;S[2]=y; S[3]=z;S[4]=0;S[5]=-x; S[6]=-y;S[7]=x;S[8]=0;
}
static inline void so3_exp(const float w[3], float R[9]){
  float th = vnorm3(w);
  if (th < 1e-12f) { mat3_identity(R); float K[9]; skew(w,K); for(int i=0;i<9;i++) R[i]+=K[i]; return; }
  float k[3] = { w[0]/th, w[1]/th, w[2]/th };
  float K[9], K2[9]; skew(k,K); mat3_mul(K,K,K2);
  mat3_identity(R);
  float s=sinf(th), c=cosf(th), a=(1.0f-c);
  for (int i=0;i<9;i++) R[i]+= s*K[i] + a*K2[i];
}
static bool solve6x6(float A[36], float b[6], float x[6]){
  float M[6][7];
  for (int r=0;r<6;r++){ for(int c=0;c<6;c++) M[r][c]=A[r*6+c]; M[r][6]=b[r]; }
  for (int k=0;k<6;k++){
    int piv=k; float maxv=fabsf(M[k][k]);
    for (int r=k+1;r<6;r++){ float v=fabsf(M[r][k]); if(v>maxv){maxv=v;piv=r;} }
    if (maxv < 1e-12f) return false;
    if (piv!=k){ for(int c=k;c<7;c++){ float t=M[k][c]; M[k][c]=M[piv][c]; M[piv][c]=t; } }
    float d=M[k][k];
    for (int c=k;c<7;c++) M[k][c]/=d;
    for (int r=k+1;r<6;r++){ float f=M[r][k]; for (int c=k;c<7;c++) M[r][c]-=f*M[k][c]; }
  }
  for (int r=5;r>=0;r--){ float s=M[r][6]; for (int c=r+1;c<6;c++) s-=M[r][c]*x[c]; x[r]=s; }
  return true;
}

// --- FK main ---
bool fk_iterative_6dof(const float H[6][3], const float p_top[6][3],
                       const float d_len[6], float R[9], float T[3],
                       int iters, float lam, float tol_r, float tol_dx)
{
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
