// forward_kinematics.cpp
// Recover roll/pitch from the 6 actuator angles by fitting the top plate pose.

#include <Arduino.h>
#include <math.h>

/* ---- Geometry provided by inverse_kinematics.cpp ---- */
// scalars
extern const float h;
extern const float d;
extern const float z0;

// arrays from IK geometry init
extern float beta_[6];   // horn azimuths (rad)
extern float Bx_[6];     // base joint x (cm)
extern float By_[6];     // base joint y (cm)
extern float px_[6];     // top plate local x (cm)
extern float py_[6];     // top plate local y (cm)

/* ---- small math ---- */
struct Vec3 { float x,y,z; };
struct Mat3 { float m[3][3]; };

static inline Vec3 v3(float x,float y,float z){ return {x,y,z}; }
static inline Vec3 add(Vec3 a, Vec3 b){ return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 sub(Vec3 a, Vec3 b){ return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 scale(Vec3 a, float s){ return {a.x*s, a.y*s, a.z*s}; }
static inline float dot(Vec3 a, Vec3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline float norm(Vec3 a){ return sqrtf(dot(a,a)); }

static inline Mat3 I(){ Mat3 R={{{1,0,0},{0,1,0},{0,0,1}}}; return R; }
static inline Vec3 mul(Mat3 R, Vec3 v){
  return {
    R.m[0][0]*v.x + R.m[0][1]*v.y + R.m[0][2]*v.z,
    R.m[1][0]*v.x + R.m[1][1]*v.y + R.m[1][2]*v.z,
    R.m[2][0]*v.x + R.m[2][1]*v.y + R.m[2][2]*v.z
  };
}
static inline Mat3 mm(Mat3 A, Mat3 B){
  Mat3 C={};
  for(int i=0;i<3;i++) for(int j=0;j<3;j++){
    float s=0; for(int k=0;k<3;k++) s+=A.m[i][k]*B.m[k][j];
    C.m[i][j]=s;
  }
  return C;
}
static inline Mat3 skew(Vec3 v){ return Mat3{{{0,-v.z,v.y},{v.z,0,-v.x},{-v.y,v.x,0}}}; }
static inline Mat3 madd(Mat3 A, Mat3 B, float s){ Mat3 C=A; for(int i=0;i<3;i++) for(int j=0;j<3;j++) C.m[i][j]+=s*B.m[i][j]; return C; }

static Mat3 so3_exp(Vec3 w){
  float th = norm(w);
  if(th < 1e-6f) return madd(I(), skew(w), 1.0f);
  Vec3 k = scale(w, 1.0f/th);
  Mat3 K=skew(k), I3=I(), K2=mm(K,K), R=I3;
  float s=sinf(th), c=cosf(th), a=1.0f-c;
  for(int i=0;i<3;i++) for(int j=0;j<3;j++){
    R.m[i][j]+=s*K.m[i][j];
    R.m[i][j]+=a*K2.m[i][j];
  }
  return R;
}

static inline void ui_times_mat_row(Vec3 ui,const Mat3&M,float out[3]){
  out[0]=ui.x*M.m[0][0]+ui.y*M.m[1][0]+ui.z*M.m[2][0];
  out[1]=ui.x*M.m[0][1]+ui.y*M.m[1][1]+ui.z*M.m[2][1];
  out[2]=ui.x*M.m[0][2]+ui.y*M.m[1][2]+ui.z*M.m[2][2];
}

/* Compute horn endpoints H[i] given base layout and alpha */
static void horn_endpoints_base(const float Bx[6], const float By[6],
                                const float beta[6], const float alpha[6],
                                float hlen, Vec3 H[6]){
  for(int i=0;i<6;i++){
    float ca=cosf(alpha[i]), sa=sinf(alpha[i]);
    float cb=cosf(beta[i]),  sb=sinf(beta[i]);
    H[i].x = Bx[i] + hlen*(ca*cb);
    H[i].y = By[i] + hlen*(ca*sb);
    H[i].z = 0.0f + hlen*(sa);
  }
}

/* 6×6 linear solver (Gauss-Jordan) */
static bool solve6(float A[6][6], float b[6], float x[6]){
  static float M[6][7];
  for(int i=0;i<6;i++){ for(int j=0;j<6;j++) M[i][j]=A[i][j]; M[i][6]=b[i]; }
  for(int i=0;i<6;i++){
    int piv=i; float maxv=fabsf(M[i][i]);
    for(int k=i+1;k<6;k++){ float v=fabsf(M[k][i]); if(v>maxv){maxv=v; piv=k;} }
    if(maxv < 1e-12f) return false;
    if(piv!=i){ for(int j=0;j<7;j++){ float t=M[i][j]; M[i][j]=M[piv][j]; M[piv][j]=t; } }
    float di=M[i][i]; for(int j=0;j<7;j++) M[i][j]/=di;
    for(int k=0;k<6;k++){ if(k==i) continue; float f=M[k][i]; for(int j=0;j<7;j++) M[k][j]-=f*M[i][j]; }
  }
  for(int i=0;i<6;i++) x[i]=M[i][6];
  return true;
}

/* Levenberg–Marquardt fit of pose (R,T) s.t. |R p_i + T - H_i| = d */
static void fk_iterative_6dof(const Vec3 H[6], const Vec3 p_local[6], float dlen,
                              Mat3 &R, Vec3 &T,
                              int iters=5, float lam=1e-3f, float tol_r=3e-5f, float tol_dx=3e-7f)
{
  static Vec3 P,v; static float L,ri; static float A[6][6], b[6], dx[6];

  for(int it=0; it<iters; ++it){
    float r2=0;
    for(int i=0;i<6;i++){
      P = add(mul(R, p_local[i]), T);
      v = sub(P, H[i]);
      L = norm(v);
      ri = L - dlen;
      r2 += ri*ri;
    }
    if (sqrtf(r2) < tol_r) break;

    for(int r=0;r<6;r++){ b[r]=0; for(int c=0;c<6;c++) A[r][c]=0; }
    for(int i=0;i<6;i++){
      P = add(mul(R, p_local[i]), T);
      v = sub(P, H[i]);
      L = norm(v);
      ri = L - dlen;

      Vec3 ui = (L>1e-12f) ? scale(v,1.0f/L) : v3(0,0,0);

      // J = [ -ui^T * R * [p]_x  |  ui^T ]
      Mat3 RK = mm(R, skew(p_local[i]));
      float jw[3]; ui_times_mat_row(ui, RK, jw);

      float Jrow[6]={-jw[0],-jw[1],-jw[2], ui.x, ui.y, ui.z};
      for(int r=0;r<6;r++){ b[r]+=Jrow[r]*ri; for(int c=0;c<6;c++) A[r][c]+=Jrow[r]*Jrow[c]; }
    }
    float lam2=lam*lam; for(int d=0;d<6;d++) A[d][d]+=lam2;
    for(int r=0;r<6;r++) b[r]=-b[r];

    if(!solve6(A,b,dx)) break;

    Mat3 dR = so3_exp( v3(dx[0],dx[1],dx[2]) );
    R = mm(R,dR); T = add(T, v3(dx[3],dx[4],dx[5]) );

    float d2=0; for(int k=0;k<6;k++) d2+=dx[k]*dx[k];
    if (sqrtf(d2)<tol_dx) break;
  }
}

static void extract_rp(const Mat3 &R, float &roll, float &pitch){
  // R = Ry(pitch)*Rx(roll)
  pitch = asinf(-R.m[2][0]);
  roll  = atan2f(-R.m[1][2], R.m[1][1]);
}

/* ---- Public API: recover roll/pitch from alpha ---- */
void fk_pose_from_alpha(const float alpha_rad[6], float& roll_deg_out, float& pitch_deg_out){
  // Horn endpoints at the base
  Vec3 H[6];
  horn_endpoints_base(Bx_, By_, beta_, alpha_rad, h, H);

  // Top plate local points
  Vec3 pTop[6];
  for(int i=0;i<6;i++) pTop[i] = v3(px_[i], py_[i], 0.0f);

  // Solve pose
  Mat3 R = I(); Vec3 T = v3(0,0,z0);
  fk_iterative_6dof(H, pTop, d, R, T, 5, 1e-3f, 3e-5f, 3e-7f);

  float r,p; extract_rp(R, r, p);
  roll_deg_out  = r * 180.0f/PI;
  pitch_deg_out = p * 180.0f/PI;
}