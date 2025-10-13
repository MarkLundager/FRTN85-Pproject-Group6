#include <Arduino.h>
#include <math.h>

/*** AVR-safe free RAM (no sbrk) ***/
static int freeRam() {
  extern char __bss_end;
  extern char* __brkval;
  char top;
  return (int) (&top - (__brkval ? __brkval : &__bss_end));
}

/*** Types ***/
struct Vec3 { float x,y,z; };
struct Mat3 { float m[3][3]; };

/*** Small math ***/
struct Lin {
  static inline Vec3 v3(float x,float y,float z){ return {x,y,z}; }
  static inline Vec3 add(Vec3 a, Vec3 b){ return {a.x+b.x, a.y+b.y, a.z+b.z}; }
  static inline Vec3 sub(Vec3 a, Vec3 b){ return {a.x-b.x, a.y-b.y, a.z-b.z}; }
  static inline Vec3 scale(Vec3 a, float s){ return {a.x*s, a.y*s, a.z*s}; }
  static inline float dot(Vec3 a, Vec3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
  static inline float norm(Vec3 a){ return sqrtf(dot(a,a)); }

  static inline Mat3 I(){ Mat3 R={{{1,0,0},{0,1,0},{0,0,1}}}; return R; }
  static inline Vec3 mul(Mat3 R, Vec3 v){
    return { R.m[0][0]*v.x + R.m[0][1]*v.y + R.m[0][2]*v.z,
             R.m[1][0]*v.x + R.m[1][1]*v.y + R.m[1][2]*v.z,
             R.m[2][0]*v.x + R.m[2][1]*v.y + R.m[2][2]*v.z };
  }
  static inline Mat3 mm(Mat3 A, Mat3 B){
    Mat3 C={};
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){
      float s=0; for(int k=0;k<3;k++) s += A.m[i][k]*B.m[k][j];
      C.m[i][j]=s;
    } return C;
  }
  static inline Mat3 madd(Mat3 A, Mat3 B, float s){ Mat3 C=A; for(int i=0;i<3;i++) for(int j=0;j<3;j++) C.m[i][j]+=s*B.m[i][j]; return C; }
  static inline Mat3 skew(Vec3 v){ return Mat3{{{ 0, -v.z, v.y},{ v.z, 0,-v.x},{-v.y, v.x, 0}}}; }
  static Mat3 so3_exp(Vec3 w){
    float th = norm(w);
    if(th < 1e-6f) return madd(I(), skew(w), 1.0f);
    Vec3 k = scale(w, 1.0f/th);
    Mat3 K = skew(k), I3 = I(), K2 = mm(K,K), R = I3;
    float s = sinf(th), c = cosf(th), a = 1.0f - c;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){ R.m[i][j]+=s*K.m[i][j]; R.m[i][j]+=a*K2.m[i][j]; }
    return R;
  }
  static inline Mat3 R_x(float r){ float c=cosf(r), s=sinf(r); return Mat3{{{1,0,0},{0,c,-s},{0,s,c}}}; }
  static inline Mat3 R_y(float p){ float c=cosf(p), s=sinf(p); return Mat3{{{ c,0,s},{0,1,0},{-s,0,c}}}; }
  static inline void ui_times_mat_row(Vec3 ui, const Mat3 &M, float out_row[3]){
    out_row[0]=ui.x*M.m[0][0]+ui.y*M.m[1][0]+ui.z*M.m[2][0];
    out_row[1]=ui.x*M.m[0][1]+ui.y*M.m[1][1]+ui.z*M.m[2][1];
    out_row[2]=ui.x*M.m[0][2]+ui.y*M.m[1][2]+ui.z*M.m[2][2];
  }
};

/*** Geometry ***/
static const float h_len = 1.64f, d_len = 17.0f, r_top = 7.8f, r_bot = 9.3f, z0 = 16.3f;
static const float pDelta = (15.25f * (float)PI / 180.0f);
static const float bDelta = (23.58f * (float)PI / 180.0f);

static float p_phi[6], b_phi[6], beta_arr[6];
static Vec3 B_base[6], p_top_arr[6];

static void build_layout(){
  for(int k=0;k<6;k++){
    float fk2 = (float)(k / 2);
    float fk1 = (float)((k + 1) / 2);
    float sgn = ((k & 1) ? -1.0f : 1.0f);
    p_phi[k]    = (2.0f*(float)PI/3.0f)*fk2 - sgn*(pDelta*0.5f) + ((float)PI/3.0f);
    b_phi[k]    = (2.0f*(float)PI/3.0f)*fk1 + sgn*(bDelta*0.5f);
    beta_arr[k] = b_phi[k] + (float)PI*0.5f*sgn;
    B_base[k]   = { r_bot*cosf(b_phi[k]), r_bot*sinf(b_phi[k]), 0.0f };
    p_top_arr[k]= { r_top*cosf(p_phi[k]), r_top*sinf(p_phi[k]), 0.0f };
  }
}

/*** IK ***/
struct IK {
  static void compute_alpha_rad(const Vec3 L[6], const float beta[6],
                                float alpha[6], uint8_t flags[6])
  {
    for(int i=0;i<6;i++){
      float lx=L[i].x, ly=L[i].y, lz=L[i].z;
      float e = 2.0f*h_len*lz;
      float f = 2.0f*h_len*(cosf(beta[i])*lx + sinf(beta[i])*ly);
      float g = (lx*lx + ly*ly + lz*lz) - (d_len*d_len - h_len*h_len);
      float A = sqrtf(e*e + f*f);
      if(A < 1e-7f){ alpha[i]=NAN; flags[i]=1; continue; }
      float val = g / (A + 1e-12f);
      if(val >  1.0f) val =  1.0f;
      if(val < -1.0f) val = -1.0f;
      alpha[i] = asinf(val) - atan2f(f, e);
      flags[i] = 0;
    }
  }
  static void horn_endpoints_base(const Vec3 B[6], const float beta[6],
                                  const float alpha[6], float h, Vec3 H[6])
  {
    for(int i=0;i<6;i++){
      float ca=cosf(alpha[i]), sa=sinf(alpha[i]);
      float cb=cosf(beta[i]),  sb=sinf(beta[i]);
      H[i].x = B[i].x + h*(ca*cb);
      H[i].y = B[i].y + h*(ca*sb);
      H[i].z = B[i].z + h*(sa);
    }
  }
};

/*** 6×6 solver (Gauss-Jordan) — uses a single 6×7 buffer ***/
static bool solve6(float A[6][6], float b[6], float x[6]){
  static float M[6][7];
  for(int i=0;i<6;i++){ for(int j=0;j<6;j++) M[i][j]=A[i][j]; M[i][6]=b[i]; }
  for(int i=0;i<6;i++){
    int piv=i; float maxv=fabsf(M[i][i]);
    for(int k=i+1;k<6;k++){ float v=fabsf(M[k][i]); if(v>maxv){maxv=v; piv=k;} }
    if(maxv < 1e-12f) return false;
    if(piv!=i){ for(int j=0;j<7;j++){ float t=M[i][j]; M[i][j]=M[piv][j]; M[piv][j]=t; } }
    float di=M[i][i];
    for(int j=0;j<7;j++) M[i][j]/=di;
    for(int k=0;k<6;k++){
      if(k==i) continue;
      float f=M[k][i];
      for(int j=0;j<7;j++) M[k][j]-=f*M[i][j];
    }
  }
  for(int i=0;i<6;i++) x[i]=M[i][6];
  return true;
}

/*** FK (Gauss–Newton) — low-RAM version (no J stored) ***/
struct FK {
  static void iterative_6dof(const Vec3 H[6], const Vec3 p_local[6], float dlen,
                             Mat3 &R, Vec3 &T,
                             int iters=6, float lam=1e-3f, float tol_r=3e-5f, float tol_dx=3e-7f)
  {
    static Vec3 P, v;              // per-leg temporaries
    static float L, ri;            // per-leg length & residual
    static float A[6][6], b[6], dx[6];

    for(int it=0; it<iters; ++it){
      // compute residual norm to early-break
      float rnorm2 = 0.0f;
      for(int i=0;i<6;i++){
        P = Lin::add(Lin::mul(R, p_local[i]), T);
        v = Lin::sub(P, H[i]);
        L = Lin::norm(v);
        ri = L - dlen;
        rnorm2 += ri*ri;
      }
      if (sqrtf(rnorm2) < tol_r) break;

      // zero normal eqs
      for(int r=0;r<6;r++){ b[r]=0; for(int c=0;c<6;c++) A[r][c]=0; }

      // accumulate A = J^T J, b = J^T r
      for(int i=0;i<6;i++){
        P = Lin::add(Lin::mul(R, p_local[i]), T);
        v = Lin::sub(P, H[i]);
        L = Lin::norm(v);
        ri = L - dlen;

        Vec3 ui = (L > 1e-12f) ? Lin::scale(v, 1.0f/L) : Lin::v3(0,0,0);

        // Jrow_w = ui^T * (R * skew(p_i))
        Mat3 RK = Lin::mm(R, Lin::skew(p_local[i]));
        float jw[3]; Lin::ui_times_mat_row(ui, RK, jw);
        float Jrow[6] = { -jw[0], -jw[1], -jw[2], ui.x, ui.y, ui.z };

        // A += Jrow^T * Jrow ; b += Jrow^T * ri
        for(int r=0;r<6;r++){
          b[r] += Jrow[r]*ri;
          for(int c=0;c<6;c++) A[r][c] += Jrow[r]*Jrow[c];
        }
      }

      // Levenberg-Marquardt damping: A += lam^2 * I
      float lam2 = lam*lam;
      for(int d=0; d<6; d++) A[d][d] += lam2;
      for(int r=0;r<6;r++) b[r] = -b[r];

      if(!solve6(A, b, dx)) break;

      Mat3 dR = Lin::so3_exp(Lin::v3(dx[0],dx[1],dx[2]));
      R = Lin::mm(R, dR);
      T = Lin::add(T, Lin::v3(dx[3],dx[4],dx[5]));

      float dnorm2=0; for(int k=0;k<6;k++) dnorm2 += dx[k]*dx[k];
      if(sqrtf(dnorm2) < tol_dx) break;
    }
  }

  static void extract_roll_pitch_from_R(const Mat3 &R, float &roll, float &pitch){
    pitch = asinf(-R.m[2][0]);
    roll  = atan2f(-R.m[1][2], R.m[1][1]);
  }
};

static float rad2deg(float r){ return r * 180.0f / (float)PI; }

/*** Test buffers ***/
static Mat3 R_bw, R_tb_des;
static Vec3 T_tb_des;
static Vec3 P_des[6];
static Vec3 Lvec[6];
static float alpha_sol[6];
static uint8_t flags_sol[6];
static Vec3 H_tip[6];

void setup(){
  pinMode(LED_BUILTIN, OUTPUT);
  Serial.begin(115200);
  delay(10);
  Serial.println(F("Boot OK"));

  build_layout();
  Serial.println(F("Layout OK"));
  Serial.print(F("Free RAM @start ≈ ")); Serial.println(freeRam());

  // Desired ~5°/5°
  float roll  = 5.0f * (float)PI / 180.0f;
  float pitch = 5.0f * (float)PI / 180.0f;

  // R_tb_des = (R_y * R_x)^T
  R_bw = Lin::mm(Lin::R_y(pitch), Lin::R_x(roll));
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) R_tb_des.m[i][j] = R_bw.m[j][i];
  T_tb_des = {0,0,z0};

  // Desired top points and leg vectors
  for(int i=0;i<6;i++) P_des[i] = Lin::add(Lin::mul(R_tb_des, p_top_arr[i]), T_tb_des);
  for(int i=0;i<6;i++) Lvec[i]  = Lin::sub(P_des[i], B_base[i]);

  // IK
  IK::compute_alpha_rad(Lvec, beta_arr, alpha_sol, flags_sol);
  bool ik_ok = true;
  for(int i=0;i<6;i++){
    if(!isfinite(alpha_sol[i])) ik_ok=false;
    float deg = rad2deg(alpha_sol[i]);
    if(deg < -90.0f || deg > 90.0f) ik_ok=false;
  }
  Serial.print(F("IK OK: ")); Serial.println(ik_ok ? F("yes") : F("NO"));
  if(!ik_ok){ Serial.println(F("Aborting before FK.")); return; }

  // Build horn tips from IK
  IK::horn_endpoints_base(B_base, beta_arr, alpha_sol, h_len, H_tip);

  // FK
  Serial.print(F("Free RAM pre-FK ≈ ")); Serial.println(freeRam());
  Serial.println(F("Entering FK..."));
  Mat3 R_fk = Lin::I();
  Vec3 T_fk = {0,0,z0};
  FK::iterative_6dof(H_tip, p_top_arr, d_len, R_fk, T_fk, 6, 1e-3f, 3e-5f, 3e-7f);

  float roll_rec, pitch_rec;
  FK::extract_roll_pitch_from_R(R_fk, roll_rec, pitch_rec);
  Serial.print(F("FK roll/pitch [deg]=("));
  Serial.print(rad2deg(roll_rec),3); Serial.print(',');
  Serial.print(rad2deg(pitch_rec),3); Serial.println(')');

  // Errors
  float max_len_err=0.0f, max_P_err=0.0f;
  for(int i=0;i<6;i++){
    Vec3 P_fk = Lin::add(Lin::mul(R_fk, p_top_arr[i]), T_fk);
    float Li = Lin::norm(Lin::sub(P_fk, H_tip[i]));
    float le = fabsf(Li - d_len); if(le > max_len_err) max_len_err = le;
    float pe = Lin::norm(Lin::sub(P_fk, P_des[i])); if(pe > max_P_err) max_P_err = pe;
  }
  Serial.print(F("Max rod-length error [cm]: ")); Serial.println(max_len_err,6);
  Serial.print(F("Max point mismatch |P_fk-P_des| [cm]: ")); Serial.println(max_P_err,6);

  Serial.print(F("Free RAM end ≈ ")); Serial.println(freeRam());
  Serial.println(F("Done."));
  digitalWrite(LED_BUILTIN, HIGH);
}

void loop(){
  static uint32_t t=0;
  if (millis() - t > 500) { t = millis(); digitalWrite(LED_BUILTIN, !digitalRead(LED_BUILTIN)); }
}