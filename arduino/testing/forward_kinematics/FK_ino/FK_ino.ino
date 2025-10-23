#include <Arduino.h>
#include <math.h>

// ---------------- Geometry (cm / radians) ----------------
static const float h_servo = 1.64f;   // servo horn length (cm)
static const float d_rod   = 17.0f;   // rod length (fixed, cm)
static const float r_top   = 7.8f;    // top platform radius (cm)
static const float r_bot   = 9.3f;    // base radius (cm)
static const float z0      = 16.3f;   // nominal height (cm)

static const float pDelta_deg = 15.25f;  // platform pair separation (deg)
static const float bDelta_deg = 23.58f;  // base pair separation (deg)
static const float pDelta = pDelta_deg * (float)M_PI / 180.0f;
static const float bDelta = bDelta_deg * (float)M_PI / 180.0f;

// Limits / tolerances
static const float SERVO_MIN_DEG = -90.0f;
static const float SERVO_MAX_DEG =  90.0f;
static const float IK_TOL   = 1e-9f;   // not used directly here, but kept for reference
static const float ROD_TOL  = 1e-3f;   // cm tolerance on |HP| == d

// ---------- Math helpers (same as before) ----------
static const float EPS = 1e-12f;

inline float deg2rad(float d){ return d * (float)M_PI / 180.0f; }

inline float vnorm3(const float v[3]) {
  return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline float vdot3(const float a[3], const float b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void vadd3(const float a[3], const float b[3], float out[3]) {
  out[0]=a[0]+b[0]; out[1]=a[1]+b[1]; out[2]=a[2]+b[2];
}

inline void vsub3(const float a[3], const float b[3], float out[3]) {
  out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2];
}

inline void mat3_identity(float R[9]) {
  R[0]=1; R[1]=0; R[2]=0;
  R[3]=0; R[4]=1; R[5]=0;
  R[6]=0; R[7]=0; R[8]=1;
}

inline void mat3_copy(const float A[9], float B[9]) {
  for (int i=0;i<9;i++) B[i]=A[i];
}

inline void mat3_mul(const float A[9], const float B[9], float C[9]) {
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

inline void mat3_vec_mul(const float A[9], const float v[3], float out[3]) {
  out[0]=A[0]*v[0]+A[1]*v[1]+A[2]*v[2];
  out[1]=A[3]*v[0]+A[4]*v[1]+A[5]*v[2];
  out[2]=A[6]*v[0]+A[7]*v[1]+A[8]*v[2];
}

inline void skew(const float v[3], float S[9]) {
  const float x=v[0], y=v[1], z=v[2];
  S[0]=0;   S[1]=-z;  S[2]=y;
  S[3]=z;   S[4]=0;   S[5]=-x;
  S[6]=-y;  S[7]=x;   S[8]=0;
}

inline void so3_exp(const float omega[3], float R_out[9]) {
  float th = vnorm3(omega);
  if (th < 1e-12f) {
    float K[9];
    skew(omega, K);
    mat3_identity(R_out);
    for (int i=0;i<9;i++) R_out[i] += K[i];
    return;
  }
  float k[3] = { omega[0]/th, omega[1]/th, omega[2]/th };
  float K[9], K2[9];
  skew(k, K);
  mat3_mul(K, K, K2);
  mat3_identity(R_out);
  float s = sinf(th), c = cosf(th), a = (1.0f - c);
  for (int i=0;i<9;i++) R_out[i] += s*K[i] + a*K2[i];
}

void build_B_p_beta(float B[6][3], float p_top[6][3], float beta[6]) {
  for (int k=0;k<6;k++){
    float fk2 = floorf(k / 2.0f);
    float fk1 = floorf((k+1) / 2.0f);
    float sgn = (k % 2 == 0) ? 1.0f : -1.0f;

    float p_phi = (2.0f*(float)M_PI/3.0f)*fk2 - sgn*(pDelta*0.5f) + (float)M_PI/3.0f;
    float b_phi = (2.0f*(float)M_PI/3.0f)*fk1 + sgn*(bDelta*0.5f);
    beta[k]     = b_phi + (float)M_PI*0.5f*sgn;

    // base anchors B
    B[k][0] = r_bot * cosf(b_phi);
    B[k][1] = r_bot * sinf(b_phi);
    B[k][2] = 0.0f;

    // top anchors p_top (in top frame)
    p_top[k][0] = r_top * cosf(p_phi);
    p_top[k][1] = r_top * sinf(p_phi);
    p_top[k][2] = 0.0f;
  }
}

// Horn endpoints in base frame: H = B + h * [cos(beta)*cos(alpha), sin(beta)*cos(alpha), sin(alpha)]
void horn_endpoints_base(const float B[6][3], const float beta[6],
                         const float alpha_rad[6], float h_len,
                         float H[6][3]) {
  for (int i=0;i<6;i++){
    float ca = cosf(alpha_rad[i]), sa = sinf(alpha_rad[i]);
    float cb = cosf(beta[i]),     sb = sinf(beta[i]);
    H[i][0] = B[i][0] + h_len * (ca * cb);
    H[i][1] = B[i][1] + h_len * (ca * sb);
    H[i][2] = B[i][2] + h_len * (sa);
  }
}

void printMatrix3(const float R[9]) {
  for (int r=0;r<3;r++){
    Serial.print(R[r*3+0], 6); Serial.print(' ');
    Serial.print(R[r*3+1], 6); Serial.print(' ');
    Serial.println(R[r*3+2], 6);
  }
}




// 6x6 solver
bool solve6x6(float A[36], float b[6], float x[6]) {
  float M[6][7];
  for (int r=0;r<6;r++){
    for (int c=0;c<6;c++) M[r][c] = A[r*6+c];
    M[r][6] = b[r];
  }
  for (int k=0;k<6;k++){
    int piv = k;
    float maxv = fabsf(M[k][k]);
    for (int r=k+1;r<6;r++){
      float v = fabsf(M[r][k]);
      if (v > maxv) { maxv=v; piv=r; }
    }
    if (maxv < 1e-12f) return false;
    if (piv != k) {
      for (int c=k;c<7;c++) { float tmp=M[k][c]; M[k][c]=M[piv][c]; M[piv][c]=tmp; }
    }
    float diag = M[k][k];
    for (int c=k;c<7;c++) M[k][c] /= diag;
    for (int r=k+1;r<6;r++){
      float f = M[r][k];
      for (int c=k;c<7;c++) M[r][c] -= f * M[k][c];
    }
  }
  for (int r=5;r>=0;r--){
    float sum = M[r][6];
    for (int c=r+1;c<6;c++) sum -= M[r][c]*x[c];
    x[r] = sum;
  }
  return true;
}

// Gauss-Newton pose solve
bool fk_iterative_6dof(
  const float H[6][3],
  const float p_top[6][3],
  const float d_len[6],
  float R[9], float T[3],
  int iters = 25,
  float lam = 1e-3f,
  float tol_r = ROD_TOL,
  float tol_dx = 1e-8f
) {
  for (int iter=0; iter<iters; ++iter) {
    float v[6][3], L[6], r[6], P[3];
    for (int i=0;i<6;i++){
      float pt[3] = { p_top[i][0], p_top[i][1], p_top[i][2] };
      mat3_vec_mul(R, pt, P);
      P[0]+=T[0]; P[1]+=T[1]; P[2]+=T[2];
      vsub3(P, H[i], v[i]);
      L[i] = fmaxf(vnorm3(v[i]), 1e-12f);
      r[i] = L[i] - d_len[i];
    }

    // Stop if residuals are all within tolerance (max norm)
    float max_abs_r = 0.0f;
    for (int i=0;i<6;i++) max_abs_r = fmaxf(max_abs_r, fabsf(r[i]));
    if (max_abs_r < tol_r) return true;

    float J[6][6];
    for (int i=0;i<6;i++){
      float ui[3] = { v[i][0]/L[i], v[i][1]/L[i], v[i][2]/L[i] };

      float S[9], RS[9], M[9];
      float pti[3] = { p_top[i][0], p_top[i][1], p_top[i][2] };
      skew(pti, S);
      mat3_mul(R, S, RS);
      for (int k=0;k<9;k++) M[k] = -RS[k];

      float col0[3] = { M[0], M[3], M[6] };
      float col1[3] = { M[1], M[4], M[7] };
      float col2[3] = { M[2], M[5], M[8] };
      J[i][0] = vdot3(ui, col0);
      J[i][1] = vdot3(ui, col1);
      J[i][2] = vdot3(ui, col2);
      J[i][3] = ui[0];
      J[i][4] = ui[1];
      J[i][5] = ui[2];
    }

    float A[36] = {0}, bvec[6] = {0};
    for (int m=0;m<6;m++){
      for (int n=0;n<6;n++){
        float s=0.0f; for (int i=0;i<6;i++) s += J[i][m]*J[i][n];
        A[m*6+n] = s + ((m==n)?(lam*lam):0.0f);
      }
      float sb=0.0f; for (int i=0;i<6;i++) sb += J[i][m]*r[i];
      bvec[m] = -sb;
    }

    float dx[6]={0};
    if (!solve6x6(A, bvec, dx)) return false;

    float dR[9], Rnew[9];
    float dx_rot[3] = { dx[0], dx[1], dx[2] };
    so3_exp(dx_rot, dR);
    mat3_mul(R, dR, Rnew);
    mat3_copy(Rnew, R);
    T[0]+=dx[3]; T[1]+=dx[4]; T[2]+=dx[5];

    float dxn = 0.0f; for (int k=0;k<6;k++) dxn += dx[k]*dx[k];
    if (sqrtf(dxn) < tol_dx) return true;
  }
  return true;
}

// Extract roll/pitch (same convention as your Python)
void extract_roll_pitch_from_R(const float R[9], float &roll, float &pitch) {
  float r20 = R[6];  // R[2,0]
  float r12 = R[5];  // R[1,2]
  float r11 = R[4];  // R[1,1]
  pitch = asinf(-r20);
  roll  = atan2f(-r12, r11);
}

// ---------- Globals for pose ----------
float R_pose[9];
float T_pose[3];

// ---------- Build geometry from your formulas ----------
void build_geometry(float H[6][3], float p_top[6][3], float d_len[6]) {
  // k = 0..5
  for (int i=0;i<6;i++){
    int k = i;

    // floor(k/2)
    float fk2 = floorf(k / 2.0f);

    // (-1)^k  -> +1 for even, -1 for odd
    float sgn = (k % 2 == 0) ? 1.0f : -1.0f;

    // p_phi = (2π/3)*floor(k/2) - ((-1)^k)*(pDelta/2) + π/3
    float p_phi = (2.0f*(float)M_PI/3.0f)*fk2 - sgn*(pDelta*0.5f) + (float)M_PI/3.0f;

    // b_phi = (2π/3)*floor((k+1)/2) + ((-1)^k)*(bDelta/2)
    float fk1 = floorf((k+1) / 2.0f);
    float b_phi = (2.0f*(float)M_PI/3.0f)*fk1 + sgn*(bDelta*0.5f);

    // beta = b_phi + (π/2)*((-1)^k)
    float beta = b_phi + (float)M_PI*0.5f*sgn;

    // Anchors in local frames
    float B[3] = { r_bot * cosf(b_phi), r_bot * sinf(b_phi), 0.0f };
    float pL[3]= { r_top * cosf(p_phi), r_top * sinf(p_phi), 0.0f };

    // Servo horn tip position in base frame: H = B + h * [cos(beta), sin(beta), 0]
    H[i][0] = B[0] + h_servo * cosf(beta);
    H[i][1] = B[1] + h_servo * sinf(beta);
    H[i][2] = 0.0f;

    // Platform anchor in top frame
    p_top[i][0] = pL[0];
    p_top[i][1] = pL[1];
    p_top[i][2] = 0.0f;

    // Rod length (fixed)
    d_len[i] = d_rod;
  }
}




void setup(){
  Serial.begin(115200);
  delay(300);

  // 1) Build geometry bits used by FK
  float B[6][3], Ptop[6][3], Beta[6];
  build_B_p_beta(B, Ptop, Beta);

  // 2) Choose servo angles alpha (deg). For a quick sanity test, zeroes work well.
  //    You can tweak any of these to test different horn orientations.
  float alpha_deg[6] = { 0, 0, 0, 0, 0, 0 };
  float alpha_rad[6];
  for (int i=0;i<6;i++) alpha_rad[i] = deg2rad(alpha_deg[i]);

  // 3) Compute horn endpoints H for these alphas
  float H[6][3];
  horn_endpoints_base(B, Beta, alpha_rad, h_servo, H);

  // 4) Rod lengths (fixed): all d_rod
  float d_len[6]; for (int i=0;i<6;i++) d_len[i] = d_rod;

  // 5) Initial pose guess (top frame wrt base): identity, at z0
  float R0[9]; mat3_identity(R0);
  float T0[3] = { 0.0f, 0.0f, z0 };

  // 6) Run FK
  bool ok = fk_iterative_6dof(H, Ptop, d_len, R0, T0, 30, 1e-3f, 1e-3f, 1e-8f);

  Serial.println(ok ? F("FK: OK") : F("FK: FAILED"));
  Serial.println(F("R ="));
  printMatrix3(R0);
  Serial.print(F("T = "));
  Serial.print(T0[0], 6); Serial.print(' ');
  Serial.print(T0[1], 6); Serial.print(' ');
  Serial.println(T0[2], 6);

  float roll, pitch;
  extract_roll_pitch_from_R(R0, roll, pitch);
  Serial.print(F("roll  = ")); Serial.println(roll, 6);
  Serial.print(F("pitch = ")); Serial.println(pitch, 6);
}

void loop(){}