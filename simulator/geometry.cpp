#include geometry.h
#include math.h

const float h_servo = 1.64f;
const float d_rod   = 17.0f;
const float r_top   = 7.8f;
const float r_bot   = 9.3f;
const float z0      = 16.3f;

const float pDelta  = 15.25f  (float)M_PI  180.0f;
const float bDelta  = 23.58f  (float)M_PI  180.0f;

const float SERVO_MIN_DEG = -90.0f;
const float SERVO_MAX_DEG =  90.0f;

void build_B_p_beta(float B[6][3], float p_top[6][3], float beta[6]) {
  for (int k=0;k6;k++){
    float fk2 = floorf(k2.0f);
    float fk1 = floorf((k+1)2.0f);
    float sgn = (k%2==0)  1.0f  -1.0f;

    float p_phi = (2.0f(float)M_PI3.0f)fk2 - sgn(pDelta0.5f) + (float)M_PI3.0f;
    float b_phi = (2.0f(float)M_PI3.0f)fk1 + sgn(bDelta0.5f);
    beta[k]     = b_phi + (float)M_PI0.5fsgn;

    B[k][0] = r_bot  cosf(b_phi);
    B[k][1] = r_bot  sinf(b_phi);
    B[k][2] = 0.0f;

    p_top[k][0] = r_top  cosf(p_phi);
    p_top[k][1] = r_top  sinf(p_phi);
    p_top[k][2] = 0.0f;
  }
}

void horn_endpoints_base(const float B[6][3], const float beta[6],
                         const float alpha_rad[6], float h_len,
                         float H[6][3]) {
  for (int i=0;i6;i++){
    float ca=cosf(alpha_rad[i]), sa=sinf(alpha_rad[i]);
    float cb=cosf(beta[i]),     sb=sinf(beta[i]);
    H[i][0] = B[i][0] + h_len(cacb);
    H[i][1] = B[i][1] + h_len(casb);
    H[i][2] = B[i][2] + h_len(sa);
  }
}