#include <Arduino.h>
#include <math.h>

// ---------- Geometry constants (cm) ----------
static const float h     = 1.64f;
static const float d     = 17.0f;
static const float r_top = 7.8f;
static const float r_bot = 9.3f;
static const float z0    = 16.3f;

// ---------- Layout angles (rad) ----------
static const float pDelta = 15.25f * PI / 180.0f;
static const float bDelta = 23.58f * PI / 180.0f;

// ---------- Precomputed geometry ----------
static float beta_[6], b_phi_[6], p_phi_[6];
static float Bx_[6], By_[6];
static float px_[6], py_[6];

static void rotXY(float roll, float pitch,
                  const float px[], const float py[], const float pz[],
                  float outx[], float outy[], float outz[]) {
  float cr = cosf(roll),  sr = sinf(roll);
  float cp = cosf(pitch), sp = sinf(pitch);
  for (int i = 0; i < 6; i++) {
    float x = px[i], y = py[i], z = pz[i];
    // R = Ry(pitch)*Rx(roll)
    float x1 =  cp*x + sp*z;
    float y1 =  sr*sp*x + cr*y - sr*cp*z;
    float z1 = -cr*sp*x + sr*y + cr*cp*z;
    outx[i] = x1; outy[i] = y1; outz[i] = z1 + z0;
  }
}

void ik_init_geometry() {
  for (int k = 0; k < 6; k++) {
    p_phi_[k] = (2*PI/3)*floorf(k/2.0f)   - powf(-1.0f,k)*(pDelta/2.0f) + PI/3.0f;
    b_phi_[k] = (2*PI/3)*floorf((k+1)/2.0f) + powf(-1.0f,k)*(bDelta/2.0f);
    beta_[k]  = b_phi_[k] + (PI/2.0f)*powf(-1.0f,k);

    Bx_[k] = r_bot * cosf(b_phi_[k]);
    By_[k] = r_bot * sinf(b_phi_[k]);
    px_[k] = r_top * cosf(p_phi_[k]);
    py_[k] = r_top * sinf(p_phi_[k]);
  }
}

// roll_deg, pitch_deg are in degrees (same convention as controller)
// Output alpha in radians
void ik_compute(float roll_deg, float pitch_deg, float alpha_rad_out[6]) {
  const float roll  = roll_deg  * PI / 180.0f;
  const float pitch = pitch_deg * PI / 180.0f;

  float pz[6] = {0,0,0,0,0,0};
  float Px[6], Py[6], Pz[6];
  rotXY(roll, pitch, px_, py_, pz, Px, Py, Pz);

  for (int i = 0; i < 6; i++) {
    float lx = Px[i] - Bx_[i];
    float ly = Py[i] - By_[i];
    float lz = Pz[i];

    float e = 2*h*lz;
    float f = 2*h*(cosf(beta_[i])*lx + sinf(beta_[i])*ly);
    float g = lx*lx + ly*ly + lz*lz - (d*d - h*h);
    float A = sqrtf(e*e + f*f);

    float val = g / (A > 1e-9f ? A : 1e-9f);
    if (val >  1.0f) val =  1.0f;
    if (val < -1.0f) val = -1.0f;

    alpha_rad_out[i] = asinf(val) - atan2f(f, e);
  }
}