#include "IK.h"
#include "geometry.h"
#include <math.h>

void compute_alpha_rad(const float v[6][3], const float beta[6],
                       float h_len, float d_fix,
                       float alpha[6], bool feasible[6])
{
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

    float phi  = atan2f(Bz, A);
    float ainc = acosf(cval);

    auto inLim = [](float a)->bool{
      float deg = a * 180.0f / (float)M_PI;
      return deg >= SERVO_MIN_DEG && deg <= SERVO_MAX_DEG;
    };

    float a1 = phi + ainc;
    float a2 = phi - ainc;

    if (inLim(a1)) { alpha[i]=a1; feasible[i]=true; }
    else if (inLim(a2)) { alpha[i]=a2; feasible[i]=true; }
    else { alpha[i]=a1; feasible[i]=false; }
  }
}