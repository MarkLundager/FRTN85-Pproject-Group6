#include <Arduino.h>
#include <math.h>
#include "controller.h"

namespace {
static const float KP_ROLL  = 3.0f;
static const float KP_PITCH = 3.0f;

static const float KI_ROLL  = 0.4f;
static const float KI_PITCH = 0.40f;
static const float ICLAMP   = 15.0f;

static const float KD_ROLL  = 1.20f;
static const float KD_PITCH = 1.20f;
static const float DERIV_LPF_ALPHA = 0.95f;

static const float DEAD_BAND_DEG = 0.3f;
static const float MAX_REF_RATE_DEG_S = 10.0f;   // if you keep slew limit
static const float REF_LIMIT_DEG = 12.0f;


static float err_roll_prev  = 0.0f;
static float err_pitch_prev = 0.0f;
static float derr_roll_filt  = 0.0f;
static float derr_pitch_filt = 0.0f;
static float ref_roll_prev = 0.0f;
static float ref_pitch_prev = 0.0f;
static float i_roll = 0.0f;
static float i_pitch = 0.0f;
}

// ---------- helpers ----------
inline float apply_deadband(float err, float threshold) {
  return (fabs(err) < threshold) ? 0.0f : err;
}

// ---------- initialization ----------
void controller_init(float initial_roll_deg, float initial_pitch_deg) {
  err_roll_prev = 0.0f;
  err_pitch_prev = 0.0f;
  derr_roll_filt = 0.0f;
  derr_pitch_filt = 0.0f;
  ref_roll_prev = initial_roll_deg;
  ref_pitch_prev = initial_pitch_deg;
  i_roll = 0.0f;
  i_pitch = 0.0f;
}

void controller_reset(float ref_roll_deg_in, float ref_pitch_deg_in) {
  controller_init(ref_roll_deg_in, ref_pitch_deg_in);
}

// ---------- main update ----------
void controller_update(float desired_roll_deg,
                       float desired_pitch_deg,
                       float measured_roll_deg,
                       float measured_pitch_deg,
                       float dt_seconds,
                       float& ref_roll_deg_out,
                       float& ref_pitch_deg_out) {
  if (dt_seconds <= 0.0f) dt_seconds = 0.0005f;

  float err_roll  = desired_roll_deg  - measured_roll_deg;
  float err_pitch = desired_pitch_deg - measured_pitch_deg;

  err_roll  = apply_deadband(err_roll,  DEAD_BAND_DEG);
  err_pitch = apply_deadband(err_pitch, DEAD_BAND_DEG);

  float derr_roll  = (err_roll  - err_roll_prev)  / dt_seconds;
  float derr_pitch = (err_pitch - err_pitch_prev) / dt_seconds;

  derr_roll_filt  = DERIV_LPF_ALPHA * derr_roll_filt  + (1.0f - DERIV_LPF_ALPHA) * derr_roll;
  derr_pitch_filt = DERIV_LPF_ALPHA * derr_pitch_filt + (1.0f - DERIV_LPF_ALPHA) * derr_pitch;

  // --- integral term with anti-windup ---
  i_roll  += KI_ROLL  * err_roll  * dt_seconds;
  i_pitch += KI_PITCH * err_pitch * dt_seconds;
  if (i_roll  >  ICLAMP) i_roll  =  ICLAMP;
  if (i_roll  < -ICLAMP) i_roll  = -ICLAMP;
  if (i_pitch >  ICLAMP) i_pitch =  ICLAMP;
  if (i_pitch < -ICLAMP) i_pitch = -ICLAMP;

  // --- PID correction (absolute output, not accumulated) ---
  float corr_roll  = KP_ROLL  * err_roll  + KD_ROLL  * derr_roll_filt + i_roll;
  float corr_pitch = KP_PITCH * err_pitch + KD_PITCH * derr_pitch_filt + i_pitch;

  // absolute reference = command + correction
  float ref_roll  = desired_roll_deg  + corr_roll;
  float ref_pitch = desired_pitch_deg + corr_pitch;

  // float ref_roll  = corr_roll;
  // float ref_pitch = corr_pitch;

  // slew-rate limit
  // float max_step = MAX_REF_RATE_DEG_S * dt_seconds;
  // float dr = ref_roll  - ref_roll_prev;
  // float dp = ref_pitch - ref_pitch_prev;
  // if (dr >  max_step) dr =  max_step;
  // if (dr < -max_step) dr = -max_step;
  // if (dp >  max_step) dp =  max_step;
  // if (dp < -max_step) dp = -max_step;
  // ref_roll  = ref_roll_prev  + dr;
  // ref_pitch = ref_pitch_prev + dp;

  // clamp
  ref_roll  = constrain(ref_roll,  -REF_LIMIT_DEG, REF_LIMIT_DEG);
  ref_pitch = constrain(ref_pitch, -REF_LIMIT_DEG, REF_LIMIT_DEG);

  // output and store for next cycle
  ref_roll_deg_out  = ref_roll;
  ref_pitch_deg_out = ref_pitch;
  err_roll_prev     = err_roll;
  err_pitch_prev    = err_pitch;
  ref_roll_prev     = ref_roll;
  ref_pitch_prev    = ref_pitch;
}