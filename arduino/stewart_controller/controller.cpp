#include <Arduino.h>
#include <math.h>

#include "controller.h"

namespace {
// ---------- PD attitude controller configuration ----------
static const float KP_ROLL  = 0.3f;   // proportional gain [deg/deg]
static const float KP_PITCH = 0.3f;
static const float KD_ROLL  = 0.05f;  // derivative gain [deg*s/deg]
static const float KD_PITCH = 0.05f;
static const float DERIV_LPF_ALPHA = 0.8f; // 0..1, closer to 1 = more smoothing
static const float REF_LIMIT_DEG    = 12.0f; // clamp controller output to safe range

static float ref_roll_deg  = 0.0f;
static float ref_pitch_deg = 0.0f;
static float err_roll_prev  = 0.0f;
static float err_pitch_prev = 0.0f;
static float derr_roll_filt  = 0.0f;
static float derr_pitch_filt = 0.0f;
}

void controller_init(float initial_roll_deg, float initial_pitch_deg) {
  ref_roll_deg  = initial_roll_deg;
  ref_pitch_deg = initial_pitch_deg;
  err_roll_prev  = 0.0f;
  err_pitch_prev = 0.0f;
  derr_roll_filt  = 0.0f;
  derr_pitch_filt = 0.0f;
}

void controller_reset(float ref_roll_deg_in, float ref_pitch_deg_in) {
  controller_init(ref_roll_deg_in, ref_pitch_deg_in);
}
static const float DEAD_BAND_DEG = 1.0f;  // ignore errors smaller than this

inline float apply_deadband(float err, float threshold) {
  return (fabs(err) < threshold) ? 0.0f : err;
}

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

  // suppress small oscillations
  err_roll  = apply_deadband(err_roll,  DEAD_BAND_DEG);
  err_pitch = apply_deadband(err_pitch, DEAD_BAND_DEG);

  float derr_roll  = (err_roll  - err_roll_prev)  / dt_seconds;
  float derr_pitch = (err_pitch - err_pitch_prev) / dt_seconds;

  derr_roll_filt  = DERIV_LPF_ALPHA * derr_roll_filt  + (1.0f - DERIV_LPF_ALPHA) * derr_roll;
  derr_pitch_filt = DERIV_LPF_ALPHA * derr_pitch_filt + (1.0f - DERIV_LPF_ALPHA) * derr_pitch;

  float ctrl_roll  = KP_ROLL  * err_roll  + KD_ROLL  * derr_roll_filt;
  float ctrl_pitch = KP_PITCH * err_pitch + KD_PITCH * derr_pitch_filt;

  ref_roll_deg  += ctrl_roll;
  ref_pitch_deg += ctrl_pitch;

  ref_roll_deg  = constrain(ref_roll_deg,  -REF_LIMIT_DEG, REF_LIMIT_DEG);
  ref_pitch_deg = constrain(ref_pitch_deg, -REF_LIMIT_DEG, REF_LIMIT_DEG);

  err_roll_prev  = err_roll;
  err_pitch_prev = err_pitch;

  ref_roll_deg_out  = ref_roll_deg;
  ref_pitch_deg_out = ref_pitch_deg;
}
