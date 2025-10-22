#pragma once

// Radians API to match your .ino
// roll_deg / pitch_deg are in DEGREES; output alpha in RADIANS.
void ik_init_geometry();  // no-op (we use geometry.h), kept for compatibility
void ik_compute(float roll_deg, float pitch_deg, float alpha_rad_out[6]);