#ifndef STEWART_CONTROLLER_CONTROLLER_H
#define STEWART_CONTROLLER_CONTROLLER_H

void controller_init(float initial_roll_deg, float initial_pitch_deg);
void controller_reset(float ref_roll_deg, float ref_pitch_deg);
void controller_update(float desired_roll_deg,
                       float desired_pitch_deg,
                       float measured_roll_deg,
                       float measured_pitch_deg,
                       float dt_seconds,
                       float& ref_roll_deg_out,
                       float& ref_pitch_deg_out);

#endif