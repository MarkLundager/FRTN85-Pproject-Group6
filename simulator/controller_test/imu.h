#pragma once

void imu_begin();
void imu_update();
void imu_get_rp(float& roll_deg, float& pitch_deg);