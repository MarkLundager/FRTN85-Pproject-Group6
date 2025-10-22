// imu.cpp
// MPU9250 (MPU6050 register map) complementary filter for roll & pitch.
// Includes axis remap for your upside-down + swapped mounting
// and dynamic accel-offset calibration (current position = 0,0).

#include <Arduino.h>
#include <Wire.h>

static const uint8_t MPU = 0x68; // AD0=GND

// Raw sensors
static float AccX, AccY, AccZ, aux;
static float GyroX, GyroY, GyroZ;

// Estimates (deg)
static float roll_est = 0.0f, pitch_est = 0.0f, yaw_est = 0.0f;

// Gyro biases
static float GyroErrorX = 0.0f, GyroErrorY = 0.0f, GyroErrorZ = 0.0f;

// Accel offsets (deg)
static float rollOffset  = 0.0f;
static float pitchOffset = 0.0f;

// Complementary filter parameters
static const float alpha = 0.96f;
static unsigned long t_prev_us = 0;

// ---------- Gyro calibration ----------
static void calibrateGyro() {
  float sx=0, sy=0, sz=0;
  const int N = 200;

  Serial.println(F("Calibrating gyro... keep still"));
  for (int i = 0; i < N; i++) {
    Wire.beginTransmission(MPU);
    Wire.write(0x43);
    Wire.endTransmission(false);
    Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
    float gx = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;
    float gy = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;
    float gz = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;

    // Apply same axis remap as runtime
    aux = gx; gx = -gy; gy = aux; gz = -gz;

    sx += gx; sy += gy; sz += gz;
    delay(5);
  }
  GyroErrorX = sx / N;
  GyroErrorY = sy / N;
  GyroErrorZ = sz / N;

  Serial.print(F("Gyro bias X=")); Serial.println(GyroErrorX, 4);
  Serial.print(F("Gyro bias Y=")); Serial.println(GyroErrorY, 4);
  Serial.print(F("Gyro bias Z=")); Serial.println(GyroErrorZ, 4);
}

// ---------- Accelerometer offset calibration ----------
static void calibrateAccelOffsets() {
  const int N = 200;
  float rollSum = 0.0f, pitchSum = 0.0f;

  Serial.println(F("Calibrating accelerometer... keep platform level and still"));

  for (int i = 0; i < N; i++) {
    // Read accelerometer registers
    Wire.beginTransmission(MPU);
    Wire.write(0x3B);
    Wire.endTransmission(false);
    Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
    float ax = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
    float ay = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
    float az = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;

    // Apply same remap (IMU mounted upside-down)
    aux = ax;
    ax = -ay;
    ay =  aux;
    az = -az;

    // Compute instantaneous roll/pitch
    float rollAcc  = atan2f(ay, sqrtf(ax*ax + az*az)) * 180.0f/PI;
    float pitchAcc = atan2f(-ax, sqrtf(ay*ay + az*az)) * 180.0f/PI;

    rollSum  += rollAcc;
    pitchSum += pitchAcc;
    delay(5);
  }

  rollOffset  = rollSum / N;
  pitchOffset = pitchSum / N;

  Serial.println(F("Accelerometer offset calibration done."));
  Serial.print(F("rollOffset = "));  Serial.println(rollOffset, 3);
  Serial.print(F("pitchOffset = ")); Serial.println(pitchOffset, 3);
}

// ---------- Initialization ----------
void imu_begin() {
  Wire.begin();

  // Wake up
  Wire.beginTransmission(MPU);
  Wire.write(0x6B);
  Wire.write(0x00);
  Wire.endTransmission(true);

  delay(100);
  calibrateGyro();
  delay(200);
  calibrateAccelOffsets();  // <-- now also calibrate accel orientation

  t_prev_us = micros();
}

// ---------- Runtime update ----------
void imu_update() {
  // --- Read accelerometer ---
  Wire.beginTransmission(MPU);
  Wire.write(0x3B);
  Wire.endTransmission(false);
  Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
  AccX = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
  AccY = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
  AccZ = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;

  // Axis remap
  aux = AccX;
  AccX = -AccY;
  AccY =  aux;
  AccZ = -AccZ;

  float accRoll  = atan2f(AccY, sqrtf(AccX*AccX + AccZ*AccZ)) * 180.0f/PI - rollOffset;
  float accPitch = atan2f(-AccX, sqrtf(AccY*AccY + AccZ*AccZ)) * 180.0f/PI - pitchOffset;

  // --- Read gyro ---
  Wire.beginTransmission(MPU);
  Wire.write(0x43);
  Wire.endTransmission(false);
  Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
  GyroX = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;
  GyroY = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;
  GyroZ = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;

  // Axis remap same as accel
  aux = GyroX;
  GyroX = -GyroY;
  GyroY =  aux;
  GyroZ = -GyroZ;

  // Remove bias
  GyroX -= GyroErrorX;
  GyroY -= GyroErrorY;
  GyroZ -= GyroErrorZ;

  // dt
  unsigned long t_now = micros();
  float dt = (t_now - t_prev_us) * 1e-6f;
  t_prev_us = t_now;
  if (dt <= 0) dt = 0.001f;

  // trust accelerometer only when close to 1 g
  float anorm = sqrtf(AccX*AccX + AccY*AccY + AccZ*AccZ);
  float w = (fabsf(anorm - 1.0f) < 0.15f) ? (1.0f - alpha) : 0.0f;

  // Complementary filter (degrees)
  roll_est  = alpha * (roll_est  + GyroX * dt) + w * accRoll;
  pitch_est = alpha * (pitch_est + GyroY * dt) + w * accPitch;
  yaw_est  += GyroZ * dt;
}

// ---------- Access ----------
void imu_get_rp(float& roll_deg, float& pitch_deg) {
  roll_deg  = roll_est;
  pitch_deg = pitch_est;
}
