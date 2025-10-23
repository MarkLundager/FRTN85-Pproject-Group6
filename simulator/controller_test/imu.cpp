




#include <Arduino.h>
#include <Wire.h>

static const uint8_t MPU = 0x68; 


static float AccX, AccY, AccZ, aux;
static float GyroX, GyroY, GyroZ;


static float roll_est = 0.0f, pitch_est = 0.0f, yaw_est = 0.0f;


static float GyroErrorX = 0.0f, GyroErrorY = 0.0f, GyroErrorZ = 0.0f;


static float rollOffset  = 0.0f;
static float pitchOffset = 0.0f;


static const float alpha = 0.96f;
static unsigned long t_prev_us = 0;


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


static void calibrateAccelOffsets() {
  const int N = 200;
  float rollSum = 0.0f, pitchSum = 0.0f;

  Serial.println(F("Calibrating accelerometer... keep platform level and still"));

  for (int i = 0; i < N; i++) {
    
    Wire.beginTransmission(MPU);
    Wire.write(0x3B);
    Wire.endTransmission(false);
    Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
    float ax = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
    float ay = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
    float az = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;

    
    aux = ax;
    ax = -ay;
    ay =  aux;
    az = -az;

    
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


void imu_begin() {
  Wire.begin();

  
  Wire.beginTransmission(MPU);
  Wire.write(0x6B);
  Wire.write(0x00);
  Wire.endTransmission(true);

  delay(100);
  calibrateGyro();
  delay(200);
  calibrateAccelOffsets();  

  t_prev_us = micros();
}


void imu_update() {
  
  Wire.beginTransmission(MPU);
  Wire.write(0x3B);
  Wire.endTransmission(false);
  Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
  AccX = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
  AccY = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;
  AccZ = (int16_t)((Wire.read()<<8) | Wire.read()) / 16384.0f;

  
  aux = AccX;
  AccX = -AccY;
  AccY =  aux;
  AccZ = -AccZ;

  float accRoll  = atan2f(AccY, sqrtf(AccX*AccX + AccZ*AccZ)) * 180.0f/PI - rollOffset;
  float accPitch = atan2f(-AccX, sqrtf(AccY*AccY + AccZ*AccZ)) * 180.0f/PI - pitchOffset;

  
  Wire.beginTransmission(MPU);
  Wire.write(0x43);
  Wire.endTransmission(false);
  Wire.requestFrom((uint8_t)MPU, (uint8_t)6, (uint8_t)true);
  GyroX = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;
  GyroY = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;
  GyroZ = (int16_t)((Wire.read()<<8) | Wire.read()) / 131.0f;

  
  aux = GyroX;
  GyroX = -GyroY;
  GyroY =  aux;
  GyroZ = -GyroZ;

  
  GyroX -= GyroErrorX;
  GyroY -= GyroErrorY;
  GyroZ -= GyroErrorZ;

  
  unsigned long t_now = micros();
  float dt = (t_now - t_prev_us) * 1e-6f;
  t_prev_us = t_now;
  if (dt <= 0) dt = 0.001f;

  
  float anorm = sqrtf(AccX*AccX + AccY*AccY + AccZ*AccZ);
  float w = (fabsf(anorm - 1.0f) < 0.15f) ? (1.0f - alpha) : 0.0f;

  
  roll_est  = alpha * (roll_est  + GyroX * dt) + w * accRoll;
  pitch_est = alpha * (pitch_est + GyroY * dt) + w * accPitch;
  yaw_est  += GyroZ * dt;
}


void imu_get_rp(float& roll_deg, float& pitch_deg) {
  roll_deg  = roll_est;
  pitch_deg = pitch_est;
}
