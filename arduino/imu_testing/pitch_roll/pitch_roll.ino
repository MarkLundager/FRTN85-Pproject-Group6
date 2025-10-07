/*
   MPU9250 Complementary Filter (Pitch & Roll)
   Corrigido para evitar offset crescente
*/

#include <Wire.h>

const int MPU = 0x68; // I2C address (AD0=GND → 0x68)

// --- Variáveis globais ---
float AccX, AccY, AccZ;
float GyroX, GyroY, GyroZ;
float roll = 0.0, pitch = 0.0, yaw = 0.0;

// Calibração
float GyroErrorX = 0.0, GyroErrorY = 0.0, GyroErrorZ = 0.0;
const float rollOffset  = 0.257;
const float pitchOffset = -1.982;

// Filtro complementar
const float alpha = 0.96f; // peso do giroscópio
unsigned long t_prev;

// === Setup ===
void setup() {
  Serial.begin(19200);
  Wire.begin();

  // Wake up MPU
  Wire.beginTransmission(MPU);
  Wire.write(0x6B);
  Wire.write(0x00);
  Wire.endTransmission(true);

  // measureFlatOffsets();
  // while (true); // parar após medir


  delay(100);
  Serial.println("Calibrar IMU... mantenha imóvel");
  calibrateGyro();
  delay(100);
  Serial.println("Calibração completa.");

  t_prev = micros();
}

// === Loop principal ===
void loop() {
  // --- Ler acelerómetro ---
  Wire.beginTransmission(MPU);
  Wire.write(0x3B);
  Wire.endTransmission(false);
  Wire.requestFrom(MPU, 6, true);
  AccX = (Wire.read() << 8 | Wire.read()) / 16384.0;
  AccY = (Wire.read() << 8 | Wire.read()) / 16384.0;
  AccZ = (Wire.read() << 8 | Wire.read()) / 16384.0;

  float accAngleX = atan2(AccY, sqrt(AccX * AccX + AccZ * AccZ)) * 180.0 / PI - rollOffset;
  float accAngleY = atan2(-AccX, sqrt(AccY * AccY + AccZ * AccZ)) * 180.0 / PI - pitchOffset;

  // --- Ler giroscópio ---
  Wire.beginTransmission(MPU);
  Wire.write(0x43);
  Wire.endTransmission(false);
  Wire.requestFrom(MPU, 6, true);
  GyroX = (Wire.read() << 8 | Wire.read()) / 131.0;
  GyroY = (Wire.read() << 8 | Wire.read()) / 131.0;
  GyroZ = (Wire.read() << 8 | Wire.read()) / 131.0;

  // Remover bias
  GyroX -= GyroErrorX;
  GyroY -= GyroErrorY;
  GyroZ -= GyroErrorZ;

  // --- Calcular dt ---
  unsigned long t_now = micros();
  float dt = (t_now - t_prev) * 1e-6f;
  t_prev = t_now;
  if (dt <= 0) dt = 0.001f;

  // --- Proteger contra aceleração linear ---
  float anorm = sqrt(AccX * AccX + AccY * AccY + AccZ * AccZ);
  bool trustAcc = fabs(anorm - 1.0f) < 0.15f;
  float w = trustAcc ? (1.0f - alpha) : 0.0f;

  // --- Filtro complementar ---
  roll  = alpha * (roll  + GyroX * dt) + w * accAngleX;
  pitch = alpha * (pitch + GyroY * dt) + w * accAngleY;
  yaw  += GyroZ * dt; // drift, sem correção

  // --- Saída ---
  Serial.print(roll); Serial.print(" / ");
  Serial.println(pitch);

  delay(5); // ~200 Hz
}

// === Calibração ===
void calibrateGyro() {
  float gyroXSum = 0, gyroYSum = 0, gyroZSum = 0;
  const int samples = 200;

  Serial.println("A calibrar giroscópio... mantenha imóvel");

  for (int i = 0; i < samples; i++) {
    Wire.beginTransmission(MPU);
    Wire.write(0x43);
    Wire.endTransmission(false);
    Wire.requestFrom(MPU, 6, true);

    GyroX = (Wire.read() << 8 | Wire.read()) / 131.0;
    GyroY = (Wire.read() << 8 | Wire.read()) / 131.0;
    GyroZ = (Wire.read() << 8 | Wire.read()) / 131.0;

    gyroXSum += GyroX;
    gyroYSum += GyroY;
    gyroZSum += GyroZ;
    delay(5);
  }

  GyroErrorX = gyroXSum / samples;
  GyroErrorY = gyroYSum / samples;
  GyroErrorZ = gyroZSum / samples;

  Serial.println("Calibração do giroscópio completa:");
  Serial.print("Bias X: "); Serial.println(GyroErrorX);
  Serial.print("Bias Y: "); Serial.println(GyroErrorY);
  Serial.print("Bias Z: "); Serial.println(GyroErrorZ);
}


// === Medir offsets de roll/pitch em superfície plana ===
void measureFlatOffsets() {
  float rollSum = 0, pitchSum = 0;
  const int samples = 200;

  Serial.println("A medir offsets de roll/pitch... mantenha imóvel e nivelado");

  for (int i = 0; i < samples; i++) {
    Wire.beginTransmission(MPU);
    Wire.write(0x3B);
    Wire.endTransmission(false);
    Wire.requestFrom(MPU, 6, true);

    AccX = (Wire.read() << 8 | Wire.read()) / 16384.0;
    AccY = (Wire.read() << 8 | Wire.read()) / 16384.0;
    AccZ = (Wire.read() << 8 | Wire.read()) / 16384.0;

    float tempRoll  = atan2(AccY, sqrt(AccX * AccX + AccZ * AccZ)) * 180.0 / PI;
    float tempPitch = atan2(-AccX, sqrt(AccY * AccY + AccZ * AccZ)) * 180.0 / PI;

    rollSum  += tempRoll;
    pitchSum += tempPitch;
    delay(5);
  }

  float rollOffset  = rollSum / samples;
  float pitchOffset = pitchSum / samples;

  Serial.println("Offsets medidos (usar como constantes no código):");
  Serial.print("rollOffset = ");  Serial.println(rollOffset, 3);
  Serial.print("pitchOffset = "); Serial.println(pitchOffset, 3);
}

