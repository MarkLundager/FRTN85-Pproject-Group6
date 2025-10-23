// stewart_controller.ino
// Controller + IK + FK + IMU, with proper headers only.

#include <Arduino.h>
#include "ik.h"
#include <Servo.h>
#include <math.h>
#include "imu.h"
#include "controller.h"
#include "geometry.h"
#include "fk.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



// ---------- Servo setup ----------
static const int SERVO_PINS[6] = {8, 9, 10, 11, 12, 13};
// Left-mounted servos are inverted (per your hardware note)
static const bool SERVO_INVERTED[6] = {false, true, false, true, false, true};

// Optional overall sign calibration (keep as in your last working build)
static const int ROLL_DIR  = -1;  // +roll tilts +Y side down
static const int PITCH_DIR = -1;  // +pitch tilts +X side down

Servo servos[6];

// Commanded pose from Serial (deg)
volatile float cmd_roll_deg  = 0.0f;
volatile float cmd_pitch_deg = 0.0f;

// --- Fusion weight: IMU primary, FK secondary ---
static const float W_IMU = 0.75f;   // 0.6..0.9 typical

// --- Keep last commanded servo angles (deg) for FK each cycle ---
static float last_alpha_deg[6] = {0,0,0,0,0,0};
static float roll_fk_deg_prev = 0.0f, pitch_fk_deg_prev = 0.0f;

// --- Loop timing ---
static unsigned long last_loop_us = 0;

// Map α [deg] to servo command [0..180]; center ~90, ±30° span → [60..120]
static inline int mapAlphaToServoDeg(float alpha_deg) {
  long m = map((long)alpha_deg, -30, 30, 60, 120);
  if (m < 0) m = 0; if (m > 180) m = 180;
  return (int)m;
}

static void parseSerialCommands() {
  if (!Serial.available()) return;
  String line = Serial.readStringUntil('\n');
  line.trim();
  if (line.length() == 0) return;

  float r=cmd_roll_deg, p=cmd_pitch_deg;
  int matched = 0;

  if (line.indexOf('r') != -1 || line.indexOf('R') != -1) {
    int ir = line.indexOf('r'); if (ir == -1) ir = line.indexOf('R');
    int ip = line.indexOf('p'); if (ip == -1) ip = line.indexOf('P');
    if (ir != -1) { int eq = line.indexOf('=', ir); if (eq != -1) { r = line.substring(eq+1).toFloat(); matched++; } }
    if (ip != -1) { int eq = line.indexOf('=', ip); if (eq != -1) { p = line.substring(eq+1).toFloat(); matched++; } }
  } else {
    char buf[64];
    line.toCharArray(buf, sizeof(buf));
    if (sscanf(buf, "%f %f", &r, &p) == 2) matched = 2;
  }

  if (matched) {
    cmd_roll_deg  = r;
    cmd_pitch_deg = p;
    Serial.print(F("[CMD] roll="));  Serial.print(cmd_roll_deg, 2);
    Serial.print(F(" pitch="));      Serial.println(cmd_pitch_deg, 2);
  }
}

// Build horn tip positions from horn angles (degrees) for FK
static inline void build_H_from_alpha_deg(const float alpha_deg[6], float H[6][3]) {
  for (int i = 0; i < 6; i++) {
    float th = alpha_deg[i] * (float)M_PI / 180.0f;
    H[i][0] = B[i][0] + H_SERVO * cosf(beta[i]) * cosf(th);
    H[i][1] = B[i][1] + H_SERVO * sinf(beta[i]) * cosf(th);
    H[i][2] =             H_SERVO *            sinf(th);
  }
}

void setup() {
  Serial.begin(115200);
  Serial.println(F("\n=== Stewart Platform Controller Boot (Controller + IK + FK + IMU) ==="));

  // Initialize geometry
  ik_init_geometry();   // your IK precompute (radians API)
  geometry_init();      // shared geometry for FK (B, p, beta, Z0, etc.)

  // Attach servos and move to neutral (flat platform)
  for (int i = 0; i < 6; i++) servos[i].attach(SERVO_PINS[i]);

  float alpha_rad_neutral[6];
  ik_compute(0.0f, 0.0f, alpha_rad_neutral);
  for (int i = 0; i < 6; i++) {
    float deg = alpha_rad_neutral[i] * 180.0f / (float)M_PI;
    last_alpha_deg[i] = deg; // seed FK
    int mapped = mapAlphaToServoDeg(deg);
    if (SERVO_INVERTED[i]) mapped = 180 - mapped;
    servos[i].write(constrain(mapped, 0, 180));
  }

  Serial.println(F("Servos set to neutral (roll=0, pitch=0)"));
  delay(1500); // wait to settle

  // Initialize IMU while level & still
  imu_begin();

  // Seed initial FK estimate (optional)
  {
    float H0[6][3]; build_H_from_alpha_deg(last_alpha_deg, H0);
    float R0[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
    float T0[3] = {0.0f, 0.0f, Z0};
    float R[3][3], T[3];
    fk_iterative_6dof(H0, p, D_ROD, R0, T0, 30, 0.10f, 1e-5f, 1e-5f, R, T);
    float r,pch; extract_roll_pitch_from_R(R, r, pch);
    roll_fk_deg_prev  = r   * 180.0f/(float)M_PI;
    pitch_fk_deg_prev = pch * 180.0f/(float)M_PI;
  }

  controller_init(cmd_roll_deg, cmd_pitch_deg);
  last_loop_us = micros();

  Serial.println(F("Ready. Send 'r=<deg> p=<deg>' or '<r> <p>'"));
}

void loop() {
  unsigned long now_us = micros();
  float dt = (now_us - last_loop_us) * 1.0e-6f;
  if (dt <= 0.0f) dt = 0.0005f; // fallback
  last_loop_us = now_us;

  // 1) Update IMU
  imu_update();
  float meas_roll_deg=0.0f, meas_pitch_deg=0.0f;
  imu_get_rp(meas_roll_deg, meas_pitch_deg);

  // 2) Handle serial commands (updates cmd_roll_deg / cmd_pitch_deg)
  parseSerialCommands();

  // 3) Fuse IMU with last-cycle FK as measured tilt for the controller
  float roll_meas_fused  = W_IMU * meas_roll_deg  + (1.0f - W_IMU) * roll_fk_deg_prev;
  float pitch_meas_fused = W_IMU * meas_pitch_deg + (1.0f - W_IMU) * pitch_fk_deg_prev;

  // 4) Run your controller (outputs absolute roll/pitch reference in deg)
  float ref_roll_deg = 0.0f, ref_pitch_deg = 0.0f;

  // Reset controller if a new command arrives (anti-windup style)
  static float last_cmd_roll = 0.0f, last_cmd_pitch = 0.0f;
  if (fabsf(cmd_roll_deg - last_cmd_roll) > 1e-3f ||
      fabsf(cmd_pitch_deg - last_cmd_pitch) > 1e-3f) {
    controller_reset(cmd_roll_deg, cmd_pitch_deg);
    last_cmd_roll  = cmd_roll_deg;
    last_cmd_pitch = cmd_pitch_deg;
  }

  controller_update(cmd_roll_deg, cmd_pitch_deg,
                    roll_meas_fused, pitch_meas_fused,
                    dt,
                    ref_roll_deg, ref_pitch_deg);

  // 5) IK → servo angles (radians out)
  float alpha_rad_cmd[6];
  ik_compute(ROLL_DIR*ref_roll_deg, PITCH_DIR*ref_pitch_deg, alpha_rad_cmd);

  // --- servo write rate limit to ~100 Hz ---
  static unsigned long lastServoWrite = 0;
  const unsigned long SERVO_WRITE_PERIOD_US = 10000; // 10 ms
  bool doWrite = (now_us - lastServoWrite) >= SERVO_WRITE_PERIOD_US;
  if (doWrite) lastServoWrite = now_us;

  if (doWrite) {
    for (int i = 0; i < 6; i++) {
      float alpha_deg = alpha_rad_cmd[i] * 180.0f / (float)M_PI;
      last_alpha_deg[i] = alpha_deg; // store for FK
      int sdeg = mapAlphaToServoDeg(alpha_deg);
      if (SERVO_INVERTED[i]) sdeg = 180 - sdeg;
      servos[i].write(sdeg);
    }
  }

  // 6) Recompute FK from last_alpha_deg for NEXT cycle’s fusion + Telemetry
  {
    float H[6][3]; build_H_from_alpha_deg(last_alpha_deg, H);
    float R_init[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
    float T_init[3] = {0.0f, 0.0f, Z0};
    float R_out[3][3], T_out[3];

    fk_iterative_6dof(H, p, D_ROD, R_init, T_init, 30, 0.10f, 1e-5f, 1e-5f, R_out, T_out);

    float roll_fk_rad, pitch_fk_rad;
    extract_roll_pitch_from_R(R_out, roll_fk_rad, pitch_fk_rad);
    roll_fk_deg_prev  = roll_fk_rad  * 180.0f / (float)M_PI;
    pitch_fk_deg_prev = pitch_fk_rad * 180.0f / (float)M_PI;

    // Telemetry at ~10 Hz
    static unsigned long lastPrint=0;
    if (millis() - lastPrint > 100) {
      lastPrint = millis();
      Serial.print(F("IMU r="));   Serial.print(meas_roll_deg, 2);
      Serial.print(F(" p="));      Serial.print(meas_pitch_deg, 2);
      Serial.print(F(" | FUSED r=")); Serial.print(roll_meas_fused, 2);
      Serial.print(F(" p="));        Serial.print(pitch_meas_fused, 2);
      Serial.print(F(" | REF r="));  Serial.print(ref_roll_deg, 2);
      Serial.print(F(" p="));        Serial.print(ref_pitch_deg, 2);
      Serial.print(F(" | FK r="));   Serial.print(roll_fk_deg_prev, 2);
      Serial.print(F(" p="));        Serial.print(pitch_fk_deg_prev, 2);
      Serial.print(F(" Z="));        Serial.println(T_out[2], 2);
    }
  }
}