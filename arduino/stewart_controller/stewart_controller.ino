// stewart_controller.ino
// Main: reads IMU roll/pitch, reads commanded roll/pitch from Serial,
// runs IK, maps to servo angles, writes servos, prints telemetry.

#include <Arduino.h>
#include <Servo.h>
#include <math.h>

#include "controller.h"

// --------- Externs from imu.cpp ---------
extern void  imu_begin();
extern void  imu_update();
extern void  imu_get_rp(float& roll_deg, float& pitch_deg);

// --------- Externs from inverse_kinematics.cpp ---------
extern void  ik_init_geometry(); // precomputes B, p, beta
extern void  ik_compute(float roll_deg, float pitch_deg, float alpha_rad_out[6]);

// ---------- Servo setup ----------
static const int SERVO_PINS[6] = {8, 9, 10, 11, 12, 13};
// Left-mounted servos are inverted (per your hardware note)
static const bool SERVO_INVERTED[6] = {false, true, false, true, false, true};

// Map α [deg] to servo command [0..180]; center ~90, ±30° span → [60..120]
static inline int mapAlphaToServoDeg(float alpha_deg) {
  long m = map((long)alpha_deg, -30, 30, 60, 120);
  if (m < 0) m = 0; if (m > 180) m = 180;
  return (int)m;
}

// Optional overall sign calibration (keep as in your last working build)
static const int ROLL_DIR  = -1;  // +roll tilts +Y side down
static const int PITCH_DIR = -1;  // +pitch tilts +X side down

Servo servos[6];

// Commanded pose from Serial (deg)
volatile float cmd_roll_deg  = 0.0f;
volatile float cmd_pitch_deg = 0.0f;
static unsigned long last_loop_us = 0;

static void parseSerialCommands() {
  // Accept lines like: "r=3.5 p=-2"  or "3.5 -2"
  if (!Serial.available()) return;
  String line = Serial.readStringUntil('\n');
  line.trim();
  if (line.length() == 0) return;

  float r=cmd_roll_deg, p=cmd_pitch_deg;
  int matched = 0;

  // Try "r=.. p=.."
  if (line.indexOf('r') != -1 || line.indexOf('R') != -1) {
    // crude parse
    int ir = line.indexOf('r');
    if (ir == -1) ir = line.indexOf('R');
    int ip = line.indexOf('p');
    if (ip == -1) ip = line.indexOf('P');
    if (ir != -1) {
      int eq = line.indexOf('=', ir);
      if (eq != -1) r = line.substring(eq+1).toFloat(), matched++;
    }
    if (ip != -1) {
      int eq = line.indexOf('=', ip);
      if (eq != -1) p = line.substring(eq+1).toFloat(), matched++;
    }
  } else {
    // Try two floats
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

void setup() {
  Serial.begin(115200);
  Serial.println(F("\n=== Stewart Platform Controller Boot ==="));

  // 1) Initialize geometry
  ik_init_geometry();

  // 2) Attach servos and move to neutral position (flat platform)
  for (int i = 0; i < 6; i++) {
    servos[i].attach(SERVO_PINS[i]);
  }

  // Compute neutral alpha and send to servos
  float alpha_rad[6];
  ik_compute(0.0f, 0.0f, alpha_rad);
  for (int i = 0; i < 6; i++) {
    float deg = alpha_rad[i] * 180.0f / PI;
    int mapped = map((int)deg, -30, 30, 60, 120);
    if (SERVO_INVERTED[i]) mapped = 180 - mapped;
    mapped = constrain(mapped, 0, 180);
    servos[i].write(mapped);
  }

  Serial.println(F("Servos set to neutral (roll=0, pitch=0)"));
  delay(1500); // wait for the platform to settle

  // 3) Initialize IMU (do this while platform is level and still)
  imu_begin();

  controller_init(cmd_roll_deg, cmd_pitch_deg);
  last_loop_us = micros();

  Serial.println(F("Stewart controller ready. Send 'r=<deg> p=<deg>' or '<r> <p>'"));
}

void loop() {
  unsigned long now_us = micros();
  float dt = (now_us - last_loop_us) * 1.0e-6f;
  if (dt <= 0.0f) {
    dt = 0.0005f; // fallback to 2 kHz in rare micros rollover race
  }
  last_loop_us = now_us;

  // 1) Update IMU fusion
  imu_update();
  float meas_roll_deg=0, meas_pitch_deg=0;
  imu_get_rp(meas_roll_deg, meas_pitch_deg);

  // 2) Handle serial commands
  parseSerialCommands();

  // Reset controller reference if a new command is provided (prevents windup)
  static float last_cmd_roll = 0.0f;
  static float last_cmd_pitch = 0.0f;
  if (fabsf(cmd_roll_deg - last_cmd_roll) > 1e-3f ||
      fabsf(cmd_pitch_deg - last_cmd_pitch) > 1e-3f) {
    controller_reset(cmd_roll_deg, cmd_pitch_deg);
    last_cmd_roll  = cmd_roll_deg;
    last_cmd_pitch = cmd_pitch_deg;
  }

  // 3) PD attitude controller
  float ref_roll_deg = 0.0f;
  float ref_pitch_deg = 0.0f;
  controller_update(cmd_roll_deg, cmd_pitch_deg,
                    meas_roll_deg, meas_pitch_deg,
                    dt,
                    ref_roll_deg, ref_pitch_deg);

  // 4) Compute IK for commanded pose (apply sign calibration)
  float alpha_rad[6];
  ik_compute(ROLL_DIR*ref_roll_deg, PITCH_DIR*ref_pitch_deg, alpha_rad);

  // 5) Map α → servo and write
  for (int i = 0; i < 6; i++) {
    float alpha_deg = alpha_rad[i] * 180.0f / PI;
    int sdeg = mapAlphaToServoDeg(alpha_deg);
    if (SERVO_INVERTED[i]) sdeg = 180 - sdeg;
    servos[i].write(sdeg);
  }

  // 6) Telemetry at ~10 Hz
  static unsigned long lastPrint=0;
  if (millis() - lastPrint > 100) {
    lastPrint = millis();
    Serial.print(F("IMU r="));   Serial.print(meas_roll_deg, 2);
    Serial.print(F(" p="));      Serial.print(meas_pitch_deg, 2);
    Serial.print(F(" | CMD r="));Serial.print(cmd_roll_deg, 2);
    Serial.print(F(" p="));      Serial.print(cmd_pitch_deg, 2);
    Serial.print(F(" | REF r="));Serial.print(ref_roll_deg, 2);
    Serial.print(F(" p="));      Serial.println(ref_pitch_deg, 2);
  }
}
