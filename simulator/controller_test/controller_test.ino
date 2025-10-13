#include <stdio.h>
#include <Arduino.h>
#include <Servo.h>

/* --------- Functions provided by your .cpp files (no headers) --------- */
extern void  imu_begin();
extern void  imu_update();
extern void  imu_get_rp(float& roll_deg, float& pitch_deg);

extern void  ik_init_geometry(); // builds geometry (B, pTop, beta, constants)
extern void  ik_compute(float roll_deg, float pitch_deg, float alpha_rad_out[6]);

extern void  fk_pose_from_alpha(const float alpha_rad[6], float& roll_deg_out, float& pitch_deg_out);

/* ------------------- Servo setup ------------------- */
static const uint8_t SERVO_PINS[6]    = {8, 9, 10, 11, 12, 13};
static const bool    SERVO_INVERTED[6]= {false, true, false, true, false, true};

// Map α[deg] to servo [0..180] (tune for your linkage)
static inline int mapAlphaToServoDeg(float alpha_deg) {
  long m = map((long)alpha_deg, -30, 30, 60, 120);
  if (m < 0) m = 0; if (m > 180) m = 180;
  return (int)m;
}

// Sign calibration between world roll/pitch and your IK convention
static const int ROLL_DIR  = -1;
static const int PITCH_DIR = -1;

static inline float clampf(float x, float lo, float hi){ return x<lo?lo:(x>hi?hi:x); }

Servo servos[6];

/* ------------------- Controller params ------------------- */
static const float IMU_TRIGGER_DEG = 0.6f;         // fire when |ΔIMU| exceeds this
static const uint32_t TRIGGER_MIN_MS = 40;         // rate limit
static const float CMD_LIMIT_DEG    = 8.0f;        // cap each corrective command
static const float K_CORR           = 1.0f;        // 0.6..1.0 typical (1=full cancel)

/* ------------------- State ------------------- */
static float last_alpha_rad[6] = {0};
static float last_fk_roll_deg=0, last_fk_pitch_deg=0;
static float last_imu_roll_deg=0, last_imu_pitch_deg=0;
static uint32_t last_trig_ms=0;

static bool manual_mode=false;
static float cmd_roll_deg=0, cmd_pitch_deg=0;

/* ------------------- Helpers ------------------- */
static void writeServosFromAlpha(const float alpha_rad[6]) {
  for (int i=0;i<6;i++) {
    float a_deg = alpha_rad[i]*180.0f/PI;
    int s = mapAlphaToServoDeg(a_deg);
    if (SERVO_INVERTED[i]) s = 180 - s;
    servos[i].write(s);
  }
}

static void parseSerialCommands() {
  if (!Serial.available()) return;
  String line = Serial.readStringUntil('\n'); line.trim();
  if (!line.length()) return;

  if (line == F("auto"))   { manual_mode=false; Serial.println(F("[MODE] auto"));   return; }
  if (line == F("manual")) { manual_mode=true;  Serial.println(F("[MODE] manual")); return; }

  float r=cmd_roll_deg, p=cmd_pitch_deg; int matched=0;
  if (line.indexOf('=')!=-1) {
    int ir=line.indexOf('r'); if(ir==-1) ir=line.indexOf('R');
    int ip=line.indexOf('p'); if(ip==-1) ip=line.indexOf('P');
    if (ir!=-1){ int eq=line.indexOf('=',ir); if(eq!=-1){ r=line.substring(eq+1).toFloat(); matched++; } }
    if (ip!=-1){ int eq=line.indexOf('=',ip); if(eq!=-1){ p=line.substring(eq+1).toFloat(); matched++; } }
  } else {
    char buf[64]; line.toCharArray(buf,sizeof(buf));
    if (sscanf(buf, "%f %f", &r, &p)==2) matched=2;
  }
  if (matched) {
    cmd_roll_deg=r; cmd_pitch_deg=p; manual_mode=true;
    Serial.print(F("[CMD manual] r=")); Serial.print(cmd_roll_deg,2);
    Serial.print(F(" p="));            Serial.println(cmd_pitch_deg,2);
  }
}

/* ------------------- Setup ------------------- */
void setup() {
  Serial.begin(115200);
  Serial.println(F("\n=== Stewart Platform Controller — Event-Triggered ==="));

  ik_init_geometry();

  for (int i=0;i<6;i++) servos[i].attach(SERVO_PINS[i]);

  // Neutral
  float alpha0[6]; ik_compute(0,0,alpha0);
  writeServosFromAlpha(alpha0);
  for(int i=0;i<6;i++) last_alpha_rad[i]=alpha0[i];

  // FK baseline from neutral
  fk_pose_from_alpha(last_alpha_rad, last_fk_roll_deg, last_fk_pitch_deg);

  // IMU
  imu_begin();
  last_imu_roll_deg=0; last_imu_pitch_deg=0;

  Serial.println(F("Ready. Send 'auto'/'manual' or 'r=.. p=..' / '<r> <p>'"));
}

/* ------------------- Loop ------------------- */
void loop() {
  // 1) IMU update
  imu_update();
  float imu_r=0, imu_p=0; imu_get_rp(imu_r, imu_p);

  // 2) Commands
  parseSerialCommands();

  // 3) AUTO: event-triggered corrections
  if (!manual_mode) {
    float d_r = imu_r - last_imu_roll_deg;
    float d_p = imu_p - last_imu_pitch_deg;
    bool trig = (fabs(d_r) >= IMU_TRIGGER_DEG) || (fabs(d_p) >= IMU_TRIGGER_DEG);
    uint32_t now=millis();

    if (trig && (now - last_trig_ms >= TRIGGER_MIN_MS)) {
      last_trig_ms = now;

      // FK from current alphas → absolute baseline
      fk_pose_from_alpha(last_alpha_rad, last_fk_roll_deg, last_fk_pitch_deg);

      // Total world tilt estimate = FK + ΔIMU
      float tot_r = last_fk_roll_deg  + d_r;
      float tot_p = last_fk_pitch_deg + d_p;

      // Opposite correction (clamped)
      float corr_r = clampf(-K_CORR*tot_r, -CMD_LIMIT_DEG, +CMD_LIMIT_DEG);
      float corr_p = clampf(-K_CORR*tot_p, -CMD_LIMIT_DEG, +CMD_LIMIT_DEG);

      // IK for the corrective command (with sign calibration)
      float alpha_new[6];
      ik_compute(ROLL_DIR*corr_r, PITCH_DIR*corr_p, alpha_new);
      writeServosFromAlpha(alpha_new);
      for(int i=0;i<6;i++) last_alpha_rad[i]=alpha_new[i];

      // reset IMU delta baseline
      last_imu_roll_deg  = imu_r;
      last_imu_pitch_deg = imu_p;

      Serial.print(F("TRIG ΔIMU=("));
      Serial.print(d_r,2); Serial.print(',');
      Serial.print(d_p,2); Serial.print(F(") FK=("));
      Serial.print(last_fk_roll_deg,2); Serial.print(',');
      Serial.print(last_fk_pitch_deg,2); Serial.print(F(") CMD=("));
      Serial.print(corr_r,2); Serial.print(',');
      Serial.print(corr_p,2); Serial.println(')');
    }
  } else {
    // MANUAL: keep applying current cmd
    float alpha[6];
    ik_compute(ROLL_DIR*cmd_roll_deg, PITCH_DIR*cmd_pitch_deg, alpha);
    writeServosFromAlpha(alpha);
    for(int i=0;i<6;i++) last_alpha_rad[i]=alpha[i];
  }

  // 4) slow telemetry
  static uint32_t t0=0;
  if (millis()-t0>120){
    t0=millis();
    Serial.print(F("IMU r=")); Serial.print(imu_r,2);
    Serial.print(F(" p="));     Serial.print(imu_p,2);
    Serial.print(F(" | mode="));Serial.println(manual_mode?F("MAN"):F("AUTO"));
  }

  delay(5);
}