#include <Servo.h>
#include <math.h>

// ---------- Geometry constants (cm) ----------
const float h = 1.64;
const float d = 17.0;
const float r_top = 7.8;
const float r_bot = 9.3;
const float z0 = 16.3;

// ---------- Layout angles (rad) ----------
const float pDelta = 15.25 * M_PI / 180.0;
const float bDelta = 23.58 * M_PI / 180.0;

// ---------- Servo objects ----------
Servo servos[6];
const int servoPins[6] = {8, 9, 10, 11, 12, 13};

// Left-mounted servos are inverted
bool servoInverted[6] = {false, true, false, true, false, true};

// ---------- Sign calibration ----------
const int ROLL_DIR  = -1;   // set to -1 if roll moves opposite
const int PITCH_DIR = -1;   // set to +1 if pitch moves opposite

// ---------- Platform geometry arrays ----------
float beta[6], b_phi[6], p_phi[6];
float Bx[6], By[6];
float px[6], py[6];

// ---------- Setup ----------
void setup() {
  Serial.begin(115200);
  for (int i = 0; i < 6; i++) servos[i].attach(servoPins[i]);

  for (int k = 0; k < 6; k++) {
    p_phi[k] = (2*M_PI/3)*floor(k/2.0) - pow(-1,k)*(pDelta/2.0) + M_PI/3;
    b_phi[k] = (2*M_PI/3)*floor((k+1)/2.0) + pow(-1,k)*(bDelta/2.0);
    beta[k]  = b_phi[k] + (M_PI/2.0)*pow(-1,k);

    Bx[k] = r_bot * cos(b_phi[k]);
    By[k] = r_bot * sin(b_phi[k]);
    px[k] = r_top * cos(p_phi[k]);
    py[k] = r_top * sin(p_phi[k]);
  }

  Serial.println("Stewart IK (with roll/pitch sign flags and mirrored servos)");
}

// ---------- Rotation helpers ----------
void rotXY(float roll, float pitch, float px[], float py[], float pz[], float outx[], float outy[], float outz[]) {
  float cr = cos(roll), sr = sin(roll);
  float cp = cos(pitch), sp = sin(pitch);
  for (int i = 0; i < 6; i++) {
    float x = px[i], y = py[i], z = pz[i];
    // R = Ry(pitch)*Rx(roll)
    float x1 = cp*x + sp*z;
    float y1 = sr*sp*x + cr*y - sr*cp*z;
    float z1 = -cr*sp*x + sr*y + cr*cp*z;
    outx[i] = x1; outy[i] = y1; outz[i] = z1 + z0;
  }
}

// ---------- IK ----------
void computeIK(float roll, float pitch, float alpha[6]) {
  float pz[6] = {0,0,0,0,0,0};
  float Px[6], Py[6], Pz[6];
  rotXY(roll, pitch, px, py, pz, Px, Py, Pz);

  for (int i = 0; i < 6; i++) {
    float lx = Px[i] - Bx[i];
    float ly = Py[i] - By[i];
    float lz = Pz[i];

    float e = 2*h*lz;
    float f = 2*h*(cos(beta[i])*lx + sin(beta[i])*ly);
    float g = lx*lx + ly*ly + lz*lz - (d*d - h*h);
    float A = sqrt(e*e + f*f);

    float val = g / A;
    if (val > 1.0) val = 1.0;
    if (val < -1.0) val = -1.0;

    alpha[i] = asin(val) - atan2(f, e);
  }
}

// ---------- Loop ----------
void loop() {
  // Commanded pose (change these)
  float roll_cmd_deg  = 0.0;
  float pitch_cmd_deg = 0.0;

  float roll  = ROLL_DIR  * roll_cmd_deg  * M_PI / 180.0f;
  float pitch = PITCH_DIR * pitch_cmd_deg * M_PI / 180.0f;

  float alpha[6];
  computeIK(roll, pitch, alpha);

  int pulse[6];
  for (int i = 0; i < 6; i++) {
    float deg = alpha[i] * 180.0f / M_PI;
    int mapped = map((int)deg, -30, 30, 60, 120);   // center around ~90
    if (servoInverted[i]) mapped = 180 - mapped;    // mirror left-mounted
    pulse[i] = constrain(mapped, 0, 180);
    servos[i].write(pulse[i]);
  }

  static unsigned long lastPrint = 0;
  if (millis() - lastPrint > 500) {
    lastPrint = millis();
    Serial.print("cmd roll=");
    Serial.print(roll_cmd_deg, 1);
    Serial.print("  cmd pitch=");
    Serial.print(pitch_cmd_deg, 1);
    Serial.print("  eff roll=");
    Serial.print(roll * 180.0 / M_PI, 2);
    Serial.print("  eff pitch=");
    Serial.print(pitch * 180.0 / M_PI, 2);
    Serial.print("  alpha(deg)=");
    for (int i = 0; i < 6; i++) { Serial.print(alpha[i] * 180.0 / M_PI, 1); Serial.print(" "); }
    Serial.print("  servo=");
    for (int i = 0; i < 6; i++) { Serial.print(pulse[i]); Serial.print(" "); }
    Serial.println();
  }

  delay(20);
}
