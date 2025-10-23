#include <Wire.h>

// --- DEVICE ADDRESSES ---
#define MPU9250_ADDR 0x68
#define AK8963_ADDR  0x0C

// --- MPU9250 REGISTERS ---
#define MPU_REG_SMPLRT_DIV 0x19
#define MPU_REG_CONFIG     0x1A
#define MPU_REG_GYRO_CFG   0x1B
#define MPU_REG_ACCEL_CFG  0x1C
#define MPU_REG_INT_PIN_CFG 0x37
#define MPU_REG_PWR_MGMT_1 0x6B
#define MPU_REG_ACCEL_XOUT_H 0x3B
#define MPU_REG_GYRO_XOUT_H  0x43

// --- AK8963 REGISTERS ---
#define AK8963_REG_ST1   0x02
#define AK8963_REG_HXL   0x03
#define AK8963_REG_ST2   0x09
#define AK8963_REG_CNTL1 0x0A  // 0x16 = 16-bit, cont. mode 2 (100 Hz)

// --- SCALES ---
#define ACCEL_SCALE 16384.0  // ±2 g
#define GYRO_SCALE  131.0    // ±250 °/s
#define MAG_SCALE   0.15     // µT/LSB (16-bit)

// --- OFFSETS ---
double accelOffsetX=0, accelOffsetY=0, accelOffsetZ=0;
double gyroXOffset=0,  gyroYOffset=0,  gyroZOffset=0;
double magMinX=0, magMaxX=0, magMinY=0, magMaxY=0, magMinZ=0, magMaxZ=0;
double magScaleX=1.0, magScaleY=1.0, magScaleZ=1.0;

// --- STATE ---
double pitch=0, roll=0, yaw=0;   // degrees
double dt=0.005;                 // seconds
unsigned long lastMicros;

// --- KALMAN (per axis) ---
struct Kalman1D {
  double angle=0, bias=0;
  double P00=0, P01=0, P10=0, P11=0;
};
double Q_angle=0.001, Q_bias=0.003, R_meas=0.03;

void kalmanUpdate(Kalman1D &kf, double gyroRate, double accelAngle, double dt, bool useAccel){
  // Predict
  double rate = gyroRate - kf.bias;
  kf.angle += dt * rate;

  kf.P00 += dt * (dt*kf.P11 - kf.P01 - kf.P10 + Q_angle);
  kf.P01 -= dt * kf.P11;
  kf.P10 -= dt * kf.P11;
  kf.P11 += Q_bias * dt;

  if(!useAccel) return; // skip update when accel is unreliable

  // Update
  double S = kf.P00 + R_meas;
  double K0 = kf.P00 / S;
  double K1 = kf.P10 / S;

  double y = accelAngle - kf.angle;
  kf.angle += K0 * y;
  kf.bias  += K1 * y;

  double P00t = kf.P00, P01t = kf.P01;
  kf.P00 -= K0 * P00t;
  kf.P01 -= K0 * P01t;
  kf.P10 -= K1 * P00t;
  kf.P11 -= K1 * P01t;
}

Kalman1D kfPitch, kfRoll;

// --- PROTOTYPES ---
void MPU9250_init();
void AK8963_init();
void calibrate_MPU9250();
void calibrate_AK8963();
void computeMagScales();
void read_MPU9250(double &ax,double &ay,double &az,double &gx,double &gy,double &gz);
bool read_AK8963(double &mx,double &my,double &mz);

void setup(){
  Serial.begin(115200);
  Wire.begin();

  MPU9250_init();
  AK8963_init();

  calibrate_MPU9250();
  calibrate_AK8963();
  computeMagScales();

  lastMicros = micros();
}

void loop(){
  double ax,ay,az,gx,gy,gz, mx,my,mz;
  read_MPU9250(ax,ay,az,gx,gy,gz);
  read_AK8963(mx,my,mz); // if false, keep previous mag

  unsigned long now = micros();
  dt = (now - lastMicros) * 1e-6;
  lastMicros = now;
  if(dt <= 0) dt = 0.001;

  // Accelerometer angles (deg)
  double accPitch = atan2(-ax, sqrt(ay*ay + az*az)) * 180.0/PI;
  double accRoll  = atan2( ay, az) * 180.0/PI;

  // Gate accel update when |a| deviates from 1 g (±0.15 g window)
  double anorm = sqrt(ax*ax + ay*ay + az*az);
  bool trustAcc = fabs(anorm - 1.0) < 0.15;

  kalmanUpdate(kfPitch, gx - gyroXOffset, accPitch, dt, trustAcc);
  kalmanUpdate(kfRoll,  gy - gyroYOffset, accRoll,  dt, trustAcc);

  pitch = kfPitch.angle;
  roll  = kfRoll.angle;

  // Magnetometer hard/soft-iron correction (µT)
  double cx = (magMinX + magMaxX)/2.0, cy = (magMinY + magMaxY)/2.0, cz = (magMinZ + magMaxZ)/2.0;
  double mxC = (mx - cx) * magScaleX;
  double myC = (my - cy) * magScaleY;
  double mzC = (mz - cz) * magScaleZ;

  // Tilt-compensated yaw (deg)
  double pr = pitch * PI/180.0, rr = roll * PI/180.0;
  double cp = cos(pr), sp = sin(pr), cr = cos(rr), sr = sin(rr);
  double Xh = mxC*cp + mzC*sp;
  double Yh = mxC*sr*sp + myC*cr - mzC*sr*cp;

  double heading = atan2(Yh, Xh);
  const double declination = 0.0928; // radians (+5°19')
  heading += declination;
  if(heading < 0) heading += 2*PI;
  if(heading > 2*PI) heading -= 2*PI;
  yaw = heading * 180.0/PI;

  Serial.print(pitch,1); Serial.print(" / ");
  Serial.print(roll,1);  Serial.print(" / ");
  Serial.println(yaw,1);

  delay(5); // ~200 Hz loop
}

// ---------- IMU ----------
void MPU9250_init(){
  // Wake
  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_PWR_MGMT_1);
  Wire.write(0x00);
  Wire.endTransmission(true);
  delay(10);

  // DLPF ~41 Hz (CONFIG=3), sample rate 1kHz/(1+4)=200 Hz
  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_CONFIG);
  Wire.write(0x03);
  Wire.endTransmission(true);

  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_SMPLRT_DIV);
  Wire.write(0x04);
  Wire.endTransmission(true);

  // Gyro ±250 dps
  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_GYRO_CFG);
  Wire.write(0x00);
  Wire.endTransmission(true);

  // Accel ±2 g
  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_ACCEL_CFG);
  Wire.write(0x00);
  Wire.endTransmission(true);

  // Bypass to access AK8963 directly
  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_INT_PIN_CFG);
  Wire.write(0x02);
  Wire.endTransmission(true);
}

void read_MPU9250(double &ax,double &ay,double &az,double &gx,double &gy,double &gz){
  Wire.beginTransmission(MPU9250_ADDR);
  Wire.write(MPU_REG_ACCEL_XOUT_H);
  Wire.endTransmission(false);
  Wire.requestFrom(MPU9250_ADDR, 14, true);
  int16_t axraw = (Wire.read()<<8)|Wire.read();
  int16_t ayraw = (Wire.read()<<8)|Wire.read();
  int16_t azraw = (Wire.read()<<8)|Wire.read();
  Wire.read(); Wire.read(); // temp
  int16_t gxraw = (Wire.read()<<8)|Wire.read();
  int16_t gyraw = (Wire.read()<<8)|Wire.read();
  int16_t gzraw = (Wire.read()<<8)|Wire.read();

  ax = axraw/ACCEL_SCALE - accelOffsetX;
  ay = ayraw/ACCEL_SCALE - accelOffsetY;
  az = azraw/ACCEL_SCALE - accelOffsetZ;
  gx = gxraw/GYRO_SCALE  - gyroXOffset;
  gy = gyraw/GYRO_SCALE  - gyroYOffset;
  gz = gzraw/GYRO_SCALE  - gyroZOffset;
}

// ---------- MAG ----------
void AK8963_init(){
  // Power down
  Wire.beginTransmission(AK8963_ADDR);
  Wire.write(AK8963_REG_CNTL1);
  Wire.write(0x00);
  Wire.endTransmission();
  delay(10);
  // 16-bit, continuous mode 2 (100 Hz)
  Wire.beginTransmission(AK8963_ADDR);
  Wire.write(AK8963_REG_CNTL1);
  Wire.write(0x16);
  Wire.endTransmission();
  delay(10);
}

bool read_AK8963(double &mx,double &my,double &mz){
  // Check data ready
  Wire.beginTransmission(AK8963_ADDR);
  Wire.write(AK8963_REG_ST1);
  Wire.endTransmission(false);
  Wire.requestFrom(AK8963_ADDR, 1, true);
  if(Wire.available()<1) return false;
  uint8_t st1 = Wire.read();
  if(!(st1 & 0x01)) return false; // not ready

  // Read 6 data + ST2
  Wire.beginTransmission(AK8963_ADDR);
  Wire.write(AK8963_REG_HXL);
  Wire.endTransmission(false);
  Wire.requestFrom(AK8963_ADDR, 7, true);
  if(Wire.available()<7) return false;

  int16_t x = (int16_t)(Wire.read() | (Wire.read()<<8)); // LSB first
  int16_t y = (int16_t)(Wire.read() | (Wire.read()<<8));
  int16_t z = (int16_t)(Wire.read() | (Wire.read()<<8));
  uint8_t st2 = Wire.read();
  if(st2 & 0x08) return false; // overflow

  mx = x * MAG_SCALE;
  my = y * MAG_SCALE;
  mz = z * MAG_SCALE;
  return true;
}

// ---------- CALIBRATION ----------
void calibrate_MPU9250(){
  Serial.println("Calibrating gyro (keep still)...");
  const int N=300;
  double ax,ay,az,gx,gy,gz;
  for(int i=0;i<N;i++){
    read_MPU9250(ax,ay,az,gx,gy,gz);
    gyroXOffset += gx;
    gyroYOffset += gy;
    gyroZOffset += gz;
    delay(5);
  }
  gyroXOffset/=N; gyroYOffset/=N; gyroZOffset/=N;
  Serial.println("Gyro bias done.");
  // If your accel is already good, leave accel offsets at 0.
}

void calibrate_AK8963(){
  Serial.println("Mag calib: rotate slowly for 10 s.");
  magMinX=magMinY=magMinZ= 1e6;
  magMaxX=magMaxY=magMaxZ=-1e6;
  unsigned long t0 = millis();
  double mx,my,mz;
  while(millis()-t0 < 10000){
    if(read_AK8963(mx,my,mz)){
      if(mx<magMinX) magMinX=mx; if(mx>magMaxX) magMaxX=mx;
      if(my<magMinY) magMinY=my; if(my>magMaxY) magMaxY=my;
      if(mz<magMinZ) magMinZ=mz; if(mz>magMaxZ) magMaxZ=mz;
    }
    delay(10);
  }
  Serial.println("Mag min/max set.");
}

void computeMagScales(){
  double rx = (magMaxX - magMinX)/2.0;
  double ry = (magMaxY - magMinY)/2.0;
  double rz = (magMaxZ - magMinZ)/2.0;
  double r_avg = (rx + ry + rz)/3.0;
  if(rx>0) magScaleX = r_avg/rx;
  if(ry>0) magScaleY = r_avg/ry;
  if(rz>0) magScaleZ = r_avg/rz;
}
