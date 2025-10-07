import processing.serial.*;
Serial myPort;
String data="";
float roll, pitch, yaw;

void setup() {
  size(960, 640, P3D);
  println(Serial.list());                 // check port index if needed
  myPort = new Serial(this, "/dev/ttyACM0", 115200);
  myPort.clear();                         // flush old bytes
  delay(200);                             // let the buffer settle
  myPort.bufferUntil('\n');
}

void draw() {
  translate(width/2, height/2, 0);
  background(233);
  textSize(22);
  text("Roll: " + int(roll) + "     Pitch: " + int(pitch) + "     Yaw: " + int(yaw),
       -200, 265);

  rotateX(radians(-pitch));
  rotateZ(radians(roll));
  rotateY(radians(yaw));

  fill(0, 76, 153); box(386, 40, 200);
}

void serialEvent(Serial p) {
  try {
    String s = p.readStringUntil('\n');
    if (s == null) return;
    s = trim(s);
    if (s.length() == 0) return;

    // Accept "/" or commas, ignore spaces/tabs
    String[] items = splitTokens(s, "/, \t");
    if (items.length == 3) {
      // Arduino now sends roll/pitch/yaw in this order
      roll  = float(items[0]);
      pitch = float(items[1]);
      yaw   = float(items[2]);
    }
  } catch (Exception e) {
    // ignore malformed line and continue
  }
}
