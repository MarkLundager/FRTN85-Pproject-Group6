#include <Servo.h>
// 8 - 0 on top
// 9 - 180 on top
// 10 - 0 on top
// 11 - 180 on top
// 12 - 0 on top
// 13 - 180 on top


Servo myServo;  // Create servo object
int servoPin = 11;  // PWM-capable pin

void setup() {
  Serial.begin(115200);
  myServo.attach(servoPin);  // attach servo signal to pin 9
  Serial.println("Servo test started.");
}

void loop() {
  // Ask the user for an angle through Serial
  if (Serial.available() > 0) {
    int angle = Serial.parseInt();  // read integer from Serial
    if (angle >= 0 && angle <= 180) {
      myServo.write(angle);  // set servo angle
      Serial.print("Moved to ");
      Serial.print(angle);
      Serial.println(" degrees.");
    } else {
      Serial.println("Please send a value between 0 and 180.");
    }

    // clear buffer
    while (Serial.available() > 0) Serial.read();
  }
}
