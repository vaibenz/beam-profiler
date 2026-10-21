// Linear stage controller for the beam profiler.
//
// Hardware: Arduino Uno, 28BYJ-48 5V stepper, ULN2003 driver board.
// The motor turns an M8 lead screw (1.25 mm pitch) that carries the camera.
// 2048 full steps per output revolution gives about 0.61 um per step.
//
// Serial protocol, 115200 baud, newline terminated:
//   H          set current position to zero
//   A <steps>  absolute move with backlash compensation
//   M <steps>  relative move (no compensation)
//   P          report position
//   R          release coils (stops the motor heating while idle)
// Every command replies "OK <position>" or "ERR <reason>".

#include <Stepper.h>

const int STEPS_PER_REV = 2048;
const int BACKLASH_STEPS = 60;   // measure for your screw and nut, see docs/hardware.md
const int RPM = 12;              // the 28BYJ-48 stalls above roughly 15 rpm

// ULN2003 IN1..IN4 on pins 8..11. The Stepper library wants them as 1, 3, 2, 4.
const int IN1 = 8, IN2 = 9, IN3 = 10, IN4 = 11;
Stepper motor(STEPS_PER_REV, IN1, IN3, IN2, IN4);

long position = 0;
String line;

void release_coils() {
  digitalWrite(IN1, LOW);
  digitalWrite(IN2, LOW);
  digitalWrite(IN3, LOW);
  digitalWrite(IN4, LOW);
}

void step_by(long n) {
  // Stepper::step blocks, and takes an int, so move in chunks.
  while (n != 0) {
    int chunk = (n > 1000) ? 1000 : (n < -1000 ? -1000 : (int)n);
    motor.step(chunk);
    n -= chunk;
    position += chunk;
  }
}

void move_absolute(long target) {
  // Always finish the move travelling in the positive direction so the nut
  // sits against the same flank of the thread every time.
  if (target < position) {
    step_by(target - BACKLASH_STEPS - position);
  }
  step_by(target - position);
}

void reply_ok() {
  Serial.print("OK ");
  Serial.println(position);
}

void handle(String cmd) {
  cmd.trim();
  if (cmd.length() == 0) return;
  char op = cmd.charAt(0);
  long arg = cmd.length() > 2 ? cmd.substring(2).toInt() : 0;
  switch (op) {
    case 'H': position = 0; reply_ok(); break;
    case 'A': move_absolute(arg); reply_ok(); break;
    case 'M': step_by(arg); reply_ok(); break;
    case 'P': reply_ok(); break;
    case 'R': release_coils(); reply_ok(); break;
    default: Serial.println("ERR unknown command");
  }
}

void setup() {
  Serial.begin(115200);
  motor.setSpeed(RPM);
  release_coils();
}

void loop() {
  while (Serial.available()) {
    char c = Serial.read();
    if (c == '\n') {
      handle(line);
      line = "";
    } else if (c != '\r') {
      line += c;
    }
  }
}
