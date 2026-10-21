"""Linear stage control over serial, matching firmware/stage_controller.

Protocol (newline terminated ASCII, 115200 baud):
    H          set the current position as zero
    A <steps>  absolute move, with backlash taken up from the negative side
    M <steps>  relative move
    P          query position
    R          release the motor coils
Replies are "OK <position_steps>" or "ERR <reason>".
"""

import time


class SerialStage:
    def __init__(self, port, mm_per_step, baud=115200, timeout=30.0):
        import serial
        self.ser = serial.Serial(port, baud, timeout=timeout)
        time.sleep(2.0)  # the Uno resets when the port opens
        self.ser.reset_input_buffer()
        self.mm_per_step = mm_per_step
        self._pos = self._cmd("P")

    def _cmd(self, line):
        self.ser.write((line + "\n").encode())
        reply = self.ser.readline().decode().strip()
        if not reply.startswith("OK"):
            raise RuntimeError(f"Stage replied {reply!r} to {line!r}")
        return int(reply.split()[1])

    def home(self):
        self._pos = self._cmd("H")

    def move_to_mm(self, z_mm):
        self._pos = self._cmd(f"A {round(z_mm / self.mm_per_step)}")

    def position_mm(self):
        return self._pos * self.mm_per_step

    def release(self):
        self._cmd("R")

    def close(self):
        self.ser.close()


class SimulatedStage:
    def __init__(self, mm_per_step=0.000610, backlash_steps=0):
        self.mm_per_step = mm_per_step
        self._pos = 0

    def home(self):
        self._pos = 0

    def move_to_mm(self, z_mm):
        self._pos = round(z_mm / self.mm_per_step)

    def position_mm(self):
        return self._pos * self.mm_per_step

    def release(self):
        pass

    def close(self):
        pass
