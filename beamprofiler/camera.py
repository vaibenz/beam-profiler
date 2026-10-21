"""Camera sources: a real USB webcam through OpenCV, and a simulated one.

Webcam caveats that matter for beam profiling:
  * Remove the lens so the beam falls directly on the bare sensor.
  * Turn off auto exposure and auto gain, or every frame has a different scale.
  * Most webcams apply a gamma curve. Run scripts/linearity.py once and pass
    the measured gamma here so pixel values are linear in intensity again.
  * The sensor sits behind a Bayer filter. For a red laser only the red pixels
    respond strongly, so the default is to take the red channel.
"""

import time

import numpy as np

from . import simulate


class WebcamSource:
    def __init__(self, index=0, exposure=None, gamma=1.0, channel="red", width=None, height=None):
        import cv2
        self.cv2 = cv2
        self.cap = cv2.VideoCapture(index)
        if not self.cap.isOpened():
            raise RuntimeError(f"Could not open camera {index}")
        if width:
            self.cap.set(cv2.CAP_PROP_FRAME_WIDTH, width)
        if height:
            self.cap.set(cv2.CAP_PROP_FRAME_HEIGHT, height)
        # 1 means manual exposure on the V4L2 backend; values are driver specific.
        self.cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 1)
        self.cap.set(cv2.CAP_PROP_GAIN, 0)
        if exposure is not None:
            self.cap.set(cv2.CAP_PROP_EXPOSURE, exposure)
        self.gamma = gamma
        self.channel = {"blue": 0, "green": 1, "red": 2}.get(channel)
        for _ in range(5):  # let exposure settle
            self.cap.read()

    def grab(self):
        ok, bgr = self.cap.read()
        if not ok:
            raise RuntimeError("Frame grab failed")
        raw = bgr[:, :, self.channel] if self.channel is not None else self.cv2.cvtColor(bgr, self.cv2.COLOR_BGR2GRAY)
        return raw

    def linearise(self, frame):
        f = np.asarray(frame, dtype=np.float64) / 255.0
        return (f ** self.gamma) * 255.0

    def close(self):
        self.cap.release()


class SimulatedSource:
    """Stands in for the webcam. Reads the stage position to decide what to render."""

    def __init__(self, stage, beam=None, cam=None, seed=1, peak_fraction=0.7):
        self.stage = stage
        self.beam = beam or simulate.BeamSpec()
        self.cam = cam or simulate.CameraSpec()
        self.rng = np.random.default_rng(seed)
        self.peak_fraction = peak_fraction
        self.gamma = 1.0
        self.shutter_open = True

    def grab(self):
        if not self.shutter_open:
            return simulate.frame(0, self.beam, self.cam, peak_fraction=0.0, rng=self.rng)
        return simulate.frame(self.stage.position_mm(), self.beam, self.cam, self.peak_fraction, self.rng)

    def linearise(self, frame):
        return np.asarray(frame, dtype=np.float64)

    def close(self):
        pass


def average(source, n=8, delay_s=0.0):
    """Average n linearised frames to beat down noise. Returns (mean_frame, saturated_fraction)."""
    frames, sat = [], 0.0
    for _ in range(n):
        raw = source.grab()
        sat = max(sat, float(np.mean(raw >= 255)))
        frames.append(source.linearise(raw))
        if delay_s:
            time.sleep(delay_s)
    return np.mean(frames, axis=0), sat
