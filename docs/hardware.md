# Hardware build

## Bill of materials (about S$50)

| Part | Notes |
|---|---|
| Arduino Uno (or clone) | runs `firmware/stage_controller` |
| 28BYJ-48 5V stepper + ULN2003 driver board | usually sold as a pair |
| M8 threaded rod, 150 mm, and two M8 nuts | 1.25 mm pitch lead screw |
| 5 mm to 8 mm flexible shaft coupler | motor to rod |
| Two 8 mm smooth rods or a drawer slide | guides the carriage so it does not rotate |
| USB webcam | any cheap UVC webcam; the lens gets removed |
| Red laser pointer, Class 2 (under 1 mW) | 650 nm |
| Plano convex lens, f = 150 to 200 mm | longer focal length gives a bigger, easier to sample waist |
| Neutral density film or two linear polarisers | attenuates the beam so the sensor does not saturate |

## Layout

```
laser ── ND / polarisers ── lens ─────── focus ───── webcam sensor on carriage
                                              ◄──── z scan ────►
```

The webcam carriage rides on the guide rails and is pushed along z by the nut
on the lead screw. Mount the laser and lens rigidly on the same base so only
the camera moves.

## Getting the numbers right

**Steps per mm.** 2048 steps per revolution and 1.25 mm per revolution gives
0.000610 mm per step, which is the default `--mm-per-step`. Check it by moving
2048 x 20 steps and measuring the carriage travel with calipers.

**Backlash.** Move +2000 steps, then −2000, then +2000 again while watching a
dial indicator or a scale under a phone camera. The lost motion on reversal is
the backlash. Set `BACKLASH_STEPS` in the firmware to a little more than it.
Absolute moves always finish travelling forwards, so it cancels out.

**Pixel pitch.** Look up the sensor in the webcam teardown or datasheet. If it
is unknown, image a ruler at a known distance with the lens still on, or
translate the beam sideways by a known stage distance and read the centroid
shift in pixels.

**Gamma.** Run `scripts/linearity.py` once and pass the gamma it reports to
`scan.py`. Gamma encoded data inflates second moment widths.

**Exposure.** Aim for a peak around 60 to 80 percent of full scale at the
waist. `scan.py` flags any saturated frame; reduce exposure or add ND and rescan.

**Expected waist.** For a beam of diameter D on a lens of focal length f,
the waist diameter is roughly 4 M² λ f / (π D). A 1.5 mm pointer beam through
f = 200 mm gives about 110 to 140 µm, which is 35 to 45 pixels at 3 µm pitch:
plenty for an accurate D4σ. The Rayleigh length is then around 10 to 15 mm, so
a ±40 mm scan covers the ISO requirement of half the points beyond 2 zR.

## Laser safety

Class 2 pointers are safe for accidental exposure only because of the blink
reflex. Never look into the beam or its reflections, keep the beam below eye
level, block it with a matte card at the end of the bench, and do not use
pointers sold as Class 3 or above.
