# Backscatter Designer

Interactive explorer for backscatter radio design parameters. Static site, no build step.

**Live:** https://matthew-marsicano.github.io/backscatter-radio-viz/

## What it shows

- **Link budget vs distance.** Desired sideband power for a d2 sweep at fixed d1, and for a co-located carrier and receiver. Sensitivity and carrier leakage overlaid.
- **Load states.** Reflection coefficient of each switch state on a Smith chart, against the ideal N-phase targets. Hover for |Γ|, angle and Z.
- **Image rejection.** IRR over phase error and amplitude mismatch of one state, with the 20 dB contour.
- **Spectrum at receiver.** Carrier leakage and every tag harmonic at fc + k·f_sw, with BLE advertising channels 37/38/39 marked. Table view included.
- **Bistatic coverage.** Received sideband power over tag position with the carrier at the origin and the receiver at d0. Sensitivity contour drawn.
- **Summary tiles.** Output frequency and BLE channel, Rx power, link margin, max range, conversion loss, IRR, sideband-to-leakage ratio.

Loads can be ideal phasors with an injected error on state 1, or arbitrary R + jX per state (Z0 = 50 Ω). All parameters live in the URL hash, so a design can be shared with **Copy link**.

## Model

- Bistatic radar equation: `Prx = Ptx Gtx Gtag² Grx λ⁴ |ck|² / ((4π)⁴ d1² d2²)`, minus switch insertion loss twice.
- The tag steps through N load states with equal dwell. Harmonic k of the stepped Γ(t):
  `ck = sin(πk/N)/(πk) · e^(−jπk/N) · Σn Γn e^(−j2πkn/N)`
- Ideal stepping keeps only k ≡ 1 (mod N). N = 2 is DSB. N ≥ 3 is SSB.
- Carrier leakage is free-space path loss over d0 minus user isolation.
- Not modeled: multipath, phase noise, polarization mismatch, near field, switch rise time.

## Run locally

```bash
python -m http.server 8000
```

Then open http://localhost:8000.

## Files

- `index.html` layout
- `style.css` theme tokens (light and dark)
- `app.js` model, controls, Plotly charts
