// Backscatter Designer: link budget, load states, spectrum and coverage for an N-phase tag.
"use strict";

const C0 = 299792458;
const Z0 = 50;
const BLE_ADV = [{ ch: 37, f: 2402 }, { ch: 38, f: 2426 }, { ch: 39, f: 2480 }];

// Parameter definitions. Slider ranges are suggestions; the number box accepts anything.
const PARAMS = {
  fc:   { label: "Carrier freq",        unit: "MHz", min: 2400, max: 2500, step: 0.5, def: 2410 },
  ptx:  { label: "Tx power",            unit: "dBm", min: -10,  max: 30,   step: 0.5, def: 20 },
  gtx:  { label: "Tx antenna gain",     unit: "dBi", min: -5,   max: 15,   step: 0.5, def: 2 },
  grx:  { label: "Rx antenna gain",     unit: "dBi", min: -5,   max: 15,   step: 0.5, def: 2 },
  sens: { label: "Rx sensitivity",      unit: "dBm", min: -110, max: -50,  step: 1,   def: -90 },
  d0:   { label: "Carrier to Rx",       unit: "m",   min: 0.1,  max: 20,   step: 0.1, def: 3 },
  iso:  { label: "Extra Tx/Rx isolation", unit: "dB", min: 0,   max: 60,   step: 1,   def: 0 },
  d1:   { label: "Carrier to tag d1",   unit: "m",   min: 0.1,  max: 10,   step: 0.1, def: 1 },
  d2:   { label: "Tag to Rx d2",        unit: "m",   min: 0.1,  max: 30,   step: 0.1, def: 3 },
  gtag: { label: "Tag antenna gain",    unit: "dBi", min: -10,  max: 6,    step: 0.5, def: 0 },
  il:   { label: "Switch IL (one way)", unit: "dB",  min: 0,    max: 6,    step: 0.1, def: 1 },
  fsw:  { label: "Switch freq f_sw",    unit: "MHz", min: 1,    max: 80,   step: 0.5, def: 16 },
  nph:  { label: "Phase states N",      unit: "",    min: 2,    max: 8,    step: 1,   def: 4 },
  gmag: { label: "|Γ| of states",       unit: "",    min: 0.1,  max: 1,    step: 0.01, def: 1 },
  perr: { label: "Phase error",         unit: "deg", min: -45,  max: 45,   step: 0.5, def: 5 },
  aerr: { label: "Amplitude mismatch",  unit: "dB",  min: -6,   max: 6,    step: 0.1, def: 0.5 },
};

const state = { mode: "phasor", z: [] };
for (const k in PARAMS) state[k] = PARAMS[k].def;

// ---------- complex helpers ----------
const cx = (re, im = 0) => ({ re, im });
const cadd = (a, b) => cx(a.re + b.re, a.im + b.im);
const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
const cexp = (th) => cx(Math.cos(th), Math.sin(th));
const cabs2 = (a) => a.re * a.re + a.im * a.im;
const db10 = (x) => 10 * Math.log10(Math.max(x, 1e-30));

const zToGamma = (z) => cdiv(cx(z.re - Z0, z.im), cx(z.re + Z0, z.im));
const gammaToZ = (g) => cmul(cx(Z0), cdiv(cx(1 + g.re, g.im), cx(1 - g.re, -g.im)));

// ---------- model ----------
function idealStates(N, mag) {
  return Array.from({ length: N }, (_, n) => { const e = cexp(2 * Math.PI * n / N); return cx(mag * e.re, mag * e.im); });
}

function phasorStates(N, mag, perrDeg, aerrDb) {
  const g = idealStates(N, mag);
  const a = Math.pow(10, -aerrDb / 20);
  const rot = cexp(perrDeg * Math.PI / 180);
  g[1] = cmul(cx(a * g[1].re, a * g[1].im), rot);
  return g;
}

function currentStates() {
  const N = state.nph;
  if (state.mode === "z") return state.z.slice(0, N).map(zToGamma);
  return phasorStates(N, state.gmag, state.perr, state.aerr);
}

// Fourier coefficient of the stepped reflection waveform at harmonic k.
function coeff(G, k) {
  const N = G.length;
  let s = cx(0);
  for (let n = 0; n < N; n++) s = cadd(s, cmul(G[n], cexp(-2 * Math.PI * k * n / N)));
  if (k === 0) return cx(s.re / N, s.im / N);
  const amp = Math.sin(Math.PI * k / N) / (Math.PI * k);
  return cmul(cx(amp * s.re, amp * s.im), cexp(-Math.PI * k / N));
}

const lambda = () => C0 / (state.fc * 1e6);

// Received power in dBm for a harmonic with |c|^2 = m, at distances d1, d2.
function prxDbm(m, d1, d2) {
  const L = lambda();
  return state.ptx + state.gtx + 2 * state.gtag + state.grx
    + 40 * Math.log10(L) - 40 * Math.log10(4 * Math.PI)
    - 20 * Math.log10(Math.max(d1, 1e-3)) - 20 * Math.log10(Math.max(d2, 1e-3))
    + db10(m) - 2 * state.il;
}

function leakDbm() {
  return state.ptx + state.gtx + state.grx + 20 * Math.log10(lambda() / (4 * Math.PI * state.d0)) - state.iso;
}

function nearestBle(fMHz) {
  const ch = Math.round((fMHz - 2402) / 2);
  if (ch < 0 || ch > 39) return null;
  const center = 2402 + 2 * ch;
  // Map RF channel index to BLE channel number.
  let num;
  if (ch === 0) num = 37; else if (ch === 12) num = 38; else if (ch === 39) num = 39;
  else if (ch < 12) num = ch - 1; else num = ch - 2;
  return { num, center, off: fMHz - center };
}

// ---------- theming ----------
function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }

function baseLayout(extra = {}) {
  const t = css("--text-secondary"), g = css("--grid");
  const axis = { color: t, gridcolor: g, zerolinecolor: g, linecolor: css("--border"), tickfont: { size: 11 }, titlefont: { size: 12 } };
  return Object.assign({
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
    font: { family: "Inter, system-ui, sans-serif", color: t, size: 12 },
    margin: { l: 56, r: 16, t: 8, b: 44 },
    hoverlabel: { bgcolor: css("--surface"), bordercolor: css("--border"), font: { color: css("--text-primary"), family: "JetBrains Mono, monospace", size: 12 } },
    legend: { orientation: "h", x: 0, y: 1.12, font: { size: 12 } },
    xaxis: Object.assign({}, axis), yaxis: Object.assign({}, axis),
  }, extra);
}
// Sequential blue. Low end recedes toward the surface in each theme.
function seqRamp() {
  const dark = getComputedStyle(document.documentElement).colorScheme.includes("dark");
  return dark
    ? [[0, "#16263b"], [0.25, "#184f95"], [0.5, "#2a78d6"], [0.75, "#6da7ec"], [1, "#cde2fb"]]
    : [[0, "#cde2fb"], [0.25, "#86b6ef"], [0.5, "#3987e5"], [0.75, "#1c5cab"], [1, "#0d366b"]];
}
const PCFG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"] };

// ---------- controls ----------
function buildControls() {
  document.querySelectorAll(".ctl[data-key]").forEach((el) => {
    const k = el.dataset.key, p = PARAMS[k];
    el.innerHTML = `
      <label for="n-${k}"><span>${p.label}</span>
        <span class="num"><input type="number" id="n-${k}" step="${p.step}" value="${state[k]}"><span class="unit">${p.unit}</span></span>
      </label>
      <input type="range" id="r-${k}" min="${p.min}" max="${p.max}" step="${p.step}" value="${state[k]}" aria-label="${p.label}">`;
    const num = el.querySelector("input[type=number]"), rng = el.querySelector("input[type=range]");
    const set = (v, src) => {
      v = parseFloat(v);
      if (!isFinite(v)) return;
      if (k === "nph") v = Math.max(2, Math.min(16, Math.round(v)));
      state[k] = v;
      if (src !== num) num.value = v;
      if (src !== rng) rng.value = v;
      if (k === "nph" || ((k === "gmag" || k === "perr" || k === "aerr") && state.mode === "phasor")) seedZ();
      schedule();
    };
    num.addEventListener("input", () => set(num.value, num));
    rng.addEventListener("input", () => set(rng.value, rng));
  });

  document.querySelectorAll(".seg button").forEach((b) => b.addEventListener("click", () => {
    state.mode = b.dataset.mode;
    document.querySelectorAll(".seg button").forEach((x) => x.classList.toggle("on", x === b));
    document.getElementById("phasor-ctl").hidden = state.mode !== "phasor";
    document.getElementById("z-ctl").hidden = state.mode !== "z";
    schedule();
  }));
}

function syncControls() {
  for (const k in PARAMS) {
    const n = document.getElementById("n-" + k), r = document.getElementById("r-" + k);
    if (n) n.value = state[k];
    if (r) r.value = state[k];
  }
  document.querySelectorAll(".seg button").forEach((x) => x.classList.toggle("on", x.dataset.mode === state.mode));
  document.getElementById("phasor-ctl").hidden = state.mode !== "phasor";
  document.getElementById("z-ctl").hidden = state.mode !== "z";
}

// Seed impedance list from the phasor model, keeping |Γ| below 1 so Z stays finite.
function seedZ() {
  const g = phasorStates(state.nph, Math.min(state.gmag, 0.98), state.perr, state.aerr);
  state.z = g.map((x) => {
    const m = Math.sqrt(cabs2(x));
    const s = m > 0.98 ? 0.98 / m : 1;
    const z = gammaToZ(cx(x.re * s, x.im * s));
    return cx(+z.re.toFixed(2), +z.im.toFixed(2));
  });
  renderZList();
}

function renderZList() {
  const host = document.getElementById("z-list");
  host.innerHTML = `<div class="zrow"><span></span><span>R (Ω)</span><span>X (Ω)</span></div>` +
    state.z.map((z, i) => `<div class="zrow"><span>${i}</span>
      <input type="number" step="any" data-i="${i}" data-p="re" value="${z.re}" aria-label="State ${i} resistance">
      <input type="number" step="any" data-i="${i}" data-p="im" value="${z.im}" aria-label="State ${i} reactance"></div>`).join("");
  host.querySelectorAll("input").forEach((inp) => inp.addEventListener("input", () => {
    const v = parseFloat(inp.value);
    if (!isFinite(v)) return;
    state.z[+inp.dataset.i][inp.dataset.p] = v;
    schedule();
  }));
}

// ---------- rendering ----------
let pending = false;
function schedule() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => { pending = false; render(); writeHash(); });
}

function tile(k, v, d = "", flag = null) {
  const f = flag ? `<div class="flag ${flag.ok ? "ok" : "no"}">${flag.ok ? "✓" : "✗"} ${flag.text}</div>` : "";
  return `<div class="tile"><div class="k">${k}</div><div class="v">${v}</div>${d ? `<div class="d">${d}</div>` : ""}${f}</div>`;
}

function render() {
  const G = currentStates();
  const c1 = coeff(G, 1), cm1 = coeff(G, -1);
  const m1 = cabs2(c1);
  const irr = db10(m1 / Math.max(cabs2(cm1), 1e-30));
  const p = prxDbm(m1, state.d1, state.d2);
  const margin = p - state.sens;
  const convLoss = -db10(m1) + 2 * state.il;
  // Co-located range (d1 = d2 = d): solve 40 log d = P(1 m) - sens.
  const p1 = prxDbm(m1, 1, 1);
  const rCo = Math.pow(10, (p1 - state.sens) / 40);
  const rD2 = Math.pow(10, (prxDbm(m1, state.d1, 1) - state.sens) / 20);
  const fOut = state.fc + state.fsw;
  const ble = nearestBle(fOut);
  const bleTxt = ble ? `ch ${ble.num} (${ble.center} MHz), ${ble.off >= 0 ? "+" : ""}${ble.off.toFixed(2)} MHz` : "outside BLE band";
  const onAdv = ble && [37, 38, 39].includes(ble.num) && Math.abs(ble.off) < 0.15;
  const sbOverLeak = p - leakDbm();

  document.getElementById("tiles").innerHTML = [
    tile("Desired sideband", `${fOut.toFixed(1)} MHz`, bleTxt, ble ? { ok: onAdv, text: onAdv ? "advertising channel" : "not on adv channel" } : null),
    tile("Rx power at d1, d2", `${p.toFixed(1)} dBm`, `${state.d1} m / ${state.d2} m`),
    tile("Link margin", `${margin >= 0 ? "+" : ""}${margin.toFixed(1)} dB`, `vs ${state.sens} dBm`, { ok: margin >= 0, text: margin >= 0 ? "closes" : "short" }),
    tile("Max d2 at d1", `${fmtDist(rD2)}`, `co-located range ${fmtDist(rCo)}`),
    tile("Conversion loss", `${convLoss.toFixed(1)} dB`, `|c₁|² ${db10(m1).toFixed(2)} dB + 2×IL`),
    tile("Image rejection", state.nph === 2 ? "0 dB" : (irr > 80 ? "> 80 dB" : `${irr.toFixed(1)} dB`), state.nph === 2 ? "DSB: image is symmetric" : "desired vs −f_sw image", state.nph === 2 ? null : { ok: irr >= 20, text: irr >= 20 ? "meets 20 dB" : "below 20 dB" }),
    tile("Sideband vs leakage", `${sbOverLeak.toFixed(1)} dBc`, `leakage ${leakDbm().toFixed(1)} dBm at ${state.fsw} MHz offset`),
  ].join("");

  drawLink(m1);
  drawSmith(G);
  drawIrr();
  drawSpec(G);
  drawMap(m1);
}

function fmtDist(d) {
  if (!isFinite(d)) return "∞";
  if (d >= 1000) return (d / 1000).toFixed(2) + " km";
  if (d >= 10) return d.toFixed(1) + " m";
  return d.toFixed(2) + " m";
}

function logspace(a, b, n) { const la = Math.log10(a), lb = Math.log10(b); return Array.from({ length: n }, (_, i) => Math.pow(10, la + (lb - la) * i / (n - 1))); }

function drawLink(m1) {
  const ds = logspace(0.1, 100, 200);
  const co = ds.map((d) => prxDbm(m1, d, d));
  const sw = ds.map((d) => prxDbm(m1, state.d1, d));
  const leak = ds.map(() => leakDbm());
  const traces = [
    { x: ds, y: sw, name: `Sweep d2, d1 = ${state.d1} m`, mode: "lines", line: { color: css("--s1"), width: 2 }, hovertemplate: "d2 %{x:.2f} m<br>%{y:.1f} dBm<extra></extra>" },
    { x: ds, y: co, name: "Co-located d1 = d2", mode: "lines", line: { color: css("--s2"), width: 2 }, hovertemplate: "d %{x:.2f} m<br>%{y:.1f} dBm<extra></extra>" },
    { x: [state.d2], y: [prxDbm(m1, state.d1, state.d2)], name: "Current", mode: "markers", marker: { color: css("--s1"), size: 10, line: { color: css("--surface"), width: 2 } }, showlegend: false, hovertemplate: "current<br>%{y:.1f} dBm<extra></extra>" },
  ];
  const lay = baseLayout({
    xaxis: Object.assign(baseLayout().xaxis, { type: "log", title: { text: "Distance (m)" } }),
    yaxis: Object.assign(baseLayout().yaxis, { title: { text: "Power (dBm)" } }),
    shapes: [
      { type: "line", xref: "paper", x0: 0, x1: 1, y0: state.sens, y1: state.sens, line: { color: css("--text-muted"), width: 1.5, dash: "dash" } },
      { type: "line", xref: "paper", x0: 0, x1: 1, y0: leak[0], y1: leak[0], line: { color: css("--s3"), width: 1.5, dash: "dot" } },
    ],
    annotations: [
      { xref: "paper", x: 1, y: state.sens, text: "sensitivity", showarrow: false, xanchor: "right", yanchor: "bottom", font: { size: 11, color: css("--text-muted") } },
      { xref: "paper", x: 1, y: leak[0], text: "carrier leakage", showarrow: false, xanchor: "right", yanchor: "bottom", font: { size: 11, color: css("--text-muted") } },
    ],
  });
  Plotly.react("p-link", traces, lay, PCFG);
}

function smithGrid() {
  const g = css("--grid"), traces = [];
  const line = (xs, ys, w = 1, col = g) => traces.push({ x: xs, y: ys, mode: "lines", line: { color: col, width: w }, hoverinfo: "skip", showlegend: false });
  const th = Array.from({ length: 181 }, (_, i) => (i / 180) * 2 * Math.PI);
  line(th.map(Math.cos), th.map(Math.sin), 1.5, css("--text-muted"));
  for (const r of [0.2, 0.5, 1, 2, 5]) {
    const c = r / (1 + r), rad = 1 / (1 + r);
    line(th.map((t) => c + rad * Math.cos(t)), th.map((t) => rad * Math.sin(t)));
  }
  for (const x of [0.2, 0.5, 1, 2, 5, -0.2, -0.5, -1, -2, -5]) {
    const xs = [], ys = [];
    for (const t of th) {
      const px = 1 + (1 / Math.abs(x)) * Math.cos(t), py = 1 / x + (1 / Math.abs(x)) * Math.sin(t);
      if (px * px + py * py <= 1.0001) { xs.push(px); ys.push(py); } else { xs.push(null); ys.push(null); }
    }
    line(xs, ys);
  }
  line([-1, 1], [0, 0]);
  return traces;
}

function drawSmith(G) {
  const N = G.length;
  const ideal = state.mode === "phasor" ? idealStates(N, state.gmag) : idealStates(N, Math.sqrt(G.reduce((a, g) => a + cabs2(g), 0) / N));
  const dc = coeff(G, 0);
  const zs = G.map(gammaToZ);
  const traces = smithGrid().concat([
    { x: ideal.map((g) => g.re), y: ideal.map((g) => g.im), name: "Ideal", mode: "markers", marker: { size: 14, color: "rgba(0,0,0,0)", line: { color: css("--text-muted"), width: 1.5 } }, hoverinfo: "skip" },
    { x: G.concat([G[0]]).map((g) => g.re), y: G.concat([G[0]]).map((g) => g.im), mode: "lines", line: { color: css("--s1"), width: 1, dash: "dot" }, hoverinfo: "skip", showlegend: false },
    {
      x: G.map((g) => g.re), y: G.map((g) => g.im), name: "States", mode: "markers+text",
      text: G.map((_, i) => String(i)), textposition: "top right", textfont: { color: css("--text-secondary"), size: 11 },
      marker: { size: 10, color: css("--s1"), line: { color: css("--surface"), width: 2 } },
      customdata: G.map((g, i) => [Math.sqrt(cabs2(g)), Math.atan2(g.im, g.re) * 180 / Math.PI, zs[i].re, zs[i].im]),
      hovertemplate: "state %{text}<br>|Γ| %{customdata[0]:.3f} ∠ %{customdata[1]:.1f}°<br>Z %{customdata[2]:.1f} %{customdata[3]:+.1f}j Ω<extra></extra>",
    },
    { x: [dc.re], y: [dc.im], name: "Mean Γ (carrier)", mode: "markers", marker: { size: 9, symbol: "x", color: css("--s2") }, hovertemplate: "mean Γ %{x:.3f} %{y:+.3f}j<extra></extra>" },
  ]);
  const ax = { range: [-1.15, 1.15], showgrid: false, zeroline: false, showticklabels: false, fixedrange: true };
  Plotly.react("p-smith", traces, baseLayout({
    margin: { l: 8, r: 8, t: 8, b: 8 },
    xaxis: Object.assign({}, ax), yaxis: Object.assign({}, ax, { scaleanchor: "x" }),
    legend: { orientation: "h", x: 0, y: 0, yanchor: "top", font: { size: 11 } },
  }), Object.assign({}, PCFG, { displayModeBar: false }));
}

function drawIrr() {
  const N = state.nph;
  const ps = Array.from({ length: 61 }, (_, i) => -30 + i);
  const as = Array.from({ length: 49 }, (_, i) => -3 + i * 0.125);
  const z = as.map((a) => ps.map((ph) => {
    const G = phasorStates(N, 1, ph, a);
    return Math.min(60, db10(cabs2(coeff(G, 1)) / Math.max(cabs2(coeff(G, -1)), 1e-30)));
  }));
  const ramp = seqRamp();
  const traces = [
    { type: "heatmap", x: ps, y: as, z, colorscale: ramp, zmin: 0, zmax: 60, colorbar: { title: { text: "dB", side: "right" }, thickness: 10, outlinewidth: 0, tickfont: { size: 10 } }, hovertemplate: "phase %{x}°<br>amp %{y:.2f} dB<br>IRR %{z:.1f} dB<extra></extra>" },
    { type: "contour", x: ps, y: as, z, contours: { start: 20, end: 20, size: 1, coloring: "none", showlabels: true, labelfont: { size: 11, color: css("--text-primary") } }, line: { color: css("--text-primary"), width: 1.5 }, showscale: false, showlegend: false, hoverinfo: "skip" },
  ];
  if (state.mode === "phasor") traces.push({ x: [state.perr], y: [state.aerr], mode: "markers", marker: { size: 11, color: css("--s2"), line: { color: css("--surface"), width: 2 } }, name: "Current", hovertemplate: "current<extra></extra>", showlegend: false });
  const lay = baseLayout({
    xaxis: Object.assign(baseLayout().xaxis, { title: { text: "Phase error (deg)" }, range: [-30, 30] }),
    yaxis: Object.assign(baseLayout().yaxis, { title: { text: "Amplitude mismatch (dB)" }, range: [-3, 3] }),
    margin: { l: 56, r: 8, t: 8, b: 44 },
  });
  if (N === 2) lay.annotations = [{ xref: "paper", yref: "paper", x: 0.5, y: 0.5, text: "N = 2 is DSB. No image rejection.", showarrow: false, bgcolor: css("--surface"), font: { color: css("--text-primary") } }];
  Plotly.react("p-irr", traces, lay, PCFG);
}

function drawSpec(G) {
  const K = Math.max(9, 2 * state.nph + 1);
  const rows = [];
  for (let k = -K; k <= K; k++) {
    if (k === 0) continue;
    const m = cabs2(coeff(G, k));
    rows.push({ k, f: state.fc + k * state.fsw, p: prxDbm(m, state.d1, state.d2), m });
  }
  const leak = leakDbm();
  const floor = Math.min(state.sens - 30, ...rows.filter((r) => r.m > 1e-12).map((r) => r.p)) - 10;
  const stems = (sel, name, color) => {
    const rs = rows.filter(sel), xs = [], ys = [];
    rs.forEach((r) => { const p = Math.max(r.p, floor); xs.push(r.f, r.f, null); ys.push(floor, p, null); });
    return [
      { x: xs, y: ys, mode: "lines", line: { color, width: 2 }, hoverinfo: "skip", showlegend: false },
      { x: rs.map((r) => r.f), y: rs.map((r) => Math.max(r.p, floor)), name, mode: "markers", marker: { color, size: 8, line: { color: css("--surface"), width: 2 } }, customdata: rs.map((r) => [r.k, r.p]), hovertemplate: "k = %{customdata[0]}<br>%{x:.1f} MHz<br>%{customdata[1]:.1f} dBm<extra></extra>" },
    ];
  };
  const traces = [
    { x: [state.fc, state.fc], y: [floor, leak], mode: "lines", line: { color: css("--s3"), width: 2 }, hoverinfo: "skip", showlegend: false },
    { x: [state.fc], y: [leak], name: "Carrier leakage", mode: "markers", marker: { color: css("--s3"), size: 8, symbol: "diamond", line: { color: css("--surface"), width: 2 } }, hovertemplate: "carrier<br>%{x:.1f} MHz<br>%{y:.1f} dBm<extra></extra>" },
    ...stems((r) => r.k === 1, "Desired (k = +1)", css("--s1")),
    ...stems((r) => r.k === -1, "Image (k = −1)", css("--s2")),
    ...stems((r) => Math.abs(r.k) > 1, "Other harmonics", css("--neutral")),
  ];
  const fmin = state.fc - K * state.fsw, fmax = state.fc + K * state.fsw;
  const shapes = BLE_ADV.map((b) => ({ type: "rect", yref: "paper", x0: b.f - 1, x1: b.f + 1, y0: 0, y1: 1, fillcolor: css("--band"), line: { width: 0 }, layer: "below" }));
  shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, y0: state.sens, y1: state.sens, line: { color: css("--text-muted"), width: 1.5, dash: "dash" } });
  shapes.push({ type: "rect", yref: "paper", x0: 2400, x1: 2483.5, y0: 0, y1: 1, fillcolor: "rgba(0,0,0,0)", line: { color: css("--grid"), width: 1, dash: "dot" }, layer: "below" });
  const ann = BLE_ADV.filter((b) => b.f > fmin && b.f < fmax).map((b) => ({ x: b.f, yref: "paper", y: 1, text: `ch${b.ch}`, showarrow: false, yanchor: "bottom", font: { size: 10, color: css("--text-muted") } }));
  Plotly.react("p-spec", traces, baseLayout({
    xaxis: Object.assign(baseLayout().xaxis, { title: { text: "Frequency (MHz)" }, range: [fmin - state.fsw * 0.5, fmax + state.fsw * 0.5] }),
    yaxis: Object.assign(baseLayout().yaxis, { title: { text: "Power (dBm)" }, range: [floor, Math.max(leak, ...rows.map((r) => r.p)) + 8] }),
    shapes, annotations: ann, margin: { l: 56, r: 16, t: 20, b: 44 }, legend: { orientation: "h", x: 0, y: -0.22 },
  }), PCFG);

  const tb = document.getElementById("spec-table");
  const all = [{ k: 0, f: state.fc, p: leak, label: "carrier leak" }].concat(rows.map((r) => ({ ...r, label: r.k === 1 ? "desired" : r.k === -1 ? "image" : "" })));
  all.sort((a, b) => a.k - b.k);
  tb.innerHTML = "<tr><th>k</th><th>Freq (MHz)</th><th>Power (dBm)</th><th></th></tr>" +
    all.map((r) => `<tr><td>${r.k}</td><td>${r.f.toFixed(2)}</td><td>${r.p < -250 ? "none" : r.p.toFixed(1)}</td><td style="text-align:left">${r.label}</td></tr>`).join("");
}

function drawMap(m1) {
  // Sensitivity contour is a Cassini oval d1*d2 = K. Size the map to fit it.
  const K = Math.pow(10, (prxDbm(m1, 1, 1) - state.sens) / 20);
  const half = Math.min(Math.max(1.25 * Math.sqrt(state.d0 * state.d0 / 4 + K), state.d0 * 0.9, 1), 500);
  const n = 121, cxm = state.d0 / 2;
  const xs = Array.from({ length: n }, (_, i) => cxm - half + 2 * half * i / (n - 1));
  const ys = Array.from({ length: n }, (_, i) => -half * 0.75 + 1.5 * half * i / (n - 1));
  const z = ys.map((y) => xs.map((x) => {
    const d1 = Math.hypot(x, y), d2 = Math.hypot(x - state.d0, y);
    return prxDbm(m1, Math.max(d1, 0.05), Math.max(d2, 0.05));
  }));
  const zmin = state.sens - 30, zmax = state.sens + 40;
  const ramp = seqRamp();
  const traces = [
    { type: "heatmap", x: xs, y: ys, z, zmin, zmax, colorscale: ramp, colorbar: { title: { text: "dBm", side: "right" }, thickness: 10, outlinewidth: 0, tickfont: { size: 10 } }, hovertemplate: "x %{x:.2f} m, y %{y:.2f} m<br>%{z:.1f} dBm<extra></extra>" },
    { type: "contour", x: xs, y: ys, z, contours: { start: state.sens, end: state.sens, size: 1, coloring: "none" }, line: { color: css("--text-primary"), width: 2 }, showscale: false, hoverinfo: "skip", name: "Sensitivity" },
    { x: [0], y: [0], mode: "markers+text", text: ["Carrier"], textposition: "bottom center", textfont: { color: css("--text-primary") }, marker: { size: 12, symbol: "triangle-up", color: css("--s3"), line: { color: css("--surface"), width: 2 } }, hoverinfo: "skip", showlegend: false },
    { x: [state.d0], y: [0], mode: "markers+text", text: ["Rx"], textposition: "bottom center", textfont: { color: css("--text-primary") }, marker: { size: 12, symbol: "square", color: css("--s2"), line: { color: css("--surface"), width: 2 } }, hoverinfo: "skip", showlegend: false },
  ];
  Plotly.react("p-map", traces, baseLayout({
    xaxis: Object.assign(baseLayout().xaxis, { title: { text: "x (m)" }, range: [xs[0], xs[n - 1]], showgrid: false }),
    yaxis: Object.assign(baseLayout().yaxis, { title: { text: "y (m)" }, scaleanchor: "x", range: [ys[0], ys[n - 1]], showgrid: false }),
    margin: { l: 56, r: 8, t: 8, b: 44 },
  }), PCFG);
}

// ---------- URL state ----------
function writeHash() {
  const q = new URLSearchParams();
  for (const k in PARAMS) if (state[k] !== PARAMS[k].def) q.set(k, state[k]);
  if (state.mode === "z") { q.set("mode", "z"); q.set("z", state.z.slice(0, state.nph).map((z) => `${z.re},${z.im}`).join(";")); }
  history.replaceState(null, "", q.toString() ? "#" + q.toString() : location.pathname);
}

function readHash() {
  for (const k in PARAMS) state[k] = PARAMS[k].def;
  state.mode = "phasor";
  const q = new URLSearchParams(location.hash.slice(1));
  for (const k in PARAMS) if (q.has(k) && isFinite(parseFloat(q.get(k)))) state[k] = parseFloat(q.get(k));
  seedZ();
  if (q.get("mode") === "z") {
    state.mode = "z";
    const zs = (q.get("z") || "").split(";").map((s) => s.split(",").map(Number));
    zs.forEach((p, i) => { if (i < state.z.length && p.length === 2 && p.every(isFinite)) state.z[i] = cx(p[0], p[1]); });
    renderZList();
  }
}

// ---------- boot ----------
function initTheme() {
  const btn = document.getElementById("theme-btn");
  let saved = null;
  try { saved = localStorage.getItem("bsd-theme"); } catch (e) {}
  if (saved) document.documentElement.dataset.theme = saved;
  btn.addEventListener("click", () => {
    const dark = getComputedStyle(document.documentElement).colorScheme.includes("dark");
    const next = dark ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    try { localStorage.setItem("bsd-theme", next); } catch (e) {}
    render();
  });
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => render());
}

buildControls();
readHash();
syncControls();
initTheme();

window.addEventListener("hashchange", () => { readHash(); syncControls(); render(); });

document.getElementById("reset-btn").addEventListener("click", () => {
  for (const k in PARAMS) state[k] = PARAMS[k].def;
  state.mode = "phasor";
  seedZ();
  syncControls();
  schedule();
});
document.getElementById("share-btn").addEventListener("click", async (e) => {
  try { await navigator.clipboard.writeText(location.href); e.target.textContent = "Copied"; }
  catch (err) { e.target.textContent = "Copy failed"; }
  setTimeout(() => (e.target.textContent = "Copy link"), 1500);
});

render();
// Re-measure legends once web fonts arrive.
if (document.fonts) document.fonts.ready.then(() => render());
