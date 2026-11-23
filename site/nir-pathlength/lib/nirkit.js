// 本書共用的顏色、計算與繪圖小工具（由第 1 章抽出；第 2 章起載入此檔）。需先載入 engine.js 與 nirdata.js。
'use strict';
/* ---------- 顏色：同一顏色在全章代表同一件事 ---------- */
const PATHS = NIR.paths;                                   // [0.25, 0.5, 1, 2, 5] mm
const PATHC = { 0.25: '#7ec8e3', 0.5: '#6fd3a0', 1: '#f2c14e', 2: '#f39060', 5: '#e8596b' };   // 光徑
const SAMPC = ['#a8a8a8', '#c9bcff', '#a28eff', '#7f62ff'];  // 水、3.32、4.87、7.45 %
const SAMPL = ['水 0%', '3.32%', '4.87%', '7.45%'];
const GOLD = '#f2c14e', OKC = '#6fd3a0', BADC = '#ff6b6b', WATERC = '#6fa8dc';
const CONC = NIR.conc, WL = NIR.wl;
const A_LIM = 1.5;
const RATIO = [1.00, 1.93, 3.95, 3.38, -0.36];             // 2258 nm 靈敏度比（以 0.25 mm 為 1），與報告一致
const IDEAL = PATHS.map(p => p / 0.25);
const pstr = p => String(p) + ' mm';
const fmt = (v, d = 2) => v.toFixed(d).replace('-', '−');

/* ---------- 由數據計算（與 R / Python 練習同一算法） ---------- */
const nMaxA = p => WL.map((_, i) => Math.max(...NIR.A[String(p)].map(s => s[i])));
const nFrac = p => nMaxA(p).filter(a => a < A_LIM).length / WL.length;
function nSegments(p) {                                    // 4 條光譜皆 A < 1.5 的連續波段
  const ok = nMaxA(p).map(a => a < A_LIM), out = [];
  for (let i = 0; i < ok.length; i++) {
    if (!ok[i]) continue;
    let j = i; while (j + 1 < ok.length && ok[j + 1]) j++;
    if (j > i) out.push([WL[i], WL[j]]);
    i = j;
  }
  return out;
}
function nFit(xs, ys) {                                    // 最小平方法直線
  const n = xs.length, mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0; xs.forEach((x, i) => { sxy += (x - mx) * (ys[i] - my); sxx += (x - mx) ** 2; });
  const slope = sxy / sxx; return { slope, icpt: my - slope * mx };
}
const B1 = NIR.band['2258']['1'];                          // 1 mm、2258 nm
const D2K = B1.d2.map(v => v * 1e3);                       // 二階導數 ×10³
const FIT1 = nFit(CONC, D2K);
const R2_1 = (() => { const my = D2K.reduce((a, b) => a + b) / 4; let r = 0, t = 0;
  D2K.forEach((y, i) => { r += (y - (FIT1.icpt + FIT1.slope * CONC[i])) ** 2; t += (y - my) ** 2; }); return 1 - r / t; })();
const LOO = CONC.map((c, j) => {                           // 留一交叉驗證
  const xs = CONC.filter((_, k) => k !== j), ys = D2K.filter((_, k) => k !== j), f = nFit(xs, ys);
  const pred = (D2K[j] - f.icpt) / f.slope; return { f, pred, err: c - pred };
});
const RMSECV = Math.sqrt(LOO.reduce((a, l) => a + l.err ** 2, 0) / 4);

/* ---------- 繪圖小工具（加 n 前綴，避免和引擎的全域名稱衝突） ---------- */
function nTitle(p, str, t0 = .2) { const e = T(p, str, { x: 110, y: 96, size: 40, weight: 700, o: 0 }); show(e, t0); return e; }
function nText(p, str, x, y, t0, { size = 28, fill = COL.chalk, anchor = 'start', weight = 500 } = {}) {
  const e = T(p, str, { x, y, size, fill, anchor, weight, o: 0 }); show(e, t0); return e;
}
function nPlot(p, box, xr, yr, o = {}) {
  const [bx, by, bw, bh] = box;
  const px = v => bx + (v - xr[0]) / (xr[1] - xr[0]) * bw;
  const py = v => by + bh - (v - yr[0]) / (yr[1] - yr[0]) * bh;
  const g = G(p, { o: 0 });
  const y0 = o.y0 != null ? py(o.y0) : by + bh;
  path(g, `M${bx} ${by}V${by + bh}`, { stroke: COL.dim, 'stroke-width': 2, fill: 'none' });
  path(g, `M${bx} ${y0}H${bx + bw}`, { stroke: COL.dim, 'stroke-width': 2, fill: 'none' });
  (o.xt || []).forEach(v => {
    path(g, `M${px(v)} ${by + bh}v8`, { stroke: COL.dim, 'stroke-width': 2 });
    T(g, o.xf ? o.xf(v) : String(v), { x: px(v), y: by + bh + 32, size: 22, fill: COL.dim, anchor: 'middle' });
  });
  (o.yt || []).forEach(v => {
    path(g, `M${bx} ${py(v)}h-8`, { stroke: COL.dim, 'stroke-width': 2 });
    T(g, o.yf ? o.yf(v) : String(v), { x: bx - 14, y: py(v) + 8, size: 22, fill: COL.dim, anchor: 'end' });
  });
  if (o.xl) T(g, o.xl, { x: bx + bw + 12, y: by + bh + 32, size: 22, fill: COL.dim });
  if (o.yl) T(g, o.yl, { x: bx - 6, y: by - 18, size: 22, fill: COL.dim });
  show(g, o.t0 ?? .2);
  return { g, px, py, box };
}
const nLineD = (P, xs, ys) => xs.map((x, i) => (i ? 'L' : 'M') + P.px(x).toFixed(1) + ' ' + P.py(ys[i]).toFixed(1)).join('');
function nCurve(P, xs, ys, color, t0, dur = 1, w = 3, extra = {}) {
  const e = path(P.g, nLineD(P, xs, ys), { stroke: color, 'stroke-width': w, fill: 'none', 'stroke-linejoin': 'round', ...extra }, { d: 0 });
  draw(e, t0, dur); return e;
}
function nDot(parent, x, y, r, fill, t0, extra = {}) {
  const e = G(parent, { x, y, s: 0 }); mk('circle', { r, fill, ...extra }, e); tw(e, { s: 1 }, t0, .4, back); return e;
}
function nBand(P, x0, x1, fill, o, t0) {                    // 直立色帶（波段標示）
  const [, by, , bh] = P.box;
  const e = mk('rect', { x: P.px(x0), y: by, width: Math.max(2, P.px(x1) - P.px(x0)), height: bh, fill }, P.g);
  put(e, { o: 0 }); tw(e, { o }, t0, .6); return e;
}
function nHLine(P, v, color, t0, dash = '10 8') {
  const [bx, , bw] = P.box;
  const e = path(P.g, `M${bx} ${P.py(v)}H${bx + bw}`, { stroke: color, 'stroke-width': 3, 'stroke-dasharray': dash }, { o: 0 });
  show(e, t0); return e;
}
function nBar(parent, x, w, yBase, yTop, fill, t0, dur = .8, extra = {}) {   // 由基線長出的長條
  const e = mk('rect', { x: x - w / 2, y: yBase, width: w, height: 0, fill, ...extra }, parent);
  const top = Math.min(yBase, yTop), hgt = Math.abs(yBase - yTop);
  tw(e, { a_y: top, a_height: hgt }, t0, dur, out); return e;
}
function nArrowHead(parent, x, y, dir, color, size = 16) {  // dir: 'r' 'd'
  const d = dir === 'r' ? `M${x} ${y}l${-size} ${-size * .6}v${size * 1.2}z` : `M${x} ${y}l${-size * .6} ${-size}h${size * 1.2}z`;
  return path(parent, d, { fill: color });
}
function nLegend(parent, items, x, y, t0, size = 24) {      // items: [[label, color]]
  const g = G(parent, { o: 0 });
  items.forEach(([lab, c], i) => {
    path(g, `M${x} ${y + i * 34 - 8}h30`, { stroke: c, 'stroke-width': 5 });
    T(g, lab, { x: x + 42, y: y + i * 34, size, fill: c });
  });
  show(g, t0); return g;
}

