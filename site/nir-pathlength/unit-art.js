// Chalk sketches for the unit headers on the contents page, one per unit, drawn in the unit's colour.
// Each entry is the inside of an SVG with viewBox 0 0 600 180; keep the left ~200 px light (the unit title sits there
// on narrow screens). Classes, styled in lib/units.css:
//   d  a stroke that draws itself in when the unit scrolls into view (--k orders the strokes)
//   t  text or a filled shape that fades in
//   a  faint chalk white instead of the unit colour;  faint  even fainter
// A sketch may add one looping motion after it is drawn: give an element a class that lib/units.css animates (swap-a,
// hop, orbit, ride, spin, tilt, breathe, ...), or add your own .seen .unit-art .NAME rule and keyframes in a <style> in index.html.
'use strict';
(() => {
  let k = 0;
  const P = (d, cls = '', extra = '') => `<path class="d ${cls}" pathLength="1" style="--k:${k++}" d="${d}" ${extra}/>`;
  const Tx = (x, y, s, size = 24, cls = '', anchor = 'start') => `<text class="t ${cls}" style="--k:${k++}" x="${x}" y="${y}" font-size="${size}" text-anchor="${anchor}">${s}</text>`;
  const C = (cx, cy, r, cls = '') => `<circle class="d ${cls}" pathLength="1" style="--k:${k++}" cx="${cx}" cy="${cy}" r="${r}"/>`;
  const reset = s => { k = 0; return s; };

  window.UNIT_ART = [
    // 1 光徑選擇：三條吸收光譜，最長的那條在頂端被壓平（飽和），虛線是 A = 1.5
    reset(P('M228 160H580', 'a') + P('M228 160V20', 'a') +
      P('M228 150C300 145 330 120 360 128S430 150 470 120S540 110 580 116') +
      P('M228 140C300 130 330 80 360 95S430 140 470 80S540 60 580 70') +
      P('M228 128C290 110 320 34 360 34H400C420 34 440 120 470 34H580', 'a') +
      P('M228 72H580', 'faint', 'stroke-dasharray="6 6"') + Tx(500, 64, 'A = 1.5', 18, 'a') + Tx(300, 28, '飽和', 20, 'a')),
  ];
})();
