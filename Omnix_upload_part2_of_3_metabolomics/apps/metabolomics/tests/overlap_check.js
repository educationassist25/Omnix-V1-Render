// Geometric no-overlap check for the rewiring network SVG (run inside the page with frame.evaluate).
// Every <text> is turned into its (possibly rotated) box; the check fails when a box
//   * intersects another text box,
//   * touches a node dot, or
//   * is crossed by any drawn line (rewired edges, reaction edges, pathway arcs, leader lines).
// Returns a list of human-readable problems (empty = clean).
module.exports = function overlapProblems() {
  const svg = document.getElementById('net');
  const M = (el, x, y) => { const m = el.getCTM(); return { x: m.a * x + m.c * y + m.e, y: m.b * x + m.d * y + m.f }; };
  const texts = [...svg.querySelectorAll('text')].map(t => {
    const b = t.getBBox(); const i = 0.6;   // glyph cells include a little line gap
    const pts = [[b.x + i, b.y + i], [b.x + b.width - i, b.y + i], [b.x + b.width - i, b.y + b.height - i], [b.x + i, b.y + b.height - i]].map(([x, y]) => M(t, x, y));
    const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
    return { name: t.textContent.slice(0, 40), pts, x0: Math.min(...xs), x1: Math.max(...xs), y0: Math.min(...ys), y1: Math.max(...ys) };
  });
  // coarse grid of text boxes so each line sample is only tested against nearby text
  const CELL = 40, grid = new Map();
  texts.forEach(t => { for (let gx = Math.floor(t.x0 / CELL); gx <= Math.floor(t.x1 / CELL); gx++) for (let gy = Math.floor(t.y0 / CELL); gy <= Math.floor(t.y1 / CELL); gy++) { const k = gx + ',' + gy; if (!grid.has(k)) grid.set(k, []); grid.get(k).push(t); } });
  const near = p => grid.get(Math.floor(p.x / CELL) + ',' + Math.floor(p.y / CELL)) || [];
  const axes = P => P.map((p, k) => { const q = P[(k + 1) % P.length]; return { x: -(q.y - p.y), y: q.x - p.x }; });
  const proj = (P, a) => { let lo = Infinity, hi = -Infinity; P.forEach(p => { const v = p.x * a.x + p.y * a.y; lo = Math.min(lo, v); hi = Math.max(hi, v); }); return [lo, hi]; };
  const polysOverlap = (A, B) => axes(A).concat(axes(B)).every(a => { const [a0, a1] = proj(A, a), [b0, b1] = proj(B, a); return a1 > b0 && b1 > a0; }) ;
  const segDist = (p, a, b) => { const dx = b.x - a.x, dy = b.y - a.y; const L = dx * dx + dy * dy; let t = L ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / L : 0; t = Math.max(0, Math.min(1, t)); return Math.hypot(p.x - a.x - t * dx, p.y - a.y - t * dy); };
  const inside = (p, P) => { let s = 0; for (let k = 0; k < P.length; k++) { const a = P[k], b = P[(k + 1) % P.length]; const c = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x); if (c !== 0) { if (s === 0) s = Math.sign(c); else if (Math.sign(c) !== s) return false; } } return true; };
  const distPoly = (p, P) => inside(p, P) ? 0 : Math.min(...P.map((a, k) => segDist(p, a, P[(k + 1) % P.length])));
  const out = [];
  for (let a = 0; a < texts.length; a++) for (let b = a + 1; b < texts.length; b++)
    if (texts[a].x1 > texts[b].x0 && texts[b].x1 > texts[a].x0 && texts[a].y1 > texts[b].y0 && texts[b].y1 > texts[a].y0 && polysOverlap(texts[a].pts, texts[b].pts)) out.push(`text/text: ${texts[a].name} | ${texts[b].name}`);
  svg.querySelectorAll('circle').forEach(c => {
    const m = c.getCTM(); const p = M(c, +c.getAttribute('cx'), +c.getAttribute('cy')); const r = (+c.getAttribute('r') + 1) * Math.hypot(m.a, m.b);
    texts.forEach(t => { if (p.x > t.x0 - r && p.x < t.x1 + r && p.y > t.y0 - r && p.y < t.y1 + r && distPoly(p, t.pts) < r) out.push(`dot/text: ${t.name}`); });
  });
  svg.querySelectorAll('path').forEach(pth => {
    if (pth.classList.contains('hit') || pth.closest('defs')) return;   // invisible hit targets, marker shapes
    const sw = parseFloat(getComputedStyle(pth).strokeWidth) || 1; const m = pth.getCTM(); const half = sw / 2 * Math.hypot(m.a, m.b);
    const L = pth.getTotalLength(); const n = Math.max(2, Math.ceil(L / 2));
    for (let k = 0; k <= n; k++) {
      const q = pth.getPointAtLength(L * k / n); const p = M(pth, q.x, q.y);
      const hit = near(p).find(t => p.x > t.x0 - half && p.x < t.x1 + half && p.y > t.y0 - half && p.y < t.y1 + half && distPoly(p, t.pts) < half);
      if (hit) { out.push(`line(${pth.getAttribute('class') || 'path'})/text: ${hit.name}`); break; }
    }
  });
  svg.querySelectorAll('line').forEach(ln => {
    const sw = parseFloat(ln.getAttribute('stroke-width')) || 1; const m = ln.getCTM(); const half = sw / 2 * Math.hypot(m.a, m.b);
    const a = M(ln, +ln.getAttribute('x1'), +ln.getAttribute('y1')), b = M(ln, +ln.getAttribute('x2'), +ln.getAttribute('y2'));
    const n = Math.max(2, Math.ceil(Math.hypot(b.x - a.x, b.y - a.y) / 2));
    for (let k = 0; k <= n; k++) { const p = { x: a.x + (b.x - a.x) * k / n, y: a.y + (b.y - a.y) * k / n };
      const hit = texts.find(t => distPoly(p, t.pts) < half); if (hit) { out.push(`line/text: ${hit.name}`); break; } }
  });
  svg.querySelectorAll('rect').forEach(r => {
    const x = +r.getAttribute('x'), y = +r.getAttribute('y'), w = +r.getAttribute('width'), h = +r.getAttribute('height');
    const R = [M(r, x, y), M(r, x + w, y), M(r, x + w, y + h), M(r, x, y + h)];
    texts.forEach(t => { if (polysOverlap(R, t.pts)) out.push(`rect/text: ${t.name}`); });
  });
  return out;
};
