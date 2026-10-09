'use strict';
/* VisRWE explorer. Data: BUNDLE (built by pipeline/09_explorer_data.py). Rendering: vis-network 9.1.9. */

// ------------------------------------------------------------------------------------------------ utilities
const $ = (id) => document.getElementById(id);
const fromB64 = (s) => Uint8Array.from(atob(s), (c) => c.charCodeAt(0));
const fmt = (v, d = 2) => (v == null || Number.isNaN(v) ? '—' : (v >= 0 ? '+' : '') + v.toFixed(d));
const fmtP = (p) => (p == null || Number.isNaN(p) ? '—' : p < 1e-3 ? p.toExponential(1) : p.toPrecision(2));
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const cssVar = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const store = {
  get(k, d) { try { const v = localStorage.getItem('visrwe.' + k); return v == null ? d : v; } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem('visrwe.' + k, v); } catch (e) { /* storage unavailable */ } },
};

async function gunzipB64(s) {
  const stream = new Blob([fromB64(s)]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Int8Array(await new Response(stream).arrayBuffer());
}

// Student t two-sided p via the regularised incomplete beta function (Numerical Recipes)
function lgamma(z) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = z, t = z + 5.5; t -= (z + 0.5) * Math.log(t); let s = 1.000000000190015;
  for (let j = 0; j < 6; j++) s += c[j] / ++y;
  return -t + Math.log(2.5066282746310005 * s / z);
}
function betacf(a, b, x) {
  const TINY = 1e-300; let c = 1, d = 1 - (a + b) * x / (a + 1); if (Math.abs(d) < TINY) d = TINY; d = 1 / d; let h = d;
  for (let m = 1; m <= 300; m++) {
    const m2 = 2 * m; let aa = m * (b - m) * x / ((a - 1 + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < TINY) d = TINY; c = 1 + aa / c; if (Math.abs(c) < TINY) c = TINY; d = 1 / d; h *= d * c;
    aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + 1 + m2));
    d = 1 + aa * d; if (Math.abs(d) < TINY) d = TINY; c = 1 + aa / c; if (Math.abs(c) < TINY) c = TINY; d = 1 / d;
    const del = d * c; h *= del; if (Math.abs(del - 1) < 3e-14) break;
  }
  return h;
}
function betai(a, b, x) {
  if (x <= 0) return 0; if (x >= 1) return 1;
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
}
const pFromR = (r, n) => (n < 4 || Number.isNaN(r) ? NaN : betai((n - 2) / 2, 0.5, (n - 2) / (n - 2 + r * r * (n - 2) / (1 - r * r + 1e-15))));
function chi2sf(x, k) { // upper tail of chi-square via the regularised gamma function
  const a = k / 2, z = x / 2; if (z <= 0) return 1;
  if (z < a + 1) { let ap = a, sum = 1 / a, del = sum; for (let n = 0; n < 500; n++) { ap += 1; del *= z / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-15) break; }
    return Math.max(0, 1 - sum * Math.exp(-z + a * Math.log(z) - lgamma(a))); }
  let b = z + 1 - a, c = 1e300, d = 1 / b, h = d;
  for (let i = 1; i < 500; i++) { const an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300; c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-15) break; }
  return Math.exp(-z + a * Math.log(z) - lgamma(a)) * h;
}
const atanh = (r) => 0.5 * Math.log((1 + r) / (1 - r));
const clampR = (r) => Math.max(-0.999, Math.min(0.999, r));

// ------------------------------------------------------------------------------------------------ data model
const D = BUNDLE, G = D.genes, NGENES = G.sym.length;
const E = {
  a: new Uint16Array(fromB64(D.edges.a).buffer), b: new Uint16Array(fromB64(D.edges.b).buffer),
  rng: new Int8Array(fromB64(D.edges.rng).buffer), r2d: new Int8Array(fromB64(D.edges.r2d).buffer), f: fromB64(D.edges.f), n: D.edges.n,
};
const ADJ = (() => { // CSR adjacency: edge indices per gene
  const deg = new Int32Array(NGENES + 1);
  for (let e = 0; e < E.n; e++) { deg[E.a[e] + 1]++; deg[E.b[e] + 1]++; }
  for (let i = 0; i < NGENES; i++) deg[i + 1] += deg[i];
  const fill = deg.slice(0, NGENES), idx = new Int32Array(2 * E.n);
  for (let e = 0; e < E.n; e++) { idx[fill[E.a[e]]++] = e; idx[fill[E.b[e]]++] = e; }
  return { off: deg, idx };
})();
const SYM = new Map(G.sym.map((s, i) => [s.toUpperCase(), i]));
const PAIR = new Map(); for (let e = 0; e < E.n; e++) PAIR.set(E.a[e] * 65536 + E.b[e], e);
const pairEdge = (i, j) => PAIR.get(Math.min(i, j) * 65536 + Math.max(i, j));
const otherEnd = (e, i) => (E.a[e] === i ? E.b[e] : E.a[e]);
const LP = D.lineage_pairs.map(([a, b, lin, r, q, rr, cls]) => ({ a, b, lin, r, q, rr, cls, e: pairEdge(a, b) }));
const LP_EDGES = new Set(LP.map((x) => x.e));
const DATASETS = ['NextGen', 'DepMap2D'];
const VEC = {};
async function vectors(ds) { if (!VEC[ds]) VEC[ds] = gunzipB64(D.vectors[ds]); return VEC[ds]; }
function geneVector(V, ds, i) {
  const n = D.screens[ds].id.length, out = new Float64Array(n);
  for (let s = 0; s < n; s++) { const v = V[i * n + s]; out[s] = v === -128 ? NaN : v / 40; }
  return out;
}

// ------------------------------------------------------------------------------------------------ state, i18n
const S = {
  ds: 'NextGen', mode: 'explore', lang: store.get('lang', (navigator.language || 'zh').startsWith('zh') ? 'zh' : 'en'),
  rmin: 0.3, pos: true, neg: true, cls: { both: true, ng: true, d2: true }, rweOnly: false, hideLib: false,
  labels: 'auto', labelSize: 12, edgeLabels: false, lpHalo: true, maxExpand: 30, modLevel: 'category',
  shown: new Set(), history: [], ghosts: [], hoverNode: null, selected: null, demean: false, hiddenLin: new Set(),
  root: null, baseLen: 0, userMoved: false, builtMode: null, // explore: centre gene, initial log length, view moved by hand, mode drawn last
};
const t = (k, ...a) => { const v = (I18N[S.lang] || I18N.en)[k]; return typeof v === 'function' ? v(...a) : (v == null ? k : v); };
const lab = (type, v) => (S.lang === 'zh' && LABELS[type] && LABELS[type].zh && LABELS[type].zh[v]) || v;
const dsName = (ds) => (ds === 'NextGen' ? t('ds_ng') : t('ds_2d'));
function applyI18n() {
  document.documentElement.lang = S.lang;
  document.querySelectorAll('[data-t]').forEach((el) => { el.textContent = t(el.dataset.t); });
  $('search').placeholder = t('search_ph');
  $('langBtn').textContent = t('lang');
  $('modeHelp').textContent = t(S.mode + '_help');
}

// ------------------------------------------------------------------------------------------------ colours
const FAM_PALETTE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#56B4E9', '#E69F00', '#8C564B', '#7F7F7F', '#BCBD22', '#17BECF', '#9467BD', '#E377C2', '#4D4D4D'];
const CAT_PALETTE = ['#8C564B', '#C49C94', '#9EDAE5', '#AEC7E8', '#C5B0D5', '#F7B6D2', '#BCBD22', '#98DF8A', '#2CA02C', '#B8860B', '#DAA520', '#DBDB8D', '#E377C2', '#1F77B4', '#17BECF',
  '#9467BD', '#C7C7C7', '#D62728', '#FF7F0E'];
let COL = {};
function readColours() { COL = { ng: cssVar('--ng'), d2: cssVar('--d2'), both: cssVar('--both'), ghost: cssVar('--ghost'), ink: cssVar('--ink'), canvas: cssVar('--canvas'), muted: cssVar('--muted') }; }
const nodeColour = (i) => (G.rwe[i] ? FAM_PALETTE[Math.max(G.fam[i], 0) % FAM_PALETTE.length] : CAT_PALETTE[Math.max(G.cat[i], 0) % CAT_PALETTE.length]);

// ------------------------------------------------------------------------------------------------ edge semantics
function edgeClass(e) { const c = E.f[e] & 3; return c === 3 ? 'both' : c === 1 ? 'ng' : c === 2 ? 'd2' : ''; }
function edgeR(e, ds = S.ds) {
  if (ds === 'NextGen') return E.rng[e] / 100;
  if (ds === 'DepMap2D') return E.r2d[e] / 100;
  return (E.f[e] & 1) ? E.rng[e] / 100 : E.r2d[e] / 100;
}
function edgeVisible(e) {
  const f = E.f[e];
  if (S.lpHalo && S.ds !== 'DepMap2D' && LP_EDGES.has(e)) return true;
  if (S.ds === 'NextGen' && !(f & 1)) return false;
  if (S.ds === 'DepMap2D' && !(f & 2)) return false;
  if (S.ds === 'compare') { const c = edgeClass(e); if (!c || !S.cls[c]) return false; }
  const r = edgeR(e);
  if (Math.abs(r) < S.rmin || (r > 0 && !S.pos) || (r < 0 && !S.neg)) return false;
  if (S.hideLib && (f & 4) && S.ds !== 'NextGen') return false;
  if (S.rweOnly && !(G.rwe[E.a[e]] && G.rwe[E.b[e]])) return false;
  return true;
}
function edgeColour(e) {
  if (S.ds === 'NextGen') return COL.ng;
  if (S.ds === 'DepMap2D') return COL.d2;
  return COL[edgeClass(e)] || COL.muted;
}
function visibleNeighbours(i) {
  const out = [];
  for (let k = ADJ.off[i]; k < ADJ.off[i + 1]; k++) { const e = ADJ.idx[k]; if (edgeVisible(e)) out.push({ j: otherEnd(e, i), e, r: edgeR(e) }); }
  return out.sort((x, y) => Math.abs(y.r) - Math.abs(x.r));
}

// ------------------------------------------------------------------------------------------------ network
let net, nodesDS, edgesDS;
function labelFor(i) {
  if (S.labels === 'none') return '';
  if (S.labels === 'all' || G.rwe[i] || i === S.hoverNode || i === S.selected) return G.sym[i];
  return S.mode === 'explore' && S.shown.size <= 120 ? G.sym[i] : '';
}
function nodeObj(i, x, y) {
  const deg = S.ds === 'DepMap2D' ? G.d2d[i] : G.dng[i], root = S.mode === 'explore' && i === S.root;
  const size = (G.rwe[i] ? 14 + Math.min(10, Math.sqrt(deg)) : 5 + Math.min(9, Math.sqrt(deg) * 0.8)) + (root ? 8 : 0);
  const col = nodeColour(i);
  const o = { id: i, label: labelFor(i), size, shape: 'dot',
    color: { background: col, border: root ? COL.ink : G.cav[i] ? '#222' : col, highlight: { background: col, border: COL.ink }, hover: { background: col, border: COL.ink } },
    borderWidth: root ? 3.5 : G.rwe[i] ? (G.cav[i] ? 2.5 : 1.5) : 0.5,
    font: { size: root ? S.labelSize + 5 : G.rwe[i] ? S.labelSize + 1 : S.labelSize - 1, color: COL.ink, strokeWidth: 3, strokeColor: COL.canvas, face: 'inherit' } };
  if (x != null) { o.x = x; o.y = y; }
  return o;
}
function edgeObj(e) {
  const r = edgeR(e), halo = S.lpHalo && LP_EDGES.has(e) && S.ds !== 'DepMap2D', light = S.mode === 'overview';
  return { id: 'e' + e, from: E.a[e], to: E.b[e], width: 0.5 + 6 * Math.max(0, Math.abs(r) - 0.25), dashes: r < 0,
    color: { color: edgeColour(e), opacity: light ? 0.32 : 0.75, highlight: edgeColour(e), hover: edgeColour(e) },
    label: S.edgeLabels ? r.toFixed(2) : '', font: { size: 9, color: COL.muted, strokeWidth: 2, strokeColor: COL.canvas, align: 'middle' },
    shadow: halo ? { enabled: true, color: '#D55E00', size: 14, x: 0, y: 0 } : false, smooth: false };
}
function createNetwork() {
  nodesDS = new vis.DataSet(); edgesDS = new vis.DataSet();
  net = new vis.Network($('net'), { nodes: nodesDS, edges: edgesDS }, {
    autoResize: true, layout: { improvedLayout: false }, physics: { enabled: false },
    interaction: { hover: true, hoverConnectedEdges: false, tooltipDelay: 0, hideEdgesOnDrag: true, multiselect: false, navigationButtons: false },
    nodes: { scaling: { label: { enabled: true, drawThreshold: 5 } } }, edges: { smooth: false, selectionWidth: 1.5 } });
  net.on('hoverNode', (p) => onHoverNode(p.node, p.pointer));
  net.on('blurNode', () => onBlurNode());
  net.on('hoverEdge', (p) => onHoverEdge(p.edge, p.pointer));
  net.on('blurEdge', () => hideTip());
  net.on('click', (p) => onClick(p));
  net.on('dragStart', (p) => { hideTip(); finishAnim(); if (!p.nodes.length) S.userMoved = true; });
  net.on('zoom', () => { hideTip(); finishAnim(); S.userMoved = true; });
  net.on('resize', () => { if (S.mode === 'explore' && S.root != null && !S.userMoved && !anim) animateTo(null, null, egoView(posMap()), 0); });
  $('net').addEventListener('mouseleave', () => { hideTip(); S.hoverNode = null; clearTimeout(ghostTimer); ghostTimer = setTimeout(clearGhosts, 300); });
}

function render() {
  hideTip(); clearGhosts(); finishAnim();
  $('net').style.display = S.mode === 'modules' ? 'none' : ''; $('matrix').style.display = S.mode === 'modules' ? 'block' : 'none';
  $('empty').style.display = S.mode === 'explore' && S.root == null ? 'flex' : 'none';
  if (S.mode === 'explore') renderExplore();
  else { nodesDS.clear(); edgesDS.clear(); S.builtMode = S.mode; if (S.mode === 'overview') renderOverview(); else renderModules(); }
  updateCounts(); updateLegend(); updateCrumbs(); updateEgoBar();
}

// ---- camera and animation: node positions and the view are driven together by one requestAnimationFrame loop (no physics, no vis focus)
const posMap = () => { const P = net.getPositions(), m = new Map(); for (const k in P) if (k[0] !== 'g') m.set(+k, P[k]); return m; };
function viewBox() { // part of the canvas not covered by the detail panel (DOM px) and the matching horizontal offset of its centre
  const el = $('net'), W = el.clientWidth || 800, H = el.clientHeight || 600, p = $('panel');
  if (p.style.display === 'block') { const pw = p.offsetWidth + 20; if (W - pw >= 0.4 * W) return { W: W - pw, H, off: -pw / 2 }; }
  return { W, H, off: 0 };
}
function egoView(pos) { // the explored gene sits at the centre of the visible canvas; every node fits around it, with room left for a hover preview
  const v = viewBox(), c = pos.get(S.root) || { x: 0, y: 0 }, pad = 110;
  let mx = 1, my = 1;
  for (const p of pos.values()) { mx = Math.max(mx, Math.abs(p.x - c.x)); my = Math.max(my, Math.abs(p.y - c.y)); }
  return { pos: { x: c.x, y: c.y }, off: v.off, scale: Math.min(1, Math.max(0.06, Math.min((v.W / 2 - 30) / (mx + pad), (v.H / 2 - 20) / (my + pad)))) };
}
function fitView(pos) { // bounding-box fit (overview)
  const v = viewBox(); let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const p of pos.values()) { x0 = Math.min(x0, p.x); x1 = Math.max(x1, p.x); y0 = Math.min(y0, p.y); y1 = Math.max(y1, p.y); }
  if (!Number.isFinite(x0)) return null;
  return { pos: { x: (x0 + x1) / 2, y: (y0 + y1) / 2 }, off: v.off, scale: Math.min(1.2, Math.max(0.02, Math.min((v.W - 60) / (x1 - x0 + 60), (v.H - 60) / (y1 - y0 + 60)))) };
}
let anim = null;
function animateTo(target, start, view, ms, done) { // target/start: Map id -> {x, y}; view: {pos, scale, off} or null (camera stays); done: called once at the end
  finishAnim();
  if (matchMedia('(prefers-reduced-motion: reduce)').matches) ms = 0;
  const B = net.body.nodes, t0 = performance.now(), me = { target, raf: 0, done }, cv = net.canvas.frame.canvas;
  // start = the canvas point now under the target's screen anchor, so a gene that is already at the anchor stays put while the scale changes
  const v0 = view && { pos: net.DOMtoCanvas({ x: cv.clientWidth / 2 + view.off, y: cv.clientHeight / 2 }), scale: net.getScale() };
  const ease = (u) => 1 - (1 - u) ** 3;
  const frame = (now) => {
    const u = ms ? Math.min(1, (now - t0) / ms) : 1, k = ease(u);
    if (target) for (const [id, b] of target) { const n = B[id]; if (!n) continue; const a = (start && start.get(id)) || b; n.x = a.x + (b.x - a.x) * k; n.y = a.y + (b.y - a.y) * k; }
    if (view) net.moveTo({ position: { x: v0.pos.x + (view.pos.x - v0.pos.x) * k, y: v0.pos.y + (view.pos.y - v0.pos.y) * k },
      scale: v0.scale + (view.scale - v0.scale) * k, offset: { x: view.off, y: 0 }, animation: false });
    net.redraw();
    if (u < 1) { me.raf = requestAnimationFrame(frame); return; }
    if (anim === me) anim = null;
    const d = me.done; me.done = null; if (d) d();
  };
  anim = me;
  if (ms) me.raf = requestAnimationFrame(frame); else frame(t0);
}
function finishAnim() { // jump a running animation to its final node positions; the camera stays where it is
  if (!anim) return;
  const a = anim; anim = null; cancelAnimationFrame(a.raf);
  if (a.target) { for (const [id, b] of a.target) { const n = net.body.nodes[id]; if (n) { n.x = b.x; n.y = b.y; } } net.redraw(); }
  const d = a.done; a.done = null; if (d) d();
}
function recentre(ms) { if (S.mode === 'explore' && S.root != null) animateTo(null, null, egoView(posMap()), ms); }
function centreOn(i, ms) { const p = net.getPositions([i])[i]; if (p) animateTo(null, null, { pos: p, scale: Math.max(net.getScale(), 0.8), off: viewBox().off }, ms); }

// ---- explore mode: the network of one gene (centre) and its called partners; partners can be expanded outwards
// S.history = expansion log [{node, kind: 'expand' | 'add', added}]; the first S.baseLen entries make the initial network of S.root
const absR = (i, j) => { const e = pairEdge(i, j); return e == null ? 0 : Math.abs(edgeR(e)); };
const strength = (i, j) => Math.max(0, Math.min(1, (absR(i, j) - 0.3) / 0.5)); // 0 at |r| = 0.3, 1 at |r| >= 0.8
const groupKey = (j) => (G.rwe[j] ? -100 + Math.max(G.fam[j], 0) : Math.max(G.cat[j], 0));
const egoOrder = (i, list) => [...list].sort((a, b) => groupKey(a) - groupKey(b) || absR(i, b) - absR(i, a)); // same colour together, stronger first
const SPACING = 40;
const topNew = (i) => visibleNeighbours(i).filter((x) => !S.shown.has(x.j)).slice(0, S.maxExpand).map((x) => x.j);
function showEgo(i, extra = []) { // entry point: search, URL, "centre on this gene", lineage-pair list
  S.mode = 'explore'; syncSeg('modeSeg', 'explore'); $('modLevelRow').style.display = 'none'; $('modeHelp').textContent = t('explore_help');
  S.root = i; S.history = extra.length ? [{ node: i, kind: 'add', added: extra }] : []; // the gene alone; the user grows the network level by level
  S.baseLen = S.history.length; S.userMoved = false;
  $('net').style.display = ''; $('matrix').style.display = 'none'; $('empty').style.display = 'none';
  selectGene(i); writeHash();
  buildEgo(true, true);
  updateLegend();
}
function renderExplore() {
  if (S.root == null) { nodesDS.clear(); edgesDS.clear(); S.shown = new Set(); S.history = []; S.builtMode = 'explore'; return; }
  buildEgo(false, true);
}
function egoReplay() { // re-derive the shown genes from the expansion log under the current dataset and filters
  const log = S.history; S.shown = new Set([S.root]); S.history = [];
  for (const h of log) {
    if (!S.shown.has(h.node)) continue;
    const added = h.kind === 'add' ? h.added.filter((j) => { const e = pairEdge(h.node, j); return !S.shown.has(j) && e != null && edgeVisible(e); }) : topNew(h.node);
    added.forEach((j) => S.shown.add(j));
    S.history.push({ node: h.node, kind: h.kind, added });
  }
}
function egoLayout() { // radial tree: root at the origin, each subtree in its own angular wedge, closer to the parent = larger |r|
  const root = S.root, parent = new Map(), kids = new Map();
  for (const h of S.history) for (const j of h.added) if (j !== root && S.shown.has(j) && !parent.has(j)) parent.set(j, h.node);
  for (const j of S.shown) if (j !== root && !parent.has(j)) parent.set(j, root);
  for (const j of parent.keys()) { const seen = new Set([j]); let q = parent.get(j); while (q !== root) { if (!S.shown.has(q) || seen.has(q)) { parent.set(j, root); break; } seen.add(q); q = parent.get(q); } }
  for (const [j, p] of parent) { if (!kids.has(p)) kids.set(p, []); kids.get(p).push(j); }
  for (const [p, list] of kids) kids.set(p, egoOrder(p, list));
  const w = new Map(), weigh = (v) => { let s = 0; for (const c of kids.get(v) || []) s += weigh(c); s = Math.max(1, 0.55 * s); w.set(v, s); return s; };
  weigh(root);
  const pos = new Map([[root, { x: 0, y: 0 }]]), wedge = new Map([[root, [-Math.PI / 2, 1.5 * Math.PI]]]);
  let level = [root], R = 0, band = 0;
  for (let depth = 1; level.length; depth++) {
    const next = [];
    for (const v of level) {
      const ch = kids.get(v) || []; if (!ch.length) continue;
      const [a0, a1] = wedge.get(v), tot = ch.reduce((s, c) => s + w.get(c), 0); let a = a0;
      for (const c of ch) { const da = (a1 - a0) * w.get(c) / tot; wedge.set(c, [a, a + da]); a += da; next.push(c); }
    }
    if (!next.length) break;
    const minW = Math.min(...next.map((c) => wedge.get(c)[1] - wedge.get(c)[0]));
    R = Math.max(depth === 1 ? 170 : R * (1 + band) + 130, SPACING / minW); band = depth === 1 ? 0.45 : 0.15;
    for (const c of next) { const [a0, a1] = wedge.get(c), th = (a0 + a1) / 2, d = R * (1 + band * (1 - strength(parent.get(c), c))); pos.set(c, { x: d * Math.cos(th), y: d * Math.sin(th) }); }
    level = next;
  }
  return { pos, parent };
}
function buildEgo(animate, forceView) { // (re)draw the explored network from the log; genes already on screen move to their new spots
  hideTip(); clearGhosts(); finishAnim();
  const old = S.builtMode === 'explore' ? net.getPositions() : {};
  egoReplay();
  const { pos, parent } = egoLayout(), start = new Map();
  for (const [i, p] of pos) {
    let s = old[i] || old['g' + i];
    if (!s) { let q = parent.get(i); while (q != null && !old[q]) q = parent.get(q); s = (q != null && old[q]) || { x: 0, y: 0 }; }
    start.set(i, animate ? { x: s.x, y: s.y } : p);
  }
  nodesDS.clear(); edgesDS.clear();
  nodesDS.add([...S.shown].map((i) => nodeObj(i, start.get(i).x, start.get(i).y)));
  edgesDS.add(inducedEdges(S.shown).map(edgeObj));
  S.builtMode = 'explore';
  animateTo(pos, start, forceView || !S.userMoved ? egoView(pos) : null, animate ? 520 : 0);
  updateCounts(); updateCrumbs(); updateEgoBar();
}
function inducedEdges(set) {
  const out = [];
  for (const i of set) for (let k = ADJ.off[i]; k < ADJ.off[i + 1]; k++) { const e = ADJ.idx[k], j = otherEnd(e, i); if (j > i && set.has(j) && edgeVisible(e)) out.push(e); }
  return out;
}
function placeAround(i, list) { // where new partners of node i go: a fan pointing away from the centre (a full outer ring for the centre itself)
  const P = posMap(), c = P.get(i), out = new Map(), n = list.length;
  if (!c || !n) return out;
  const root = P.get(S.root) || { x: 0, y: 0 }, occ = [...P.values()], items = egoOrder(i, list);
  let ang, dist;
  if (i === S.root) { // first level: the same ring as egoLayout; later rounds: an outer ring
    let maxD = 0; for (const p of occ) maxD = Math.max(maxD, Math.hypot(p.x - c.x, p.y - c.y));
    const first = maxD < 1, R = Math.max(first ? 170 : maxD + 110, n * SPACING / (2 * Math.PI)), band = first ? 0.45 : 0.2;
    ang = (k) => -Math.PI / 2 + 2 * Math.PI * (k + 0.5) / n; dist = (j) => R * (1 + band * (1 - strength(i, j)));
  } else { // larger fans alternate between two staggered rows, which halves the arc they need
    const rows = n > 6 ? 2 : 1, th = Math.atan2(c.y - root.y, c.x - root.x), d0 = Math.max(100 + 1.5 * n, (n - 1) * SPACING / rows / Math.PI);
    const span = n > 1 ? Math.min(Math.PI, (n - 1) * SPACING / rows / d0) : 0;
    ang = (k) => th + (n > 1 ? span * (k / (n - 1) - 0.5) : 0); dist = (j, k) => d0 * (1 + 0.25 * (1 - strength(i, j))) + (rows === 2 && k % 2 ? 46 : 0);
  }
  items.forEach((j, k) => { // step outwards along the ray until the spot is free
    const a = ang(k); let d = dist(j, k), p;
    for (let step = 0; step < 40; step++) { p = { x: c.x + d * Math.cos(a), y: c.y + d * Math.sin(a) }; if (!occ.some((q) => Math.hypot(q.x - p.x, q.y - p.y) < 30)) break; d += 16; }
    occ.push(p); out.set(j, p);
  });
  return out;
}
function addGenes(list, target, from, start) { // add genes at their final spots, growing out of `from`; the centre gene stays centred
  if (!list.length) return;
  finishAnim();
  const P = posMap(), c = P.get(from) || P.get(S.root) || { x: 0, y: 0 };
  start = start || new Map(list.map((j) => [j, { x: c.x, y: c.y }]));
  list.forEach((j) => S.shown.add(j));
  nodesDS.add(list.map((j) => nodeObj(j, start.get(j).x, start.get(j).y)));
  const have = new Set(edgesDS.getIds()), add = [];
  for (const j of list) for (let k = ADJ.off[j]; k < ADJ.off[j + 1]; k++) { const e = ADJ.idx[k]; if (S.shown.has(otherEnd(e, j)) && edgeVisible(e) && !have.has('e' + e)) { add.push(edgeObj(e)); have.add('e' + e); } }
  edgesDS.add(add);
  const all = posMap(); for (const [j, p] of target) all.set(j, p);
  animateTo(target, start, S.userMoved ? null : egoView(all), 480);
  refreshLabels(); updateCounts(); updateEgoBar();
}
function expand(i) {
  if (S.mode !== 'explore' || !S.shown.has(i)) { showEgo(i); return; }
  const add = topNew(i);
  if (!add.length) return;
  S.history.push({ node: i, kind: 'expand', added: add });
  addGenes(add, placeAround(i, add), i);
  updateCrumbs();
}
const isExpanded = (i) => S.history.slice(S.baseLen).some((h) => h.node === i && h.added.some((j) => S.shown.has(j)));
function retract(list, into) { // fold genes back into the node they grew from, then drop them; the view refits what is left
  finishAnim(); clearGhosts();
  const rm = list.filter((j) => S.shown.has(j) && j !== S.root);
  if (!rm.length) return;
  const P = posMap(), c = P.get(into) || P.get(S.root) || { x: 0, y: 0 };
  rm.forEach((j) => S.shown.delete(j));
  const rest = new Map([...P].filter(([j]) => S.shown.has(j)));
  animateTo(new Map(rm.map((j) => [j, { x: c.x, y: c.y }])), new Map(rm.map((j) => [j, P.get(j) || c])), S.userMoved ? null : egoView(rest), 380, () => {
    nodesDS.remove(rm);
    edgesDS.remove(edgesDS.getIds().filter((id) => { const e = +String(id).slice(1); return String(id)[0] === 'e' && (!S.shown.has(E.a[e]) || !S.shown.has(E.b[e])); }));
    refreshLabels(); updateCounts();
  });
  updateCrumbs(); updateEgoBar();
}
function collapse(i) { // fold back everything grown out of node i (its expansions and their descendants); for the centre gene: all of it
  let gone;
  if (i === S.root) {
    const keep = new Set([S.root, ...S.history.slice(0, S.baseLen).flatMap((h) => h.added)]);
    S.history = S.history.slice(0, S.baseLen);
    gone = [...S.shown].filter((j) => !keep.has(j));
  } else {
    const drop = new Set([i]), out = new Set();
    for (let grew = true; grew;) { grew = false; for (const h of S.history.slice(S.baseLen)) if (drop.has(h.node)) for (const j of h.added) if (!out.has(j)) { out.add(j); drop.add(j); grew = true; } }
    S.history = S.history.filter((h, k) => k < S.baseLen || !drop.has(h.node));
    const still = new Set(S.history.flatMap((h) => h.added));
    gone = [...out].filter((j) => !still.has(j));
  }
  retract(gone, i);
}
function undo() {
  if (S.mode !== 'explore' || S.history.length <= S.baseLen) return;
  const h = S.history.pop(), still = new Set(S.history.flatMap((x) => x.added));
  retract(h.added.filter((j) => !still.has(j)), h.node);
}
function resetEgo() { if (S.root == null) return; S.history = S.history.slice(0, S.baseLen); S.userMoved = false; buildEgo(true, true); }

// ghost preview on hover: the partners an expansion would add, drawn where they would land
let ghostTimer = null;
function showGhosts(i) {
  clearGhosts();
  const nb = visibleNeighbours(i).filter((x) => !S.shown.has(x.j)).slice(0, S.maxExpand);
  if (!nb.length) return;
  finishAnim();
  const spot = placeAround(i, nb.map((x) => x.j)); if (!spot.size) return;
  S.ghosts = nb.map((x) => x.j); S.ghostOwner = i;
  nodesDS.add(nb.map((x) => { const p = spot.get(x.j);
    return { id: 'g' + x.j, x: p.x, y: p.y, label: S.labels === 'none' ? '' : G.sym[x.j], size: 5, shape: 'dot',
      color: { background: COL.ghost, border: COL.ghost, highlight: { background: COL.ghost, border: COL.ink }, hover: { background: COL.ghost, border: COL.ink } },
      font: { size: S.labelSize - 2, color: COL.muted, strokeWidth: 2, strokeColor: COL.canvas } }; }));
  edgesDS.add(nb.map((x) => ({ id: 'ge' + x.e, from: i, to: 'g' + x.j, dashes: [4, 4], width: 1, color: { color: COL.ghost, opacity: 0.8 }, smooth: false })));
}
function clearGhosts() {
  if (!S.ghosts.length) return;
  edgesDS.remove(edgesDS.getIds().filter((id) => String(id).startsWith('ge')));
  nodesDS.remove(S.ghosts.map((j) => 'g' + j));
  S.ghosts = [];
}
function addGhost(j) { // click on a grey node: add only that partner, exactly where it was previewed
  const from = S.ghostOwner, p = net.getPositions(['g' + j])['g' + j];
  clearGhosts();
  if (from == null || !p || S.shown.has(j)) return;
  S.history.push({ node: from, kind: 'add', added: [j] });
  const spot = new Map([[j, { x: p.x, y: p.y }]]);
  addGenes([j], spot, from, spot);
  updateCrumbs();
}

// ---- overview mode
const LAYOUT = new Map(D.layout.idx.map((g, k) => [g, [D.layout.x[k], D.layout.y[k]]]));
function renderOverview() {
  const nodes = new Set(), edges = [];
  for (let e = 0; e < E.n; e++) {
    const a = E.a[e], b = E.b[e];
    if (!(G.rwe[a] || G.rwe[b]) || !LAYOUT.has(a) || !LAYOUT.has(b) || !edgeVisible(e)) continue;
    edges.push(e); nodes.add(a); nodes.add(b);
  }
  nodesDS.add([...nodes].map((i) => { const p = LAYOUT.get(i); return nodeObj(i, p[0], p[1]); }));
  edgesDS.add(edges.map(edgeObj));
  S.shown = nodes;
}

// ---- modules mode: heat map of module over-representation among each RWE gene's top partners
// level 'category' = 16 broad categories (all cells); level 'pathway' = Reactome / CORUM modules called in any dataset (top columns shown)
const MOD = D.modules;
const fdrShade = (fdr, colour) => { // intensity ~ -log10 FDR, saturating at 1e-6; transparent when not significant at 0.25
  if (fdr == null || fdr >= 0.25) return 'transparent';
  const a = Math.min(1, -Math.log10(Math.max(fdr, 1e-12)) / 6);
  return colour + Math.round(40 + 215 * a).toString(16).padStart(2, '0');
};
function renderModules() {
  const L = MOD[S.modLevel], names = L.names;
  const bit = S.ds === 'NextGen' ? 1 : S.ds === 'DepMap2D' ? 2 : 3;
  const isCalled = (m) => (S.ds === 'compare' ? ((m[7] & 3) === 3 ? S.cls.both : (m[7] & 1) ? S.cls.ng : (m[7] & 2) ? S.cls.d2 : false) : (m[7] & bit));
  const called = L.rows.filter(isCalled);
  let cols = [...new Set(called.map((m) => m[1]))];
  if (S.modLevel === 'pathway') { const cnt = {}; called.forEach((m) => { cnt[m[1]] = (cnt[m[1]] || 0) + 1; }); cols = cols.sort((a, b) => cnt[b] - cnt[a]).slice(0, 45); }
  else cols.sort((a, b) => a - b);
  const colSet = new Set(cols), genes = [...new Set(called.filter((m) => colSet.has(m[1])).map((m) => m[0]))];
  const cell = new Map(L.rows.map((m) => [m[0] * 4096 + m[1], m]));
  const firstCol = (i) => { for (let k = 0; k < cols.length; k++) { const m = cell.get(i * 4096 + cols[k]); if (m && isCalled(m)) return k; } return 99; };
  genes.sort((a, b) => firstCol(a) - firstCol(b) || G.sym[a].localeCompare(G.sym[b]));
  const short = (n) => n.replace(/^Reactome: /, '').replace(/ R-HSA-\d+$/, '').replace(/^CORUM: /, 'CORUM · ').slice(0, 58);
  const cw = 18, ch = 15, lw = 96, th = S.modLevel === 'pathway' ? 300 : 190, W = lw + cols.length * cw + 260, H = th + genes.length * ch + 40;
  let svg = `<svg width="${W}" height="${H}" font-size="11" fill="${COL.ink}">`;
  cols.forEach((c, k) => { const x = lw + k * cw + cw / 2;
    svg += `<text transform="translate(${x + 4} ${th - 6}) rotate(-55)" font-size="10.5">${esc(S.modLevel === 'category' ? lab('cat', names[c]) : short(names[c]))}<title>${esc(names[c])} (${L.sizes[c]})</title></text>`; });
  genes.forEach((i, r) => {
    const y = th + r * ch;
    svg += `<text class="rowlab" data-g="${i}" x="${lw - 6}" y="${y + ch - 4}" text-anchor="end" fill="${FAM_PALETTE[Math.max(G.fam[i], 0) % FAM_PALETTE.length]}">${esc(G.sym[i])}${G.lib[i] ? ' †' : ''}</text>`;
    cols.forEach((c, k) => {
      const m = cell.get(i * 4096 + c), x = lw + k * cw;
      svg += `<rect x="${x}" y="${y}" width="${cw - 1}" height="${ch - 1}" fill="${COL.canvas}" stroke="${COL.muted}" stroke-opacity="0.15"/>`;
      if (!m) return;
      const tip = `${G.sym[i]} · ${names[c]}\nNextGen: ${m[5]}/100 top partners, FDR ${fmtP(m[2])}\nDepMap 2D: ${m[6]}/100, FDR ${fmtP(m[3])}\nCas12a 2D (library check): FDR ${fmtP(m[4])}`;
      if (S.ds === 'compare') {
        svg += `<path d="M${x} ${y}h${cw - 1}L${x} ${y + ch - 1}z" fill="${fdrShade(m[2], COL.ng)}"/><path d="M${x + cw - 1} ${y}v${ch - 1}h${-(cw - 1)}z" fill="${fdrShade(m[3], COL.d2)}"/>`;
      } else {
        svg += `<rect x="${x}" y="${y}" width="${cw - 1}" height="${ch - 1}" fill="${fdrShade(S.ds === 'NextGen' ? m[2] : m[3], S.ds === 'NextGen' ? COL.ng : COL.d2)}"/>`;
      }
      if (isCalled(m)) svg += `<rect x="${x + 0.5}" y="${y + 0.5}" width="${cw - 2}" height="${ch - 2}" fill="none" stroke="${COL.ink}" stroke-width="0.9"/>`;
      if (m[7] & 4) svg += `<circle cx="${x + cw - 5}" cy="${y + 4}" r="1.6" fill="${COL.ink}"/>`;
      svg += `<rect x="${x}" y="${y}" width="${cw - 1}" height="${ch - 1}" fill="transparent"><title>${esc(tip)}</title></rect>`;
    });
  });
  const zh = S.lang === 'zh';
  svg += `<text x="${lw}" y="${th + genes.length * ch + 22}" fill="${COL.muted}">${zh ? '颜色深浅 = −log10 FDR（超几何，前 100 个正相关伙伴）；黑框 = 当前数据集 FDR < 0.05；右上角圆点 = Cas12a 2D 也显著（文库检查）；† = 2D 文库敏感基因' : 'shade = −log10 FDR (hypergeometric, top 100 positive partners); outline = FDR < 0.05 in the current dataset; dot = also called in Cas12a 2D (library check); † = 2D library-sensitive gene'}${S.ds === 'compare' ? (zh ? '；左上 = NextGen，右下 = DepMap 2D' : '; upper-left = NextGen, lower-right = DepMap 2D') : ''}${S.modLevel === 'pathway' ? (zh ? `；显示被最多 RWE 基因富集的 ${cols.length} 个模块` : `; showing the ${cols.length} modules called for most RWE genes`) : ''}</text></svg>`;
  $('matrix').innerHTML = svg;
  $('matrix').querySelectorAll('.rowlab').forEach((el) => el.addEventListener('click', () => selectGene(+el.dataset.g)));
  S.shown = new Set(genes);
  S.matrixCells = called.filter((m) => colSet.has(m[1])).length;
}

// ------------------------------------------------------------------------------------------------ labels, counts, legend, crumbs
function refreshLabels() { nodesDS.update(nodesDS.getIds().filter((id) => typeof id === 'number').map((i) => ({ id: i, label: labelFor(i) }))); }
function updateCounts() {
  $('counts').textContent = S.mode === 'modules' ? t('counts_modules', S.shown.size, S.matrixCells || 0) : t('counts', nodesDS.length - S.ghosts.length, edgesDS.length - S.ghosts.length);
}
function updateLegend() {
  const ln = (c, k, dash) => `<div><span class="ln${dash ? ' dash' : ''}" style="border-color:${c}"></span>${t(k)}</div>`;
  let h = '';
  if (S.ds === 'NextGen') h += ln(COL.ng, 'edge_ng');
  if (S.ds === 'DepMap2D') h += ln(COL.d2, 'edge_2d');
  if (S.ds === 'compare') h += ln(COL.both, 'c_both') + ln(COL.ng, 'c_ng') + ln(COL.d2, 'c_2d');
  h += ln(COL.muted, 'edge_neg', true);
  h += `<div class="sec" style="margin-top:8px">${t('legend_nodes')}</div><div>${t('node_rwe')}</div><div>${t('node_partner')}</div>`;
  if (S.mode === 'explore') h += `<div><span class="sw" style="background:${COL.ghost}"></span>${t('node_ghost')}</div><div class="hint">${t('node_root')}</div><div class="hint">${t('shift_hint')}</div>`;
  $('legend').innerHTML = h;
  $('clsRow').style.display = S.ds === 'compare' ? '' : 'none';
}
function updateCrumbs() {
  const nodes = S.mode === 'explore' && S.root != null && S.history.length ? [...new Set([S.root, ...S.history.map((h) => h.node)])] : [];
  $('crumbs').innerHTML = nodes.length ? `<span class="hint">${t('history')}:</span>` + nodes.map((i) => `<span data-g="${i}">${esc(G.sym[i])}</span>`).join('') : '';
  $('crumbs').querySelectorAll('span[data-g]').forEach((el) => el.addEventListener('click', () => focusGene(+el.dataset.g)));
  $('undoBtn').disabled = !(S.mode === 'explore' && S.history.length > S.baseLen);
}
function updateEgoBar() { // one line above the canvas: whose network this is and how many of its partners are on screen
  const bar = $('egoBar'), on = S.mode === 'explore' && S.root != null;
  bar.style.display = on ? 'block' : 'none';
  if (!on) return;
  const g = `<b>${esc(G.sym[S.root])}</b>`, ds = S.ds === 'compare' ? t('ds_cmp') : dsName(S.ds), nb = visibleNeighbours(S.root), shown = nb.filter((x) => S.shown.has(x.j)).length;
  bar.innerHTML = nb.length ? t('ego_bar', g, ds, shown, nb.length) + (shown < nb.length ? ` · <span class="k">${t(shown ? 'ego_more' : 'ego_start')}</span>` : '')
    : t('ego_none', g, ds) + (G.infam[S.root] ? '' : `<br><span class="k">${t('not_dependency')}</span>`);
}
function restyle() { // colours, fonts and labels changed; positions stay
  if (S.mode === 'modules') { render(); return; }
  nodesDS.update(nodesDS.getIds().filter((id) => typeof id === 'number').map((i) => nodeObj(i)));
  edgesDS.update(edgesDS.getIds().filter((id) => String(id)[0] === 'e').map((id) => edgeObj(+String(id).slice(1))));
  updateLegend();
}
const syncSeg = (id, v) => document.querySelectorAll(`#${id} button`).forEach((b) => b.classList.toggle('on', b.dataset.v === v));
function writeHash() { // shareable link: the centre gene in explore, otherwise the selected gene
  const g = S.mode === 'explore' ? S.root : S.selected;
  try { history.replaceState(null, '', '#' + (g != null ? 'gene=' + encodeURIComponent(G.sym[g]) + '&' : '') + 'ds=' + S.ds); } catch (e) { /* file:// */ }
}

// ------------------------------------------------------------------------------------------------ hover tooltip
function showTip(html, pointer) {
  const tip = $('tip'); tip.innerHTML = html; tip.style.display = 'block';
  const rect = $('net').getBoundingClientRect(), x = rect.left + (pointer?.DOM?.x ?? 0), y = rect.top + (pointer?.DOM?.y ?? 0), w = tip.offsetWidth, h = tip.offsetHeight;
  tip.style.left = Math.min(window.innerWidth - w - 8, x + 14) + 'px'; tip.style.top = Math.min(window.innerHeight - h - 8, y + 14) + 'px';
}
function hideTip() { $('tip').style.display = 'none'; }
function geneTip(i, ghost) {
  const cls = (c) => (c >= 0 ? lab('cls', D.dict.cls[c]) : '—');
  const pct = (v) => (v == null ? '—' : Math.round(v * 100) + '%');
  return `<b>${esc(G.sym[i])}</b>${G.rwe[i] ? ' <span class="pill both">RWE</span>' : ''}${G.lib[i] ? ' <span class="pill warn">2D lib</span>' : ''}<br>
    <span class="k">${G.rwe[i] ? t('family') + ': ' + esc(D.dict.fam[G.fam[i]] || '—') : t('category') + ': ' + esc(lab('cat', D.dict.cat[G.cat[i]] || ''))}</span><br>
    <table class="kv"><tr><th></th><th class="ng">NextGen</th><th class="d2">DepMap 2D</th></tr>
    <tr><td>${t('class')}</td><td>${cls(G.cng[i])}</td><td>${cls(G.c2d[i])}</td></tr>
    <tr><td>${t('mean')}</td><td>${fmt(G.mng[i])}</td><td>${fmt(G.m2d[i])}</td></tr>
    <tr><td>${t('frac_dep')}</td><td>${pct(G.fng[i])}</td><td>${pct(G.f2d[i])}</td></tr>
    <tr><td>${t('n_partners')}</td><td>${G.infam[i] ? G.dng[i] : '—'}</td><td>${G.infam[i] ? G.d2d[i] : '—'}</td></tr></table>
    ${ghost ? '<span class="k">' + (S.lang === 'zh' ? '点击加入网络' : 'click to add') + '</span>'
    : S.mode !== 'explore' || !S.shown.has(i) ? '' : isExpanded(i) ? `<span class="k">${t('tip_collapse')}</span>`
    : visibleNeighbours(i).some((x) => !S.shown.has(x.j)) ? `<span class="k">${t('tip_expand')}</span>` : ''}`;
}
function edgeTip(e) {
  const call = (bit) => (E.f[e] & bit ? `<span class="pill ${bit === 1 ? 'ng' : 'd2'}">${t('called')}</span>` : t('not_called'));
  return `<b>${esc(G.sym[E.a[e]])} — ${esc(G.sym[E.b[e]])}</b><table class="kv"><tr><th></th><th>r</th><th></th></tr>
    <tr><td class="ng">NextGen</td><td>${fmt(E.rng[e] / 100)}</td><td>${call(1)}</td></tr>
    <tr><td class="d2">DepMap 2D</td><td>${fmt(E.r2d[e] / 100)}</td><td>${call(2)}</td></tr></table>${LP_EDGES.has(e) ? `<span class="pill warn">${t('lineage_pair')}</span>` : ''}`;
}
function onHoverNode(id, ev) {
  if (typeof id === 'string' && id.startsWith('g')) { clearTimeout(ghostTimer); showTip(geneTip(+id.slice(1), true), ev); return; }
  if (typeof id === 'string') return;
  S.hoverNode = id; showTip(geneTip(id, false), ev);
  if (S.mode === 'explore') { clearTimeout(ghostTimer); ghostTimer = setTimeout(() => (isExpanded(id) ? clearGhosts() : showGhosts(id)), 120); }
}
function onBlurNode() {
  hideTip(); S.hoverNode = null;
  if (S.mode === 'explore') { clearTimeout(ghostTimer); ghostTimer = setTimeout(clearGhosts, 700); }
}
function onHoverEdge(id, pointer) {
  const s = String(id), e = s.startsWith('ge') ? +s.slice(2) : s.startsWith('e') ? +s.slice(1) : null;
  if (e != null) showTip(edgeTip(e), pointer);
}
function onClick(p) {
  hideTip();
  let id = p.nodes[0];
  if (S.mode === 'explore' && typeof id === 'number' && !S.shown.has(id)) { finishAnim(); id = net.getNodeAt(p.pointer.DOM); } // clicked a gene that is folding away: use what is underneath
  if (id != null) {
    if (typeof id === 'string' && id.startsWith('g')) { addGhost(+id.slice(1)); return; }
    if (typeof id === 'string') return;
    if (S.mode === 'explore') {
      clearGhosts(); selectGene(id); if (p.event?.srcEvent?.shiftKey || isExpanded(id)) collapse(id); else expand(id);
      if (S.shown.has(id)) showTip(geneTip(id, false), p.pointer); // the pointer is still on the gene: say what the next click does
      return;
    }
    selectGene(id); if (S.mode === 'overview') centreOn(id, 450);
    return;
  }
  const eid = p.edges[0];
  if (eid != null && String(eid).startsWith('e')) { showEdge(+String(eid).slice(1)); return; }
  closePanel();
}

// ------------------------------------------------------------------------------------------------ gene card
function selectGene(i) {
  S.selected = i; refreshLabels();
  if (S.mode !== 'explore') writeHash();
  showGene(i);
}
function focusGene(i) { // breadcrumb: look at a gene that is already on screen (the centre gene re-centres the whole network)
  if (S.mode !== 'explore' || !S.shown.has(i)) { showEgo(i); return; }
  selectGene(i);
  if (i === S.root) { S.userMoved = false; recentre(450); } else { S.userMoved = true; centreOn(i, 450); }
}
function partnerRows(i, ds) {
  const rows = [];
  for (let k = ADJ.off[i]; k < ADJ.off[i + 1]; k++) { const e = ADJ.idx[k]; if (E.f[e] & (ds === 'NextGen' ? 1 : 2)) rows.push({ e, j: otherEnd(e, i), r: edgeR(e, ds) }); }
  return rows.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
}
function showGene(i) {
  const c = D.cards[i], cls = (k) => (k >= 0 ? lab('cls', D.dict.cls[k]) : '—'), pct = (v) => (v == null ? '—' : Math.round(v * 100) + '%');
  let h = `<button class="close" onclick="closePanel()">✕</button><h2>${esc(G.sym[i])}</h2><div class="sub">${esc(G.key[i])}${G.rwe[i] ? ' · RWE' : ''}</div>`;
  if (G.rwe[i]) h += `<div class="sub">${t('family')}: ${esc(D.dict.fam[G.fam[i]] || '—')} · ${t('role')}: ${esc(lab('role', D.dict.role[G.role[i]] || '—'))}${c && c.mods ? ' · ' + t('mods') + ': ' + esc(c.mods) : ''}</div>`;
  else h += `<div class="sub">${t('category')}: ${esc(lab('cat', D.dict.cat[G.cat[i]] || ''))}</div>`;
  const isRoot = S.mode === 'explore' && S.root === i, onScreen = S.mode === 'explore' && S.shown.has(i);
  h += `<div class="btnrow">${isRoot ? '' : `<button class="iconbtn" id="egoBtn">${t('centre')}</button>`}${onScreen ? `<button class="iconbtn" id="expBtn">${t('expand')}</button>` : ''}</div>`;
  if (c && c.caveat) h += `<div class="note">${t('caveat')}: ${esc(c.caveat)}</div>`;
  h += `<table class="kv"><tr><th></th><th class="ng">NextGen (n=${D.meta.n.NextGen})</th><th class="d2">DepMap 2D (n=${D.meta.n.DepMap2D})</th></tr>
    <tr><td>${t('class')}</td><td>${cls(G.cng[i])}</td><td>${cls(G.c2d[i])}</td></tr><tr><td>${t('mean')}</td><td>${fmt(G.mng[i])}</td><td>${fmt(G.m2d[i])}</td></tr>
    <tr><td>${t('frac_dep')}</td><td>${pct(G.fng[i])}</td><td>${pct(G.f2d[i])}</td></tr><tr><td>${t('n_partners')}</td><td>${G.infam[i] ? G.dng[i] : '—'}</td><td>${G.infam[i] ? G.d2d[i] : '—'}</td></tr></table>`;
  if (!G.infam[i]) h += `<div class="note">${t('not_dependency')}</div>`;
  if (G.lib[i]) h += `<div class="note">${t('lib_flag')}</div>`;
  if (c) {
    const extra = [];
    if (c.org) extra.push(`<tr><td>${t('organoid')}</td><td colspan="2">Δ=${fmt(c.org[0])}, FDR=${fmtP(c.org[1])}${c.org[2] ? ' · ' + t('essential_dropped') : ''}${c.org[3] ? ' · <b>' + (c.org[0] < 0 ? t('more_dep_org') : t('less_dep_org')) + '</b>' : ''}</td></tr>`);
    if (c.lin_ng || c.lin_2d) extra.push(`<tr><td>${t('top_lineage')}</td><td class="ng">${c.lin_ng ? esc(lab('lineage', c.lin_ng[0])) + ' ' + fmt(c.lin_ng[1]) + ' (FDR ' + fmtP(c.lin_ng[2]) + ')' : '—'}</td><td class="d2">${c.lin_2d ? esc(lab('lineage', c.lin_2d[0])) + ' ' + fmt(c.lin_2d[1]) + ' (FDR ' + fmtP(c.lin_2d[2]) + ')' : '—'}</td></tr>`);
    if (c.bio) extra.push(`<tr><td>${t('biomarkers')}</td><td colspan="2">${esc(c.bio)}</td></tr>`);
    if (c.libdiff != null) extra.push(`<tr><td>${t('library')}</td><td colspan="2">${fmt(c.libdiff)}</td></tr>`);
    if (extra.length) h += `<h3>${t('gene_card')}</h3><table class="kv">${extra.join('')}</table>`;
    const mods = (ds) => (c.mod && c.mod[ds] && c.mod[ds].length ? c.mod[ds].map(([name, k, fdr]) => `<span class="pill ${ds === 'NextGen' ? 'ng' : 'd2'}" title="${esc(name)}">${esc(name.replace(/^Reactome: /, '').replace(/ R-HSA-\d+$/, ''))} (${k}/100, FDR ${fmtP(fdr)})</span>`).join('') : '—');
    if (c.mod && (c.mod.NextGen || c.mod.DepMap2D)) h += `<h3>${t('modules')}</h3><div>NextGen: ${mods('NextGen')}</div><div>DepMap 2D: ${mods('DepMap2D')}</div>`;
    if (c.claim) h += `<h3>${t('literature')}</h3><div class="note">${esc(c.claim)}…<br><b>${esc(c.conc || '')}</b></div>`;
  }
  h += `<h3>${t('strip')}</h3><div class="plot" id="strip">…</div>`;
  if (G.infam[i]) {
    h += `<h3>${t('partners')} <span class="hint">(${t('select_edge')})</span></h3><div class="tabs"><button data-ds="NextGen" class="on">NextGen (${G.dng[i]})</button><button data-ds="DepMap2D">DepMap 2D (${G.d2d[i]})</button></div><div id="plist"></div>`;
  }
  openPanel(h);
  const list = (ds) => { const rows = partnerRows(i, ds).slice(0, 200);
    $('plist').innerHTML = `<table class="list"><tr><th>${S.lang === 'zh' ? '伙伴' : 'partner'}</th><th class="ng">r NextGen</th><th class="d2">r 2D</th><th>${t('cls_label')}</th></tr>` +
      rows.map((x) => `<tr class="pick" data-e="${x.e}"><td>${esc(G.sym[x.j])}${G.rwe[x.j] ? ' ·RWE' : ''}</td><td>${fmt(E.rng[x.e] / 100)}</td><td>${fmt(E.r2d[x.e] / 100)}</td><td>${classPill(x.e)}</td></tr>`).join('') + '</table>';
    $('plist').querySelectorAll('tr[data-e]').forEach((tr) => tr.addEventListener('click', () => showEdge(+tr.dataset.e))); };
  if (G.infam[i]) { list('NextGen'); $('panel').querySelectorAll('.tabs button').forEach((b) => b.addEventListener('click', () => { $('panel').querySelectorAll('.tabs button').forEach((x) => x.classList.toggle('on', x === b)); list(b.dataset.ds); })); }
  $('egoBtn')?.addEventListener('click', () => showEgo(i));
  $('expBtn')?.addEventListener('click', () => { $('expBtn').disabled = true; expand(i); });
  drawStrip(i);
}
function classPill(e) { const c = edgeClass(e); return c ? `<span class="pill ${c}">${t('c_' + c)}</span>` : ''; }

async function drawStrip(i) {
  if (typeof DecompressionStream === 'undefined') { $('strip').textContent = t('no_gz'); return; }
  const V = { NextGen: await vectors('NextGen'), DepMap2D: await vectors('DepMap2D') };
  if (!$('strip')) return;
  const L = D.dict.lineage, rows = [];
  L.forEach((name, li) => DATASETS.forEach((ds) => {
    const v = geneVector(V[ds], ds, i), lin = D.screens[ds].lin, pts = [];
    for (let s = 0; s < v.length; s++) if (lin[s] === li && !Number.isNaN(v[s])) pts.push(v[s]);
    if (pts.length) rows.push({ li, ds, pts });
  }));
  const all = rows.flatMap((r) => r.pts); if (!all.length) { $('strip').textContent = '—'; return; }
  let lo = Math.min(...all), hi = Math.max(...all); const pad = (hi - lo) * 0.05 || 0.1; lo -= pad; hi += pad;
  const W = 390, lw = 112, rh = 9, H = rows.length * rh + 26, X = (v) => lw + (v - lo) / (hi - lo) * (W - lw - 8);
  let s = `<svg viewBox="0 0 ${W} ${H}" font-size="10" fill="${COL.ink}">`;
  [Math.ceil(lo), 0, Math.floor(hi)].filter((v, k, a) => v >= lo && v <= hi && a.indexOf(v) === k).forEach((v) => { s += `<line x1="${X(v)}" x2="${X(v)}" y1="2" y2="${H - 18}" stroke="${COL.muted}" stroke-width="0.5" stroke-dasharray="${v === 0 ? '' : '2 2'}"/><text x="${X(v)}" y="${H - 6}" text-anchor="middle" fill="${COL.muted}">${v}</text>`; });
  let prev = -1;
  rows.forEach((r, k) => {
    const y = 6 + k * rh, col = r.ds === 'NextGen' ? COL.ng : COL.d2;
    if (r.li !== prev) { s += `<text x="${lw - 6}" y="${y + 6}" text-anchor="end" fill="${D.dict.lcol[r.li]}">${esc(lab('lineage', L[r.li]))}</text>`; prev = r.li; }
    r.pts.forEach((v) => { s += `<circle cx="${X(v).toFixed(1)}" cy="${(y + 2 + Math.random() * 5).toFixed(1)}" r="1.6" fill="${col}" fill-opacity="0.55"/>`; });
    const m = r.pts.reduce((a, b) => a + b, 0) / r.pts.length;
    s += `<line x1="${X(m)}" x2="${X(m)}" y1="${y}" y2="${y + 9}" stroke="${col}" stroke-width="2"/>`;
  });
  s += `<text x="${lw}" y="${H - 6}" fill="${COL.ng}">● NextGen</text><text x="${lw + 70}" y="${H - 6}" fill="${COL.d2}">● DepMap 2D</text></svg>`;
  $('strip').innerHTML = s;
}

// ------------------------------------------------------------------------------------------------ edge card (lineage decomposition)
function corrStats(x, y, idx) {
  let n = 0, sx = 0, sy = 0; const xs = [], ys = [];
  for (const s of idx) if (!Number.isNaN(x[s]) && !Number.isNaN(y[s])) { xs.push(x[s]); ys.push(y[s]); sx += x[s]; sy += y[s]; n++; }
  if (n < 3) return { n, r: NaN, xs, ys, mx: sx / (n || 1), my: sy / (n || 1) };
  const mx = sx / n, my = sy / n; let sxx = 0, syy = 0, sxy = 0;
  for (let k = 0; k < n; k++) { const a = xs[k] - mx, b = ys[k] - my; sxx += a * a; syy += b * b; sxy += a * b; }
  return { n, r: sxy / Math.sqrt(sxx * syy), xs, ys, mx, my, sxx, syy, sxy };
}
function decompose(x, y, ds) {
  const lin = D.screens[ds].lin, all = [...Array(lin.length).keys()];
  const per = D.dict.lineage.map((_, li) => corrStats(x, y, all.filter((s) => lin[s] === li)));
  const tot = corrStats(x, y, all);
  let sxx = 0, syy = 0, sxy = 0; const means = [];
  per.forEach((p) => { if (p.n >= 3) { sxx += p.sxx; syy += p.syy; sxy += p.sxy; means.push([p.mx, p.my]); } });
  const within = sxy / Math.sqrt(sxx * syy);
  let between = NaN; if (means.length >= 4) { const mx = means.reduce((a, b) => a + b[0], 0) / means.length, my = means.reduce((a, b) => a + b[1], 0) / means.length;
    let a = 0, b = 0, c = 0; means.forEach(([u, v]) => { a += (u - mx) ** 2; b += (v - my) ** 2; c += (u - mx) * (v - my); }); between = c / Math.sqrt(a * b); }
  const use = per.filter((p) => p.n >= 6 && Number.isFinite(p.r));
  let pQ = NaN; if (use.length >= 2) { const w = use.map((p) => p.n - 3), z = use.map((p) => atanh(clampR(p.r))), sw = w.reduce((a, b) => a + b, 0);
    const zb = w.reduce((a, b, k) => a + b * z[k], 0) / sw; pQ = chi2sf(w.reduce((a, b, k) => a + b * (z[k] - zb) ** 2, 0), use.length - 1); }
  return { per, tot, within, between, pQ };
}
async function showEdge(e) {
  const a = E.a[e], b = E.b[e], lp = LP.find((x) => x.e === e);
  let h = `<button class="close" onclick="closePanel()">✕</button><h2>${esc(G.sym[a])} — ${esc(G.sym[b])}</h2><div class="sub">${t('edge_card')} · ${classPill(e) || t('not_called')}</div><div id="edgeBody">…</div>`;
  openPanel(h);
  if (typeof DecompressionStream === 'undefined') { $('edgeBody').textContent = t('no_gz'); return; }
  const V = { NextGen: await vectors('NextGen'), DepMap2D: await vectors('DepMap2D'), Cas12a2D: await vectors('Cas12a2D') };
  const dec = {}, vec = {};
  const c12 = corrStats(geneVector(V.Cas12a2D, 'Cas12a2D', a), geneVector(V.Cas12a2D, 'Cas12a2D', b), [...Array(D.screens.Cas12a2D.id.length).keys()]);
  DATASETS.forEach((ds) => { vec[ds] = [geneVector(V[ds], ds, a), geneVector(V[ds], ds, b)]; dec[ds] = decompose(vec[ds][0], vec[ds][1], ds); });
  const rr = { NextGen: E.rng[e] / 100, DepMap2D: E.r2d[e] / 100 };
  const nn = { NextGen: dec.NextGen.tot.n, DepMap2D: dec.DepMap2D.tot.n };
  const zd = (atanh(clampR(rr.NextGen)) - atanh(clampR(rr.DepMap2D))) / Math.sqrt(1 / (nn.NextGen - 3) + 1 / (nn.DepMap2D - 3));
  const pdiff = 2 * (1 - normCdf(Math.abs(zd)));
  let body = `<table class="kv"><tr><th></th><th>r</th><th>n</th><th>p</th><th></th></tr>` + DATASETS.map((ds) => `<tr><td class="${ds === 'NextGen' ? 'ng' : 'd2'}">${dsName(ds)}</td><td>${fmt(rr[ds])}</td><td>${nn[ds]}</td><td>${fmtP(pFromR(rr[ds], nn[ds]))}</td><td>${E.f[e] & (ds === 'NextGen' ? 1 : 2) ? '<span class="pill ' + (ds === 'NextGen' ? 'ng' : 'd2') + '">' + t('called') + '</span>' : t('not_called')}</td></tr>`).join('') +
    `<tr><td class="k" title="${esc(t('c12_note'))}">Cas12a 2D ⓘ</td><td>${fmt(c12.r)}</td><td>${c12.n}</td><td>${fmtP(pFromR(c12.r, c12.n))}</td><td class="k">${t('library_check')}</td></tr>` +
    `<tr><td colspan="5">${t('diff_p')}: ${fmtP(pdiff)}</td></tr></table>`;
  if (E.f[e] & 4) body += `<div class="note">${t('lib_flag')}</div>`;
  if (lp) body += `<div class="note"><b>${t('lineage_pair')}</b>: ${esc(lab('lineage', lp.lin))} r=${fmt(lp.r)}, FDR=${fmtP(lp.q)}; ${S.lang === 'zh' ? '其余谱系' : 'other lineages'} r=${fmt(lp.rr)} (${esc(lp.cls)})</div>`;
  body += `<h3>${t('decomposition')}</h3><table class="kv"><tr><th></th><th class="ng">NextGen</th><th class="d2">DepMap 2D</th></tr>
    <tr><td>${t('pooled')}</td><td>${fmt(dec.NextGen.tot.r)}</td><td>${fmt(dec.DepMap2D.tot.r)}</td></tr>
    <tr><td>${t('within')}</td><td>${fmt(dec.NextGen.within)}</td><td>${fmt(dec.DepMap2D.within)}</td></tr>
    <tr><td>${t('between')}</td><td>${fmt(dec.NextGen.between)}</td><td>${fmt(dec.DepMap2D.between)}</td></tr>
    <tr><td>${t('het')}</td><td>${fmtP(dec.NextGen.pQ)}</td><td>${fmtP(dec.DepMap2D.pQ)}</td></tr></table>
    <div class="hint">${t('forest_note')}</div><div class="plot">${forestSVG(dec, lp)}</div>
    <h3>${t('scatter')}</h3><label class="row"><input type="checkbox" id="demean" ${S.demean ? 'checked' : ''}> ${t('demean')}</label>
    <div class="plot" id="scatters"></div><div class="lg" id="linLegend"></div>`;
  $('edgeBody').innerHTML = body;
  const redraw = () => { $('scatters').innerHTML = DATASETS.map((ds) => scatterSVG(vec[ds][0], vec[ds][1], ds, G.sym[a], G.sym[b])).join('');
    $('linLegend').innerHTML = D.dict.lineage.map((L, li) => `<span data-li="${li}" class="${S.hiddenLin.has(li) ? 'off' : ''}"><span class="sw" style="background:${D.dict.lcol[li]}"></span>${esc(lab('lineage', L))}</span>`).join('') + ' <span>● 3D ▲ 2D</span>';
    $('linLegend').querySelectorAll('span[data-li]').forEach((el) => el.addEventListener('click', () => { const li = +el.dataset.li; S.hiddenLin.has(li) ? S.hiddenLin.delete(li) : S.hiddenLin.add(li); redraw(); })); };
  $('demean').addEventListener('change', () => { S.demean = $('demean').checked; redraw(); });
  redraw();
}
function normCdf(z) { const t_ = 1 / (1 + 0.2316419 * z), d = 0.3989423 * Math.exp(-z * z / 2); const p = d * t_ * (0.3193815 + t_ * (-0.3565638 + t_ * (1.781478 + t_ * (-1.821256 + t_ * 1.330274)))); return 1 - p; }
function forestSVG(dec, lp) {
  const L = D.dict.lineage, rows = L.map((_, li) => li).filter((li) => DATASETS.some((ds) => dec[ds].per[li].n >= 6));
  const lw = 104, pw = 190, tw = 92, W = lw + pw + tw, rh = 20, top = 16, H = top + rows.length * rh + 22, X = (v) => lw + (Math.max(-1, Math.min(1, v)) + 1) / 2 * pw;
  const ci = (r, n) => { const z = atanh(clampR(r)), se = 1 / Math.sqrt(Math.max(n - 3, 1)); return [Math.tanh(z - 1.96 * se), Math.tanh(z + 1.96 * se)]; };
  let s = `<svg viewBox="0 0 ${W} ${H}" font-size="10" fill="${COL.ink}"><text x="${lw + pw + 6}" y="10" fill="${COL.ng}">NG</text><text x="${lw + pw + 48}" y="10" fill="${COL.d2}">2D</text>`;
  s += `<line x1="${X(0)}" x2="${X(0)}" y1="${top - 4}" y2="${H - 18}" stroke="${COL.muted}" stroke-width="0.6"/>`;
  DATASETS.forEach((ds) => { const r = dec[ds].tot.r; if (Number.isFinite(r)) s += `<line x1="${X(r)}" x2="${X(r)}" y1="${top - 4}" y2="${H - 18}" stroke="${ds === 'NextGen' ? COL.ng : COL.d2}" stroke-dasharray="3 3"/>`; });
  rows.forEach((li, k) => {
    const y = top + k * rh + rh / 2;
    if (lp && lp.lin === L[li]) s += `<rect x="0" y="${y - rh / 2}" width="${W}" height="${rh}" fill="#D55E00" opacity="0.12"/>`;
    s += `<text x="${lw - 6}" y="${y + 3}" text-anchor="end" fill="${D.dict.lcol[li]}">${esc(lab('lineage', L[li]))}</text>`;
    DATASETS.forEach((ds, d) => {
      const p = dec[ds].per[li], col = ds === 'NextGen' ? COL.ng : COL.d2, yy = y + (d ? 4 : -4);
      if (p.n >= 6 && Number.isFinite(p.r)) { const [lo, hi] = ci(p.r, p.n);
        s += `<line x1="${X(lo)}" x2="${X(hi)}" y1="${yy}" y2="${yy}" stroke="${col}" stroke-width="1.4"/><circle cx="${X(p.r)}" cy="${yy}" r="3.2" fill="${d ? COL.canvas : col}" stroke="${col}" stroke-width="1.4"><title>${dsName(ds)} · ${L[li]} · n=${p.n} · r=${fmt(p.r)}</title></circle>`; }
      s += `<text x="${lw + pw + 6 + d * 42}" y="${y + 3}" fill="${col}">${p.n >= 6 ? fmt(p.r) : '—'}</text>`;
    });
  });
  [-1, -0.5, 0, 0.5, 1].forEach((v) => { s += `<text x="${X(v)}" y="${H - 6}" text-anchor="middle" fill="${COL.muted}">${v}</text>`; });
  return s + '</svg>';
}
function scatterSVG(x, y, ds, xl, yl) {
  const lin = D.screens[ds].lin, type = D.screens[ds].type, pts = [];
  const mean = {}; if (S.demean) D.dict.lineage.forEach((_, li) => { const st = corrStats(x, y, [...Array(lin.length).keys()].filter((s) => lin[s] === li)); mean[li] = [st.mx, st.my]; });
  for (let s = 0; s < x.length; s++) { if (Number.isNaN(x[s]) || Number.isNaN(y[s]) || S.hiddenLin.has(lin[s])) continue; const m = S.demean ? mean[lin[s]] : [0, 0]; pts.push([x[s] - m[0], y[s] - m[1], lin[s], type[s], s]); }
  const W = 200, H = 200, ml = 30, mb = 26, mt = 16, mr = 6;
  if (!pts.length) return `<svg viewBox="0 0 ${W} ${H}" style="width:49%;display:inline-block"></svg>`;
  let x0 = Math.min(...pts.map((p) => p[0])), x1 = Math.max(...pts.map((p) => p[0])), y0 = Math.min(...pts.map((p) => p[1])), y1 = Math.max(...pts.map((p) => p[1]));
  const px = (x1 - x0) * 0.06 || 0.1, py = (y1 - y0) * 0.06 || 0.1; x0 -= px; x1 += px; y0 -= py; y1 += py;
  const X = (v) => ml + (v - x0) / (x1 - x0) * (W - ml - mr), Y = (v) => H - mb - (v - y0) / (y1 - y0) * (H - mt - mb);
  const col = ds === 'NextGen' ? COL.ng : COL.d2;
  let s = `<svg viewBox="0 0 ${W} ${H}" style="width:49%;display:inline-block" font-size="9" fill="${COL.ink}"><text x="${ml}" y="10" fill="${col}" font-weight="600">${dsName(ds)} (n=${pts.length})</text>`;
  s += `<rect x="${ml}" y="${mt}" width="${W - ml - mr}" height="${H - mt - mb}" fill="none" stroke="${COL.muted}" stroke-width="0.5"/>`;
  // least-squares line through the shown points
  const n = pts.length, mx = pts.reduce((a, p) => a + p[0], 0) / n, my = pts.reduce((a, p) => a + p[1], 0) / n; let sxx = 0, sxy = 0, syy = 0;
  pts.forEach((p) => { sxx += (p[0] - mx) ** 2; sxy += (p[0] - mx) * (p[1] - my); syy += (p[1] - my) ** 2; });
  if (sxx > 0) { const b = sxy / sxx; s += `<line x1="${X(x0)}" x2="${X(x1)}" y1="${Y(my + b * (x0 - mx))}" y2="${Y(my + b * (x1 - mx))}" stroke="${col}" stroke-width="1.2" stroke-dasharray="4 3"/>`; }
  const scr = D.screens[ds];
  pts.forEach((p) => { const c = D.dict.lcol[p[2]] || '#888', cx = X(p[0]).toFixed(1), cy = Y(p[1]).toFixed(1),
    tip = `<title>${scr.id[p[4]]} · ${scr.model[p[4]]} · ${scr.sub[p[4]]} · ${scr.type[p[4]]} · ${scr.lib[p[4]]}\n${xl}=${p[0].toFixed(2)}  ${yl}=${p[1].toFixed(2)}</title>`;
    s += String(p[3])[0] === '3' ? `<circle cx="${cx}" cy="${cy}" r="2.3" fill="${c}" fill-opacity="0.75">${tip}</circle>` : `<path d="M${cx} ${(+cy - 2.8).toFixed(1)}l2.5 4.6h-5z" fill="${c}" fill-opacity="0.75">${tip}</path>`; });
  s += `<text x="${(ml + W - mr) / 2}" y="${H - 4}" text-anchor="middle">${esc(xl)}</text><text transform="translate(9 ${(mt + H - mb) / 2}) rotate(-90)" text-anchor="middle">${esc(yl)}</text>`;
  s += `<text x="${W - mr}" y="${mt + 10}" text-anchor="end" fill="${col}">r=${fmt(sxy / Math.sqrt(sxx * syy))}</text>`;
  return s + '</svg>';
}

// ------------------------------------------------------------------------------------------------ panel, about
function openPanel(html) { const p = $('panel'), was = p.style.display === 'block'; p.innerHTML = html; p.style.display = 'block'; p.scrollTop = 0; if (!was) panelMoved(); }
function closePanel() { const was = $('panel').style.display === 'block'; $('panel').style.display = 'none'; S.selected = null; refreshLabels(); if (was) panelMoved(); }
function panelMoved() { if (S.mode === 'explore' && !S.userMoved && !anim) recentre(350); } // keep the centre gene in the middle of the uncovered canvas
window.closePanel = closePanel;
function aboutHTML() {
  const zh = S.lang === 'zh';
  const cite = `Zhao S. VisRWE: dependencies and co-dependencies of RNA-modification writers, readers and erasers in DepMap NextGen organoid and DepMap 2D CRISPR screens. GitHub repository https://github.com/shuyang-zhao/VisRWE_DepMap (2026).`;
  const paper = 'Neiswender, Maffa, Brenan et al. A dependency map enhanced with next-generation 3D cancer models. Nature (2026). doi:10.1038/s41586-026-10843-7';
  return `<button class="close" onclick="document.getElementById('about').style.display='none'" style="float:right;border:0;background:none;font-size:18px;cursor:pointer">✕</button>
  <h2>VisRWE · DepMap</h2>
  <p><b>${zh ? '作者' : 'Author'}</b>: Shuyang Zhao, Department of Chemistry, University of California, Riverside.</p>
  <p>${zh ? '本页面用 DepMap NextGen 论文公开的数据与统计口径（筛选级 Chronos gene effect；依赖基因之间的 Pearson 相关；每个结果表一次 Benjamini–Hochberg 校正；显著 = FDR < 0.005 且 |r| > 0.3），对 217 个 RNA 修饰 writer / reader / eraser（RWE）基因及其共依赖伙伴，并排展示 NextGen（147 个类器官 / 神经球筛选，Cas12a 文库）与 DepMap 2D（同谱系的 597 个传统 2D 筛选）的结果。'
    : 'This page shows, side by side, the dependencies and co-dependencies of 217 RNA-modification writer / reader / eraser (RWE) genes in DepMap NextGen (147 organoid / neurosphere screens, Cas12a library) and DepMap 2D (597 traditional 2D screens of the same lineages), using the data and statistical conventions of the NextGen paper (screen-level Chronos gene effect; Pearson correlation among dependency genes; one Benjamini–Hochberg correction per results table; call = FDR < 0.005 and |r| > 0.3).'}</p>
  <h3>${zh ? '如何引用' : 'How to cite'}</h3><pre>${esc(cite)}</pre>
  <p>${zh ? '同时请引用原始数据：' : 'Please also cite the original data:'}</p><pre>${esc(paper)}
Tsherniak A, et al. Defining a cancer dependency map. Cell 170, 564–576 (2017).
Dempster JM, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biology 22, 343 (2021).
DepMap portal, Broad Institute: https://depmap.org (file set "NextGen Model Manuscript 2026").</pre>
  <h3>${zh ? '说明' : 'Notes'}</h3><ul>
  <li>${zh ? 'NextGen 全部使用 Cas12a（Humagne-CD）文库；2D 对照以 Cas9 为主（Avana / KY），少数基因的 2D 效应受文库影响，已标记。' : 'All NextGen screens use the Cas12a (Humagne-CD) library; the 2D comparator is mostly Cas9 (Avana / KY); genes whose 2D effect is library-driven are flagged.'}</li>
  <li>${zh ? '共依赖不等于功能相互作用；许多 RWE 基因的依赖来自核糖体 / 线粒体翻译等通用机器。' : 'Co-dependency is not interaction; many RWE dependencies reflect general ribosome / mitochondrial-translation machinery.'}</li>
  <li>${zh ? '谱系分解与谱系限定关联仅用于描述，发现仍以合并模型为准。' : 'Lineage decomposition and lineage-restricted pairs are descriptive; discovery is on the pooled data.'}</li></ul>
  <p>${zh ? '代码与数据表' : 'Code and tables'}: <a href="https://github.com/shuyang-zhao/VisRWE_DepMap" target="_blank" rel="noopener">github.com/shuyang-zhao/VisRWE_DepMap</a> · MIT License · ${zh ? '数据版本' : 'build'} ${esc(D.meta.built)} · v${esc(D.meta.version)}</p>
  <p class="hint">${zh ? 'DepMap 数据由 Broad Institute 以 CC BY 4.0 发布；本页面内嵌的效应向量为量化后的派生数据，仅供浏览。' : 'DepMap data are released by the Broad Institute under CC BY 4.0; the effect vectors embedded here are quantised derived data for browsing only.'}</p>`;
}

// ------------------------------------------------------------------------------------------------ controls
function refilter() { // dataset or filters changed: the explored network is re-derived (partners follow the new calls), other views redrawn
  if (S.mode === 'explore' && S.root != null) { buildEgo(true, false); updateLegend(); } else render();
}
function setDataset(ds) { S.ds = ds; syncSeg('dsSeg', ds); writeHash(); refilter(); }
function setMode(m) {
  S.mode = m; syncSeg('modeSeg', m);
  $('modLevelRow').style.display = m === 'modules' ? '' : 'none';
  $('modeHelp').textContent = t(m + '_help');
  render();
  if (m === 'overview') { const v = fitView(posMap()); if (v) animateTo(null, null, v, 0); }
  writeHash();
}
function onFilter() {
  S.rmin = +$('rmin').value; $('rminVal').textContent = S.rmin.toFixed(2);
  S.pos = $('pos').checked; S.neg = $('neg').checked; S.rweOnly = $('rweOnly').checked; S.hideLib = $('hideLib').checked;
  S.cls = { both: $('cBoth').checked, ng: $('cNg').checked, d2: $('c2d').checked };
  refilter();
}
function searchGene() {
  const q = $('search').value.trim(), i = SYM.get(q.toUpperCase());
  $('searchMsg').textContent = q && i == null ? t('not_found', q) : '';
  if (i != null) { showEgo(i); $('search').blur(); }
}
function lineagePairList() {
  $('lpList').innerHTML = `<table class="list">` + LP.slice(0, 60).map((x, k) => `<tr class="pick" data-k="${k}"><td>${esc(G.sym[x.a])}–${esc(G.sym[x.b])}</td><td><span class="sw" style="background:${D.dict.lcol[D.dict.lineage.indexOf(x.lin)]}"></span>${esc(lab('lineage', x.lin))}</td><td>${fmt(x.r)}</td></tr>`).join('') + '</table>';
  $('lpList').querySelectorAll('tr[data-k]').forEach((tr) => tr.addEventListener('click', () => {
    const x = LP[+tr.dataset.k]; if (S.ds === 'DepMap2D') { S.ds = 'NextGen'; syncSeg('dsSeg', S.ds); }
    showEgo(x.a, [x.b]); if (x.e != null) showEdge(x.e);
  }));
}
function bind() {
  document.querySelectorAll('#dsSeg button').forEach((b) => b.addEventListener('click', () => setDataset(b.dataset.v)));
  document.querySelectorAll('#modeSeg button').forEach((b) => b.addEventListener('click', () => setMode(b.dataset.v)));
  document.querySelectorAll('#labelSeg button').forEach((b) => b.addEventListener('click', () => { S.labels = b.dataset.v; document.querySelectorAll('#labelSeg button').forEach((x) => x.classList.toggle('on', x === b)); refreshLabels(); }));
  document.querySelectorAll('#modLevel button').forEach((b) => b.addEventListener('click', () => { S.modLevel = b.dataset.v; document.querySelectorAll('#modLevel button').forEach((x) => x.classList.toggle('on', x === b)); if (S.mode === 'modules') render(); }));
  ['rmin', 'pos', 'neg', 'rweOnly', 'hideLib', 'cBoth', 'cNg', 'c2d'].forEach((id) => $(id).addEventListener('change', onFilter));
  $('maxExp').addEventListener('input', () => { S.maxExpand = +$('maxExp').value; $('maxExpVal').textContent = S.maxExpand; });
  $('maxExp').addEventListener('change', () => { if (S.mode === 'explore' && S.root != null) buildEgo(true, false); });
  $('rmin').addEventListener('input', () => { $('rminVal').textContent = (+$('rmin').value).toFixed(2); });
  $('labelSize').addEventListener('change', () => { S.labelSize = +$('labelSize').value; restyle(); });
  $('edgeLabels').addEventListener('change', () => { S.edgeLabels = $('edgeLabels').checked; edgesDS.update(edgesDS.getIds().filter((id) => String(id).startsWith('e')).map((id) => edgeObj(+String(id).slice(1)))); });
  $('lpHalo').addEventListener('change', () => { S.lpHalo = $('lpHalo').checked; refilter(); });
  $('undoBtn').addEventListener('click', undo);
  $('resetBtn').addEventListener('click', () => { if (S.mode === 'explore') resetEgo(); else { closePanel(); setMode(S.mode); } });
  $('fitBtn').addEventListener('click', () => { finishAnim(); if (S.mode === 'explore') { S.userMoved = false; recentre(400); } else { const v = fitView(posMap()); if (v) animateTo(null, null, v, 400); } });
  $('layoutBtn').addEventListener('click', () => { if (S.mode === 'explore' && S.root != null) { S.userMoved = false; buildEgo(true, true); } });
  $('pngBtn').addEventListener('click', () => { const c = net.canvas.frame.canvas, a = document.createElement('a'); a.href = c.toDataURL('image/png'); a.download = `visrwe_${S.mode}_${S.ds}.png`; a.click(); });
  $('themeBtn').addEventListener('click', () => { const r = document.documentElement, dark = r.dataset.theme ? r.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    r.dataset.theme = dark ? 'light' : 'dark'; store.set('theme', r.dataset.theme); readColours(); restyle(); });
  $('langBtn').addEventListener('click', () => { S.lang = S.lang === 'zh' ? 'en' : 'zh'; store.set('lang', S.lang); applyI18n(); updateLegend(); updateCounts(); updateCrumbs(); updateEgoBar(); lineagePairList(); if ($('about').style.display === 'flex') $('aboutBox').innerHTML = aboutHTML(); closePanel(); });
  $('aboutBtn').addEventListener('click', () => { $('aboutBox').innerHTML = aboutHTML(); $('about').style.display = 'flex'; });
  $('about').addEventListener('click', (ev) => { if (ev.target === $('about')) $('about').style.display = 'none'; });
  $('search').addEventListener('change', searchGene);
  $('search').addEventListener('keydown', (ev) => { if (ev.key === 'Enter') searchGene(); });
  $('empty').querySelectorAll('[data-g]').forEach((el) => el.addEventListener('click', () => { const i = SYM.get(el.dataset.g); if (i != null) showEgo(i); }));
  $('emptyOverview').addEventListener('click', () => setMode('overview'));
  document.addEventListener('keydown', (ev) => { if (ev.key === 'Escape') { closePanel(); $('about').style.display = 'none'; }
    if ((ev.ctrlKey || ev.metaKey) && ev.key === 'z' && S.mode === 'explore' && !/^(INPUT|TEXTAREA)$/.test(ev.target.tagName)) undo(); });
}

// ------------------------------------------------------------------------------------------------ start
function start() {
  const theme = store.get('theme', ''); if (theme) document.documentElement.dataset.theme = theme;
  readColours(); applyI18n(); bind();
  $('geneList').innerHTML = G.sym.map((s, i) => (G.rwe[i] ? `<option value="${s}">RWE</option>` : `<option value="${s}">`)).join('');
  document.querySelector('#labelSeg button[data-v="auto"]').classList.add('on');
  document.querySelector('#modLevel button[data-v="category"]').classList.add('on');
  $('modLevelRow').style.display = 'none';
  $('rminVal').textContent = S.rmin.toFixed(2); $('maxExpVal').textContent = S.maxExpand;
  createNetwork(); lineagePairList();
  const h = new URLSearchParams(location.hash.slice(1));
  if (['NextGen', 'DepMap2D', 'compare'].includes(h.get('ds'))) S.ds = h.get('ds');
  syncSeg('dsSeg', S.ds); syncSeg('modeSeg', S.mode);
  render();
  $('loading').style.display = 'none';
  const g = h.get('gene') && SYM.get(h.get('gene').toUpperCase()); if (g != null) showEgo(g); // otherwise: empty canvas asking for a gene
  if (typeof DecompressionStream !== 'undefined') { vectors('NextGen'); vectors('DepMap2D'); }
  window.__ready = true;
}
start();
