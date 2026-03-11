/* Oscilloscope Viewer (Web) */

const state = {
  files: [],
  eventsIndex: {},
  events: {},
  eventOrder: [],
  currentEvent: null,
  decim: 20,
  theme: 'dark',
  tool: 'zoom',
  stackVert: true,
  view: {},
  scales: {},
  sgFiltered: {},
  hiddenChannels: new Set(),
  channelMenuChannels: [1, 2, 3, 4, 5, 6, 7, 8],
};

const STORAGE_KEY = 'oscope_web_viewer_session_v2';
let saveTimer = null;

function byId(id){ return document.getElementById(id); }
function clamp(x, a, b){ return Math.max(a, Math.min(b, x)); }
function linspace(a, b, n){ const arr=new Array(n); const d=(b-a)/(n-1); for(let i=0;i<n;i++) arr[i]=a+i*d; return arr; }
function mean(arr){ let s=0; for(const v of arr) s+=v; return s/(arr.length || 1); }
function arrMin(arr){ let m=Infinity; for(const v of arr){ if(v<m) m=v; } return m; }
function arrMax(arr){ let m=-Infinity; for(const v of arr){ if(v>m) m=v; } return m; }

function showBusy(text){
  const el = byId('busy');
  const t = byId('busy-text');
  if(t) t.textContent = text || 'Working…';
  if(el) el.classList.remove('hidden');
}

function hideBusy(){
  const el = byId('busy');
  if(el) el.classList.add('hidden');
}

function saveSessionDebounced(){
  clearTimeout(saveTimer);
  saveTimer = setTimeout(saveSessionNow, 400);
}

function saveSessionNow(){
  try{
    const snapshot = {
      currentEvent: state.currentEvent,
      decim: state.decim,
      theme: state.theme,
      tool: state.tool,
      stackVert: state.stackVert,
      view: state.view,
      hiddenChannels: Array.from(state.hiddenChannels),
      ui: {
        sgChan: byId('sg-chan')?.value || null,
        sgWidth: byId('sg-width')?.value || null,
      },
    };
    localStorage.setItem(STORAGE_KEY, JSON.stringify(snapshot));
  } catch (err){
    console.warn('Session save failed', err);
  }
}

function loadSession(){
  try{
    const raw = localStorage.getItem(STORAGE_KEY);
    if(!raw) return;
    const snap = JSON.parse(raw);
    if(!snap || typeof snap !== 'object') return;

    state.currentEvent = snap.currentEvent || null;
    state.decim = Number.isFinite(snap.decim) ? snap.decim : state.decim;
    state.theme = snap.theme || state.theme;
    state.tool = snap.tool || state.tool;
    state.stackVert = (typeof snap.stackVert === 'boolean') ? snap.stackVert : state.stackVert;
    state.view = snap.view || {};
    state.hiddenChannels = new Set(Array.isArray(snap.hiddenChannels) ? snap.hiddenChannels.map(Number) : []);

    if(byId('decim')) byId('decim').value = String(state.decim);
    if(byId('tool-select')) byId('tool-select').value = state.tool;
    if(byId('stack-vert')) byId('stack-vert').checked = state.stackVert;

    const ui = snap.ui || {};
    if(byId('sg-chan') && ui.sgChan) byId('sg-chan').value = ui.sgChan;
    if(byId('sg-width') && ui.sgWidth) byId('sg-width').value = ui.sgWidth;
  } catch (err){
    console.warn('Session load failed', err);
  }
}

window.addEventListener('beforeunload', saveSessionNow);

const compute = (()=>{
  const WORKERS_OK = (typeof window !== 'undefined' && typeof Worker !== 'undefined' && (location.protocol === 'http:' || location.protocol === 'https:'));
  let worker = null;
  let nextId = 1;
  const pending = new Map();

  function ensure(){
    if(!WORKERS_OK) return null;
    if(worker) return worker;
    try{
      worker = new Worker('worker.js');
      worker.onmessage = (ev)=>{
        const msg = ev.data || {};
        const cb = pending.get(msg.id);
        if(!cb) return;
        pending.delete(msg.id);
        if(msg.ok) cb.resolve(msg.result);
        else cb.reject(new Error(msg.error || 'Worker error'));
      };
    } catch (err){
      console.warn('Worker init failed; using main thread fallback.', err);
      worker = null;
    }
    return worker;
  }

  function call(type, payload, transfer){
    if(!ensure()){
      if(type === 'sg'){
        return Promise.resolve({ y: savgol(payload.y, payload.window, payload.poly || 2) });
      }
      if(type === 'parse_trc'){
        const out = TrcReader.parseBuffer(payload.buf);
        return Promise.resolve({ t: out.t, y: out.y });
      }
      return Promise.reject(new Error('Unknown task: ' + type));
    }

    const id = nextId++;
    return new Promise((resolve, reject)=>{
      pending.set(id, { resolve, reject });
      try{
        worker.postMessage({ id, type, ...payload }, transfer || []);
      } catch (err){
        pending.delete(id);
        reject(err);
      }
    });
  }

  return {
    sg: (y, window, poly=2)=> call('sg', { y, window, poly }),
    parseTrc: (buf)=> call('parse_trc', { buf }, [buf]),
  };
})();

function parseFileName(name){
  const m1 = name.match(/^C(\d+)[-_](\d+)/i);
  if(m1) return { ch: parseInt(m1[1], 10), event: m1[2] };
  const m2 = name.match(/^C(\d+).*-(\d+)\.(?:trc|TRC)$/);
  if(m2) return { ch: parseInt(m2[1], 10), event: m2[2] };
  return null;
}

async function parseCSV(file){
  const text = await file.text();
  const lines = text.split(/\r?\n/).filter(Boolean);
  let start = 0;
  if(lines[0] && lines[0].toLowerCase().includes('time')) start = 1;

  const t = [];
  const y = [];
  for(let i=start;i<lines.length;i++){
    const parts = lines[i].split(/[\t,;\s]+/).filter(Boolean);
    if(parts.length < 2) continue;
    const tt = parseFloat(parts[0]);
    const yy = parseFloat(parts[1]);
    if(Number.isFinite(tt) && Number.isFinite(yy)){
      t.push(tt);
      y.push(yy);
    }
  }
  return { t, y };
}

class TrcReader {
  static findWavedescOffset(buf){
    const bytes = new Uint8Array(buf, 0, Math.min(4096, buf.byteLength));
    const pat = new TextEncoder().encode('WAVEDESC');
    for(let i=0;i<=bytes.length-pat.length;i++){
      let ok = true;
      for(let j=0;j<pat.length;j++){
        if(bytes[i+j] !== pat[j]){ ok = false; break; }
      }
      if(ok) return i;
    }
    return -1;
  }

  static readString(view, pos, len){
    const arr = [];
    for(let i=0;i<len;i++){
      const b = view.getUint8(pos+i);
      if(b === 0) break;
      arr.push(b);
    }
    return new TextDecoder('ascii').decode(new Uint8Array(arr));
  }

  static detectEndian(view, base){
    const be = view.getUint16(base+34, false);
    const le = view.getUint16(base+34, true);
    if(le===1 || be===256) return true;
    if(le===0 || be===0) return false;
    return true;
  }

  static parseBuffer(buf){
    const base = this.findWavedescOffset(buf);
    if(base < 0) throw new Error('WAVEDESC not found');

    const view = new DataView(buf);
    const little = this.detectEndian(view, base);
    const getI32 = (off)=> view.getInt32(base+off, little);
    const getF32 = (off)=> view.getFloat32(base+off, little);
    const getF64 = (off)=> view.getFloat64(base+off, little);

    const commTypeBE = view.getUint16(base+32, false);
    const commTypeLE = view.getUint16(base+32, true);
    const is16 = (commTypeLE===1 || commTypeBE===256 || (commTypeLE!==0 && commTypeBE!==0));

    const lWAVEDESC = getI32(36);
    const lUSERTEXT = getI32(40);
    const lTRIGTIME = getI32(48);
    const lRISTIME = getI32(52);
    const lWAVE1 = getI32(60);

    const VERTICAL_GAIN = getF32(156);
    const VERTICAL_OFFSET = getF32(160);
    const HORIZ_INTERVAL = getF32(176);
    const HORIZ_OFFSET = getF64(180);

    const dataOff = base + lWAVEDESC + lUSERTEXT + lTRIGTIME + lRISTIME;
    const n = lWAVE1 / (is16 ? 2 : 1);

    const yraw = new Array(n);
    if(is16){
      for(let i=0;i<n;i++) yraw[i] = view.getInt16(dataOff + 2*i, little);
    } else {
      for(let i=0;i<n;i++) yraw[i] = view.getInt8(dataOff + i);
    }

    const y = yraw.map(v => VERTICAL_GAIN * v - VERTICAL_OFFSET);
    const t = new Array(n);
    for(let i=0;i<n;i++) t[i] = i * HORIZ_INTERVAL + HORIZ_OFFSET;

    return { t, y, d: { TEMPLATE_NAME: this.readString(view, base+16, 16) } };
  }
}

function savgol(y, window, poly){
  window = (window|0);
  if(window < 5) window = 5;
  if(window % 2 === 0) window += 1;

  const n = y.length;
  const half = (window - 1) / 2;
  const out = new Array(n).fill(0);

  function weights(m){
    const X = [];
    for(let i=-m;i<=m;i++) X.push([1, i, i*i]);

    const A = [[0,0,0],[0,0,0],[0,0,0]];
    for(const row of X){
      for(let i=0;i<3;i++){
        for(let j=0;j<3;j++) A[i][j] += row[i] * row[j];
      }
    }

    const [a,b,c,d,e,f,g,h,k] = [A[0][0],A[0][1],A[0][2],A[1][0],A[1][1],A[1][2],A[2][0],A[2][1],A[2][2]];
    const A_ = e*k-f*h;
    const B_ = -(d*k-f*g);
    const C_ = d*h-e*g;
    const D_ = -(b*k-c*h);
    const E_ = a*k-c*g;
    const F_ = -(a*h-b*g);
    const G_ = b*f-c*e;
    const H_ = -(a*f-b*e);
    const K_ = a*e-b*d;
    const det = a*A_ + b*B_ + c*C_;
    const inv = [[A_/det, D_/det, G_/det], [B_/det, E_/det, H_/det], [C_/det, F_/det, K_/det]];
    const c0 = [inv[0][0], inv[1][0], inv[2][0]];

    const W = [];
    for(let i=-m;i<=m;i++) W.push(c0[0] + c0[1]*i + c0[2]*i*i);
    return W;
  }

  const W = weights(half);
  for(let i=0;i<n;i++){
    let acc=0;
    let ws=0;
    for(let k=-half;k<=half;k++){
      let idx = i + k;
      if(idx<0) idx = 0;
      if(idx>=n) idx = n-1;
      const w = W[k+half];
      acc += w * y[idx];
      ws += w;
    }
    out[i] = acc / ws;
  }

  return out;
}

function getPlotColors(){
  const dark = state.theme !== 'light';
  return {
    bg: dark ? '#000' : '#fff',
    axes: dark ? '#666' : '#666',
    raw: dark ? '#ffffff' : '#000000',
    sg: dark ? '#5fd35f' : '#2ca02c',
    grid: dark ? '#222' : '#ddd',
    text: dark ? '#bbb' : '#333',
  };
}

function getPlotPadding(){
  return { padL: 55, padR: 12, padT: 10, padB: 28 };
}

function niceTicks(min, max, maxTicks=6){
  const span = Math.max(1e-20, max - min);
  const rough = span / Math.max(1, maxTicks);
  const pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
  const candidates = [1,2,5,10].map(m => m*pow10);
  let step = candidates[0];
  for(const c of candidates){
    if(Math.abs(c - rough) < Math.abs(step - rough)) step = c;
  }
  const niceMin = Math.ceil(min/step) * step;
  const niceMax = Math.floor(max/step) * step;
  const ticks = [];
  for(let v=niceMin; v<=niceMax + 0.5*step; v+=step) ticks.push(v);
  return { ticks, step };
}

function fmtTick(v){
  const av = Math.abs(v);
  if((av > 0 && av < 1e-3) || av >= 1e4) return v.toExponential(1);
  const d = av >= 1 ? 2 : 4;
  return v.toFixed(d).replace(/\.0+$/, '').replace(/\.$/, '');
}

function setCanvasSize(canvas, widthPx, heightPx){
  const dpr = window.devicePixelRatio || 1;
  widthPx = Math.max(200, Math.floor(widthPx));
  heightPx = Math.max(80, Math.floor(heightPx));
  canvas.style.width = widthPx + 'px';
  canvas.style.height = heightPx + 'px';
  canvas.width = Math.round(widthPx * dpr);
  canvas.height = Math.round(heightPx * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}

function createCanvasWrap(ch){
  const wrap = document.createElement('div');
  wrap.className = 'plot-item';

  const header = document.createElement('div');
  header.className = 'plot-header';
  header.innerHTML = `<strong>Channel ${ch}</strong>`;

  const canvas = document.createElement('canvas');
  canvas.className = 'plot-canvas';
  canvas.dataset.ch = String(ch);

  attachCanvasInteraction(canvas, ch);

  wrap.appendChild(header);
  wrap.appendChild(canvas);
  return { wrap, canvas };
}

function attachCanvasInteraction(canvas, ch){
  let active = false;
  let x0 = 0;
  let y0 = 0;
  let lastX = 0;
  let lastY = 0;

  const getScale = ()=> state.scales[`${state.currentEvent}:ch${ch}`];

  const toT = (xp, scale)=>{
    const { tmin, tmax, padL, padR, width } = scale;
    const L = padL ?? 55;
    const R = padR ?? 12;
    return tmin + clamp((xp - L) / (width - L - R), 0, 1) * (tmax - tmin);
  };

  const toV = (yp, scale)=>{
    const { ymin, ymax, padT, padB, height } = scale;
    const T = padT ?? 10;
    const B = padB ?? 28;
    return ymax - clamp((yp - T) / (height - T - B), 0, 1) * (ymax - ymin);
  };

  const pointerPos = (e)=>{
    const r = canvas.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  };

  canvas.addEventListener('pointerdown', (e)=>{
    canvas.setPointerCapture(e.pointerId);
    const p = pointerPos(e);
    active = true;
    x0 = p.x;
    y0 = p.y;
    lastX = p.x;
    lastY = p.y;
  });

  canvas.addEventListener('pointermove', (e)=>{
    if(!active) return;
    if(state.tool !== 'pan') return;

    const scale = getScale();
    if(!scale) return;

    const p = pointerPos(e);
    const dx = p.x - lastX;
    const dy = p.y - lastY;

    const w = scale.width - scale.padL - scale.padR;
    const h = scale.height - scale.padT - scale.padB;
    if(w <= 0 || h <= 0) return;

    const dt = dx / w * (scale.tmax - scale.tmin);
    const dv = dy / h * (scale.ymax - scale.ymin);

    const evk = `${state.currentEvent}:ch${ch}`;
    const vr = state.view[evk] || { tmin: scale.tmin, tmax: scale.tmax, ymin: scale.ymin, ymax: scale.ymax };

    if(!e.shiftKey){
      vr.tmin -= dt;
      vr.tmax -= dt;
    }
    if(e.shiftKey){
      vr.ymin += dv;
      vr.ymax += dv;
    }

    state.view[evk] = vr;
    lastX = p.x;
    lastY = p.y;
    plotWaveforms();
  });

  canvas.addEventListener('pointerup', (e)=>{
    if(!active) return;
    active = false;

    if(state.tool !== 'zoom'){
      saveSessionDebounced();
      return;
    }

    const scale = getScale();
    if(!scale) return;

    const p = pointerPos(e);
    const tx0 = Math.min(x0, p.x);
    const tx1 = Math.max(x0, p.x);
    const ty0 = Math.min(y0, p.y);
    const ty1 = Math.max(y0, p.y);

    if(Math.abs(tx1 - tx0) < 3 || Math.abs(ty1 - ty0) < 3){
      saveSessionDebounced();
      return;
    }

    const tA = toT(tx0, scale);
    const tB = toT(tx1, scale);
    const vA = toV(ty1, scale);
    const vB = toV(ty0, scale);

    if(!(tB > tA) || !(vB > vA)){
      saveSessionDebounced();
      return;
    }

    const evk = `${state.currentEvent}:ch${ch}`;
    state.view[evk] = { tmin: tA, tmax: tB, ymin: vA, ymax: vB };
    plotWaveforms();
    saveSessionDebounced();
  });

  canvas.addEventListener('wheel', (e)=>{
    const scale = getScale();
    if(!scale) return;
    e.preventDefault();

    const r = canvas.getBoundingClientRect();
    const x = e.clientX - r.left;
    const y = e.clientY - r.top;

    const factor = (e.deltaY > 0) ? 1.1 : 0.9;
    const tCur = toT(x, scale);
    const vCur = toV(y, scale);

    const evk = `${state.currentEvent}:ch${ch}`;
    const vr = state.view[evk] || { tmin: scale.tmin, tmax: scale.tmax, ymin: scale.ymin, ymax: scale.ymax };

    const doX = !e.shiftKey;
    const doY = e.shiftKey || e.ctrlKey;

    if(doX){
      const span = vr.tmax - vr.tmin;
      const ns = span * factor;
      const a = (tCur - vr.tmin) / span;
      vr.tmin = tCur - a * ns;
      vr.tmax = vr.tmin + ns;
    }

    if(doY){
      const span = vr.ymax - vr.ymin;
      const ns = span * factor;
      const a = (vCur - vr.ymin) / span;
      vr.ymin = vCur - a * ns;
      vr.ymax = vr.ymin + ns;
    }

    state.view[evk] = vr;
    plotWaveforms();
    saveSessionDebounced();
  }, { passive: false });

  canvas.addEventListener('dblclick', ()=>{
    const evk = `${state.currentEvent}:ch${ch}`;
    delete state.view[evk];
    plotWaveforms();
    saveSessionDebounced();
  });
}

function updateCanvasCursors(){
  document.querySelectorAll('.plot-canvas').forEach((canvas)=>{
    canvas.classList.toggle('pan', state.tool === 'pan');
    canvas.classList.toggle('zoom', state.tool === 'zoom');
  });
}

function visibleChannelsForEvent(ev){
  const obj = state.events[ev];
  if(!obj || !obj.channels) return [];
  const all = Object.keys(obj.channels).map(Number).sort((a,b)=>a-b);
  return all.filter(ch => !state.hiddenChannels.has(ch));
}

function formatChannelVisibilityIndicator(){
  const hidden = Array.from(state.hiddenChannels).sort((a,b)=>a-b);
  if(hidden.length === 0) return '';
  return ' • ' + hidden.join(',');
}

function updateChannelVisibilitySummary(){
  const summary = byId('channel-menu-summary');
  if(!summary) return;
  const suffix = formatChannelVisibilityIndicator();
  summary.textContent = 'Channels' + suffix;
  if(!suffix) summary.title = 'Toggle which channel plots are shown';
  else summary.title = `Hidden channels: ${suffix.replace(' • ', '')}`;
}

function rebuildChannelVisibilityMenu(channels){
  const nextChannels = (channels && channels.length)
    ? Array.from(new Set(channels.map(Number))).sort((a,b)=>a-b)
    : [1, 2, 3, 4, 5, 6, 7, 8];
  state.channelMenuChannels = nextChannels;
  const allowed = new Set(nextChannels);
  state.hiddenChannels = new Set(Array.from(state.hiddenChannels).filter(ch => allowed.has(ch)));

  const list = byId('channel-toggle-list');
  if(!list) return;
  list.innerHTML = '';

  for(const ch of nextChannels){
    const row = document.createElement('label');
    row.className = 'channel-toggle-row';

    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = !state.hiddenChannels.has(ch);
    cb.addEventListener('change', ()=>{
      if(cb.checked) state.hiddenChannels.delete(ch);
      else state.hiddenChannels.add(ch);
      updateChannelVisibilitySummary();
      plotWaveforms();
      saveSessionDebounced();
    });

    const txt = document.createElement('span');
    txt.textContent = `Channel ${ch}`;

    row.appendChild(cb);
    row.appendChild(txt);
    list.appendChild(row);
  }

  updateChannelVisibilitySummary();
}

function setHiddenChannels(channels){
  state.hiddenChannels = new Set(Array.from(channels || []).map(Number));
  rebuildChannelVisibilityMenu(state.channelMenuChannels);
  plotWaveforms();
  saveSessionDebounced();
}

function showAllChannels(){
  setHiddenChannels(new Set());
}

function hideAllChannels(){
  setHiddenChannels(new Set(state.channelMenuChannels));
}

async function loadEvent(ev){
  if(!ev || !state.eventsIndex[ev]) return;
  if(state.events[ev] && state.events[ev].channels) return;

  const channels = state.eventsIndex[ev].channels;
  const obj = { channels: {} };

  for(const ch of Object.keys(channels)){
    const file = channels[ch];
    let t = [];
    let y = [];

    if(/\.(trc)$/i.test(file.name)){
      try{
        const buf = await file.arrayBuffer();
        try{
          const out = await compute.parseTrc(buf);
          t = out.t;
          y = out.y;
        } catch (_workerErr){
          const out = TrcReader.parseBuffer(buf);
          t = out.t;
          y = out.y;
        }
      } catch (err){
        console.warn('TRC parse failed', file.name, err);
        continue;
      }
    } else {
      const out = await parseCSV(file);
      t = out.t;
      y = out.y;
    }

    const ym = mean(y);
    obj.channels[parseInt(ch, 10)] = { t, y: y.map(v => v - ym) };
  }

  state.events[ev] = obj;
}

async function getEventChannelData(ev, ch){
  if(!state.events[ev]) await loadEvent(ev);
  const evt = state.events[ev];
  if(!evt) return null;
  return evt.channels[ch] || null;
}

function refreshEventSelect(){
  const sel = byId('event-select');
  sel.innerHTML = '';
  for(const ev of state.eventOrder){
    const opt = document.createElement('option');
    opt.value = ev;
    opt.textContent = ev;
    sel.appendChild(opt);
  }
  if(state.currentEvent && state.eventOrder.includes(state.currentEvent)){
    sel.value = state.currentEvent;
  } else if(state.eventOrder.length){
    state.currentEvent = state.eventOrder[0];
    sel.value = state.currentEvent;
  }
}

function refreshChannelCombos(){
  const ev = state.currentEvent;
  const channelsIdx = ev && state.eventsIndex[ev] ? state.eventsIndex[ev].channels : {};
  const arr = Object.keys(channelsIdx).map(Number).sort((a,b)=>a-b);

  const sgChan = byId('sg-chan');
  sgChan.innerHTML = '';
  for(const ch of arr){
    const opt = document.createElement('option');
    opt.value = String(ch);
    opt.textContent = String(ch);
    sgChan.appendChild(opt);
  }

  if(arr.length){
    if(!arr.includes(parseInt(sgChan.value, 10))) sgChan.value = String(arr[0]);
  }

  rebuildChannelVisibilityMenu(arr);
}

function plotWaveforms(){
  const ev = state.currentEvent;
  if(!ev) return;

  const obj = state.events[ev];
  if(!obj) return;

  const container = byId('plots');
  container.innerHTML = '';

  const channels = visibleChannelsForEvent(ev);
  if(channels.length === 0){
    const empty = document.createElement('div');
    empty.className = 'plot-item';
    empty.textContent = 'All channels are hidden for this event. Use Channels → Show All.';
    container.appendChild(empty);
    return;
  }

  state.scales = {};

  const mainEl = document.querySelector('main');
  const viewportH = window.innerHeight || document.documentElement.clientHeight || 800;
  const availH = Math.max(140, (mainEl ? mainEl.clientHeight : (viewportH - (document.querySelector('header')?.offsetHeight || 120))) - 8);
  const containerW = container.clientWidth || (mainEl?.clientWidth || 1000);

  const gapY = 4;
  const overhead = 30;
  let perCanvasH = state.stackVert
    ? Math.max(70, Math.floor((availH - (channels.length-1)*gapY - channels.length*overhead) / Math.max(1, channels.length)))
    : 260;

  if(state.stackVert) perCanvasH = Math.max(70, Math.floor(perCanvasH * 0.95));

  for(const ch of channels){
    const { wrap, canvas } = createCanvasWrap(ch);
    container.appendChild(wrap);

    const g = setCanvasSize(canvas, containerW - 24, perCanvasH);
    const rect = canvas.getBoundingClientRect();
    const w = Math.floor(rect.width);
    const h = Math.floor(rect.height);

    const { padL, padR, padT, padB } = getPlotPadding();
    const colors = getPlotColors();

    g.clearRect(0, 0, w, h);
    g.fillStyle = colors.bg;
    g.fillRect(0, 0, w, h);

    const data = obj.channels[ch];
    if(!data || !data.t || data.t.length === 0) continue;

    const t = data.t;
    const y = data.y;

    let tmin = t[0];
    let tmax = t[t.length-1];
    let ymin = arrMin(y);
    let ymax = arrMax(y);

    const sgKey = `${ev}:ch${ch}`;
    const sgv = state.sgFiltered[sgKey];
    if(sgv && Array.isArray(sgv.y) && sgv.y.length){
      const sgMin = arrMin(sgv.y);
      const sgMax = arrMax(sgv.y);
      if(sgMin < ymin) ymin = sgMin;
      if(sgMax > ymax) ymax = sgMax;
    }

    const evk = `${ev}:ch${ch}`;
    const vr = state.view[evk];
    if(vr){
      if(Number.isFinite(vr.tmin) && Number.isFinite(vr.tmax) && vr.tmax > vr.tmin){
        tmin = vr.tmin;
        tmax = vr.tmax;
      }
      if(Number.isFinite(vr.ymin) && Number.isFinite(vr.ymax) && vr.ymax > vr.ymin){
        ymin = vr.ymin;
        ymax = vr.ymax;
      }
    }

    if(ymax === ymin){
      ymax += 1;
      ymin -= 1;
    }

    state.scales[evk] = { tmin, tmax, ymin, ymax, padL, padR, padT, padB, width: w, height: h };

    const tx = (tt)=> padL + (tt - tmin)/(tmax-tmin) * (w - padL - padR);
    const ty = (vv)=> (h - padB) - (vv - ymin)/(ymax-ymin) * (h - padT - padB);

    const xTicks = niceTicks(tmin, tmax, 6).ticks;
    const yTicks = niceTicks(ymin, ymax, 5).ticks;

    g.strokeStyle = colors.grid;
    g.lineWidth = 1;
    g.setLineDash([2,4]);
    g.beginPath();
    for(const xt of xTicks){
      const X = tx(xt);
      g.moveTo(X, padT);
      g.lineTo(X, h - padB);
    }
    for(const yt of yTicks){
      const Y = ty(yt);
      g.moveTo(padL, Y);
      g.lineTo(w - padR, Y);
    }
    g.stroke();
    g.setLineDash([]);

    g.strokeStyle = colors.axes;
    g.lineWidth = 1;
    g.beginPath();
    g.moveTo(padL, padT);
    g.lineTo(padL, h-padB);
    g.lineTo(w-padR, h-padB);
    g.stroke();

    g.fillStyle = colors.text;
    g.font = '12px sans-serif';
    g.textAlign = 'center';
    g.textBaseline = 'top';
    for(const xt of xTicks){ g.fillText(fmtTick(xt), tx(xt), h - padB + 2); }

    g.textAlign = 'right';
    g.textBaseline = 'middle';
    for(const yt of yTicks){ g.fillText(fmtTick(yt), padL - 6, ty(yt)); }

    g.textAlign = 'center';
    g.textBaseline = 'bottom';
    g.fillText('t (s)', (padL + (w - padR))/2, h - 4);

    g.save();
    g.translate(12, padT + (h - padT - padB)/2);
    g.rotate(-Math.PI/2);
    g.textAlign = 'center';
    g.textBaseline = 'top';
    g.fillText('y', 0, 0);
    g.restore();

    const dec = Math.max(1, state.decim|0);
    g.strokeStyle = colors.raw;
    g.lineWidth = 1.2;
    g.beginPath();
    for(let i=0;i<t.length;i+=dec){
      const X = tx(t[i]);
      const Y = ty(y[i]);
      if(i===0) g.moveTo(X,Y);
      else g.lineTo(X,Y);
    }
    g.stroke();

    if(sgv && Array.isArray(sgv.t) && Array.isArray(sgv.y) && sgv.t.length){
      g.strokeStyle = colors.sg;
      g.lineWidth = 1.4;
      g.beginPath();
      for(let i=0;i<sgv.t.length;i+=dec){
        const X = tx(sgv.t[i]);
        const Y = ty(sgv.y[i]);
        if(i===0) g.moveTo(X,Y);
        else g.lineTo(X,Y);
      }
      g.stroke();
    }
  }

  updateCanvasCursors();
}

function applyTheme(){
  if(state.theme === 'light'){
    document.documentElement.style.setProperty('--bg', '#f0f0f0');
    document.documentElement.style.setProperty('--fg', '#222');
    document.documentElement.style.setProperty('--panel', '#ffffff');
    document.documentElement.style.setProperty('--btn-bg', '#e0e0e0');
    document.documentElement.style.setProperty('--border', '#cfcfcf');
    document.documentElement.style.setProperty('--accent', '#0d6efd');
    if(byId('theme-btn')) byId('theme-btn').textContent = 'Dark Mode';
  } else {
    document.documentElement.style.setProperty('--bg', '#1e1e1e');
    document.documentElement.style.setProperty('--fg', '#f0f0f0');
    document.documentElement.style.setProperty('--panel', '#2b2b2b');
    document.documentElement.style.setProperty('--btn-bg', '#333333');
    document.documentElement.style.setProperty('--border', '#444');
    document.documentElement.style.setProperty('--accent', '#00aaff');
    if(byId('theme-btn')) byId('theme-btn').textContent = 'Light Mode';
  }
}

async function changeEvent(ev){
  state.currentEvent = ev;
  await loadEvent(ev);
  refreshChannelCombos();
  plotWaveforms();
  saveSessionDebounced();
}

function initUI(){
  byId('file-input').addEventListener('change', async (e)=>{
    state.files = Array.from(e.target.files || []).filter(f => /\.(csv|txt|trc)$/i.test(f.name));
    state.eventsIndex = {};
    state.events = {};
    state.eventOrder = [];
    state.currentEvent = null;
    state.sgFiltered = {};

    for(const file of state.files){
      const meta = parseFileName(file.name);
      if(!meta) continue;
      if(!state.eventsIndex[meta.event]) state.eventsIndex[meta.event] = { channels: {} };
      state.eventsIndex[meta.event].channels[meta.ch] = file;
    }

    state.eventOrder = Object.keys(state.eventsIndex).sort((a,b)=>parseInt(a,10)-parseInt(b,10));

    refreshEventSelect();
    if(state.currentEvent){
      showBusy('Loading event…');
      try{
        await loadEvent(state.currentEvent);
      } finally {
        hideBusy();
      }
    }
    refreshChannelCombos();
    plotWaveforms();
    saveSessionDebounced();
  });

  byId('event-select').addEventListener('change', async (e)=>{
    await changeEvent(e.target.value);
  });

  byId('prev-btn').addEventListener('click', async ()=>{
    const arr = state.eventOrder;
    const idx = arr.indexOf(state.currentEvent);
    if(idx > 0){
      byId('event-select').value = arr[idx-1];
      await changeEvent(arr[idx-1]);
    }
  });

  byId('next-btn').addEventListener('click', async ()=>{
    const arr = state.eventOrder;
    const idx = arr.indexOf(state.currentEvent);
    if(idx >= 0 && idx < arr.length - 1){
      byId('event-select').value = arr[idx+1];
      await changeEvent(arr[idx+1]);
    }
  });

  byId('decim').addEventListener('input', (e)=>{
    state.decim = clamp(parseInt(e.target.value, 10) || 1, 1, 200);
    plotWaveforms();
    saveSessionDebounced();
  });

  byId('stack-vert').addEventListener('change', (e)=>{
    state.stackVert = !!e.target.checked;
    plotWaveforms();
    saveSessionDebounced();
  });

  byId('tool-select').addEventListener('change', (e)=>{
    state.tool = e.target.value;
    updateCanvasCursors();
    saveSessionDebounced();
  });

  byId('reset-view-btn').addEventListener('click', ()=>{
    const ev = state.currentEvent;
    if(!ev) return;
    const keys = Object.keys(state.view).filter(k => k.startsWith(`${ev}:ch`));
    for(const k of keys) delete state.view[k];
    plotWaveforms();
    saveSessionDebounced();
  });

  byId('autoscale-y-btn').addEventListener('click', ()=>{
    const ev = state.currentEvent;
    if(!ev || !state.events[ev]) return;

    const obj = state.events[ev];
    for(const chStr of Object.keys(obj.channels)){
      const ch = parseInt(chStr, 10);
      if(state.hiddenChannels.has(ch)) continue;
      const evk = `${ev}:ch${ch}`;
      const scale = state.scales[evk];
      if(!scale) continue;

      const data = obj.channels[ch];
      const idx = data.t.map((tt, i)=> (tt >= scale.tmin && tt <= scale.tmax) ? i : -1).filter(i => i >= 0);
      if(idx.length < 2) continue;

      const yy = idx.map(i => data.y[i]);
      const sgKey = `${ev}:ch${ch}`;
      const sgv = state.sgFiltered[sgKey];
      if(sgv && sgv.y && sgv.y.length){
        const sgy = idx.map(i => sgv.y[i]).filter(v => Number.isFinite(v));
        yy.push(...sgy);
      }

      let ymin = arrMin(yy);
      let ymax = arrMax(yy);
      if(ymax === ymin){ ymax += 1; ymin -= 1; }

      state.view[evk] = { ...(state.view[evk] || {}), ymin, ymax, tmin: scale.tmin, tmax: scale.tmax };
    }

    plotWaveforms();
    saveSessionDebounced();
  });

  byId('theme-btn').addEventListener('click', ()=>{
    state.theme = (state.theme === 'dark') ? 'light' : 'dark';
    applyTheme();
    plotWaveforms();
    saveSessionDebounced();
  });

  byId('show-all-channels').addEventListener('click', (e)=>{
    e.preventDefault();
    showAllChannels();
  });

  byId('hide-all-channels').addEventListener('click', (e)=>{
    e.preventDefault();
    hideAllChannels();
  });

  byId('sg-btn').addEventListener('click', async ()=>{
    const ev = state.currentEvent;
    if(!ev) return;

    const ch = parseInt(byId('sg-chan').value, 10);
    const data = await getEventChannelData(ev, ch);
    if(!data) return;

    const w = parseInt(byId('sg-width').value, 10) || 200;

    showBusy('Filtering…');
    try{
      const res = await compute.sg(data.y, w, 2);
      const yy = (res && res.y) ? res.y : savgol(data.y, w, 2);
      state.sgFiltered[`${ev}:ch${ch}`] = { t: [...data.t], y: yy };
    } finally {
      hideBusy();
    }

    plotWaveforms();
    saveSessionDebounced();
  });

  byId('sg-clear-btn').addEventListener('click', ()=>{
    const ev = state.currentEvent;
    if(!ev) return;
    const ch = parseInt(byId('sg-chan').value, 10);
    delete state.sgFiltered[`${ev}:ch${ch}`];
    plotWaveforms();
    saveSessionDebounced();
  });
  byId('sg-chan').addEventListener('change', saveSessionDebounced);
  byId('sg-width').addEventListener('change', saveSessionDebounced);

  byId('help-btn').addEventListener('click', ()=>{
    alert('Load CSV or TRC data (e.g. C1-00012.trc). This viewer is waveform-only (fitting removed).\n\nTools: Zoom (drag rectangle), Pan (drag; hold Shift for vertical). Scroll wheel zooms X; hold Shift/Ctrl to zoom Y. Double-click resets a plot. Reset View clears channel views for current event.\n\nUse the Channels menu to show/hide plots by channel.');
  });

  document.addEventListener('keydown', (e)=>{
    if(e.target && (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || e.target.isContentEditable)) return;

    const k = e.key.toLowerCase();
    if(k === ','){
      byId('prev-btn').click();
      e.preventDefault();
    } else if(k === '.'){
      byId('next-btn').click();
      e.preventDefault();
    } else if(k === 's'){
      byId('sg-btn').click();
      e.preventDefault();
    } else if(k === '0'){
      showAllChannels();
      e.preventDefault();
    } else if(k >= '1' && k <= '8'){
      const chNum = parseInt(k, 10);
      const sg = byId('sg-chan');
      const has = Array.from(sg.options).some(o => parseInt(o.value, 10) === chNum);
      if(has){
        sg.value = String(chNum);
        e.preventDefault();
      }
    }
  });

  window.addEventListener('resize', ()=> plotWaveforms());
}

function startup(){
  loadSession();
  applyTheme();
  rebuildChannelVisibilityMenu();
  initUI();
  updateCanvasCursors();
}

startup();
