/* Oscilloscope Palantir (Web) */

// Simple state
const state = {
  files: [],            // File objects from input
  eventsIndex: {},      // {eventKey: {channels: {ch: File}}} (no parsing yet)
  events: {},           // loaded data: {eventKey: {channels: {ch:{t,y}}}}
  eventOrder: [],       // [eventKey]
  currentEvent: null,
  decim: 20,
  results: {},          // fit results keyed like evt_<id>_<fit>_ch<ch>
  fitRegion: {},        // per-channel [t0,t1]
  sgFiltered: {},
  theme: 'dark',
  featureFilter: null,  // array of filtered event keys
  scales: {},           // per-ev:ch -> {tmin,tmax,ymin,ymax}
  fitResolution: 2000,
  datasetId: null,
};

// Compute worker wrapper to offload DSP/fits off the main thread
const compute = (function(){
  let w=null; let nextId=1; const pending=new Map();
  function ensure(){
    if(!w){
      try { w = new Worker('worker.js'); }
      catch (e) { console.warn('Worker init failed, using main-thread fallback', e); w=null; }
      if(w){
        w.onmessage = (ev)=>{
          const {id, ok, result, error} = ev.data||{};
          const p = pending.get(id); if(!p) return; pending.delete(id);
          if(ok) p.resolve(result); else p.reject(new Error(error||'Worker error'));
        };
      }
    }
    return w;
  }
  function call(type, payload){
    if(!ensure()){
      // Fallback path executes synchronously on main thread
      return new Promise((resolve, reject)=>{
        try{
          if(type==='sg'){
            const {y, window, poly} = payload; resolve({y: savgol(y, window, poly||2)});
          } else if(type==='fit'){
            const {fitName, t, y, p0} = payload; const model = FIT_LIB[fitName]; if(!model||!model.func) throw new Error('No model'); resolve({params: lsqFit(model, t, y, p0)});
          } else if(type==='low_pass_max'){
            const {t, y, width, invert} = payload; const wOdd = (width%2===0)? width+1: width; const yproc = invert? y.map(v=>-v): y.slice(); const ys = savgol(yproc, wOdd, 2); const idx=argMax(ys); const t_at=t[Math.max(0, Math.min(idx, t.length-1))]; const val = invert? -ys[idx]: ys[idx]; resolve({value: val, t_at, width: wOdd});
          } else { reject(new Error('Unknown compute call: '+type)); }
        } catch (e){ reject(e); }
      });
    }
    const id = nextId++;
    return new Promise((resolve, reject)=>{ pending.set(id, {resolve, reject}); try{ w.postMessage({id, type, ...payload}); } catch(e){ pending.delete(id); reject(e);} });
  }
  return {
    fit: (fitName, t, y, p0)=> call('fit', {fitName, t, y, p0}),
    sg:  (y, width, poly=2)=> call('sg', {y, window: width, poly}),
    lowPassMax: (t, y, width, invert)=> call('low_pass_max', {t, y, width, invert}),
  };
})();

// Utilities
function byId(id){ return document.getElementById(id); }
function clamp(x,a,b){ return Math.max(a, Math.min(b, x)); }
function linspace(a,b,n){ const arr=new Array(n); const d=(b-a)/(n-1); for(let i=0;i<n;i++)arr[i]=a+i*d; return arr; }
function mean(arr){ let s=0; for(const v of arr) s+=v; return s/arr.length; }
function stdev(arr){ const m=mean(arr); let s=0; for(const v of arr){ const d=v-m; s+=d*d;} return Math.sqrt(s/arr.length); }
function arrMax(arr){ let m=-Infinity; for(let i=0;i<arr.length;i++){ const v=arr[i]; if(v>m) m=v; } return m; }
function arrMin(arr){ let m=Infinity; for(let i=0;i<arr.length;i++){ const v=arr[i]; if(v<m) m=v; } return m; }
function argMax(arr){ let m=-Infinity, idx=-1; for(let i=0;i<arr.length;i++){ const v=arr[i]; if(v>m){ m=v; idx=i; } } return idx; }
function arrAbsMax(arr){ let m=-Infinity; for(let i=0;i<arr.length;i++){ const v=Math.abs(arr[i]); if(v>m) m=v; } return m; }

// Parse CSV file (time, value) with or without header
async function parseCSV(file){
  const text = await file.text();
  const lines = text.split(/\r?\n/).filter(Boolean);
  let start=0;
  if(lines[0].toLowerCase().includes('time')) start=1;
  const t=[]; const y=[];
  for(let i=start;i<lines.length;i++){
    const parts = lines[i].split(/[\t,;\s]+/).filter(Boolean);
    if(parts.length<2) continue;
    const tt = parseFloat(parts[0]);
    const yy = parseFloat(parts[1]);
    if(Number.isFinite(tt) && Number.isFinite(yy)){
      t.push(tt); y.push(yy);
    }
  }
  return {t,y};
}

// Parse LeCroy .trc file into {t,y,d}
class TrcReader {
  static findWavedescOffset(buf){
    const bytes = new Uint8Array(buf, 0, Math.min(4096, buf.byteLength));
    const pat = new TextEncoder().encode('WAVEDESC');
    for(let i=0;i<=bytes.length-pat.length;i++){
      let ok=true; for(let j=0;j<pat.length;j++){ if(bytes[i+j]!==pat[j]){ ok=false; break; } }
      if(ok) return i;
    }
    return -1;
  }
  static readString(view, pos, len){
    const arr=[]; for(let i=0;i<len;i++){ const b=view.getUint8(pos+i); if(b===0) break; arr.push(b); }
    return new TextDecoder('ascii').decode(new Uint8Array(arr));
  }
  static detectEndian(view, base){
    // COMM_ORDER at offset 34: 0=big, 1=little
    const be = view.getUint16(base+34, false);
    const le = view.getUint16(base+34, true);
    if(le===1 || be===256) return true; // little-endian
    if(le===0 || be===0) return false; // big-endian
    // default to little
    return true;
  }
  static async parseFile(file){
    const buf = await file.arrayBuffer();
    const base = this.findWavedescOffset(buf);
    if(base<0) throw new Error('WAVEDESC not found');
    const view = new DataView(buf);
    // Template name (debug only)
    const template = this.readString(view, base+16, 16);
    // Sample format at 32: 0=byte, 1=word
    const commTypeBE = view.getUint16(base+32, false);
    const commTypeLE = view.getUint16(base+32, true);
    // If either equals 1 or 256, treat as 16-bit
    const is16 = (commTypeLE===1 || commTypeBE===256 || (commTypeLE!==0 && commTypeBE!==0));
    const little = this.detectEndian(view, base);
    const getI32 = (off)=> view.getInt32(base+off, little);
    const getI16 = (off)=> view.getInt16(base+off, little);
    const getF32 = (off)=> view.getFloat32(base+off, little);
    const getF64 = (off)=> view.getFloat64(base+off, little);
    // lengths
    const lWAVEDESC = getI32(36);
    const lUSERTEXT = getI32(40);
    const lTRIGTIME = getI32(48);
    const lRISTIME  = getI32(52);
    const lWAVE1    = getI32(60);
    // waveform info
    const VERTICAL_GAIN   = getF32(156);
    const VERTICAL_OFFSET = getF32(160);
    const HORIZ_INTERVAL  = getF32(176);
    const HORIZ_OFFSET    = getF64(180);
    // sample array start
    const dataOff = base + lWAVEDESC + lUSERTEXT + lTRIGTIME + lRISTIME;
    const n = lWAVE1 / (is16? 2: 1);
    const yraw = new Array(n);
    if(is16){ for(let i=0;i<n;i++){ yraw[i]=view.getInt16(dataOff + 2*i, little); } }
    else{ for(let i=0;i<n;i++){ yraw[i]=view.getInt8(dataOff + i); } }
    // scale
    const y = yraw.map(v=> VERTICAL_GAIN * v - VERTICAL_OFFSET);
    const t = new Array(n); for(let i=0;i<n;i++){ t[i] = (i+1)*HORIZ_INTERVAL + HORIZ_OFFSET; }
    return {t,y,d:{ TEMPLATE_NAME: template }};
  }
}

// File grouping by name: supports CSV and TRC naming
function parseFileName(name){
  // examples: C1-00012.csv, C4-5.csv, C1-00012.trc, C4_xxx-5.trc
  const m1 = name.match(/^C(\d+)[-_](\d+)/i);
  if(m1) return { ch: parseInt(m1[1],10), event: m1[2] };
  const m2 = name.match(/^C(\d+).*-(\d+)\.(?:trc|TRC)$/);
  if(m2) return { ch: parseInt(m2[1],10), event: m2[2] };
  return null;
}

// Simple canvas plotting
function createOrGetCanvas(container, ch){
  const id = `plot-ch${ch}`;
  let wrap = document.getElementById(`${id}-wrap`);
  if(!wrap){
    wrap = document.createElement('div'); wrap.className='plot-item'; wrap.id=`${id}-wrap`;
    const header = document.createElement('div'); header.className='plot-header'; header.innerHTML = `<strong>Channel ${ch}</strong>`;
    const canvas = document.createElement('canvas'); canvas.className='plot-canvas'; canvas.id=id; const cw = Math.max(800, container.clientWidth-40 || 1200); canvas.width = cw; canvas.height=250;
    wrap.appendChild(header); wrap.appendChild(canvas); container.appendChild(wrap);

    // Selection handlers
    let isDown=false; let x0=null; let x1=null; const pad=40;
    const onDown = (e)=>{ isDown=true; const rect=canvas.getBoundingClientRect(); x0 = e.clientX - rect.left; };
    const onMove = (e)=>{ if(!isDown) return; const rect=canvas.getBoundingClientRect(); x1 = e.clientX - rect.left; };
    const onUp   = (e)=>{
      if(!isDown) return; isDown=false; const rect=canvas.getBoundingClientRect(); x1 = e.clientX - rect.left;
      if(x0==null || x1==null) return; const ev = state.currentEvent; if(!ev) return;
      const scale = state.scales[`${ev}:ch${ch}`]; if(!scale) return; const {tmin,tmax}=scale; const w=canvas.width;
      const toT = (xp)=> tmin + clamp((xp - pad)/(w-2*pad),0,1) * (tmax - tmin);
      const ta = toT(Math.min(x0,x1)); const tb = toT(Math.max(x0,x1));
      if(Math.abs(tb-ta) <= 0){ return; }
      state.fitRegion[`${ev}:ch${ch}`] = [ta,tb];
      plotWaveforms();
    };
    canvas.addEventListener('mousedown', onDown);
    canvas.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
  }
  return document.getElementById(id);
}

function getPlotColors(){
  const dark = state.theme !== 'light';
  return {
    bg: dark? '#000' : '#fff',
    axes: dark? '#666' : '#666',
    base: dark? '#ffffff' : '#000000',
    fit: '#ff3030',
    sg: dark? '#5fd35f' : '#2ca02c',
  };
}

function plotWaveforms(){
  const ev = state.currentEvent; if(!ev) return;
  const obj = state.events[ev]; if(!obj) return;
  const container = byId('plots');
  container.innerHTML='';
  const channels = Object.keys(obj.channels).map(Number).sort((a,b)=>a-b);
  state.scales = state.scales || {};
  for(const ch of channels){
    const cv = createOrGetCanvas(container, ch);
    const g = cv.getContext('2d'); const w=cv.width, h=cv.height; const pad=40;
    // theme-aware background and colors
    const colors = getPlotColors();
    g.clearRect(0,0,w,h); g.fillStyle=colors.bg; g.fillRect(0,0,w,h);
    const data = obj.channels[ch]; const t=data.t, y=data.y; if(!t||t.length===0) continue;
    let tmin=t[0], tmax=t[t.length-1]; let ymin=Infinity, ymax=-Infinity; for(const v of y){ if(v<ymin) ymin=v; if(v>ymax) ymax=v; }
    // include fit overlays/markers in y-range so they are visible
    for(const [k,rec] of Object.entries(state.results)){
      const m = k.match(new RegExp(`^evt_${ev}_([^_]+)_ch${ch}(?:_\\d+)?$`));
      if(!m) continue; const fitName=m[1]; const model=FIT_LIB[fitName];
      const region = rec.fit_region || [t[0], t[t.length-1]]; const t0r=Math.max(region[0], t[0]); const t1r=Math.min(region[1], t[t.length-1]);
      if(!(t1r>t0r)) continue;
      if(fitName==='low_pass_max'){
        if(Number.isFinite(rec.value)) { if(rec.value<ymin) ymin=rec.value; if(rec.value>ymax) ymax=rec.value; }
        continue;
      }
      if(!model || !model.func) continue;
      const N= Math.min(300, Math.max(50, Math.floor((t1r-t0r)/(t[1]-t[0]) / (state.decim||1))));
      const ts = linspace(t0r, t1r, N);
      let yf; try{ yf = model.func(ts, ...rec.params); } catch(e){ yf=null; }
      if(!yf) continue; if(rec.inverted) yf = yf.map(v=>-v);
      for(let i=0;i<yf.length;i++){ const v=yf[i]; if(Number.isFinite(v)){ if(v<ymin) ymin=v; if(v>ymax) ymax=v; } }
    }
    if(ymax===ymin){ ymax+=1; ymin-=1; }
    state.scales[`${ev}:ch${ch}`] = {tmin,tmax,ymin,ymax};
    const tx=(tt)=> pad + (tt - tmin)/(tmax-tmin)*(w-2*pad);
    const ty=(vv)=> h-pad - (vv - ymin)/(ymax-ymin)*(h-2*pad);
    // axes
    g.strokeStyle=colors.axes; g.lineWidth=1; g.beginPath(); g.moveTo(pad,pad); g.lineTo(pad,h-pad); g.lineTo(w-pad,h-pad); g.stroke();
    // waveform
    const dec = Math.max(1, state.decim|0);
    g.strokeStyle = colors.base; g.lineWidth=1.2; g.beginPath(); for(let i=0;i<t.length;i+=dec){ const X=tx(t[i]); const Y=ty(y[i]); if(i===0) g.moveTo(X,Y); else g.lineTo(X,Y);} g.stroke();
    // SG overlay
    const sgKey = `${ev}:ch${ch}`; const sgv = state.sgFiltered[sgKey]; if(sgv){ const {t:tt,y:yy}=sgv; g.strokeStyle=colors.sg; g.beginPath(); for(let i=0;i<tt.length;i+=dec){const X=tx(tt[i]); const Y=ty(yy[i]); if(i===0) g.moveTo(X,Y); else g.lineTo(X,Y);} g.stroke(); }
    // Region overlay
    const reg = state.fitRegion[`${ev}:ch${ch}`];
    if(reg && reg.length===2){ g.strokeStyle='#ff5050'; g.setLineDash([5,4]); g.beginPath(); g.moveTo(tx(reg[0]), pad); g.lineTo(tx(reg[0]), h-pad); g.moveTo(tx(reg[1]), pad); g.lineTo(tx(reg[1]), h-pad); g.stroke(); g.setLineDash([]); }
    // Fit overlays/markers for this event+channel
    for(const [k,rec] of Object.entries(state.results)){
      const m = k.match(new RegExp(`^evt_${ev}_([^_]+)_ch${ch}(?:_\\d+)?$`));
      if(!m) continue; const fitName = m[1]; const model = FIT_LIB[fitName]; const invert=!!rec.inverted;
      const region = rec.fit_region || [t[0], t[t.length-1]];
      if(fitName==='low_pass_max'){
        if(Number.isFinite(rec.t_at) && Number.isFinite(rec.value)){
          const X=tx(rec.t_at); const Y=ty(rec.value); const colors=getPlotColors();
          const color = '#2ca02c';
          g.fillStyle = color; g.beginPath(); g.arc(X, Y, 4, 0, Math.PI*2); g.fill();
        }
        continue;
      }
      if(!model || !model.func) continue;
      // Draw on a dense time grid within the region, independent of waveform sample indices
      const a = Math.max(region[0], t[0]); const b = Math.min(region[1], t[t.length-1]);
      if(!(b>a)) continue;
      const dt = (t.length>1)? (t[1]-t[0]) : (b-a);
      let N = Math.floor((b-a)/Math.max(dt, 1e-12));
      N = Math.min(1000, Math.max(100, N*2));
      const tt = linspace(a, b, N);
      let yy; let err=null; try { yy = model.func(tt, ...rec.params); } catch(e) { yy = null; err=e; }
      if(!yy || yy.length!==tt.length) {
        if(state.debugOverlay){ g.fillStyle='#ffcc00'; g.font='10px sans-serif'; g.fillText(`${fitName}: eval failed${err? ' ('+String(err)+')':''}`, pad+6, pad+12); }
        continue;
      }
      const yplot = invert? yy.map(v=>-v): yy;
      const decFit = Math.max(1, Math.floor(tt.length/800));
      g.strokeStyle=colors.fit; g.lineWidth=2.8; g.beginPath();
      for(let i=0;i<tt.length;i+=decFit){ const X=tx(tt[i]); const Y=ty(yplot[i]); if(i===0) g.moveTo(X,Y); else g.lineTo(X,Y);} g.stroke();
      if(state.debugOverlay){
        const idxs=[0, Math.floor((tt.length-1)/2), tt.length-1]; g.fillStyle='#ffcc00';
        for(const ii of idxs){ const X=tx(tt[ii]); const Y=ty(yplot[ii]); if(Number.isFinite(X)&&Number.isFinite(Y)){ g.beginPath(); g.arc(X,Y,2.5,0,Math.PI*2); g.fill(); } }
        let yminF=Infinity,ymaxF=-Infinity; for(const v of yplot){ if(Number.isFinite(v)){ if(v<yminF) yminF=v; if(v>ymaxF) ymaxF=v; } }
        g.fillText(`${fitName}: min=${yminF.toExponential(2)} max=${ymaxF.toExponential(2)}`, pad+6, pad+24);
      }
    }
}
}

// Savitzky–Golay filter (odd window, poly=2)
function savgol(y, window, poly){
  window = (window|0); if(window<5) window=5; if(window%2===0) window+=1;
  const n = y.length; const half = (window-1)/2; const out = new Array(n).fill(0);
  // precompute convolution weights for poly=2 symmetric
  // Using least squares for a quadratic; for simplicity, reuse central weights from known formula
  // For robustness across window sizes we compute weights by solving normal equations
  function weights(m){
    // positions -m..m
    const X = [];
    for(let i=-m;i<=m;i++){ X.push([1,i,i*i]); }
    // normal matrix A = X^T X; target basis e0=[1,0,0] to estimate y at 0
    const A = [[0,0,0],[0,0,0],[0,0,0]];
    for(const row of X){ for(let i=0;i<3;i++) for(let j=0;j<3;j++) A[i][j]+=row[i]*row[j]; }
    // Solve A * c = e0 for c (3x3)
    function inv3(M){
      const [a,b,c,d,e,f,g,h,k]=[M[0][0],M[0][1],M[0][2],M[1][0],M[1][1],M[1][2],M[2][0],M[2][1],M[2][2]];
      const A_ = e*k-f*h, B_ = -(d*k-f*g), C_ = d*h-e*g;
      const D_ = -(b*k-c*h), E_ = a*k-c*g, F_ = -(a*h-b*g);
      const G_ = b*f-c*e, H_ = -(a*f-b*e), K_ = a*e-b*d;
      const det = a*A_ + b*B_ + c*C_;
      return [[A_/det,D_/det,G_/det],[B_/det,E_/det,H_/det],[C_/det,F_/det,K_/det]];
    }
    const invA = inv3(A);
    const e0 = [1,0,0];
    const c = [invA[0][0]*e0[0]+invA[0][1]*e0[1]+invA[0][2]*e0[2], invA[1][0]*e0[0]+invA[1][1]*e0[1]+invA[1][2]*e0[2], invA[2][0]*e0[0]+invA[2][1]*e0[1]+invA[2][2]*e0[2]];
    // convolution weights w_i are dot(c,[1,i,i^2]) for i=-m..m
    const W=[]; for(let i=-m;i<=m;i++){ W.push(c[0]+c[1]*i+c[2]*i*i); }
    return W;
  }
  const W = weights(half);
  for(let i=0;i<n;i++){
    let acc=0, wsum=0;
    for(let k=-half;k<=half;k++){
      let idx = i+k; if(idx<0) idx=0; if(idx>=n) idx=n-1;
      const w=W[k+half]; acc += w*y[idx]; wsum+=w;
    }
    out[i]=acc/wsum;
  }
  return out;
}

// Fit function library (ported)
const FIT_LIB = {
  linear: {
    names: ['m','b'],
    func: (T, m, b) => T.map(t => m*t + b)
  },
  QD3Fit: {
    names: ['t0','q','v','C'],
    func: (T, t0,q,v,C) => {
      v = v*100; const PUT_l=20.0, PUT_d=0.90; const tau1=0.13e-3, tau2=0.15e-3;
      const dt1=PUT_d/v, dt2=(PUT_l-2*PUT_d)/v; const signal=new Array(T.length).fill(0);
      for(let i=0;i<T.length;i++){
        const t=T[i];
        if(t>t0 && t<t0+dt1){ signal[i]=(q/dt1)*(t-t0); }
        else if(t>=t0+dt1 && t<t0+dt1+dt2){ signal[i]= q*Math.exp(-(t-(t0+dt1))/tau1); }
        else if(t>=t0+dt1+dt2 && t<t0+2*dt1+dt2){
          // need qq at end of idx2; approximate with boundary value
          const qq = q*Math.exp(-(dt2)/tau1);
          signal[i] = qq - (q/dt1)*(t-(t0+dt1+dt2));
        } else if(t>=t0+2*dt1+dt2){
          const qq2 = q*Math.exp(-(dt2)/tau1) - (q/dt1)*( (2*dt1+dt2) - (dt1+dt2) );
          signal[i] = qq2 * Math.exp(-(t-(t0+2*dt1+dt2))/tau2);
        }
      }
      return signal.map(v=>v+C);
    }
  },
  QDMFit: {
    names: ['t0','q','v','C'],
    func: (T, t0,q,v,C) => {
      v=v*100; const PUT_l=15.3, PUT_d=2.0; const tau1=0.4e-3, tau2=0.4e-3;
      const dt1=PUT_d/v, dt2=(PUT_l-2*PUT_d)/v; const signal=new Array(T.length).fill(0);
      for(let i=0;i<T.length;i++){
        const t=T[i];
        if(t>t0 && t<t0+dt1){ signal[i]=(q/dt1)*(t-t0); }
        else if(t>t0+dt1 && t<t0+dt1+dt2){ signal[i]= q*Math.exp(-(t-(t0+dt1))/tau1); }
        else if(t>t0+dt1+dt2 && t<t0+2*dt1+dt2){
          const qq = q*Math.exp(-(dt2)/tau1);
          signal[i] = qq - (q/dt1)*(t-(t0+dt1+dt2));
        } else if(t>t0+2*dt1+dt2){
          const qq2 = q*Math.exp(-(dt2)/tau1) - (q/dt1)*( (2*dt1+dt2) - (dt1+dt2) );
          signal[i] = qq2 * Math.exp(-(t-(t0+2*dt1+dt2))/tau2);
        }
      }
      return signal.map(v=>v+C);
    }
  },
  CSA_pulse: {
    names: ['t0','C0','C1','C2','T0','T1','T2','C'],
    func: (x,t0,C0,C1,C2,T0,T1,T2,C) => {
      // Ensure positive time constants to avoid NaN/Inf
      const e = 1e-12;
      const T0p = Math.max(Math.abs(T0), e);
      const T1p = Math.max(Math.abs(T1), e);
      const T2p = Math.max(Math.abs(T2), e);
      const out = new Array(x.length);
      for(let i=0;i<x.length;i++){
        const xx = x[i];
        const left  = (xx<=t0?1:0);
        const right = (xx>=t0?1:0);
        const gauss = Math.exp(-((xx-t0)*(xx-t0))/(T0p*T0p));
        const rise  = (1.0-Math.exp(-(xx-t0)/T1p));
        const decay = Math.exp(-(xx-t0)/T2p);
        const val = (C0 - left*C1*gauss + right*( C2*rise*decay - C1 )) + C;
        out[i] = Number.isFinite(val) ? val : (C0 + C);
      }
      return out;
    }
  },
  skew_gaussian: {
    names: ['A','xi','omega','alpha','C'],
    func: (x,A,xi,omega,alpha,C) => {
      const out=new Array(x.length);
      const invsqrt2pi=1/Math.sqrt(2*Math.PI);
      const op = Math.max(Math.abs(omega), 1e-12);
      for(let i=0;i<x.length;i++){
        const t=(x[i]-xi)/op; const phi=invsqrt2pi*Math.exp(-0.5*t*t); const Phi=0.5*(1+erf(alpha*t/Math.sqrt(2)));
        const val = A*2/op*phi*Phi + C; out[i] = Number.isFinite(val)? val: C;
      }
      return out;
    }
  },
  gaussian: {
    names: ['A','mu','sigma','C'],
    func: (x,A,mu,sigma,C) => x.map(xx => A*Math.exp(-0.5*((xx-mu)/sigma)**2)+C)
  }
  ,
  // Alternative CSA pulse with different parameter naming and numerically-stable eval
  CSA_pulse_alt: {
    // t0: transition time; b0: baseline; a_neg: negative Gaussian lobe amplitude; a_pos: positive lobe amplitude
    // tg: Gaussian width; tr: rise time; td: decay time; off: constant offset
    names: ['t0','b0','a_neg','a_pos','tg','tr','td','off'],
    func: (x, t0, b0, a_neg, a_pos, tg, tr, td, off) => {
      const e = 1e-12;
      const tgP = Math.max(Math.abs(tg), e);
      const trP = Math.max(Math.abs(tr), e);
      const tdP = Math.max(Math.abs(td), e);
      const out = new Array(x.length);
      for(let i=0;i<x.length;i++){
        const dt = x[i] - t0;
        const left  = dt <= 0 ? 1 : 0;
        const right = dt >= 0 ? 1 : 0;
        const gauss = Math.exp(-(dt*dt)/(tgP*tgP));
        const rise  = (1.0 - Math.exp(-dt/trP));
        const decay = Math.exp(-dt/tdP);
        // Keep structure similar to original but with explicit variables
        let val = b0 - left * a_neg * gauss + right * (a_pos * rise * decay - a_neg) + off;
        if(!Number.isFinite(val)) val = b0 + off;
        out[i] = val;
      }
      return out;
    }
  },
  // Alternative Skew-Gaussian with different naming and explicit safeguards
  skew_gaussian_alt: {
    // amp: amplitude; loc: center; scale: width (>0); shape: skew; off: offset
    names: ['amp','loc','scale','shape','off'],
    func: (x, amp, loc, scale, shape, off) => {
      const out = new Array(x.length);
      const invsqrt2pi = 1/Math.sqrt(2*Math.PI);
      const s = Math.max(Math.abs(scale), 1e-12);
      for(let i=0;i<x.length;i++){
        const t = (x[i] - loc) / s;
        const phi = invsqrt2pi * Math.exp(-0.5*t*t);
        const Phi = 0.5 * (1 + erf(shape * t / Math.sqrt(2)));
        let v = amp * 2/s * phi * Phi + off;
        if(!Number.isFinite(v)) v = off;
        out[i] = v;
      }
      return out;
    }
  },
  // Measurement (no curve): SG low-pass then take max (invert => min)
  low_pass_max: {
    names: ['width'],
  },
};

// erf implementation
function erf(x){
  // Abramowitz-Stegun approximation
  const sign = x<0?-1:1; x=Math.abs(x);
  const a1=0.254829592,a2=-0.284496736,a3=1.421413741,a4=-1.453152027,a5=1.061405429,p=0.3275911;
  const t=1/(1+p*x); const y=1-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*Math.exp(-x*x);
  return sign*y;
}

// Derived metrics for QD3/QDM
function qd3Cal(q, v, rho=7500){
  const c = 1e-12 / 20 * 0.95;
  const charge = c * q;
  const mass = Math.abs(2 * 2.2e6 * charge / (v*v));
  const radius = Math.pow((3.0 * mass) / (4.0 * Math.PI * rho), 1/3);
  return {charge, mass, radius};
}
function qdmCal(q, v, rho=7500){
  const c = 1e-12 / 20 * 0.95 * 2.35;
  const charge = c * q;
  const mass = Math.abs(2 * 2.2e6 * charge / (v*v));
  const radius = Math.pow((3.0 * mass) / (4.0 * Math.PI * rho), 1/3);
  return {charge, mass, radius};
}
function attachDerivedMetrics(fitName, rec){
  try{
    if(!rec || !Array.isArray(rec.params)) return;
    if(fitName==='QD3Fit'){
      const t0=rec.params[0], q=rec.params[1], v=rec.params[2];
      if(Number.isFinite(q) && Number.isFinite(v) && v!==0){ Object.assign(rec, qd3Cal(q, v)); }
    } else if(fitName==='QDMFit'){
      const t0=rec.params[0], q=rec.params[1], v=rec.params[2];
      if(Number.isFinite(q) && Number.isFinite(v) && v!==0){ Object.assign(rec, qdmCal(q, v)); }
    }
  } catch(_){}
}

// Simple Nelder-Mead for least squares
function nelderMead(f, x0, opts={}){
  const maxIter = opts.maxIter||800;
  const ftol = opts.ftol||1e-8;
  const alpha=1, gamma=2, rho=0.5, sigma=0.5; // standard params
  const n=x0.length;
  let simplex = new Array(n+1);
  simplex[0]=x0.slice();
  const step = opts.step||1e-2;
  for(let i=0;i<n;i++){
    const xi=x0.slice(); xi[i]=xi[i]!==0? xi[i]*(1+step): step; simplex[i+1]=xi;
  }
  function sort(){ simplex.sort((a,b)=>f(a)-f(b)); }
  sort();
  let it=0;
  let prevBest = f(simplex[0]);
  while(it++<maxIter){
    sort();
    const best=simplex[0], worst=simplex[n], secondWorst=simplex[n-1];
    const fbest = f(best), fworst = f(worst);
    if(Math.abs(fworst - fbest) < ftol*(1+Math.abs(fbest))){ break; }
    // centroid of all but worst
    const centroid=new Array(n).fill(0);
    for(let i=0;i<n;i++) for(let j=0;j<n;j++) centroid[j]+=simplex[i][j];
    for(let j=0;j<n;j++) centroid[j]/=n;
    // reflection
    const xr=centroid.map((c,j)=>c+alpha*(c-worst[j]));
    if(f(xr)<f(secondWorst)){ simplex[n]=xr; continue; }
    if(f(xr)<f(best)){ // expansion
      const xe=centroid.map((c,j)=>c+gamma*(xr[j]-c));
      simplex[n]= f(xe)<f(xr)? xe: xr; continue;
    }
    // contraction
    const xc=centroid.map((c,j)=>c+rho*(worst[j]-c));
    if(f(xc)<f(worst)){ simplex[n]=xc; continue; }
    // shrink
    for(let i=1;i<n+1;i++) simplex[i]=simplex[0].map((b,j)=> b + sigma*(simplex[i][j]-b));
  }
  sort();
  return {x: simplex[0], f: f(simplex[0]), iters: it};
}

function lsqFit(model, t, y, p0){
  const f = (p)=>{
    try{
      const yy=model.func(t, ...p);
      if(!yy || yy.length!==y.length) return 1e24;
      let s=0; for(let i=0;i<y.length;i++){ const d=yy[i]-y[i]; if(!Number.isFinite(d)) return 1e24; s+=d*d; }
      return Number.isFinite(s)? s: 1e24;
    }catch(e){
      return 1e24;
    }
  };
  // faster: limit iterations and allow early stop
  const res = nelderMead(f, p0, {maxIter: 600, step: 1e-2, ftol: 1e-8});
  return res.x;
}

// Subsample time series for faster fitting
function subsample(t, y, maxPts=2000){
  const n=t.length; if(n<=maxPts) return {t,y};
  const step = Math.ceil(n/maxPts);
  const tt=new Array(Math.ceil(n/step)); const yy=new Array(tt.length);
  let j=0; for(let i=0;i<n;i+=step){ tt[j]=t[i]; yy[j]=y[i]; j++; }
  return {t:tt, y:yy};
}

// Guess parameters (basic heuristics)
function guessParams(fitName, names, t, y){
  const g={};
  const amp = arrMax(y)-arrMin(y);
  // Improve t0 and v for QD3/QDM by simple derivative heuristics
  if(fitName==='QD3Fit' || fitName==='QDMFit'){
    let rise=0, fall=0; if(t.length>2){
      // use a light smoothing length to reduce noise in derivative
      const n=t.length; const half=Math.min(51, Math.max(5, Math.floor(n/50))); const win=half%2? half: half+1;
      const ys = savgol(y, win, 2);
      const diff=new Array(n-1); for(let i=0;i<n-1;i++) diff[i]=ys[i+1]-ys[i];
      rise = argMax(diff); let minD=Infinity, idxMin=-1; for(let i=0;i<diff.length;i++){ if(diff[i]<minD){ minD=diff[i]; idxMin=i; } } fall=idxMin;
      const dt = Math.max(1e-12, Math.abs(t[rise]-t[fall]));
      const det = (fitName==='QD3Fit')? 0.19: 0.133; // detector constants from desktop heuristics
      g.v = det / dt;
      g.t0 = t[Math.max(0, Math.min(fall, t.length-1))];
    }
  }
  if(fitName==='CSA_pulse'){
    const n=t.length; if(n>3){
      // Use smoothed derivative to locate transition near t0
      const half=Math.min(51, Math.max(5, Math.floor(n/50))); const win=half%2? half: half+1;
      const ys = savgol(y, win, 2);
      const diff=new Array(n-1); for(let i=0;i<n-1;i++) diff[i]=ys[i+1]-ys[i];
      const idx = argMax(diff.map(v=>Math.abs(v)));
      g.t0 = t[Math.max(0, Math.min(idx, n-1))];
    }
    const w = Math.max(1e-9, t[t.length-1]-t[0]);
    const preIdx = Math.max(0, Math.floor(0.1*t.length));
    const base = mean(y.slice(0, preIdx||1));
    g.C  = mean(y);
    g.C0 = base;
    const A = amp || (arrAbsMax(y) || 1e-3);
    g.C1 = 0.5*A;
    g.C2 = A;
    g.T0 = w/20;
    g.T1 = w/10;
    g.T2 = w/5;
  }
  if(fitName==='CSA_pulse_alt'){
    const n=t.length; if(n>3){
      const half=Math.min(51, Math.max(5, Math.floor(n/50))); const win=half%2? half: half+1;
      const ys = savgol(y, win, 2);
      const diff=new Array(n-1); for(let i=0;i<n-1;i++) diff[i]=ys[i+1]-ys[i];
      const idx = argMax(diff.map(v=>Math.abs(v)));
      g.t0 = t[Math.max(0, Math.min(idx, n-1))];
    }
    const w = Math.max(1e-9, t[t.length-1]-t[0]);
    const preIdx = Math.max(1, Math.floor(0.1*t.length));
    const base = mean(y.slice(0, preIdx));
    g.off = 0.0;     // keep offset near 0 since we detrend
    g.b0 = base;
    const A = amp || (arrAbsMax(y) || 1e-3);
    g.a_neg = 0.5*A;
    g.a_pos = A;
    g.tg = w/20; g.tr = w/10; g.td = w/5;
  }
  if(names.includes('t0') && g.t0===undefined) g.t0 = t[0];
  if(names.includes('q')) g.q = amp;
  if(names.includes('v') && g.v===undefined) g.v = 1.0;
  if(names.includes('C')) g.C = mean(y);
  if(fitName==='skew_gaussian_alt'){
    const i=argMax(y);
    g.amp = amp || (arrAbsMax(y) || 1e-3);
    g.loc = t[i];
    g.scale = Math.max( (t[t.length-1]-t[0])/6, 1e-6 );
    g.shape = 0.0;
    g.off = 0.0; // detrended
  }
  if(fitName==='gaussian'){
    const i=argMax(y);
    g.A = amp; g.mu=t[i]; g.sigma = Math.max(1e-9, (t[t.length-1]-t[0])/10); g.C = mean(y);
  }
  if(fitName==='skew_gaussian'){
    const area=y.reduce((a,v,i)=> a + (i? (y[i-1]+v)*(t[i]-t[i-1])*0.5:0), 0);
    const i=argMax(y);
    g.A=area; g.xi=t[i]; g.omega=(t[t.length-1]-t[0])/6; g.alpha=0; g.C=mean(y);
  }
  if(fitName==='low_pass_max'){
    const w = parseInt(byId('sg-width').value, 10) || 200;
    return [w % 2 === 0 ? w+1 : w];
  }
  return names.map(n => (n in g)? g[n]: 0);
}

// File loading
byId('file-input').addEventListener('change', async (e)=>{
  state.files = Array.from(e.target.files||[]).filter(f=>/\.(csv|txt|trc)$/i.test(f.name));
  state.eventsIndex = {}; state.events={}; state.eventOrder=[]; state.currentEvent=null; state.results={}; state.sgFiltered={}; state.featureFilter=null;
  // derive dataset id from filenames for scoping
  state.datasetId = computeDatasetId(state.files.map(f=> f.webkitRelativePath || f.name).sort());
  // Index only (do not parse now)
  for(const file of state.files){
    const meta = parseFileName(file.name);
    if(!meta) continue; const key = meta.event; const ch = meta.ch;
    if(!state.eventsIndex[key]) state.eventsIndex[key] = {channels:{}};
    state.eventsIndex[key].channels[ch] = file;
  }
  state.eventOrder = Object.keys(state.eventsIndex).sort((a,b)=>parseInt(a)-parseInt(b));
  if(state.eventOrder.length>0){
    state.currentEvent=state.eventOrder[0];
    await loadEvent(state.currentEvent);
  }
  refreshEventSelect(); refreshChannelCombos(); plotWaveforms();
});

// Busy overlay helpers
function showBusy(text){ const el=byId('busy'); const t=byId('busy-text'); if(t) t.textContent = text || 'Working…'; if(el) el.classList.remove('hidden'); }
function updateBusy(text){ const t=byId('busy-text'); if(t) t.textContent = text; }
function hideBusy(){ const el=byId('busy'); if(el) el.classList.add('hidden'); }

function computeDatasetId(names){
  // FNV-1a 32-bit hash of concatenated names
  let h = 0x811c9dc5;
  for(const name of names){
    for(let i=0;i<name.length;i++){
      h ^= name.charCodeAt(i);
      h = (h + ((h<<1) + (h<<4) + (h<<7) + (h<<8) + (h<<24))) >>> 0;
    }
  }
  return 'ds_' + h.toString(16).padStart(8,'0');
}

async function loadEvent(ev){
  if(!ev || !state.eventsIndex[ev]) return;
  // If already loaded, skip
  if(state.events[ev] && state.events[ev].channels) return;
  const channels = state.eventsIndex[ev].channels;
  const obj = {channels:{}};
  for(const ch of Object.keys(channels)){
    const file = channels[ch]; let t=[], y=[];
    if(/\.(trc)$/i.test(file.name)){
      try{ const out = await TrcReader.parseFile(file); t=out.t; y=out.y; }
      catch(err){ console.warn('TRC parse failed', file.name, err); continue; }
    } else {
      const out = await parseCSV(file); t=out.t; y=out.y;
    }
    const ym = mean(y); const yc = y.map(v=>v-ym);
    obj.channels[parseInt(ch,10)] = {t, y: yc};
  }
  state.events[ev] = obj;
}

async function getEventChannelData(ev, ch){
  if(!state.events[ev]) await loadEvent(ev);
  const evt = state.events[ev]; if(!evt) return null;
  return evt.channels[ch] || null;
}

function refreshEventSelect(){
  const sel = byId('event-select'); sel.innerHTML='';
  const keys = state.featureFilter? state.featureFilter: state.eventOrder;
  for(const k of keys){ const opt=document.createElement('option'); opt.value=k; opt.textContent=k; sel.appendChild(opt);}  
  sel.value = state.currentEvent||'';
}

function refreshChannelCombos(){
  const chs = new Set();
  const ev = state.currentEvent; if(!ev) return;
  const channelsIdx = state.eventsIndex[ev]? state.eventsIndex[ev].channels: {};
  Object.keys(channelsIdx).forEach(k=>chs.add(parseInt(k)));
  const arr = Array.from(chs).sort((a,b)=>a-b);
  const sgc = byId('sg-chan'); sgc.innerHTML=''; arr.forEach(ch=>{ const o=document.createElement('option'); o.value=ch; o.textContent=String(ch); sgc.appendChild(o);});
  const fc = byId('fit-chan'); fc.innerHTML=''; arr.forEach(ch=>{ const o=document.createElement('option'); o.value=ch; o.textContent=String(ch); fc.appendChild(o);});
}

// Event selection
byId('event-select').addEventListener('change', async (e)=>{ state.currentEvent=e.target.value; await loadEvent(state.currentEvent); refreshChannelCombos(); plotWaveforms(); });
byId('prev-btn').addEventListener('click', async ()=>{ const arr=(state.featureFilter? state.featureFilter: state.eventOrder); const idx=arr.indexOf(state.currentEvent); if(idx>0){ state.currentEvent=arr[idx-1]; byId('event-select').value=state.currentEvent; await loadEvent(state.currentEvent); refreshChannelCombos(); plotWaveforms(); }});
byId('next-btn').addEventListener('click', async ()=>{ const arr=(state.featureFilter? state.featureFilter: state.eventOrder); const idx=arr.indexOf(state.currentEvent); if(idx<arr.length-1){ state.currentEvent=arr[idx+1]; byId('event-select').value=state.currentEvent; await loadEvent(state.currentEvent); refreshChannelCombos(); plotWaveforms(); }});

// Decimation
byId('decim').addEventListener('input', (e)=>{ state.decim=clamp(parseInt(e.target.value,10)||1,1,100); plotWaveforms(); });
byId('fit-res').addEventListener('input', (e)=>{ const v=parseInt(e.target.value,10); if(Number.isFinite(v)) state.fitResolution = clamp(v, 200, 10000); });

// Fit select
const fitSel = byId('fit-select'); Object.keys(FIT_LIB).forEach(k=>{ const o=document.createElement('option'); o.value=k; o.textContent=k; fitSel.appendChild(o); });

// SG filter
byId('sg-btn').addEventListener('click', async ()=>{
  const ev=state.currentEvent; if(!ev) return; const ch=parseInt(byId('sg-chan').value,10);
  const data = await getEventChannelData(ev, ch); if(!data) return; const {t,y}=data;
  const w=parseInt(byId('sg-width').value,10)||200;
  showBusy('Filtering…');
  try{
    const res = await compute.sg(y, w, 2);
    const yy = res && res.y ? res.y : savgol(y, w, 2);
    state.sgFiltered[`${ev}:ch${ch}`]={t:[...t],y:yy};
  } finally {
    hideBusy();
  }
  plotWaveforms();
});
byId('sg-clear-btn').addEventListener('click', ()=>{ const ev=state.currentEvent; if(!ev) return; const ch=parseInt(byId('sg-chan').value,10); delete state.sgFiltered[`${ev}:ch${ch}`]; plotWaveforms(); });

// Dynamic fit
byId('run-fit-btn').addEventListener('click', async ()=>{
  const ev=state.currentEvent; if(!ev) return; const fname=byId('fit-select').value; const model=FIT_LIB[fname];
  const ch=parseInt(byId('fit-chan').value,10); const invert=byId('fit-invert').checked;
  const dta = await getEventChannelData(ev, ch); if(!dta) return; let {t,y}=dta; if(fname!=='low_pass_max' && invert) y = y.map(v=>-v);
  // selection window: use saved region if present
  const reg = state.fitRegion[`${ev}:ch${ch}`];
  if(reg && reg.length===2){ const mask=t.map((tt,i)=> tt>=reg[0] && tt<=reg[1]? i: -1).filter(i=>i>=0); if(mask.length>5){ t=mask.map(i=>t[i]); y=mask.map(i=>y[i]); } }
  // subsample for fitting speed
  const ss = subsample(t, y, state.fitResolution||2000); const tf=ss.t, yf=ss.y;
  const keyBase = `evt_${ev}_${fname}_ch${ch}`; let idx=1; let key=keyBase; while(state.results[key]){ idx+=1; key=`${keyBase}_${idx}`; }
  const fitRegion = (reg && reg.length===2)? [reg[0], reg[1]]: [t[0], t[t.length-1]];
  showBusy('Fitting…');
  try{
    if(fname==='low_pass_max'){
      const width = (guessParams(fname, (model && model.names)||['width'], t, y)[0])|0;
      const res = await compute.lowPassMax(tf, yf, width, invert);
      state.results[key] = { params: [res.width], param_names: ['width'], fit_region: fitRegion, inverted: invert, value: res.value, t_at: res.t_at, dataset_id: state.datasetId };
    } else {
      const p0 = guessParams(fname, model.names, t, y);
      const res = await compute.fit(fname, tf, yf, p0);
      const popt = res && res.params ? res.params : lsqFit(model, tf, yf, p0);
      const rec = { params: popt, param_names: model.names, fit_region: fitRegion, inverted: invert, dataset_id: state.datasetId };
      attachDerivedMetrics(fname, rec);
      state.results[key] = rec;
    }
  } catch(e){ console.error('Fit failed', e); alert('Fit failed: '+e.message); hideBusy(); return; }
  // replot; overlay handled in plotWaveforms
  plotWaveforms();
  hideBusy();
});

// Clear fit(s) for current event/fit/channel
byId('clear-fit-btn').addEventListener('click', ()=>{
  const ev=state.currentEvent; if(!ev) return; const fname=byId('fit-select').value; const ch=parseInt(byId('fit-chan').value,10);
  const pref=`evt_${ev}_${fname}_ch${ch}`; let removed=0; for(const k of Object.keys(state.results)){ if(k===pref || k.startsWith(pref+'_')){ delete state.results[k]; removed++; }} alert(`Removed ${removed} fit(s).`); plotWaveforms();
});

// Adjust Fit dialog
byId('adjust-btn').addEventListener('click', ()=>{
  const ev=state.currentEvent; if(!ev) return; const fname=byId('fit-select').value; const ch=parseInt(byId('fit-chan').value,10);
  const keys = Object.keys(state.results).filter(k=> new RegExp(`^evt_${ev}_${fname}_ch${ch}(?:_\\d+)?$`).test(k));
  if(keys.length===0){ alert('No fits found to adjust for this event/channel/fit.'); return; }
  openAdjustDialog(ev, ch, fname, keys);
});

function openAdjustDialog(ev, ch, fitName, keys){
  const modal = byId('adjust'); const sel=byId('adjust-fit-select'); const form=byId('adjust-params');
  sel.innerHTML=''; keys.forEach(k=>{ const o=document.createElement('option'); o.value=k; o.textContent=k; sel.appendChild(o); });
  const snapshot = new Map(); keys.forEach(k=> snapshot.set(k, { params: [...state.results[k].params] }));
  function buildFormForKey(key){
    form.innerHTML=''; const rec=state.results[key]; const names=rec.param_names; const vals=rec.params;
    for(let i=0;i<names.length;i++){
      const lab=document.createElement('label'); lab.textContent=names[i];
      const inp=document.createElement('input'); inp.type='number'; inp.step='any'; inp.value = String(vals[i]); inp.dataset.pname=names[i];
      const wrap=document.createElement('div'); wrap.className='modal-row'; wrap.appendChild(lab); wrap.appendChild(inp); form.appendChild(wrap);
    }
  }
  function readParams(){ const rec=state.results[sel.value]; const names=rec.param_names; const inputs=form.querySelectorAll('input'); const arr=[]; inputs.forEach(inp=>{ const v=parseFloat(inp.value); arr.push(Number.isFinite(v)? v: 0); }); return arr; }
  function plotPreview(){ const key=sel.value; const rec=state.results[key]; const model=FIT_LIB[fitName]; if(!model) return; const newParams=readParams(); rec.params=[...newParams]; plotWaveforms(); }
  async function doRefit(){
    const key=sel.value; const rec=state.results[key]; const model=FIT_LIB[fitName]; if(!model) return;
    const data = await getEventChannelData(ev, ch); if(!data) return; let {t,y}=data; if(rec.inverted) y=y.map(v=>-v);
    const region = rec.fit_region; if(region && region.length===2){ const mask=t.map((tt,i)=> tt>=region[0] && tt<=region[1]? i: -1).filter(i=>i>=0); if(mask.length>5){ t=mask.map(i=>t[i]); y=mask.map(i=>y[i]); } }
    const ss=subsample(t,y,state.fitResolution||2000); const p0=readParams(); const popt=lsqFit(model, ss.t, ss.y, p0);
    rec.params=[...popt]; // update
    attachDerivedMetrics(fitName, rec);
    buildFormForKey(key); // refresh inputs to reflect refit
    plotWaveforms();
  }
  // init
  buildFormForKey(keys[0]);
  modal.classList.remove('hidden');
  sel.onchange = ()=> buildFormForKey(sel.value);
  byId('adjust-close').onclick = ()=>{ cancelAdjust(); };
  byId('adjust-cancel').onclick = ()=>{ cancelAdjust(); };
  byId('adjust-plot').onclick = ()=>{ plotPreview(); };
  byId('adjust-refit').onclick = async ()=>{ showBusy('Refitting…'); try{ await doRefit(); } finally { hideBusy(); } };
  byId('adjust-apply').onclick = ()=>{ // persist current form values and recompute derived fields
    const key=sel.value; const rec=state.results[key]; rec.params = readParams(); attachDerivedMetrics(fitName, rec); modal.classList.add('hidden'); plotWaveforms(); cleanup(); };
  function cancelAdjust(){ // restore snapshot
    snapshot.forEach((v,k)=>{ if(state.results[k]) state.results[k].params = [...v.params]; });
    modal.classList.add('hidden'); plotWaveforms(); cleanup();
  }
  function cleanup(){ sel.onchange=null; byId('adjust-close').onclick=null; byId('adjust-cancel').onclick=null; byId('adjust-plot').onclick=null; byId('adjust-refit').onclick=null; byId('adjust-apply').onclick=null; }
}

// Batch fit (basic): fit across a range using current fit and channel
byId('batch-fit-btn').addEventListener('click', async ()=>{
  const fname=byId('fit-select').value; const model=FIT_LIB[fname]; const ch=parseInt(byId('fit-chan').value,10); const invert=byId('fit-invert').checked;
  const evKeys = state.featureFilter? state.featureFilter: state.eventOrder;
  if(evKeys.length===0) return;
  // prompt simple range and window
  const start = prompt('Start event (key)', evKeys[0]); if(!start) return;
  const end = prompt('End event (key)', evKeys[evKeys.length-1]); if(!end) return;
  const t0s = prompt('Time start (s)', ''); const t1s = prompt('Time end (s)', '');
  const t0 = parseFloat(t0s), t1=parseFloat(t1s); const added=[];
  const keys = evKeys.slice(evKeys.indexOf(start), evKeys.indexOf(end)+1);
  showBusy('Batch fitting…');
  for(const ev of keys){ const data = await getEventChannelData(ev, ch); if(!data) continue; let {t,y}=data; let t0w=t0, t1w=t1; // prefer per-channel saved region
    const reg = state.fitRegion[`${ev}:ch${ch}`]; if(reg && reg.length===2){ t0w = reg[0]; t1w = reg[1]; }
    if(isFinite(t0w) && isFinite(t1w)){ const mask=t.map((tt,i)=> tt>=t0w && tt<=t1w? i: -1).filter(i=>i>=0); if(mask.length>5){ t=mask.map(i=>t[i]); y=mask.map(i=>y[i]); } }
    const ss=subsample(t,y,state.fitResolution||2000);
    try{
      const base=`evt_${ev}_${fname}_ch${ch}`; let i=1; let k=base; while(state.results[k]){ i++; k=`${base}_${i}`; }
      if(fname==='low_pass_max'){
        const width = (guessParams(fname, (model && model.names)||['width'], t, y)[0])|0;
        const res = await compute.lowPassMax(ss.t, ss.y, width, invert);
        state.results[k]={ params:[res.width], param_names:['width'], fit_region:[t[0],t[t.length-1]], inverted: invert, value: res.value, t_at: res.t_at, dataset_id: state.datasetId };
      } else {
        const yy = invert? ss.y.map(v=>-v): ss.y;
        const p0=guessParams(fname, model.names, ss.t, yy);
        const fres = await compute.fit(fname, ss.t, yy, p0);
        const popt = fres && fres.params ? fres.params : lsqFit(model, ss.t, yy, p0);
        const rec = {params:popt,param_names:model.names,fit_region:[t[0],t[t.length-1]],inverted:invert, dataset_id: state.datasetId};
        attachDerivedMetrics(fname, rec);
        state.results[k]=rec;
      }
      added.push(k);
    }catch(e){/*skip*/}
    // free memory for non-current events
    if(ev!==state.currentEvent){ delete state.events[ev]; }
    updateBusy(`Batch fitting… ${added.length}/${keys.length}`);
  }
  alert(`Batch complete. Added ${added.length} fits.`);
  plotWaveforms();
  hideBusy();
});

// Feature scan
byId('scan-btn').addEventListener('click', async ()=>{
  const evKeys=state.eventOrder; if(evKeys.length===0) return; const start=prompt('Start event (key)', evKeys[0]); if(!start) return; const end=prompt('End event (key)', evKeys[evKeys.length-1]); if(!end) return; const ch=parseInt(prompt('Channel (1-4)','1'),10); const thr=parseFloat(prompt('Threshold × std','5'));
  const useAbs = confirm('Use absolute value for peak? Ok=Yes, Cancel=No');
  const keys = evKeys.slice(evKeys.indexOf(start), evKeys.indexOf(end)+1);
  const matches=[];
  for(const ev of keys){ const data = await getEventChannelData(ev, ch); if(!data) continue; const {y}=data; const sd=stdev(y); if(sd<=0) { if(ev!==state.currentEvent){ delete state.events[ev]; } continue; } const peak = useAbs? arrAbsMax(y) : arrMax(y); if(peak>thr*sd) matches.push(ev); if(ev!==state.currentEvent){ delete state.events[ev]; } }
  if(matches.length===0){ alert('No matches.'); return; }
  // preview & filter
  const apply = confirm(`Matches (${matches.length}). Apply as navigation filter?`);
  if(apply){ state.featureFilter = matches; refreshEventSelect(); if(!matches.includes(state.currentEvent)){ state.currentEvent=matches[0]; byId('event-select').value=state.currentEvent; } plotWaveforms(); }
});
byId('clear-filter-btn').addEventListener('click', ()=>{ state.featureFilter=null; refreshEventSelect(); plotWaveforms(); });

// Export results
byId('export-btn').addEventListener('click', ()=>{
  const data = JSON.stringify(state.results, null, 2);
  const blob = new Blob([data], {type:'application/json'});
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download='fit_results.json'; a.click(); URL.revokeObjectURL(a.href);
});
byId('clear-results-btn').addEventListener('click', ()=>{ if(confirm('Clear all results?')){ state.results={}; alert('Cleared.'); }});

// Import results (JSON)
byId('import-btn').addEventListener('click', ()=>{ byId('import-json-input').click(); });
byId('import-json-input').addEventListener('change', async (e)=>{
  const f = e.target.files && e.target.files[0]; if(!f) return;
  try{
    const text = await f.text();
    const obj = JSON.parse(text);
    if(typeof obj !== 'object' || Array.isArray(obj)) throw new Error('Invalid JSON structure');
    const onlyDs = true; // default to safe merge with dataset guard
    let added=0, skipped=0, renamed=0;
    for(const [key, rec] of Object.entries(obj)){
      if(!/^evt_/.test(key) || typeof rec!=='object'){ skipped++; continue; }
      // dataset guard: allow import if dataset_id matches or absent; else allow but keep dataset_id
      const dsOk = (!rec.dataset_id) || (rec.dataset_id===state.datasetId);
      if(!dsOk){ /* keep but count */ }
      let newKey = key;
      if(state.results[newKey]){ let i=1; const base=newKey; while(state.results[newKey]){ i++; newKey = base+`_${i}`; } if(newKey!==key) renamed++; }
      state.results[newKey] = rec;
      added++;
    }
    alert(`Imported ${added} result(s). ${renamed? `Renamed ${renamed} duplicate key(s). `:''}${skipped? `Skipped ${skipped} invalid entr${skipped>1?'ies':'y'}.`:''}`);
    plotWaveforms();
  } catch (err){ alert('Import failed: '+String(err&&err.message||err)); }
  e.target.value='';
});

// HDF5 export/import (requires h5wasm vendor bundle). Falls back to CSV if unavailable on export.
async function ensureH5Wasm(){
  if(window.h5wasmReady) return true;
  // Prefer ESM build if available
  const esmPaths = [
    './vendor/h5wasm/dist/esm/hdf5_hl.js',
    'vendor/h5wasm/dist/esm/hdf5_hl.js'
  ];
  for(const p of esmPaths){
    try{
      const mod = await import(p);
      if(mod && mod.ready){ await mod.ready; window.h5wasmMod = mod; window.h5wasmReady = true; return true; }
    }catch(e){ /* try next */ }
  }
  // Fallback to classic script loader for alternative builds
  const scripts = [
    'vendor/h5wasm/hdf5_hl.js',
    'vendor/h5wasm/h5wasm.js',
    'vendor/h5wasm/hdf5_wasm.js'
  ];
  for(const src of scripts){
    try{
      await new Promise((res, rej)=>{ const s=document.createElement('script'); s.src=src; s.onload=()=>res(); s.onerror=()=>rej(new Error('load failed')); document.head.appendChild(s); });
      window.h5wasmMod = window.h5wasm || window.h5wasmDefault || null; // best-effort
      window.h5wasmReady = !!window.h5wasmMod;
      return window.h5wasmReady;
    } catch(e){ /* try next */ }
  }
  return false;
}

function rowsToColumns(rows){
  if(rows.length===0) return {};
  const cols = Array.from(rows.reduce((s,r)=>{ Object.keys(r).forEach(k=>s.add(k)); return s; }, new Set()));
  const out={}; cols.forEach(c=> out[c] = rows.map(r=> r[c] ?? null));
  return out;
}

function columnsToRows(cols){
  const keys = Object.keys(cols||{});
  if(keys.length===0) return [];
  const n = Math.max(...keys.map(k=> (cols[k]||[]).length));
  const rows = new Array(n).fill(null).map(()=>({}));
  for(const k of keys){ const arr = cols[k]||[]; for(let i=0;i<n;i++){ const v = arr[i]; if(v!=null) rows[i][k]=v; } }
  return rows;
}

async function exportHDF5(){
  const ok = await ensureH5Wasm(); if(!ok) throw new Error('h5wasm not available');
  const mod = window.h5wasmMod; const File = mod.File || (mod.h5wasm && mod.h5wasm.File);
  const Module = mod.Module || (mod.h5wasm && mod.h5wasm.Module);
  const FS = mod.FS || (mod.h5wasm && mod.h5wasm.FS) || (Module && Module.FS);
  if(!File || !Module || !FS) throw new Error('h5wasm API not found');
  const rows = buildResultsRows(); if(rows.length===0) throw new Error('No fits to export');
  const filename = `/fit_results_${Date.now()}.h5`;
  const f = new File(filename, 'w');
  try{
    const root = f; const grp = root.create_group('results');
    const cols = rowsToColumns(rows);
    for(const [name, data] of Object.entries(cols)){
      // Normalize types: numbers/strings; booleans to 0/1 for consistency
      const arr = data.map(v=> typeof v==='boolean'? (v?1:0): v);
      grp.create_dataset({ name, data: arr });
    }
    f.flush();
  } finally {
    f.close();
  }
  const bytes = FS.readFile(filename);
  const blob = new Blob([bytes], {type:'application/x-hdf5'});
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download='fit_results.h5'; a.click(); URL.revokeObjectURL(a.href);
}

function mergeRowsIntoResults(rows){
  let added=0, renamed=0;
  for(const row of rows){
    const event = row.event; const ch = row.channel; const fit = row.fit_type; if(event==null || ch==null || !fit) continue;
    const fitName = String(fit);
    const paramNames = (FIT_LIB[fitName] && FIT_LIB[fitName].names) ? FIT_LIB[fitName].names : Object.keys(row).filter(k=> !['event','channel','fit_type','inverted','dataset_id','value','t_at','charge','mass','radius'].includes(k));
    const params = paramNames.map(n=> row[n] ?? 0);
    const keyBase = `evt_${event}_${fitName}_ch${ch}`; let key=keyBase; let i=1; while(state.results[key]){ i++; key = `${keyBase}_${i}`; }
    if(key!==keyBase) renamed++;
    const rec = { params, param_names: paramNames, fit_region: null, inverted: !!row.inverted, dataset_id: row.dataset_id || state.datasetId };
    if(row.value!=null) rec.value = row.value; if(row.t_at!=null) rec.t_at=row.t_at; if(row.charge!=null) rec.charge=row.charge; if(row.mass!=null) rec.mass=row.mass; if(row.radius!=null) rec.radius=row.radius;
    // Compute derived metrics if missing
    attachDerivedMetrics(fitName, rec);
    state.results[key] = rec; added++;
  }
  return {added, renamed};
}

async function importHDF5File(file){
  const ok = await ensureH5Wasm(); if(!ok) throw new Error('h5wasm not available');
  const mod = window.h5wasmMod; const File = mod.File || (mod.h5wasm && mod.h5wasm.File);
  const Module = mod.Module || (mod.h5wasm && mod.h5wasm.Module);
  const FS = mod.FS || (mod.h5wasm && mod.h5wasm.FS) || (Module && Module.FS);
  if(!File || !Module || !FS) throw new Error('h5wasm API not found');
  const path = `/import_${Date.now()}.h5`;
  const buf = new Uint8Array(await file.arrayBuffer());
  FS.writeFile(path, buf);
  let f=null;
  try{
    f = new File(path, 'r');
    const root = f; const grp = root.get('/results'); if(!grp) throw new Error('results group not found');
    const names = grp.keys(); const cols = {};
    for(const name of names){
      const ds = grp.get(name); if(!ds || ds.type!=='Dataset') continue;
      const val = ds.value; // typed array or list
      cols[name] = Array.isArray(val) ? val : Array.from(val||[]);
    }
    const rows = columnsToRows(cols);
    const {added, renamed} = mergeRowsIntoResults(rows);
    alert(`Imported ${added} row(s) from HDF5.${renamed? ` Renamed ${renamed} duplicate key(s).`:''}`);
    plotWaveforms();
  } finally {
    if(f) f.close();
  }
}

function buildResultsRows(){
  const rows=[];
  const re = /^evt_(\w+)_([^_]+)_ch(\d+)(?:_\d+)?$/;
  for(const [k,rec] of Object.entries(state.results)){
    const m = k.match(re); if(!m) continue;
    const [_, event, fit, ch] = m;
    const row = { event, channel: parseInt(ch,10), fit_type: fit, inverted: !!rec.inverted, dataset_id: rec.dataset_id||state.datasetId||'' };
    if(Array.isArray(rec.param_names) && Array.isArray(rec.params)){
      for(let i=0;i<rec.param_names.length;i++) row[rec.param_names[i]] = rec.params[i];
    }
    if(rec.value!=null) row.value = rec.value;
    if(rec.t_at!=null) row.t_at = rec.t_at;
    if(rec.charge!=null) row.charge = rec.charge;
    if(rec.mass!=null) row.mass = rec.mass;
    if(rec.radius!=null) row.radius = rec.radius;
    rows.push(row);
  }
  return rows;
}

function downloadCSVFromRows(rows, filename){
  if(rows.length===0){ alert('No fits to export.'); return; }
  const cols = Array.from(rows.reduce((s,r)=>{ Object.keys(r).forEach(k=>s.add(k)); return s; }, new Set()));
  const header = cols.join(',');
  const lines = rows.map(r=> cols.map(c=>{
    const v = r[c]; if(v==null) return '';
    if(typeof v==='number') return String(v);
    const t = String(v).replace(/"/g,'""'); return '"'+t+'"';
  }).join(','));
  const csv = [header, ...lines].join('\n');
  const blob = new Blob([csv], {type:'text/csv'});
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download=filename; a.click(); URL.revokeObjectURL(a.href);
}

byId('export-h5-btn').addEventListener('click', async ()=>{
  try{
    await exportHDF5();
  } catch(e){
    console.warn('HDF5 export failed', e);
    if(confirm('HDF5 export failed. Download CSV fallback instead?')){
      const rows = buildResultsRows(); downloadCSVFromRows(rows, 'fit_results.csv');
    } else {
      alert('Ensure h5wasm ESM bundle is served and try again.');
    }
  }
});

byId('import-h5-btn').addEventListener('click', async ()=>{ byId('import-h5-input').click(); });

byId('import-h5-input').addEventListener('change', async (e)=>{
  const f = e.target.files && e.target.files[0]; if(!f) return;
  try { await importHDF5File(f); } catch (err){ alert('HDF5 import failed: '+String(err&&err.message||err)); }
  e.target.value='';
});

// Fit info panel (simple)
byId('show-info-btn').addEventListener('click', ()=>{
  const panel=byId('fit-info'); panel.classList.remove('hidden');
  const table=byId('fit-table'); const summary=byId('fit-summary');
  const rows=[]; const evt=state.currentEvent;
  const onlyDs = byId('fit-info-only-ds').checked;
  const curDs = state.datasetId;
  for(const [k,rec] of Object.entries(state.results)){
    const m = k.match(/^evt_(\w+)_([^_]+)_ch(\d+)/);
    if(!m) continue; const event=m[1], fit=m[2], ch=m[3];
    if(onlyDs && curDs && rec.dataset_id && rec.dataset_id !== curDs) continue;
    const params = rec.params.map(v=>Number(v).toExponential(3)).join(', ');
    const inv = rec.inverted? 'Y':'N';
    const val = (rec.value!=null && isFinite(rec.value))? Number(rec.value).toExponential(3): '';
    const tat = (rec.t_at!=null && isFinite(rec.t_at))? Number(rec.t_at).toExponential(3): '';
    const q = (rec.charge!=null && isFinite(rec.charge))? Number(rec.charge).toExponential(3): '';
    const mkg = (rec.mass!=null && isFinite(rec.mass))? Number(rec.mass).toExponential(3): '';
    const r = (rec.radius!=null && isFinite(rec.radius))? Number(rec.radius).toExponential(3): '';
    rows.push(`<tr><td>${event}</td><td>${ch}</td><td>${fit}</td><td>${inv}</td><td>${val}</td><td>${tat}</td><td>${q}</td><td>${mkg}</td><td>${r}</td><td>${params}</td></tr>`);
  }
  table.innerHTML = `<table><thead><tr><th>Event</th><th>Ch</th><th>Fit</th><th>Inv</th><th>Value</th><th>t_at</th><th>Charge</th><th>Mass</th><th>Radius</th><th>Params</th></tr></thead><tbody>${rows.join('')}</tbody></table>`;
  summary.textContent = evt? `Current: ${evt}`: '';
});
byId('fit-info-close').addEventListener('click', ()=> byId('fit-info').classList.add('hidden'));
byId('fit-info-only-ds').addEventListener('change', ()=>{ if(!byId('fit-info').classList.contains('hidden')) byId('show-info-btn').click(); });

// Theme toggle
byId('theme-btn').addEventListener('click', ()=>{
  if(state.theme==='dark'){
    state.theme='light';
    document.documentElement.style.setProperty('--bg','#f0f0f0');
    document.documentElement.style.setProperty('--fg','#222');
    document.documentElement.style.setProperty('--panel','#ffffff');
    document.documentElement.style.setProperty('--btn-bg','#e0e0e0');
  } else {
    state.theme='dark';
    document.documentElement.style.setProperty('--bg','#1e1e1e');
    document.documentElement.style.setProperty('--fg','#f0f0f0');
    document.documentElement.style.setProperty('--panel','#2b2b2b');
    document.documentElement.style.setProperty('--btn-bg','#333333');
  }
  plotWaveforms();
});

// Help
byId('help-btn').addEventListener('click', ()=>{
  alert('Load CSV data (C<ch>-<event>.csv). Use SG filter, run fits, batch, and feature scan. Results export as JSON.');
});

// Keyboard shortcuts
document.addEventListener('keydown', (e)=>{
  if(e.target && (e.target.tagName==='INPUT' || e.target.tagName==='SELECT' || e.target.isContentEditable)) return;
  const k = e.key.toLowerCase();
  const fitSelEl = byId('fit-select');
  const setFit = (name)=>{ const opt = Array.from(fitSelEl.options).find(o=>o.value.toLowerCase()===name.toLowerCase()); if(opt){ fitSelEl.value=opt.value; e.preventDefault(); }};
  if(k==='q') setFit('QD3Fit');
  else if(k==='d') setFit('QDMFit');
  else if(k==='c') setFit('CSA_pulse');
  else if(k==='w') setFit('skew_gaussian');
  else if(k==='g') setFit('gaussian');
  else if(k==='x') setFit('low_pass_max');
  else if(k==='i'){ const cb=byId('fit-invert'); cb.checked=!cb.checked; e.preventDefault(); }
  else if(k==='s'){ byId('sg-btn').click(); e.preventDefault(); }
  else if(k==='a'){ byId('adjust-btn').click(); e.preventDefault(); }
  else if(k==='r'){ byId('clear-fit-btn').click(); e.preventDefault(); }
  else if(k===',' ){ byId('prev-btn').click(); e.preventDefault(); }
  else if(k==='.') { byId('next-btn').click(); e.preventDefault(); }
});

// Clear all fits for current event/channel
byId('clear-chan-btn').addEventListener('click', ()=>{
  const ev=state.currentEvent; if(!ev) return; const ch=parseInt(byId('fit-chan').value,10);
  let removed=0; const re = new RegExp(`^evt_${ev}_.+_ch${ch}(?:_\\d+)?$`);
  for(const k of Object.keys(state.results)){
    if(re.test(k)){ delete state.results[k]; removed++; }
  }
  alert(`Removed ${removed} fit(s) for ch${ch}.`);
  plotWaveforms();
});
