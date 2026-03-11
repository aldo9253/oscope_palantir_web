/* Compute Worker for Oscilloscope Viewer (Web)
 * Handles SG filtering and TRC parsing off the main thread.
 */

function savgol(y, window, poly){
  window |= 0;
  if(window < 5) window = 5;
  if(window % 2 === 0) window += 1;

  const n = y.length;
  const m = (window - 1) >> 1;
  const out = new Float64Array(n);

  function weights(half){
    const X = [];
    for(let i=-half;i<=half;i++) X.push([1, i, i*i]);

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
    for(let i=-half;i<=half;i++) W.push(c0[0] + c0[1]*i + c0[2]*i*i);
    return W;
  }

  const W = weights(m);
  for(let i=0;i<n;i++){
    let acc = 0;
    let ws = 0;
    for(let k=-m;k<=m;k++){
      let idx = i + k;
      if(idx < 0) idx = 0;
      if(idx >= n) idx = n - 1;
      const w = W[k + m];
      acc += w * y[idx];
      ws += w;
    }
    out[i] = acc / ws;
  }

  return out;
}

function handleSG(msg){
  const { y, window, poly } = msg;
  const yy = savgol(y, window | 0, poly | 0);
  return { y: Array.from(yy) };
}

function findWavedescOffset(buf){
  const bytes = new Uint8Array(buf, 0, Math.min(8192, buf.byteLength));
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

function detectEndian(view, base){
  const be = view.getUint16(base+34, false);
  const le = view.getUint16(base+34, true);
  if(le === 1 || be === 256) return true;
  if(le === 0 || be === 0) return false;
  return true;
}

function handleParseTrc(msg){
  const { buf } = msg;
  if(!buf) throw new Error('No buffer');

  const base = findWavedescOffset(buf);
  if(base < 0) throw new Error('WAVEDESC not found');

  const view = new DataView(buf);
  const little = detectEndian(view, base);
  const getI32 = (off)=> view.getInt32(base+off, little);
  const getF32 = (off)=> view.getFloat32(base+off, little);
  const getF64 = (off)=> view.getFloat64(base+off, little);

  const commTypeLE = view.getUint16(base+32, true);
  const commTypeBE = view.getUint16(base+32, false);
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
  const bytesPer = is16 ? 2 : 1;
  const nFromLen = lWAVE1 > 0 ? Math.floor(lWAVE1 / bytesPer) : 0;
  const maxN = Math.max(0, Math.floor((buf.byteLength - dataOff) / bytesPer));
  const n = Math.max(0, Math.min(nFromLen || maxN, maxN));
  if(n <= 0) throw new Error('No waveform samples');

  const yraw = new Int32Array(n);
  if(is16){
    for(let i=0;i<n;i++) yraw[i] = view.getInt16(dataOff + 2*i, little);
  } else {
    for(let i=0;i<n;i++) yraw[i] = view.getInt8(dataOff + i);
  }

  const y = new Float64Array(n);
  const t = new Float64Array(n);
  for(let i=0;i<n;i++){
    y[i] = VERTICAL_GAIN * yraw[i] - VERTICAL_OFFSET;
    t[i] = i * HORIZ_INTERVAL + HORIZ_OFFSET;
  }

  return { t: Array.from(t), y: Array.from(y), meta: { endian: little ? 'LE' : 'BE', comm_type: is16 ? 'int16' : 'int8' } };
}

self.onmessage = (e)=>{
  const { id, type } = e.data || {};
  try{
    let result;
    if(type === 'sg') result = handleSG(e.data);
    else if(type === 'parse_trc') result = handleParseTrc(e.data);
    else throw new Error('Unknown task type: ' + type);
    self.postMessage({ id, ok: true, result });
  } catch (err){
    self.postMessage({ id, ok: false, error: String(err && err.message || err) });
  }
};
