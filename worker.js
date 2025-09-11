/* Compute Worker for Oscilloscope Palantir (Web)
 * Handles heavy DSP and fitting off the main thread.
 */

// Utilities
function mean(arr){ let s=0; for(let i=0;i<arr.length;i++) s+=arr[i]; return s/(arr.length||1); }
function arrMax(arr){ let m=-Infinity, idx=-1; for(let i=0;i<arr.length;i++){ const v=arr[i]; if(v>m){ m=v; idx=i; } } return {m, idx}; }
function arrMin(arr){ let m=Infinity, idx=-1; for(let i=0;i<arr.length;i++){ const v=arr[i]; if(v<m){ m=v; idx=i; } } return {m, idx}; }
function erf(x){ const sign=x<0?-1:1; x=Math.abs(x); const a1=0.254829592,a2=-0.284496736,a3=1.421413741,a4=-1.453152027,a5=1.061405429,p=0.3275911; const t=1/(1+p*x); const y=1-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*Math.exp(-x*x); return sign*y; }

// Savitzky–Golay filter (odd window, poly=2)
function savgol(y, window, poly){
  window|=0; if(window<5) window=5; if(window%2===0) window+=1; const n=y.length; const m=(window-1)>>1; const out=new Float64Array(n);
  function weights(m){
    const X=[]; for(let i=-m;i<=m;i++){ X.push([1,i,i*i]); }
    const A=[[0,0,0],[0,0,0],[0,0,0]]; for(const r of X){ for(let i=0;i<3;i++) for(let j=0;j<3;j++) A[i][j]+=r[i]*r[j]; }
    const [a,b,c,d,e,f,g,h,k] = [A[0][0],A[0][1],A[0][2],A[1][0],A[1][1],A[1][2],A[2][0],A[2][1],A[2][2]];
    const A_=e*k-f*h, B_=-(d*k-f*g), C_=d*h-e*g, D_=-(b*k-c*h), E_=a*k-c*g, F_=-(a*h-b*g), G_=b*f-c*e, H_=-(a*f-b*e), K_=a*e-b*d; const det=a*A_+b*B_+c*C_;
    const inv=[[A_/det,D_/det,G_/det],[B_/det,E_/det,H_/det],[C_/det,F_/det,K_/det]]; const c0=[inv[0][0],inv[1][0],inv[2][0]];
    const W=[]; for(let i=-m;i<=m;i++){ W.push(c0[0]+c0[1]*i+c0[2]*i*i); } return W;
  }
  const W=weights(m);
  for(let i=0;i<n;i++){
    let acc=0,ws=0; for(let k=-m;k<=m;k++){ let idx=i+k; if(idx<0) idx=0; if(idx>=n) idx=n-1; const w=W[k+m]; acc+=w*y[idx]; ws+=w; }
    out[i]=acc/ws; }
  return out;
}

// Fit function library
const FIT_LIB = {
  linear: { names:['m','b'], func:(T,m,b)=> T.map(t=>m*t+b) },
  QD3Fit: { names:['t0','q','v','C'], func:(T,t0,q,v,C)=>{ v=v*100; const PUT_l=20.0,PUT_d=0.90; const tau1=0.13e-3,tau2=0.15e-3; const dt1=PUT_d/v, dt2=(PUT_l-2*PUT_d)/v; const out=new Float64Array(T.length); for(let i=0;i<T.length;i++){ const t=T[i]; let s=0; if(t>t0 && t<t0+dt1){ s=(q/dt1)*(t-t0);} else if(t>=t0+dt1 && t<t0+dt1+dt2){ s=q*Math.exp(-(t-(t0+dt1))/tau1);} else if(t>=t0+dt1+dt2 && t<t0+2*dt1+dt2){ const qq=q*Math.exp(-dt2/tau1); s= qq - (q/dt1)*(t-(t0+dt1+dt2)); } else if(t>=t0+2*dt1+dt2){ const qq2 = q*Math.exp(-dt2/tau1) - (q/dt1)*( (2*dt1+dt2) - (dt1+dt2) ); s = qq2 * Math.exp(-(t-(t0+2*dt1+dt2))/tau2);} out[i]=s+C;} return out; }},
  QDMFit: { names:['t0','q','v','C'], func:(T,t0,q,v,C)=>{ v=v*100; const PUT_l=15.3,PUT_d=2.0; const tau1=0.4e-3,tau2=0.4e-3; const dt1=PUT_d/v, dt2=(PUT_l-2*PUT_d)/v; const out=new Float64Array(T.length); for(let i=0;i<T.length;i++){ const t=T[i]; let s=0; if(t>t0 && t<t0+dt1){ s=(q/dt1)*(t-t0);} else if(t>t0+dt1 && t<t0+dt1+dt2){ s=q*Math.exp(-(t-(t0+dt1))/tau1);} else if(t>t0+dt1+dt2 && t<t0+2*dt1+dt2){ const qq=q*Math.exp(-dt2/tau1); s= qq - (q/dt1)*(t-(t0+dt1+dt2)); } else if(t>t0+2*dt1+dt2){ const qq2 = q*Math.exp(-dt2/tau1) - (q/dt1)*( (2*dt1+dt2) - (dt1+dt2) ); s = qq2 * Math.exp(-(t-(t0+2*dt1+dt2))/tau2);} out[i]=s+C;} return out; }},
  CSA_pulse: { names:['t0','C0','C1','C2','T0','T1','T2','C'], func:(x,t0,C0,C1,C2,T0,T1,T2,C)=>{ const e=1e-12; const T0p=Math.max(Math.abs(T0),e),T1p=Math.max(Math.abs(T1),e),T2p=Math.max(Math.abs(T2),e); const out=new Float64Array(x.length); for(let i=0;i<x.length;i++){ const xx=x[i]; const left=(xx<=t0?1:0); const right=(xx>=t0?1:0); const gauss=Math.exp(-((xx-t0)*(xx-t0))/(T0p*T0p)); const rise=(1.0-Math.exp(-(xx-t0)/T1p)); const decay=Math.exp(-(xx-t0)/T2p); let val=(C0 - left*C1*gauss + right*( C2*rise*decay - C1 )) + C; if(!Number.isFinite(val)) val=C0+C; out[i]=val; } return out; }},
  skew_gaussian: { names:['A','xi','omega','alpha','C'], func:(x,A,xi,omega,alpha,C)=>{ const out=new Float64Array(x.length); const invsqrt2pi=1/Math.sqrt(2*Math.PI); const op=Math.max(Math.abs(omega),1e-12); for(let i=0;i<x.length;i++){ const t=(x[i]-xi)/op; const phi=invsqrt2pi*Math.exp(-0.5*t*t); const Phi=0.5*(1+erf(alpha*t/Math.sqrt(2))); out[i]=A*2/op*phi*Phi + C; } return out; }},
  gaussian: { names:['A','mu','sigma','C'], func:(x,A,mu,sigma,C)=>{ const out=new Float64Array(x.length); for(let i=0;i<x.length;i++){ const t=(x[i]-mu)/sigma; out[i]=A*Math.exp(-0.5*t*t)+C; } return out; }},
};

// Nelder–Mead least squares
function nelderMead(f, x0, opts={}){
  const maxIter=opts.maxIter||800, ftol=opts.ftol||1e-8; const alpha=1,gamma=2,rho=0.5,sigma=0.5; const n=x0.length;
  let simplex=new Array(n+1); simplex[0]=x0.slice(); const step=opts.step||1e-2; for(let i=0;i<n;i++){ const xi=x0.slice(); xi[i]=xi[i]!==0? xi[i]*(1+step): step; simplex[i+1]=xi; }
  function sort(){ simplex.sort((a,b)=>f(a)-f(b)); }
  sort(); let it=0; while(it++<maxIter){ sort(); const best=simplex[0], worst=simplex[n], secondWorst=simplex[n-1]; const fbest=f(best), fworst=f(worst); if(Math.abs(fworst-fbest)<ftol*(1+Math.abs(fbest))) break; const centroid=new Array(n).fill(0); for(let i=0;i<n;i++) for(let j=0;j<n;j++) centroid[j]+=simplex[i][j]; for(let j=0;j<n;j++) centroid[j]/=n; const xr=centroid.map((c,j)=>c+alpha*(c-worst[j])); if(f(xr)<f(secondWorst)){ simplex[n]=xr; continue; } if(f(xr)<f(best)){ const xe=centroid.map((c,j)=>c+gamma*(xr[j]-c)); simplex[n]= f(xe)<f(xr)? xe: xr; continue; } const xc=centroid.map((c,j)=>c+rho*(worst[j]-c)); if(f(xc)<f(worst)){ simplex[n]=xc; continue; } for(let i=1;i<n+1;i++) simplex[i]=simplex[0].map((b,j)=> b + sigma*(simplex[i][j]-b)); }
  sort(); return {x: simplex[0], f: f(simplex[0]), iters: it};
}

function lsqFit(model, t, y, p0){
  const f=(p)=>{ try{ const yy=model.func(t, ...p); if(!yy || yy.length!==y.length) return 1e24; let s=0; for(let i=0;i<y.length;i++){ const d=yy[i]-y[i]; if(!Number.isFinite(d)) return 1e24; s+=d*d; } return Number.isFinite(s)? s: 1e24; } catch{ return 1e24; } };
  const res=nelderMead(f, p0, {maxIter:600, step:1e-2, ftol:1e-8});
  return res.x;
}

function handleFit(msg){
  const {fitName, t, y, p0} = msg; const model = FIT_LIB[fitName];
  if(!model || !model.func) throw new Error('Unknown or non-functional fit: '+fitName);
  const params = lsqFit(model, t, y, p0);
  return {params};
}

function handleSG(msg){ const {y, window, poly} = msg; const yy = savgol(y, window|0, poly|0); return {y: Array.from(yy)}; }

function handleLowPassMax(msg){
  const {t, y, width, invert} = msg; const w = Math.max(5, (width|0)|0) | 1; // odd
  const yproc = invert? y.map(v=>-v): y.slice();
  const ys = savgol(yproc, w, 2);
  const {m:peak, idx} = arrMax(ys);
  const t_at = t[Math.max(0, Math.min(idx, t.length-1))];
  const val = invert? -peak: peak;
  return {value: val, t_at, width: w};
}

self.onmessage = (e)=>{
  const {id, type} = e.data || {};
  try{
    let result;
    switch(type){
      case 'fit':        result = handleFit(e.data); break;
      case 'sg':         result = handleSG(e.data); break;
      case 'low_pass_max': result = handleLowPassMax(e.data); break;
      case 'parse_trc':  result = handleParseTrc(e.data); break;
      default: throw new Error('Unknown task type: '+type);
    }
    self.postMessage({id, ok:true, result});
  } catch (err){
    self.postMessage({id, ok:false, error: String(err && err.message || err)});
  }
};

// TRC parsing in worker
function findWavedescOffset(buf){
  const bytes = new Uint8Array(buf, 0, Math.min(8192, buf.byteLength));
  const pat = new TextEncoder().encode('WAVEDESC');
  for(let i=0;i<=bytes.length-pat.length;i++){
    let ok=true; for(let j=0;j<pat.length;j++){ if(bytes[i+j]!==pat[j]){ ok=false; break; } }
    if(ok) return i;
  }
  return -1;
}
function detectEndian(view, base){
  const be = view.getUint16(base+34, false);
  const le = view.getUint16(base+34, true);
  if(le===1 || be===256) return true;
  if(le===0 || be===0) return false;
  return true;
}
function handleParseTrc(msg){
  const {buf} = msg; if(!buf) throw new Error('No buffer');
  const base = findWavedescOffset(buf); if(base<0) throw new Error('WAVEDESC not found');
  const view = new DataView(buf);
  const little = detectEndian(view, base);
  const getI16 = (off)=> view.getInt16(base+off, little);
  const getU16 = (off)=> view.getUint16(base+off, little);
  const getI32 = (off)=> view.getInt32(base+off, little);
  const getF32 = (off)=> view.getFloat32(base+off, little);
  const getF64 = (off)=> view.getFloat64(base+off, little);
  // Types and lengths
  const commTypeLE = view.getUint16(base+32, true);
  const commTypeBE = view.getUint16(base+32, false);
  const is16 = (commTypeLE===1 || commTypeBE===256 || (commTypeLE!==0 && commTypeBE!==0));
  const lWAVEDESC = getI32(36);
  const lUSERTEXT = getI32(40);
  const lTRIGTIME = getI32(48);
  const lRISTIME  = getI32(52);
  const lWAVE1    = getI32(60);
  // Scales
  const VERTICAL_GAIN   = getF32(156);
  const VERTICAL_OFFSET = getF32(160);
  const HORIZ_INTERVAL  = getF32(176);
  const HORIZ_OFFSET    = getF64(180);
  // Data start and count
  const dataOff = base + lWAVEDESC + lUSERTEXT + lTRIGTIME + lRISTIME;
  const bytesPer = is16? 2: 1;
  const nFromLen = lWAVE1>0? Math.floor(lWAVE1/bytesPer): 0;
  const maxN = Math.max(0, Math.floor((buf.byteLength - dataOff)/bytesPer));
  const n = Math.max(0, Math.min(nFromLen || maxN, maxN));
  if(n<=0) throw new Error('No waveform samples');
  // Read samples
  const yraw = new Int32Array(n);
  if(is16){ for(let i=0;i<n;i++){ yraw[i]=view.getInt16(dataOff + 2*i, little); } }
  else { for(let i=0;i<n;i++){ yraw[i]=view.getInt8(dataOff + i); } }
  const y = new Float64Array(n); for(let i=0;i<n;i++){ y[i] = VERTICAL_GAIN * yraw[i] - VERTICAL_OFFSET; }
  const t = new Float64Array(n); for(let i=0;i<n;i++){ t[i] = (i)*HORIZ_INTERVAL + HORIZ_OFFSET; }
  return { t: Array.from(t), y: Array.from(y), meta: { endian: little? 'LE':'BE', comm_type: is16? 'int16':'int8' } };
}
