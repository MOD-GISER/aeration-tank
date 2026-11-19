(function(){
'use strict';
const D=window.TANK_DATA, C=window.TANK_CONFIG, $=id=>document.getElementById(id);
const rad=d=>d*Math.PI/180, deg=r=>r*180/Math.PI;
// ---------------------------------------------------------------- coordinates
const G=D.georef; let mode=C.defaults.scaleToSurvey?'similarity':'rigid'; let zEx=C.defaults.zExag;
const REF={E:G.similarity.cd[0],N:G.similarity.cd[1]};
function tf(){const t=G[mode];const th=rad(t.theta_deg);return{s:t.s,c:Math.cos(th),sn:Math.sin(th),cs:t.cs,cd:t.cd};}
let T=tf();
function loc2en(x,y){const dx=x-T.cs[0],dy=y-T.cs[1];return[T.cd[0]+T.s*(dx*T.c-dy*T.sn),T.cd[1]+T.s*(dx*T.sn+dy*T.c)];}
function en2loc(E,N){const u=(E-T.cd[0])/T.s,v=(N-T.cd[1])/T.s;return[T.cs[0]+u*T.c+v*T.sn,T.cs[1]-u*T.sn+v*T.c];}
function loc2w(x,y,z){const e=loc2en(x,y);return[e[0]-REF.E,e[1]-REF.N,z*zEx];}
function utm2ll(E,N){const a=6378137,f=1/298.257223563,k0=0.9996,E0=500000,lam0=rad(27);
 const n=f/(2-f),A=a/(1+n)*(1+n*n/4+Math.pow(n,4)/64);
 const b=[0,n/2-2*n*n/3+37*Math.pow(n,3)/96-Math.pow(n,4)/360,n*n/48+Math.pow(n,3)/15-437*Math.pow(n,4)/1440,17*Math.pow(n,3)/480-37*Math.pow(n,4)/840,4397*Math.pow(n,4)/161280];
 const d=[0,2*n-2*n*n/3-2*Math.pow(n,3)+116*Math.pow(n,4)/45,7*n*n/3-8*Math.pow(n,3)/5-227*Math.pow(n,4)/45,56*Math.pow(n,3)/15-136*Math.pow(n,4)/35,4279*Math.pow(n,4)/630];
 const xi=N/(k0*A),eta=(E-E0)/(k0*A);let x2=xi,e2=eta;
 for(let j=1;j<=4;j++){x2-=b[j]*Math.sin(2*j*xi)*Math.cosh(2*j*eta);e2-=b[j]*Math.cos(2*j*xi)*Math.sinh(2*j*eta);}
 const chi=Math.asin(Math.sin(x2)/Math.cosh(e2));let phi=chi;for(let j=1;j<=4;j++)phi+=d[j]*Math.sin(2*j*chi);
 return[deg(phi),deg(lam0+Math.atan2(Math.sinh(e2),Math.cos(x2)))];}
window.__utm2ll=utm2ll;
// ---------------------------------------------------------------- math
const M={
 persp(f,a,n,fa){const t=1/Math.tan(f/2);return[t/a,0,0,0,0,t,0,0,0,0,(fa+n)/(n-fa),-1,0,0,2*fa*n/(n-fa),0];},
 ortho(w,h,n,f){return[1/w,0,0,0,0,1/h,0,0,0,0,-2/(f-n),0,0,0,-(f+n)/(f-n),1];},
 look(e,c,u){let zx=e[0]-c[0],zy=e[1]-c[1],zz=e[2]-c[2];let l=Math.hypot(zx,zy,zz);zx/=l;zy/=l;zz/=l;
  let xx=u[1]*zz-u[2]*zy,xy=u[2]*zx-u[0]*zz,xz=u[0]*zy-u[1]*zx;l=Math.hypot(xx,xy,xz);xx/=l;xy/=l;xz/=l;
  const yx=zy*xz-zz*xy,yy=zz*xx-zx*xz,yz=zx*xy-zy*xx;
  return[xx,yx,zx,0,xy,yy,zy,0,xz,yz,zz,0,-(xx*e[0]+xy*e[1]+xz*e[2]),-(yx*e[0]+yy*e[1]+yz*e[2]),-(zx*e[0]+zy*e[1]+zz*e[2]),1];},
 mul(a,b){const o=new Array(16);for(let i=0;i<4;i++)for(let j=0;j<4;j++){let s=0;for(let k=0;k<4;k++)s+=a[k*4+j]*b[i*4+k];o[i*4+j]=s;}return o;},
 inv(m){const i=new Array(16),a=m;
  i[0]=a[5]*a[10]*a[15]-a[5]*a[11]*a[14]-a[9]*a[6]*a[15]+a[9]*a[7]*a[14]+a[13]*a[6]*a[11]-a[13]*a[7]*a[10];
  i[4]=-a[4]*a[10]*a[15]+a[4]*a[11]*a[14]+a[8]*a[6]*a[15]-a[8]*a[7]*a[14]-a[12]*a[6]*a[11]+a[12]*a[7]*a[10];
  i[8]=a[4]*a[9]*a[15]-a[4]*a[11]*a[13]-a[8]*a[5]*a[15]+a[8]*a[7]*a[13]+a[12]*a[5]*a[11]-a[12]*a[7]*a[9];
  i[12]=-a[4]*a[9]*a[14]+a[4]*a[10]*a[13]+a[8]*a[5]*a[14]-a[8]*a[6]*a[13]-a[12]*a[5]*a[10]+a[12]*a[6]*a[9];
  i[1]=-a[1]*a[10]*a[15]+a[1]*a[11]*a[14]+a[9]*a[2]*a[15]-a[9]*a[3]*a[14]-a[13]*a[2]*a[11]+a[13]*a[3]*a[10];
  i[5]=a[0]*a[10]*a[15]-a[0]*a[11]*a[14]-a[8]*a[2]*a[15]+a[8]*a[3]*a[14]+a[12]*a[2]*a[11]-a[12]*a[3]*a[10];
  i[9]=-a[0]*a[9]*a[15]+a[0]*a[11]*a[13]+a[8]*a[1]*a[15]-a[8]*a[3]*a[13]-a[12]*a[1]*a[11]+a[12]*a[3]*a[9];
  i[13]=a[0]*a[9]*a[14]-a[0]*a[10]*a[13]-a[8]*a[1]*a[14]+a[8]*a[2]*a[13]+a[12]*a[1]*a[10]-a[12]*a[2]*a[9];
  i[2]=a[1]*a[6]*a[15]-a[1]*a[7]*a[14]-a[5]*a[2]*a[15]+a[5]*a[3]*a[14]+a[13]*a[2]*a[7]-a[13]*a[3]*a[6];
  i[6]=-a[0]*a[6]*a[15]+a[0]*a[7]*a[14]+a[4]*a[2]*a[15]-a[4]*a[3]*a[14]-a[12]*a[2]*a[7]+a[12]*a[3]*a[6];
  i[10]=a[0]*a[5]*a[15]-a[0]*a[7]*a[13]-a[4]*a[1]*a[15]+a[4]*a[3]*a[13]+a[12]*a[1]*a[7]-a[12]*a[3]*a[5];
  i[14]=-a[0]*a[5]*a[14]+a[0]*a[6]*a[13]+a[4]*a[1]*a[14]-a[4]*a[2]*a[13]-a[12]*a[1]*a[6]+a[12]*a[2]*a[5];
  i[3]=-a[1]*a[6]*a[11]+a[1]*a[7]*a[10]+a[5]*a[2]*a[11]-a[5]*a[3]*a[10]-a[9]*a[2]*a[7]+a[9]*a[3]*a[6];
  i[7]=a[0]*a[6]*a[11]-a[0]*a[7]*a[10]-a[4]*a[2]*a[11]+a[4]*a[3]*a[10]+a[8]*a[2]*a[7]-a[8]*a[3]*a[6];
  i[11]=-a[0]*a[5]*a[11]+a[0]*a[7]*a[9]+a[4]*a[1]*a[11]-a[4]*a[3]*a[9]-a[8]*a[1]*a[7]+a[8]*a[3]*a[5];
  i[15]=a[0]*a[5]*a[10]-a[0]*a[6]*a[9]-a[4]*a[1]*a[10]+a[4]*a[2]*a[9]+a[8]*a[1]*a[6]-a[8]*a[2]*a[5];
  let d=a[0]*i[0]+a[1]*i[4]+a[2]*i[8]+a[3]*i[12];d=1/d;return i.map(x=>x*d);},
 tp(m,x,y,z,w){return[m[0]*x+m[4]*y+m[8]*z+m[12]*w,m[1]*x+m[5]*y+m[9]*z+m[13]*w,m[2]*x+m[6]*y+m[10]*z+m[14]*w,m[3]*x+m[7]*y+m[11]*z+m[15]*w];}
};
// ---------------------------------------------------------------- data decode
function b64(s){const bin=atob(s),u=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)u[i]=bin.charCodeAt(i);return u.buffer;}
const EL=D.elements.map((e,k)=>{const v=new Float32Array(b64(e.v));const i=e.it===16?new Uint16Array(b64(e.i)):new Uint32Array(b64(e.i));return{id:k,n:e.n,l:e.l,g:e.g,s:e.s,v,i,nt:e.nt};});
// ---------------------------------------------------------------- gl
const cv=$('c');const gl=cv.getContext('webgl',{antialias:true,alpha:true,preserveDrawingBuffer:false})||cv.getContext('experimental-webgl');
if(!gl){document.body.innerHTML='<p style="padding:20px;font:16px sans-serif">WebGL is not available in this browser.</p>';return;}
function sh(t,src){const s=gl.createShader(t);gl.shaderSource(s,src);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw new Error(gl.getShaderInfoLog(s));return s;}
function prog(vs,fs){const p=gl.createProgram();gl.attachShader(p,sh(gl.VERTEX_SHADER,vs));gl.attachShader(p,sh(gl.FRAGMENT_SHADER,fs));gl.linkProgram(p);if(!gl.getProgramParameter(p,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(p));return p;}
const P1=prog('attribute vec3 p;attribute vec3 n;uniform mat4 vp;varying vec3 vn;varying float vz;void main(){vn=n;vz=p.z;gl_Position=vp*vec4(p,1.0);}',
 'precision mediump float;varying vec3 vn;varying float vz;uniform vec3 col;uniform float clipz;uniform float tint;void main(){if(vz>clipz)discard;vec3 n=normalize(vn);vec3 l1=normalize(vec3(-0.45,-0.35,0.82));vec3 l2=normalize(vec3(0.6,0.5,0.25));float d=abs(dot(n,l1))*0.62+abs(dot(n,l2))*0.22;vec3 c=col*(0.50+d*0.8);c=mix(c,vec3(1.0,0.85,0.1),tint);gl_FragColor=vec4(c,1.0);}');
const P2=prog('attribute vec3 p;uniform mat4 vp;void main(){gl_Position=vp*vec4(p,1.0);}','precision mediump float;uniform vec3 col;void main(){gl_FragColor=vec4(col,1.0);}');
const U1={vp:gl.getUniformLocation(P1,'vp'),col:gl.getUniformLocation(P1,'col'),clipz:gl.getUniformLocation(P1,'clipz'),tint:gl.getUniformLocation(P1,'tint')};
const U2={vp:gl.getUniformLocation(P2,'vp'),col:gl.getUniformLocation(P2,'col')};
function hex(h){return[parseInt(h.substr(1,2),16)/255,parseInt(h.substr(3,2),16)/255,parseInt(h.substr(5,2),16)/255];}
// per element world buffers (rebuilt when transform / exaggeration change)
function rebuild(){
 T=tf();
 for(const e of EL){
  const nt=e.i.length/3;const pos=new Float32Array(nt*9),nor=new Float32Array(nt*9);
  const wv=new Float32Array(e.v.length);
  for(let k=0;k<e.v.length;k+=3){const w=loc2w(e.v[k],e.v[k+1],e.v[k+2]);wv[k]=w[0];wv[k+1]=w[1];wv[k+2]=w[2];}
  let mn=[1e9,1e9,1e9],mx=[-1e9,-1e9,-1e9];
  for(let t=0;t<nt;t++){
   const a=e.i[t*3]*3,b=e.i[t*3+1]*3,c=e.i[t*3+2]*3;
   const ux=wv[b]-wv[a],uy=wv[b+1]-wv[a+1],uz=wv[b+2]-wv[a+2],vx=wv[c]-wv[a],vy=wv[c+1]-wv[a+1],vz=wv[c+2]-wv[a+2];
   let nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;const l=Math.hypot(nx,ny,nz)||1;nx/=l;ny/=l;nz/=l;
   const idx=[a,b,c];
   for(let q=0;q<3;q++){const s=idx[q];for(let d=0;d<3;d++){const val=wv[s+d];pos[t*9+q*3+d]=val;if(val<mn[d])mn[d]=val;if(val>mx[d])mx[d]=val;}nor[t*9+q*3]=nx;nor[t*9+q*3+1]=ny;nor[t*9+q*3+2]=nz;}
  }
  e.pos=pos;e.mn=mn;e.mx=mx;e.nt=nt;
  if(!e.bp){e.bp=gl.createBuffer();e.bn=gl.createBuffer();}
  gl.bindBuffer(gl.ARRAY_BUFFER,e.bp);gl.bufferData(gl.ARRAY_BUFFER,pos,gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER,e.bn);gl.bufferData(gl.ARRAY_BUFFER,nor,gl.STATIC_DRAW);
 }
 buildLines();
 let lo=[1e9,1e9,1e9],hi=[-1e9,-1e9,-1e9];for(const e of EL)for(let d=0;d<3;d++){lo[d]=Math.min(lo[d],e.mn[d]);hi[d]=Math.max(hi[d],e.mx[d]);}
 bounds={lo,hi};
}
let bounds=null,lineBuf=null,lineN=0,ctlBuf=null,ctlN=0,gridLabels=[],ctlLabels=[];
function buildLines(){
 const z0=1.15*zEx,L=[];gridLabels=[];
 for(const g of D.grids){
  const a=g.dir==='v'?loc2w(g.pos,g.a,1.15):loc2w(g.a,g.pos,1.15),b=g.dir==='v'?loc2w(g.pos,g.b,1.15):loc2w(g.b,g.pos,1.15);
  L.push(a[0],a[1],z0+0.05,b[0],b[1],z0+0.05);gridLabels.push({n:g.name,p:a},{n:g.name,p:b});
 }
 lineBuf=lineBuf||gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,lineBuf);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(L),gl.STATIC_DRAW);lineN=L.length/3;
 const K=[];ctlLabels=[];
 for(const p of G.points){const w=loc2w(p.x,p.y,1.15);K.push(w[0],w[1],w[2],w[0],w[1],w[2]+12*zEx/Math.max(zEx,1)+4);ctlLabels.push({n:p.name,p:[w[0],w[1],w[2]+12*zEx/Math.max(zEx,1)+4]});}
 ctlBuf=ctlBuf||gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,ctlBuf);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(K),gl.STATIC_DRAW);ctlN=K.length/3;
}
// ---------------------------------------------------------------- state / UI
const vis={};for(const k in C.layers)vis[k]=C.layers[k].visible;let showGrid=C.grids.visible,showCtl=C.control.visible;
const cam={tx:0,ty:0,tz:0,dist:300,az:-60,el:35,ortho:C.defaults.ortho};
let cutOn=C.defaults.cutOn,cutZ=9;
function fit(){const b=bounds;cam.tx=(b.lo[0]+b.hi[0])/2;cam.ty=(b.lo[1]+b.hi[1])/2;cam.tz=(b.lo[2]+b.hi[2])/2;cam.dist=Math.hypot(b.hi[0]-b.lo[0],b.hi[1]-b.lo[1])*0.95;}
function setView(v){if(v==='plan'){cam.el=89.9;cam.az=-90;cam.ortho=true;}else if(v==='3d'){cam.el=32;cam.az=-62;cam.ortho=false;}fit();syncUI();}
function matrices(){
 const w=cv.clientWidth,h=cv.clientHeight,asp=w/h;const el=rad(cam.el),az=rad(cam.az);
 const eye=[cam.tx+cam.dist*Math.cos(el)*Math.cos(az),cam.ty+cam.dist*Math.cos(el)*Math.sin(az),cam.tz+cam.dist*Math.sin(el)];
 const V=M.look(eye,[cam.tx,cam.ty,cam.tz],[0,0,1]);
 const P=cam.ortho?M.ortho(cam.dist*0.5*asp,cam.dist*0.5,cam.dist*0.01,cam.dist*8):M.persp(rad(38),asp,cam.dist*0.01,cam.dist*10);
 return{V,P,VP:M.mul(P,V),eye};}
function project(p,m){const q=M.tp(m.VP,p[0],p[1],p[2],1);if(q[3]<=0)return null;return[(q[0]/q[3]*0.5+0.5)*cv.clientWidth,(1-(q[1]/q[3]*0.5+0.5))*cv.clientHeight];}
// ---------------------------------------------------------------- picking
function ray(px,py,m){const inv=M.inv(m.VP);const nx=px/cv.clientWidth*2-1,ny=1-py/cv.clientHeight*2;
 const a=M.tp(inv,nx,ny,-1,1),b=M.tp(inv,nx,ny,1,1);const A=[a[0]/a[3],a[1]/a[3],a[2]/a[3]],B=[b[0]/b[3],b[1]/b[3],b[2]/b[3]];
 let d=[B[0]-A[0],B[1]-A[1],B[2]-A[2]];const l=Math.hypot(d[0],d[1],d[2]);d=d.map(x=>x/l);return{o:A,d};}
function hitBox(r,mn,mx){let t0=0,t1=1e12;for(let i=0;i<3;i++){const inv=1/(r.d[i]||1e-12);let a=(mn[i]-r.o[i])*inv,b=(mx[i]-r.o[i])*inv;if(a>b){const t=a;a=b;b=t;}t0=Math.max(t0,a);t1=Math.min(t1,b);if(t0>t1)return false;}return true;}
function pick(px,py){const m=matrices();const r=ray(px,py,m);let best=1e12,bel=null;const cz=cutOn?cutZ*zEx:1e9;
 for(const e of EL){if(!vis[e.l])continue;if(!hitBox(r,e.mn,e.mx))continue;const P=e.pos;
  for(let t=0;t<e.nt;t++){const o=t*9;
   const e1x=P[o+3]-P[o],e1y=P[o+4]-P[o+1],e1z=P[o+5]-P[o+2],e2x=P[o+6]-P[o],e2y=P[o+7]-P[o+1],e2z=P[o+8]-P[o+2];
   const hx=r.d[1]*e2z-r.d[2]*e2y,hy=r.d[2]*e2x-r.d[0]*e2z,hz=r.d[0]*e2y-r.d[1]*e2x;const a=e1x*hx+e1y*hy+e1z*hz;if(Math.abs(a)<1e-12)continue;
   const f=1/a,sx=r.o[0]-P[o],sy=r.o[1]-P[o+1],sz=r.o[2]-P[o+2];const u=f*(sx*hx+sy*hy+sz*hz);if(u<0||u>1)continue;
   const qx=sy*e1z-sz*e1y,qy=sz*e1x-sx*e1z,qz=sx*e1y-sy*e1x;const v=f*(r.d[0]*qx+r.d[1]*qy+r.d[2]*qz);if(v<0||u+v>1)continue;
   const tt=f*(e2x*qx+e2y*qy+e2z*qz);if(tt>1e-6&&tt<best){const zz=r.o[2]+r.d[2]*tt;if(zz<=cz){best=tt;bel=e;}}}}
 if(!bel)return null;return{w:[r.o[0]+r.d[0]*best,r.o[1]+r.d[1]*best,r.o[2]+r.d[2]*best],el:bel};}
function w2info(w){const E=w[0]+REF.E,N=w[1]+REF.N,Z=w[2]/zEx;const l=en2loc(E,N);const ll=utm2ll(E,N);return{E,N,Z,lat:ll[0],lon:ll[1],lx:l[0],ly:l[1]};}
// ---------------------------------------------------------------- popup
const dp=C.decimals;let picks=[];const pop=$('pop');
function showPick(h){const i=w2info(h.w);picks.push({w:h.w,i,el:h.el});if(picks.length>2)picks.shift();selEl=h.el;renderPop();}
let selEl=null;
function renderPop(){const p=picks[picks.length-1];if(!p){pop.style.display='none';return;}const i=p.i;let extra='';
 if(picks.length===2){const a=picks[0].i,b=i;const dh=Math.hypot(b.E-a.E,b.N-a.N),dz=b.Z-a.Z;extra='<div class="sep"></div><div class="row"><span>From previous point</span></div><div class="row"><span>Horizontal</span><b>'+dh.toFixed(3)+' m</b></div><div class="row"><span>ΔZ</span><b>'+(dz>=0?'+':'')+dz.toFixed(3)+' m</b></div><div class="row"><span>Slope dist.</span><b>'+Math.hypot(dh,dz).toFixed(3)+' m</b></div>';}
 pop.innerHTML='<div class="ph"><div><div class="pn">'+esc(p.el.n)+'</div><div class="pl">'+esc(C.layers[p.el.l].name)+'</div></div><button id="px" title="Close">×</button></div>'+
 '<div class="row"><span>E</span><b>'+i.E.toFixed(dp.en)+'</b></div><div class="row"><span>N</span><b>'+i.N.toFixed(dp.en)+'</b></div><div class="row"><span>Z</span><b>'+i.Z.toFixed(dp.z)+'</b></div>'+
 '<div class="row s"><span>Lat</span><b>'+i.lat.toFixed(dp.ll)+'</b></div><div class="row s"><span>Lon</span><b>'+i.lon.toFixed(dp.ll)+'</b></div><div class="row s"><span>Local x,y</span><b>'+i.lx.toFixed(3)+', '+i.ly.toFixed(3)+'</b></div>'+extra+
 '<div class="btns"><button id="pc">Copy E,N,Z</button><button id="pm">Clear points</button></div><div class="pe">WGS 84 / UTM 35N (EPSG:32635)</div>';
 pop.style.display='block';
 $('px').onclick=()=>{picks=[];selEl=null;renderPop();};
 $('pm').onclick=()=>{picks=[];selEl=null;renderPop();};
 $('pc').onclick=()=>{const t=i.E.toFixed(dp.en)+', '+i.N.toFixed(dp.en)+', '+i.Z.toFixed(dp.z);copy(t);$('pc').textContent='Copied';setTimeout(()=>{const b=$('pc');if(b)b.textContent='Copy E,N,Z';},1200);};}
function copy(t){try{navigator.clipboard.writeText(t);}catch(e){const a=document.createElement('textarea');a.value=t;document.body.appendChild(a);a.select();try{document.execCommand('copy');}catch(_){}document.body.removeChild(a);}}
function esc(s){return String(s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));}
function placePop(m){const p=picks[picks.length-1];if(!p){return;}const s=project(p.w,m);if(!s){pop.style.visibility='hidden';return;}pop.style.visibility='visible';
 const W=cv.clientWidth,H=cv.clientHeight,pw=pop.offsetWidth,ph=pop.offsetHeight;let x=s[0]+16,y=s[1]-ph/2;if(x+pw>W-8)x=s[0]-pw-16;x=Math.max(8,x);y=Math.max(8,Math.min(H-ph-8,y));pop.style.left=x+'px';pop.style.top=y+'px';}
// ---------------------------------------------------------------- draw
const labs=$('labs');let labEls=[];
function ensureLabels(n){while(labEls.length<n){const d=document.createElement('div');d.className='lab';labs.appendChild(d);labEls.push(d);}}
function draw(){
 const W=cv.clientWidth,H=cv.clientHeight,dpr=Math.min(window.devicePixelRatio||1,2);if(cv.width!==Math.round(W*dpr)||cv.height!==Math.round(H*dpr)){cv.width=Math.round(W*dpr);cv.height=Math.round(H*dpr);}
 gl.viewport(0,0,cv.width,cv.height);gl.clearColor(0,0,0,0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.enable(gl.DEPTH_TEST);
 const m=matrices();
 gl.useProgram(P1);gl.uniformMatrix4fv(U1.vp,false,new Float32Array(m.VP));gl.uniform1f(U1.clipz,cutOn?cutZ*zEx:1e9);
 const a=gl.getAttribLocation(P1,'p'),n=gl.getAttribLocation(P1,'n');gl.enableVertexAttribArray(a);gl.enableVertexAttribArray(n);
 for(const e of EL){if(!vis[e.l])continue;const c=hex(C.layers[e.l].color);gl.uniform3fv(U1.col,c);gl.uniform1f(U1.tint,e===selEl?0.45:0);
  gl.bindBuffer(gl.ARRAY_BUFFER,e.bp);gl.vertexAttribPointer(a,3,gl.FLOAT,false,0,0);gl.bindBuffer(gl.ARRAY_BUFFER,e.bn);gl.vertexAttribPointer(n,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.TRIANGLES,0,e.nt*3);}
 gl.disableVertexAttribArray(n);
 gl.useProgram(P2);gl.uniformMatrix4fv(U2.vp,false,new Float32Array(m.VP));const b=gl.getAttribLocation(P2,'p');gl.enableVertexAttribArray(b);
 if(showGrid){gl.uniform3fv(U2.col,hex(C.grids.color));gl.bindBuffer(gl.ARRAY_BUFFER,lineBuf);gl.vertexAttribPointer(b,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.LINES,0,lineN);}
 if(showCtl){gl.uniform3fv(U2.col,hex(C.control.color));gl.bindBuffer(gl.ARRAY_BUFFER,ctlBuf);gl.vertexAttribPointer(b,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.LINES,0,ctlN);}
 gl.disableVertexAttribArray(b);
 // html labels
 let k=0;const want=(showGrid?gridLabels.length:0)+(showCtl?ctlLabels.length:0)+picks.length;ensureLabels(want);
 const put=(txt,p,cls)=>{const s=project(p,m);const d=labEls[k++];if(!s){d.style.display='none';return;}d.style.display='block';d.className='lab '+cls;d.textContent=txt;d.style.transform='translate('+s[0].toFixed(1)+'px,'+s[1].toFixed(1)+'px) translate(-50%,-50%)';};
 if(showGrid)for(const g of gridLabels)put(g.n,g.p,'g');
 if(showCtl)for(const g of ctlLabels)put(g.n,g.p,'c');
 picks.forEach((p,j)=>put(j===picks.length-1?'':'',p.w,'m'+(j===picks.length-1?'1':'0')));
 for(;k<labEls.length;k++)labEls[k].style.display='none';
 placePop(m);
}
let raf=0;function req(){if(!raf)raf=requestAnimationFrame(()=>{raf=0;draw();});}
// ---------------------------------------------------------------- input
let drag=null,pts=new Map(),moved=0;
cv.addEventListener('contextmenu',e=>e.preventDefault());
cv.addEventListener('pointerdown',e=>{cv.setPointerCapture(e.pointerId);pts.set(e.pointerId,{x:e.clientX,y:e.clientY});moved=0;
 drag={x:e.clientX,y:e.clientY,btn:e.button,shift:e.shiftKey,d0:pinchDist(),c0:pinchCtr()};});
cv.addEventListener('pointermove',e=>{
 if(!pts.has(e.pointerId)){hover(e);return;}
 const prev=pts.get(e.pointerId);pts.set(e.pointerId,{x:e.clientX,y:e.clientY});
 if(pts.size>=2){const d=pinchDist(),c=pinchCtr();if(drag&&drag.d0){cam.dist=Math.max(2,Math.min(2000,cam.dist*drag.d0/d));panBy(c.x-drag.c0.x,c.y-drag.c0.y);drag.d0=d;drag.c0=c;moved=99;req();}return;}
 const dx=e.clientX-prev.x,dy=e.clientY-prev.y;moved+=Math.abs(dx)+Math.abs(dy);if(!drag)return;
 if(drag.btn===2||drag.shift||drag.btn===1)panBy(dx,dy);else{cam.az-=dx*0.4;cam.el=Math.max(2,Math.min(89.9,cam.el+dy*0.4));if(cam.ortho&&cam.el<89){cam.ortho=false;syncUI();}}
 req();});
function pinchDist(){if(pts.size<2)return 0;const a=[...pts.values()];return Math.hypot(a[0].x-a[1].x,a[0].y-a[1].y);}
function pinchCtr(){const a=[...pts.values()];return a.length<2?{x:0,y:0}:{x:(a[0].x+a[1].x)/2,y:(a[0].y+a[1].y)/2};}
function panBy(dx,dy){const az=rad(cam.az),el=rad(cam.el);const s=cam.dist/cv.clientHeight*(cam.ortho?1:0.69);
 const rx=-Math.sin(az),ry=Math.cos(az);const fx=-Math.cos(az)*Math.sin(el),fy=-Math.sin(az)*Math.sin(el);
 cam.tx+=(-dx*rx+dy*fx)*s;cam.ty+=(-dx*ry+dy*fy)*s;}
cv.addEventListener('pointerup',e=>{const was=pts.size;pts.delete(e.pointerId);
 if(was===1&&moved<5&&e.button===0){const r=cv.getBoundingClientRect();const h=pick(e.clientX-r.left,e.clientY-r.top);if(h)showPick(h);else{picks=[];selEl=null;renderPop();}req();}
 if(pts.size===0)drag=null;});
cv.addEventListener('pointercancel',e=>{pts.delete(e.pointerId);drag=null;});
cv.addEventListener('wheel',e=>{e.preventDefault();cam.dist=Math.max(2,Math.min(2000,cam.dist*Math.exp(e.deltaY*0.0012)));req();},{passive:false});
cv.addEventListener('dblclick',e=>{const r=cv.getBoundingClientRect();const h=pick(e.clientX-r.left,e.clientY-r.top);if(h){cam.tx=h.w[0];cam.ty=h.w[1];cam.tz=h.w[2];req();}});
let hv=0;function hover(e){const now=performance.now();if(now-hv<70)return;hv=now;const r=cv.getBoundingClientRect();const h=pick(e.clientX-r.left,e.clientY-r.top);const el=$('hov');
 if(!h){el.textContent='E – , N – , Z –';return;}const i=w2info(h.w);el.textContent='E '+i.E.toFixed(3)+'   N '+i.N.toFixed(3)+'   Z '+i.Z.toFixed(3)+'   ·   '+h.el.n;}
// ---------------------------------------------------------------- panel
function buildPanel(){
 const L=$('layers');L.innerHTML='';
 for(const k in C.layers){const d=document.createElement('label');d.className='ly';d.innerHTML='<input type="checkbox" '+(vis[k]?'checked':'')+'><i style="background:'+C.layers[k].color+'"></i><span>'+C.layers[k].name+'</span>';d.querySelector('input').onchange=e=>{vis[k]=e.target.checked;if(selEl&&!vis[selEl.l]){selEl=null;}req();};L.appendChild(d);}
 const add=(label,color,get,set)=>{const d=document.createElement('label');d.className='ly';d.innerHTML='<input type="checkbox" '+(get()?'checked':'')+'><i style="background:'+color+'"></i><span>'+label+'</span>';d.querySelector('input').onchange=e=>{set(e.target.checked);req();};L.appendChild(d);};
 add('Grids 1–10 / A–F',C.grids.color,()=>showGrid,v=>showGrid=v);add('Control points P01, P21, P36, P49',C.control.color,()=>showCtl,v=>showCtl=v);
 const t=G[mode];const rows=G.points.map((p,i)=>{const r=G[mode].residuals[i];return'<tr><td>'+p.name+'</td><td>'+p.E.toFixed(3)+'</td><td>'+p.N.toFixed(3)+'</td><td>'+(Math.hypot(r.dE,r.dN)*100).toFixed(1)+'</td></tr>';}).join('');
 $('ctl').innerHTML='<table><tr><th></th><th>E</th><th>N</th><th>fit cm</th></tr>'+rows+'</table>';
 $('geo').textContent='Rotation '+t.theta_deg.toFixed(4)+'° · scale '+t.s.toFixed(6)+' · max residual '+(t.max_resid*100).toFixed(1)+' cm';
}
function syncUI(){$('sc').checked=mode==='similarity';$('ze').value=zEx;$('zev').textContent=zEx.toFixed(1)+'×';$('cut').checked=cutOn;$('cz').value=cutZ;$('czv').textContent=cutZ.toFixed(2);$('orth').checked=cam.ortho;buildPanel();req();}
$('sc').onchange=e=>{mode=e.target.checked?'similarity':'rigid';rebuild();picks=[];selEl=null;renderPop();buildPanel();req();};
$('ze').oninput=e=>{zEx=parseFloat(e.target.value);$('zev').textContent=zEx.toFixed(1)+'×';rebuild();picks=[];renderPop();req();};
$('cut').onchange=e=>{cutOn=e.target.checked;req();};$('cz').oninput=e=>{cutZ=parseFloat(e.target.value);$('czv').textContent=cutZ.toFixed(2);req();};
$('orth').onchange=e=>{cam.ortho=e.target.checked;req();};
$('b3').onclick=()=>setView('3d');$('bp').onclick=()=>setView('plan');$('bf').onclick=()=>{fit();req();};
$('go').onclick=()=>{const E=parseFloat($('gE').value),N=parseFloat($('gN').value);if(isNaN(E)||isNaN(N)){$('gmsg').textContent='Enter E and N';return;}
 const l=en2loc(E,N);const w=loc2w(l[0],l[1],0);const m={VP:null};
 // drop a vertical ray from above to find the surface height
 let best=null;const o=[w[0],w[1],50*zEx];const r={o,d:[0,0,-1]};let bt=1e12,bel=null;
 for(const e of EL){if(!vis[e.l])continue;if(!hitBox(r,e.mn,e.mx))continue;const P=e.pos;for(let t=0;t<e.nt;t++){const q=t*9;
  const ax=P[q],ay=P[q+1],bx=P[q+3],by=P[q+4],cx=P[q+6],cy=P[q+7];const d=(by-cy)*(ax-cx)+(cx-bx)*(ay-cy);if(Math.abs(d)<1e-12)continue;
  const l1=((by-cy)*(w[0]-cx)+(cx-bx)*(w[1]-cy))/d,l2=((cy-ay)*(w[0]-cx)+(ax-cx)*(w[1]-cy))/d,l3=1-l1-l2;if(l1<-1e-6||l2<-1e-6||l3<-1e-6)continue;
  const z=l1*P[q+2]+l2*P[q+5]+l3*P[q+8];const tt=o[2]-z;if(tt>0&&tt<bt){bt=tt;bel=e;best=z;}}}
 if(bel===null){$('gmsg').textContent='Point is outside the model footprint';cam.tx=w[0];cam.ty=w[1];req();return;}
 $('gmsg').textContent='';const h={w:[w[0],w[1],best],el:bel};picks=[];showPick(h);cam.tx=w[0];cam.ty=w[1];cam.tz=best;cam.dist=Math.min(cam.dist,120);req();};
$('tg').onclick=()=>$('panel').classList.toggle('hide');
$('team').textContent=C.team;document.title=C.team+' – '+C.title;$('ttl').textContent=C.title;$('sub').textContent=C.subtitle;

window.addEventListener('resize',req);window.addEventListener('keydown',e=>{if(e.key==='Escape'){picks=[];selEl=null;renderPop();req();}});
if(window.innerWidth<700)$('panel').classList.add('hide');
rebuild();fit();syncUI();
window.__viewer={cam,pick,w2info,loc2en,en2loc,req,draw,get EL(){return EL;},get mode(){return mode;},setMode(m){mode=m;rebuild();syncUI();},showPick,setZ(z){zEx=z;rebuild();}};
})();
