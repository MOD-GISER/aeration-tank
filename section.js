// Section tool: draw a line on the model -> true-scale vertical section (cut faces are computed exactly from the mesh).
(function(){
'use strict';
function init(){
const V=window.__viewer,D=V.D,C=V.C,$=id=>document.getElementById(id);
const EL=V.EL;
const rad=d=>d*Math.PI/180,deg=r=>r*180/Math.PI;
const PAD={l:62,b:38,t:44,r:18};
const S={A:null,B:null,u:[1,0],v:[0,1],L:0,depth:40,zx:1,k:1,s0:0,z0:0,open:false,cut:[],hov:null,pins:[],capA:null,capturing:false,visSig:''};
// ------------------------------------------------------------ DOM
const st=document.createElement('style');st.textContent=`
#secp{position:fixed;right:0;top:0;bottom:0;width:50%;display:none;flex-direction:column;background:var(--bg1);border-left:1px solid var(--line);z-index:4;color:var(--ink)}
body.sec #secp{display:flex}
body.sec #c,body.sec #labs{right:auto;width:50%}
body.sec #top{right:calc(50% + 12px)}
body.sec #panel{right:calc(50% + 12px)}
body.sec #pop{max-width:46vw}
#secp .sh{display:flex;flex-wrap:wrap;gap:6px;align-items:center;padding:8px 10px;border-bottom:1px solid var(--line);background:var(--panel)}
#secp .sh b#sect{font-size:13px}#secp .sh .sp{flex:1}
#secp .sh label{display:flex;align-items:center;gap:6px;font-size:12px;color:var(--mut)}
#secp .sh input[type=range]{width:90px}
#secp .sh button{padding:4px 10px;font-size:12px}
#secw{position:relative;flex:1;min-height:0;background:#fff;touch-action:none}
#secw canvas{position:absolute;left:0;top:0;width:100%;height:100%;display:block}
#seco{cursor:crosshair}
#secr{padding:7px 10px;border-top:1px solid var(--line);background:var(--panel);font-size:12px;font-variant-numeric:tabular-nums;min-height:34px;line-height:1.45}
#secr b{font-weight:650}
@media (max-width:700px){#secp{left:0;right:0;width:auto;top:42%}body.sec #c,body.sec #labs{width:100%;bottom:58%;height:42%}body.sec #top{right:8px}body.sec #panel{right:8px;max-height:30%}#secp .sh input[type=range]{width:60px}}
`;document.head.appendChild(st);
const root=document.createElement('div');root.id='secp';
root.innerHTML='<div class="sh"><b id="sect">Section</b><span class="sp"></span><label>Z× <select id="secz"><option value="1">1 (true)</option><option value="2">2</option><option value="5">5</option><option value="10">10</option></select></label><label>Depth <input type="range" id="secd" min="0" max="300" step="1"><b id="secdv"></b></label><button id="secfl" title="Look the other way">Flip</button><button id="secfit">Fit</button><button id="secpng">PNG</button><button id="secx" title="Close">×</button></div><div id="secw"><canvas id="secg"></canvas><canvas id="seco"></canvas></div><div id="secr"></div>';
document.body.appendChild(root);
const gc=$('secg'),oc=$('seco'),ctx=oc.getContext('2d');
// panel UI (left panel)
$('secui').innerHTML='<button id="secdraw" style="width:100%">Draw section: pick 2 points</button>'+
 '<div class="small" style="margin:6px 0 2px">or type the two end points (E, N)</div>'+
 '<div class="gr"><input id="sAE" placeholder="Start E" inputmode="decimal"><input id="sAN" placeholder="Start N" inputmode="decimal"></div>'+
 '<div class="gr" style="margin-top:6px"><input id="sBE" placeholder="End E" inputmode="decimal"><input id="sBN" placeholder="End N" inputmode="decimal"></div>'+
 '<button id="secgo" style="width:100%;margin-top:6px">Create section</button>'+
 '<div class="gr" style="margin-top:8px"><select id="secgrid" style="flex:1;min-width:0;font:inherit;color:var(--ink);background:transparent;border:1px solid var(--line);border-radius:8px;padding:6px"></select><button id="secgg">Along axis</button></div>'+
 '<div class="small" id="secmsg"></div>';
$('secgrid').innerHTML=D.grids.map((g,i)=>'<option value="'+i+'">Axis '+g.name+(g.dir==='v'?' (vertical)':' (horizontal)')+'</option>').join('');
const msg=t=>$('secmsg').textContent=t||'';
// ------------------------------------------------------------ geometry helpers
let zmin=1e9,zmax=-1e9;for(const e of EL){for(let k=2;k<e.v.length;k+=3){if(e.v[k]<zmin)zmin=e.v[k];if(e.v[k]>zmax)zmax=e.v[k];}
 let a=[1e9,1e9,1e9],b=[-1e9,-1e9,-1e9];for(let k=0;k<e.v.length;k+=3)for(let d=0;d<3;d++){if(e.v[k+d]<a[d])a[d]=e.v[k+d];if(e.v[k+d]>b[d])b[d]=e.v[k+d];}e.lb=[a,b];}
function setLine(A,B){const dx=B[0]-A[0],dy=B[1]-A[1],L=Math.hypot(dx,dy);S.A=A;S.B=B;S.L=L;S.u=[dx/L,dy/L];S.v=[-dy/L,dx/L];}
function bearing(vx,vy){const o=V.loc2en(0,0),p=V.loc2en(vx,vy);let b=deg(Math.atan2(p[0]-o[0],p[1]-o[1]));if(b<0)b+=360;
 const c=['N','NE','E','SE','S','SW','W','NW'][Math.round(b/45)%8];return b.toFixed(1)+'° ('+c+')';}
// ------------------------------------------------------------ exact cut of the meshes by the vertical plane
function cutAll(){const out=[];const A=S.A,u=S.u,v=S.v;
 for(const e of EL){if(!V.vis[e.l])continue;
  const a=e.lb[0],b=e.lb[1];let dmin=1e9,dmax=-1e9;for(const x of [a[0],b[0]])for(const y of [a[1],b[1]]){const d=(x-A[0])*v[0]+(y-A[1])*v[1];if(d<dmin)dmin=d;if(d>dmax)dmax=d;}
  if(dmax<0||dmin>=0)continue;
  const P=e.v,I=e.i;let f0=0;
  for(let pi=0;pi<e.pr.length;pi++){const nf=e.pr[pi];const segs=[];
   for(let t=f0;t<f0+nf;t++){const ix=[I[t*3]*3,I[t*3+1]*3,I[t*3+2]*3];const d=[0,0,0],ps=[0,0,0];
    for(let q=0;q<3;q++){d[q]=(P[ix[q]]-A[0])*v[0]+(P[ix[q]+1]-A[1])*v[1];ps[q]=d[q]>=0?1:0;}
    if(ps[0]===ps[1]&&ps[1]===ps[2])continue;
    const pts=[];
    for(let q=0;q<3;q++){const q2=(q+1)%3;if(ps[q]===ps[q2])continue;let i0=ix[q],i1=ix[q2],d0=d[q],d1=d[q2];
     // canonical order so shared edges give bit-identical points
     if(P[i0]>P[i1]||(P[i0]===P[i1]&&(P[i0+1]>P[i1+1]||(P[i0+1]===P[i1+1]&&P[i0+2]>P[i1+2])))){let tmp=i0;i0=i1;i1=tmp;tmp=d0;d0=d1;d1=tmp;}
     const tt=d0/(d0-d1);const x=P[i0]+tt*(P[i1]-P[i0]),y=P[i0+1]+tt*(P[i1+1]-P[i0+1]),z=P[i0+2]+tt*(P[i1+2]-P[i0+2]);
     pts.push([(x-A[0])*u[0]+(y-A[1])*u[1],z]);}
    if(pts.length===2)segs.push(pts);}
   f0+=nf;
   // chain segments into loops
   const key=p=>Math.round(p[0]*1e4)+'_'+Math.round(p[1]*1e4),adj=new Map();
   segs.forEach((sg,i)=>{for(const p of sg){const k=key(p);if(!adj.has(k))adj.set(k,[]);adj.get(k).push(i);}});
   const used=new Uint8Array(segs.length),loops=[];
   for(let i=0;i<segs.length;i++){if(used[i])continue;used[i]=1;const loop=[segs[i][0],segs[i][1]];const k0=key(segs[i][0]);let cur=key(segs[i][1]);
    while(cur!==k0){const nx=(adj.get(cur)||[]).find(j=>!used[j]);if(nx===undefined)break;used[nx]=1;const sg=segs[nx];const first=key(sg[0])===cur;loop.push(first?sg[1]:sg[0]);cur=key(first?sg[1]:sg[0]);}
    if(loop.length>=3)loops.push(loop);}
   if(loops.length)out.push({e,loops});}
 }
 return out;}
// ------------------------------------------------------------ GL (beyond view)
const gl=gc.getContext('webgl',{antialias:true,alpha:false,preserveDrawingBuffer:false});
let GP=null,UL=null,built=false;
function sh(t,src){const s=gl.createShader(t);gl.shaderSource(s,src);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw new Error(gl.getShaderInfoLog(s));return s;}
function glInit(){if(GP||!gl)return;
 GP=gl.createProgram();
 gl.attachShader(GP,sh(gl.VERTEX_SHADER,'attribute vec3 p;attribute vec3 n;uniform vec2 A;uniform vec2 U;uniform vec2 Vd;uniform vec4 M;uniform vec4 R;uniform float Zx;varying vec3 vn;void main(){vec2 q=p.xy-A;float s=dot(q,U);float d=dot(q,Vd);float px=(s-M.x)*M.z+R.x;float py=R.w-((p.z-M.y)*M.z*Zx+R.y);gl_Position=vec4(px/R.z*2.0-1.0,1.0-py/R.w*2.0,d/M.w*2.0-1.0,1.0);vn=vec3(dot(n.xy,Vd),0.0,n.z);}'));
 gl.attachShader(GP,sh(gl.FRAGMENT_SHADER,'precision mediump float;varying vec3 vn;uniform vec3 col;uniform float isline;void main(){if(isline>0.5){gl_FragColor=vec4(col,1.0);return;}vec3 n=normalize(vn);float l=0.58+0.30*abs(n.x)+0.12*abs(n.z);vec3 c=col*l;c=mix(c,vec3(1.0),0.28);gl_FragColor=vec4(c,1.0);}'));
 gl.linkProgram(GP);if(!gl.getProgramParameter(GP,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(GP));
 UL={};for(const k of ["A","U","Vd","M","R","Zx","col","isline"])UL[k]=gl.getUniformLocation(GP,k);}
function hex(h){return[parseInt(h.substr(1,2),16)/255,parseInt(h.substr(3,2),16)/255,parseInt(h.substr(5,2),16)/255];}
function buildBuffers(){if(built||!gl)return;glInit();
 for(const e of EL){const nt=e.i.length/3,pos=new Float32Array(nt*9),nor=new Float32Array(nt*9);
  for(let t=0;t<nt;t++){const ia=e.i[t*3]*3,ib=e.i[t*3+1]*3,ic=e.i[t*3+2]*3,idx=[ia,ib,ic];
   const ux=e.v[ib]-e.v[ia],uy=e.v[ib+1]-e.v[ia+1],uz=e.v[ib+2]-e.v[ia+2],vx=e.v[ic]-e.v[ia],vy=e.v[ic+1]-e.v[ia+1],vz=e.v[ic+2]-e.v[ia+2];
   let nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;const l=Math.hypot(nx,ny,nz)||1;nx/=l;ny/=l;nz/=l;
   for(let q=0;q<3;q++)for(let d=0;d<3;d++){pos[t*9+q*3+d]=e.v[idx[q]+d];nor[t*9+q*3+d]=d===0?nx:d===1?ny:nz;}}
  e.sb={p:gl.createBuffer(),n:gl.createBuffer(),e:gl.createBuffer(),nt,ne:e.ed.length/3};
  gl.bindBuffer(gl.ARRAY_BUFFER,e.sb.p);gl.bufferData(gl.ARRAY_BUFFER,pos,gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER,e.sb.n);gl.bufferData(gl.ARRAY_BUFFER,nor,gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER,e.sb.e);gl.bufferData(gl.ARRAY_BUFFER,e.ed,gl.STATIC_DRAW);}
 built=true;}
// ------------------------------------------------------------ view mapping
let W=0,H=0,dpr=1;
const sx=s=>PAD.l+(s-S.s0)*S.k,sy=z=>H-PAD.b-(z-S.z0)*S.k*S.zx,inv=(px,py)=>[(px-PAD.l)/S.k+S.s0,(H-PAD.b-py)/(S.k*S.zx)+S.z0];
function size(){const r=$('secw').getBoundingClientRect();dpr=Math.min(window.devicePixelRatio||1,2);W=Math.max(100,r.width);H=Math.max(100,r.height);
 for(const c of [gc,oc]){const w=Math.round(W*dpr),h=Math.round(H*dpr);if(c.width!==w||c.height!==h){c.width=w;c.height=h;}}}
function fitView(){size();const aw=W-PAD.l-PAD.r,ah=H-PAD.t-PAD.b,zr=zmax-zmin+0.8;S.k=Math.min(aw/S.L,ah/(zr*S.zx));
 S.s0=-(aw/S.k-S.L)/2;S.z0=zmin-0.4-(ah/(S.k*S.zx)-zr)/2;}
// ------------------------------------------------------------ drawing
function niceStep(pxPerM,minPx){const c=[0.05,0.1,0.25,0.5,1,2,5,10,20,25,50,100];for(const s of c)if(s*pxPerM>=minPx)return s;return 100;}
function drawGL(){if(!gl)return;buildBuffers();gl.viewport(0,0,gc.width,gc.height);gl.clearColor(1,1,1,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
 if(S.depth<=0)return;gl.enable(gl.DEPTH_TEST);gl.useProgram(GP);
 gl.uniform2f(UL.A,S.A[0],S.A[1]);gl.uniform2f(UL.U,S.u[0],S.u[1]);gl.uniform2f(UL.Vd,S.v[0],S.v[1]);gl.uniform4f(UL.M,S.s0,S.z0,S.k,S.depth);gl.uniform4f(UL.R,PAD.l,PAD.b,W,H);gl.uniform1f(UL.Zx,S.zx);
 const a=gl.getAttribLocation(GP,'p'),n=gl.getAttribLocation(GP,'n');gl.enableVertexAttribArray(a);
 gl.enable(gl.POLYGON_OFFSET_FILL);gl.polygonOffset(1,1);
 gl.enableVertexAttribArray(n);gl.uniform1f(UL.isline,0);
 for(const e of EL){if(!V.vis[e.l])continue;gl.uniform3fv(UL.col,hex(C.layers[e.l].color));
  gl.bindBuffer(gl.ARRAY_BUFFER,e.sb.p);gl.vertexAttribPointer(a,3,gl.FLOAT,false,0,0);gl.bindBuffer(gl.ARRAY_BUFFER,e.sb.n);gl.vertexAttribPointer(n,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.TRIANGLES,0,e.sb.nt*3);}
 gl.disable(gl.POLYGON_OFFSET_FILL);gl.disableVertexAttribArray(n);gl.vertexAttrib3f(n,0,0,1);gl.uniform1f(UL.isline,1);gl.uniform3f(UL.col,0.25,0.29,0.33);
 for(const e of EL){if(!V.vis[e.l]||!e.sb.ne)continue;gl.bindBuffer(gl.ARRAY_BUFFER,e.sb.e);gl.vertexAttribPointer(a,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.LINES,0,e.sb.ne);}
 gl.disableVertexAttribArray(a);}
let hatch=null;
function hatchPat(){if(hatch)return hatch;const c=document.createElement('canvas');c.width=c.height=8;const g=c.getContext('2d');g.strokeStyle='rgba(255,255,255,.55)';g.lineWidth=1;g.beginPath();g.moveTo(-1,9);g.lineTo(9,-1);g.moveTo(-1,1);g.lineTo(1,-1);g.moveTo(7,9);g.lineTo(9,7);g.stroke();hatch=ctx.createPattern(c,'repeat');return hatch;}
function dark(h,f){const c=hex(h);return'rgb('+c.map(x=>Math.round(x*255*f)).join(',')+')';}
function drawOverlay(){ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,W,H);
 const ink='#1a2733',mut='#5b6b7a',grid='#e3e8ee';ctx.font='11px system-ui,sans-serif';
 // ruler grid
 const zs=niceStep(S.k*S.zx,46),ss=niceStep(S.k,70);ctx.lineWidth=1;ctx.strokeStyle=grid;ctx.fillStyle=mut;ctx.textBaseline='middle';ctx.textAlign='right';
 const zlo=S.z0-0.01,zhi=S.z0+(H-PAD.b)/(S.k*S.zx);
 for(let z=Math.ceil(zlo/zs)*zs;z<=zhi;z+=zs){const y=Math.round(sy(z))+0.5;if(y<PAD.t-4||y>H-PAD.b)continue;ctx.beginPath();ctx.moveTo(PAD.l,y);ctx.lineTo(W-PAD.r,y);ctx.stroke();ctx.fillText((z>=0?'+':'')+z.toFixed(zs<0.1?2:zs<1?2:1),PAD.l-6,y);}
 ctx.textAlign='center';ctx.textBaseline='top';
 const shi=S.s0+(W-PAD.l)/S.k;
 for(let s=Math.ceil(S.s0/ss)*ss;s<=shi;s+=ss){const x=Math.round(sx(s))+0.5;if(x<PAD.l||x>W-PAD.r)continue;ctx.strokeStyle=grid;ctx.beginPath();ctx.moveTo(x,PAD.t);ctx.lineTo(x,H-PAD.b);ctx.stroke();ctx.fillText(s.toFixed(ss<1?2:ss<5?1:0),x,H-PAD.b+6);}
 // clip drawing area
 ctx.save();ctx.beginPath();ctx.rect(PAD.l,PAD.t,W-PAD.l-PAD.r,H-PAD.t-PAD.b);ctx.clip();
 // level +1.15
 {const y=sy(1.15);ctx.setLineDash([6,4]);ctx.strokeStyle='#8a6d00';ctx.beginPath();ctx.moveTo(PAD.l,y);ctx.lineTo(W-PAD.r,y);ctx.stroke();ctx.setLineDash([]);ctx.fillStyle='#8a6d00';ctx.textAlign='right';ctx.textBaseline='bottom';ctx.fillText('+1.15 raft top',W-PAD.r-4,y-2);}
 // grid lines crossing the section
 for(const g of D.grids){const A=S.A,B=S.B;let t,y0,y1,s;
  if(g.dir==='v'){if(Math.abs(B[0]-A[0])<1e-9)continue;t=(g.pos-A[0])/(B[0]-A[0]);const yy=A[1]+t*(B[1]-A[1]);if(t<0||t>1||yy<g.a||yy>g.b)continue;}
  else{if(Math.abs(B[1]-A[1])<1e-9)continue;t=(g.pos-A[1])/(B[1]-A[1]);const xx=A[0]+t*(B[0]-A[0]);if(t<0||t>1||xx<g.a||xx>g.b)continue;}
  s=t*S.L;const x=sx(s);ctx.setLineDash([10,4,2,4]);ctx.strokeStyle='#d8334a';ctx.beginPath();ctx.moveTo(x,PAD.t);ctx.lineTo(x,H-PAD.b);ctx.stroke();ctx.setLineDash([]);
  ctx.fillStyle='#d8334a';ctx.beginPath();ctx.arc(x,PAD.t+9,9,0,7);ctx.fill();ctx.fillStyle='#fff';ctx.textAlign='center';ctx.textBaseline='middle';ctx.font='600 11px system-ui,sans-serif';ctx.fillText(g.name,x,PAD.t+9);ctx.font='11px system-ui,sans-serif';}
 // cut faces
 ctx.lineJoin='round';
 for(const c of S.cut){const col=C.layers[c.e.l].color;const path=new Path2D();
  for(const lp of c.loops){lp.forEach((p,i)=>{const x=sx(p[0]),y=sy(p[1]);if(i)path.lineTo(x,y);else path.moveTo(x,y);});path.closePath();}
  ctx.fillStyle=dark(col,0.62);ctx.fill(path,'evenodd');ctx.save();ctx.clip(path,'evenodd');ctx.fillStyle=hatchPat();ctx.fillRect(0,0,W,H);ctx.restore();
  ctx.strokeStyle='#0d141b';ctx.lineWidth=1.6;ctx.stroke(path);}
 // pins + hover
 const mark=(p,col)=>{const x=sx(p.s),y=sy(p.z);ctx.strokeStyle=col;ctx.lineWidth=2;
  if(p.type==='Endpoint'){ctx.strokeRect(x-6,y-6,12,12);}else if(p.type==='Midpoint'){ctx.beginPath();ctx.moveTo(x,y-7);ctx.lineTo(x-7,y+6);ctx.lineTo(x+7,y+6);ctx.closePath();ctx.stroke();}
  else if(p.type==='On line'){ctx.beginPath();ctx.moveTo(x-6,y-6);ctx.lineTo(x+6,y+6);ctx.moveTo(x+6,y-6);ctx.lineTo(x-6,y+6);ctx.stroke();}
  else{ctx.beginPath();ctx.arc(x,y,4,0,7);ctx.stroke();}};
 if(S.pins.length===2){const a=S.pins[0],b=S.pins[1];ctx.strokeStyle='#1f6feb';ctx.lineWidth=1.5;ctx.setLineDash([5,3]);ctx.beginPath();ctx.moveTo(sx(a.s),sy(a.z));ctx.lineTo(sx(b.s),sy(a.z));ctx.lineTo(sx(b.s),sy(b.z));ctx.stroke();ctx.setLineDash([]);
  ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(sx(a.s),sy(a.z));ctx.lineTo(sx(b.s),sy(b.z));ctx.stroke();}
 S.pins.forEach(p=>mark(p,'#1f6feb'));if(S.hov)mark(S.hov,'#17b26a');
 ctx.restore();
 // frame + title
 ctx.strokeStyle='#9aa8b5';ctx.lineWidth=1;ctx.strokeRect(PAD.l+0.5,PAD.t+0.5,W-PAD.l-PAD.r,H-PAD.t-PAD.b);
 ctx.fillStyle=ink;ctx.textAlign='left';ctx.textBaseline='top';ctx.font='600 13px system-ui,sans-serif';ctx.fillText('Section A – B',PAD.l,8);
 ctx.font='11px system-ui,sans-serif';ctx.fillStyle=mut;const eA=V.loc2en(S.A[0],S.A[1]),eB=V.loc2en(S.B[0],S.B[1]);
 ctx.fillText('A  E '+eA[0].toFixed(3)+'  N '+eA[1].toFixed(3)+'      B  E '+eB[0].toFixed(3)+'  N '+eB[1].toFixed(3),PAD.l,24);
 ctx.textAlign='right';ctx.fillText('A',PAD.l+2,PAD.t-10);ctx.fillText('B',W-PAD.r,PAD.t-10);
 ctx.textAlign='left';ctx.fillText('Level (m)',4,PAD.t-10);if(S.zx!==1){ctx.fillStyle='#c0392b';ctx.font='600 12px system-ui,sans-serif';ctx.fillText('Vertical scale ×'+S.zx+' (not true scale)',PAD.l+170,8);}
}
let raf=0;function redraw(){if(raf)return;raf=requestAnimationFrame(()=>{raf=0;if(!S.open)return;drawGL();drawOverlay();});}
// ------------------------------------------------------------ snapping inside the section view
function candidates(px,py,tol){const out=[];const m=[S.s0,S.z0];
 const P=(s,z)=>[sx(s),sy(z)];
 for(const c of S.cut)for(const lp of c.loops){const n=lp.length;
  for(let i=0;i<n;i++){const a=lp[i],b=lp[(i+1)%n];const A=P(a[0],a[1]),B=P(b[0],b[1]);
   if(Math.max(A[0],B[0])<px-tol||Math.min(A[0],B[0])>px+tol||Math.max(A[1],B[1])<py-tol||Math.min(A[1],B[1])>py+tol)continue;
   const dA=Math.hypot(A[0]-px,A[1]-py);if(dA<tol)out.push({k:0,d:dA,s:a[0],z:a[1],type:'Endpoint',el:c.e});
   const mx=(A[0]+B[0])/2,my=(A[1]+B[1])/2,dM=Math.hypot(mx-px,my-py);if(dM<tol)out.push({k:1,d:dM,s:(a[0]+b[0])/2,z:(a[1]+b[1])/2,type:'Midpoint',el:c.e});
   const dx=B[0]-A[0],dy=B[1]-A[1],L2=dx*dx+dy*dy;if(L2>1e-9){let t=((px-A[0])*dx+(py-A[1])*dy)/L2;t=Math.max(0,Math.min(1,t));const d=Math.hypot(A[0]+t*dx-px,A[1]+t*dy-py);
    if(d<tol)out.push({k:2,d,s:a[0]+(b[0]-a[0])*t,z:a[1]+(b[1]-a[1])*t,type:'On line',el:c.e});}}}
 return out;}
function inside(c,s,z){let ins=false;for(const lp of c.loops){const n=lp.length;let w=false;for(let i=0,j=n-1;i<n;j=i++){const a=lp[i],b=lp[j];if((a[1]>z)!==(b[1]>z)&&s<(b[0]-a[0])*(z-a[1])/(b[1]-a[1])+a[0])w=!w;}if(w)ins=!ins;}return ins;}
function pointAt(px,py,touch){const tol=touch?24:12;const cs=candidates(px,py,tol);let best=null;
 if(V.snapOn&&cs.length){cs.sort((a,b)=>a.k-b.k||a.d-b.d);best=cs[0];}
 const g=inv(px,py);if(!best){let el=null;for(const c of S.cut)if(inside(c,g[0],g[1])){el=c.e;break;}best={s:g[0],z:g[1],type:'Free',el};}
 return best;}
function fmtPoint(p){const x=S.A[0]+S.u[0]*p.s,y=S.A[1]+S.u[1]*p.s,en=V.loc2en(x,y);return{E:en[0],N:en[1]};}
function readout(){const p=S.hov;let h='';
 if(p){const q=fmtPoint(p);h+='Station <b>'+p.s.toFixed(3)+'</b> · Z <b>'+(p.z>=0?'+':'')+p.z.toFixed(3)+'</b> · E <b>'+q.E.toFixed(3)+'</b> · N <b>'+q.N.toFixed(3)+'</b>'+(p.type!=='Free'?' · '+p.type:'')+(p.el?' · '+p.el.n:'');}
 else h+='Move over the section. Snap works on the cut outline.';
 if(S.pins.length===2){const a=S.pins[0],b=S.pins[1],ds=b.s-a.s,dz=b.z-a.z,d=Math.hypot(ds,dz);h+='<br>Dimension: horizontal <b>'+Math.abs(ds).toFixed(3)+'</b> · vertical <b>'+Math.abs(dz).toFixed(3)+'</b> · distance <b>'+d.toFixed(3)+'</b> m'+(Math.abs(ds)>1e-6?' · slope <b>'+(dz/ds*100).toFixed(2)+'%</b> (1:'+(Math.abs(dz)>1e-9?Math.abs(ds/dz).toFixed(2):'∞')+')':'');}
 else if(S.pins.length===1)h+='<br>Click a second point to measure.';
 $('secr').innerHTML=h;}
// ------------------------------------------------------------ open / close
function recompute(){S.cut=cutAll();S.visSig=JSON.stringify(V.vis);}
function openSection(A,B){if(Math.hypot(B[0]-A[0],B[1]-A[1])<0.05){msg('The two points are too close.');return;}
 setLine(A,B);S.pins=[];S.hov=null;S.open=true;document.body.classList.add('sec');window.dispatchEvent(new Event('resize'));
 V.setSecLine(A,B);recompute();$('secd').value=S.depth;$('secdv').textContent=S.depth+' m';
 $('sect').textContent='Section · looking '+bearing(S.v[0],S.v[1])+' · length '+S.L.toFixed(3)+' m';
 setTimeout(()=>{fitView();readout();redraw();},30);msg('');}
function closeSection(){S.open=false;document.body.classList.remove('sec');V.setSecLine(null,null);window.dispatchEvent(new Event('resize'));}
$('secx').onclick=closeSection;
$('secfit').onclick=()=>{fitView();redraw();};
$('secfl').onclick=()=>{const A=S.B,B=S.A;openSection(A,B);};
$('secz').onchange=e=>{S.zx=parseFloat(e.target.value);fitView();redraw();};
$('secd').oninput=e=>{S.depth=parseFloat(e.target.value);$('secdv').textContent=S.depth+' m';redraw();};
$('secpng').onclick=()=>{drawGL();const c=document.createElement('canvas');c.width=gc.width;c.height=gc.height;const g=c.getContext('2d');g.fillStyle='#fff';g.fillRect(0,0,c.width,c.height);g.drawImage(gc,0,0);g.drawImage(oc,0,0);
 c.toBlob(b=>{if(!b)return;const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='section.png';document.body.appendChild(a);a.click();setTimeout(()=>{URL.revokeObjectURL(a.href);a.remove();},500);});};
// pan / zoom / hover / pin
const pts=new Map();let mv=0,last=null;
oc.addEventListener('pointerdown',e=>{oc.setPointerCapture(e.pointerId);pts.set(e.pointerId,{x:e.offsetX,y:e.offsetY});mv=0;last=null;});
oc.addEventListener('pointermove',e=>{const x=e.offsetX,y=e.offsetY;
 if(!pts.has(e.pointerId)){S.hov=pointAt(x,y,false);readout();redraw();return;}
 const prev=pts.get(e.pointerId);pts.set(e.pointerId,{x,y});
 if(pts.size>=2){const a=[...pts.values()];const d=Math.hypot(a[0].x-a[1].x,a[0].y-a[1].y),cx=(a[0].x+a[1].x)/2,cy=(a[0].y+a[1].y)/2;
  if(last){zoomAt(cx,cy,d/last.d);S.s0-=(cx-last.cx)/S.k;S.z0+=(cy-last.cy)/(S.k*S.zx);}last={d,cx,cy};mv=99;redraw();return;}
 const dx=x-prev.x,dy=y-prev.y;mv+=Math.abs(dx)+Math.abs(dy);if(mv>4){S.s0-=dx/S.k;S.z0+=dy/(S.k*S.zx);S.hov=null;redraw();}});
oc.addEventListener('pointerup',e=>{const was=pts.size;pts.delete(e.pointerId);if(pts.size<2)last=null;
 if(was===1&&mv<5){const p=pointAt(e.offsetX,e.offsetY,e.pointerType==='touch');S.hov=p;if(S.pins.length>=2)S.pins=[];S.pins.push(p);readout();redraw();}});
oc.addEventListener('pointercancel',e=>{pts.delete(e.pointerId);last=null;});
function zoomAt(px,py,f){const g=inv(px,py);S.k=Math.max(0.2,Math.min(4000,S.k*f));S.s0=g[0]-(px-PAD.l)/S.k;S.z0=g[1]-(H-PAD.b-py)/(S.k*S.zx);}
oc.addEventListener('wheel',e=>{e.preventDefault();zoomAt(e.offsetX,e.offsetY,Math.exp(-e.deltaY*0.0015));redraw();},{passive:false});
oc.addEventListener('dblclick',()=>{fitView();redraw();});
oc.addEventListener('pointerleave',()=>{if(pts.size===0){S.hov=null;readout();redraw();}});
let rzT=0;window.addEventListener('resize',()=>{if(!S.open)return;clearTimeout(rzT);rzT=setTimeout(()=>{fitView();readout();redraw();},120);});
// redraw when layers change in the main view
setInterval(()=>{if(!S.open)return;const sig=JSON.stringify(V.vis);if(sig!==S.visSig){recompute();readout();redraw();}},300);
// ------------------------------------------------------------ creating the line
function stopCapture(){S.capturing=false;S.capA=null;V.hooks.capture=null;V.hooks.move=null;$('bsec').classList.remove('on');$('secdraw').textContent='Draw section: pick 2 points';V.setSecLine(S.open?S.A:null,S.open?S.B:null);}
function startCapture(){S.capturing=true;S.capA=null;$('bsec').classList.add('on');$('secdraw').textContent='Pick the start point… (tap again to cancel)';$('panel').classList.contains('hide')||0;
 V.hooks.capture=info=>{const p=[info.i.lx,info.i.ly];
  if(!S.capA){S.capA=p;$('secdraw').textContent='Pick the end point…';V.setSecLine(p,null,p);return;}
  const A=S.capA;stopCapture();openSection(A,p);};
 V.hooks.move=info=>{if(S.capA&&info)V.setSecLine(S.capA,null,[info.i.lx,info.i.ly]);};
 msg('Tap the model for the start point, then for the end point. Snap picks exact corners.');}
function toggleCapture(){if(S.capturing){stopCapture();msg('');}else startCapture();}
$('secdraw').onclick=toggleCapture;$('bsec').onclick=toggleCapture;
$('secgo').onclick=()=>{const v=['sAE','sAN','sBE','sBN'].map(i=>parseFloat($(i).value));if(v.some(isNaN)){msg('Enter all four numbers (E, N for both points).');return;}
 stopCapture();openSection(V.en2loc(v[0],v[1]),V.en2loc(v[2],v[3]));};
$('secgg').onclick=()=>{const g=D.grids[parseInt($('secgrid').value)];stopCapture();
 if(g.dir==='v')openSection([g.pos,g.a],[g.pos,g.b]);else openSection([g.a,g.pos],[g.b,g.pos]);};
window.__sec={S,cutAll,openSection,recompute,redraw,fitView,size};
}
if(window.__viewer)init();else document.addEventListener('tankready',init);
})();
