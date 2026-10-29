import {simulate,defaults,model,derivative,spectrum} from './dist/engine.js';
import fs from 'node:fs';
const assert=(v,msg)=>{if(!v)throw Error(msg);};
let r=simulate(),half=simulate({dt:.25}),zero=simulate({drive:0}),closed=simulate({gamma:0,relax:0}),tests=[];
let metric=(name,value,tolerance)=>{assert(value<tolerance,name+': '+value);tests.push({name,value,tolerance,passed:true});};
let d=r.m.d,L=d*d;
metric('Trace preservation',Math.max(...r.frames.map(f=>Math.abs(Array.from({length:d},(_,i)=>f.rho[i*d+i]).reduce((a,b)=>a+b,0)-1))),1e-10);
metric('Hermiticity',Math.max(...r.frames.flatMap(f=>Array.from({length:L},(_,k)=>{let i=Math.floor(k/d),j=k%d;return Math.max(Math.abs(f.rho[k]-f.rho[j*d+i]),Math.abs(f.rho[L+k]+f.rho[L+j*d+i]));}))),1e-10);
metric('Half-step convergence, density matrix',Math.max(...r.frames.flatMap((f,i)=>f.rho.map((v,k)=>Math.abs(v-half.frames[i].rho[k])))),2e-6);
metric('Zero-field ground-state stationarity',Math.max(...zero.frames.map(f=>Math.abs(f.excited))),1e-12);
metric('Closed-system purity error',Math.max(...closed.frames.map(f=>Math.abs(f.purity-1))),2e-6);
let m=model(),sp=spectrum(m),k=0;metric('Oscillator-strength sum rule',Math.abs(sp.lines.reduce((s,l)=>s+l.strength,0)-m.g.reduce((s,g)=>s+g.mu[k]**2,0)),1e-10);
// Analytic single-site decay and optical dephasing: exact generator check.
let one=model({n:2,j:0,drive:0}),y=new Float64Array(18);y[4]=1;let dy=derivative(one,0,y);metric('Analytic excited-state decay generator',Math.abs(dy[4]+defaults.relax)+Math.abs(dy[0]-defaults.relax),1e-14);
fs.writeFileSync('dist/validation.json',JSON.stringify({tests,model:'helix-v0.1',status:'passed'},null,2));fs.writeFileSync('/tmp/helix-reference-input.json',JSON.stringify({params:r.m.p,H:Array.from(r.m.H),g:r.m.g,frames:r.frames}));console.log(JSON.stringify(tests,null,2));
