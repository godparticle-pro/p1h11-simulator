import {solve,geometryScan,compare} from './research-engine.js';
onmessage=({data})=>{try{if(data.type==='compare'){let a=solve({...data.p,pol:'x'}),b=solve({...data.p,pol:'y'});postMessage({type:'compare',data:compare(a,b)});}else {let run=solve(data.p);postMessage({type:'run',run,scan:geometryScan(run.m.p)});}}catch(e){postMessage({type:'error',message:e.message});}};
