import json,numpy as np
from scipy.integrate import solve_ivp
runs=json.load(open('/tmp/p1h11-reference.json'));error=0.;min_eigen=0.
for a in runs:
 p=a['p'];d=p['n']+1;L=d*d;H0=np.array(a['H']).reshape(d,d).astype(complex);mu=np.array([g['mu'] for g in a['g']]);r=np.array([g['r'] for g in a['g']]);e={'x':[1,0,0],'y':[0,1,0],'plus':[1/np.sqrt(2),1j/np.sqrt(2),0],'minus':[1/np.sqrt(2),-1j/np.sqrt(2),0]}[p['pol']]
 x=r[:,0]-p['offset']*p['waist'];rr=np.hypot(x,r[:,1]);ell=p['ell'];l=abs(ell);amp=np.exp(-(rr/p['waist'])**2)*(np.sqrt(2)*rr/p['waist'])**l/np.sqrt(2 if l==2 else 1);phase=ell*np.arctan2(r[:,1],x)+2*np.pi*r[:,2]/(1239.841984/p['carrier']);v=-p['drive']*(mu@np.array(e))*amp*np.exp(1j*phase)
 mask=(np.arange(d)>0).astype(float);rate=(mask[:,None]+mask[None,:])/2;deph=p['gamma']*rate;np.fill_diagonal(deph,0)
 def rhs(t,y):
  R=y.reshape(d,d);H=H0.copy();pulse=np.exp(-.5*((t-p['center'])/p['sigma'])**2);H[1:,0]=v*pulse;H[0,1:]=np.conj(v)*pulse;out=-1j*(H@R-R@H)/.6582119569-deph*R-p['relax']*rate*R;out[0,0]+=p['relax']*np.trace(R[1:,1:]);return out.ravel()
 ts=[f['t'] for f in a['frames']];y0=np.zeros(d*d,complex);y0[0]=1;sol=solve_ivp(rhs,[0,p['end']],y0,t_eval=ts,method='DOP853',rtol=1e-10,atol=1e-12)
 assert sol.success
 raw=np.array([f['rho'] for f in a['frames']]);R=(raw[:,:L]+1j*raw[:,L:]).reshape(-1,d,d);error=max(error,float(np.max(np.abs(R-sol.y.T.reshape(-1,d,d)))));min_eigen=min(min_eigen,float(np.linalg.eigvalsh(R).min()))
assert error<2e-6,error
assert min_eigen>-2e-8,min_eigen
path='dist/research-validation.json';report=json.load(open(path));report['tests'] += [{'name':'Independent SciPy DOP853 across 10 configurations','value':error,'tolerance':2e-6,'passed':True},{'name':'Positive semidefinite density matrices across 10 configurations','value':max(0,-min_eigen),'tolerance':2e-8,'passed':True}];json.dump(report,open(path,'w'),indent=2);print('Independent maximum error:',error,'; minimum rho eigenvalue:',min_eigen)
