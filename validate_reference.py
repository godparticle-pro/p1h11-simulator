import json,numpy as np
from scipy.integrate import solve_ivp
p='dist/validation.json'
a=json.load(open('/tmp/helix-reference-input.json'));q=a['params'];n=q['n'];d=n+1;H0=np.array(a['H']).reshape(d,d);mu=np.array([x['mu'] for x in a['g']]);hbar=.6582119569
rho0=np.zeros((d,d),complex);rho0[0,0]=1
mask=(np.arange(d)>0).astype(float);rate=(mask[:,None]+mask[None,:])/2
D=q['gamma']*rate;np.fill_diagonal(D,0)
def rhs(t,y):
 r=y.reshape(d,d);H=H0.copy();v=-q['drive']*np.exp(-.5*((t-q['center'])/q['sigma'])**2)*mu[:,0];H[0,1:]=v;H[1:,0]=v
 out=-1j*(H@r-r@H)/hbar-D*r-q['relax']*rate*r;out[0,0]+=q['relax']*np.trace(r[1:,1:]);return out.ravel()
t=np.array([f['t'] for f in a['frames']]);sol=solve_ivp(rhs,[0,q['end']],rho0.ravel(),t_eval=t,rtol=1e-10,atol=1e-12,method='DOP853')
r=np.array([f['rho'] for f in a['frames']]);R=(r[:,:d*d]+1j*r[:,d*d:]).reshape(-1,d,d)
err=float(np.max(np.abs(R-sol.y.T.reshape(-1,d,d))));mine=float(min(np.linalg.eigvalsh(x).min() for x in R));assert err<2e-6;assert mine>-2e-8
report=json.load(open(p));report['tests'] +=[{'name':'Independent SciPy DOP853 agreement','value':err,'tolerance':2e-6,'passed':True},{'name':'Density-matrix positivity violation','value':max(0,-mine),'tolerance':2e-8,'passed':True}];report['scope']='Default preset, dt=0.5 fs; numerical correctness does not establish material validity.';json.dump(report,open(p,'w'),indent=2);print(json.dumps(report,indent=2))
