"""Census version 9, item A2: UNTRUSTED numerical exploration (no decision rests on it). The second-order descent test of
Proposition 2.10 for an exact code above the size limit: the extra flexes (numerical kernel of M^T M minus the rotations)
and, for a few flex directions f, the linear program  max eta  s.t.  l_e(w) + Q_e(f) <= -eta  (e in E),  w tangent,
|w| <= 10  (HiGHS interior point). eta > 0 would give a descent direction; by LP duality eta = -max over the nonnegative
stresses of total 1 of H_omega(f)/2 (up to the box).   python3 -I explore9_descent.py D N DISC
"""
import sys, os, time, json
CERT=os.environ.get('QV_CERTIFY', os.path.expanduser('~/quantavane-spherical-codes/certify'))
sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/pylib')); sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/flintlib'))
sys.path.insert(0, CERT+'/census'); sys.path.insert(0, CERT+'/v5')
import numpy as np, scipy.sparse as sps
from scipy.optimize import linprog
import core
d, n, disc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
SNAP=os.environ.get('QV_SNAP', os.path.expanduser('~/cohn-table-dspace-153543/files'))
C=core.build(SNAP,d,n,disc); V,E=C['V'],C['E']; nV=len(V); m=len(E)
PQ=C['PQ']; Y=np.array([[float(x) for x in row] for row in PQ])[:nV]
c=float(np.mean([Y[a]@Y[b] for a,b in E]))
G=Y@Y.T; np.fill_diagonal(G,-2); 
print('c',c,'next below c', np.sort(G[G< c-1e-9])[-3:])
# sparse rigidity matrix
rows=[];cols=[];vals=[]
for a in range(nV):
    for k in range(d): rows.append(a); cols.append(a*d+k); vals.append(Y[a,k])
for e,(a,b) in enumerate(E):
    for k in range(d):
        rows+= [nV+e, nV+e]; cols += [a*d+k, b*d+k]; vals += [Y[b,k], Y[a,k]]
M=sps.csr_matrix((vals,(rows,cols)),shape=(nV+m,nV*d))
t=time.time(); ev,Vec=np.linalg.eigh((M.T@M).toarray()); print('eigh',round(time.time()-t,1))
g=d*(d-1)//2; null=Vec[:,ev<1e-9]; print('null',null.shape[1],'g',g, 'ev gap', ev[null.shape[1]-1], ev[null.shape[1]])
R=[]
for s_ in range(d):
    for t_ in range(s_+1,d):
        A=np.zeros((d,d)); A[s_,t_]=1; A[t_,s_]=-1; R.append((Y@A.T).ravel())
Q_,_=np.linalg.qr(np.array(R).T); F=null-Q_@(Q_.T@null); U,S2,_=np.linalg.svd(F,full_matrices=False); k=int((S2>1e-6).sum()); F=U[:,:k]
print('extra flexes',k)
def Qe(f):
    f=f.reshape(nV,d); return np.array([f[a]@f[b]-c/2*(f[a]@f[a]+f[b]@f[b]) for a,b in E])
# LP for each basis flex and random combos
Ml=M[nV:]   # contact rows: l_e(w)
T=M[:nV]    # tangency rows
best=None
rng=np.random.default_rng(1)
dirs=[F[:,i] for i in range(k)]+[F@rng.standard_normal(k) for _ in range(2)]
for it,f in enumerate(dirs):
    f=f/np.abs(f).max()
    q=Qe(f)
    # variables: w (nV*d), eta ; minimize -eta ; l_e(w) + eta <= -q_e ; T w = 0 ; |w|<=10
    A_ub=sps.hstack([Ml, sps.csr_matrix(np.ones((m,1)))]).tocsr()
    b_ub=-q
    A_eq=sps.hstack([T, sps.csr_matrix((nV,1))]).tocsr(); b_eq=np.zeros(nV)
    cobj=np.zeros(nV*d+1); cobj[-1]=-1
    bounds=[(-10,10)]*(nV*d)+[(None,1.0)]
    t=time.time(); res=linprog(cobj,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=bounds,method='highs-ipm', options={'time_limit': 600.0})
    eta=-res.fun if res.status==0 else None
    print('dir',it,'LP status',res.status,'eta',eta,'Q range',q.min(),q.max(),'time',round(time.time()-t,1))
    if eta is not None and eta>0 and (best is None or eta>best[0]): best=(eta,f,res.x[:-1])
if best:
    eta,f,w=best
    f=f.reshape(nV,d); w=w.reshape(nV,d)
    Gn=G.copy()
    for s in [1e-1,5e-2,2e-2,1e-2,5e-3,2e-3,1e-3]:
        Z=Y+s*f+s*s*w; Z/=np.linalg.norm(Z,axis=1,keepdims=True); H=Z@Z.T; np.fill_diagonal(H,-2)
        print('s',s,'max cos',H.max(),'improvement',c-H.max())
    pass
