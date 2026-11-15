"""Census version 9, item A2: UNTRUSTED numerical exploration (no decision rests on it). Does an exact code whose extra
flexes carry H_omega = 0 lie in a continuous family of configurations with the same contacts? From the exact core, the
extra flexes are computed numerically; the second-order terms Q_e(f) = <f_a,f_b> - c(|f_a|^2+|f_b|^2)/2 are evaluated on
the contacts; then, for steps h along a flex, Gauss-Newton (sparse normal equations) projects x + h f back onto
{|x_a| = 1, <x_a,x_b> = c on the contacts}, and the residual and the largest non-contact inner product are reported.
   python3 -I explore9_family.py D N DISC
"""
import sys, os, time
CERT=os.environ.get('QV_CERTIFY', os.path.expanduser('~/quantavane-spherical-codes/certify'))
sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/pylib')); sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/flintlib'))
sys.path.insert(0, CERT+'/census'); sys.path.insert(0, CERT+'/v5')
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spla
import core
d, n, disc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
SNAP=os.environ.get('QV_SNAP', os.path.expanduser('~/cohn-table-dspace-153543/files'))
C=core.build(SNAP,d,n,disc); V,E=C['V'],C['E']; nV=len(V); m=len(E)
Y=np.array([[float(x) for x in row] for row in C['PQ']])[:nV]
c=float(np.mean([Y[a]@Y[b] for a,b in E]))
Ea=np.array([a for a,b in E]); Eb=np.array([b for a,b in E])
def jac(X):
    rows=[];cols=[];vals=[]
    for a in range(nV):
        for k in range(d): rows.append(a); cols.append(a*d+k); vals.append(2*X[a,k])
    for e,(a,b) in enumerate(E):
        for k in range(d):
            rows+= [nV+e, nV+e]; cols += [a*d+k, b*d+k]; vals += [X[b,k], X[a,k]]
    return sps.csr_matrix((vals,(rows,cols)),shape=(nV+m,nV*d))
def res(X):
    return np.concatenate([(X*X).sum(1)-1, (X[Ea]*X[Eb]).sum(1)-c])
J=jac(Y)
t0=time.time(); ev,Vec=np.linalg.eigh((J.T@J).toarray()); print('eigh', round(time.time()-t0,1), flush=True)
null=Vec[:,ev<1e-9]
R=[]
for s_ in range(d):
    for t_ in range(s_+1,d):
        A_=np.zeros((d,d)); A_[s_,t_]=1; A_[t_,s_]=-1; R.append((Y@A_.T).ravel())
Q_,_=np.linalg.qr(np.array(R).T); F=null-Q_@(Q_.T@null); U2,S2,_=np.linalg.svd(F,full_matrices=False); k=int((S2>1e-6).sum()); F=U2[:,:k]
print('extra flexes', k, flush=True)
G0=Y@Y.T; np.fill_diagonal(G0,-2); nc=np.ones_like(G0,bool); nc[Ea,Eb]=False; nc[Eb,Ea]=False; np.fill_diagonal(nc,False)
print('largest non-contact inner product at t=0', G0[nc].max(), flush=True)
for i in range(k):
    f=F[:,i].reshape(nV,d)
    Qe=(f[Ea]*f[Eb]).sum(1)-c/2*((f[Ea]**2).sum(1)+(f[Eb]**2).sum(1))
    print('flex', i, 'max |Q_e(f)| over the contacts', np.abs(Qe).max(), 'max |f|', np.abs(f).max(), flush=True)
f=F[:,0].reshape(nV,d)/np.abs(F[:,0]).max()
for h in (1e-3, 1e-2, 3e-2, 1e-1):
    X=Y+h*f
    X/=np.linalg.norm(X,axis=1,keepdims=True)
    for it in range(12):
        r=res(X); nr=np.abs(r).max()
        if nr<1e-14: break
        Jx=jac(X)
        dx=spla.lsqr(Jx, -r, atol=1e-16, btol=1e-16, iter_lim=20000)[0]
        X=X+dx.reshape(nV,d)
    r=res(X); G=X@X.T; np.fill_diagonal(G,-2)
    # distance to the orbit of Y (Procrustes)
    U_,_,Vt_=np.linalg.svd(X.T@Y); Xr=X@(U_@Vt_)
    print('h', h, 'residual', np.abs(r).max(), 'iterations', it, 'max non-contact', G[nc].max(), 'max cos', G.max(),
          'distance to the orbit', np.linalg.norm(Xr-Y), flush=True)
