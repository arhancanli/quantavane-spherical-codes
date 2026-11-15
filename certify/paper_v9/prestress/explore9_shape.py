"""Census version 9, item A2: UNTRUSTED numerical exploration (no decision rests on it). The shape of the extra flexes
of an exact code: the points they move, and whether a flex is a rotation of the moved points.
   python3 -I explore9_shape.py D N DISC
"""
import sys, os
CERT=os.environ.get('QV_CERTIFY', os.path.expanduser('~/quantavane-spherical-codes/certify'))
sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/pylib')); sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/flintlib'))
sys.path.insert(0, CERT+'/census'); sys.path.insert(0, CERT+'/v5')
import numpy as np, scipy.sparse as sps
import core
d, n, disc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
SNAP=os.environ.get('QV_SNAP', os.path.expanduser('~/cohn-table-dspace-153543/files'))
C=core.build(SNAP,d,n,disc); V,E=C['V'],C['E']; nV=len(V); m=len(E)
Y=np.array([[float(x) for x in row] for row in C['PQ']])[:nV]
c=float(np.mean([Y[a]@Y[b] for a,b in E]))
G=Y@Y.T; vals=np.round(G[np.triu_indices(nV,1)],9); print('inner products', sorted(set(vals))[-6:])
rows=[];cols=[];vv=[]
for a in range(nV):
    for k in range(d): rows.append(a); cols.append(a*d+k); vv.append(Y[a,k])
for e,(a,b) in enumerate(E):
    for k in range(d):
        rows+= [nV+e, nV+e]; cols += [a*d+k, b*d+k]; vv += [Y[b,k], Y[a,k]]
M=sps.csr_matrix((vv,(rows,cols)),shape=(nV+m,nV*d))
ev,Vec=np.linalg.eigh((M.T@M).toarray()); null=Vec[:,ev<1e-9]
R=[]
for s_ in range(d):
    for t_ in range(s_+1,d):
        A_=np.zeros((d,d)); A_[s_,t_]=1; A_[t_,s_]=-1; R.append((Y@A_.T).ravel())
Q_,_=np.linalg.qr(np.array(R).T); F=null-Q_@(Q_.T@null); U2,S2,_=np.linalg.svd(F,full_matrices=False); k=int((S2>1e-6).sum()); F=U2[:,:k]
print('extra flexes', k)
for i in range(k):
    f=F[:,i].reshape(nV,d); nr=np.linalg.norm(f,axis=1)
    supp=np.where(nr>1e-8*nr.max())[0]
    print('flex', i, 'support size', len(supp), 'of', nV, 'norm levels', np.unique(np.round(nr[supp]/nr.max(),6))[:8])
# a common support: points moved by some extra flex
nr=np.linalg.norm(F.reshape(nV,d,k),axis=(1,2)); S=np.where(nr>1e-8*nr.max())[0]
print('union support', len(S))
# is the motion of S a rotation fixing the complement? least squares for skew A with A y_a = f_a on S, for each flex
for i in range(k):
    f=F[:,i].reshape(nV,d)
    # unknowns: the g entries of A (upper triangle)
    idx=[(s,t) for s in range(d) for t in range(s+1,d)]
    Mx=np.zeros((len(S)*d, len(idx))); bx=f[S].ravel()
    for j,(s,t) in enumerate(idx):
        A_=np.zeros((d,d)); A_[s,t]=1; A_[t,s]=-1
        Mx[:,j]=(Y[S]@A_.T).ravel()
    sol,res,_,_=np.linalg.lstsq(Mx,bx,rcond=None)
    print('flex', i, 'rotation of the support: residual', np.linalg.norm(Mx@sol-bx)/np.linalg.norm(bx))
# the complement and the contacts between S and the rest
T=np.setdiff1d(np.arange(nV),S); cross=sum(1 for a,b in E if (a in set(S))!=(b in set(S)))
print('fixed points', len(T), 'contacts across', cross, 'contacts inside S', sum(1 for a,b in E if a in set(S) and b in set(S)))
print('rank of the support points', np.linalg.matrix_rank(Y[S],1e-9), 'rank of the fixed points', np.linalg.matrix_rank(Y[T],1e-9) if len(T) else 0)
