"""Census version 9, item A2: UNTRUSTED numerical exploration (no decision rests on it). The second-order form H_omega on
the extra flexes of an exact code above the size limit, for a basis of the equilibrium stresses constant on the classes of
contacts (by numbers of common neighbours, as run_above.py).   python3 -I explore9_classes.py D N DISC
"""
import sys, os, time
CERT=os.environ.get('QV_CERTIFY', os.path.expanduser('~/quantavane-spherical-codes/certify'))
sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/pylib')); sys.path.insert(0, os.path.expanduser('~/Library/Caches/quantavane-census/flintlib'))
sys.path.insert(0, CERT+'/census'); sys.path.insert(0, CERT+'/v5'); sys.path.insert(0, CERT+'/paper_v8/above')
import numpy as np, scipy.sparse as sps, scipy.linalg as sl
import core, run_above as RA
d, n, disc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
SNAP=os.environ.get('QV_SNAP', os.path.expanduser('~/cohn-table-dspace-153543/files'))
C=core.build(SNAP,d,n,disc); V,E=C['V'],C['E']; nV=len(V); m=len(E)
Y=np.array([[float(x) for x in row] for row in C['PQ']])[:nV]
c=float(np.mean([Y[a]@Y[b] for a,b in E]))
key=RA.neighbour_classes(E,nV); classes=sorted(set(key)); ci={k:i for i,k in enumerate(classes)}
print('classes', len(classes), [sum(1 for x in key if x==k) for k in classes])
# equilibrium: for each point a, sum_b w_ab y_b - lam_a y_a = 0  -> unknowns: class weights (nc) and lam (nV)
nc=len(classes)
A=np.zeros((nV*d, nc+nV))
for e,(a,b) in enumerate(E):
    A[a*d:(a+1)*d, ci[key[e]]]+=Y[b]; A[b*d:(b+1)*d, ci[key[e]]]+=Y[a]
for a in range(nV): A[a*d:(a+1)*d, nc+a]=-Y[a]
U,S,Vt=np.linalg.svd(A, full_matrices=True); ns=Vt[(S>1e-9).sum():]
print('class stress space dim', ns.shape[0], 'singular values tail', S[-4:])
basis=[v[:nc] for v in ns]
# flexes
rows=[];cols=[];vals=[]
for a in range(nV):
    for k in range(d): rows.append(a); cols.append(a*d+k); vals.append(Y[a,k])
for e,(a,b) in enumerate(E):
    for k in range(d):
        rows+= [nV+e, nV+e]; cols += [a*d+k, b*d+k]; vals += [Y[b,k], Y[a,k]]
M=sps.csr_matrix((vals,(rows,cols)),shape=(nV+m,nV*d))
ev,Vec=np.linalg.eigh((M.T@M).toarray()); null=Vec[:,ev<1e-9]
R=[]
for s_ in range(d):
    for t_ in range(s_+1,d):
        A_=np.zeros((d,d)); A_[s_,t_]=1; A_[t_,s_]=-1; R.append((Y@A_.T).ravel())
Q_,_=np.linalg.qr(np.array(R).T); F=null-Q_@(Q_.T@null); U2,S2,_=np.linalg.svd(F,full_matrices=False); k=int((S2>1e-6).sum()); F=U2[:,:k]
print('extra flexes', k)
def Hmat(w):
    H=np.zeros((k,k))
    for i in range(k):
        fi=F[:,i].reshape(nV,d)
        for j in range(k):
            fj=F[:,j].reshape(nV,d)
            H[i,j]=sum(w[ci[key[e]]]*(fi[a]@fj[b]+fi[b]@fj[a]-c*fi[a]@fj[a]-c*fi[b]@fj[b]) for e,(a,b) in enumerate(E))
    return H
for t,bv in enumerate(basis):
    print('basis', t, 'class weights', np.round(bv/np.abs(bv).max(),4), 'H eig', np.linalg.eigvalsh(Hmat(bv)))
