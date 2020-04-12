import numpy as np
import numpy.linalg as la
import helm3d_dir as h3
import fmm3dpy as fmm3d

x = np.loadtxt('../geometries/sphere_192_o03.go3')

norder = int(x[0])
npatches = int(x[1])
npols = int((norder+1)*(norder+2)/2)
npts = npatches*npols 

# setup geometry in the correct format
norders = norder*np.ones(npatches)
iptype = np.ones(npatches)
srcvals = x[2::].reshape(12,npts)
ixyzs = np.arange(npatches+1)*npols+1

# convert values to coefs
srccoefs = h3.surf_vals_to_coefs(norders,ixyzs,iptype,srcvals[0:9,:])

wts = h3.get_qwts(norders,ixyzs,iptype,srcvals)

print("error in area of sphere = ",str(sum(wts)-4*np.pi))

# using two sources in exterior
# note: code will currently break if one source is used (fmm interace issue)
#
xyz_out = np.array([[3.17,-0.03,3.15],[6.13,-4.1,2.22]]).transpose()
zk = 1.1 + 1j*0
c = np.array([1 + 1j*0,1+1.1j])
out = fmm3d.h3ddir(zk=zk,sources=xyz_out,targets=srcvals[0:3,:],charges=c,pgt=1)
rhs = out.pottarg

alpha = 3.0
beta = 0
zpars = np.array([zk,alpha,beta],dtype=complex)
eps = 0.51e-6

nifds,nrfds,nzfds = h3.helm_comb_dir_fds_mem(norders,ixyzs,iptype,srccoefs,srcvals,
  eps,zpars)

ifds,rfds,zfds = h3.helm_comb_dir_fds_init(norders,ixyzs,iptype,srccoefs,srcvals,eps,
  zpars,nifds,nrfds,nzfds)

nent = npts**2
col_ptr = np.arange(npts+1)*npts+1
row_ind = np.tile(np.arange(npts)+1,npts)

xmat = h3.helm_comb_dir_fds_matgen(norders,ixyzs,iptype,srccoefs,srcvals,
  eps,zpars,ifds,rfds,zfds,col_ptr,row_ind)

xmat = xmat.reshape(npts,npts).transpose()
sigma = la.solve(xmat,rhs)

xyz_in = np.array([[0.17,-0.03,0.15],[0.13,-0.1,0.22]]).transpose()
ntarg = np.shape(xyz_in)[1]

ipatch_id = -1*np.ones(2)
uvs_targ = np.zeros((2,ntarg)) 

pot_comp = h3.lpcomp_helm_comb_dir(norders,ixyzs,iptype,srccoefs,srcvals,
   xyz_in,ipatch_id,uvs_targ,eps,zpars,sigma)

out = fmm3d.h3ddir(zk=zk,sources=xyz_out,targets=xyz_in,charges=c,pgt=1)
pot_ex = out.pottarg
erra = la.norm(pot_ex-pot_comp)
print("error in solution = ",str(erra))