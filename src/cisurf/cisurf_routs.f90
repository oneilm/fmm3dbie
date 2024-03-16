!
!  This file contains the routines required to for computing high
!  order meshes from input flat/quadratic triangulations/
!  quadrilaterizations.
!
!  Recall that there are 4 levels of meshes
!  1. User Mesh (U), stores adjacency information, and points
!     defining the input surface. The user mesh is the one that
!     defines \sigma(x) for constructing the level set
!
!     \sigma(x) = \sum_{j=1}^{N_{U}} \sigma_{j} \exp{-|x-c_{j}|^2/2\sigma_{0}^2}/
!                  \sum_{j=1}^{N_{U}} \exp{-|x-c_{j}|^2/2\sigma_{0}^2} \, ,
!     where c_{j} are the centroids of the elements on the user mesh, 
!     N_{U} is the number of elements on the user mesh, 
!     \sigma_{j} = r_{j}/\lambda, 
!     and \sigma_{0} = c_{\sigma_{0}} \max_{j} r_{j}
!     where r_{j} is the bounding radius of the smallest sphere
!     centered at the triangle centroid which completely contains
!     the triangle,
!     and c_{\sigma_{0}}, \lambda are user defined parameters.
!
!     Recommended values for \lambda, and c_{\sigma_{0}} are 5 and 2\sqrt{5}
!     respectively.
!
!     The following information regarding the user mesh is required: 
!      * Points defining the geometry, (nv, verts(3,nv))
!      * Elements defined via their point ids (ne, numesh, umesh_tri_ind(numesh), 
!                                             umesh_tri_se(2,ne)
!      * type of elements in the mesh (ielemtypes)
!      * centroid and radii (umesh_cms(3,ne), umesh_rads(ne) = D_{j}/2)
!      * Normals at the points on the surface (rnormals(3,nv))
!     
!     
!  2. Quadrature Mesh (Q), a refinement of triangles/quads on U
!     such that i
!       \phi(x) = \int_{U} \nabla n(y) \cdot 
!          \nabla_{y} (\erf(|x-y|/sqrt{2} \sigma(x))/4\pi |x-y|) dS 
!       \phi(x) = \int_{U} K(x,y, n(y)) dS
!     is well approximated by the quadrature rule on the mesh Q, i.e. 
!     if y_{j}, n_{j}, and w_{j} are the nodes, normals and weights on
!     the mesh Q, j=1,2, \ldots nquad then
!        |\phi(x) - \sum_{j=1}^{nquad} K(x,y_{j},n_{j}) w_{j}| < \eps_{mach}
!     and a similar result holds for its gradient, for all points
!     x \in R^{3}
!
!     The following information regarding the quadrature mesh is required: 
!      * qpts, qnormals, qwts: discretization nodes, normals, and weights
!                              for computing \phi, and \nabla \phi accurately
!
!  3. Control mesh (C), the mesh which serves as the domain of the charts in
!     the atlas for discretizing \phi = 1/2. Let $\Gamma$ denote 
!     the $\phi = 1/2$ level set.
!     Let $H$ be the pseudonormal vector field constructed on C.
!     Then $\Gamma$ is discretized as $C + hH$ where h is parameterized on 
!     the control mesh and computed on high order discretization rule at 
!     the basis elements.
!
!     The following information regarding the control mesh is required:
!     * Points defining the geometry, (ncv, verts(3,ncv))
!     * Elements defined via their point ids (nce, ncmesh, cmesh_tri_ind(ncmesh), 
!                                             cmesh_tri_se(2,nce)
!     * type of elements in the mesh (ielemtypes(nce))
!     * Normals at the points on the control mesh (rcnormals(3,ncv))
!     * High order discretization nodes (npts, xyzs(3,npts), H(3,npts), h(npts), 
!       hinit(npts)),
!       where, H(3,npts) are the pseudonormals constructed from rcnormals, 
!       and the connectivity information, h is the height function from 
!       the skeleton mesh, that takes us to the \phi = 1/2 surface, and
!       hinit is an initialization for this height function
!
!  4. Discretized surface (S), go3 format of high order phi=1/2 surface
!     with a given U, Q, and S.
!
!  Q is a function of U, and a user defined high order quadrature rule on the
!  elements. There are several options for constructing Q from U, 
!  * Q = U with the high order rule
!  * Q is a uniform refinement of U with the high order rule
!  * Q is determined adaptively starting at U to resolve \phi, and  \nabla \phi
!
!  C in theory can be chosen arbitrarily as long as it is a surface close
!  enough to the \phi=1/2 surface which is watertight. There are
!  several options for choosing C, 
!  * C = U
!  * C is a uniform refinement of U
!  * C is refined adaptively with U as an initial guess such that the certain
!    metrics on Phi=1/2 surface with that given C are well-resolved
! 
!  This file contains the following routines
!    compute_sigma: 
!      Given U, \lambda, c_{\sigma_{0}} compute \sigma
!    compute_qmesh:
!      Determine qmesh given U, \lambda, c_{\sigma_{0}} with 
!      input quadrature order, number of uniform refinements, and 
!      subsequent adaptive refinement if requested. (The adaptive
!      part of this routine might have to be written in C, since
!      the number of elements required are not apriori known)
!    compute_pseudonormals_verts: 
!      Given points, element info, normals, compute psuedonormals
!      at all vertices on the control mesh.
!    compute_pseudonormals_pts:
!      Given pseudonormals at the vertices on the mesh, evaluate
!      the psuedonormals at given interior discretization points
!    compute_phi:
!      Given U, Q, xyz, compute phi, and grad phi 
!    compute_level_set: Given U, Q, C, compute S
!      
!
!  
      subroutine compute_sigma(ncm, cms, rads, c0, rlam, nt, targs, &
         sigma, grad_sigma)
!
!  Let x_{i} = targs_{i}, y_{j} = cms_{j}, and r_{j} = rads_{j}.
!  This routine computes
!  \sigma(x) = \sum_{j=1}^{ncm} \sigma_{j} \exp{-|x-y_{j}|^2/(2\sigma_{0}^2)}/
!               \sum_{j=1}^{ncm} \exp{-|x-y_{j}|^2/(2\sigma_{0}^2)}
!
!  with \sigma_{j} = r_{j}/ rlam, and 
!  \sigma_{0} = c_{0} \max_{j} r_{j}, and its gradient 
!  at x_{i}. This routine uses a primitive version of a gauss transform
!  where the sources are bin sort on boxes of size 6 \sqrt{2} \sigma_{0}, 
!  and then doing a direct calculation with all nearest neighbors
!  
!  Input arguments
!    - ncm: integer
!        number of centroids
!    - cms: real *8 (3,ncm)
!        xyz coordinates of centroids
!    - rads: real *8 (ncm)
!        radii associated with the centroids
!    - c_{0}: real *8
!        scaling parameter for computing \sigma_{0}
!    - rlam: real *8
!        scaling parameter for defining \sigma_{j}
!    - nt: integer
!        number of targets
!    - targs: real *8 (3,nt)
!        xyz coordinates of the targets
!
!  Output arguments:
!    - sigma: real *8 (nt)
!        sigma at the target locations
!    - grad_sigma: real *8(3,nt)
!        gradient of sigma at target locations
!     
      implicit none
      integer, intent(in) :: ncm, nt
      real *8, intent(in) :: c0, rlam
      real *8, intent(in) :: cms(3,ncm), rads(ncm), targs(3,nt)

      real *8, intent(out) :: sigma(nt), grad_sigma(3,nt)

!
!  Temporary variables
!
      real *8, allocatable :: rnum(:), rnum_grad(:,:)
      real *8, allocatable :: rden(:), rden_grad(:,:)
      real *8, allocatable :: charges(:)
      real *8 radmax, sigma0

      integer i,j,k,l


!
!  Compute sigma0, \sigma_{j} which is stored in charges,
!  and initialize the numerator and denominator of 
!  the sigma evaluator
!
      radmax = maxval(rads)
      sigma0 = c0*radmax

      allocate(charges(ncm))
      allocate(rnum(nt), rnum_grad(3,nt), rden(nt), rden_grad(3,nt))
      
!$OMP PARALLEL DO DEAFULT(SHARED) 
      do i=1,ncm
         charges(i) = rads(i)/rlam
      enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nt
        rnum(i) = 0
        rnum_grad(1,i) = 0
        rnum_grad(2,i) = 0
        rnum_grad(3,i) = 0

        rden(i) = 0
        rden_grad(1,i) = 0
        rden_grad(2,i) = 0
        rden_grad(3,i) = 0
      enddo
!$OMP END PARALLEL DO

!
!  Determine bounding box
!
     

!
!
!  Bin sort sources and targets
!
      return
      end

      
