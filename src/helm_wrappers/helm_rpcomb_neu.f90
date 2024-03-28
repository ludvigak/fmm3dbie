

      subroutine getnearquad_helm_rpcomb_neu(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = S_{k}[\rho] + i*alpha*D_{k}[S_{i|k|}[\rho]]    -   (1)
!
!  and returns quantities related to evaluating du/dn on surface
!  at the surface discretization nodes
!
!  If values at other nodes is desired then the solution
!  should be reinterpolated to those nodes
!
!  On imposing the boundary condition, we get the following operator
!
!  du/dn = -I/2 + S_{k}' + i \alpha (D_{k}'-D_{i|k|}') S_{i|k|}
!    - i \alpha I/4 + i \alpha (S_{i|k|}')^2
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!  
!  NOTES:
!    - wnear must be of size (4,nquad) as 4 different layer
!      potentials are returned
!      * the first kernel is S_{k}'
!      * the second kernel is S_{i|k|}
!      * the third kernel is S_{i|k|}'
!      * the fourth kernel is D_{k}'-D_{i|k|}'
! 
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev 
!                     nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!        zpars(1) = k 
!        zpars(2) = alpha
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: complex *16(4,nquad)
!        The desired near field quadrature
!        wnear(1,:) - stores the quadrature corrections for S_{k}'
!        wnear(2,:) - stores the quadrature correction for S_{i|k|}        
!        wnear(3,:) - stores the quadrature correction for S_{i|k|}'
!        wnear(4,:) - stores the quadrature correction for 
!                     D_{k}'-D_{i|k|}'
!               
!
  
      implicit none 
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(2)
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer, intent(in) :: nquad
      complex *16, intent(out) :: wnear(4,nquad)
      
      complex *16 zpars_tmp(3)
      integer ipars(2)
      real *8 dpars(1)
      

      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)


      complex *16 alpha, beta, ima, zk
      integer i, j, ndi, ndd, ndz, idim

      integer ipv
      complex *16, allocatable :: wneartmp(:)
      
      integer ndtarg, ntarg

      procedure (), pointer :: fker
      external h3d_sprime, h3d_slp, h3d_dprime_diff

      data ima/(0.0d0,1.0d0)/

      ndz=1
      ndd=1
      ndi=2
      ndtarg = 12
      ntarg = npts
      zk = zpars(1)

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
      do i=1,nquad
        do idim = 1,4
          wnear(idim,i) = 0
        enddo
      enddo
!$OMP END PARALLEL DO

      allocate(wneartmp(nquad))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      if (iquadtype.eq.1) then
        ipv=1

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => h3d_sprime 
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(1,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO

        zpars_tmp(1) = ima*abs(zk)
        fker => h3d_slp
        ipv = 1
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(2,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO

        fker => h3d_sprime
        ipv = 1
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(3,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO

        zpars_tmp(1) = zk
        zpars_tmp(2) = ima*abs(zk)
        ndz = 2
        fker => h3d_dprime_diff
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(4,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO
      
      endif

      return
      end subroutine getnearquad_helm_rpcomb_neu
!
!
!
!
!
      subroutine lpcomp_helm_rpcomb_neu_addsub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, &
        ixyzso, srcover, whtsover, lwork, work, ndim, sigma, pot)

!
!  This subroutine evaluates the neumann data corresponding to
!  the following integral representation:
!  
!  u = S_{k}[\sigma]+i*alpha*D_{k}[S_{i|k|}[\sigma]]
!
!  du/dn = S_{k}'[\sigma] + i\alpha S_{i|k|}'^2 [\sigma] + 
!    i \alpha (D_{k}' - D_{i|k|}') S_{i|k|}
!
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!        
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested
!    - ndd: integer
!        number of real parameters defining the kernel/
!        integral representation (unused in this routine)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation (must be two)
!    - zpars: real *8(ndz)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
!        * zpars(2) = alpha
!    - ndi: integer
!        number of integer parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ipars: real *8(ndi)
!        integer parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction, must be 4
!    - wnear: complex *16(nker, nquad)
!        precomputed quadrature corrections            
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - lwork: integer
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - ndim: integer
!        number of densities per point on the surface,
!        must be 1 for this routine
!    - sigma: complex *16(npts)
!        The density sigma above                                        
!
!  Output arguments:
!    - pot: complex *16 (npts)
!        du/dn corresponding to representation
!
  
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, intent(in) :: eps
      
      integer, intent(in) :: ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)

      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      integer, intent(in) :: nker
      complex *16, intent(in) :: wnear(nker,nquad)

      integer, intent(in) :: nptso
      integer, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      integer, intent(in) :: ndim
    
      complex *16, intent(in) :: sigma(npts)
    
      complex *16, intent(out) :: pot(npts)

      real *8, allocatable :: sources(:,:), srctmp(:,:)
      complex *16, allocatable :: charges(:), dipvec(:,:), sigmaover(:)
      complex *16, allocatable :: sigma_aux(:)

      integer npols
      integer ns, ntarg
      integer ifcharge, ifdipole
      integer ifpgh, ifpghtarg
    
 

      integer i, j, jpatch, jquadstart, jstart
      complex *16 pottmp, gradtmp(3)
      complex *16, allocatable :: pot_aux(:), grad_aux(:,:)
      complex *16, allocatable :: phi1(:), phi2(:)

!$    real *8 omp_get_wtime      

      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:), dtmp2(:,:)
      real *8 thresh
  
      real *8 over4pi
      integer nss, l, ier
      complex *16 ima, ztmp

      integer nd, ntarg0, nmax

      real *8 done, pi, rr
      integer nddtmp, nditmp, ndztmp
      complex *16 zpars2(2)
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (nd=1,ntarg0=1)

      ns = nptso
      ntarg = npts
      done = 1
      pi = atan(done)*4

      ifpgh = 0
      ifpghtarg = 2
      
      allocate(sources(3,ns), srctmp(3,npts))
      allocate(charges(ns), dipvec(3,ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(npts), grad_aux(3,npts))
      allocate(phi1(npts), phi2(npts))
      allocate(sigma_aux(npts))
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), ctmp2(nmax), dtmp2(3,nmax))
! 
!       oversample density
!
      call oversample_fun_surf(2, npatches, norders, ixyzs, iptype, & 
        npts, sigma, novers, ixyzso, ns, sigmaover)



!
!  extract source and target info
!
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i) = sigmaover(i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        srctmp(1,i) = srcvals(1,i)
        srctmp(2,i) = srcvals(2,i)
        srctmp(3,i) = srcvals(3,i)
        pot_aux(i) = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
      enddo

!$OMP END PARALLEL DO

! 
!  Compute S_{k}'
      ier = 0
      call hfmm3d_t_c_g(eps, zpars(1), ns, sources, charges, npts, &
        srctmp, pot_aux, grad_aux, ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot(i) = grad_aux(1,i)*srcvals(10,i) + &
          grad_aux(2,i)*srcvals(11,i) + &
          grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    
      
!        compute threshold for ignoring local computation
      call get_fmm_thresh(12, ns, srcover, 12, npts, srcvals, thresh)
!
!       Add near field precomputed contribution
!
    

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(1,jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!
      ifcharge = 1
      ifdipole = 0 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(ctmp2, l, jstart, nss, pottmp, gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss) = charges(l)
          enddo
        enddo

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call h3ddirectcg(nd, zpars(1), srctmp2, ctmp2, nss, &
          srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        pot(i) = pot(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i) - gradtmp(3)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      

!
!
!  Now handle the computation of S_{i|k|}[\sigma] and 
!  S_{i|k|}'[\sigma] and hold them in separate arrays phi1 and phi2
!
!
      ztmp = ima*abs(zpars(1))
      ier = 0
      call hfmm3d_t_c_g(eps, ztmp, ns, sources, charges, npts, &
        srctmp, pot_aux, grad_aux, ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        phi1(i) = pot_aux(i)
        phi2(i) = grad_aux(1,i)*srcvals(10,i) + &
          grad_aux(2,i)*srcvals(11,i) + &
          grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    

!
!  Add near field precomputed contribution
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, pottmp, npols, l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            phi1(i) = phi1(i) + &
              wnear(2,jquadstart+l-1)*sigma(jstart+l-1)
            phi2(i) = phi2(i) + &
              wnear(3,jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(ctmp2, nss, l, pottmp, gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss)= charges(l) 
          enddo
        enddo
        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0
        call h3ddirectcg(nd, ztmp, srctmp2, ctmp2, nss, srctmp(1,i), &
          ntarg0, pottmp, gradtmp, thresh)
        phi1(i) = phi1(i) - pottmp  
        phi2(i) = phi2(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i) - gradtmp(3)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
!
!  End of computing phi1 = S_{i|k|}[\sigma]
!  and phi2 = S_{i|k|}'[\sigma]
!

!
!  Now compoute S_{i|k|}'[\phi2] and add i\alpha S_{i|k|}'[phi2]
!  to pot
!
      call oversample_fun_surf(2, npatches, norders, ixyzs, iptype, & 
        npts, phi2, novers, ixyzso, ns, sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        charges(i) = sigmaover(i)*whtsover(i)*over4pi 
      enddo
!$OMP END PARALLEL DO

      call hfmm3d_t_c_g(eps, ztmp, ns, sources, charges, npts, &
        srctmp, pot_aux, grad_aux, ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot_aux(i) = grad_aux(1,i)*srcvals(10,i) + &
          grad_aux(2,i)*srcvals(11,i) + &
          grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    

!
!       Add near field precomputed contribution
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot_aux(i) = pot_aux(i) + &
              wnear(3,+jquadstart+l-1)*phi2(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(ctmp2, nss, l, jstart, pottmp, gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss) = charges(l) 
          enddo
        enddo
        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call h3ddirectcg(nd, ztmp, srctmp2, ctmp2, nss, srctmp(1,i), &
          ntarg0, pottmp, gradtmp, thresh)
        pot_aux(i) = pot_aux(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i) - gradtmp(3)*srcvals(12,i)
        pot(i) = pot(i) + ima*zpars(2)*pot_aux(i)
      enddo
!$OMP END PARALLEL DO      

!
! End of adding i \alpha S_{ik}'^2[\sigma] to pot
!
!

!
!  Begin computation of D_{k}'[\phi1] 
!  we will not handle to subtraction of the near correction
!  until compouting D_{i|k|}'[\phi1] and take care of
!  the total subtraction together
!
!  the array phi2 is no longer neeeded, so we will reuse it
!  for temporary storage
!

      call oversample_fun_surf(2, npatches, norders, ixyzs, iptype, & 
        npts, phi1, novers, ixyzso, ns, sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*over4pi 
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*over4pi 
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*over4pi 
      enddo
!$OMP END PARALLEL DO

      call hfmm3d_t_d_g(eps, zpars(1), ns, sources, dipvec, npts, &
        srctmp, pot_aux, grad_aux, ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        phi2(i) = grad_aux(1,i)*srcvals(10,i) + &
          grad_aux(2,i)*srcvals(11,i) + grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO
      

      call hfmm3d_t_d_g(eps, ztmp, ns, sources, dipvec, npts, &
        srctmp, pot_aux, grad_aux, ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot_aux(i) = phi2(i) - grad_aux(1,i)*srcvals(10,i) - &
          grad_aux(2,i)*srcvals(11,i) - grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO
      

!
!       Add near field precomputed contribution
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, pottmp, npols, l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot_aux(i) = pot_aux(i) + &
              wnear(4,jquadstart+l-1)*phi1(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

      zpars2(1) = zpars(1)
      zpars2(2) = ztmp
      nddtmp = 0
      nditmp = 0
      ndztmp = 2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, l, jpatch, pottmp, rr)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            pottmp = 0
            rr = sqrt((srcover(1,l) - srcvals(1,i))**2 + &
              (srcover(2,l) - srcvals(2,i))**2 + &
              (srcover(3,l) - srcvals(3,i))**2)
            if(rr.lt.thresh) goto 1311
            call h3d_dprime_diff(srcover(1,l), 12, srcvals(1,i), &
              nddtmp, dpars, ndztmp, zpars2, nditmp, ipars, pottmp)
            pot_aux(i) = pot_aux(i) - pottmp*sigmaover(l)*whtsover(l)
 1311       continue            
          enddo
        enddo
        pot(i) = pot(i)+ ima*zpars(2)*pot_aux(i)
      enddo
!$OMP END PARALLEL DO      

!
! End of adding i \alpha (D_{k}' - D_{ik}')S_{ik}'[\sigma] to pot
!


      return
      end subroutine lpcomp_helm_rpcomb_neu_addsub
!
!
!
!
!

      subroutine helm_rpcomb_neu_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, eps_gmres, niter, errs, rres, soln, siksoln)
!
!
!  This subroutine solves the Helmholtz Neumann problem
!  on the exterior of an object where the potential
!  is represented as a right preconditioned 
!  combined field integral representation.
!
!
!  Representation:
!    u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  Boundary condition:
!    u'=f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: complex *16(npts)
!        Neumann data
!    - eps_gmres: real *8
!        gmres tolerance requested
!      
!
!  output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: complex *16(npts)
!        density which solves the neumann problem \rho
!    - siksoln: complx *16(npts)
!        sik[\rho] which can be used for far field
!        computations later
!				 


      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      integer, intent(in) :: ifinout
      complex *16, intent(in) :: zpars(2)
      complex *16, intent(in) :: rhs(npts)
      integer, intent(in) :: numit


      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer, intent(out) :: niter

      complex *16, intent(out) :: soln(npts), siksoln(npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg


      real *8 eps2
      integer nover, npolso, nptso
      integer nnz, nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)

      complex *16, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer i, j, jstart

      integer ipars(2)
      real *8 t1, t2
!$    real *8 omp_get_wtime

      integer nker
      real *8 rfac, rfac0
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over
      complex *16 ima
  

      ima=(0.0d0,1.0d0)

!
!        setup targets as on surface discretization points
! 
      ndtarg = 12  
      ntarg = npts
      allocate(targs(ndtarg,npts), uvs_targ(2,ntarg), ipatch_id(ntarg))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(:,i) = srcvals(:,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)
    

      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, & 
        srccoefs, cms, rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      call findnearmem(cms, npatches, rad_near, ndtarg, targs, npts, &
         nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, targs, npts, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,&
         iquad)

      ikerorder = -1
!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
       rads, npts, srccoefs, ndtarg, npts, targs, ikerorder, zpars(1), &
       nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)
!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      nker = 4
      allocate(wnear(4,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        do j=1,4
          wnear(j,i)=0
        enddo
      enddo
!$OMP END PARALLEL DO    

      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      

      call getnearquad_helm_rpcomb_neu(npatches, norders, &
       ixyzs, iptype, npts, srccoefs, srcvals, &
       eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
       iquad, rfac0, nquad, wnear)

      call cpu_time(t2)
!$      t2 = omp_get_wtime()     
      
      print *, "done generating near quadrature, now starting gmres"

      call helm_rpcomb_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
        errs, rres, soln, siksoln)
 

      return
      end
!
!
!
!
!
      subroutine helm_rpcomb_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
        errs, rres, soln, siksoln)
      !
!
!  This subroutine solves the Helmholtz Neumann problem
!  on the exterior of an object where the potential
!  is represented as a right preconditioned 
!  combined field integral representation.
!
!
!  Representation:
!    u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  Boundary condition:
!    u'=f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  Input arguments:
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: complex *16(npts)
!        Neumann data
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction
!    - wnear: real *8(nker, nquad)
!        precomputed quadrature corrections 
!        wnear(1,:) - stores the quadrature corrections for S_{k}'
!        wnear(2,:) - stores the quadrature correction for S_{i|k|}        
!        wnear(3,:) - stores the quadrature correction for S_{i|k|}'
!        wnear(4,:) - stores the quadrature correction for 
!                     D_{k}'-D_{i|k|}'
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - eps_gmres: real *8
!        gmres tolerance requested
!      
!
!  output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: complex *16(npts)
!        density which solves the neumann problem \rho
!    - siksoln: complx *16(npts)
!        sik[\rho] which can be used for far field
!        computations later
!				 


      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      integer, intent(in) :: ifinout
      complex *16, intent(in) :: zpars(2)
      complex *16, intent(in) :: rhs(npts)
      integer, intent(in) :: numit

      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer, intent(out) :: niter
      
      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)

      integer, intent(in) :: nker
      complex *16, intent(in) :: wnear(nker,nquad)

      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches), ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
  
      complex *16, intent(out) :: soln(npts), siksoln(npts)

      complex *16 zid

      procedure (), pointer :: fker
      external lpcomp_helm_rpcomb_neu_addsub

      integer ndd, ndi, ndz, lwork, ndim
      real *8 dpars, work
      integer ipars

      integer ndtarg
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      real *8, allocatable :: wts(:)

      complex *16 zpars_use(3)
      !

      zid=((-1)**(ifinout)/2.0d0-ima*zpars(2)/4)
      fker => lpcomp_helm_rpcomb_neu_addsub

      ndd = 0
      ndi = 0
      ndz = 2

      lwork = 0
      ndim = 1
      allocate(wts(npts))

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
      srcvals, wts)
!

      call zgmres_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, wts, &
            eps, ndd, dpars, ndz, zpars, ndi, ipars, &
            nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
            nptso, ixyzso, srcover, whtsover, lwork, work, &
            ndim, fker, zid, rhs, numit, eps_gmres, niter, errs, &
            rres, soln)

!   Now compute sik(soln)            

      zpars_use(1) = ima*zpars(1)
      zpars_use(2) = 1
      zpars_use(3) = 0
      
      ndtarg = 12
      call lpcomp_helm_comb_dir_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
        eps, zpars_use, nnz, row_ptr, col_ind, iquad, nquad, &
        wnear(2,1:nquad), soln, novers, nptso, ixyzso, srcover, &
        whtsover, siksoln)


      return
      end
!
!
!      
!
      subroutine helm_rpcomb_neu_solver_memest(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,zpars,numit,rmem)
!
!
!  This subroutine is the memory estimation routine for
!  the Helmholtz Neumann problem
!  solver where the potential
!  is represented as a right preconditioned 
!  combined field integral representation.
!
!
!  Representation:
!    u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  Boundary condition:
!    u'=f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  Input:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - numit: integer
!        max number of gmres iterations
!      
!
!  output
!    - rmem: double precision
!        estimated memory required by code in GB. Note that
!        this is meant to serve as an estimate only. 
!        The exact memory usage might be between (0.75,1.25)*rmem
!        The memory estimate may not be reliable for a
!        very small number of points
! 
!
!
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      complex *16 zpars(3)
      real *8 rmem

      integer *8 lmem8,bigint
      real *8 rmemfmm
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)

      real *8, allocatable :: srcover(:,:),wover(:),sources(:,:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars,ifcharge,ifdipole,ifpgh,ifpghtarg
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over,iper

!
!
!       gmres variables
!
      integer numit,k,l
      complex *16 temp
      
      bigint = numit+1
      bigint = bigint*npts*2
      lmem8 = lmem8 + bigint



      lmem8 = lmem8 + numit*(numit+5)*2 + npts*2


      done = 1
      pi = atan(done)*4


!
!
!        setup targets as on surface discretization points
! 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))
      lmem8 = lmem8 + ndtarg*npts + 3*ntarg 

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
       ipatch_id,uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      lmem8 = lmem8 + 5*npatches

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, &
             col_ind)

      allocate(iquad(nnz+1)) 
      lmem8 = lmem8 + npts+1+ 2*nnz
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind, &
              iquad)

      ikerorder = 0


!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))
      lmem8 = lmem8 + 2*npatches + 1

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
         rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1), &
         nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))
      lmem8 = lmem8 + 15*npts_over

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over, &
             srcover,wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      lmem8 = lmem8 + nquad*2*4
      rmem = lmem8*8/1024/1024/1024

      ifcharge = 1
      ifdipole = 0
      if(abs(zpars(2)).gt.1.0d-16) ifdipole = 1

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,npts_over))
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts_over
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      

      iper = 0
      rmemfmm = 0
      call hfmm3d_memest(1,eps,zpars(1),npts_over,sources,ifcharge, &
        ifdipole,iper,ifpgh,npts,targs,ifpghtarg,rmemfmm)
      rmem = rmem + rmemfmm
      
!
      return
      end
!
!
!
!
!
!
!
      subroutine lpcomp_helm_rpcomb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars,sigma,sigma1,pot)
!
!
!  This subroutine evaluates the dirichlet data for
!  the right preconditioned
!  combined field representation.
!
!  Representation:
!    u = S_{k}[\sigma]+i*alpha*D_{k}[S_{i|k|}[\sigma]]
!
!
!  Input:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!        target info, the first three coordinates must be
!        the xyz components of the target
!    - ipatch_id: integer(ntarg)
!        patch on which target is on if on-surface, =-1/0 if 
!        target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates on patch, if target is on surface
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - sigma: complex *16(npts)
!        density sigma above
!    - sigma1: complex *16 (npts)
!        sik(sigma)
!
!  Output arguments:
!    - pot: complex *16(ntarg)
!        potential at the target locations
!				 
!
      implicit none
      integer, intent(in) :: npatches
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg),eps
      complex *16, intent(in) :: zpars(2),sigma(npts),sigma1(npts)
      complex *16, intent(out) :: pot(ntarg)

      complex *16, allocatable :: pottmp2(:)
      complex *16 zpars_tmp(3),ima
      integer, allocatable :: ipatch_id_src(:)
      real *8, allocatable :: uvs_src(:,:)
      integer i,ndtarg0

      data ima/(0.0d0,1.0d0)/


      allocate(pottmp2(ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pottmp2(i) = 0
        pot(i) = 0
      enddo
!$OMP END PARALLEL DO
      
      
!
!  compute S_{k} [sigma]
!
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = 1.0d0 
      zpars_tmp(3) = 0

      
      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars_tmp,sigma,pot)
      
      
!
!  compute D_{k} [pottmp]
!
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = 0
      zpars_tmp(3) = 1.0d0

      
      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars_tmp,sigma1,pottmp2)


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pot(i) = pot(i) + pottmp2(i)*zpars(2)*ima
      enddo
!$OMP END PARALLEL DO

     
      return
      end subroutine lpcomp_helm_rpcomb_dir

