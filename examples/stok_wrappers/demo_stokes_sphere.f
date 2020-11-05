      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      
      real *8 ts(2), rres
      real *8, allocatable :: rfacs(:,:), errs(:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),novers(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: strac(:,:),du(:,:)
      real *8, allocatable :: wnear_s(:),wnear_d(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)

      real *8 xyz_out(3),xyz_in(3,10),stracmat(3,3),smat(3,3), dmat(3,3)
      real *8 velgrad(3,3), vel(3), pre, tractemp(3)
      real *8 sigout(3), uin(3), uintest(3,10), dpars(2), st1(3), du1(3)
      real *8 udir(3), uneu(3,10), uavecomp(3), uavetest(3)
      real *8 st2(3), du2(3), uconst(3)
      real *8 v(3), omega(3), r0(3), udiff(3,10), udiff2(3,10)      
      real *8, allocatable :: uval(:,:), tracval(:,:) 
      complex * 16 zpars,ztmp

      call prini(6,13)

      done = 1
      pi = atan(done)*4


c
c  stokes demo on a sphere
c
c  You can change the following three parameters
c
c   ipars(1) decides the refinement on the sphere
c   # of triangles = 12*4**(ipars(1))   
c
c   eps - precision requested for evaluating layer potential
c
c   norder = order of discretization used on each triangle.
c     if norder = p, the (p+1)th order Vioreanu
c     Rokhlin nodes are used in the discretization
c
c

      igeomtype = 1
      ipars(1) = 2
      norder = 3
      eps = 0.51d-5

      npatches = 12*(4**ipars(1))


      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1,1) = 0.17d0
      xyz_in(2,1) = 0.23d0
      xyz_in(3,1) = -0.11d0

      do i = 2,10
         xyz_in(1,i) = 0.2d0*cos(100.0d0*i+5.0d0)
         xyz_in(2,i) = 0.2d0*cos(211.0d0*i+5.0d0)
         xyz_in(3,i) = 0.2d0*cos(357.0d0*i+5.0d0)
      enddo


      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols
      call prinf('n=*',npts,1)
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(targs(3,npts))
      ifplot = 0



      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(ixyzso(npatches+1),novers(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
c
c
c   test green's identity by verfiying that
c    u = S[trac] - D[u] for a known solution
c

      allocate(uval(3,npts),tracval(3,npts))

      sigout(1) = 1.1d0
      sigout(2) = -0.27d0
      sigout(3) = .31d0
      
      do i=1,npts
         call st3d_slp_vec(9,xyz_out,3,srcvals(1,i),0,dpars,0,zpars,0,
     1        ipars,smat)
         call st3d_strac_vec(9,xyz_out,12,srcvals(1,i),0,dpars,0,zpars,
     1        0,ipars,stracmat)
         uval(1,i) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
         uval(2,i) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
         uval(3,i) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)
         tracval(1,i) = stracmat(1,1)*sigout(1)+ stracmat(1,2)*sigout(2)
     1        + stracmat(1,3)*sigout(3)
         tracval(2,i) = stracmat(2,1)*sigout(1)+ stracmat(2,2)*sigout(2)
     1        + stracmat(2,3)*sigout(3)
         tracval(3,i) = stracmat(3,1)*sigout(1)+ stracmat(3,2)*sigout(2)
     1        + stracmat(3,3)*sigout(3)
         
      enddo

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)


      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO      

c
c    find near quadrature correction interactions
c
      call findnearmem(cms,npatches,rad_near,12,srcvals,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,12,srcvals,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = 0 

c
c    estimate oversampling for far-field, and oversample geometry
c

      ztmp = 0

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,12,npts,srcvals,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)

      
c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear_s(6*nquad),wnear_d(6*nquad))

      wnear_s = 0
      wnear_d = 0
      
      iquadtype = 1

c
c
c   compute near quadrature correction for stokeslet
c

      dpars(1) = 1
      dpars(2) = 0

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_stok_comb_vel(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,12,npts,srcvals,
     1     ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     1     iquad,rfac0,nquad,wnear_s)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time for stokeslet=*',t2-t1,1)
      s1 = (npts+0.0d0)/(t2-t1)
      call prin2('quadrutre generation speed in pps=*',s1,1)
c
c   Apply layer potential to evaluate vecloity
c
      allocate(strac(3,npts))
      strac = 0
      call lpcomp_stok_comb_vel_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,12,npts,srcvals,eps,
     2   dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear_s,
     3   tracval,novers,npts_over,ixyzso,srcover,wover,strac)
     
      dpars(1) = 0
      dpars(2) = 1
      

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_stok_comb_vel(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,12,npts,srcvals,
     1     ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     1     iquad,rfac0,nquad,wnear_d)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time for stresslet=*',t2-t1,1)
      s1 = (npts+0.0d0)/(t2-t1)
      call prin2('quadrutre generation speed in pps=*',s1,1)

      allocate(du(3,npts))
      du = 0
      call lpcomp_stok_comb_vel_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,12,npts,srcvals,eps,
     2   dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear_d,
     3   uval,novers,npts_over,ixyzso,srcover,wover,du)


      erra = 0
      ra = 0
      do i=1,npts
        do j=1,3
          uin(j) = (strac(j,i) - du(j,i))/2/pi
          erra = erra + abs(uin(j)-uval(j,i))**2*wts(i)
          ra = ra + abs(uval(j,i))**2*wts(i)
        enddo
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in green identity=*',erra,1)
      
      stop
      end


      


      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

