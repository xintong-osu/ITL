c-----------------------------------------------------------------------
      subroutine exact(uu,vv,xx,yy,n,time,visc,u0,v0)
c
c     This routine creates initial conditions for an exact solution
c     to the Navier-Stokes equations based on the paper of Walsh [1],
c     with an additional translational velocity (u0,v0).
c     
c     The computational domain is [0,2pi]^2 with doubly-periodic 
c     boundary conditions.
c     
c     Walsh's solution consists of an array of vortices determined 
c     as a linear combinations of eigenfunctions of having form:
c     
c         cos(pi m x)cos(pi n y), cos(pi m x)sin(pi n y)
c         sin(pi m x)cos(pi n y), sin(pi m x)sin(pi n y)
c     
c     and
c
c         cos(pi k x)cos(pi l y), cos(pi k x)sin(pi l y)
c         sin(pi k x)cos(pi l y), sin(pi k x)sin(pi l y)
c     
c     While there are constraints on admissible (m,n),(k,l)
c     pairings, Walsh shows that there is a large class of
c     possible pairings that give rise to very complex vortex
c     patterns.
c     
c     Walsh's solution applies either to unsteady Stokes or 
c     unsteady Navier-Stokes.  The solution is a non-translating
c     decaying array of vortices that decays at the rate 
c
c          exp ( -4 pi^2 (m^2+n^2) visc time ),
c
c     with (m^2+n^2) = (k^2+l^2). A nearly stationary state may
c     be obtained by taking the viscosity to be extremely small,
c     so the effective decay is negligible.   This limit, however,
c     leads to an unstable state, thus diminsishing the value of 
c     Walsh's solution as a high-Reynolds number test case.
c
c     It is possible to extend Walsh's solution to a stable convectively-
c     dominated case by simulating an array of vortices that translate
c     at arbitrary speed by adding a constant to the initial velocity field.  
c     This approach provides a good test for convection-diffusion dynamics.
c     
c     The approach can also be extended to incompressible MHD with unit
c     magnetic Prandtl number Pm.
c     
c [1] Owen Walsh, "Eddy Solutions of the Navier-Stokes Equations,"
c     in The Navier-Stokes Equations II - Theory and Numerical Methods,
c     Proceedings, Oberwolfach 1991, J.G. Heywood, K. Masuda,
c     R. Rautmann,  S.A. Solonnikov, Eds., Springer-Verlag, pp. 306--309
c     (1992).
c
c     2/23/02; 6/2/09;  pff
c
c
      include 'SIZE'
      include 'INPUT'
c
      real uu(n),vv(n),xx(n),yy(n)
c
      real cpsi(2,5), a(2,5)
      save cpsi     , a

c     data a / .4,.45 , .4,.2 , -.2,-.1 , .2,.05, -.09,-.1 / ! See eddy.m
c     data cpsi / 0,65 , 16,63 , 25,60 , 33,56 , 39,52 /     ! See squares.f
c     data cpsi / 0,85 , 13,84 , 36,77 , 40,75 , 51,68 /


c     This data from Walsh's Figure 1 [1]:

      data a / -.2,-.2, .25,0.,   0,0  ,  0,0  ,  0,0  /
      data cpsi / 0, 5 ,  3, 4 ,  0,0  ,  0,0  ,  0,0  /

      one   = 1.
      pi    = 4.*atan(one)

      aa    = cpsi(2,1)**2
      arg   = -visc*time*aa  ! domain is [0:2pi]
      e     = exp(arg)
c
c     ux = psi_y,  uy = -psi_x
c
      do i=1,n
         x = xx(i) - u0*time
         y = yy(i) - v0*time

         sx = sin(cpsi(2,1)*x)
         cx = cos(cpsi(2,1)*x)
         sy = sin(cpsi(2,1)*y)
         cy = cos(cpsi(2,1)*y)
         u  =  a(1,1)*cpsi(2,1)*cy 
         v  =  a(2,1)*cpsi(2,1)*sx

         do k=2,5
            s1x = sin(cpsi(1,k)*x)
            c1x = cos(cpsi(1,k)*x)
            s2x = sin(cpsi(2,k)*x)
            c2x = cos(cpsi(2,k)*x)

            s1y = sin(cpsi(1,k)*y)
            c1y = cos(cpsi(1,k)*y)
            s2y = sin(cpsi(2,k)*y)
            c2y = cos(cpsi(2,k)*y)
            
            c1  = cpsi(1,k)
            c2  = cpsi(2,k)

            if (k.eq.2) u = u + a(1,k)*s1x*c2y*c2
            if (k.eq.2) v = v - a(1,k)*c1x*s2y*c1
            if (k.eq.2) u = u - a(2,k)*s2x*c1y*c1
            if (k.eq.2) v = v + a(2,k)*c2x*s1y*c2

            if (k.eq.3) u = u - a(1,k)*s1x*c2y*c2
            if (k.eq.3) v = v + a(1,k)*c1x*s2y*c1
            if (k.eq.3) u = u - a(2,k)*c2x*c1y*c1
            if (k.eq.3) v = v - a(2,k)*s2x*s1y*c2

            if (k.eq.4) u = u + a(1,k)*c1x*c2y*c2
            if (k.eq.4) v = v + a(1,k)*s1x*s2y*c1
            if (k.eq.4) u = u + a(2,k)*c2x*c1y*c1
            if (k.eq.4) v = v + a(2,k)*s2x*s1y*c2

            if (k.eq.5) u = u - a(1,k)*s1x*c2y*c2
            if (k.eq.5) v = v + a(1,k)*c1x*s2y*c1
            if (k.eq.5) u = u - a(2,k)*s2x*c1y*c1
            if (k.eq.5) v = v + a(2,k)*c2x*s1y*c2
         enddo
         uu(i) = u*e + u0
         vv(i) = v*e + v0
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C
      udiff =0.
      utrans=0.
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C
      ffx = 0.0
      ffy = 0.0
      ffz = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C
      qvol   = 0.0
      source = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine userchk	
      include 'SIZE'  
      include 'TOTAL' 
c
      common /exacu/ ue(lx1,ly1,lz1,lelt),ve(lx1,ly1,lz1,lelt)
      common /exacd/ ud(lx1,ly1,lz1,lelt),vd(lx1,ly1,lz1,lelt)

      call itl() ! calls ITL

      ifield = 1  ! for outpost

      n    = nx1*ny1*nz1*nelv
      visc = param(2)
      u0   = param(96)
      v0   = param(97)
      call exact  (ue,ve,xm1,ym1,n,time,visc,u0,v0)
      if (istep.eq.0     ) call outpost(ue,ve,vx,pr,t,'   ')

      call sub3   (ud,ue,vx,n)
      call sub3   (vd,ve,vy,n)
      if (istep.eq.nsteps) call outpost(ud,vd,vx,pr,t,'   ')

      umx = glamax(vx,n)
      vmx = glamax(vy,n)
      uex = glamax(ue,n)
      vex = glamax(ve,n)
      udx = glamax(ud,n)
      vdx = glamax(vd,n)

      if (nid.eq.0) then
          write(6,11) istep,time,udx,umx,uex,u0,'  X err'
          write(6,11) istep,time,vdx,vmx,vex,v0,'  Y err'
   11     format(i5,1p5e14.6,a7)
      endif


      if (istep.le.5) then        !  Reset velocity to eliminate 
         call copy (vx,ue,n)      !  start-up contributions to
         call copy (vy,ve,n)      !  temporal-accuracy behavior.
      endif

      return
      end
c-----------------------------------------------------------------------
c
c calls ITL
c
      subroutine itl
      include 'SIZE'  
      include 'TOTAL' 

      integer first
      integer b

      ! ADD-BY-LEETEN 07/10/2011-BEGIN
	integer nvpb 	! #voxels/block
	integer bo	! block offset

        integer rv_vec_id ! ID of the random variable for the vector orientation
        save rv_vec_id
        data rv_vec_id /0/

        ! ADD-BY-LEETEN 07/22/2011-BEGIN
        integer rv_vecm_id ! ID of the random variable for the vector orientation
        save rv_vecm_id
        data rv_vecm_id /0/
        ! ADD-BY-LEETEN 07/22/2011-END

        integer time_step 
        save time_step 
        data time_step /0/

        integer time_step_mod

      ! ADD-BY-LEETEN 07/10/2011-END

c nelv = number of elements (blocks) per process
c nx1, ny1, nz1 = number of grid points per block in x,y,z
c vx, vy, vz = velocity data
c xm1, xm2, xm2 = geometry data
c velocity and geometry data ordered in (x,y,z,element) order, x changing fastest

c first time step
      save first
      data first /1/
      if (first.eq.1) then
          first = 0
          call ITL_begin()

          call ITL_add_random_field(nelv, 3, rf_id)
          call ITL_bind_random_field(rf_id)

          ! vector
          call ITL_add_random_variable(rv_vec_id)
          call ITL_bind_random_variable(rv_vec_id)
          ! MOD-BY-LEETEN 07/22/2011-FROM:
          ! call ITL_random_varable_as_vector3(1, 2, 3, 1) ! 1 mean using the vector orientation
          ! TO:
          ! MOD-BY-LEETEN 07/23/2011-FROM:
          ! call ITL_random_varable_as_vector2(1, 2, "dir")
          ! TO:
          call ITL_random_varable_as_vector3(1, 2, 3, "dir")
          ! MOD-BY-LEETEN 07/23/2011-END
          ! MOD-BY-LEETEN 07/22/2011-END

          ! ADD-BY-LEETEN 07/22/2011-BEGIN
          call ITL_add_random_variable(rv_vecm_id)
          call ITL_bind_random_variable(rv_vecm_id)
          ! MOD-BY-LEETEN 07/23/2011-FROM:
          ! call ITL_random_varable_as_vector2(1, 2, "abs")
          ! TO:
          call ITL_random_varable_as_vector3(1, 2, 3, "abs")
          ! MOD-BY-LEETEN 07/23/2011-END
          ! ADD-BY-LEETEN 07/22/2011-END

          do b = 1, nelv
             nvpb = nx1 * ny1 * nz1 
             bo = 1 + (b - 1) * nvpb

             ! specify the block size and geometry
             call ITL_bind_block(b)
             call ITL_block_size3(nx1, ny1, nz1)

             call ITL_geom_rect_dim_coord(1, xm1, bo, 1)
             call ITL_geom_rect_dim_coord(2, ym1, bo, nx1)
             call ITL_geom_rect_dim_coord(3, zm1, bo, nx1 * ny1)
             call ITL_dump_bound_block_geom_2tmp()
          enddo
      endif

c every 100 time steps
      time_step = time_step + 1
      
      time_step_mod = modulo(time_step, 100)
      if( time_step_mod.eq.1 ) then 
           if( 0.eq.1 ) then ! MOD-BY-LEETEN 07/22/2011-FROM:
	          do b = 1, nelv
	             nvpb = nx1 * ny1 * nz1 
	             bo = 1 + (b - 1) * nvpb
	
	             ! specfiy the data
	             call ITL_bind_block(b)
	             call ITL_bind_data_component(1) ! specify the U component
	             call ITL_data_source(vx, bo, 1) 
	             call ITL_bind_data_component(2) ! specify the V component 
	             call ITL_data_source(vy, bo, 1) 
	             call ITL_bind_data_component(3) ! specify the W component 
	             call ITL_data_source(vz, bo, 1) 
	
                     ! dump the feature vector
                     ! call ITL_dump_bound_block_feature_vector_2tmp(rv_vec_id)

                     ! compute and dump the entropy
                     ! call ITL_dump_bound_block_global_entropy_2tmp(rv_vec_id)

                     ! compute and dump the entropy
                     ! call ITL_dump_bound_block_global_entropy_2tmp(rv_vecm_id)
	          enddo

         else ! MOD-BY-LEETEN 07/22/2011-TO:

          do b = 1, nelv
             nvpb = nx1 * ny1 * nz1 
             bo = 1 + (b - 1) * nvpb

             ! specfiy the data
             call ITL_bind_block(b)
             call ITL_bind_data_component(1) ! specify the U component
             call ITL_data_source(vx, bo, 1) 
             call ITL_bind_data_component(2) ! specify the V component 
             call ITL_data_source(vy, bo, 1) 
             call ITL_bind_data_component(3) ! specify the W component 
             call ITL_data_source(vz, bo, 1) 
          enddo

          call ITL_use_domain_range(rv_vecm_id) ! obtain the range over the entire domain

          do b = 1, nelv
             ! specfiy the data
             call ITL_bind_block(b)

             ! dump the feature vector
             call ITL_dump_bound_block_feature_vector_2tmp(rv_vec_id)

             ! compute and dump the entropy
             call ITL_dump_bound_block_global_entropy_2tmp(rv_vec_id)

             ! compute and dump the entropy
             call ITL_dump_bound_block_global_entropy_2tmp(rv_vecm_id)
          enddo

          endif ! MOD-BY-LEETEN 07/22/2011-END
      endif

c last time step
      if (lastep.eq.1) then
         call ITL_end()
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      ux=0.0
      uy=0.0
      uz=0.0
      temp=0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      common /exacu/ ue(lx1,ly1,lz1,lelt),ve(lx1,ly1,lz1,lelt)
      common /exacd/ ud(lx1,ly1,lz1,lelt),vd(lx1,ly1,lz1,lelt)

      integer icalld
      save    icalld
      data    icalld  /0/

      n = nx1*ny1*nz1*nelv
      if (icalld.eq.0) then
         icalld = icalld + 1
         time = 0.
         u0   = param(96)
         v0   = param(97)
         call exact (ue,ve,xm1,ym1,n,time,visc,u0,v0)
      endif

      ie = gllel(ieg)
      ux=ue(ix,iy,iz,ie)
      uy=ve(ix,iy,iz,ie)
      uz=0.0
      temp=0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'
      integer e

      one   = 1.
      twopi = 8.*atan(one)

      do e=1,nelv   !  Rescale mesh to [0,2pi]^2
      do i=1,4      !  Assumes original domain in .rea file on [0,1]
         xc(i,e) = twopi*xc(i,e)
         yc(i,e) = twopi*yc(i,e)
      enddo 
      enddo 

      param(66) = 4
      param(67) = 4

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      return
      end
c-----------------------------------------------------------------------
c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
