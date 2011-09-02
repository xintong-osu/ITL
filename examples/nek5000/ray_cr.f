c-----------------------------------------------------------------------
c
c     A short example of Rayleigh-Benard convection.
c
c     This case specific to the Ra_c determination, with a very short
c     domain (Lx := 2pi/k_c) and using only 3 elements in the x-direction
c     to capture a single roll-pair at the prescribed critical wavenumber, k_c.
c
c     Ra_c is determined by running two or more cases at Ra > Ra_c and computing
c     the kinetic energy levels, Ek, once the flow has reached a steady state.
c     For two or more values of Ra close to but greater than Ra_c, one can
c     linearly extrapolate back to Ek=0 and thus estimate Ra_c.   This is
c     done in userchk.
c
c     Parameters are set in routine rayleigh_const.
c
c     For the given nondimensionalization, set rho==1 (parameter p1 in .rea
c     file) and visocity (p2) to be the desired Prandtl number.
c
c     Rayleigh number is set as Ra=p76 and fundamental wavenumber kc=p75.
c
c     The buoyancy is ffy = Ra Pr T, where T is determined by 
c     boundary and initial conditions.
c
c     For the Dirichlet-Dirichlet (dd) case, the critical Rayleigh number is 
c     around 1707 for wavenumber k_c = 3.117 (Chandrasekhar 1961 "Hydrodynamic
c     and Hydromagnetic Stability", Table III on p. 43)
c
c     GEOMETRY:
c
c     There are two primary cases ray_##.box and ray_9.box.  The former 
c     specifies 3 elements in x, with ##=dd, dn, or nn, according to the
c     choice of Dirichlet-Direchlet, Dirchlet-Neumann, or Neumann-Neumann
c     boundary conditions. The latter specifies 9 elements in x. All cases 
c     have just a single element in y, corresponding to an Nth-order spectral
c     method.  For the Ra_c cases, the length in x is rescaled in usrdat2()
c     to [0:2pi/k_c].   The ray_9 case is designed to illustrate a typical
c     Rayleigh-Benard simulation, while the 3-element cases are set up to
c     determine Ra_c as quickly as possible.  The Ra_c estimates are typically
c     good to 5 digits, which is about as much as can be expected given that
c     k_c is specified to only 4 digits.
c
c
c     NOTES:
c
c     A time trace of volume-avearged kinetic energy Ek vs t is output to
c     the logfile. 
c
c     For long domains, be careful about selecting an even number of elements 
c     in x as it appears that the RB system likes to lock onto the grid spacing
c     and give a number of rolls that matches the number of elements, if the
c     elements have order-unity aspect ratio, as in the present case.
c     Thus, in this case, a 9 element mesh is likely to be more faithful
c     to the linear stability theory, at least for modest polynomial orders
c     of lx1=8.
c     
c     It appears that one cannot realize Courant conditions of CFL ~ 0.5
c     with these cases because of the explicit Boussinesq treatment.
c     The given value dt=.02 is stable with lx1=8.
c
c     Use 'grep rayl *.log' to get critical Rayleigh number estimates
c     for each of the different boundary condition cases.
c
c
c-----------------------------------------------------------------------
      subroutine rayleigh_const

      include 'SIZE'
      include 'INPUT'

      common /rayleigh_r/ rapr,ta2pr
      common /rayleigh_c/ Ra,Ra1,Ra2,Ra3,Prnum,Ta2,Ek1,Ek2,Ek3,ck

      Prnum  = param(2)	! Pr=1

      ck     = param(75)
      Ra     = param(76)   ! Ball-park estimate for Ra_c
      Ta2    = param(77)   ! Taylor number squared (=0 for this case)

      Ra3 = Ra
      Ra2 = Ra3 + 10
      Ra1 = Ra2 + 10

      Ra  = Ra1       ! Start with largest Rayleigh number

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      udiff  = 0
      utrans = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      common /rayleigh_r/ rapr,ta2pr

      buoy = temp*rapr

      if (if3d) then
         ffx  =   uy*Ta2Pr
         ffy  = - ux*Ta2Pr
         ffz  = buoy
      elseif (ifaxis) then
         ffx  = -buoy
         ffy  =  0.
      else
         ffx  = 0.
         ffy  = buoy
      endif
c     write(6,*) ffy,temp,rapr,'ray',ieg

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol   = 0.0
      source = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'NEKUSE'
      common /rayleigh_r/ rapr,ta2pr

      ux=0.
      uy=0.
      uz=0.

      temp=0.  !     Temp = 0 on top, 1 on bottom

      if (if3d) then
         temp = 1-z
      elseif (ifaxis) then  !      domain is on interval x in [-1,0]
         temp = 1.+x
      else                  ! 2D:  domain is on interval y in [0,1]
         temp = 1.-y
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer idum
      save    idum
      data    idum /99/

      ran = 2.e7*(ieg+x*sin(y)) + 1.e6*ix*iy + 1.e7*ix 
      ran = 1.e9*sin(ran)
      ran = 1.e9*sin(ran)
      ran = cos(ran)
      ran = ran1(idum)
      amp = .005

      temp = 1-y + ran*amp*(1-y)*y*x  ! 2D 

      ux=0.0
      uy=0.0
      uz=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

      common /rayleigh_c/ Ra,Ra1,Ra2,Ra3,Prnum,Ta2,Ek1,Ek2,Ek3,ck

      call rayleigh_const


      param(66) = 4
      param(67) = 4

      one = 1.
      pi  = 4.*atan(one)

      x0 = 0.
      x1 = 2*pi/ck
      zero = 0.
      call rescale_x(xm1,x0,x1)

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      common /rayleigh_r/ rapr,ta2pr
      common /rayleigh_c/ Ra,Ra1,Ra2,Ra3,Prnum,Ta2,Ek1,Ek2,Ek3,ck

      real Ek0,Ek,t12,t23
      save Ek0,Ek,t12,t23
      integer kstep,icase,i12,i23
      save    kstep,icase,i12,i23
      data    kstep,icase /0,0/

      ifxyo = .true.  ! For VisIt
      if (istep.gt.iostep) ifxyo = .false.

      rapr    = ra*prnum
      ta2pr   = ta2*prnum

      n    = nx1*ny1*nz1*nelv
      Ek0  = Ek
      Ek   = .5*(glsc3(vx,vx,bm1,n)+ glsc3(vy,vy,bm1,n))/volvm1

      sigma = 1.e-4
      de    = abs(Ek-Ek0)/dt

      if (nid.eq.0.and.(mod(istep,100).eq.0.or.istep.lt.100))
     $   write(6,6) kstep,istep,time,Ek,de,ra,icase
    6 format(2i7,1p4e13.5,i3,' ekde')

      kstep = kstep+1
      if (kstep.gt.1000.and.de.lt.sigma*Ek) then

         if (icase.eq.0) then  ! First pass

            Ek1   = Ek
            icase = icase+1    ! Set up 2nd pass
            kstep = 0
            Ra    = Ra2
            t12   = time
            i12   = istep

         elseif (icase.eq.1) then  ! 2nd pass

            Ek2   = Ek
            icase = icase+1    ! Set up 3rd pass
            kstep = 0
            Ra    = Ra3
            t23   = time
            i23   = istep


         elseif (icase.eq.2) then  ! Make estimate and exit
            Ek3 = Ek
            r12 = Ra1 - (Ra2-Ra1)*Ek1/(Ek2-Ek1)
            r23 = Ra2 - (Ra3-Ra2)*Ek2/(Ek3-Ek2)

            if (r23.gt.600.) Ra_c = 657.511 ! Values from Chandrashekar
            if (r23.gt.1000) Ra_c = 1100.65 
            if (r23.gt.1500) Ra_c = 1707.76 
            
            if (nid.eq.0) write(6,1) r12,r23,Ra_c,ck,ra1,ra2,ra3
  1         format(3f9.3,f9.4,3f8.1,' rayleigh')

            if (nid.eq.0) write(6,2) i12,i23,t12,t23
            call ITL_end() ! ADD-BY-LEETEN 09/01/2011
  2         format(2i8,2f9.4,' transitions')
            call exitti('Done in userchk.$',istep)

         endif
      endif

	! ADD-BY-LEETEN 08/05/2011-BEGIN
      	call itl()                ! calls ITL
	! ADD-BY-LEETEN 08/05/2011-END

      return
      end
c-----------------------------------------------------------------------

! ADD-BY-LEETEN 07/22/2011-BEGIN
c
c calls ITL
c
      subroutine itl
      include 'SIZE'  
      include 'TOTAL' 

      integer first
      integer b

	integer nvpb 	! #voxels/block
	integer bo	! block offset

		! ADD-BY-LEETEN 07/22/2011-BEGIN
        integer rf_id           ! ID of the random field
        integer rv_t_id         ! ID of the random variable for temperature
        integer rv_vec_id       ! ID of the random variable for vector

        save rv_t_id
        data rv_t_id /0/

        save rv_vec_id
        data rv_vec_id /0/

        integer time_step
        integer time_step_mod
        save time_step
        data time_step /0/
     
		! ADD-BY-LEETEN 07/22/2011-END

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

	  ! MOD-BY-LEETEN 08/05/2011-FROM:
          	! call ITL_begin()
	  ! TO:
	  ! read the session's name
	ierr = 0
	IF(NID.EQ.0) THEN
        	OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD',ERR=24)
	        READ(8,10) SESSION
        	READ(8,10) PATH
 10		FORMAT(A132)
		CLOSE(UNIT=8)
        	GOTO 23
 24     	ierr = 1
 23   	ENDIF

	call err_chk(ierr,' Cannot open SESSION.NAME!$')

	! send the session's name to all process and then call ITL
	call bcast(SESSION,132*CSIZE)

	! initialize ITL by sending the name of this test
	call ITL_begin(132, SESSION) 
	! MOD-BY-LEETEN 08/05/2011-END

          ! create a random field of nelv blocks and 4 data components, which will be the temperature and the 3D vectors
          call ITL_add_random_field(nelv, 4, rf_id)
          call ITL_bind_random_field(rf_id)

          ! ADD-BY-LEETEN 08/12/2011-BEGIN
          call ITL_set_local2global_mapping(lglel, 1) ! 1: the IDs are 1-based
          ! ADD-BY-LEETEN 08/12/2011-END

          ! temperature
          call ITL_add_random_variable(rv_t_id)
          call ITL_bind_random_variable(rv_t_id)
          ! ADD-BY-LEETEN 08/06/2011-BEGIN
          call ITL_rv_name("temperature ")
          ! ADD-BY-LEETEN 08/06/2011-END

#if 0 
! MOD-BY-LEETEN 08/05/2011-FROM:
          ! MOD-BY-LEETEN 07/22/2011-FROM:
          ! call ITL_random_varable_as_scalar(1)
          ! TO:
	  call ITL_random_varable_as_scalar(1, "abs")
          ! MOD-BY-LEETEN 07/22/2011-END
          call ITL_set_random_variable_range(-2.0, +2.0)
#else
! MOD-BY-LEETEN 08/05/2011-TO:
	  call ITL_random_varable_as_scalar(1, "raw")
          ! MOD-BY-LEETEN 07/22/2011-END
          call ITL_set_random_variable_range(0.0, 1.0)
#endif
! MOD-BY-LEETEN 08/05/2011-END

	  ! ADD-BY-LEETEN 07/31/2011-BEGIN
          call ITL_set_n_bins(16)
	  ! ADD-BY-LEETEN 07/31/2011-END

          ! vector
          call ITL_add_random_variable(rv_vec_id)
          call ITL_bind_random_variable(rv_vec_id)
          ! ADD-BY-LEETEN 08/06/2011-BEGIN
          call ITL_rv_name("vec ")
          ! ADD-BY-LEETEN 08/06/2011-END

          ! MOD-BY-LEETEN 07/22/2011-FROM: 
          ! call ITL_random_varable_as_vector3(2, 3, 4, 1) ! 1 mean using the vector orientation as the random variable
          ! TO:
          call ITL_random_varable_as_vector2(2, 3, "dir")
          ! MOD-BY-LEETEN 07/22/2011-END

		  ! ADD-BY-LEETEN 07/31/2011-BEGIN
          call ITL_set_n_bins(16)
		  ! ADD-BY-LEETEN 07/31/2011-END

         do b = 1, nelv
            call ITL_bind_block(b)
            call ITL_block_size3(nx1, ny1, nz1)

         ! ADD-BY-LEETEN 08/06/2011-BEGIN
         enddo

         call ITL_bind_data_component(1)
         call ITL_data_name('temperature ')
         call ITL_bind_data_component(2)
         call ITL_data_name('u ')
         call ITL_bind_data_component(3)
         call ITL_data_name('v ')
         call ITL_bind_data_component(4)
         call ITL_data_name('w ')

         ! Create the NetCDF.
         ! In order to finish the define mode,
         ! the NetCDF file must be created after
         ! 1. The #blocks and #data components are known
         ! 2. The dimension of the data block are given
         ! Besides, for accuracy and readability, 
         ! it is recommended to set the name of the data components
         ! and random variable.
         call ITL_nc_create()

         do b = 1, nelv
            call ITL_bind_block(b)
         ! ADD-BY-LEETEN 08/06/2011-END

            nvpb = nx1 * ny1 * nz1 
            bo = 1 + (b - 1) * nvpb

            ! specify and dump the block geometry
            call ITL_geom_rect_dim_coord(1, xm1, bo, 1)
            call ITL_geom_rect_dim_coord(2, ym1, bo, nx1)
            call ITL_geom_rect_dim_coord(3, zm1, bo, nx1 * ny1)
            ! DEL-BY-LEETEN 08/06/2011  call ITL_dump_bound_block_geom_2tmp()
         enddo

         ! ADD-BY-LEETEN 08/06/2011-BEGIN
         call ITL_nc_wr_geom()
         ! ADD-BY-LEETEN 08/06/2011-END
      endif

c every 10 time step
      time_step = time_step + 1

      ! MOD-BY-LEETEN 08/07/2011-FROM:
        ! time_step_mod = modulo(time_step, 10)
        ! if( time_step_mod.eq.1 ) then
      ! MOD-BY-LEETEN 08/07/2011-TO:
      ! for 
      if( time_step.le.100 ) then
      ! MOD-BY-LEETEN 08/07/2011-END
	! ADD-BY-LEETEN 08/05/2011-BEGIN
	call ITL_set_time_stamp(time_step)
	! ADD-BY-LEETEN 08/05/2011-END

         do b = 1, nelv
            call ITL_bind_block(b)

            nvpb = nx1 * ny1 * nz1 
            bo = 1 + (b - 1) * nvpb

            call ITL_bind_data_component(1) ! specify the temperature
            ! MOD-BY-LEETEN 08/05/2011-FROM: 
		!call ITL_data_source(temp, bo, 1)
            ! TO:
            call ITL_data_source(t, bo, 1)
            ! MOD-BY-LEETEN 08/05/2011-END
            call ITL_bind_data_component(2) ! specify the U component
            call ITL_data_source(vx, bo, 1) 
            call ITL_bind_data_component(3) ! specify the V component 
            call ITL_data_source(vy, bo, 1) 
            call ITL_bind_data_component(4) ! specify the W component 
            call ITL_data_source(vz, bo, 1) 
         enddo

         ! ADD-BY-LEETEN 08/06/2011-BEGIN
         call ITL_nc_wr_data()
         call ITL_nc_wr_rv(rv_t_id)
         call ITL_nc_wr_rv(rv_vec_id)
         ! ADD-BY-LEETEN 08/06/2011-END

         do b = 1, nelv
            call ITL_bind_block(b)

            ! compute and dump the feature vector/entropy for temperature

            ! DEL-BY-LEETEN 08/06/2011  call ITL_dump_bound_block_feature_vector_2tmp(rv_t_id) 
            call ITL_dump_bound_block_global_entropy_2tmp(rv_t_id)

            ! compute and dump the feature vector/entropy for the 3D vector field
            ! DEL-BY-LEETEN 08/06/2011  call ITL_dump_bound_block_feature_vector_2tmp(rv_vec_id)
            call ITL_dump_bound_block_global_entropy_2tmp(rv_vec_id)
         enddo
      endif

c last time step
      if (lastep.eq.1) then
          call ITL_end()
      endif

      return
      end
! ADD-BY-LEETEN 07/22/2011-END

c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
