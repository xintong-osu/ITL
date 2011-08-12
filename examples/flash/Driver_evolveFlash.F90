!!****if* source/Driver/DriverMain/Split/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This routine implements the Strang splitting scheme for time
!! advancement. A single step in the this driver 
!! includes two sweeps, the first one in order XYZ, and
!! the second one in order ZYX. This driver works with directionally
!! split operators only. The routine also controls the regridding of
!! the mesh if necessary and the simulation output.
!!
!!  
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_myPE or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif


subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_myPE, dr_numProcs, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_simGeneration, dr_fSweepDir, dr_rSweepDir,&
       dr_nstep, dr_dtOld, dr_dtNew, dr_restart, dr_elapsedWCTime, &
       dr_redshiftInitial, dr_redshiftFinal, dr_redshift, dr_redshiftOld, &
       dr_useRedshift
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Diffuse_interface, ONLY : Diffuse
  use Particles_interface, ONLY : Particles_advance, Particles_dump

!--- for ITL, add extra grid interface functions
!  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
!    Grid_getListOfBlocks, Grid_updateRefinement
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement, Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_getCellCoords, Grid_getBlkData
!---

  use Hydro_interface, ONLY : Hydro
  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface, ONLY :IO_output,IO_outputFinal
  use Cosmology_interface, ONLY : Cosmology_redshiftHydro, &
    Cosmology_solveFriedmannEqn, Cosmology_getRedshift

  implicit none

#include "constants.h"
#include "Flash.h"

!--- declarations for ITL
  ! ADD-BY-LEETEN 07/20/2011-BEGIN
  integer rf_id 
  integer rv_id
  ! ADD-BY-LEETEN 07/20/2011-END
  ! ADD-BY-LEETEN 08/12/2011-BEGIN
  integer time_step
  integer time_step_mod
  save time_step
  data time_step /0/
  ! ADD-BY-LEETEN 08/12/2011-END
  integer :: i, b, lb, nx, ny, nz
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  real, pointer :: solnData(:,:,:,:)
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: xs
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: ys
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: zs
  integer :: start(MDIM)
  real, allocatable :: data(:,:,:)
  integer :: size(3)
!---

!--- Matt's example
!integer :: iBlk, lb, lnblocks, imin, imax, jmin, jmax, kmin, kmax, &
!isize, jsize, ksize, i, j, k
!integer :: blkList(MAXBLOCKS)
!integer :: blkLimits(2, MDIM)
!integer :: blkLimitsGC(2, MDIM)
!real, dimension(GRID_ILO_GC:GRID_IHI_GC) :: xCoord
!real, dimension(GRID_JLO_GC:GRID_JHI_GC) :: yCoord
!real, dimension(GRID_KLO_GC:GRID_KHI_GC) :: zCoord
!real :: dens, x, y, z
!real, pointer :: solnData(:,:,:,:)
!---

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(4,2) :: strBuff
  character(len=15) :: numToStr

  logical :: gridChanged
  logical :: endRun !Should we end our run on this iteration?
  logical :: shortenedDt !Is the last timestep being shortened to reach dr_tmax?

  endRun = .false.

  call Logfile_stamp(dr_myPE, 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

!!******************************************************************************
!! Start of Evolution Loop
!!******************************************************************************

  do dr_nstep = dr_nbegin, dr_nend


     call dr_shortenLastDt(dr_dt, dr_simTime, dr_tmax, shortenedDt, 2)
     !!Step forward in time. See bottom of loop for time step calculation.
     
     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)
     if (dr_myPE == MASTER_PE) then

        write (numToStr(1:), '(I10)') dr_nstep
        write (strBuff(1,1), "(A)") "n"
        write (strBuff(1,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_simTime
        write (strBuff(2,1), "(A)") "t"
        write (strBuff(2,2), "(A)") trim(adjustl(numToStr))

        if (.not. dr_useRedshift) then

           write (numToStr(1:), "(1PE12.6)") dr_dt
           write (strBuff(3,1), "(A)") "dt"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           call Logfile_stamp(dr_myPE, strBuff(1:3,:), 3, 2, "step")

        else

           write (numToStr(1:), "(F8.3)") dr_redshift
           write (strBuff(3,1), "(A)") "z"
           write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))

           write (numToStr(1:), "(1PE12.6)") dr_dt
           write (strBuff(4,1), "(A)") "dt"
           write (strBuff(4,2), "(A)") trim(adjustl(NumToStr))
           
           call Logfile_stamp(dr_myPE, strBuff, 4, 2, "step")

        endif

     end if

     !!--------------------------------------------------------------------
     !!- Start Physics Sequence
     !!--------------------------------------------------------------------
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0) print*, 'going into Hydro/MHD'
#endif

     call Timers_start("cosmology")
     call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
     call Timers_stop("cosmology")

     dr_simTime = dr_simTime + dr_dt
     dr_simGeneration = 0

     call Timers_start("hydro")
#ifdef DEBUG_DRIVER
      if (dr_myPE == 0) print*,'going into hydro'
#endif
     call Hydro(dr_myPE, dr_numProcs, blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, dr_fSweepDir)

     call Timers_stop("hydro")

     
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Hydro/MHD timestep'
#endif


     call Diffuse(blockCount, blockList, dr_dt, pass=1)
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Diffuse '
#endif

     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt, pass=1)
     call Timers_stop("sourceTerms")
#ifdef DEBUG_DRIVER
      if (dr_myPE == 0) print*,'done source terms'
    if (dr_myPE == 0)   prin*, 'return from Drivers_sourceTerms '
#endif
     call Timers_start("Particles_advance")
     call Particles_advance(dr_dtOld, dr_dt)
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Particles_advance '
#endif
     call Timers_stop("Particles_advance")
   
     call Gravity_potentialListOfBlocks(blockCount,blockList)
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Gravity_potential '
#endif

     call Timers_start("cosmology")
     call Cosmology_redshiftHydro(dr_myPE, dr_numprocs, blockCount, blockList)
     call Timers_stop("cosmology")

#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from redshiftHydro '
#endif
!!******************************************************************************
!!Second "half-step" of the evolution loop
!!******************************************************************************

#ifdef DEBUG_DRIVER
     if (dr_myPE == 0) print*, 'going into Hydro/MHD'
#endif
     call Timers_start("cosmology")
     call Cosmology_solveFriedmannEqn(dr_simTime, dr_dt)
     call Timers_stop("cosmology")

     dr_simTime = dr_simTime + dr_dt
     dr_simGeneration = 0
     call Timers_start("hydro")
#ifdef DEBUG_DRIVER
      if (dr_myPE == 0) print*,'going into hydro'
#endif
     call Hydro(dr_myPE, dr_numProcs, blockCount, blockList, &
                dr_simTime, dr_dt, dr_dtOld, dr_rSweepDir)
     call Timers_stop("hydro")
  
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Hydro/MHD timestep'
#endif

     call Diffuse(blockCount, blockList, dr_dt, pass=2)

     call Timers_start("sourceTerms")
     call Driver_sourceTerms(blockCount, blockList, dr_dt)
     call Timers_stop("sourceTerms")

#ifdef DEBUG_DRIVER
      if (dr_myPE == 0) print*,'done source terms'
    if (dr_myPE == 0)   print*, 'return from Drivers_sourceTerms '
#endif
     call Timers_start("Particles_advance")
     call Particles_advance(dr_dt, dr_dt)
     call Timers_stop("Particles_advance")
     
#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Particles_advance '
#endif
     call Gravity_potentialListOfBlocks(blockCount,blockList)

#ifdef DEBUG_DRIVER
     if (dr_myPE == 0)  print*, 'return from Gravity_potential '
#endif
     call Timers_start("cosmology")
     call Cosmology_redshiftHydro(dr_myPE, dr_numprocs, blockCount, blockList)
     call Timers_stop("cosmology")

     !--------------------------------------------------------------------
     !- End Physics Sequence -- Start Simulation Bookkeeping
     !--------------------------------------------------------------------

     !output a plotfile before the grid changes
     !call Timers_start("IO_output")
     !call IO_output(dr_myPE, dr_numProcs, dr_simTime, &
     !     dr_dt, dr_nstep+1, dr_nbegin, endRun, PLOTFILE_AND_PARTICLEFILE)
     !call Timers_stop("IO_output")


     !!if (itemp_limit) .eq. 1) call Hydro_timstepPrecompute()

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement(dr_myPE, dr_nstep, dr_simTime, gridChanged)
     call Timers_stop("Grid_updateRefinement")

     if (gridChanged) dr_simGeneration = dr_simGeneration + 1

     dr_dtOld = dr_dt                     ! backup needed old 
     ! calculate new
     
     call Timers_start("compute dt")
     call Driver_computeDt(dr_myPE, dr_numProcs,  &
                         dr_nbegin, dr_nstep, &
                         dr_simTime, dr_dtOld, dr_dtNew)
     call Timers_stop("compute dt")
     dr_dt = dr_dtNew                                        ! store new
     

     !!-----------------------------------------------------------------
     !! Output for current step in evolution
     !!-----------------------------------------------------------------

     call Timers_start("IO_output")
     call IO_output(dr_myPE, dr_numProcs, dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun, CHECKPOINT_FILE_ONLY)
     call Timers_stop("IO_output")

!!*****************************************************************************
!!  Evolution Loop -- check termination conditions
!!*****************************************************************************


     !Exit if this step was handled specially as the last step
     if(shortenedDt) exit
     !Exit if a .dump_restart or .kill was found during the last step
     if(endRun) exit

     !call Particles_dump(dr_myPE, blockCount, blockList, dr_nstep, dr_simTime, dr_dt)


     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift  
     !!        (also called redshiftFinal) 
     !!  (iii) the wall clock time is greater than the maximum 
     !!        (wall_clock_time_max)

     !!Update redshift from Driver's POV.  Need this for exit condition. -PR
     !!old redshift needed for accurate restarts.
     dr_redshiftOld = dr_redshift
     call Cosmology_getRedshift(dr_redshift)
     
     if (dr_simTime >= dr_tmax) then
        if(dr_myPE == MASTER_PE) then
           print *, "exiting: reached max SimTime"
        end if
        exit
     end if
     
     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_elapsedWCTime >  dr_wallClockTimeLimit) then
        if(dr_myPE == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        end if
        exit
     end if

     if (dr_redshift < dr_redshiftfinal .and. dr_useRedshift) then
        if(dr_mype == MASTER_PE) then
           print *, "exiting: reached redshiftfinal"
        end if
        exit
     end if

!--- call ITL here
     ! assuming all blocks are same number of cells
     ! get size and allocate data space
     call Grid_getBlkIndexLimits(blockList(1), blkLimits, blkLimitsGC)
     nx = blkLimits(HIGH, IAXIS) - blkLimits(LOW, IAXIS) + 1
     ny = blkLimits(HIGH, JAXIS) - blkLimits(LOW, JAXIS) + 1
     nz = blkLimits(HIGH, KAXIS) - blkLimits(LOW, KAXIS) + 1
     start(1) = 1
     start(2) = 1
     start(3) = 1
     size(1) = nx
     size(2) = ny
     size(3) = nz
     allocate(data(nx, ny, nz))

     if( time_step.lt.1 ) then	  ! ADD-BY-LEETEN 08/12/2011
     call ITL_add_random_field(blockCount, 1, rf_id)
     call ITL_bind_random_field(rf_id)

     ! scalar
     call ITL_add_random_variable(rv_id)
     call ITL_bind_random_variable(rv_id)

     call ITL_rv_name("d ")	  ! ADD-BY-LEETEN 08/12/2011
     ! MOD-BY-LEETEN 07/20/2011-FROM:
     ! call ITL_random_varable_as_scalar(1) ! 1 mean using the vector orientation
     ! TO:
     call ITL_random_varable_as_scalar(1, "raw") 
     ! ADD-BY-LEETEN 07/31/2011-BEGIN
     call ITL_set_n_bins(16) 
     ! ADD-BY-LEETEN 07/31/2011-END
     ! MOD-BY-LEETEN 07/20/2011-END

     do b = 1, blockCount
        lb = blockList(b)
        call Grid_getBlkData(lb, CENTER, PDEN_VAR, INTERIOR, start, data, &
             size)
        call Grid_getCellCoords(IAXIS, lb, CENTER, .false., xs, nx)
        call Grid_getCellCoords(JAXIS, lb, CENTER, .false., ys, ny)
        call Grid_getCellCoords(KAXIS, lb, CENTER, .false., zs, nz)

	! specify the block size and geometry
	call ITL_bind_block(b)
	call ITL_block_size3(nx, ny, nz)

	call ITL_geom_rect_dim_coord(1, xs, 1, 1)
	call ITL_geom_rect_dim_coord(2, ys, 1, 1)
	call ITL_geom_rect_dim_coord(3, zs, 1, 1)
	! DEL-BY-LEETEN 08/12/2011-BEGIN	call ITL_dump_bound_block_geom_2tmp()
     ! ADD-BY-LEETEN 08/12/2011-BEGIN
     enddo
     call ITL_bind_data_component(1) 	! specify the data component
     call ITL_data_name("d ")

     call ITL_nc_create()
     call ITL_nc_wr_geom()
     endif

     time_step = time_step + 1;
     if( time_step.le.100 ) then
	call ITL_set_time_stamp(time_step)
        do b = 1, blockCount
           call ITL_bind_block(b)
     ! ADD-BY-LEETEN 08/12/2011-END

        ! specify the data
	call ITL_bind_data_component(1) 	! specify the data component
	call ITL_data_source(data, 1, 1) 	! specify the array to the data
     ! ADD-BY-LEETEN 07/20/2011-BEGIN
     enddo

     ! Find the range of the entire domain 
     call ITL_use_domain_range(rv_id)

     ! ADD-BY-LEETEN 08/12/2011-BEGIN
     call ITL_nc_wr_data()
     call ITL_nc_wr_rv(rv_id)
     ! ADD-BY-LEETEN 08/12/2011-END

     do b = 1, blockCount
	call ITL_bind_block(b)
     ! ADD-BY-LEETEN 07/20/2011-END

	call ITL_dump_bound_block_feature_vector_2tmp(rv_id)

        ! compute and dump the entropy
	call ITL_dump_bound_block_local_entropy3_2tmp(rv_id, REAL(4), REAL(4), REAL(4))
     enddo
     endif	     ! ADD-BY-LEETEN 08/12/2011

     deallocate(data)
!---

!--- Matt's example
!call Grid_getListOfBlocks(LEAF, blkList, lnblocks)

!do iBlk = 1, lnblocks

! lb = blkList(iBlk)

! call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
! imin = blkLimits(LOW,  IAXIS)
! imax = blkLimits(HIGH, IAXIS)
! jmin = blkLimits(LOW,  JAXIS)
! jmax = blkLimits(HIGH, JAXIS)
! kmin = blkLimits(LOW,  KAXIS)
! kmax = blkLimits(HIGH, KAXIS)

! isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
! jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
! ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

! call Grid_getBlkPtr(lb, solnData)

! call grid_GetCellCoords(IAXIS, lb, CENTER, .true., xCoord, isize)
! call grid_GetCellCoords(JAXIS, lb, CENTER, .true., yCoord, jsize)
! call grid_GetCellCoords(KAXIS, lb, CENTER, .true., zCoord, ksize)

! do k = kmin, kmax
! do j = jmin, jmax
! do i = imin, imax

!   x = xCoord(i)
!   y = yCoord(j)
!   z = zCoord(k)
!   dens = solnData(PDEN_VAR,i,j,k)

!   if (dens > 5.E-30) then
!      print*, "The location of this cell is", x,y,z, "and its density is", dens
!   endif

! end do
! end do
! end do

! call Grid_releaseBlkPtr(lb, solnData)
!end do
!---

  enddo
  !The value of dr_nstep after the loop is (dr_nend + 1) if the loop iterated for
  !the maximum number of times.  However, we need to retain the value that
  !dr_nstep had during the last loop iteration, otherwise the number for nstep
  !that will be stored in a final checkpoint file will be wrong.
  dr_nstep = min(dr_nstep,dr_nend)

!!******************************************************************************
!! End of Evolution Loop
!!******************************************************************************

  call Timers_stop("evolution")

  call Logfile_stamp(dr_myPE, 'Exiting evolution loop' , '[Driver_evolveFlash]')

  !if a file termination, this may already be done.
  if(.NOT.endRun) call IO_outputFinal(dr_myPE, dr_numProcs)

  call Timers_getSummary(dr_myPE, max(0,dr_nstep-dr_nbegin+1))


  call Logfile_stamp(dr_myPE, "FLASH run complete.", "LOGFILE_END")

  call Logfile_close(dr_myPE)


  return
  
end subroutine Driver_evolveFlash



