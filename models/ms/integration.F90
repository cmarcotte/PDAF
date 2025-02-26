!$Id$
!BOP
!
! !ROUTINE: integration  --- integration routine for the Lorenz96 model
!
! !INTERFACE:
SUBROUTINE integration(time, nsteps)

! !DESCRIPTION:
! Routine to perform model integration with the Lorenz96 model.
!
! For simplicity of the implementation with PDAF,
! the time stepping is separated into a single routine.
! This allows to simply put the assimilation routine
! assimilation\_pdaf() in between the main program and
! the integration routine. If the time stepping is not
! available as a separate routine, a different
! implementation style is required.
!
! !REVISION HISTORY:
! 2009-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &            ! Timing
       ONLY: timeit
  USE mod_model, &        ! Model variables
       ONLY: x, dt, dim_state
#ifdef USE_PDAF
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype, incremental, model_error
#endif
  USE output_netcdf, &    ! NetCDF output
       ONLY: write_netcdf, close_netcdf

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(inout) :: time   ! Model time
  INTEGER, INTENT(in) :: nsteps ! Number of time steps to be performed
!EOP

! local variables
  INTEGER :: step               ! Time step counter
  REAL, ALLOCATABLE :: x1(:), x2(:), x3(:), x4(:) ! Temporary arrays for RK4
  REAL, ALLOCATABLE :: x_tmp(:)

#ifdef USE_PDAF
  EXTERNAL :: distribute_stateinc_pdaf ! Routine to add state increment for IAU
#endif


! **********************
! *** Initialization ***
! **********************

  ! Allocate temporary arrays for RK4
  ALLOCATE(x1(dim_state))
  ALLOCATE(x2(dim_state))
  ALLOCATE(x3(dim_state))
  ALLOCATE(x4(dim_state))
  ALLOCATE(x_tmp(dim_state))


! *********************************
! *** Perform model integration ***
! *********************************

  CALL timeit(5, 'new')

! *** time stepping loop ***
  integrate: DO step = 1, nsteps

#ifdef USE_PDAF
     ! For incremental updating (SEEK, SEIK, and LSEIK)
     IF (incremental == 1 &
          .AND. (filtertype==0 .OR. filtertype == 1 .OR. filtertype == 3)) THEN
        CALL PDAF_incremental(nsteps, distribute_stateinc_pdaf)
     END IF
#endif

! *** model time step - RK4 ***

     ! Intermediate steps
     CALL mitsch_dxdt(dim_state, x, x1, time)
     x1 = dt * x1

     x_tmp = x + x1/2.0
     CALL mitsch_dxdt(dim_state, x_tmp, x2, time)
     x2 = dt * x2

     x_tmp = x + x2/2.0
     CALL mitsch_dxdt(dim_state, x_tmp, x3, time)
     x3 = dt * x3

     x_tmp = x + x3
     CALL mitsch_dxdt(dim_state, x_tmp, x4, time)
     x4 = dt * x4

     ! New value of x
     x = x + x1/6.0 + x2/3.0 + x3/3.0 + x4/6.0

     ! Increment time
     time = time + dt


#ifdef USE_PDAF
     ! *** PDAF: Add model error ***
     IF (model_error) CALL add_model_noise(dt, dim_state, x)
#endif

#ifndef USE_PDAF
     ! Write NetCDF output
     CALL write_netcdf(step, time, dim_state, x)
#endif

  END DO integrate

#ifndef USE_PDAF
  ! Close NetCDF file
  CALL close_netcdf()
#endif

  DEALLOCATE(x1, x2, x3, x4, x_tmp)

  CALL timeit(5, 'old')

END SUBROUTINE integration

! !ROUTINE: mitsch_dxdt  --- compute dx/dt for Lorenz96 model
!
! !INTERFACE:
SUBROUTINE mitsch_dxdt(dim_state, x, dxdt, time)

! !DESCRIPTION:
! This function computes the time derivate for the Mitchell-Schaeffer model
!
! !REVISION HISTORY:
! 2009-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: forcing
!  #ifdef USE_PDAF
!  USE mod_assimilation, & ! Variables for assimilation
!      ONLY: filtertype, incremental, model_error
!  #endif
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_state     ! Dimension of model state
  REAL, INTENT(in)  :: x(dim_state)    ! Model state
  REAL, INTENT(out) :: dxdt(dim_state) ! Time derivate
  REAL, INTENT(in)  :: time            ! Model time
!EOP

! local variables
  INTEGER :: i, im1, ip1, np, nx, nv
  REAL :: du, dv, lu, d, ti, tu, to, tc, u, v, ug, k, H, Is
  REAL, external :: sigmoid, stimulus, envelope

! constant parameters
  np = 7						! number of parameters
  nx = ((dim_state - np)/2)		! number of grid points
  nv = 2						! number of fields

! state is laid out as:
! x(1:nx) = u,
! x(nx+1:nv*nx) = v,
! x(nv*nx+1:dim_state) = \(d, k, ti, tu, to, tc, ug\)
 d  = x(dim_state - np + 1)
 k  = x(dim_state - np + 2)
 ti = x(dim_state - np + 3)
 tu = x(dim_state - np + 4)
 to = x(dim_state - np + 5)
 tc = x(dim_state - np + 6)
 ug = x(dim_state - np + 7)

! Compute derivative
  DO i = 1, nx
     ! Cyclic boundary conditions:
     im1 = i - 1
     IF(im1 < 1) im1 = 1	!nx
     ip1 = i + 1
     IF(ip1 > nx) ip1 = nx	!1
     ! select local i => u,v => H, Is
     u = x(i)
     v = x(i + nx)
     H = sigmoid(u - ug, k)
     Is = stimulus(i, 10, time, 100.0)
     ! Laplacian
     lu = x(ip1) + x(im1) - 2*u
     ! change in u, v
     du = v*u*u*(1 - u)/ti - u/tu
     dv = (1 - H)*(1 - v)/to - H * v/tc
     ! Compute derivative
     dxdt(i) = d*lu + du + forcing * Is
     dxdt(i + nx) = dv
  END DO
! parameter evolution is 0, identically
  DO i = nv*nx+1, dim_state
  	dxdt(i) = 0
  END DO

#ifdef USE_PDAF
 	! *** PDAF: Add parameter error ***
  	CALL add_model_noise(forcing, 7, dxdt((dim_state - 7 + 1):dim_state))
#endif
END SUBROUTINE mitsch_dxdt

! function for smoothed sigmoid
! !ROUTINE: sigmoid  --- compute smoothed sigmoid
REAL FUNCTION sigmoid(x, k) RESULT(s)
  REAL, INTENT(in) :: x, k
  s = 0.5*(1.0 + TANH(k*x))
END FUNCTION sigmoid

! function for periodic forcing stimulus
! !ROUTINE: stimulus  --- compute forcing stimulus
REAL FUNCTION stimulus(i, l, t, p) RESULT(s)
	INTEGER, INTENT(in) :: i, l
  	REAL, INTENT(in) :: t, p
  	s = EXP(-ABS(REAL(i)/REAL(l)))
    s = s * envelope(t/p)
END FUNCTION stimulus

! function for spectrally smoothed envelope
! !ROUTINE: envelope  --- compute smoothed envelope
REAL FUNCTION envelope(a) RESULT(s)
	REAL, INTENT(in) :: a
  	s = COS(2.0*PI*a)**500
END FUNCTION envelope
