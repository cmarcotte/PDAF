!$Id$
!BOP
!
! !ROUTINE: initialize  --- initialize the Lorenz96 model
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the Lorenz 96 model.
! The model is introduced in E. N. Lorenz, Predictability: A problem
! partly solved. Proc. Seminar on Predictability, Vol. 1, ECMWF,
! Reading, UK, 1-18 and also described in E. N. Lorenz and K. A. Emanuel
! Optimal Sites for Supplementary Weather Observations: Simulation
! with a Small Model, J. Atm. Sci. 55 (1998) 399-414.
!
! Here the model is implemented using a 4th order Runge-Kutta scheme.
! Because this model is typically used with small dimensions, it is
! implemented here ithout parallelization.
!
! !REVISION HISTORY:
! 2009-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: dim_state, dt, step_null, step_final, x, forcing
  USE mod_modeltime, &    ! Time information for model integration
       ONLY: time, total_steps
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, mype_world, init_parallel, finalize_parallel
  USE mod_memcount, &     ! Counting allocated memory
       ONLY: memcount, memcount_ini, memcount_get
  USE parser, &           ! Parse command lines
       ONLY: handle, parse
  USE output_netcdf, &    ! NetCDF output
       ONLY: file_state, delt_write, init_netcdf

  IMPLICIT NONE
  INTEGER :: i
! !ARGUMENTS:
!EOP


! *** Model specifications ***
  dim_state   = (2*512)+7     ! State dimension
  forcing     = 0.0    ! Forcing parameter
  dt          = 0.05   ! Size of time step
  step_null   = 0      ! initial time step
  total_steps = 5000   ! number of time steps

  file_state = 'state.nc'  ! Name of output file for state trajectory
  delt_write = 1           ! Output interval in time steps (0 for no output)


! *** Parse command line options ***
  handle = 'dim_state'               ! state dimension of dummy model
  CALL parse(handle, dim_state)
  handle = 'forcing'                 ! Forcing parameter
  CALL parse(handle, forcing)
  handle = 'dt'                      ! Time step size
  CALL parse(handle, dt)
  handle = 'step_null'               ! Initial time step
  CALL parse(handle, step_null)
  handle = 'total_steps'             ! total number of time steps
  CALL parse(handle, total_steps)
  handle = 'file_state'              ! Name of output file
  CALL parse(handle, file_state)
  handle = 'delt_write'              ! Output interval
  CALL parse(handle, delt_write)

! *** Define final time step ***
  step_final = step_null + total_steps

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL'
     WRITE (*, '(22x,a)') 'MODEL: Mitchell-Schaeffer'
     WRITE (*, '(5x, a, i7)') &
          'Global model dimension:', dim_state
     WRITE (*, '(10x,a,es10.2)') 'Forcing parameter:', forcing
     WRITE (*, '(17x, a, i7, a, i7, a1)') &
          'Time steps:', total_steps, '  (final step:', step_final, ')'
     IF (step_null /= 0) WRITE (*, '(10x, a, i7)') 'Initial time step:', step_null
     WRITE (*, '(13x, a, f10.3)') 'Time step size:', dt
#ifndef USE_PDAF
     IF (delt_write > 0) THEN
        WRITE (*, '(5x,a,i4,a)') &
             'Write output after each ',delt_write,' time steps'
        WRITE (*, '(10x,a,a)') 'Output file: ',TRIM(file_state)
     ELSE
        WRITE (*, '(5x,a)') 'Do not write model output!'
     END IF
#endif
  END IF screen2


! *** Initialize dimensions and fields with domain decompsition

  ! Allocate a model field
  ALLOCATE(x(dim_state))
  ! count allocated memory
  CALL memcount(1, 'r', dim_state)

  ! Initialize model field to rest state
  x(1:((dim_state - 7)/2)) = 0.0
  x(((1+(dim_state - 7)/2)):(dim_state - 7)) = 1.0

#ifdef USE_PDAF
	! initialize the model parameter deviations (std.dev = 0.05)
	! print to stdout the model parameters
	CALL add_model_noise(forcing, 7, x((dim_state - 7 + 1):dim_state))
#else
  ! Initialize model parameters in state vector to zero
	x(dim_state - 7 + 1 : dim_state) = 0.0
#endif
  x(dim_state - 7 + 1) = x(dim_state - 7 + 1) + 0.1		! diffusive coupling
  x(dim_state - 7 + 2) = x(dim_state - 7 + 2) + 10.0		! sigmoidal strength
  x(dim_state - 7 + 3) = x(dim_state - 7 + 3) + 0.3		! ti
  x(dim_state - 7 + 4) = x(dim_state - 7 + 4) + 6.0		! tu
  x(dim_state - 7 + 5) = x(dim_state - 7 + 5) + 120.0		! to
  x(dim_state - 7 + 6) = x(dim_state - 7 + 6) + 150.0		! tc
  x(dim_state - 7 + 7) = x(dim_state - 7 + 7) + 0.13		! ug

  ! Initialize model time
  time = REAL(step_null) * dt

  ! Initialize netcdf output
#ifndef USE_PDAF
  CALL init_netcdf(step_null, time, dt, forcing, dim_state, x)
#endif

END SUBROUTINE initialize
