! Copyright (c) 2004-2024 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$
!BOP
!
! !ROUTINE: PDAFomi_put_state_nonlin_nondiagR --- Interface to PDAF for global filters
!
! !INTERFACE:
SUBROUTINE PDAFomi_put_state_nonlin_nondiagR(collect_state_pdaf, &
     init_dim_obs_pdafomi, obs_op_pdafomi, likelihood_pdafomi, prepoststep_pdaf, &
     outflag)

! !DESCRIPTION:
! Interface routine called from the model during the 
! forecast of each ensemble state to transfer data
! from the model to PDAF and to perform the analysis
! step.
!
! This routine provides the simplified interface
! where names of user-provided subroutines are
! fixed. It simply calls the routine with the
! full interface using pre-defined routine names.
!
! The routine supports all global filters.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2024-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                ! User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdafomi, &  ! Initialize dimension of observation vector
       obs_op_pdafomi, &               ! Observation operator
       likelihood_pdafomi              ! Compute likelihood
  EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize observation vector
       PDAFomi_init_obsvar_cb, &       ! Initialize mean observation error variance
       PDAFomi_init_obscovar_cb, &     ! Initialize mean observation error variance
       PDAFomi_add_obs_error_cb        ! Add observation error covariance matrix

! !CALLING SEQUENCE:
! Called by: model code  
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_nonlin_nondiagR -- START'

  IF (TRIM(filterstr) == 'NETF') THEN
     CALL PDAF_put_state_netf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, likelihood_pdafomi, outflag)
  ELSEIF (TRIM(filterstr) == 'PF') THEN
     CALL PDAF_put_state_pf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, likelihood_pdafomi, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_nonlin_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_nonlin_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_nonlin_nondiagR
