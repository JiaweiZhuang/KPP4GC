SUBROUTINE onestep(all_c_before, all_phy, all_photol, n_sample, all_c_after)
! Start from many initial conditions and integrate a single time step

use GCKPP_Model
USE GCKPP_Initialize, ONLY: Initialize
USE gckpp_Precision, ONLY: sp, dp

IMPLICIT NONE

! =========================
! -- Variable declaration -
! =========================

! For KPP solver
REAL(dp) :: TIN, TOUT ! DT is already defined in gckpp_Global
INTEGER  :: ICNTRL(20), ISTATUS(20), IERR
REAL(dp) :: RCNTRL(20), RSTATE(20)

! For internal use
INTEGER, PARAMETER :: JVN_ = 130 ! number of JVals; defined in Headers/CMN_FJX_MOD.F
INTEGER, PARAMETER :: NSPEC_ = 240 ! nspec already defined in gckpp_Parameters.f90
integer :: II ! for looping through data

! For interfacing with f2py
! use lower case to be consistent with Python signature
INTEGER :: n_sample ! total number of samples to loop over
REAL*8 :: all_c_before(NSPEC_, n_sample) ! NSPEC is already defined by KPP,
REAL*8 :: all_c_after(NSPEC_, n_sample) ! but cannot be used in f2py declaration
REAL*8 :: all_phy(4, n_sample)
REAL*8 :: all_photol(JVN_, n_sample)

! below are important F2PY signatures, not comments!

!f2py intent(in) :: n_sample, all_c_before, all_phy, all_photol
!f2py intent(out) :: all_c_after
!f2py depend(n_sample) :: all_c_after, all_phy, all_photol

! =========================
! --- KPP configuration ---
! =========================

!%%%%% TIMESTEPS %%%%%
DT        = 30d0 * 60d0 ! in seconds
TIN       = 0d0
TOUT      = TIN + DT
 
!%%%%% CONVERGENCE CRITERIA () %%%%%
! ATOL(NVAR), RTOL(NVAR) are defined in gckpp_Global.F90
ATOL      = 1e-2_dp    
RTOL      = 0.5e-2_dp 

!%%%%% SOLVER OPTIONS %%%%%
! taken from flexchem_mod.F90
ICNTRL    = 0 ! Zero all slots of IC
ICNTRL(1) = 1
ICNTRL(2) = 0
ICNTRL(3) = 4
ICNTRL(7) = 1
RCNTRL    = 0.0_dp ! Zero all slots of RCNTRL
!! RCNTRL(3) = HSAVE_KPP(I,J,L) ! Only do one time step. Don't have this.

! ========================= 
! Pass data to KPP and solve
! =========================

! Don't care about heterogenous reactions for now
HET(:,:) = 0.0_dp

! Loop over different initial conditions
DO II = 1, n_sample

  ! -- pass input data --

  ! initial concentrations
  VAR(1:NVAR) = all_c_before(1:NVAR,II)
  FIX(:) = all_c_before(NVAR+1:NSPEC_,II)

  ! to compute themo reaction rate
  TEMP = all_phy(1,II)
  PRESS = all_phy(2,II)
  NUMDEN = all_phy(3,II)
  H2O = all_phy(4,II)

  ! photolysis reaction rate
  PHOTOL(1:JVN_) = all_photol(:,II)

  CALL Update_RCONST()


  ! -- solve --
  ! Call the KPP integrator
  CALL Integrate( TIN,    TOUT,    ICNTRL,      &
                 RCNTRL, ISTATUS, RSTATE, IERR )

  ! Print grid box indices to screen if integrate failed
  IF ( IERR < 0 ) THEN
     WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', II
  ENDIF
  
  ! -- save result --

  all_c_after(1:NVAR,II) = VAR(1:NVAR)
  all_c_after(NVAR+1:NSPEC_,II) = FIX(:)

END DO ! End loop over samples

END SUBROUTINE