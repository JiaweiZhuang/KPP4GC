PROGRAM gckpp_Check_RCONST

use netcdf
use GCKPP_Model
USE GCKPP_Initialize, ONLY: Initialize
USE gckpp_Precision, ONLY: sp, dp

IMPLICIT NONE

! =========================
! -- Variable declaration -
! =========================

! For IO
character (len = *), parameter :: INPUT_FILE_NAME = "./KPP_fields.nc"
character (len = *), parameter :: OUTPUT_FILE_NAME = "./RCONST.nc"
logical :: is_write = .true.

integer :: ncid, varid ! can reuse them for different variables
integer :: x_dimid, y_dimid, z_dimid, n_dimid, dimids(4) ! for NetCDF write

! For data
integer :: II, JJ, LL ! for looping through data
integer, parameter :: IIPAR=72, JJPAR=46, LLPAR=72
! NREACT are already defined by KPP
REAL*8 :: ALL_RCONST_COMPUTED(NREACT, IIPAR, JJPAR, LLPAR)
REAL*8 :: ALL_PHY(4, IIPAR, JJPAR, LLPAR)

! For timing
real :: timeStart, timeEnd

! =========================
! ---- Read input data ----
! =========================

WRITE (6,*) 'Read input data'

! Open the file with read-only access
call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )

! Read the data
call check( nf90_inq_varid(ncid, "PHY", varid) )
call check( nf90_get_var(ncid, varid, ALL_PHY) )

! Close the file, freeing all resources.
call check( nf90_close(ncid) )

CALL Initialize() ! don't seem to do anything

! --- timing starts here ---
call cpu_time(time=timeStart)

! Don't care about additional rates for now
PHOTOL(:) = 0.0_dp
HET(:,:) = 0.0_dp

DO LL=1,LLPAR
! LL = 1 ! save time, just do one layer

WRITE (6,*) 'Loop through layer', LL

DO JJ=1,JJPAR
DO II=1,IIPAR

! pass input data
TEMP = ALL_PHY(1,II,JJ,LL)
PRESS = ALL_PHY(2,II,JJ,LL)
NUMDEN = ALL_PHY(3,II,JJ,LL)
H2O = ALL_PHY(4,II,JJ,LL)

! Compute RCONST
CALL Update_RCONST()

! record output results
ALL_RCONST_COMPUTED(:,II,JJ,LL) = RCONST(:)

ENDDO ! FOR II
ENDDO ! FOR JJ

ENDDO ! FOR LL

call cpu_time(time=timeEnd)
! --- timing ends here ---
  
WRITE(6,*) "Time for RCONST computation", timeEnd-timeStart, "seconds"

! ========================= 
! ----  Write results  ----
! =========================

IF (is_write) THEN

WRITE (6,*) 'Write results'

call check( nf90_create(OUTPUT_FILE_NAME, NF90_CLOBBER, ncid) )

call check( nf90_def_dim(ncid, "lon", IIPAR, x_dimid) )
call check( nf90_def_dim(ncid, "lat", JJPAR, y_dimid) )
call check( nf90_def_dim(ncid, "lev", LLPAR, z_dimid) )
call check( nf90_def_dim(ncid, "nreact", NREACT, n_dimid) )

dimids =  (/n_dimid, x_dimid, y_dimid, z_dimid/) 
call check( nf90_def_var(ncid, "RCONST", NF90_DOUBLE, dimids, varid) )
call check( nf90_enddef(ncid) )
call check( nf90_put_var(ncid, varid, ALL_RCONST_COMPUTED) )

call check( nf90_close(ncid) )

ENDIF ! is_write

! ========================= 
! - NetCDF util function  -
! =========================

contains

  subroutine check(status)
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

END PROGRAM gckpp_Check_RCONST
