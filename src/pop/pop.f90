!Module that contains population based functions

MODULE pop
  USE util
  IMPLICIT NONE

  CONTAINS

!=====================================================================
!                      SUBPROGRAMS 
!=====================================================================

!---------------------------------------------------------------------
! Initialize species concentrations
!---------------------------------------------------------------------
! Variables
!	conc	:	1D real(8), list of populations
!	ID_list	:	1D char(8), list of species	

  SUBROUTINE pop_init(conc,ID_list)
    IMPLICIT NONE
    REAl(KIND=8), PARAMETER :: ppb2mpc=2.46D10
    !inout
    CHARACTER(LEN=8),DIMENSION(0:), INTENT(IN) :: ID_list
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: conc
    !internal
    CHARACTER(LEN=20) :: fname
    CHARACTER(LEN=8) :: ID
    REAL(KIND=8) :: ppb
    INTEGER :: i,nline,idx
    LOGICAL :: ex

    fname='pop0.txt'

    !read in ppb 
    INQUIRE(file=fname,EXIST=ex)
    IF (.NOT. ex) THEN
      WRITE(*,*) "@pop:pop_init - failed to open pop0.txt"
      STOP
    END IF

    conc = (/ (0.0D0, i=0, SIZE(conc)-1) /)
    
    !assign to ID
    nline = file_numlines(fname) 
    OPEN(unit=4,file=fname,status='old',access='sequential')
    READ(4,*) 
    READ(4,*)
    DO i=0,nline-3
      READ(4,*) ID,ppb 
      idx = ID2idx(ID,ID_list) 
      IF (idx .GE. 0) conc(idx) = ppb*ppb2mpc
    END DO
    CLOSE(unit=4)

  END SUBROUTINE pop_init

!---------------------------------------------------------------------
END MODULE pop
