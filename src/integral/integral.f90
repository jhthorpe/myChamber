!Module for integration of populations

MODULE integral
  USE gas
  USE photo
  IMPLICIT NONE

  CONTAINS

!=====================================================================
!                       SUBPROGRAMS
!=====================================================================

!---------------------------------------------------------------------
! Integration control scheme 
!---------------------------------------------------------------------
! Variables
!	
! NOTES
!    This is hardcoded to deal with constant Actinic Flux, and should be fixed later
 
  SUBROUTINE integrate(ID_list,nID,gas_rxns,gas_nrxn,gas_coef,photo_rxns,photo_nrxn,&
                       photo_QY,photo_CS,photo_AF,photo_ID,conc,options)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8),DIMENSION(0:), INTENT(IN) :: ID_list
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: gas_coef,photo_QY,photo_CS,photo_AF
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: conc
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: options
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: gas_rxns,photo_rxns
    INTEGER, INTENT(IN) :: nID,photo_ID,gas_nrxn,photo_nrxn
    !internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: nconc
    CHARACTER(LEN=20) :: fname
    CHARACTER(LEN=8) :: d2
    REAL(KIND=8) :: time, tstart, tstep, tlen, T 
    INTEGER :: fID
    INTEGER :: i,j,k

    WRITE(*,*) "Starting integration"
    WRITE(*,*) "populations in ppb"

    fname = 'pop.out' 
    fID = 42 
    tstart = options(4)
    tlen = options(5)
    tstep = options(6)
    T = options(7)
    time = tstart
    ALLOCATE(nconc(0:nID-1))
    nconc = conc

    OPEN(unit=fID,file=fname,status='replace',access='sequential')
    d2='time'
    WRITE(fID,*) d2,ID_list(4:nID-1) 
    DO WHILE (time .LT. tstart+tlen)
      WRITE(fID,*) time, conc(4:)/2.46D10
      CALL integral_step(tstep,ID_list,nID,photo_ID,conc,nconc,gas_rxns,gas_nrxn,gas_coef,&
                         photo_rxns,photo_nrxn,photo_QY,photo_CS,photo_AF,T,options)
      time = time + tstep
    END DO

    CLOSE(unit=fID,status='keep')
    DEALLOCATE(nconc)

  END SUBROUTINE integrate

!---------------------------------------------------------------------
! One timestep of integral 
!---------------------------------------------------------------------
! Variables

  SUBROUTINE integral_step(tstep,ID_list,nID,photo_ID,conc,nconc,gas_rxns,gas_nrxn,gas_coef,photo_rxns,&
                           photo_nrxn,photo_QY,photo_CS,photo_AF,T,options)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8),DIMENSION(0:), INTENT(IN) :: ID_list
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: gas_coef,photo_QY,photo_CS,photo_AF
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: conc,nconc
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: options
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: gas_rxns,photo_rxns
    REAL(KIND=8), INTENT(IN) :: tstep,T 
    INTEGER, INTENT(IN) :: nID,photo_ID,gas_nrxn,photo_nrxn
    !internal
    INTEGER :: i,j,k     

    !go through gas phase reactions
    
    DO i=0,gas_nrxn-1
      CALL gas_react(tstep,ID_list,conc,nconc,gas_rxns(i,:),gas_coef(i,:),T)
    END DO

    !go through photolysis reactions
    DO i=0,photo_nrxn-1
      CALL photo_react(tstep,ID_list,conc,nconc,photo_rxns(i,:),photo_QY(i,:),&
                       photo_CS(i,:),photo_AF(i,:))
    END DO

    conc = nconc

  END SUBROUTINE integral_step

!---------------------------------------------------------------------

END MODULE integral
