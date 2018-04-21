!Module for auxilary functions for gas phase kinetics reactions

MODULE gas
  USE util
  IMPLICIT NONE

  CONTAINS
!=====================================================================
!                       SUBPROGRAMS
!=====================================================================

!---------------------------------------------------------------------
! Read gas phase input 
!---------------------------------------------------------------------
! Variables
!	ID_list		:	1D char(4), ID's of chemical species
!	gas_rxns	:	2D int, reactions matrix
!	gas_coef	:	2D dp, coef of general rate eq
!	gas_nrxn	:	int, number of gas phase reactions
!	gas_nID		:	int, number of gas phase species

  SUBROUTINE gas_read(ID_list,gas_rxns,gas_coef,gas_nrxn,nID)
    IMPLICIT NONE

    !inout
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: gas_rxns
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: gas_coef
    INTEGER, INTENT(INOUT) :: gas_nrxn,nID

    !internal
    CHARACTER(LEN=20) :: fname
    CHARACTER(LEN=8) :: dummy
    INTEGER :: nlines,nrxn,fID,mID,mrxn,gas_nID
    INTEGER :: i,j,k,l

    WRITE(*,*) 
    WRITE(*,*) "Reading gas phase data..."

    !Initialize data
    ALLOCATE(gas_coef(0:1,0:7))
    ALLOCATE(gas_rxns(0:1,0:5))
    ALLOCATE(ID_list(0:3))
    gas_nrxn = 0
    dummy = 'X'
    ID_list(0) = dummy
    dummy = 'N2'
    ID_list(1) = dummy
    gas_nID = 2
    mrxn = 2
    mID = 4

    !Read gas phase data
    fname = 'gas.txt'
    nlines = file_numlines(fname)
    fID = 2
    OPEN(unit=fID,file=fname,status='old',access='sequential')
    READ(fID,*)
    READ(fID,*)
    DO i=0,nlines-3
      gas_nrxn = gas_nrxn + 1
      IF (gas_nrxn .GT. mrxn) THEN
        CALL int4_2Dgrow1(gas_rxns)
        CALL real8_2Dgrow1(gas_coef)
        mrxn = mrxn*2
      END IF 
      CALL gas_readline(fID,ID_list,gas_rxns(i,:),gas_coef(i,:),gas_nID,mID)   
    END DO
    CLOSE(unit=fID)
    nID = nID + gas_nID

    !write to ID file
    OPEN(unit=3,file='ID.txt',status='replace',access='sequential')
    DO i=0,nID-1
      WRITE(3,*) i,ID_list(i)
    END DO
    CLOSE(unit=3,status='keep')

    WRITE(*,*) "ID's stored in ID.txt"

  END SUBROUTINE gas_read
    
!---------------------------------------------------------------------
! Read gas.txt data and add to datasets
!---------------------------------------------------------------------
! WARNING - assuming array index from 0
! Variables
!	ID_list		:	1D char(4), ID's of chemical species
!	gas_rxns	:	12D int, reactions matrix
!	gas_coef	:	1D dp, coef of general rate eq
!	gas_nID		:	int, number of gas phase species 
!	mID		:	int, max number of species

  SUBROUTINE gas_readline(fID,ID_list,gas_rxns,gas_coef,gas_nID,mID)
    IMPLICIT NONE

    !inout
    INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: gas_rxns
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: gas_coef
    INTEGER, INTENT(INOUT) :: gas_nID,mID
    INTEGER, INTENT(IN) :: fID
    
    !internal
    CHARACTER(LEN=8),DIMENSION(0:6) :: ID
    CHARACTER(LEN=2) :: dum2
    CHARACTER(LEN=1) :: dum1
    REAL(KIND=8) :: a,b,c,d,e,f,g,h
    INTEGER :: idx,i

    !read in variables
    READ(fID,*) ID(0),dum1,ID(1),dum1,ID(2),dum1,ID(3),dum1,ID(4),dum1,ID(5), &
                dum1,dum2,a,dum2,b,dum2,c,dum2,d,dum2,e,dum2,f,dum2,g,dum2,ID(6)

    DO i=0,5 

    !find ID and index of species 
      idx = ID2idx(ID(i),ID_list)

      !add to ID_list
      IF (idx .LT. 0) THEN
        !Check if we're about to overflow
        IF (gas_nID + 1 .GT. mID) THEN
          CALL chr8_1Dgrow(ID_list)
          mID = mID*2
        END IF
        ID_list(gas_nID) = ID(i)
        gas_rxns(i) = gas_nID
        gas_nID = gas_nID + 1
      ELSE
        gas_rxns(i) = idx
      END IF

    END DO

    !get ID number of h, and add coefs to gas_coef matrix
    idx = ID2idx(ID(6),ID_list)
    IF (idx .LT. 0) THEN
      IF (gas_nID + 1 .GT. mID) THEN
        CALL chr8_1Dgrow(ID_list)
        mID = mID*2
      END IF
      ID_list(gas_nID) = ID(4)
      gas_coef = [a,b,c,d,e,f,g,gas_nID*1.0D0]
      gas_nID = gas_nID + 1
    ELSE
      gas_coef = [a,b,c,d,e,f,g,idx*1.0D0]
    END IF

  END SUBROUTINE gas_readline

!---------------------------------------------------------------------
! Calculates rate of a reaction from general rate eq 
!---------------------------------------------------------------------
! WARNING - assume index from zero
! Variables
!	A		:	1D real(8), coefficients of the reaction 
!	conc		:	1D real(8), concentrations of species	
!	T		:	real(8), temperature in K

  REAL(KIND=8) FUNCTION gas_k(A,T,conc)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: A,conc
    REAL(KIND=8), INTENT(IN) :: T
    INTEGER :: idx
    idx = NINT(A(7))
    IF (idx .EQ. 0) THEN
      gas_k = A(0)*DEXP(A(1)*T/A(2))*DEXP(A(3)/T)*(A(4)*T/A(5))**A(6)
    ELSE
      gas_k = A(0)*DEXP(A(1)*T/A(2))*DEXP(A(3)/T)*(A(4)*T/A(5))**A(6)*conc(idx)
    END IF 
  END FUNCTION gas_k

!---------------------------------------------------------------------
END MODULE gas
