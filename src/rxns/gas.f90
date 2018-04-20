!Module for auxilary functions for gas phase kinetics reactions

MODULE gas
  USE util
  IMPLICIT NONE

  CONTAINS
!=====================================================================
!                       SUBROUTINES 
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

  SUBROUTINE gas_read(ID_list,gas_rxns,gas_coef,gas_nrxn,gas_nID)
    IMPLICIT NONE

    !inout
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: gas_rxns
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: gas_coef
    INTEGER, INTENT(INOUT) :: gas_nrxn,gas_nID

    !internal
    CHARACTER(LEN=20) :: fname
    INTEGER :: nlines,nrxn,nID,fID,mID,mrxn
    INTEGER :: i,j,k,l

    !Initialize data
    ALLOCATE(gas_coef(0:1,0:7))
    ALLOCATE(gas_rxns(0:1,0:3))
    ALLOCATE(ID_list(0:3))
    gas_nrxn = 0
    gas_nID = 0
    mrxn = 2
    mID = 4

    !Read gas phase data
    fname = 'gas.txt'
    nlines = file_numlines(fname)
    fID = 2
    OPEN(unit=fID,file=fname,status='old',access='sequential')
    READ(fID,*)
    READ(fID,*)
    DO i=0,nlines-4
      gas_nrxn = gas_nrxn + 1
      IF (gas_nrxn .GT. mrxn) THEN
        CALL int4_2Dgrow1(gas_rxns)
        CALL real8_2Dgrow1(gas_coef)
        mrxn = mrxn*2
      END IF 
      CALL gas_readline(fID,ID_list,gas_rxns(i,:),gas_coef(i,:),gas_nID,mID)   
    END DO
    CLOSE(unit=fID)

    WRITE(*,*) ID_list
    WRITE(*,*) 
    DO i=0,gas_nID-1
      WRITE(*,*) i,ID_list(i)
    END DO

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
    CHARACTER(LEN=8),DIMENSION(0:4) :: ID
    CHARACTER(LEN=2) :: dum2
    CHARACTER(LEN=1) :: dum1
    REAL(KIND=8) :: a,b,c,d,e,f,g,h
    INTEGER :: idx,i

    !read in variables
    READ(fID,*) ID(0),dum1,ID(1),dum1,ID(2),dum1,ID(3), &
                dum1,dum2,a,dum2,b,dum2,c,dum2,d,dum2,e,dum2,f,dum2,g,dum2,ID(4)

    DO i=0,3 

    !find ID and index of species 
      idx = gas_ID2idx(ID(i),ID_list)

      !add to ID_list
      IF (idx .LT. 0) THEN
        !Check if we're about to overflow
        IF (gas_nID + 1 .GT. mID) THEN
          CALL chr8_1Dgrow(ID_list)
          mID = mID*2
        END IF
        !add to list
        WRITE(*,*) ID(i),"not found in ID_list, adding"
        ID_list(gas_nID) = ID(i)
        gas_nID = gas_nID + 1
      ELSE
        WRITE(*,*) ID(i),"found in list, skipping"
      END IF

      !Add species to gas_rxns matrix
      gas_rxns(i) = idx

    END DO

    !get ID number of h, and add coefs to gas_coef matrix
    idx = gas_ID2idx(ID(4),ID_list)
    IF (idx .LT. 0) THEN
      IF (gas_nID + 1 .GT. mID) THEN
        CALL chr8_1Dgrow(ID_list)
        mID = mID*2
      END IF
      WRITE(*,*) ID(4),"not found in list, adding"
      ID_list(gas_nID) = ID(4)
      gas_nID = gas_nID + 1
    ELSE 
      WRITE(*,*) ID(4), "found in list, skipping"
    END IF
    gas_coef = [a,b,c,d,e,f,g,idx*1.0D0]

  END SUBROUTINE gas_readline

!---------------------------------------------------------------------
! finds index of ID from list, -1 if not there
!---------------------------------------------------------------------
! WARNING - assume index from zero
! Variables
!	ID	:	char(8), ID of chemical
!	ID_list	:	1D int, list of chemical ID's

  INTEGER FUNCTION gas_ID2idx(ID,ID_list)
    IMPLICIT NONE
    CHARACTER(LEN=8), DIMENSION(0:), INTENT(IN) :: ID_list
    CHARACTER(LEN=8), INTENT(IN) :: ID
    INTEGER :: i,idx
    idx = -1
    DO i=0,SIZE(ID_list)-1
      IF (ID_list(i) .EQ. ID) THEN
        idx = i
        EXIT
      END IF
    END DO
    gas_ID2idx = idx
  END FUNCTION gas_ID2idx

!---------------------------------------------------------------------
END MODULE gas
