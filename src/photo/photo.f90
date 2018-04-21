! Module for photolysis reactions
MODULE photo
  USE util
  IMPLICIT NONE

  CONTAINS

!=====================================================================
!                       SUBPROGRAMS
!=====================================================================

!---------------------------------------------------------------------
! Read photolysis input 
!---------------------------------------------------------------------
! Variables
!       ID_list         :       1D char(8), ID's of chemical species
!	photo_rxns	:	2D int, reactions matrix
!	photo_QY	:	2D real(8), array of quantum yields
!	photo_CS	:	2D real(8), array of cross section
!	photo_AF	:	2D real(8), array of actinic flux	
!	photo_nrxn	:	int, number of photolysis reactions
!	nID		:	int, number of unique species
! NOTES
!    Here, I will use 1 nm as the integer values between 0 and 1000 nm wavelengths
!    There are more efficient ways to do this, but I am lazy 
  SUBROUTINE photo_read(ID_list,photo_rxns,photo_QY,photo_CS,photo_AF,photo_nrxn,nID,photo_ID)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: photo_QY, photo_CS, photo_AF 
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: photo_rxns
    INTEGER, INTENT(INOUT) :: photo_nrxn,nID,photo_ID
    !internal
    CHARACTER(LEN=20) :: fname
    CHARACTER(LEN=8) :: dummy
    INTEGER :: mID,mrxn,nlines,fID,idx
    INTEGER :: i,j,k

    WRITE(*,*)
    WRITE(*,*) "Reading photolysis data..."

    photo_ID = nID

    !Initialize data
    ALLOCATE(photo_rxns(0:1,0:5))
    ALLOCATE(photo_QY(0:1,0:1000))
    ALLOCATE(photo_CS(0:1,0:1000))
    ALLOCATE(photo_AF(0:1,0:1000))
    mID = SIZE(ID_list)
    mrxn = 2
    photo_nrxn = 0
    
    !add in X and hv if they don't already exist
    IF (nID + 2 .GT. mID) THEN
      CALL chr8_1Dgrow(ID_list)
      mID = mID*2
    END IF
    dummy = 'X'
    idx = ID2idx(dummy,ID_list)
    IF (idx .LT. 0) THEN
      ID_list(nID) = dummy
      nID = nID + 1
    END IF
    dummy = 'hv'
    idx = ID2idx(dummy,ID_list)
    IF (idx .LT. 0) THEN
      ID_list(nID) = dummy
      nID = nID + 1 
    END IF

    !read photo.txt data to get reaction list 
    fname = 'photo.txt'
    nlines = file_numlines(fname)
    fID = 1
    OPEN(unit=fID,file=fname,status='old',access='sequential')
    READ(fID,*)
    READ(fID,*) 
    DO i=0,nlines-3
      photo_nrxn = photo_nrxn + 1 
      IF (photo_nrxn .GT. mrxn) THEN
        CALL real8_2Dgrow1(photo_QY)
        CALL real8_2Dgrow1(photo_CS)
        CALL real8_2Dgrow1(photo_AF)
        CALL int4_2Dgrow1(photo_rxns)
        mrxn = mrxn*2
      END IF
      CALL photo_readline(fID,ID_list,photo_rxns(i,:),nID,mID)
    END DO
    CLOSE(unit=fID)

    !append data to ID.txt
    OPEN(unit=1,file='ID.txt',access='APPEND',status='old')
    DO i=photo_ID,nID-1
      WRITE(1,*) i,ID_list(i)
    END DO
    CLOSE(unit=1)
    WRITE(*,*) "photolysis ID appended to ID.txt"

    !get QY,CS,and AF data
    CALL photo_getQY(ID_list,photo_rxns,photo_nrxn,photo_QY)

  END SUBROUTINE photo_read

!---------------------------------------------------------------------
! Read photo.txt data and add to rxn matrix
!---------------------------------------------------------------------
! WARNING - assuming array index from 0
! Variables
!	fID		:	int, readfile ID
!	ID_list		:	1D char(8), ID list
!	photo_rxns	:	1D int, line of reactions
!	nID		:	int, number of unique species
!	mID		:	int, maxID of ID_list

  SUBROUTINE photo_readline(fID,ID_list,photo_rxns,nID,mID)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    INTEGER, DIMENSION(0:), INTENT(INOUT) :: photo_rxns
    INTEGER, INTENT(INOUT) :: nID, mID
    INTEGER, INTENT(IN) :: fID
    !internal
    CHARACTER(LEN=8),DIMENSION(0:5) :: ID
    CHARACTER(LEN=1) :: dum
    INTEGER :: idx
    INTEGER :: i,j,k
    !read in variables
    READ(fID,*) ID(0),dum,ID(1),dum,ID(2),dum,ID(3),dum,ID(4),dum,ID(5)
    !find ID of species
    DO i=0,5
      idx = ID2idx(ID(i),ID_list) 
      IF (idx .LT. 0) THEN
        !check for overflow
        IF (nID+1 .GT. mID) THEN
          CALL chr8_1Dgrow(ID_list)
          mID = mID*2
        END IF
        ID_list(nID) = ID(i)
        photo_rxns(i) = nID
        nID = nID + 1
      ELSE
        photo_rxns(i) = idx
      END IF
    END DO
  END SUBROUTINE photo_readline

!---------------------------------------------------------------------
! Read in Quantum Yield Data
!---------------------------------------------------------------------
! Values
!	ID_list		:	1D char(8), list of ID's
!	photo_rxns	:	2D int, list of photolysis reactions
!	photo_nrxn	:	int, number of photolysis reactions 
!	photo_QY	:	2D real(8), Quantum Yield matrix

  SUBROUTINE photo_getQY(ID_list,photo_rxns,photo_nrxn,photo_QY)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:),INTENT(INOUT) :: photo_QY
    CHARACTER(LEN=8), DIMENSION(0:), INTENT(IN) :: ID_list
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: photo_rxns
    INTEGER, INTENT(IN) :: photo_nrxn
    !internal
    CHARACTER(LEN=20) :: fname
    INTEGER :: fID,val
    INTEGER :: i,j,k

    WRITE(*,*) "Reading data from QY.txt"
  
    fname = 'QY.txt'
    fID = 1
    
    !go through each reaction
    OPEN(unit=fID,file=fname,status='old',access='sequential',action='READ')
    DO i=0,photo_nrxn-1
      val = photo_goto(fID,ID_list,photo_rxns(i,:))
      IF (val .LT. 0) THEN
        WRITE(*,*) "failed to find this reaction in QY.txt"
        WRITE(*,*) photo_rxns(i,:)
        CLOSE(unit=fID)
        STOP
      END IF
      
      REWIND(unit=fID)   
    END DO 
    CLOSE(unit=fID)
    
  END SUBROUTINE photo_getQY

!---------------------------------------------------------------------
! Goto start of reaction in file
!---------------------------------------------------------------------
! Values
!	fID		:	int, unit id of file
!	ID_list		:	1D char(8), list of species ID
!	photo_rxns	:	1D int, list of species in reaction
!	rID		:	1D char(8), list of read IDs
!	pID		:	1D char(8), list of photo ID

  INTEGER FUNCTION photo_goto(fID,ID_list,photo_rxns)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8),DIMENSION(0:),INTENT(IN) :: ID_list
    INTEGER, DIMENSION(0:), INTENT(IN) :: photo_rxns
    INTEGER, INTENT(IN) :: fID
    !internal
    CHARACTER(LEN=8), DIMENSION(0:5) :: rID,pID
    CHARACTER(LEN=8) :: str
    INTEGER :: val,io,eq
    INTEGER :: i,j,k
    LOGICAL :: found,equal

    val = -1
    found = .FALSE.
  
    !get the characters we are looking for
    str = 'notamol'
    DO i=0,SIZE(photo_rxns)-1
      pID(i) = ID_list(photo_rxns(i))
    END DO 

    !goto line 
    DO WHILE(.NOT. found)
      READ(fID,*,iostat=io) str
      IF (io .NE. 0) STOP
  
      !we have a potential match, with correct format
      IF (str .EQ. pID(0)) THEN
        BACKSPACE(unit=fID)
        READ(fID,*) rID

        !determine if they are equal
        eq = 1
        DO j=0,SIZE(pID)-1
          IF (pID(j) .NE. rID(j)) eq = eq * 0
        END DO
        IF (eq .EQ. 1) THEN
          val = 0
          found = .TRUE.
        END IF
     
      END IF

    END DO
    
    photo_goto = val

  END FUNCTION photo_goto

!---------------------------------------------------------------------

END MODULE photo
