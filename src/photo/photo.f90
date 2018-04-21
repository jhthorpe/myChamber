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
  SUBROUTINE photo_read(ID_list,photo_rxns,photo_QY,photo_CS,photo_AF,photo_nrxn,nID)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: photo_QY, photo_CS, photo_AF 
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: photo_rxns
    INTEGER, INTENT(INOUT) :: photo_nrxn,nID
    !internal
    CHARACTER(LEN=20) :: fname
    CHARACTER(LEN=8) :: dummy
    INTEGER :: mID,mrxn,nlines,fID,start_ID,end_ID,idx
    INTEGER :: i,j,k

    WRITE(*,*)
    WRITE(*,*) "Reading photolysis data"

    start_ID = nID
    end_ID = -1

    !Initialize data
    ALLOCATE(photo_rxns(0:1,0:6))
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
    fID = 5
    OPEN(unit=fID,file=fname,status='old',access='sequential')
    READ(fID,*)
    READ(fID,*) 
    DO i=0,nlines-3
      READ(fID,*)
    END DO
    CLOSE(unit=fID)
  
    !check below is correct
    !end_ID = start_ID + nlines-3

    !DO i=start_ID,end_ID
    !END DO
    

  END SUBROUTINE photo_read
  

END MODULE photo
