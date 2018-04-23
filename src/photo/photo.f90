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
!
!    Note that actinic flux is currently considered a constant, this should be fixed

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
    ALLOCATE(photo_AF(0:0,0:1000))
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
        CALL int4_2Dgrow1(photo_rxns)
        mrxn = mrxn*2
      END IF
      CALL photo_readline(fID,ID_list,photo_rxns(i,:),nID,mID)
    END DO
    CLOSE(unit=fID)

    !zero arrays
    CALL real8_2Dzero(photo_QY)
    CALL real8_2Dzero(photo_CS)
    CALL real8_2Dzero(photo_AF)

    !append data to ID.txt
    OPEN(unit=1,file='ID.txt',access='APPEND',status='old')
    DO i=photo_ID,nID-1
      WRITE(1,*) i,ID_list(i)
    END DO
    CLOSE(unit=1)
    WRITE(*,*) "photolysis ID appended to ID.txt"

    !get QY,CS for each reaction
    fname = 'QY.txt'
    CALL photo_getdata(ID_list,photo_rxns,photo_nrxn,photo_QY,fname)
    fname = 'CS.txt'
    CALL photo_getdata(ID_list,photo_rxns,photo_nrxn,photo_CS,fname)
    !get AF for the system
    fname = 'AF.txt'
    CALL photo_getAF(photo_AF,fname)

    WRITE(*,*) 

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
! Read in Data from designated file 
!---------------------------------------------------------------------
! Values
!	ID_list		:	1D char(8), list of ID's
!	photo_rxns	:	2D int, list of photolysis reactions
!	photo_nrxn	:	int, number of photolysis reactions 
!	photo_QY	:	2D real(8), Quantum Yield matrix

  SUBROUTINE photo_getdata(ID_list,photo_rxns,photo_nrxn,photo_data,fname)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:),INTENT(INOUT) :: photo_data
    CHARACTER(LEN=8), DIMENSION(0:), INTENT(IN) :: ID_list
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: photo_rxns
    CHARACTER(LEN=20),INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: photo_nrxn
    !internal
    INTEGER :: fID,val,io
    INTEGER :: i,j,k

    WRITE(*,*) "Reading data from ", fname
    fID = 1
    
    !go through each reaction
    OPEN(unit=fID,file=fname,status='old',access='sequential',action='READ')
    DO i=0,photo_nrxn-1

      !goto start of reaction
      val = photo_goto(fID,ID_list,photo_rxns(i,:))
      IF (val .LT. 0) THEN
        WRITE(*,*) "failed to find this reaction in ",fname
        WRITE(*,*) photo_rxns(i,:)
        CLOSE(unit=fID)
        STOP
      END IF

      !get data
      CALL photo_filldata(fID,photo_data(i,:))

      REWIND(unit=fID)   
    END DO 
    CLOSE(unit=fID)
    
  END SUBROUTINE photo_getdata

!---------------------------------------------------------------------
! Goto start of reaction in given file
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

    !search for first match
    DO WHILE(.NOT. found)
      READ(fID,*,iostat=io) str
      IF (io .NE. 0) RETURN
  
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
!  Fill in photo data from designated file until hit io fail
!---------------------------------------------------------------------
! Values
!	fID		:	int, unit number of file
!	photo_data	:	1D real(8), data to fill

  SUBROUTINE photo_filldata(fID,photo_data)
    IMPLICIT NONE
    !inout
    REAL(KIND=8),DIMENSION(0:),INTENT(INOUT) :: photo_data
    INTEGER, INTENT(IN) :: fID
    !internal
    REAL(KIND=8) :: peak,width,val,fact
    INTEGER :: io
    INTEGER :: i,j,k
    io = 0
    !read factor
    READ(fID,*,iostat=io) fact
    IF (io .NE. 0) STOP "need factor in photolysis data"
    DO WHILE (io .EQ. 0) 
      READ(fID,*,iostat=io) peak,width,val
      !if good, fill in data
      IF (io .EQ. 0) THEN
        DO i=NINT(peak)-NINT(width),NINT(peak)+NINT(width)
          photo_data(i) = val/fact
        END DO 
      ELSE
        RETURN
      END IF 
    END DO
    BACKSPACE(unit=fID)
  END SUBROUTINE photo_filldata

!---------------------------------------------------------------------
! Get the Actinic Flux for our system
!---------------------------------------------------------------------
! WARNING - currently considereing only constant AF
! Values
!	photo_AF	:	2D real(8), actinic flux per timestep
!	fname		:	char(20), filename

  SUBROUTINE photo_getAF(photo_AF,fname)
    IMPLICIT NONE
    !inout
    REAL(KIND=8),DIMENSION(0:,0:), INTENT(INOUT) :: photo_AF
    CHARACTER(LEN=20), INTENT(IN) :: fname
    !internal
    INTEGER :: fID
    INTEGER :: i
    fID = 1
    WRITE(*,*) "Reading data from ",fname
    OPEN(unit=fID,file=fname,access='sequential',status='old')
    READ(fID,*) 
    READ(fID,*) 
    READ(fID,*) 
    READ(fID,*) 
    CALL photo_filldata(fID,photo_AF(0,:)) 
    CLOSE(unit=fID)
  END SUBROUTINE photo_getAF

!---------------------------------------------------------------------
! Perform one step of photolysis reactions 
!---------------------------------------------------------------------
! WARNING - currently considereing only constant AF

  SUBROUTINE photo_react(tstep,ID_list,conc,nconc,photo_rxns,photo_QY,&
                         photo_CS,photo_AF)
    IMPLICIT NONE
    !inout
    CHARACTER(LEN=8),DIMENSION(0:), INTENT(IN) :: ID_list
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: conc,nconc
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: photo_QY,photo_CS,photo_AF
    INTEGER, DIMENSION(0:), INTENT(IN) :: photo_rxns
    REAL(KIND=8), INTENT(IN) :: tstep
    !internal
    REAL(KIND=8) :: rate,coef
    INTEGER :: i,j,k

    !k
    coef = photo_k(photo_QY,photo_CS,photo_AF)
    rate = coef

    !rate
    DO i=0,2
      IF (photo_rxns(i) .GT. 3 .OR. photo_rxns(i) .EQ. 1) rate = rate*conc(photo_rxns(i))
    END DO

     !reactants
     DO i=0,2
       IF (photo_rxns(i) .GT. 3 .OR. photo_rxns(i) .EQ. 1) THEN
         nconc(photo_rxns(i)) = conc(photo_rxns(i)) - rate*tstep
       END IF
     END DO
     !products
     DO i=3,5
       IF (photo_rxns(i) .GT. 3 .OR. photo_rxns(i) .EQ. 1) THEN
         nconc(photo_rxns(i)) = conc(photo_rxns(i)) + rate*tstep
       END IF
     END DO
 
     !check for zeros
     DO i=0,5
       IF (nconc(photo_rxns(i)) .LT. 0.0D0) THEN
!         WRITE(*,*) "WARNING - ", ID_list(photo_rxns(i)), " has a concentration of", nconc(photo_rxns(i)), ": zeroing"
!         nconc(photo_rxns(i)) = 0.0D0
       END IF
     END DO

  END SUBROUTINE photo_react
 
!---------------------------------------------------------------------
! Get photolysis rate constant 
!---------------------------------------------------------------------
! WARNING - currently considereing only constant AF
! WARNING - currently assuming QY,CS,and AF are equally indexed
  REAL(KIND=8) FUNCTION photo_k(QY,CS,AF)
    IMPLICIT NONE
    !inout
    REAL(KIND=8),DIMENSION(0:),INTENT(IN) :: QY,CS,AF
    !internal
    REAL(KIND=8) :: val
    INTEGER :: i,j,k
    val = 0.0D0
    DO i=0,SIZE(QY)-1 
      val = val + QY(i)*CS(i)*AF(i)
    END DO
    photo_k = val
  END FUNCTION photo_k
  
!---------------------------------------------------------------------
END MODULE photo
