! Module that contains useful utility functions

MODULE util
  IMPLICIT NONE

  CONTAINS

!---------------------------------------------------------------------
! Double size of 1st dimension of 2D-int(kind=4) array 
!---------------------------------------------------------------------
! WARNING - currently assumes that array indexs from 0
! Variables
!	A	:	2D int(4), array to be reallocated

  SUBROUTINE int4_2Dgrow1(A)
    IMPLICIT NONE
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: A
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: B
    INTEGER :: stat,n,m,i,j

    !get shape of A
    n = SIZE(A(:,0))
    m = SIZE(A(0,:))
    
    !check if allocatable
    IF (.NOT. ALLOCATED(A)) THEN
      WRITE(*,*) "@util:int4_2Dreallocate - array not allocated"
      STOP
    ELSE
      ALLOCATE(B(0:n-1,0:m-1),STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not deallocate array"
        STOP
      END IF

      !temp copy
      B(:,:) = A(:,:)

      DEALLOCATE(A,STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not deallocate array"
        STOP
      END IF

      !resize
      ALLOCATE(A(0:2*n-1,0:m-1),STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not allocate array"
        STOP
      END IF
       
      !move the data
      DO j=0,m-1
        DO i=0,n-1
          A(i,j) = B(i,j)
        END DO 
      END DO
      DEALLOCATE(B)

    END IF
    
  END SUBROUTINE int4_2Dgrow1

!---------------------------------------------------------------------
! Double length of 1D-chr(8) array 
!---------------------------------------------------------------------
! WARNING - currently assumes that array indexs from 0
! Variables
!	A	:	2D chr(8), array to be reallocated

  SUBROUTINE chr8_1Dgrow(A)
    IMPLICIT NONE
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: A
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: B
    INTEGER :: stat,n,i

    !get shape of A
    n = SIZE(A(:))
    
    !check if allocatable
    IF (.NOT. ALLOCATED(A)) THEN
      WRITE(*,*) "@util:int4_2Dreallocate - array not allocated"
      STOP
    ELSE
      ALLOCATE(B(0:n-1),STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not deallocate array"
        STOP
      END IF

      !temp copy
      B(:) = A(:)

      DEALLOCATE(A,STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not deallocate array"
        STOP
      END IF

      !resize
      ALLOCATE(A(0:2*n-1),STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not allocate array"
        STOP
      END IF
       
      !move the data
      DO i=0,n-1
        A(i) = B(i)
      END DO 
      DEALLOCATE(B)

    END IF
    
  END SUBROUTINE chr8_1Dgrow

!---------------------------------------------------------------------
! Double size of 1st dimension of 2D-real(kind=8) array 
!---------------------------------------------------------------------
! WARNING - currently assumes that array indexs from 0
! Variables
!	A	:	2D real(8), array to be reallocated

  SUBROUTINE real8_2Dgrow1(A)
    IMPLICIT NONE
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: A
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: B
    INTEGER :: stat,n,m,i,j

    !get shape of A
    n = SIZE(A(:,0))
    m = SIZE(A(0,:))
    
    !check if allocatable
    IF (.NOT. ALLOCATED(A)) THEN
      WRITE(*,*) "@util:int4_2Dreallocate - array not allocated"
      STOP
    ELSE
      ALLOCATE(B(0:n-1,0:m-1),STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not deallocate array"
        STOP
      END IF

      !temp copy
      B(:,:) = A(:,:)

      DEALLOCATE(A,STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not deallocate array"
        STOP
      END IF

      !resize
      ALLOCATE(A(0:2*n-1,0:m-1),STAT=stat)
      IF (stat .NE. 0) THEN
        WRITE(*,*) "util:int4_2Dreallocate - could not allocate array"
        STOP
      END IF
       
      !move the data
      DO j=0,m-1
        DO i=0,n-1
          A(i,j) = B(i,j)
        END DO 
      END DO
      DEALLOCATE(B)

    END IF
    
  END SUBROUTINE real8_2Dgrow1

!---------------------------------------------------------------------
! Get number of lines in a file
!---------------------------------------------------------------------
! Variables
!	fname	:	char(20), file name to open	

  INTEGER FUNCTION file_numlines(fname)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: fname
    INTEGER :: nlines,io
    LOGICAL :: ex
    INQUIRE(file=fname,EXIST=ex)
    IF (.NOT. ex) THEN
      WRITE(*,*) "@util:file_numlines -",fname,"does not exist"
      STOP
    END IF
    nlines=0
    io=0
    OPEN(unit=1,file=fname,status='old',access='sequential')
    DO WHILE (io .EQ. 0)
      READ(1,*,iostat=io)
      IF (io .EQ. 0) nlines = nlines + 1
    END DO
    CLOSE(unit=1)
    file_numlines = nlines
  END FUNCTION file_numlines

!---------------------------------------------------------------------
! finds index of ID from list, -1 if not there
!---------------------------------------------------------------------
! WARNING - assume index from zero
! Variables
!       ID      :       char(8), ID of chemical
!       ID_list :       1D int, list of chemical ID's

  INTEGER FUNCTION ID2idx(ID,ID_list)
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
    ID2idx = idx
  END FUNCTION ID2idx

!---------------------------------------------------------------------
END MODULE util
