! Module that contains useful utility functions

MODULE util
  IMPLICIT NONE

  CONTAINS

!---------------------------------------------------------------------
! Change size of int(kind=4) array 
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

END MODULE util
