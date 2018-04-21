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
    INTEGER :: mID,mrxn,nlines,fID
    INTEGER :: i,j,k

    
    

  END SUBROUTINE photo_read
  

END MODULE photo
