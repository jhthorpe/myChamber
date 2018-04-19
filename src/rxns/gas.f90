!Module for auxilary functions for gas phase kinetics reactions

MODULE gas
  USE util
  IMPLICIT NONE

  CONTAINS
!=====================================================================
!                       SUBROUTINES 

!---------------------------------------------------------------------
! Read gas phase input 
! Variables
!	ID_list		:	1D char(4), ID's of chemical species
!	gas_rxns	:	2D int, reactions matrix
!	gas_coef	:	2D dp, coef of general rate eq
!	gas_num		:	int, number of gas phase reactions

  SUBROUTINE gas_read(ID_list,gas_rxns,gas_coef,gas_num)
    IMPLICIT NONE
    !inout
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: gas_rxns
    INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ID_list
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: gas_coef
    INTEGER, INTENT(INOUT) :: gas_num

    !internal
    INTEGER :: i,j,k,l

    ALLOCATE(gas_rxns(0:7,0:3))
    WRITE(*,*) SHAPE(gas_rxns)

    CALL int4_2Dgrow1(gas_rxns)
    WRITE(*,*) SHAPE(gas_rxns)


  END SUBROUTINE gas_read
    
  

END MODULE gas
