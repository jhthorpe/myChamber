PROGRAM test_gas
  USE gas
  IMPLICIT NONE

  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: gas_rxns
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:)  :: ID_list
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: gas_coef
  INTEGER  :: gas_num

  CALL gas_read(ID_list,gas_rxns,gas_coef,gas_num)

END PROGRAM test_gas
