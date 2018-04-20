! Controling function for myChamber

PROGRAM myChamber
  USE gas
  IMPLICIT NONE
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: gas_rxns
  CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)  :: ID_list
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: gas_coef
  INTEGER  :: gas_nrxn, gas_nID
  LOGICAL :: ex

  WRITE(*,*) "Starting myChamber"

  CALL EXECUTE_COMMAND_LINE('parse_chamber')
  INQUIRE(file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from parse_chamber"
    STOP
  END IF

  CALL gas_read(ID_list,gas_rxns,gas_coef,gas_nrxn,gas_nID)

END PROGRAM myChamber
