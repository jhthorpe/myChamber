! Controling function for myChamber

PROGRAM myChamber
  USE gas 
  USE integral
  USE pop
  IMPLICIT NONE
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: gas_rxns
  CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)  :: ID_list
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: gas_coef
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: conc
  INTEGER  :: gas_nrxn, nID
  LOGICAL :: ex
  INTEGER :: i

  WRITE(*,*) "Starting myChamber"

  CALL EXECUTE_COMMAND_LINE('parse_chamber')
  INQUIRE(file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from parse_chamber"
    STOP
  END IF

  !read in gas phase reaction data
  CALL gas_read(ID_list,gas_rxns,gas_coef,gas_nrxn,nID)

  !read in photolysis reaction data

  !initialize concentrations
  ALLOCATE(conc(0:nID-1))
  CALL pop_init(conc,ID_list)

  !integrate rate equations

END PROGRAM myChamber
