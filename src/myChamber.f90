! Controling function for myChamber

PROGRAM myChamber
  USE gas 
  USE integral
  USE pop
  USE photo
  IMPLICIT NONE
  INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: gas_rxns,photo_rxns
  CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)  :: ID_list
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)  :: gas_coef,photo_QY,photo_CS,photo_AF
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: conc
  REAl(KIND=8), DIMENSION(0:7) :: options
  INTEGER  :: gas_nrxn, photo_nrxn,nID,photo_ID
  LOGICAL :: ex
  INTEGER :: i

  WRITE(*,*) "Starting myChamber"

  CALL EXECUTE_COMMAND_LINE('parse_chamber')
  INQUIRE(file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from parse_chamber"
    STOP
  END IF


  nID = 0
  !read in gas phase reaction data
  CALL gas_read(ID_list,gas_rxns,gas_coef,gas_nrxn,nID)
  !read in photolysis reaction data
  CALL photo_read(ID_list,photo_rxns,photo_QY,photo_CS,photo_AF,photo_nrxn,nID,photo_ID)

  WRITE(*,*) "=============================="


  !initialize concentrations
  ALLOCATE(conc(0:nID-1))
  CALL pop_init(conc,ID_list)

  !read in options
  OPEN(unit=1,file='jobinfo',access='sequential',status='old')
  READ(1,*) options
  CLOSE(unit=1)

  !integrate rate equations
  CALL integrate(ID_list,nID,gas_rxns,gas_nrxn,gas_coef,photo_rxns,photo_nrxn,photo_QY,photo_CS,photo_AF,photo_ID,&
                 conc,options)
  WRITE(*,*) 
  WRITE(*,*) "Simulation complete"
  WRITE(*,*) 

END PROGRAM myChamber

