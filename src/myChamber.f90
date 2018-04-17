! Controling function for myChamber

PROGRAM myChamber
  IMPLICIT NONE
  LOGICAL :: ex

  WRITE(*,*) "Starting myChamber"

  CALL EXECUTE_COMMAND_LINE('parse_chamber')
  INQUIRE(file='error',EXIST=ex)
  IF (ex) THEN
    WRITE(*,*) "Error from parse_chamber"
    STOP
  END IF

END PROGRAM myChamber
