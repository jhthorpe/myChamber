!	Parse input for myChamber 

!=====================================================================
!                       MAIN 

PROGRAM parser

  IMPLICIT NONE

  !See options.txt for simulation options

  !Variables
  REAL(KIND=8), DIMENSION(0:7) :: options
  CHARACTER(LEN=4), DIMENSION(0:7) :: opt_name

  !internal
  CHARACTER(LEN=20),DIMENSION(0:1) :: line
  CHARACTER(LEN=20) :: str
  REAL(KIND=8) :: val,val1,val2,val3
  INTEGER :: i

  !defaults
  options = [0.0,0.0,0.0,0.0,0.0,1.0,0.5,298.15]
  opt_name = ['SUN ','LAT ','LON ','DNUM','TSTR','TLEN','TSTP','TEMP']

  WRITE(*,*) 
  WRITE(*,*) "parse called" 
  WRITE(*,*) "Input Parameters"
 
  !get parameters
  OPEN (unit=1,file="input.dat",status="old",access="sequential")

  !1) SIMULATION parameters
  READ(1,*)
  READ(1,*)
  READ(1,*)
  READ(1,*) line    
  options(0) = getSUN(line(1))
  READ(1,*) str, val
  options(1) = getLAT(val)
  READ(1,*) str, val
  options(2) = getLON(val)
  READ(1,*) str, val
  options(3) = getDNUM(val)
  READ(1,*) str, val1,val2,val3
  options(4) = getTSTR(val1,val2,val3)
  READ(1,*) str, val1,val2,val3
  options(5) = getTLEN(val1,val2,val3)
  READ(1,*) str, val1,val2,val3
  options(6) = getTSTP(val1,val2,val3)
  READ(1,*) 
  READ(1,*)
  READ(1,*) str, val
  options(7) = getTEMP(val)
  
  
  CLOSE (unit=1)

  !print parameters
  WRITE(*,*)
  WRITE(*,*) "Paremeters"
  DO i=0,SIZE(options)-1
    WRITE(*,*) opt_name(i),":",options(i)
  END DO
  WRITE(*,*) 
  WRITE(*,*) "=============================="

  !Write to jobinfo file
  OPEN(unit=2,file='jobinfo',status='replace',access='sequential')
  WRITE(2,*) options
  CLOSE(unit=2,status='keep')

  CONTAINS
  
!=====================================================================
!                       FUNCTIONS 

!---------------------------------------------------------------------
! Get SUN treatment
  REAL(KIND=8) FUNCTION getSUN(chr)
    IMPLICIT NONE
    CHARACTER(LEN=20),INTENT(IN) :: chr
    getSUN = 0.0D0
    IF (chr .EQ. 'NONE') THEN
      getSUN = 0.0D0
    ELSE IF (chr .EQ. 'INTERNAL') THEN
      getSUN = 1.0D0
    ELSE IF (chr .EQ. 'EXTERNAL') THEN
      getSUN = 2.0D0
    ELSE
      WRITE(*,*) "WARNING - defaulting to SUN= NONE"
    END IF
  END FUNCTION getSUN

!---------------------------------------------------------------------
! Get Lattitude
  REAL(KIND=8) FUNCTION getLAT(val)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: val
    REAL(KIND=8) :: lat
    lat = val
    IF (val .LT. 0.0D0 .OR. val .GT. 180.0D0) THEN
      WRITE(*,*) "WARNING - defaulting to LAT= 0.0"
      lat = 0.0D0
    END IF
    getLAT = lat
  END FUNCTION getLAT
!---------------------------------------------------------------------
! Get Longitude
  REAL(KIND=8) FUNCTION getLON(val)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: val
    REAL(KIND=8) :: lon
    lon = val
    IF (val .LT. -180.0D0 .OR. val .GT. 180.0D0) THEN
      WRITE(*,*) "WARNING - defaulting to LON= 0.0"
      lon = 0.0D0
    END IF
    getLON = lon
  END FUNCTION getLON
!---------------------------------------------------------------------
! Get Day Number 
  REAL(KIND=8) FUNCTION getDNUM(val)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: val
    REAL(KIND=8) :: DNUM
    DNUM = NINT(val)*1.0D0
    IF (val .LT. 0 .OR. val .GT. 365) THEN
      WRITE(*,*) "WARNING - defaulting to DNUM= 0"
      DNUM = 0.0D0
    END IF
    getDNUM = DNUM
  END FUNCTION getDNUM
!---------------------------------------------------------------------
! Get Start Time 
  REAL(KIND=8) FUNCTION getTSTR(hr,mn,sc)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in) :: hr,mn,sc
    REAL(KIND=8) :: time
    time = hr*360.0D0
    time = mn*60.0D0 + time
    time = sc + time
    IF (time .GT. 86400.0D0) THEN
      WRITE(*,*) "WARNING - defaulting to TSTR= 00 00 00"
      time = 0.0D0
    END IF
    getTSTR = time
  END FUNCTION getTSTR
!---------------------------------------------------------------------
! Get Simulation length
  REAL(KIND=8) FUNCTION getTLEN(hr,mn,sc)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in) :: hr,mn,sc
    REAL(KIND=8) :: time
    time = hr*360.0D0
    time = mn*60.0D0 + time
    time = sc + time
    getTLEN = time
  END FUNCTION getTLEN
!---------------------------------------------------------------------
! Get time step 
  REAL(KIND=8) FUNCTION getTSTP(hr,mn,sc)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in) :: hr,mn,sc
    REAL(KIND=8) :: time
    time = hr*360.0D0
    time = mn*60.0D0 + time
    time = sc + time
    getTSTP = time
  END FUNCTION getTSTP
!---------------------------------------------------------------------
! Get temperatures 
  REAL(KIND=8) FUNCTION getTEMP(val)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(in) :: val 
    REAL(KIND=8) :: temp
    temp = 298.15
    IF (val .LT. 200.0 .OR. val .GT. 400.0) THEN
      WRITE(*,*) "WARNING- defaulting to TEMP= 298.15 K" 
    ELSE
      temp = val 
    END IF
    getTEMP = temp
  END FUNCTION getTEMP


END PROGRAM parser
