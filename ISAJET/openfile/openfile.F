      SUBROUTINE OPENFILES
      CHARACTER*100 CHROOT  
      CHARACTER*100 FILNAM
      LOGICAL EXISTS
      CHROOT=' '
      CALL GETENVF('ALICE_ROOT',CHROOT)
      LNROOT = LNBLNK(CHROOT)
      IF(LNROOT.LE.0) THEN
         FILNAM='myjob.par'
      ELSE
         FILNAM=CHROOT(1:LNROOT)//'/ISAJET/data/myjob.par'
      ENDIF

      INQUIRE(FILE=FILNAM,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
         PRINT*,'**********************************'
         PRINT*,'*        I S A J E T             *'
         PRINT*,'*        -----------             *'
         PRINT*,'*   File myjob.par   not found   *'
         PRINT*,'*         Program STOP           *'
         PRINT*,'*   Check CERN_ROOT environment  *'
         PRINT*,'*           variable             *'
         PRINT*,'**********************************'
         STOP
      ENDIF

      OPEN(51, FILE=FILNAM
     &, STATUS="OLD",FORM = "FORMATTED") 
      OPEN(52, FILE="myjob.dat", STATUS="OLD",FORM = "UNFORMATTED") 

      OPEN(53, FILE="myjob.lis", STATUS="OLD",FORM = "FORMATTED") 


      FILNAM=CHROOT(1:LNROOT)//'/ISAJET/data/decay.cpp'
      INQUIRE(FILE=FILNAM,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
         PRINT*,'**********************************'
         PRINT*,'*        I S A J E T             *'
         PRINT*,'*        -----------             *'
         PRINT*,'*   File decay.cpp   not found   *'
         PRINT*,'*         Program STOP           *'
         PRINT*,'*   Check CERN_ROOT environment  *'
         PRINT*,'*           variable             *'
         PRINT*,'**********************************'
         STOP
      ENDIF

      OPEN(54, FILE=FILNAM
     &, STATUS="OLD",FORM = "FORMATTED") 
      RETURN 
      END
