*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*
      INTEGER FUNCTION CHDIRF (CHNAME)
C
C CERN PROGLIB# Z265    CHDIRF          .VERSION KERNVAX  2.36  910306
C ORIG.  31/07/90, Federico Carminati
C MODIF. 25/02/91, Jamie Shiers,  to work in detached mode
C
C-    Simulate UNIX system call

      INCLUDE '($LNMDEF)'
      STRUCTURE /ITMLST/
         UNION
            MAP
            INTEGER*2 BUFFER_LENGTH
            INTEGER*2 ITEM_CODE
            INTEGER*4 BUFFER_ADDRESS
            INTEGER*4 RETURN_LENGTH_ADDRESS
            ENDMAP
            MAP
            INTEGER*4 END_LIST /0/
            ENDMAP
         END UNION
      END STRUCTURE

      RECORD /ITMLST/ LNM_LIST(2)

      COMMON /SLATE/ ISTAT, ISLATE(39)

      CHARACTER    CHNAME*(*), CHSTRI*255, CHLOGN*127
      LOGICAL      FIRST
      INTEGER*2    LIST2(2)
      INTEGER      LIST4(4), SYS$SETDDIR, SYS$TRNLNM, SYS$CRELNM
      EQUIVALENCE (LIST2, LIST4)
      DATA         LENRET/0/

      CHSTRI=CHNAME
      FIRST=.TRUE.
      LENCHN=INDEX(CHNAME,' ')-1
      IF (LENCHN.LE.0) LENCHN=LEN(CHNAME)
*
*     Cater for users who like the <> characters as directory delimiters
*
      CALL CTRANS('<','[',CHSTRI,1,LENCHN)
      CALL CTRANS('>',']',CHSTRI,1,LENCHN)
  10  IF (CHSTRI(1:1).NE.'[') THEN
        IF (INDEX(CHSTRI,':[') .NE. 0) THEN
*
*---        A directory and a device are specified
*
*         ISTAT=LIB$SET_LOGICAL('SYS$DISK',CHSTRI(1:INDEX(CHSTRI,':')),
*    +                          'LNM$PROCESS')

          LNM_LIST(1).BUFFER_LENGTH           = INDEX(CHSTRI,':')
          LNM_LIST(1).ITEM_CODE               = LNM$_STRING
          LNM_LIST(1).BUFFER_ADDRESS          = %LOC(CHSTRI)
          LNM_LIST(1).RETURN_LENGTH_ADDRESS   = %LOC(LENRET)

          LNM_LIST(2).END_LIST                = 0

          ISTAT=SYS$CRELNM(,'LNM$PROCESS','SYS$DISK',,LNM_LIST)
          IF(.NOT.ISTAT) GO TO 20
          CHSTRI=CHSTRI(INDEX(CHSTRI,'['):)
        ELSEIF (CHSTRI(LENCHN:LENCHN).EQ.':') THEN
*
*---        A device only is specified
*
*         ISTAT=LIB$SET_LOGICAL('SYS$DISK',CHSTRI(:LENCHN),
*    +                          'LNM$PROCESS')

          LNM_LIST(1).BUFFER_LENGTH           = LENCHN
          LNM_LIST(1).ITEM_CODE               = LNM$_STRING
          LNM_LIST(1).BUFFER_ADDRESS          = %LOC(CHSTRI)
          LNM_LIST(1).RETURN_LENGTH_ADDRESS   = %LOC(LENRET)

          LNM_LIST(2).END_LIST                = 0

          ISTAT=SYS$CRELNM(,'LNM$PROCESS','SYS$DISK',,LNM_LIST)
          GO TO 20
        ELSEIF (FIRST) THEN
          FIRST=.FALSE.
          LIST2(1)=LEN(CHLOGN)
          LIST2(2)=LNM$_STRING
          LIST4(2)=%LOC(CHLOGN)
          LIST4(3)=%LOC(LENCH)
          LIST4(4)=0
          ISTAT = SYS$TRNLNM(LNM$M_CASE_BLIND,'LNM$FILE_DEV',
     +           CHSTRI(1:INDEX(CHSTRI,' ')-1),,LIST2)
          IF(.NOT.ISTAT) GO TO 20
          CHSTRI=CHLOGN(:LENCH)
          GO TO 10
        ELSE
          ISTAT=0
          GO TO 20
        END IF
      END IF
      ISTAT = SYS$SETDDIR(CHSTRI, %VAL(0), %VAL(0))
  20  IF (ISTAT)  THEN
        ISTAT = 0
        IRET = 0
      ELSE
        IRET = -1
      ENDIF
      CHDIRF = IRET
      RETURN
      END
