*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzophr55.inc
*
#if defined(CERNLIB_QMIBMD)
C--       Output random with Fortran calls

      CALL ADWRTKD (LUNX,LQ(LWR),NWR,JREC,ISW)
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 2
          IQUEST(3) = ISW
          IQUEST(4) = LUNX
          CALL ZTELL (19,0)
        ENDIF
      GO TO 999
#elif defined(CERNLIB_QMIRTD)
C--       Output random with Fortran calls

      CALL ADWRTKD (LUNX,LQ(LWR),NWR,JREC,ISW)
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 2
          IQUEST(3) = ISW
          IQUEST(4) = LUNX
          CALL ZTELL (19,0)
        ENDIF
      GO TO 999
#elif 1
C--       Output random with Fortran calls

      WRITE (LUNX,REC=JREC,IOSTAT=ISW)  (LQ(LWR+J-1),J=1,NWR)
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 2
          IQUEST(3) = ISW
          IQUEST(4) = LUNX
          CALL ZTELL (19,0)
        ENDIF
      GO TO 999

#endif
