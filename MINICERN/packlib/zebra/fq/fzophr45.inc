*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:10  mclareni
* Zebra
*
*
*
* fzophr45.inc
*
#if defined(CERNLIB_QMIBMD)
C--       Output sequential with Fortran calls
C--       For AUTODBL machines
      IF (IACMOX.EQ.0)             GO TO 48
      CALL ADWRTK (LUNX, LQ(LWR), NWR, ISW)
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 1
          IQUEST(3) = ISW
          IQUEST(4) = LUNX
          CALL ZTELL (19,0)
        ENDIF
      GO TO 999
#elif defined(CERNLIB_QMIRTD)
C--       Output sequential with Fortran calls
C--       For AUTODBL machines
  48  CALL ADWRTK (LUNX, LQ(LWR), NWR, ISW)
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 1
          IQUEST(3) = ISW
          IQUEST(4) = LUNX
          CALL ZTELL (19,0)
        ENDIF
      GO TO 999
#elif 1
   48 WRITE (LUNX,IOSTAT=ISW)  (LQ(LWR+J-1),J=1,NWR)
      IF (ISW.NE.0)  THEN
          IQUEST(1) = 19
          IQUEST(2) = 1
          IQUEST(3) = ISW
          IQUEST(4) = LUNX
          CALL ZTELL (19,0)
        ENDIF
      GO TO 999
#endif
