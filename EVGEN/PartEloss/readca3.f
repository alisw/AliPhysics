      SUBROUTINE ELOSS_READCA3
      REAL*8           xx(400), da(30), ca(30,260), rrr(30)
      COMMON /data/    xx, da, ca, rrr
*
      CHARACTER*100 CHROOT  
      CHARACTER*100 FILNAM
      CHROOT=' '
      CALL GETENVF('ALICE_ROOT',CHROOT)
      LNROOT = LNBLNK(CHROOT)

      IF(LNROOT.LE.0) THEN
         FILNAM='DATAcca3.dat'
      ELSE
         FILNAM=CHROOT(1:LNROOT)//'/FAST/eloss/DATAcca3.dat'
      ENDIF

      OPEN(UNIT=20,FILE=FILNAM,STATUS='OLD',ERR=90)
      nn = 1
 100  read (20,*,end=110) xx(nn), ca(1,nn), ca(2,nn), ca(3,nn),
     +     ca(4,nn), ca(5,nn), ca(6,nn), ca(7,nn), ca(8,nn),
     +     ca(9,nn), ca(10,nn), ca(11,nn), ca(12,nn), ca(13,nn),
     +     ca(14,nn), ca(15,nn), ca(16,nn), ca(17,nn), ca(18,nn),
     +     ca(19,nn), ca(20,nn), ca(21,nn), ca(22,nn), ca(23,nn),
     +     ca(24,nn), ca(25,nn), ca(26,nn), ca(27,nn), ca(28,nn),
     +     ca(29,nn), ca(30,nn)
*         print*, 0.005*(nn-1), ca(1,nn), ca(2,nn), ca(3,nn),
*     +     ca(4,nn), ca(5,nn), ca(6,nn), ca(7,nn), ca(8,nn),
*     +     ca(9,nn), ca(10,nn), ca(11,nn), ca(12,nn), ca(13,nn),
*     +     ca(14,nn), ca(15,nn), ca(16,nn), ca(17,nn), ca(18,nn),
*     +     ca(19,nn), ca(20,nn), ca(21,nn), ca(22,nn), ca(23,nn),
*     +     ca(24,nn), ca(25,nn), ca(26,nn), ca(27,nn), ca(28,nn),
*     +     ca(29,nn), ca(30,nn)
         nn = nn + 1
         goto 100
 110     continue
      close(20)
*
      IF(LNROOT.LE.0) THEN
         FILNAM='DATAdca3.dat'
      ELSE
         FILNAM=CHROOT(1:LNROOT)//'/FAST/eloss/DATAdca3.dat'
      ENDIF

      OPEN(UNIT=21,FILE=FILNAM,STATUS='OLD',ERR=90)
      nn = 1
 101  read (21,*,end=111) rrr(nn), da(nn)
*         print*, rrr(nn), da(nn)
         nn = nn + 1
         goto 101
 111     continue
      close(21)
*
      goto 888
 90   PRINT*, 'input - output error'
 888  continue
      RETURN 
      END











