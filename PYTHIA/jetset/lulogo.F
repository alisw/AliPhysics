 
C********************************************************************* 
 
      SUBROUTINE LULOGO 
 
C...Purpose: to write logo for JETSET and PYTHIA programs. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /LUDAT1/ 
      SAVE /PYPARS/ 
      CHARACTER MONTH(12)*3, LOGO(48)*32, REFER(22)*36, LINE*79, 
     &VERS*1, SUBV*3, DATE*2, YEAR*4 
 
C...Data on months, logo, titles, and references. 
      DATA MONTH/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', 
     &'Oct','Nov','Dec'/ 
      DATA (LOGO(J),J=1,10)/ 
     &'PPP  Y   Y TTTTT H   H III   A  ', 
     &'P  P  Y Y    T   H   H  I   A A ', 
     &'PPP    Y     T   HHHHH  I  AAAAA', 
     &'P      Y     T   H   H  I  A   A', 
     &'P      Y     T   H   H III A   A', 
     &'JJJJ EEEE TTTTT  SSS  EEEE TTTTT', 
     &'   J E      T   S     E      T  ', 
     &'   J EEE    T    SSS  EEE    T  ', 
     &'J  J E      T       S E      T  ', 
     &' JJ  EEEE   T    SSS  EEEE   T  '/ 
      DATA (LOGO(J),J=11,29)/ 
     &'            *......*            ', 
     &'       *:::!!:::::::::::*       ', 
     &'    *::::::!!::::::::::::::*    ', 
     &'  *::::::::!!::::::::::::::::*  ', 
     &' *:::::::::!!:::::::::::::::::* ', 
     &' *:::::::::!!:::::::::::::::::* ', 
     &'  *::::::::!!::::::::::::::::*! ', 
     &'    *::::::!!::::::::::::::* !! ', 
     &'    !! *:::!!:::::::::::*    !! ', 
     &'    !!     !* -><- *         !! ', 
     &'    !!     !!                !! ', 
     &'    !!     !!                !! ', 
     &'    !!                       !! ', 
     &'    !!        ep             !! ', 
     &'    !!                       !! ', 
     &'    !!                 pp    !! ', 
     &'    !!   e+e-                !! ', 
     &'    !!                       !! ', 
     &'    !!                          '/ 
      DATA (LOGO(J),J=30,48)/ 
     &'Welcome to the Lund Monte Carlo!', 
     &'                                ', 
     &'  This is PYTHIA version x.xxx  ', 
     &'Last date of change: xx xxx 199x', 
     &'                                ', 
     &'  This is JETSET version x.xxx  ', 
     &'Last date of change: xx xxx 199x', 
     &'                                ', 
     &'          Main author:          ', 
     &'       Torbjorn Sjostrand       ', 
     &' Dept. of theoretical physics 2 ', 
     &'       University of Lund       ', 
     &'         Solvegatan 14A         ', 
     &'      S-223 62 Lund, Sweden     ', 
     &'   phone: +46 - 46 - 222 48 16  ', 
     &'   E-mail: torbjorn@thep.lu.se  ', 
     &'                                ', 
     &'  Copyright Torbjorn Sjostrand  ', 
     &'     and CERN, Geneva 1993      '/ 
      DATA (REFER(J),J=1,6)/ 
     &'The latest program versions and docu',
     &'mentation is found on WWW address   ',
     &'http://thep.lu.se/tf2/staff/torbjorn',
     &'/Welcome.html                       ',
     &'                                    ',
     &'                                    '/
      DATA (REFER(J),J=7,22)/ 
     &'When you cite these programs, priori', 
     &'ty should always be given to the    ', 
     &'latest published description. Curren', 
     &'tly this is                         ', 
     &'T. Sjostrand, Computer Physics Commu', 
     &'n. 82 (1994) 74.                    ', 
     &'The most recent long description (un', 
     &'published) is                       ', 
     &'T. Sjostrand, LU TP 95-20 and CERN-T',
     &'H.7112/93 (revised August 1995).    ', 
     &'Also remember that the programs, to ', 
     &'a large extent, represent original  ', 
     &'physics research. Other publications', 
     &' of special relevance to your       ', 
     &'studies may therefore deserve separa', 
     &'te mention.                         '/ 
 
C...Check if PYTHIA linked. 
      IF(MSTP(183)/10.NE.199) THEN 
        LOGO(32)=' Warning: PYTHIA is not loaded! ' 
        LOGO(33)='Did you remember to link PYDATA?' 
      ELSE 
        WRITE(VERS,'(I1)') MSTP(181) 
        LOGO(32)(26:26)=VERS 
        WRITE(SUBV,'(I3)') MSTP(182) 
        LOGO(32)(28:30)=SUBV 
        WRITE(DATE,'(I2)') MSTP(185) 
        LOGO(33)(22:23)=DATE 
        LOGO(33)(25:27)=MONTH(MSTP(184)) 
        WRITE(YEAR,'(I4)') MSTP(183) 
        LOGO(33)(29:32)=YEAR 
      ENDIF 
 
C...Check if JETSET linked. 
      IF(MSTU(183)/10.NE.199) THEN 
        LOGO(35)='  Error: JETSET is not loaded!  ' 
        LOGO(36)='Did you remember to link LUDATA?' 
      ELSE 
        WRITE(VERS,'(I1)') MSTU(181) 
        LOGO(35)(26:26)=VERS 
        WRITE(SUBV,'(I3)') MSTU(182) 
        LOGO(35)(28:30)=SUBV 
        WRITE(DATE,'(I2)') MSTU(185) 
        LOGO(36)(22:23)=DATE 
        LOGO(36)(25:27)=MONTH(MSTU(184)) 
        WRITE(YEAR,'(I4)') MSTU(183) 
        LOGO(36)(29:32)=YEAR 
      ENDIF 
 
C...Loop over lines in header. Define page feed and side borders. 
      DO 100 ILIN=1,48 
      LINE=' ' 
      IF(ILIN.EQ.1) THEN 
        LINE(1:1)='1' 
      ELSE 
        LINE(2:3)='**' 
        LINE(78:79)='**' 
      ENDIF 
 
C...Separator lines and logos. 
      IF(ILIN.EQ.2.OR.ILIN.EQ.3.OR.ILIN.EQ.47.OR.ILIN.EQ.48) THEN 
        LINE(4:77)='***********************************************'// 
     &  '***************************' 
      ELSEIF(ILIN.GE.6.AND.ILIN.LE.10) THEN 
        LINE(6:37)=LOGO(ILIN-5) 
        LINE(44:75)=LOGO(ILIN) 
      ELSEIF(ILIN.GE.13.AND.ILIN.LE.31) THEN 
        LINE(6:37)=LOGO(ILIN-2) 
        LINE(44:75)=LOGO(ILIN+17) 
      ELSEIF(ILIN.GE.34.AND.ILIN.LE.44) THEN 
        LINE(5:40)=REFER(2*ILIN-67) 
        LINE(41:76)=REFER(2*ILIN-66) 
      ENDIF 
 
C...Write lines to appropriate unit. 
      IF(MSTU(183)/10.EQ.199) THEN 
        WRITE(MSTU(11),'(A79)') LINE 
      ELSE 
        WRITE(*,'(A79)') LINE 
      ENDIF 
  100 CONTINUE 
 
C...Check that matching subversions are linked. 
      IF(MSTU(183)/10.EQ.199.AND.MSTP(183)/10.EQ.199) THEN 
        IF(MSTU(182).LT.MSTP(186)) WRITE(MSTU(11), 
     &  '(/'' Warning: JETSET subversion too old for PYTHIA''/)') 
        IF(MSTP(182).LT.MSTU(186)) WRITE(MSTU(11), 
     &  '(/'' Warning: PYTHIA subversion too old for JETSET''/)') 
      ENDIF 
 
      RETURN 
      END 
