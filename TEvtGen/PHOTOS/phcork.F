


      SUBROUTINE PHCORK(MODCOR)
      implicit none
C.----------------------------------------------------------------------
C.
C.    PHCORK: corrects kinmatics of subbranch needed if host program
C.            produces events with the shaky momentum conservation
C
C.    Input Parameters:   Common /PHOEVT/, MODCOR
C.                        MODCOR >0 type of action
C.                               =1 no action
C.                               =2 corrects energy from mass
C.                               =3 corrects mass from energy
C.                               =4 corrects energy from mass for 
C.                                  particles up to .4 GeV mass, 
C.                                  for heavier ones corrects mass,
C.                               =0 execution mode 
C.
C.    Output Parameters:  corrected /PHOEVT/
C.
C.    Author(s):  P.Golonka, Z. Was               Created at:  01/02/99
C.						  Modified  :  08/02/99
C.----------------------------------------------------------------------
      INTEGER NMXPHO
      PARAMETER (NMXPHO=10000)
      
      REAL*8 M,P2,PX,PY,PZ,E,EN,MCUT
      INTEGER MODCOR,MODOP,I,IEV,IPRINT
      INTEGER IDPHO,ISTPHO,JDAPHO,JMOPHO,NEVPHO,NPHO
      REAL*8 PPHO,VPHO
      COMMON/PHOEVT/NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
     &JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)

      INTEGER PHLUN
      COMMON/PHOLUN/PHLUN

      COMMON /PHNUM/ IEV
      SAVE MODOP
      DATA MODOP  /0/
      SAVE IPRINT
      DATA IPRINT /0/
      SAVE MCUT
      IF (MODCOR.NE.0) THEN 
C       INITIALIZATION
        MODOP=MODCOR

        WRITE(PHLUN,*) 'Message from PHCORK(MODCOR):: initialization'
        IF     (MODOP.EQ.1) THEN
          WRITE(PHLUN,*) 'MODOP=1 -- no corrections on event: DEFAULT' 
        ELSEIF (MODOP.EQ.2) THEN
          WRITE(PHLUN,*) 'MODOP=2 -- corrects Energy from mass'
        ELSEIF (MODOP.EQ.3) THEN
          WRITE(PHLUN,*) 'MODOP=3 -- corrects mass from Energy'
        ELSEIF (MODOP.EQ.4) THEN
          WRITE(PHLUN,*) 'MODOP=4 -- corrects Energy from mass to Mcut'
          WRITE(PHLUN,*) '           and mass from  energy above  Mcut '
          MCUT=0.4
          WRITE(PHLUN,*) 'Mcut=',MCUT,'GeV'
        ELSE
          WRITE(PHLUN,*) 'PHCORK wrong MODCOR=',MODCOR
          STOP
        ENDIF
        RETURN
      ENDIF

      IF (MODOP.EQ.0.AND.MODCOR.EQ.0) THEN
        WRITE(PHLUN,*) 'PHCORK lack of initialization'
        STOP
      ENDIF

C execution mode
C ==============
C ============== 

     
        PX=0
        PY=0
        PZ=0
        E =0

      IF    (MODOP.EQ.1) THEN
C     -----------------------
C       In this case we do nothing
        RETURN
      ELSEIF(MODOP.EQ.2) THEN
C     -----------------------
CC      lets loop thru all daughters and correct their energies 
CC      according to E^2=p^2+m^2

       DO I=3,NPHO
         
         PX=PX+PPHO(1,I)
         PY=PY+PPHO(2,I)
         PZ=PZ+PPHO(3,I)

         P2=PPHO(1,I)**2+PPHO(2,I)**2+PPHO(3,I)**2

         EN=SQRT( PPHO(5,I)**2 + P2)
         
         IF (IPRINT.EQ.1)
     &   WRITE(PHLUN,*) 'CORRECTING ENERGY OF ',I,':',
     &        PPHO(4,I),'=>',EN

         PPHO(4,I)=EN
         E = E+PPHO(4,I)

       ENDDO
      
      ELSEIF(MODOP.EQ.3) THEN
C     -----------------------

CC      lets loop thru all daughters and correct their masses 
CC      according to E^2=p^2+m^2

       DO I=3,NPHO
         
         PX=PX+PPHO(1,I)
         PY=PY+PPHO(2,I)
         PZ=PZ+PPHO(3,I)
         E = E+PPHO(4,I)

         P2=PPHO(1,I)**2+PPHO(2,I)**2+PPHO(3,I)**2

         M=SQRT(ABS( PPHO(4,I)**2 - P2))

         IF (IPRINT.EQ.1)
     &   WRITE(PHLUN,*) 'CORRECTING MASS OF ',I,':',
     &        PPHO(5,I),'=>',M

         PPHO(5,I)=M

       ENDDO
      

      ELSEIF(MODOP.EQ.4) THEN
C     -----------------------
            
CC      lets loop thru all daughters and correct their masses 
CC      or energies according to E^2=p^2+m^2

       DO I=3,NPHO
         
         PX=PX+PPHO(1,I)
         PY=PY+PPHO(2,I)
         PZ=PZ+PPHO(3,I)

         P2=PPHO(1,I)**2+PPHO(2,I)**2+PPHO(3,I)**2

         M=SQRT(ABS( PPHO(4,I)**2 - P2))

         IF (M.GT.MCUT) THEN
          IF (IPRINT.EQ.1)
     &    WRITE(PHLUN,*) 'CORRECTING MASS OF ',I,':',
     &         PPHO(5,I),'=>',M
          PPHO(5,I)=M
          E = E+PPHO(4,I)
         ELSE

          EN=SQRT( PPHO(5,I)**2 + P2)

         IF (IPRINT.EQ.1)
     &    WRITE(PHLUN,*) 'CORRECTING ENERGY OF ',I,':',
     &        PPHO(4,I),'=>',EN

          PPHO(4,I)=EN
          E = E+PPHO(4,I)
         ENDIF

       ENDDO
      ENDIF
C     -----

       IF (IPRINT.EQ.1) THEN
        WRITE(PHLUN,*) 'CORRECTING MOTHER'
        WRITE(PHLUN,*) 'PX:',PPHO(1,1),'=>',PX-PPHO(1,2)
        WRITE(PHLUN,*) 'PY:',PPHO(2,1),'=>',PY-PPHO(2,2)
        WRITE(PHLUN,*) 'PZ:',PPHO(3,1),'=>',PZ-PPHO(3,2)
        WRITE(PHLUN,*) ' E:',PPHO(4,1),'=>',E-PPHO(4,2)
       ENDIF

       PPHO(1,1)=PX-PPHO(1,2)
       PPHO(2,1)=PY-PPHO(2,2)
       PPHO(3,1)=PZ-PPHO(3,2)
       PPHO(4,1)=E -PPHO(4,2)

       P2=PPHO(1,1)**2+PPHO(2,1)**2+PPHO(3,1)**2

       IF (PPHO(4,1)**2.GT.P2) THEN
          M=SQRT( PPHO(4,1)**2 - P2 )
          IF (IPRINT.EQ.1)
     &    WRITE(PHLUN,*) ' M:',PPHO(5,1),'=>',M
          PPHO(5,1)=M
       ENDIF

      CALL PHLUPA(25)

      END
