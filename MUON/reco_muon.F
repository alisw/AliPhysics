*  25/5/99
**  Authors J.P. Cussonneau & P. Lautridou 

************** commentaires **********************

* bon muon = tous les hits de la trace (pas forcement 28) proviennent de muons
* ghost = bonne trace dans laquelle certains hits proviennent délectrons, de muons ou sont ambigus
 
* ALPHATOP
* alpha to P

* EFF
* efficacite

* HHIT 
* h au vertex

* HTOP
* h to P

* INDEXMAX
* nombre de candidats a ordonner (a partir de l'ímpulsion)

* INDEXTAB
* pour recuperer les candidats

* ISTAT
* =1 si bon muon, =2 si ghost, =0 autrement

* ITCHECK
* =1 si bonne trace, =0 autrement

* IT_LIST
* permet a partir du numero de la trace de retrouver la numero du hit

* IT_NP
* compte le nombre de plans touches par trace

* ITRACK
* permet de retrouver le numero de la trace a partir du numero de hit

* ITTROUGH
* pour une trace et une station donnee, dit si la trace est passee dans la chambre

* IVERTEX
* =0 si point (0,0) du vertex impose, sinon coordonnes libres (pour le fit)

* JCAN 
* numero du hit pour une station et un candidat

* JCANTYP
* nombre de hits par trace au niveau des stations 4 et 5 (3 ou 4)

* JJOUT
* numero de hit associe a une chambre et a une trace

* NMUONALL 
* nombre de bons muons trouves

* NERR
* = nombre de fois ou lón ná pas trouve la bonne trace dans la station

* NERRALL
* nombre de cas ou pas de hit trouve par station

* NGHOSTALL
* nombre de fantomes trouves

* NRES
* nombre de resonances dans lácceptance

* NRESF
* nombre de resonances trouvees

* NTRACFALL
* nombre de traces totales trouvees

* NTRMUALL
* nombre total de muons dans lácceptance


****************************************************************
      SUBROUTINE reconstmuon(IFIT,IDEBUGC,NEV,IDRES,IREADGEANT)
****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/DEBEVT/IDEBUG
      common/dstation/idstation

*      idsation = 8
      idstation = 20

      IDEBUG=IDEBUGC

** Read events          
      CALL RECO_READEVT(NEV,IDRES,IREADGEANT)


** Trackfinding         
      CALL RECO_TRACKF(IDRES,IREADGEANT)

** Precision fit  
      IF (IFIT.EQ.1) THEN
         CALL RECO_PRECISION
      ENDIF

** Calculate         
      CALL RECO_SUM  

      END

********************************************
      BLOCK DATA
********************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NBSTATION=5)
* --
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
      
* ZPLANE = position de la premiere chambre de chaque station
* exemple (pour la premiere station) : 
*           zch (dans AliMUONv0) = 540 (position du centre des 2 chambres) 
*           implique :
*           ZPLANE = -530 si DZ_PL = 20cm

*     dstation = 8cm :
*         DATA DZ_PL/8.,8.,8.,8.,8./   
*         DATA ZPLANE/-511.,-686.0,-971.,-1245.,-1445./ 

*     dstation=20cm :
         DATA DZ_PL/20.,20.,20.,20.,20./     
         DATA ZPLANE/-518.,-680.,-965.,-1239.,-1439./ 
*      end if

*      DATA ZCOIL,ZMAGEND/-825.0,-1125./  ! Constant field 3 Tm
      DATA ZCOIL,ZMAGEND/-805.0,-1233./ ! CCC magn. field map M.B 
**      DATA ZCOIL,ZMAGEND/-795.1,-1242.9/ ! CCC magn. field map M.B 

* --
      END

********************************************
      SUBROUTINE cutpxz(spxzcut)
********************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/ACUTPXZ/ACUTPXZ

      ACUTPXZ = SPXZCUT

      END

********************************************
      SUBROUTINE sigmacut(ssigcut)
********************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP

      SIGCUT = SSIGCUT

      END

********************************************
      SUBROUTINE xpreci(sxprec)
********************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP

      XPREC = SXPREC
 
      END

********************************************
      SUBROUTINE ypreci(syprec)
********************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP

      YPREC = SYPREC

      END

*************************************************************************
      SUBROUTINE RECO_INIT(seff,sb0,sbl3)
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CALL TRACKF_INIT(seff,sb0,sbl3)

      CALL PREC_INIT 

      RETURN
      END

*************************************************************************
      SUBROUTINE TRACKF_INIT(seff,sb0,sbl3)
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
**
      PARAMETER(NBSTATION=5,NTRMAX=500) 
**      
      COMMON/REVENT/IEVBKGI,NBKGMAX,MAXUPSEV
**      
      COMMON/MAGNET/BL3,B0
**      
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
**
      COMMON/FILED/FILERES,FILEBKG,FILEOUT,FILEMIN
**
      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP
**
      COMMON/TRACKSUM/NRES(5),NRESF,NTRMUALL,NMUONALL,NGHOSTALL,
     &     NTRACKFALL,NERRALL(NBSTATION),IR               
**      
      COMMON/ACUTPXZ/ACUTPXZ
**
      COMMON/DEBEVT/IDEBUG
**
      CALL HIST_CREATE
*
      EFF = SEFF
      B0 = SB0
      BL3 = SBL3

      PXZCUT = ACUTPXZ 
      AMAGLEN = ZMAGEND-ZCOIL
      ZM = AMAGLEN/2.+ZCOIL
      ALPHATOP = 0.01*0.3*B0*ABS(AMAGLEN)
      HTOP = ALPHATOP*ZM
      HCUT = ABS(HTOP)/PXZCUT

      print*,'TRACK_INIT hcut= ',hcut
      print*,'TRACK_INIT eff = ',eff
      print*,'TRACK_INIT b0 = ',b0
      print*,'TRACK_INIT bl3 = ',bl3
      print*,'TRACK_INIT sigmacut = ',sigcut
      print*,'TRACK_INIT cutpxz = ',pxzcut
      print*,'TRACK_INIT xprec = ',xprec
      print*,'TRACK_INIT yprec = ',yprec

      EFF2 = EFF**2               ! PROBA. DEUX CHAMBRES TOUCHES
      EFF1 = EFF2+2.*EFF*(1.-EFF) ! PROBA. AU MOINS UNE CHAMBRE TOUCHE
** Used only for stations 4 & 5
      PHIPREC   = SQRT(2.)*XPREC/DZ_PL(5) ! PHI = (OZ , PROJ. DANS XOZ)
      ALAMPREC  = SQRT(2.)*YPREC/DZ_PL(5) ! LAM = (OM , PROJ. DANS XOZ)

      DO I = 1,5
         NRES(I) = 0
      ENDDO
      NRESF = 0 
      NTRMUALL = 0 
      NMUONALL = 0
      NGHOSTALL = 0
      NTRACKFALL = 0
      DO I = 1,NBSTATION
         NERRALL(I) = 0
      ENDDO   
      IR = 0 
*
      RETURN
      END

*************************************************************************
      SUBROUTINE PREC_INIT
*************************************************************************
*  
*
*************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)      
*            
      PARAMETER(NPLANE=10,NBSTATION=5) 
*           
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
*      
      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
*           
      COMMON/PRECCUT/PCUT,PTCUT,CHI2CUT
*
      COMMON/PRECSUM/NRESF1,NMUONALL1,NGHOSTALL1,NTRACKFALL1 
*          
      COMMON/DEBEVT/IDEBUG
*

      DATA THICK/0.03D0/   ! X/X0=3% 
**      DATA THICK/0.02D0/   ! X/X0=2% chambre

**      DATA B0,BL3/10.,2.0/ ! Champ magnetique dans le dipole et dans L3 en kgauss
      DATA B0,BL3/7.,2.0/ ! Magnetic field in the dipole & L3 in kgauss
      DATA ZMAGS/805.0D0/,ZMAGE/1233.0D0/,ZABS/503.D0/ ! CCC not used when
                                                       ! magn. field map 
*      DATA ZMAGS/825.0D0/,ZMAGE/1125.0D0/,ZABS/503.D0/
      DATA XMAG/190.0D0/

**      DATA XPREC/0.0100D0/,YPREC/0.144337D0/ ! CCC
      DATA XPREC/0.0100D0/,YPREC/0.2D0/

* Input parameters       
      CONST = 0.299792458D-3*B0*(ZMAGE-ZMAGS)
      J = 0
      DO I = 1,5
         J = J+1
         ZPLANEP(J) = -ZPLANE(I)
         J = J+1
         ZPLANEP(J) = -ZPLANE(I)+DZ_PL(I)  
      ENDDO 

      do i=1,10
         print*,'zplanep(',i,')=',ZPLANEP(I)
      end do

      PCUT = 3.        ! Coupure en PXZ muon (GeV/c)
      PTCUT = 0.       ! Coupure en Pt muon (GeV/c)
 
      CHI2CUT = 1.E4  ! Coupure sur le CHI2 du fit

      X01 = 18.8     ! C (cm)
      X02 = 10.397   ! Concrete (cm)
      X03 = 0.56     ! Plomb (cm)
      X04 = 47.26    ! Polyethylene (cm)
      X05 = 0.35     ! W (cm) 

** Calcul des parametres pour la correction de Branson de l'absorbeur     
      ANBP = (315.**3-90.**3)/X01 +(467.**3-315.**3)/X02+
     &       (472.**3-467.**3)/X03+(477.**3-472.**3)/X04+
     &       (482.**3-477.**3)/X03+(487.**3-482.**3)/X04+
     &       (492.**3-487.**3)/X03+(497.**3-492.**3)/X04+
     &       (502.**3-497.**3)/X03
      ADBP = (315.**2-90.**2)/X01 +(467.**2-315.**2)/X02+
     &       (472.**2-467.**2)/X03+(477.**2-472.**2)/X04+
     &       (482.**2-477.**2)/X03+(487.**2-482.**2)/X04+
     &       (492.**2-487.**2)/X03+(497.**2-492.**2)/X04+
     &       (502.**2-497.**2)/X03
      ZBP1 = 2./3.*ANBP/ADBP
      ANBP = (315.**3-90.**3)/X01 +(467.**3-315.**3)/X02+
     &       (503.**3-467.**3)/X05
      ADBP = (315.**2-90.**2)/X01 +(467.**2-315.**2)/X02+
     &       (503.**2-467.**2)/X05
      ZBP2 = 2./3.*ANBP/ADBP
*      
      IF (IDEBUG.GE.1) THEN
         PRINT *,' PREC_INIT B0 (kgauss)',B0,' BL3 (kgauss)',BL3
         PRINT *,' PREC_INIT ZMAGE (cm)',ZMAGE,' ZMAGS (cm)',ZMAGS,
     &        ' XMAG (cm)',XMAG
         PRINT *,' PREC_INIT ZABS (cm)',ZABS,' ZBP1 (cm)',ZBP1,
     &        ' ZBP2 (cm)',ZBP2
         PRINT *,' PREC_INIT Radiation length absorber X01 (cm)',X01,
     &        ' X02 (cm)',X02
         PRINT *,' PREC_INIT X03(cm)',X03,' X04 (cm)',X04,' X05 (cm)',
     &           X05
         PRINT *,' PREC_INIT Radiation length chamber THICK (%)',
     &           THICK*100.
         PRINT *,' PREC_INIT XPREC (cm)',XPREC,' YPREC (cm)',YPREC
         PRINT *,' PREC_INIT Coupure en Pxz (GeV/c): ',PCUT
         PRINT *,' PREC_INIT Coupure en Pt (GeV/c): ',PTCUT
         PRINT *,' PREC_INIT Coupure en CHI2 : ',CHI2CUT
      ENDIF

*      PAUSE

      NRESF1 = 0 
      NMUONALL1 = 0
      NGHOSTALL1 = 0
      NTRACKFALL1 = 0
           
*  Magnetic Field Map GC
**      OPEN (UNIT=40,FILE='data/field02.dat',
**     &      STATUS='UNKNOWN')

      CALL INITFIELD
                  
**      CLOSE(40)

      RETURN
      END
      
*************************************************************************
      SUBROUTINE RECO_READEVT(NEV,IDRES,IREADGEANT)
*************************************************************************
*
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      PARAMETER(NTRMAX=500) 

      PARAMETER (NBSTATION=5,MAXIDG=20000,MAXHITTOT=20000,
     &           MAXHITCH=10000,MAXHIT=1000,NBCHAMBER=10)

      COMMON/RHITG/ITYPG(MAXIDG),XTRG(MAXIDG),YTRG(MAXIDG),
     &             PTOTG(MAXIDG),IDG(MAXIDG),IZCH(MAXIDG),
     &             PVERT1G(MAXIDG),PVERT2G(MAXIDG),PVERT3G(MAXIDG),
     &     ZVERTG(MAXIDG),NHITTOT1,CX(MAXIDG),CY(MAXIDG),CZ(MAXIDG),
     &     XGEANT(MAXIDG),YGEANT(MAXIDG),CLSIZE1(MAXIDG),CLSIZE2(MAXIDG)

      DIMENSION TYPG(MAXIDG),ZCH(MAXIDG)

      REAL*4 R1,R2
      DATA R1,R2/0.,1./

      IF (IREADGEANT.eq.1) THEN  ! GEANT hits

         CALL TRACKF_READ_GEANT(ITYPG,XTRG,YTRG,PTOTG,IDG,IZCH,PVERT1G,
     &        PVERT2G,PVERT3G,ZVERTG,NHITTOT1,CX,CY,CZ,IEVR,NEV,
     &        XGEANT,YGEANT,CLSIZE1,CLSIZE2)
      ELSE ! reconstructed hits
         CALL TRACKF_READ_SPOINT(ITYPG,XTRG,YTRG,PTOTG,IDG,IZCH,PVERT1G,
     &        PVERT2G,PVERT3G,ZVERTG,NHITTOT1,CX,CY,CZ,IEVR,NEV,
     &        XGEANT,YGEANT,CLSIZE1,CLSIZE2)
      ENDIF
      
      do i=1,NHITTOT1
         TYPG(i)=ITYPG(i)
         call chfill(100,sngl(typg(i)),R1,R2)
         call chfill(101,sngl(ygeant(i)),R1,R2)
         call chfill(102,sngl(xgeant(i)),R1,R2)
         ZCH(i)=IZCH(i)
         call chfill(103,sngl(zch(i)),R1,R2)
         call chfill(104,sngl(ptotg(i)),R1,R2)
         call chfill(105,sngl(pvert2g(i)),R1,R2)
         call chfill(106,sngl(pvert1g(i)),R1,R2)
         call chfill(107,sngl(pvert3g(i)),R1,R2)
         call chfill(108,sngl(zvertg(i)),R1,R2)
         call chfill(109,sngl(ytrg(i)),R1,R2)
         call chfill(110,sngl(xtrg(i)),R1,R2)
      enddo

      do i=1,NHITTOT1
         CALL CHFILL (999,SNGL(PTOTG(I)),R1,R2)
      enddo

      CALL TRACKF_STAT(IDRES,IREADGEANT)

      RETURN
      END

*************************************************************************
      SUBROUTINE TRACKF_STAT(IDRES,IREADGEANT)
*************************************************************************
*  Associate hits between two chambers inside a station
*  Simulate spatial resolution and chamber efficiency
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (NBSTATION=5,MAXIDG=20000,MAXHITTOT=20000,
     &           MAXHITCH=10000,MAXHIT=1000,NBCHAMBER=10)
*      
      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP
*
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
*
* HITS GEANT initiaux par chambre
      COMMON/RHITG/ITYPG(MAXIDG),XTRG(MAXIDG),YTRG(MAXIDG),
     &             PTOTG(MAXIDG),IDG(MAXIDG),IZCH(MAXIDG),
     &             PVERT1G(MAXIDG),PVERT2G(MAXIDG),PVERT3G(MAXIDG),
     &     ZVERTG(MAXIDG),NHITTOT1,CX(MAXIDG),CY(MAXIDG),CZ(MAXIDG),
     &     XGEANT(MAXIDG),YGEANT(MAXIDG),CLSIZE1(MAXIDG),CLSIZE2(MAXIDG)

* HITS GEANT associes par station
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

* 
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
*
      COMMON/DEBEVT/IDEBUG
      common/dstation/idstation
*
      DIMENSION RMIN(NBCHAMBER),RMAX1(NBCHAMBER)
      DIMENSION XMA(NBCHAMBER,MAXHITCH),YMA(NBCHAMBER,MAXHITCH),
     &          IMARK(NBCHAMBER,MAXHITCH)
* 
      DIMENSION IEFFI(MAXHITTOT)
      DIMENSION IH(NBCHAMBER,MAXHIT) 
      DIMENSION NHIT(NBCHAMBER)
      DIMENSION DXMAX(NBSTATION),DYMAX(NBSTATION),VRES(2,5)
  
      REAL*4 RNDM,RN,RN1,RN2,R1,R2

* Chambre 10 deg.      
      DATA RMAX1/91.5,91.5,122.5,122.5,158.3,158.3,260.,260.,260.,260./
* Zone de recherche entre deux plans d'une station 

*      if (idstation.eq.8) then
*         DATA DXMAX/1.,1.,1.2,2.4,2.4/ ! dz_ch = 8 cm
*      else if (idstation.eq.20) then
         DATA DXMAX/1.5,1.5,3.,6.,6./ ! dz_ch = 20cm
*      end if

      DATA DYMAX/0.22,0.22,0.22,0.22,0.22/  ! CCC Upsilon dz_ch = 20 cm

      DATA R1,R2/0.,1./
       
       
      ICH = 0
      DO IZ=1,5
         ICH = ICH+1
         RMIN(ICH) =  ABS(ZPLANE(IZ)*TAN(2.*ACOS(-1.)/180))
         IF (IZ.GT.2) RMIN(ICH) = 30.
         ICH = ICH+1
         RMIN(ICH) =  ABS(ZPLANE(IZ)*TAN(2.*ACOS(-1.)/180))
         IF (IZ.GT.2) RMIN(ICH) = 30.
      ENDDO   

*  Initialisations 
      DO ICH = 1,10
         NHIT(ICH) = 0
      ENDDO
         
*  1 ere boucle de lecture des hits initiaux 

      
      IF (IREADGEANT.EQ.1) THEN
         DO I = 1,2
            DO J = 1,5
               VRES(I,J) = 0.
            ENDDO
         ENDDO    
         IMU = 0 
         DO I = 1,NHITTOT1
            ICH = IZCH(I)
            IZ = INT(FLOAT(ICH+1)/2.)
            IMOD = MOD(ICH,2)

            IF (IMOD.NE.0.AND.IZ.LE.5) THEN
               CALL CHFILL2(1000+IZ,SNGL(XGEANT(I)),SNGL(YGEANT(I)),R2)
            ENDIF   
            IF (ICH.EQ.9) THEN
               ISTAK = IDG(I)
               ISTAK = MOD(ISTAK,30000)
               ISTAK = MOD(ISTAK,10000)

               IF ((ITYPG(I).EQ.5.OR.ITYPG(I).EQ.6).AND.
     &              ISTAK.EQ.IDRES) THEN ! upsilon
                  IMU = IMU+1
                  VRES(1,IMU) = XGEANT(I)
                  VRES(2,IMU) = YGEANT(I)
               ENDIF
            ENDIF   
         ENDDO
      ENDIF   


      DO I = 1,NHITTOT1 ! Boucle sur les hits GEANT de toutes les ch.
      
**        IF (ITYPG(I).NE.5.AND.ITYPG(I).NE.6) GOTO 1 ! CCC

         ICH = IZCH(I)
         
         IF (ICH.GT.10) GO TO 1
         
         IF (IREADGEANT.EQ.1) THEN ! GEANT hits

            IF (ICH.EQ.9.OR.ICH.EQ.10) THEN
               DNUM = 999.
               DO IM = 1,IMU
                  DNU = SQRT((XGEANT(I)-VRES(1,IM))**2+
     &                 (YGEANT(I)-VRES(2,IM))**2)
                  IF (DNU.LT.DNUM) DNUM = DNU
               ENDDO
               IF (DNUM.GT.50.) GO TO 1 ! discard hits far from MUONS
            ENDIF

            CALL RANNOR(RN1,RN2)          ! CCC
            X = XGEANT(I)
            Y = YGEANT(I)
*            X = XGEANT(I) + RN1 * XPREC
*            Y = YGEANT(I) + RN2 * YPREC
*     efficacite des chambres        
            IEFFI(I) = 1
            RN = RNDM()              
            IF (RN.GT.EFF) IEFFI(I) = 0
**            IF (ICH.EQ.9.OR.ICH.EQ.10) THEN 
**               PRINT *,' HIT GEANT',' ICH=',ICH,' I=',I,' X =',X,' Y=',
**     &              Y,' IDG=',IDG(I)
**            ENDIF   

            IF (ITYPG(I).EQ.5.OR.ITYPG(I).EQ.6) THEN
               ISTAK = IDG(I)
               ISTAK = MOD(ISTAK,30000)
               ISTAK = MOD(ISTAK,10000)

               IF (ISTAK.EQ.IDRES) then
                  dx=xgeant(i)-x
                  dy=ygeant(i)-y
                  IZ = INT(FLOAT(IZCH(i)+1)/2.)
                  ichx=110+IZ
                  ichy=120+IZ
                  call chfill(ichx,sngl(dy),R1,R2)
                  call chfill(ichy,sngl(dx),R1,R2)
                  if (iz.eq.1) call chfill(116,sngl(dy),R1,R2)	
                  if (iz.eq.2) call chfill(117,sngl(dy),R1,R2)	
               end if
            end if

         ELSE  ! reconstructed hits 

            IEFFI(I) = 1

            X = XTRG(I)
            Y = YTRG(I)

*     étude des hits geant avec un seul fichier
*            CALL RANNOR(RN1,RN2)          ! CCCC
*            X = XGEANT(I) + RN1 * XPREC
*            Y = YGEANT(I) + RN2 * YPREC

            IF (ITYPG(I).EQ.5.OR.ITYPG(I).EQ.6) THEN
               ISTAK = IDG(I)
               ISTAK = MOD(ISTAK,30000)
               ISTAK = MOD(ISTAK,10000)

               IF (ISTAK.EQ.IDRES) then
                  dx=xgeant(i)-x
                  dy=ygeant(i)-y
                  IZ = INT(FLOAT(IZCH(i)+1)/2.)
                  ichx=110+IZ
                  ichy=120+IZ
                  call chfill(ichx,sngl(dy),R1,R2)
                  call chfill(ichy,sngl(dx),R1,R2)
                  if (iz.eq.1) call chfill(116,sngl(dy),R1,R2)	
                  if (iz.eq.2) call chfill(117,sngl(dy),R1,R2)	
               end if
            end if

         ENDIF

         R = SQRT(X**2+Y**2)
**         IF (R.LT.RMIN(ICH).OR.R.GT.RMAX1(ICH)) then
**            if (ich.le.10) then
**               print*,'* chambre ',ich,' * hit ',i
**               print*,'ityp=',itypg(i),' x=',X,' y=',Y
**               print*,'R=',R,' RMIN=',RMIN(ICH),' RMAX1=',RMAX1(ICH)
**            endif
**            GO TO 1        ! CCC   
**         end if

 	 NHIT(ICH) = NHIT(ICH)+1
	 IH(ICH,NHIT(ICH)) = I
	 XMA(ICH,NHIT(ICH)) = X
	 YMA(ICH,NHIT(ICH)) = Y
         IMARK(ICH,NHIT(ICH)) = 0

  1      CONTINUE
      ENDDO

* Association des hits entre chambres d'une station
      II = 0            ! nombre de hits GEANT par station
      DO ICH1 = 1,10,2  ! loop on chamber
         IZ = INT(FLOAT(ICH1+1)/2.)
         JHIT(IZ) = 0
         ICH2 = ICH1+1

         DO I1 = 1,NHIT(ICH1) ! loop on hits in 1st chamber
            II = II+1
            IFIND = 0 
            I = IH(ICH1,I1)

            ITYP(II) = ITYPG(I) 
            XTR(II) = XTRG(I) 
            YTR(II) = YTRG(I)
            PTOT(II) = PTOTG(I)
            ID(II) = IDG(I)
            IZST(II) = IZ
            PVERT1(II) = PVERT1G(I)  
            PVERT2(II) = PVERT2G(I)  
            PVERT3(II) = PVERT3G(I)
            ZVERT(II) = ZVERTG(I)
 
            IF (IEFFI(I).EQ.1) THEN
               X1 = XMA(ICH1,I1)
               Y1 = YMA(ICH1,I1)
               ID1 = IDG(I)
               XEXT1 = (ZPLANE(IZ)-DZ_PL(IZ))/ZPLANE(IZ)*X1
               YEXT1 = (ZPLANE(IZ)-DZ_PL(IZ))/ZPLANE(IZ)*Y1
    
               DO I2 = 1,NHIT(ICH2)  ! loop on hits in 2nd chamber
                  J = IH(ICH2,I2)

                  IF (IEFFI(J).EQ.1) THEN
                     X2 = XMA(ICH2,I2)
                     Y2 = YMA(ICH2,I2)
                     ID2 = IDG(J)
                     DX = X2-XEXT1
                     DY = Y2-YEXT1

                     IF (ID1.EQ.ID2.AND.
     &                    (ITYP(II).EQ.5.OR.ITYP(II).EQ.6)) THEN
                        CALL CHFILL(70+IZ,SNGL(DX),R1,R2)
                        CALL CHFILL(80+IZ,SNGL(DY),R1,R2)
                     ENDIF   
                     DX = ABS(DX) 
                     DY = ABS(DY)

                     IF (DX.LT.DXMAX(IZ).AND.DY.LT.(SIGCUT*DYMAX(IZ)) ! CCC
     &                   ) THEN
                        IFIND = 1
                        IMARK(ICH2,I2) = 1
                        JHIT(IZ) =  JHIT(IZ)+1
                        XM(IZ,JHIT(IZ)) = X1
                        YM(IZ,JHIT(IZ)) = Y1
                        IZM(IZ,JHIT(IZ)) = 1
                        PHM(IZ,JHIT(IZ)) = -ATAN((X2-X1)/DZ_PL(IZ))
                        ALM(IZ,JHIT(IZ)) = ATAN((Y2-Y1)/DZ_PL(IZ)*
     &                                     COS(PHM(IZ,JHIT(IZ)))) 
                        XMR(IZ,JHIT(IZ),1) = X1
                        YMR(IZ,JHIT(IZ),1) = Y1
                        XMR(IZ,JHIT(IZ),2) = X2
                        YMR(IZ,JHIT(IZ),2) = Y2
                        IP(IZ,JHIT(IZ)) = II
                  
                        ISTAK = ID2
                        ISTAK = MOD(ISTAK,30000)
                        ISTAK = MOD(ISTAK,10000)

                        IF ((ITYPG(J).EQ.5.OR.ITYPG(J).EQ.6).AND.
     &                       ISTAK.EQ.IDRES) THEN ! upsilon or J/psi
                          
                           ITYP(II) = ITYPG(J) 
                           XTR(II) = XTRG(J) 
                           YTR(II) = YTRG(J)
                           PTOT(II) = PTOTG(J)
                           ID(II) = IDG(J)
                           PVERT1(II) = PVERT1G(J)  
                           PVERT2(II) = PVERT2G(J)  
                           PVERT3(II) = PVERT3G(J)
                           ZVERT(II) = ZVERTG(J)

                        ENDIF   

                     ENDIF
                  ENDIF 
               ENDDO  ! loop on hits in 2nd chamber

               IF (IFIND.EQ.0) THEN ! No possible association
                  JHIT(IZ) =  JHIT(IZ)+1
                  XM(IZ,JHIT(IZ)) = X1
                  YM(IZ,JHIT(IZ)) = Y1
                  IZM(IZ,JHIT(IZ)) = 1
                  IP(IZ,JHIT(IZ)) = II
                  PHM(IZ,JHIT(IZ)) = 10.
                  ALM(IZ,JHIT(IZ)) = 10.
                  XMR(IZ,JHIT(IZ),1) = X1
                  YMR(IZ,JHIT(IZ),1) = Y1
                  XMR(IZ,JHIT(IZ),2) = 0.
                  YMR(IZ,JHIT(IZ),2) = 0.
               ENDIF
            ENDIF
         ENDDO  ! end loop on hits in 1st chamber
      ENDDO  ! end loop on chamber 


* On conserve les HITS (x,y) de la 2nde chambre des stations      

      DO ICH = 2,10,2 ! Loop on 2nd chambers
         IZ = INT(FLOAT(ICH+1)/2.)
         DO I = 1,NHIT(ICH) ! Loop on hits
            J = IH(ICH,I)

            IF (IMARK(ICH,I).EQ.0) THEN ! hit not already associated
  
               II = II+1
            
               ITYP(II) = ITYPG(J) 
               XTR(II) = XTRG(J) 
               YTR(II) = YTRG(J)
               PTOT(II) = PTOTG(J)
               ID(II) = IDG(J)
               IZST(II) = IZ
               PVERT1(II) = PVERT1G(J)  
               PVERT2(II) = PVERT2G(J)  
               PVERT3(II) = PVERT3G(J)
               ZVERT(II) = ZVERTG(I)  
                
               IF (IEFFI(J).EQ.1) THEN
                  JHIT(IZ) =  JHIT(IZ)+1
                  XM(IZ,JHIT(IZ)) = XMA(ICH,I)
                  YM(IZ,JHIT(IZ)) = YMA(ICH,I) 
                  IZM(IZ,JHIT(IZ)) = 2
                  PHM(IZ,JHIT(IZ)) = 10.
                  ALM(IZ,JHIT(IZ)) = 10.
                  IP(IZ,JHIT(IZ)) = II  
                  XMR(IZ,JHIT(IZ),1) = 1000.
                  YMR(IZ,JHIT(IZ),1) = 1000.   
                  XMR(IZ,JHIT(IZ),2) = XMA(ICH,I)
                  YMR(IZ,JHIT(IZ),2) = YMA(ICH,I)                    
               ENDIF 
            ENDIF
         ENDDO ! End loop on hits
      ENDDO  ! End loop on 2nd chambers


      NHITTOT = II ! total number of hits in stations
*      

      IF (IDEBUG.GE.2) THEN
         PRINT *,'TRACKF_STAT nb hits:',NHITTOT 
      ENDIF     

**      DO IZ = 5,5
**         PRINT *,' IZ=',IZ,' JHIT(IZ)=',JHIT(IZ)
**         DO J = 1,JHIT(IZ)
**            II = IP(IZ,J)
**            PRINT *,' HIT ASS.',' IZ=',IZ,' II=',II,' X =',
**     &           XM(IZ,J),' Y=',YM(IZ,J),' ID=',ID(II)
**            PRINT *,' HIT ASS.',' IZ=',IZ,' II=',II,' X1 =',
**     &           XMR(IZ,J,1),' Y1=',YMR(IZ,J,1),' ID=',ID(II)
**            PRINT *,' HIT ASS.',' IZ=',IZ,' II=',II,' X2 =',
**     &           XMR(IZ,J,2),' Y2=',YMR(IZ,J,2),' ID=',ID(II)
**         ENDDO
**      ENDDO 

**    DO IZ = 1,NBSTATION
**       PRINT *,' IZ=',IZ,' JHIT(IZ)=',JHIT(IZ)
**       DO J = 1,JHIT(IZ)
**          PRINT *,' PHM(IZ,J)=',PHM(IZ,J),' ALM(IZ,J)=',ALM(IZ,J)
**       ENDDO
**     ENDDO 
      
      RETURN
      END
*************************************************************************
      SUBROUTINE TRACKF_STAT_NEW(IDRES)
*************************************************************************
*  Associate hits between two chambers inside a station
*  Simulate spatial resolution and chamber efficiency
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (NBSTATION=5,MAXIDG=20000,MAXHITTOT=20000,
     &           MAXHITCH=10000,MAXHIT=1000,NBCHAMBER=10)
*      
      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP
*
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
*
* HITS GEANT initiaux par chambre
      COMMON/RHITG/ITYPG(MAXIDG),XTRG(MAXIDG),YTRG(MAXIDG),
     &             PTOTG(MAXIDG),IDG(MAXIDG),IZCH(MAXIDG),
     &             PVERT1G(MAXIDG),PVERT2G(MAXIDG),PVERT3G(MAXIDG),
     &     ZVERTG(MAXIDG),NHITTOT1,CX(MAXIDG),CY(MAXIDG),CZ(MAXIDG),
     &     XGEANT(MAXIDG),YGEANT(MAXIDG),CLSIZE1(MAXIDG),CLSIZE2(MAXIDG)

* HITS GEANT associes par station
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

* 
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
*
      COMMON/DEBEVT/IDEBUG
*
      DIMENSION RMIN(NBCHAMBER),RMAX1(NBCHAMBER)
      DIMENSION XMA(NBCHAMBER,MAXHITCH),YMA(NBCHAMBER,MAXHITCH),
     &          IMARK(NBCHAMBER,MAXHITCH)
* 
      DIMENSION IEFFI(MAXHITTOT)
      DIMENSION IH(NBCHAMBER,MAXHIT) 
      DIMENSION NHIT(NBCHAMBER)
      DIMENSION DXMAX(NBSTATION),DYMAX(NBSTATION),I2C(1000)
  
      DIMENSION DIST(2),NMUON(2),NHITMUON(2,5),NMUONGOOD(2)

      REAL*4 RNDM,RN,RN1,RN2,R1,R2

* Chambre 10 deg.      
      DATA RMAX1/91.5,91.5,122.5,122.5,158.3,158.3,260.,260.,260.,260./
* Zone de recherche entre deux plans d'une station 
**      DATA DXMAX/2.,1.5,2.5,3.,3./
      DATA DXMAX/1.5,1.5,3.,3.,3./
      DATA DYMAX/0.22,0.22,0.21,0.21,0.21/
      DATA R1,R2/0.,1./

      ICH = 0
      DO IZ=1,5
         ICH = ICH+1
         RMIN(ICH) =  ABS(ZPLANE(IZ)*TAN(2.*ACOS(-1.)/180))
         IF (IZ.GT.2) RMIN(ICH) = 30.
         ICH = ICH+1
         RMIN(ICH) =  ABS(ZPLANE(IZ)*TAN(2.*ACOS(-1.)/180))
         IF (IZ.GT.2) RMIN(ICH) = 30.
      ENDDO   

*    Initialisations 
      DO ICH = 1,10
         NHIT(ICH) = 0
      ENDDO
         
      DO NCH = 1,10
         DO J=1,2
            DIST(J)=999.
            NMUON(J)=0
         ENDDO
         DO I = 1,NHITTOT1
            IF (IZCH(I).EQ.NCH) THEN
               ISTAK = IDG(I)
               ISTAK = MOD(ISTAK,30000)
               ISTAK = MOD(ISTAK,10000)
               IF (ISTAK.EQ.IDRES.AND.IDG(I).EQ.50116) THEN
                  DISTMIN=(XTRG(I)-XGEANT(I))**2+(YTRG(I)-YGEANT(I))**2
                  IF (DISTMIN.LT.DIST(1)) THEN
                     DIST(1)=DISTMIN
                     NMUONGOOD(1)=I
                  ENDIF
                  NMUON(1)=NMUON(1)+1
                  NHITMUON(1,NMUON(1))=I
               ELSE 
                IF (ISTAK.EQ.IDRES.AND.IDG(I).EQ.70116) THEN
                   DISTMIN=(XTRG(I)-XGEANT(I))**2+(YTRG(I)-YGEANT(I))**2
                   IF (DISTMIN.LT.DIST(2)) THEN
                        DIST(2)=DISTMIN
                        NMUONGOOD(2)=I
                   ENDIF
                   NMUON(2)=NMUON(2)+1
                   NHITMUON(2,NMUON(2))=I
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         DO J=1,2
            IF (NMUON(J).GE.2) THEN
             print*,'j=',j,' nmuon=',nmuon(j)
             print*,'chambre',nch
             DO K=1,NMUON(J)
               IF (NHITMUON(J,K).NE.NMUONGOOD(J)) IDG(NHITMUON(J,K))=999 ! flag les mauvais hits MUONS
             ENDDO
            ENDIF
         ENDDO
      ENDDO
      

* 1 ere boucle Lecture des hits initiaux 

      DO I = 1,NHITTOT1 ! Boucle sur les hits GEANT de toutes les ch.
      
         ICH = IZCH(I)
         
         X = XTRG(I)
         Y = YTRG(I)

         R = SQRT(X**2+Y**2)
         IF (R.LT.RMIN(ICH).OR.R.GT.RMAX1(ICH)) then
            if (ich.le.10.and.i.le.28) then
               print*,'****** chambre ',ich,' ****** hit ',i
               print*,'ityp=',itypg(i)
               print*,'x=',XTRG(I),' y=',YTRG(I)
               print*,'R=',R,' RMIN=',RMIN(ICH),' RMAX1=',RMAX1(ICH)
           endif
            GO TO 1   
         end if

         IEFFI(I) = 1

 	 NHIT(ICH) = NHIT(ICH)+1
	 IH(ICH,NHIT(ICH)) = I
	 XMA(ICH,NHIT(ICH)) = XTRG(I)
	 YMA(ICH,NHIT(ICH)) = YTRG(I)
         IMARK(ICH,NHIT(ICH)) = 0

         print*,' XTRG(I)=', XTRG(I),' YTRG(I)=', YTRG(I),' IDG(I)=',
     &        IDG(I),' ICH=',ICH
	 
  1      CONTINUE
      ENDDO
      
      
* Association des hits entre chambres d'une station
      II = 0            ! nombre de hits GEANT par station
      DO ICH1 = 1,10,2
         IZ = INT(FLOAT(ICH1+1)/2.)
         JHIT(IZ) = 0
         ICH2 = ICH1+1

         DO I1 = 1,NHIT(ICH1)
            II = II+1
            IFIND = 0 
            I = IH(ICH1,I1)

            ITYP(II) = ITYPG(I) 
            XTR(II) = XTRG(I) 
            YTR(II) = YTRG(I)
            PTOT(II) = PTOTG(I)
            ID(II) = IDG(I)
            IZST(II) = IZ
            PVERT1(II) = PVERT1G(I)  
            PVERT2(II) = PVERT2G(I)  
            PVERT3(II) = PVERT3G(I)
            ZVERT(II) = ZVERTG(I)
 
            IF (IEFFI(I).EQ.1) THEN
               X1 = XMA(ICH1,I1)
               Y1 = YMA(ICH1,I1)
               ID1 = IDG(I)
               XEXT1 = (ZPLANE(IZ)-DZ_PL(IZ))/ZPLANE(IZ)*X1
               YEXT1 = (ZPLANE(IZ)-DZ_PL(IZ))/ZPLANE(IZ)*Y1    
               KC = 0
               PRINT *,'***** DEBUT RECHERCHE',' ID1=',ID1,' ich1=',ICH1
               PRINT *,' XTR(II)=', XTR(II),' YTR(II)=', YTR(II) 
               PRINT *,'  ITYP(II)=', ITYP(II)
               DO I2 = 1,NHIT(ICH2)
                  J = IH(ICH2,I2)
                  IF (IEFFI(J).EQ.1) THEN
                     X2 = XMA(ICH2,I2)
                     DX = X2-XEXT1
                     DX = ABS(DX) 
                     IF (DX.LT.DXMAX(IZ)) THEN
                        KC = KC + 1
                        I2C(KC) = I2
                        ID2 = IDG(J)
                        print *,' DX=',DX,' KC=',KC,' ID2=',ID2 
                     ENDIF
                  ENDIF
               ENDDO   
               DYOLD = 999.
               I2FIND = 0
               DO IKC = 1,KC
                  I2 = I2C(IKC)
                  Y2 = YMA(ICH2,I2)
                  DY = Y2-YEXT1
                  DY = ABS(DY)
                  J = IH(ICH2,I2)
                  ID2 = IDG(J)
                  IF (DY.LT.DYOLD.AND.DY.LT.(SIGCUT*DYMAX(IZ))) THEN
                     DYOLD = DY
                     I2FIND = I2
                     PRINT *,' ID2=',ID2,' DY=',DY
                  ENDIF
               ENDDO   
               IF (I2FIND.GT.0) THEN
                  I2 = I2FIND
                  J = IH(ICH2,I2)
                  ID2 = IDG(J)
                  IFIND = 1
                  IMARK(ICH2,I2) = 1
                  JHIT(IZ) =  JHIT(IZ)+1
                  X2 = XMA(ICH2,I2)
                  Y2 = YMA(ICH2,I2)
                  XM(IZ,JHIT(IZ)) = X1
                  YM(IZ,JHIT(IZ)) = Y1
                  IZM(IZ,JHIT(IZ)) = 1
                  PHM(IZ,JHIT(IZ)) = -ATAN((X2-X1)/DZ_PL(IZ))
                  ALM(IZ,JHIT(IZ)) = ATAN((Y2-Y1)/DZ_PL(IZ)*
     &                 COS(PHM(IZ,JHIT(IZ)))) 
                  XMR(IZ,JHIT(IZ),1) = X1
                  YMR(IZ,JHIT(IZ),1) = Y1
                  XMR(IZ,JHIT(IZ),2) = X2
                  YMR(IZ,JHIT(IZ),2) = Y2
                  IP(IZ,JHIT(IZ)) = II
                  
                  ISTAK = ID2
                  ISTAK = MOD(ISTAK,30000)
                  ISTAK = MOD(ISTAK,10000)
* test
                  IF (ISTAK.EQ.IDRES.AND.ID1.NE.999) THEN
                          
                     ITYP(II) = ITYPG(J) 
                     PTOT(II) = PTOTG(J)
                     XTR(II) = XTRG(I) 
                     YTR(II) = YTRG(I)
                     ID(II) = IDG(J)
                     PVERT1(II) = PVERT1G(J)  
                     PVERT2(II) = PVERT2G(J)  
                     PVERT3(II) = PVERT3G(J)
                     ZVERT(II) = ZVERTG(J)

                  ENDIF  
               ENDIF   
  
                           
          
               IF (IFIND.EQ.0) THEN
                  JHIT(IZ) =  JHIT(IZ)+1
                  XM(IZ,JHIT(IZ)) = X1
                  YM(IZ,JHIT(IZ)) = Y1
                  IZM(IZ,JHIT(IZ)) = 1
                  IP(IZ,JHIT(IZ)) = II
                  PHM(IZ,JHIT(IZ)) = 10.
                  ALM(IZ,JHIT(IZ)) = 10.
                  XMR(IZ,JHIT(IZ),1) = X1
                  YMR(IZ,JHIT(IZ),1) = Y1
                  XMR(IZ,JHIT(IZ),2) = 0.
                  YMR(IZ,JHIT(IZ),2) = 0.
               ENDIF
               CALL CHFILL2(1000+IZ,SNGL(X1),SNGL(Y1),R2)
            ENDIF
         ENDDO
      ENDDO  

* On conserve les HITS de la 2nde chambre des stations      

      DO ICH = 2,10,2 
         IZ = INT(FLOAT(ICH+1)/2.)
         DO I = 1,NHIT(ICH)
            J = IH(ICH,I)

            IF (IMARK(ICH,I).EQ.0) THEN
**               print *,' ich=',ich,' i=',i,' j=',j
  
               II = II+1
            
               ITYP(II) = ITYPG(J) 
               XTR(II) = XTRG(J) 
               YTR(II) = YTRG(J)
               PTOT(II) = PTOTG(J)
               ID(II) = IDG(J)
               IZST(II) = IZ
               PVERT1(II) = PVERT1G(J)  
               PVERT2(II) = PVERT2G(J)  
               PVERT3(II) = PVERT3G(J)
               ZVERT(II) = ZVERTG(I)  
                
               IF (IEFFI(J).EQ.1) THEN
                  JHIT(IZ) =  JHIT(IZ)+1
**                  XM(IZ,JHIT(IZ)) = ZPLANE(IZ)/(ZPLANE(IZ)-DZ_PL)
**     &                              *XMA(ICH,I)
**                  YM(IZ,JHIT(IZ)) = ZPLANE(IZ)/(ZPLANE(IZ)-DZ_PL)
**     &                              *YMA(ICH,I) 
                  XM(IZ,JHIT(IZ)) = XMA(ICH,I)
                  YM(IZ,JHIT(IZ)) = YMA(ICH,I) 
                  IZM(IZ,JHIT(IZ)) = 2
                  PHM(IZ,JHIT(IZ)) = 10.
                  ALM(IZ,JHIT(IZ)) = 10.
                  IP(IZ,JHIT(IZ)) = II  
                  XMR(IZ,JHIT(IZ),1) = 1000.
                  YMR(IZ,JHIT(IZ),1) = 1000.   
                  XMR(IZ,JHIT(IZ),2) = XMA(ICH,I)
                  YMR(IZ,JHIT(IZ),2) = YMA(ICH,I)                    
               ENDIF 
            ENDIF
         ENDDO
      ENDDO


      NHITTOT = II
*      
      IF (IDEBUG.GE.2) THEN
         PRINT *,'TRACKF_MICRO nb hits:',NHITTOT 
      ENDIF     
      DO IZ = 1,NBSTATION
         PRINT *,' IZ=',IZ,' JHIT(IZ)=',JHIT(IZ)
         DO J = 1,JHIT(IZ)
            II = IP(IZ,J)
            PRINT *,' ID(II)=',ID(II)
            PRINT *,' XMR(IZ,J,1)=', XMR(IZ,J,1),
     &            ' YMR(IZ,J,1)=', YMR(IZ,J,1)   
            PRINT *,' XMR(IZ,J,2)=', XMR(IZ,J,2),
     &            ' YMR(IZ,J,2)=', YMR(IZ,J,2)   
        ENDDO
      ENDDO 
      
      RETURN
      END

*************************************************************************
      SUBROUTINE RECO_TRACKF(IDRES,IREADGEANT)
*************************************************************************
*
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
* 
      PARAMETER (MAXIDG=20000)
* 
      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5,MAXCAN=1000,NBCHAMBER=10,NTRMAX=500)
*      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
*      
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

* 
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
*      
      COMMON/VERIFGEANT/ITTROUGH(MAXTRK,NBSTATION),
     &                  IT_LIST(MAXTRK),IT_NP(MAXTRK),ITCHECK(MAXTRK),
     &                  ITRACK(MAXHITTOT)
*     
      COMMON/HCHHIT/HHIT(MAXHITCH),INDEXTAB(MAXHITCH),INDEXMAX
* 
      COMMON/PRECIS/EEXM(NBSTATION),EEYM(NBSTATION),EEPH(NBSTATION),
     &     EEAL(NBSTATION)
*
      COMMON/MEASUR/XMES(NBSTATION),YMES(NBSTATION),IZMES(NBSTATION),
     &              PHMES(NBSTATION),ALMES(NBSTATION),MPOS(NBSTATION),
     &              MANG(NBSTATION)
*
      COMMON/PLANE/XPL(NBSTATION,2),YPL(NBSTATION,2),
     &             PHPL(NBSTATION),ALPL(NBSTATION),CHI2PL
*
      COMMON/CANDIDAT/JCAN(NBSTATION,MAXCAN),JCANTYP(MAXCAN),
     &          EEX(MAXCAN),EEY(MAXCAN),EEP(MAXCAN),EEA(MAXCAN)
*     
      COMMON /VERTEX/ERRV,IVERTEX
*
      COMMON/TRACKSUM/NRES(5),NRESF,NTRMUALL,NMUONALL,NGHOSTALL,
     &     NTRACKFALL,NERRALL(NBSTATION),IR              
*            
      COMMON/TRACKFOUT/IEVOUT,NTREVT,JJOUT(NBCHAMBER,NTRMAX),
     &               ISTAT(NTRMAX),PXZOUT(NTRMAX),TPHIOUT(NTRMAX),
     &               TALAMOUT(NTRMAX),XVERTOUT(NTRMAX),YVERTOUT(NTRMAX),
     &               CHI2OUT(NTRMAX),
     &               XMESOUT(NBCHAMBER,NTRMAX),YMESOUT(NBCHAMBER,NTRMAX)  
     &              ,PXVOUT(NTRMAX),PYVOUT(NTRMAX),PZVOUT(NTRMAX)
*
      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP
*
      COMMON/DEBEVT/IDEBUG
      common/dstation/idstation

      DIMENSION PEST(5),PSTEP(5),NERR(NBSTATION),IFIND2(10),NMU45(2)
      DIMENSION NMU345(2),NMUONF(2)
*
      REAL*4 R2
      DATA R2/1./

** GEANT informations   
      DO I = 1,MAXTRK
         IT_LIST(I)= 0  ! ID of the GEANT tracks
      ENDDO
      DO I = 1,MAXHITTOT
         ITRACK(I) = 0  ! Track number associated to hit number I
      ENDDO

      NTRTOT = 0
*  BOUCLE SUR LES HITS
      DO IH = 1,NHITTOT
         DO IT = 1,NTRTOT            
            IF (ID(IH).EQ.IT_LIST(IT))  THEN
               IT1 = IT
               GOTO 4
            ENDIF   
         ENDDO
         NTRTOT = NTRTOT +1     ! NB de trace GEANT
         IT_LIST(NTRTOT) = ID(IH) ! ID DE LA TRACE NTRTOT
         IT_NP(NTRTOT) = 0      ! NOMBRE DE HITS PAR TRACE
         DO II=1,NBSTATION
            ITTROUGH(NTRTOT,II) = 0
         ENDDO
         IT1 = NTRTOT
 4       IT_NP(IT1) = IT_NP(IT1)+1 ! Number of crossed stations per track
         ITTROUGH(IT1,IZST(IH))=IH ! =IH si la trace IT touche le plan IZST
         ITRACK(IH) = IT1
      ENDDO
      IF (IDEBUG.GE.2) THEN
         PRINT *,'RECO_TRACKF nb total de trace GEANT:',NTRTOT
      ENDIF 

      NTRPART=0
      NTRMU = 0
      
      DO IT = 1,NTRTOT
         ITCHECK(IT) = 0
         IF ((ITTROUGH(IT,1)*ITTROUGH(IT,2)*ITTROUGH(IT,3)*
     &        ITTROUGH(IT,4)*ITTROUGH(IT,5)).NE.0)
     &        THEN  ! track crossing all stations
            NTRPART=NTRPART+1
            IH = ITTROUGH(IT,NBSTATION)
            IF (ITYP(IH).EQ.5.OR.ITYP(IH).EQ.6) THEN
               ISTAK = ID(IH)
               ISTAK = MOD(ISTAK,30000)
               ISTAK = MOD(ISTAK,10000)
*     test
               pt=sqrt(pvert1(ih)**2+pvert2(ih)**2)
               thet=datan2(pt,pvert3(ih))*180/3.1416
               pp=sqrt(pt**2+pvert3(ih)**2)
               
               IF (ISTAK.EQ.IDRES) THEN ! psi or upsilon
                  NTRMU = NTRMU+1
                  NTRMUALL = NTRMUALL+1
                  ITCHECK(IT) = 1
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      IF (IDEBUG.GE.2) THEN
         PRINT *,'RECO_TRACKF nb of part. GEANT crossing 5 st.:',
     &        NTRPART
         PRINT *,'RECO_TRACKF nb of muons/res. GEANT crossing 5 st.:',
     &        NTRMU
      ENDIF
      
**         CALL H_ACCEPTANCE(5)
**         CALL H_ACCEPTANCE(4)
**         PAUSE
         
      NCAN = 0 

*  Recherche 5 -> 4
      CALL ORDONNE_HIT(5,HCUT)

      DO IH = 1,INDEXMAX
         JJ = INDEXTAB(IH)
         X1 = XM(5,JJ)-(ZPLANE(5)-ZPLANE(4))*TAN(PHM(5,JJ))
         Y1 = YM(5,JJ)+(ZPLANE(5)-ZPLANE(4))*TAN(ALM(5,JJ))
     &                                            /COS(PHM(5,JJ))
         X2 = XM(5,JJ)-(ZPLANE(5)-ZPLANE(4)+DZ_PL(4))*TAN(PHM(5,JJ))
         Y2 = YM(5,JJ)+(ZPLANE(5)-ZPLANE(4)+DZ_PL(4))*TAN(ALM(5,JJ))
     &                                            /COS(PHM(5,JJ))

*  Domaine de recherche dans la st. 4
         HCONST = 0.0136*SQRT(0.06)/HTOP  ! -> X/X0=0.03 % / chamber
         EPH2 = 2.0*PHIPREC**2  + (HCONST*HHIT(JJ))**2  
         EAL2 = 2.0*ALAMPREC**2 + (HCONST*HHIT(JJ))**2
         EXM12 = (PHIPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4))**2 + XPREC**2
         EYM12 = (ALAMPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4))**2 + YPREC**2  
         EXM22 = (PHIPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4)+DZ_PL(4))**2 + XPREC**2
         EYM22 = (ALAMPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4)+DZ_PL(4))**2 + YPREC**2  
 
         EPH = SIGCUT*SQRT(EPH2)
         EAL = SIGCUT*SQRT(EAL2)
         EX1 = SIGCUT*SQRT(EXM12)
         EY1 = SIGCUT*SQRT(EYM12)
         EX2 = SIGCUT*SQRT(EXM22)
         EY2 = SIGCUT*SQRT(EYM22)
         
**         P2 = (HTOP/HHIT(JJ))**2

**         EPH=SIGCUT*SQRT(9.14D-7+1.2D-3/P2)
**         EAL=SIGCUT*SQRT(1.84D-4)
**         EX1=SIGCUT*SQRT(1.95D-2+6.37/P2)
**         EY1=SIGCUT*SQRT(3.89+151./P2)
**         EX2=EX1
**         EY2=EY2
 
* renvoie le num de hit de 4 le plus pres dans le domaine de recherche

         CALL DISTMIN4(X1,Y1,PHM(5,JJ),ALM(5,JJ),4,EX1,EY1,EPH,EAL,
     &        IFIND,IFIND2)
         P1 = PTOT(IP(5,JJ))
         CALL CHECK_HISTO4(11,4,IFIND,5,JJ,X1,Y1,PHM(5,JJ),ALM(5,JJ),P1,
     &                     EX1,EY1,EPH,EAL)   
         IF (IFIND.GT.0) THEN 
            CALL STOCK_CANDIDAT(5,JJ,4,IFIND,IFIND2,EX1,EY1,EPH,EAL,
     &           NCAN,1)
         ELSE
            CALL DISTMIN2(X1,Y1,X2,Y2,4,EX1,EY1,EX2,EY2,IFIND,IFIND2)
            CALL CHECK_HISTO2(0,4,IFIND,5,JJ,X1,Y1,X2,Y2,P1,EX1,EY1,
     &           EX2,EY2)   
            IF (IFIND.GT.0) THEN
               CALL STOCK_CANDIDAT(5,JJ,4,IFIND,IFIND2,EX1,EY1,EPH,EAL,
     &              NCAN,2)
            ENDIF 
         ENDIF   
      ENDDO
         
*   Recherche 4 -> 5

      CALL ORDONNE_HIT(4,HCUT)

      DO IH = 1,INDEXMAX
         JJ = INDEXTAB(IH)
         X1 = XM(4,JJ)-(ZPLANE(4)-ZPLANE(5))*TAN(PHM(4,JJ))
         Y1 = YM(4,JJ)+(ZPLANE(4)-ZPLANE(5))*TAN(ALM(4,JJ))
     &                                            /COS(PHM(4,JJ))
         X2 = XM(4,JJ)-(ZPLANE(4)-ZPLANE(5)+DZ_PL(5))*TAN(PHM(4,JJ))
         Y2 = YM(4,JJ)+(ZPLANE(4)-ZPLANE(5)+DZ_PL(5))*TAN(ALM(4,JJ))
     &                                            /COS(PHM(4,JJ))
*   Domaine de recherche dans la st. 5
         HCONST = 0.0136*SQRT(0.06)/HTOP  ! -> X/X0=0.03 / chamber
**         EPH2 = 2.0*PHIPREC**2  + (HCONST*HHIT(JJ))**2   
**         EAL2 = 2.0*ALAMPREC**2 + (HCONST*HHIT(JJ))**2
**         EX12 = 2.0*(XPREC/SQRT(2.))**2    
**     &          + (ZPLANE(5)-ZPLANE(4))**2/2.*EPH2
**         EY12 = 2.0*YPREC**2 + (ZPLANE(5)-ZPLANE(4))**2/2.*EAL2
**         EX22 = 2.0*(XPREC/SQRT(2.))**2    
**     &          + (ZPLANE(5)-ZPLANE(4)-DZ_PL(5))**2/2.*EPH2
**         EY22 = 2.0*YPREC**2 + (ZPLANE(5)-ZPLANE(4)-DZ_PL(5))**2/2.
**     &        *EAL2
         EPH2 = 2.0*PHIPREC**2  + (HCONST*HHIT(JJ))**2  
         EAL2 = 2.0*ALAMPREC**2 + (HCONST*HHIT(JJ))**2
         EXM12 = (PHIPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4))**2 + XPREC**2
         EYM12 = (ALAMPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4))**2 + YPREC**2  
         EXM22 = (PHIPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4)-DZ_PL(5))**2 + XPREC**2
         EYM22 = (ALAMPREC**2+(HCONST*HHIT(JJ))**2)*
     &        (ZPLANE(5)-ZPLANE(4)-DZ_PL(5))**2 + YPREC**2  

         EPH = SIGCUT*SQRT(EPH2)
         EAL = SIGCUT*SQRT(EAL2)
         EX1 = SIGCUT*SQRT(EXM12)
         EY1 = SIGCUT*SQRT(EYM12)
         EX2 = SIGCUT*SQRT(EXM22)
         EY2 = SIGCUT*SQRT(EYM22)


**         P2 = (HTOP/HHIT(JJ))**2

**         EPH=SIGCUT*SQRT(9.14D-7+1.2D-3/P2)
**         EAL=SIGCUT*SQRT(1.84D-4)
**         EX1=SIGCUT*SQRT(1.95D-2+6.37/P2)
**         EY1=SIGCUT*SQRT(3.89+151./P2)
**         EX2=EX1
**         EY2=EY2
*   Renvoie le num de hit de 5 le plus pres dans le domaine de recherche
         CALL DISTMIN2(X1,Y1,X2,Y2,5,EX1,EY1,EX2,EY2,IFIND,IFIND2)
         P1 = PTOT(IP(4,JJ))
         CALL CHECK_HISTO2(0,5,IFIND,4,JJ,X1,Y1,X2,Y2,P1,EX1,EY1,
     &        EX2,EY2)   
         IF (IFIND.GT.0) THEN
            DO ICAN=1,NCAN 
               IF (IFIND.EQ.JCAN(5,ICAN).AND.JJ.EQ.JCAN(4,ICAN)) GOTO 40 
**               IF (JJ.EQ.JCAN(4,ICAN).AND.
**     &              ABS(XM(5,JCAN(5,ICAN))-XM(5,IFIND)).LT.XPREC/10.)
**     &            GO TO 40 ! elimine les doubles comptages de traces ccc
**               IF (IFIND.EQ.JCAN(5,ICAN).AND.
**     &              ABS(XM(4,JCAN(4,ICAN))-XM(4,JJ)).LT.XPREC/10.)
**     &            GO TO 40 ! elimine les doubles comptages de traces ccc

               DIST1 = SQRT(((XM(5,JCAN(5,ICAN))-XM(5,IFIND))
     &              /(0.1*XPREC))**2+((YM(5,JCAN(5,ICAN))-YM(5,IFIND))
     &              /(0.1*YPREC))**2)
               DIST2 = SQRT(((XM(4,JCAN(4,ICAN))-XM(4,JJ))
     &              /(0.1*XPREC))**2+((YM(4,JCAN(4,ICAN))-YM(4,JJ))
     &              /(0.1*YPREC))**2)
               IF (DIST1.LT.2..AND.DIST2.LT.2.)
     &            GO TO 40 ! elimine les doubles comptages de traces ccc
            ENDDO   
            CALL STOCK_CANDIDAT(4,JJ,5,IFIND,IFIND2,EX1,EY1,EPH,EAL,NCAN
     &           ,3)
         ENDIF
40       CONTINUE           
      ENDDO


      NMU45(1) = 0
      NMU45(2) = 0
      DO ICAN = 1,NCAN
         JJ4 = JCAN(4,ICAN)  
         JJ5 = JCAN(5,ICAN)  
         IF (JJ4.GT.0.AND.JJ5.GT.0) THEN
            ID4 = ID(IP(4,JJ4))
            ID5 = ID(IP(5,JJ5))
            IT = ITRACK(IP(5,JJ5))
            IF (ITCHECK(IT).EQ.1) THEN
               IF (ID4.EQ.ID5) THEN
                  IF (ITYP(IP(5,JJ5)).EQ.5) NMU45(1) = 1
                  IF (ITYP(IP(5,JJ5)).EQ.6) NMU45(2) = 1
               ENDIF   
            ENDIF
         ENDIF   
      ENDDO
      IF (NMU45(1).GE.1.AND.NMU45(2).EQ.1) NRES(1) = NRES(1)+1   

      IF (IDEBUG.GE.2) THEN         
         PRINT *,'RECO_TRACKF nb candidat recherche 4->5 et 5->4 :'
     &        ,NCAN
         PRINT *,'RECO_TRACKF nb of good muons 4->5 et 5->4 :'
     &        ,(NMU45(1)+NMU45(2))
        
      ENDIF


      NMU345(1) = 0
      NMU345(2) = 0
*
* -- Boucle sur les candidats (4,5) NCAN
*         
      DO I = 1,NBSTATION
         NERR(I) = 0
      ENDDO
      NMUONF(1) = 0
      NMUONF(2) = 0
      NGHOST = 0
      NTRACKF = 0
     
**      GO TO 125 ! CCC

      DO ICAN = 1,NCAN
         JJ1 = 0
         JJ2 = 0
         JJ3 = 0
         DO ICH = 1,NBSTATION
            MPOS(ICH) = 0
            MANG(ICH) = 0
         ENDDO
         MPOS(4) = 1
         MPOS(5) = 1
         MANG(4) = 1
         MANG(5) = 1
         IF (JCANTYP(ICAN).EQ.2) MANG(4)=0
         IF (JCANTYP(ICAN).EQ.3) MANG(5)=0
         EEXM(5) = EEX(ICAN)/SIGCUT
         EEYM(5) = EEY(ICAN)/SIGCUT
         EEPH(5) = EEP(ICAN)/SIGCUT
         EEAL(5) = EEA(ICAN)/SIGCUT
         EEXM(4) = EEX(ICAN)/SIGCUT
         EEYM(4) = EEY(ICAN)/SIGCUT
         EEPH(4) = EEP(ICAN)/SIGCUT
         EEAL(4) = EEA(ICAN)/SIGCUT
         JJ5 = JCAN(5,ICAN)
         JJ4 = JCAN(4,ICAN)
         P = PTOT(IP(5,JJ5))
         IF (IZM(4,JJ4).EQ.1) THEN
            ZPL4 = ZPLANE(4)
         ELSE
            ZPL4 = ZPLANE(4)-DZ_PL(4)
         ENDIF   
         IF (IZM(5,JJ5).EQ.1) THEN
            ZPL5 = ZPLANE(5)
         ELSE
            ZPL5 = ZPLANE(5)-DZ_PL(5)
         ENDIF   
         TPHEST = (XM(5,JJ5) - XM(4,JJ4))/(ZPL5-ZPL4)
         PHEST = ATAN(TPHEST)
         TALEST = -(YM(5,JJ5) - YM(4,JJ4))*COS(PHEST)
     &               /(ZPL5-ZPL4)
         PXZEST = -HTOP/(XM(4,JJ4) - ZPL4*TPHEST)
         PEST(1) = 1.0/PXZEST
         PEST(2) = TPHEST - ALPHATOP/PXZEST ! PHI emission au vertex
         PEST(3) = TALEST       ! tan(lambda)     !h=zm*ang deviation alpha
         PEST(4) = 0.0          !alpha=qbl/pxz
         PEST(5) = 0.0
         PSTEP(1) = 0.003       ! =d(1/p)=delta(p)/p**2
         PSTEP(2) = 0.001       ! 0.5 degre
         PSTEP(3) = 0.001       ! 0.5 degre
         PSTEP(4) = 0.0
         PSTEP(5) = 1.0
         XMES(4) = XM(4,JJ4)
         YMES(4) = YM(4,JJ4)
         IZMES(4) = IZM(4,JJ4)
         PHMES(4) = PHM(4,JJ4)
         ALMES(4) = ALM(4,JJ4)
         XMES(5) = XM(5,JJ5)
         YMES(5) = YM(5,JJ5)
         IZMES(5) = IZM(5,JJ5)
         PHMES(5) = PHM(5,JJ5)
         ALMES(5) = ALM(5,JJ5)

         IVERTEX = 0 ! Vertex = (0,0,0)     

* -- Fit 4,5,V pour la recherche dans 3

         CALL TRACKF_FIT(IVERTEX,PEST,PSTEP,PXZINV45,TPHI45,TALAM45,
     &        XVERT,YVERT)
         
* -- Recherche dans la station 3
         
         P2 = (1.0D0 + TALAM45**2)/(PXZINV45**2) ! P**2

         if (idstation.eq.8) then !  DZ_CH = 8 CM
            EPH=SIGCUT*SQRT(3.6D-6+0.011/P2)
            EAL=SIGCUT*SQRT(6.85D-4)
            EXM=SIGCUT*SQRT(0.034+446./P2)
            EYM=SIGCUT*SQRT(0.049+354./P2)
         else if (idstation.eq.20) then ! DZ_CH = 20 CM
            EPH=SIGCUT*SQRT(4.1D-7+0.015/P2)
            EAL=SIGCUT*SQRT(1.04D-4)
            EXM=SIGCUT*SQRT(0.0+459./P2)
            EYM=SIGCUT*SQRT(0.042+345./P2)
         end if

         EEXM(3) = EXM/SIGCUT
         EEYM(3) = EYM/SIGCUT
         EEPH(3) = EPH/SIGCUT
         EEAL(3) = EAL/SIGCUT
**         DO IEX = 4,5
**            EEXM(IEX) = EEXM(3) 
**            EEYM(IEX) = EEYM(3)
**            EEPH(IEX) = EEPH(3)
**            EEAL(IEX) = EEAL(3)
**         ENDDO   
         X1 = XPL(3,1)
         Y1 = YPL(3,1) 
         X2 = XPL(3,2)
         Y2 = YPL(3,2) 
         PHI1 = PHPL(3)
         ALAM1 = ALPL(3)

         CALL DISTMIN4(X1,Y1,PHI1,ALAM1,3,EXM,EYM,EPH,EAL,IPA3,IFIND2)
         
         P1 = PTOT(IP(5,JJ5))

         CALL CHECK_HISTO4(61,3,IPA3,5,JJ5,X1,Y1,PHI1,ALAM1,P1,
     &        EXM,EYM,EPH,EAL)
         
         IF (IPA3.NE.0) THEN
            JJ3 = IPA3
            XMES(3) = XM(3,JJ3)
            YMES(3) = YM(3,JJ3)
            IZMES(3) = IZM(3,JJ3)
            PHMES(3) = PHM(3,JJ3)
            ALMES(3) = ALM(3,JJ3)
            MPOS(3) = 1
            MANG(3) = 1
            GO TO 124
         ELSE
            CALL DISTMIN2(X1,Y1,X2,Y2,3,EXM,EYM,0.D0,0.D0,IP3,IFIND2)
            CALL CHECK_HISTO2(0,3,IP3,5,JJ5,X1,Y1,X2,Y2,P1,EXM,EYM,
     &           0.D0,0.D0)   
         ENDIF
         IF (IP3.NE.0) THEN
            JJ3 = IP3
            XMES(3) = XM(3,JJ3)
            YMES(3) = YM(3,JJ3)
            IZMES(3) = IZM(3,JJ3)
            MPOS(3) = 1
            MANG(3) = 0
         ELSE
            GO TO 123
         ENDIF     
         
         
*     -- Fit 3,4,5 pour la recherche dans 2 
         
 124     CONTINUE
         
         IF (JJ5.GT.0.AND.JJ4.GT.0.AND.JJ3.GT.0) THEN           
            ID4 = ID(IP(4,JJ4))
            ID5 = ID(IP(5,JJ5))
            ID3 = ID(IP(3,JJ3))
            IT = ITRACK(IP(5,JJ5))
            IF (ITCHECK(IT).EQ.1) THEN
               IF ((ID5.EQ.ID3).AND.(ID5.EQ.ID4)) THEN
                  IF (ITYP(IP(5,JJ5)).EQ.5) NMU345(1) = 1
                  IF (ITYP(IP(5,JJ5)).EQ.6) NMU345(2) = 1
               ENDIF   
            ENDIF
         ENDIF   

         IVERTEX = 1 ! CCC    

         PEST(1) = PXZINV45
         PEST(2) = TPHI45
         PEST(3) = TALAM45
         PEST(4) = 0.
         PEST(5) = 0.
         PSTEP(1) = 0.003
         PSTEP(2) = 0.001
         PSTEP(3) = 0.001
         PSTEP(4) = 1.
         PSTEP(5) = 1. 

         CALL TRACKF_FIT(IVERTEX,PEST,PSTEP,PXZINV345,TPHI345,
     &        TALAM345,XVERT,YVERT)
         
*     -- Recherche dans la st. 2
         
         P2 = (1.0D0 + TALAM345**2)/(PXZINV345**2) ! P**2

         if (idstation.eq.8) then !  DZ_CH = 8 CM
            EPH=SIGCUT*SQRT(3.63D-6+4.8D-3/P2)
            EAL=SIGCUT*SQRT(6.49D-4)
            EXM=SIGCUT*SQRT(1.64D-2+821./P2)
            EYM=SIGCUT*SQRT(4.83D-2+866./P2)
         else if (idstation.eq.20) then ! DZ_CH = 20 CM
            EPH=SIGCUT*SQRT(5.78D-7+5.97D-3/P2)
            EAL=SIGCUT*SQRT(1.D-4)
            EXM=SIGCUT*SQRT(0.+1453./P2)
            EYM=SIGCUT*SQRT(2.78D-2+504./P2)
         end if

         EEXM(2) = EXM/SIGCUT
         EEYM(2) = EYM/SIGCUT
         EEPH(2) = EPH/SIGCUT
         EEAL(2) = EAL/SIGCUT
**         DO IEX = 3,5
**            EEXM(IEX) = EEXM(2) 
**            EEYM(IEX) = EEYM(2)
**            EEPH(IEX) = EEPH(2)
**            EEAL(IEX) = EEAL(2)
**         ENDDO   

         X1 = XPL(2,1)
         Y1 = YPL(2,1) 
         PHI1 = PHPL(2)
         ALAM1 = ALPL(2)
         CALL DISTMIN4(X1,Y1,PHI1,ALAM1,2,EXM,EYM,EPH,EAL,IPA2,IFIND2)
         P1 = PTOT(IP(5,JJ5))
         CALL CHECK_HISTO4(21,2,IPA2,5,JJ5,X1,Y1,PHI1,ALAM1,P1,
     &        EXM,EYM,EPH,EAL)   
*     -- Recherche dans la st. 1

         if (idstation.eq.8) then !  DZ_CH = 8 CM
            EPH=SIGCUT*SQRT(3.54D-6+4.49D-3/P2)
            EAL=SIGCUT*SQRT(6.14D-4)
            EXM=SIGCUT*SQRT(6.43D-2+986./P2)  
            EYM=SIGCUT*SQRT(4.66D-2+1444./P2) 
         else if (idstation.eq.20) then ! DZ_CH = 20 CM
            EPH=SIGCUT*SQRT(4.96D-7+5.87D-3/P2)
            EAL=SIGCUT*SQRT(1.D-4)
            EXM=SIGCUT*SQRT(6.98D-4+1467./P2)  
            EYM=SIGCUT*SQRT(5.22D-2+1013./P2)  
         end if

         EEXM(1) = EXM/SIGCUT
         EEYM(1) = EYM/SIGCUT
         EEPH(1) = EPH/SIGCUT
         EEAL(1) = EAL/SIGCUT
**         DO IEX = 2,5
**            EEXM(IEX) = EEXM(1) 
**            EEYM(IEX) = EEYM(1)
**            EEPH(IEX) = EEPH(1)
**            EEAL(IEX) = EEAL(1)
**         ENDDO   

         X1 = XPL(1,1)
         Y1 = YPL(1,1) 
         PHI1 = PHPL(1)
         ALAM1 = ALPL(1)
         CALL DISTMIN4(X1,Y1,PHI1,ALAM1,1,EXM,EYM,EPH,EAL,IPA1,IFIND2)
         CALL CHECK_HISTO4(31,1,IPA1,5,JJ5,X1,Y1,PHI1,ALAM1,P1,
     &        EXM,EYM,EPH,EAL)   
*     -- A partir de P+A dans la st. 2     
         IPA2PA1 = 0
         IF (IPA2.GT.0) THEN
            PEST(1) = PXZINV345
            PEST(2) = TPHI345
            PEST(3) = TALAM345
            PEST(4) = 0.
            PEST(5) = 0.
            PSTEP(1) = 0.003
            PSTEP(2) = 0.001
            PSTEP(3) = 0.001
            PSTEP(4) = 1.
            PSTEP(5) = 1. 
            XMES(2) = XM(2,IPA2)
            YMES(2) = YM(2,IPA2)
            IZMES(2) = IZM(2,IPA2)
            PHMES(2) = PHM(2,IPA2)
            ALMES(2) = ALM(2,IPA2)
            MPOS(2) = 1
            MANG(2) = 1
            IVERTEX = 1 
*     -- Fit V,2,3,4,5 pour la recherche dans 1 
            
            CALL TRACKF_FIT(IVERTEX,PEST,PSTEP,PXZINV,TPHI,TALAM,
     &           XVERT,YVERT)
            
*     !!! ATTENTION aux erreurs
            
            P2 = (1.0D0 + TALAM**2)/(PXZINV**2)

            if (idstation.eq.8) then !  DZ_CH = 8 CM
               EPH=SIGCUT*SQRT(3.15D-6+9.21D-3/P2)
               EAL=SIGCUT*SQRT(5.94D-4)
               EXM=SIGCUT*SQRT(4.16D-2+182./P2)
               EYM=SIGCUT*SQRT(2.13)
            else if (idstation.eq.20) then ! DZ_CH = 20 CM
               EPH=SIGCUT*SQRT(9.58D-7+8.93D-3/P2)
               EAL=SIGCUT*SQRT(1.D-4)
               EXM=SIGCUT*SQRT(1.92D-2+103.3/P2)
               EYM=SIGCUT*SQRT(6.3D-2+81.1/P2)
            end if

            EEXM(1) = EXM/SIGCUT
            EEYM(1) = EYM/SIGCUT
            EEPH(1) = EPH/SIGCUT
            EEAL(1) = EAL/SIGCUT
**            DO IEX = 2,5
**               EEXM(IEX) = EEXM(1) 
**               EEYM(IEX) = EEYM(1)
**               EEPH(IEX) = EEPH(1)
**               EEAL(IEX) = EEAL(1)
**            ENDDO   

            X1 = XPL(1,1)
            Y1 = YPL(1,1) 
            X2 = XPL(1,2)
            Y2 = YPL(1,2) 
            PHI1 = PHPL(1)
            ALAM1 = ALPL(1)
            CALL DISTMIN4(X1,Y1,PHI1,ALAM1,1,EXM,EYM,EPH,EAL,
     &           IPA2PA1,IFIND2)
            CALL CHECK_HISTO4(41,1,IPA2PA1,5,JJ5,X1,Y1,PHI1,ALAM1,
     &           P1,EXM,EYM,EPH,EAL)
            IF (IPA2PA1.GT.0) THEN
               JJ2 = IPA2
               JJ1 = IPA2PA1

               
               XMES(1) = XM(1,JJ1)
               YMES(1) = YM(1,JJ1)
               IZMES(1) = IZM(1,JJ1)
               PHMES(1) = PHM(1,JJ1)
               ALMES(1) = ALM(1,JJ1)
               MPOS(1) = 1
               MANG(1) = 1
               GOTO 123
            ELSE
               CALL DISTMIN2(X1,Y1,X2,Y2,1,EXM,EYM,0.D0,0.D0,IPA2P1,
     &              IFIND2)
               CALL CHECK_HISTO2(0,1,IPA2P1,5,JJ5,X1,Y1,X2,Y2,P1,
     &              EXM,EYM,0.D0,0.D0)   
            ENDIF 
         ENDIF   
*     --  A partir de P+A dans la st. 1     
         IPA1PA2 = 0
         IF (IPA1.GT.0) THEN
            PEST(1) = PXZINV345
            PEST(2) = TPHI345
            PEST(3) = TALAM345
            PEST(4) = 0.
            PEST(5) = 0.
            PSTEP(1) = 0.003
            PSTEP(2) = 0.001
            PSTEP(3) = 0.001
            PSTEP(4) = 1.
            PSTEP(5) = 1. 
            XMES(1) = XM(1,IPA1)
            YMES(1) = YM(1,IPA1)
            IZMES(1) = IZM(1,IPA1)
            PHMES(1) = PHM(1,IPA1)
            ALMES(1) = ALM(1,IPA1)
            MPOS(1) = 1
            MANG(1) = 1
            MPOS(2) = 0
            MANG(2) = 0
            IVERTEX = 1 
*     -- Fit V,1,3,4,5 pour la recherche dans 2 
            
            CALL TRACKF_FIT(IVERTEX,PEST,PSTEP,PXZINV,TPHI,TALAM,
     &           XVERT,YVERT) 
            
*     !!! ATTENTION aux erreurs
            
            P2 = (1.0D0 + TALAM**2)/(PXZINV**2)

            if (idstation.eq.8) then !  DZ_CH = 8 CM
               EPH=SIGCUT*SQRT(3.15D-6+9.21D-3/P2)
               EAL=SIGCUT*SQRT(5.94D-4)
               EXM=SIGCUT*SQRT(4.16D-2+182./P2)
               EYM=SIGCUT*SQRT(2.13)
            else if (idstation.eq.20) then ! DZ_CH = 20 CM
               EPH=SIGCUT*SQRT(9.58D-7+8.93D-3/P2)
               EAL=SIGCUT*SQRT(1.D-4)
               EXM=SIGCUT*SQRT(4.0D-2+11.4/P2)
               EYM=SIGCUT*SQRT(5.5D-2+44.35/P2)
            end if

            EEXM(2) = EXM/SIGCUT
            EEYM(2) = EYM/SIGCUT
            EEPH(2) = EPH/SIGCUT
            EEAL(2) = EAL/SIGCUT
**            DO IEX = 1,5
**               EEXM(IEX) = EEXM(2) 
**               EEYM(IEX) = EEYM(2)
**               EEPH(IEX) = EEPH(2)
**               EEAL(IEX) = EEAL(2)
**            ENDDO   

            X1 = XPL(2,1)
            Y1 = YPL(2,1) 
            X2 = XPL(2,2)
            Y2 = YPL(2,2) 
            PHI1 = PHPL(2)
            ALAM1 = ALPL(2)
            
            CALL DISTMIN4(X1,Y1,PHI1,ALAM1,2,EXM,EYM,EPH,EAL,
     &           IPA1PA2,IFIND2)
            
            CALL CHECK_HISTO4(51,2,IPA1PA2,5,JJ5,X1,Y1,PHI1,ALAM1,
     &           P1,EXM,EYM,EPH,EAL)
            IF (IPA1PA2.GT.0) THEN
               JJ1 = IPA1
               JJ2 = IPA1PA2
               XMES(2) = XM(2,JJ2)
               YMES(2) = YM(2,JJ2)
               IZMES(2) = IZM(2,JJ2)
               PHMES(2) = PHM(2,JJ2)
               ALMES(2) = ALM(2,JJ2)
               MPOS(2) = 1
               MANG(2) = 1
               GOTO 123
            ELSE
               CALL DISTMIN2(X1,Y1,X2,Y2,2,EXM,EYM,0.D0,0.D0,IPA1P2,
     &              IFIND2)
               CALL CHECK_HISTO2(0,2,IPA1P2,5,JJ5,X1,Y1,X2,Y2,P1,
     &              EXM,EYM,0.D0,0.D0) 
            ENDIF 
         ENDIF
*     -- Selection par type de candidat            
         IF (IPA1PA2.EQ.0.OR.IPA2PA1.EQ.0) THEN
            IF (IPA2.GT.0.AND.IPA2P1.GT.0) THEN
               JJ2 = IPA2
               JJ1 = IPA2P1
               XMES(1) = XM(1,JJ1)
               YMES(1) = YM(1,JJ1)
               IZMES(1) = IZM(1,JJ1)
               MPOS(1) = 1
               MANG(1) = 0
               MPOS(2) = 1
               MANG(2) = 1
            ELSEIF(IPA1.GT.0.AND.IPA1P2.GT.0) THEN
               JJ1 = IPA1
               JJ2 = IPA1P2
               XMES(2) = XM(2,JJ2)
               YMES(2) = YM(2,JJ2)
               IZMES(2) = IZM(2,JJ2)
               MPOS(1) = 1
               MANG(1) = 1
               MPOS(2) = 1
               MANG(2) = 0
            ENDIF
         ENDIF
         
 123     CONTINUE 
***   
         NTRACKFOLD = NTRACKF
         NMUONFOLD = NMUONF(1)+NMUONF(2)
         NGHOSTOLD = NGHOST
         CALL CHECK_MUON(JJ1,JJ2,JJ3,JJ4,JJ5,NTRACKF,NMUONF,
     &        NGHOST,NERR,PXZINV,TPHI,TALAM,XVERT,YVERT)
         
         IF (NTRACKF.GT.NTRACKFOLD) THEN
            ISTAT(NTRACKF) = 0 
            IF ((NMUONF(1)+NMUONF(2)).GT.NMUONFOLD) ISTAT(NTRACKF) = 1 ! Good muon
            IF (NGHOST.GT.NGHOSTOLD) ISTAT(NTRACKF) = 2 ! ghost track
            PXZOUT(NTRACKF) = 1./PXZINV
            TPHIOUT(NTRACKF) = TPHI  
            TALAMOUT(NTRACKF) = TALAM 
            XVERTOUT(NTRACKF) = XVERT 
            YVERTOUT(NTRACKF) = YVERT
            PXVOUT(NTRACKF) = PVERT1(IP(5,JJ5))
            PYVOUT(NTRACKF) = PVERT2(IP(5,JJ5))
            PZVOUT(NTRACKF) = PVERT3(IP(5,JJ5))
            CHI2OUT(NTRACKF) = CHI2PL
         ENDIF
***   
      ENDDO                     ! end loop on candidats NCAN
      
 125  IF (NMU345(1).EQ.1.AND.NMU345(2).EQ.1) NRES(2) = NRES(2) + 1

      IF (IDEBUG.GE.2) THEN         
         PRINT *,'RECO_TRACKF nb of good muons 3 4 5 :'
     &        ,(NMU345(1)+NMU345(2))
        
      ENDIF
     
      IF (IDEBUG.GE.2) THEN
         PRINT *,'RECO_TRACKF nb of track/evt     :',NTRACKF        
         PRINT *,'RECO_TRACKF nb of good muon/evt :',NMUONF(1)+NMUONF(2)         
         PRINT *,'RECO_TRACKF nb of ghost/evt     :',NGHOST         
         DO I = 1,4
            PRINT *,'RECO_TRACKF nb of error in st',I,'/evt:',NERR(I)
         ENDDO 
      ENDIF           
      IF (NTRMU.GE.2) NRES(5) = NRES(5)+1
      IF ((NMUONF(1)+NMUONF(2)).EQ.2) NRESF = NRESF+1
      NMUONALL = NMUONALL+NMUONF(1)+NMUONF(2)
      NGHOSTALL = NGHOSTALL+NGHOST
      NTRACKFALL = NTRACKFALL+NTRACKF
      DO I = 1,4
         NERRALL(I) = NERRALL(I)+NERR(I)
      ENDDO  
      
      CALL TRACKF_OUT(IR,NTRACKF) 
***   
      RETURN
      END

**********************************************************************         
      SUBROUTINE TRACKF_OUT(IR,NTRACKF)
**********************************************************************         
*   
*  
**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5)
      PARAMETER (NBCHAMBER=10,NTRMAX=500)
  
      COMMON/RECOUT/JJO(NBCHAMBER,NTRMAX),XMESO(NBCHAMBER,NTRMAX),
     &              YMESO(NBCHAMBER,NTRMAX)

      COMMON/TRACKFOUT/IEVOUT,NTREVT,JJOUT(NBCHAMBER,NTRMAX),
     &               ISTAT(NTRMAX),PXZOUT(NTRMAX),TPHIOUT(NTRMAX),
     &               TALAMOUT(NTRMAX),XVERTOUT(NTRMAX),YVERTOUT(NTRMAX),
     &               CHI2OUT(NTRMAX),
     &               XMESOUT(NBCHAMBER,NTRMAX),YMESOUT(NBCHAMBER,NTRMAX)
     &              ,PXVOUT(NTRMAX),PYVOUT(NTRMAX),PZVOUT(NTRMAX)
**
      IEVOUT = IR  
      NTREVT = NTRACKF 
      IF (NTREVT.GT.0) THEN
         DO ITR = 1,NTREVT
            ICH = 0
            DO IST = 1,NBSTATION
               DO ILOOP = 1,2 
                  ICH = ICH + 1 
                  JJOUT(ICH,ITR) = JJO(ICH,ITR)
                  IF (JJOUT(ICH,ITR).GT.0) THEN
                      XMESOUT(ICH,ITR) = XMESO(ICH,ITR)
                      YMESOUT(ICH,ITR) = YMESO(ICH,ITR)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF


      RETURN
      END
**********************************************************************          
      SUBROUTINE CHECK_MUON(JJ1,JJ2,JJ3,JJ4,JJ5,
     &                  NTRACKF,NMUONF,NGHOST,NERR,PXZINV,TPHI,TALAM,
     &                  XVERT,YVERT)  
**********************************************************************          
*   Check muon track candidate using GEANT informations
*
**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
* --
      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5)
      PARAMETER (NBCHAMBER=10,NTRMAX=500)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
      
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

      COMMON/VERIFGEANT/ITTROUGH(MAXTRK,NBSTATION),
     &                  IT_LIST(MAXTRK),IT_NP(MAXTRK),ITCHECK(MAXTRK),
     &                  ITRACK(MAXHITTOT)
     
      COMMON/MEASUR/XMES(NBSTATION),YMES(NBSTATION),IZMES(NBSTATION),
     &              PHMES(NBSTATION),ALMES(NBSTATION),MPOS(NBSTATION),
     &              MANG(NBSTATION)

      COMMON/RECOUT/JJO(NBCHAMBER,NTRMAX),XMESO(NBCHAMBER,NTRMAX),
     &              YMESO(NBCHAMBER,NTRMAX)

      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)      

      COMMON/MAGNET/BL3,B0

      REAL*4 R1,R2
      DATA R1,R2/0.,1./
 
      DIMENSION NERR(NBSTATION),JJK(NBSTATION),NMUONF(2) 
      
*      print*,' *** appel de la subroutine check_muon ***'

      IF (JJ1.GT.0.AND.JJ2.GT.0.AND.JJ3.GT.0.AND.JJ4.GT.0
     &     .AND.JJ5.GT.0) THEN
         IDFIND = ID(IP(5,JJ5))
         IT = ITRACK(IP(5,JJ5))
         JJK(1) = JJ1
         JJK(2) = JJ2
         JJK(3) = JJ3
         JJK(4) = JJ4
         JJK(5) = JJ5
         NTRACKF = NTRACKF+1
         DO I = 1,NBCHAMBER
            JJO(I,NTRACKF) = 0
         ENDDO
         ICH1 = -1
         DO IST = 1,5
            ICH1 = ICH1 + 2  
            IF (XMR(IST,JJK(IST),1).LT.999.) THEN
               JJO(ICH1,NTRACKF) = JJK(IST)
               XMESO(ICH1,NTRACKF) = XMR(IST,JJK(IST),1)
               YMESO(ICH1,NTRACKF) = YMR(IST,JJK(IST),1)
               ICH2 = ICH1 + 1  
               IF (MANG(IST).EQ.1) THEN
                  JJO(ICH2,NTRACKF) = JJK(IST)
                  XMESO(ICH2,NTRACKF) = XMR(IST,JJK(IST),2)
                  YMESO(ICH2,NTRACKF) = YMR(IST,JJK(IST),2)
               ENDIF
            ELSE
               ICH2 = ICH1+1
               JJO(ICH2,NTRACKF) = JJK(IST)
               XMESO(ICH2,NTRACKF) = XMR(IST,JJK(IST),2)
               YMESO(ICH2,NTRACKF) = YMR(IST,JJK(IST),2)
            ENDIF
         ENDDO
         
         CALL CHFILL(700,SNGL(XVERT),R1,R2)
         CALL CHFILL(701,SNGL(YVERT),R1,R2)
         
         IF (ITCHECK(IT).EQ.1) THEN
            
            
            IF (ID(IP(1,JJ1)).EQ.IDFIND .AND.
     &           ID(IP(2,JJ2)).EQ.IDFIND .AND.
     &           ID(IP(3,JJ3)).EQ.IDFIND .AND.
     &           ID(IP(4,JJ4)).EQ.IDFIND) THEN

               IF (ITYP(IP(5,JJ5)).EQ.5) NMUONF(1) = 1 ! Bon muon
               IF (ITYP(IP(5,JJ5)).EQ.6) NMUONF(2) = 1 ! Bon muon
               
               PT = SQRT(PVERT1(IP(5,JJ5))**2+PVERT2(IP(5,JJ5))**2)
               PZ =  PVERT3(IP(5,JJ5))
               E = SQRT(PT**2+PZ**2+0.1056**2)
               YRAP = 0.5*DLOG((E+PZ)/(E-PZ))
               PHIM = DATAN2(PVERT2(IP(5,JJ5)),PVERT1(IP(5,JJ5)))
               CALL CHFILL(801,SNGL(YRAP),R1,R2)
               CALL CHFILL(901,SNGL(PT),R1,R2)
               CALL CHFILL(911,SNGL(PHIM),R1,R2)
            ELSE
               
               NGHOST = NGHOST+1 ! ghost
               
               PT = SQRT(PVERT1(IP(5,JJ5))**2+PVERT2(IP(5,JJ5))**2)
               PZ =  PVERT3(IP(5,JJ5))
               E = SQRT(PT**2+PZ**2+0.1056**2)
               YRAP = 0.5*DLOG((E+PZ)/(E-PZ))
               PHIM = DATAN2(PVERT2(IP(5,JJ5)),PVERT1(IP(5,JJ5)))
               CALL CHFILL(802,SNGL(YRAP),R1,R2)
               CALL CHFILL(902,SNGL(PT),R1,R2)
               CALL CHFILL(912,SNGL(PHIM),R1,R2)
            ENDIF
         ENDIF    
      ENDIF
      IF (ITCHECK(ITRACK(IP(5,JJ5))).EQ.1) THEN
         IF (JJ3.EQ.0) NERR(3) = NERR(3)+1
         IF (JJ3.NE.0) THEN

            IF (JJ1.EQ.0) NERR(1) = NERR(1)+1
            IF (JJ2.EQ.0) NERR(2) = NERR(2)+1
            IF (JJ1.EQ.0) print*,'hit not found stations 1', NERR(1) 
            IF (JJ2.EQ.0) print*,'hit not found stations 2', NERR(2) 

         ENDIF   
         PT = SQRT(PVERT1(IP(5,JJ5))**2+PVERT2(IP(5,JJ5))**2)
         PZ =  PVERT3(IP(5,JJ5))
         E = SQRT(PT**2+PZ**2+0.1056**2)
         YRAP = 0.5*DLOG((E+PZ)/(E-PZ))
         PHIM = DATAN2(PVERT2(IP(5,JJ5)),PVERT1(IP(5,JJ5)))
         CALL CHFILL(800,SNGL(YRAP),R1,R2)
         CALL CHFILL(900,SNGL(PT),R1,R2)
         CALL CHFILL(910,SNGL(PHIM),R1,R2)
      ENDIF
      
      RETURN
      END


*************************************************************************
      SUBROUTINE RECO_SUM
*************************************************************************
*
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5)
      PARAMETER (NBCHAMBER=10,NTRMAX=500)
*
      COMMON/TRACKFOUT/IEVOUT,NTREVT,JJOUT(NBCHAMBER,NTRMAX),
     &               ISTAT(NTRMAX),PXZOUT(NTRMAX),TPHIOUT(NTRMAX),
     &               TALAMOUT(NTRMAX),XVERTOUT(NTRMAX),YVERTOUT(NTRMAX),
     &               CHI2OUT(NTRMAX), 
     &               XMESOUT(NBCHAMBER,NTRMAX),YMESOUT(NBCHAMBER,NTRMAX)
     &              ,PXVOUT(NTRMAX),PYVOUT(NTRMAX),PZVOUT(NTRMAX)
*
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
*      
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT
*  
      COMMON/RECOUT/JJO(NBCHAMBER,NTRMAX),XMESO(NBCHAMBER,NTRMAX),
     &              YMESO(NBCHAMBER,NTRMAX)
*
      COMMON/TRACKSUM/NRES(5),NRESF,NTRMUALL,NMUONALL,NGHOSTALL,
     &     NTRACKFALL,NERRALL(NBSTATION),IR                
*
      COMMON/PRECSUM/NRESF1,NMUONALL1,NGHOSTALL1,NTRACKFALL1          
*
      REAL*4 PXR,PYR,PZR,ZVR,CHI2R,PXV,PYV,PZV
      COMMON/PAWCR4/IEVR,NTRACKR,ISTATR(NTRMAX),ISIGNR(NTRMAX),
     &              PXR(NTRMAX),PYR(NTRMAX),PZR(NTRMAX),ZVR(NTRMAX),
     &              CHI2R(NTRMAX),PXV(NTRMAX),PYV(NTRMAX),PZV(NTRMAX)
*
      COMMON/DEBEVT/IDEBUG
*      
      DIMENSION ISEL(NTRMAX)

*  
      CALL RECO_SELECT(ISEL)
*
      NMUF = 0
      NGHF = 0

      DO ITR = 1,NTREVT
         IF (ISEL(ITR).EQ.1) THEN
            NTRACKFALL1 =  NTRACKFALL1 + 1
            IF (ISTAT(ITR).EQ.1) THEN
               NMUF = NMUF + 1
               NMUONALL1 = NMUONALL1 + 1
            ELSEIF (ISTAT(ITR).EQ.2) THEN 
               NGHF = NGHF + 1
               NGHOSTALL1 = NGHOSTALL1 + 1
            ENDIF
         ENDIF   
      ENDDO
      
      IF (NMUF.GE.2) NRESF1 = NRESF1 + 1 
*
      NTRACKR = 0
      DO ITR = 1,NTREVT
         IF (ISEL(ITR).EQ.1) THEN
            NTRACKR = NTRACKR + 1
            ISTATR(NTRACKR) = ISTAT(ITR)
            ISIGNR(NTRACKR) = 1 
            IF (PXZOUT(ITR).LT.0.) ISIGNR(NTRACKR) = -1
            PXZ = ABS(PXZOUT(ITR))
            PHI = ATAN(TPHIOUT(ITR)) 
            ALAM = ATAN(TALAMOUT(ITR)) 
            PYR(NTRACKR) = PXZ*SIN(PHI)
            PXR(NTRACKR) = PXZ*TAN(ALAM)
            PZR(NTRACKR) = PXZ*COS(PHI)
            ZVR(NTRACKR) = 0.
            CHI2R(NTRACKR) =  CHI2OUT(ITR)
            PXV(NTRACKR) = PXVOUT(ITR)
            PYV(NTRACKR) = PYVOUT(ITR) 
            PZV(NTRACKR) = PZVOUT(ITR)   
         ENDIF   
      ENDDO     

      CALL CHFNT(IEVR,NTRACKR,ISTATR,ISIGNR,
     &     PXR,PYR,PZR,ZVR,CHI2R,PXV,PYV,PZV)

      IF (IDEBUG.GE.2) THEN 
         PRINT *,'RECO_SUM evt number :',IEVOUT
         PRINT *,'RECO_SUM nb of track /evt     :',NTRACKR        
         PRINT *,'RECO_SUM nb of good muon /evt :',NMUF         
         PRINT *,'RECO_SUM nb of ghost /evt     :',NGHF         
         IF (NTRACKR.GT.0) THEN
            DO ITR = 1,NTRACKR
               PRINT *,'RECO_SUM track number :',ITR
               PRINT *,'RECO_SUM CHI2OUT :',CHI2R(ITR)
               PRINT *,' PX GEANT= ', PXV(ITR),' PX RECONS= ',PYR(ITR)
               PRINT *,' PY GEANT= ', PYV(ITR),' PY RECONS= ',PXR(ITR)
               PRINT *,' PZ GEANT= ', PZV(ITR),' PZ RECONS= ',PZR(ITR)
               PXZV = SQRT( PYV(ITR)**2+PZV(ITR)**2)
               PXZR = SQRT( PXR(ITR)**2+PZR(ITR)**2) 
               PRINT *,' PXZ GEANT= ', PXZV,' PXZ RECONS= ',PXZR
               FIV=ATAN2(PYV(ITR),PZV(ITR))
               ALAMV=ATAN2(PXV(ITR),PXZV)
               FIR=ATAN2(PXR(ITR),PZR(ITR))
               ALAMR=ATAN2(PYR(ITR),PXZR)
**               PRINT *,' PHI GEANT= ',FIV,' PXZ RECONS= ',FIR
**               PRINT *,' ALAM GEANT= ',ALAMV,' ALAM RECONS= ',ALAMR
            ENDDO
         ENDIF
      ENDIF 

      RETURN
      END

*************************************************************************
      SUBROUTINE RECO_TERM
*************************************************************************
*
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER(NBSTATION=5)
*
      COMMON/TRACKSUM/NRES(5),NRESF,NTRMUALL,NMUONALL,NGHOSTALL,
     &     NTRACKFALL,NERRALL(NBSTATION),IR  
*              
      COMMON/PRECSUM/NRESF1,NMUONALL1,NGHOSTALL1,NTRACKFALL1          
*
      COMMON/DEBEVT/IDEBUG
*
      CHARACTER*50 FILEBKG,FILERES,FILEOUT,FILEMIN
*
      IF (IDEBUG.GE.1) THEN
         PRINT *,'    '
         PRINT *,'RECO_TERM ***** SUMMARY TRACK-FINDING *****'
         PRINT *,'RECO_TERM nb of resonances :',NRES(5)
         PRINT *,'RECO_TERM nb of resonances  45 :',NRES(1)
         PRINT *,'RECO_TERM nb of resonances 345 :',NRES(2)
         PRINT *,'RECO_TERM nb of resonances found :',NRESF
         PRINT *,'RECO_TERM nb of muon track       :',NTRMUALL
         PRINT *,'RECO_TERM nb of track found      :',NTRACKFALL
         PRINT *,'RECO_TERM nb of muon track found :',NMUONALL
         PRINT *,'RECO_TERM nb of ghost found      :',NGHOSTALL
         DO I = 1,4
            PRINT *,'RECO_TERM nb of error in st',I,':',NERRALL(I)
         ENDDO

         PRINT *,'    '
         PRINT *,'RECO_TERM ***** SUMMARY PRECISION *****'
         PRINT *,'RECO_TERM nb of resonances found :',NRESF1
         PRINT *,'RECO_TERM nb of track found      :',NTRACKFALL1
         PRINT *,'RECO_TERM nb of muon track found :',NMUONALL1
         PRINT *,'RECO_TERM nb of ghost found      :',NGHOSTALL1
      ENDIF
* 
      CALL HIST_CLOSED

      RETURN
      END

*************************************************************************
      SUBROUTINE RECO_PRECISION
*************************************************************************
*
*
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5,MAXCAN=1000,NTRMAX=500)
      PARAMETER (NPLANE=10) 
*
      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
*      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
*      
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT
*
      COMMON/TRACKFOUT/IEVOUT,NTREVT,JJOUT(NPLANE,NTRMAX),
     &               ISTAT(NTRMAX),PXZOUT(NTRMAX),TPHIOUT(NTRMAX),
     &               TALAMOUT(NTRMAX),XVERTOUT(NTRMAX),YVERTOUT(NTRMAX),
     &               CHI2OUT(NTRMAX),
     &               XMESOUT(NPLANE,NTRMAX),YMESOUT(NPLANE,NTRMAX)
     &              ,PXVOUT(NTRMAX),PYVOUT(NTRMAX),PZVOUT(NTRMAX)

*
      COMMON/MEAS/LPLANE(NPLANE),XMP(NPLANE),YMP(NPLANE)
*      
      COMMON/FCNOUT/PXZEA,ALAMEA,PHIEA,XEA,YEA,NPLU,CHI2
*
      COMMON/PRECCUT/PCUT,PTCUT,CHI2CUT
*
      DIMENSION PARMU(MAXCAN,NPLANE,2),LPLANEMU(MAXCAN,NPLANE)  
*
      IF (NTREVT.EQ.0) RETURN

      DO ITR = 1,NTREVT
         ICH = 0
         DO IST = 1,NBSTATION
            DO ILOOP = 1,2 
               ICH = ICH + 1 
**               print *,' ich=',ich
               IF (JJOUT(ICH,ITR).GT.0) THEN
                  LPLANEMU(ITR,ICH) = 1
                  PARMU(ITR,ICH,1) = XMESOUT(ICH,ITR)
                  PARMU(ITR,ICH,2) = YMESOUT(ICH,ITR)
**                  print *,' x,y ', PARMU(ITR,ICH,1),PARMU(ITR,ICH,2) 
               ELSE
                  LPLANEMU(ITR,ICH) = 0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
*
      NTRACK = 0
      DO ICAN = 1,NTREVT
         DO ICH = 1,NPLANE
            LPLANE(ICH) = LPLANEMU(ICAN,ICH)
            XMP(ICH) = PARMU(ICAN,ICH,1)
            YMP(ICH) = PARMU(ICAN,ICH,2)
         ENDDO 

         IF (LPLANE(1).GT.0) THEN
            X1 = XMP(1)
            Y1 = YMP(1)
            IPL1 = 1
         ELSE
            X1 = XMP(2)
            Y1 = YMP(2)
            IPL1 = 2
         ENDIF
         IF (LPLANE(3).GT.0) THEN
            X2 = XMP(3)
            Y2 = YMP(3)
            IPL2 = 3
         ELSE
            X2 = XMP(4)
            Y2 = YMP(4)
            IPL2 = 4
         ENDIF
         IF (LPLANE(7).GT.0) THEN
            X3 = XMP(7)
            IPL3 = 7
         ELSE
            X3 = XMP(8)
            IPL3 = 8
         ENDIF
         IF (LPLANE(9).GT.0) THEN
            X4 = XMP(9)
            IPL4 = 9
         ELSE
            X4 = XMP(10)
            IPL4 = 10
         ENDIF
   
         PHIAV = DATAN2((X2-X1),(ZPLANEP(IPL2)-ZPLANEP(IPL1)))          
         PHIAP = DATAN2((X4-X3),(ZPLANEP(IPL4)-ZPLANEP(IPL3)))
         DPHI = (PHIAP-PHIAV)
         ASIGN = 1.
         IF (DPHI.LT.0.) ASIGN = -1. ! CCC
         PXZ = CONST/DABS(DPHI)
** Cuts PXZ           
         IF (PXZ.LT.PCUT) GO TO 2       
            
         PXZINVI = ASIGN/PXZ ! CCC
**         PXZINVI = 1./PXZOUT(ICAN) ! CCC
**         PXZINVI = -1./49. ! CCC
         PHII = PHIAV
         ALAMI = DATAN2((Y2-Y1),DSQRT((X2-X1)**2
     &            +(ZPLANEP(IPL2)-ZPLANEP(IPL1))**2))
         XVR = X1
         YVR = Y1
**         print *,' avant prec_fit pxzi phii alami x y',1./ PXZINVI,
**     &         PHII, ALAMI ,XVR,YVR             
**         PRINT *,' X1 X2 X3 X4',X1,X2,X3,X4
**         PRINT *,' Z1 Z2 Z3 Z4',ZPLANEP(IPL1),ZPLANEP(IPL2),
**     &          ZPLANEP(IPL3),ZPLANEP(IPL4)
**         PRINT *,' CONST= ',CONST 


         IF (CHI2OUT(ICAN).GT.CHI2CUT) GO TO 2

** Fit des traces apres l'absorbeur           
         CALL PREC_FIT (PXZINVI,PHII,ALAMI,XVR,YVR,
     &          PXZINVF,PHIF,ALAMF,XVERTF,YVERTF,EPXZINV,EPHI,EALAM,
     &          EXVERT,EYVERT)
     
** Correction de Branson       
         CALL BRANSON(PXZEA,PHIEA,ALAMEA,XEA,YEA)
         
         PXZ1 = DABS(PXZEA)
         PX1 = PXZ1*DSIN(PHIEA)
         PY1 = PXZ1*DTAN(ALAMEA)
         PT1 = DSQRT(PX1**2+PY1**2)
** Cuts PT
         IF (PT1.LT.PTCUT) GO TO 2
** Cuts CHI2
         IF ((CHI2/FLOAT(2*NPLU-5)).GT.CHI2CUT) GO TO 2

         NTRACK = NTRACK + 1
         DO ICH = 1,NPLANE           
            JJOUT(ICH,NTRACK) =  JJOUT(ICH,ICAN)
            XMESOUT(ICH,NTRACK) =  XMESOUT(ICH,ICAN)
            YMESOUT(ICH,NTRACK) =  YMESOUT(ICH,ICAN)
         ENDDO
         ISTAT(NTRACK) = ISTAT(ICAN)
         PXZOUT(NTRACK) = PXZEA
         TPHIOUT(NTRACK) = DTAN(PHIEA)
         TALAMOUT(NTRACK) = DTAN(ALAMEA) 
         XVERTOUT(NTRACK) = XEA        
         YVERTOUT(NTRACK) = YEA 
         CHI2OUT(NTRACK) =  CHI2/FLOAT(2*NPLU-5) 


**         print *,' reco_precision pxz tphi talam xvert yvert chi2',
**     &       PXZEA,PHIEA,ALAMEA,
**     &       XEA,YEA,CHI2/FLOAT(2*NPLU-5)  
  2      CONTINUE 
      ENDDO
      NTREVT = NTRACK 

      RETURN
      END


************************************************************************
      SUBROUTINE ORDONNE_HIT(ICH,HCUT)
**************************************************************************
*    
*  Sort hits in station ICH according to the "impact parameter" HHIT
*
**************************************************************************
 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
     
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)      
      
      COMMON/VERIFGEANT/ITTROUGH(MAXTRK,NBSTATION),
     &                  IT_LIST(MAXTRK),IT_NP(MAXTRK),ITCHECK(MAXTRK),
     &                  ITRACK(MAXHITTOT)
     
      COMMON/HCHHIT/HHIT(MAXHITCH),INDEXTAB(MAXHITCH),INDEXMAX
*
      COMMON/DEBEVT/IDEBUG
* 
      REAL*4 H4(MAXHITCH) 
* tri des impulsion par ordre decroissant
* le tab INDEXTAB contient les j ordonnes
* INDEXMAX est l indice max du tableau = NBHIT si pas de contrainte

      JJ=0
      INDEXMAX = 0
* boucle sur le nombre de hits candidats de la station

      DO J=1,JHIT(ICH)

          IF (PHM(ICH,J).LT.6.3) THEN        !2pi=6.3 radian
            JJ=JJ+1
* calcul du h dans XOY a z=0
            HHIT(J)=ABS(XM(ICH,J)-ZPLANE(ICH)*PHM(ICH,J))
* cut en Pxz
            IF (HHIT(J).LT.HCUT) THEN
               INDEXMAX=INDEXMAX+1
               INDEXTAB(INDEXMAX)=J
            ELSEIF(ITCHECK(ITRACK(IP(ICH,J))).EQ.1) THEN
               IF (IDEBUG.GE.2) THEN
                  PRINT *,'ORDONNE_HIT rejet muon/res in st.',ICH,
     &                 ' h=',HHIT(J)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      NBHIT=JHIT(ICH)

      DO I = 1,NBHIT
         H4(I) = SNGL(HHIT(I))
      ENDDO   

      CALL SORTZV(H4,INDEXTAB,INDEXMAX,-1,0,1)

      DO I = 1,NBHIT
         HHIT(I) = DBLE(H4(I))
      ENDDO   
**      PRINT *,'ORDONNE st. numero',ICH
**      PRINT *,'ORDONNE nb de hits initiaux dans st.',ICH,':',JHIT(ICH)
**      PRINT *,'ORDONNE nb de hits avec mes. angulaire:',JJ
**      PRINT *,'ORDONNE nb de hits avec mes. ang. et cut en Pxz:',INDEXMAX
      IF (IDEBUG.GE.2) THEN 
         PRINT *,'ORDONNE_HIT nb de hits accepte dans st.',ICH,':',
     &           INDEXMAX
      ENDIF

      RETURN
      END
***********************************************************************************
      SUBROUTINE DISTMIN4(X1,Y1,PHI1,ALAM1,ICH,EX,EY,EPHI,ELAM,IFIND,
     &     IFIND2)
***********************************************************************************
*    Find the nearest hit in station ICH in the (X,Y,lambda,phi) phase space
*
***********************************************************************************
 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      PARAMETER (MAXHITCH=10000,NBSTATION=5)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
      DIMENSION IFIND2(10)

      IFIND = 0
      DO I = 1,10
         IFIND2(I) = 0
      ENDDO
     
      DISTMIN=4.
      NF = 0  
      DO I=1,JHIT(ICH)
         IF (PHM(ICH,I).LE.6.3) THEN ! vector measurement
            IF (ABS(PHI1-PHM(ICH,I)) .LT. EPHI .AND.
     &         ABS(ALAM1-ALM(ICH,I)) .LT. ELAM .AND.
     &         ABS(X1-XM(ICH,I)) .LT. EX .AND.
     &         ABS(Y1-YM(ICH,I)) .LT. EY) THEN
               DIST = ((PHI1-PHM(ICH,I))/EPHI)**2 +
     &                ((ALAM1-ALM(ICH,I))/ELAM)**2 +
     &                ((X1-XM(ICH,I))/EX)**2 +
     &                ((Y1-YM(ICH,I))/EY)**2
               NF = NF+1
               IF (NF.LE.10) IFIND2(NF) = I 
               IF (DIST .LT. DISTMIN) THEN
                  DISTMIN = DIST
                  IFIND = I
               ENDIF
            ENDIF
         ENDIF
      ENDDO
       
      RETURN
      END
***********************************************************************************
      SUBROUTINE DISTMIN2(X1,Y1,X2,Y2,ICH,EX1,EY1,EX2,EY2,IFIND,IFIND2)
***********************************************************************************
*    Find the nearest hit in station ICH in the (X,Y) space
*
***********************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      PARAMETER (MAXHITCH=10000,NBSTATION=5)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
      DIMENSION IFIND2(10)

      IFIND = 0
      DO I = 1,10
         IFIND2(I) = 0
      ENDDO

      DISTMIN=2.
      NF = 0
      DO I=1,JHIT(ICH)
         IF (IZM(ICH,I).EQ.1) THEN ! 1st chamber
            X = X1
            Y = Y1
         ELSE ! 2nd chamber
            X = X2
            Y = Y2
         ENDIF    
         EX = EX1
         EY = EY1
         IF (ICH.EQ.4.OR.ICH.EQ.5) THEN
            IF (IZM(ICH,I).EQ.1) THEN
               EX = EX1
               EY = EY1
            ELSE
               EX = EX2
               EY = EY2
            ENDIF    
         ENDIF    
         IF (ABS(X-XM(ICH,I)) .LT. EX .AND.
     &      ABS(Y-YM(ICH,I)) .LT. EY) THEN
            DIST = ((X-XM(ICH,I))/EX)**2 +
     &             ((Y-YM(ICH,I))/EY)**2
            NF = NF+1
            IF (NF.LE.10) IFIND2(NF) = I
            IF (DIST .LT. DISTMIN) THEN
               DISTMIN = DIST
               IFIND = I
            ENDIF
         ENDIF
      ENDDO
       
      RETURN
      END
********************************************************************************
      SUBROUTINE H_ACCEPTANCE(ICH)
********************************************************************************
* Etude de l'acceptance des resonnances en  fonction du H 
* dans la station ICH
*
*  INPUT :    ICH
*  OUTPUT :  Histo #1
********************************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
      PARAMETER (MAXHITCH=10000,MAXTRK=50000,MAXHITTOT=20000,
     &           NBSTATION=5)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
     
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

 
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
      
      COMMON/VERIFGEANT/ITTROUGH(MAXTRK,NBSTATION),
     &                  IT_LIST(MAXTRK),IT_NP(MAXTRK),ITCHECK(MAXTRK),
     &                  ITRACK(MAXHITTOT)
     
      COMMON/HCHHIT/HHIT(MAXHITCH),INDEXTAB(MAXHITCH),INDEXMAX

      REAL*4 R1,R2
      DATA R1,R2/0.,1./
 
      NMUONI = 0
      DO J = 1,JHIT(ICH)
         IF (ITYP(IP(ICH,J)).EQ.5.OR.ITYP(IP(ICH,J)).EQ.6) THEN
            ISTAK = ID(IP(ICH,J))
            ISTAK = MOD(ISTAK,30000)
            ISTAK = MOD(ISTAK,10000)
            IF (ISTAK.EQ.0) THEN
**              PRINT *,'ACCEPT. id du muon dans st.',ICH,':',ITYP(IP(ICH,J))
               NMUONI = NMUONI+1
            ENDIF
         ENDIF
      ENDDO
*      PRINT *,'ACCEPT. nb de muons/res total dans st.',ICH,':',NMUONI
*      pause

      DO IH = 1,500
         HCUT = IH
*   Sort hits in st. z
         CALL ORDONNE_HIT(ICH,HCUT)
         NMUON = 0
         DO IND = 1,INDEXMAX
            IIND = IP(ICH,INDEXTAB(IND))
            IDPART = ITYP(IIND)
            ISTAK = ID(IIND)
            ISTAK = MOD(ISTAK,30000)
            ISTAK = MOD(ISTAK,10000)
**            PRINT *,' IDPART=',IDPART,' ISTAK=',ISTAK
            IF (IDPART.EQ.5.OR.IDPART.EQ.6.AND.ISTAK.EQ.0) THEN
               NMUON = NMUON+1
            ENDIF
         ENDDO
         IF (NMUON.EQ.2.AND.NMUONI.EQ.2) THEN
            CALL CHFILL(ICH*100,SNGL(HCUT),R1,R2)
         ENDIF
      ENDDO
 
      RETURN
      END

********************************************************************************
       SUBROUTINE OLDFOLLOW(ZSTR,PEST,IFLAG,XPL,YPL,PHPL,ALPL)
********************************************************************************
*   Calculate the particle trajectory in the spectrometer and
*   (XPL,YPL,PHPL,ALPL)
*   for the 5 stations.
*       
********************************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
       PARAMETER(NBSTATION=5) 
*
       DIMENSION XPL(NBSTATION,2),YPL(NBSTATION,2),PHPL(NBSTATION),
     &           ALPL(NBSTATION),PEST(NBSTATION)

       COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
       
       COMMON /MEASUR/XMES(NBSTATION),YMES(NBSTATION),IZMES(NBSTATION),
     &                PHMES(NBSTATION),ALMES(NBSTATION),MPOS(NBSTATION),
     &                MANG(NBSTATION)
       COMMON /MAGNET/BL3,B0

       LOGICAL LFLAG, LFLAG1

       XSTR = PEST(4)
       YSTR = PEST(5)
       PXZINV = PEST(1)
       TPHI = PEST(2)
       PHI = ATAN(TPHI)
       TALAM = PEST(3)
       PXZ = 1.0/PXZINV
       PY = ABS(PXZ)*TALAM
       PX = -ABS(PXZ)*SIN(PHI)
       PZ = -ABS(PXZ)*COS(PHI)
       PXY = SQRT(PX**2 + PY**2)
       FI=ATAN2(PY,PX)
       SINFI = SIN(FI)
       COSFI = COS(FI)
       TTHET = PZ/PXY
       RS = PXY*(100.0/(0.299792458*BL3))
       IF(PXZINV.LT.0.0) RS = -RS
*       XC = XSTR + RS*SIN(FI)
*       YC = YSTR - RS*COS(FI)
       PX0 = PX
       PY0 = PY
       LFLAG = .TRUE.
       LFLAG1 = .TRUE.
*       PRINT *, XC,YC,RS,FI,TTHET,PXY,PZ

       DO J = 1,5

          IF (IFLAG.EQ.3 .OR. MPOS(J).EQ.1) THEN
             IF(ZPLANE(J) .GT. ZCOIL) THEN
*                DFI =  (ZPLANE(J)-ZSTR)/(TTHET*RS)
*                FIN = FI - DFI
*                XPL(J,1) = XC - RS*SIN(FIN)
*                YPL(J,1) = YC + RS*COS(FIN)
                 DFR =  (ZPLANE(J)-ZSTR)/TTHET
                 XPL(J,1) = XSTR + DFR*COSFI + 0.5D0*DFR*DFR*SINFI/RS
                 YPL(J,1) = YSTR + DFR*SINFI - 0.5D0*DFR*DFR*COSFI/RS
                 DFR2 =  (ZPLANE(J)-DZ_PL(J)-ZSTR)/TTHET
                 XPL(J,2) = XSTR + DFR2*COSFI + 0.5D0*DFR2*DFR2*SINFI/RS
                 YPL(J,2) = YSTR + DFR2*SINFI - 0.5D0*DFR2*DFR2*COSFI/RS
                IF (IFLAG.EQ.3 .OR. MANG(J).EQ.1) THEN
*                   PX=PXY*COS(FIN)
*                   PY=PXY*SIN(FIN)
                   PX = PX0 + DFR * (PY0 - 0.5D0*PX0*DFR/RS) / RS
                   PY = PY0 - DFR * (PX0 + 0.5D0*PY0*DFR/RS) / RS
                   PHPL(J)=ATAN(PX/PZ)
                   ALPL(J)=ATAN(PY/SQRT(PX**2+PZ**2))
                ENDIF
              ELSE 
                IF( LFLAG) THEN
*                   DFI =  (ZCOIL-ZSTR)/(TTHET*RS)
*                   FIN = FI - DFI
*                   XCOIL = XC - RS*SIN(FIN)
*                   YCOIL = YC + RS*COS(FIN)
                   DFR =  (ZCOIL-ZSTR)/TTHET
                   XCOIL = XSTR + DFR*COSFI + 0.5D0*DFR*DFR*SINFI/RS
                   YCOIL = YSTR + DFR*SINFI - 0.5D0*DFR*DFR*COSFI/RS
*                   PX=PXY*COS(FIN)
*                   PY=PXY*SIN(FIN)
                   PX = PX0 + DFR * (PY0 - 0.5D0*PX0*DFR/RS) / RS
                   PY = PY0 - DFR * (PX0 + 0.5D0*PY0*DFR/RS) / RS
                   PXZ = SQRT(PX**2 + PZ**2)
                   PHI=ATAN(PX/PZ)
                   TALAM = PY/PXZ
                   ALAM = ATAN(TALAM)
                   RD = PXZ*(100.0/(0.299792458*B0))
                   IF(PXZINV.LT.0.0) RD = -RD
                   ZC = ZCOIL - RD*SIN(PHI)
                   XC = XCOIL + RD*COS(PHI)
                   IF(ABS(ZMAGEND-ZC).GT.ABS(RD)) STOP 'FOLLOW'
                   LFLAG = .FALSE.
                ENDIF
                IF(ZPLANE(J) .GT. ZMAGEND) THEN 
                  FIN = ASIN((ZPLANE(J) - ZC)/RD)
                  XPL(J,1)= XC - RD*COS(FIN)
                  YPL(J,1)= YCOIL - RD*(FIN - PHI)*TALAM
                  FIN2 = ASIN((ZPLANE(J)-DZ_PL(J) - ZC)/RD)
                  XPL(J,2)= XC - RD*COS(FIN2)
                  YPL(J,2)= YCOIL - RD*(FIN2 - PHI)*TALAM
                  PHPL(J)=FIN
                  ALPL(J)=ALAM
                ELSE
                  IF (LFLAG1) THEN
                    FIN = ASIN((ZMAGEND - ZC)/RD)
                    XMAGEND = XC - RD*COS(FIN)
                    YMAGEND = YCOIL - RD*(FIN - PHI)*TALAM
                    TPHI = TAN(FIN)
                    CPHI = COS(FIN)
                    LFLAG1 = .FALSE.
                  ENDIF
                  XPL(J,1) = XMAGEND + (ZPLANE(J)-ZMAGEND)*TPHI
                  YPL(J,1) = YMAGEND - (ZPLANE(J)-ZMAGEND)*TALAM/CPHI
                  XPL(J,2) = XMAGEND + (ZPLANE(J)-DZ_PL(J)-ZMAGEND)*TPHI
                  YPL(J,2) = YMAGEND - (ZPLANE(J)-DZ_PL(J)-ZMAGEND)*
     &                      TALAM/CPHI
                  PHPL(J)=FIN
                  ALPL(J)=ALAM
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       RETURN
       END
********************************************************************************
       SUBROUTINE FOLLOW(ZSTR,PEST,IFLAG,XPL,YPL,PHPL,ALPL)
********************************************************************************
*   Calculate the particle trajectory in the spectrometer 
*   (XPL,YPL,PHPL,ALPL)
*   for the 5 stations.
*       Runge Kutta
********************************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
       PARAMETER(NBSTATION=5,NPLANE=10) 
*
       DIMENSION XPL(NBSTATION,2),YPL(NBSTATION,2),PHPL(NBSTATION),
     &           ALPL(NBSTATION),PEST(NBSTATION)

       COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)

       COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
       
       COMMON /MEASUR/XMES(NBSTATION),YMES(NBSTATION),IZMES(NBSTATION),
     &                PHMES(NBSTATION),ALMES(NBSTATION),MPOS(NBSTATION),
     &                MANG(NBSTATION)
 
       
       DIMENSION VECT(7),VOUT(7)
       
       STEP = 6.                ! 1 cm
       NSTEPMAX = 5000
                
       ASIGN = 1.
       IF (PEST(1).LT.0.) ASIGN = -1.
       TPHI = -1.*PEST(2)
       PHI = DATAN(TPHI)
       TALAM = PEST(3)
       ALAM = DATAN(TALAM)
       PXZ = DABS(1./PEST(1))
       
       PX = PXZ*DSIN(PHI)
       PY = PXZ*DTAN(ALAM)
       PZ = PXZ*DCOS(PHI)
       PTOT = PXZ/DCOS(ALAM)
       
       VECT(1) = PEST(4)
       VECT(2) = PEST(5)
       VECT(3) = 0.
       VECT(4) = PX/PTOT
       VECT(5) = PY/PTOT
       VECT(6) = PZ/PTOT
       VECT(7) = PTOT
              
       Z = VECT(3)
       NSTEP = 0
*
** Runge Kutta  
**      PRINT *,' AV GRKUTA ASIGN',ASIGN,' THET',THET
       ISTOLD = 0 
       DO ICH = 1,NPLANE

          IST = INT(FLOAT(ICH+1)/2.)


          DO WHILE (Z.GE.0..AND.Z.LT.ABS(ZPLANEP(ICH))
     &         .AND.NSTEP.LE.NSTEPMAX)
             NSTEP = NSTEP+1 
**          WRITE(6,*) NSTEP,(VECT(I),I=1,7)
**             CALL RECO_GRKUTA (ASIGN,STEP,VECT,VOUT) ! CCC
             CALL RECO_GHELIX (ASIGN,STEP,VECT,VOUT) 
             DO I = 1,7
                VECT(I) = VOUT(I)
             ENDDO   
             Z = VECT(3)
          ENDDO
          IF (IST.NE.ISTOLD) THEN
             IPCH = 1
          ELSE
             IPCH = 2
          ENDIF   
          XPL(IST,IPCH) = VECT(1)-(Z-ABS(ZPLANEP(ICH)))*VECT(4)/VECT(6)
          YPL(IST,IPCH) = VECT(2)-(Z-ABS(ZPLANEP(ICH)))*VECT(5)/VECT(6)
          IF (IPCH.EQ.2) THEN
             DX = XPL(IST,2)-XPL(IST,1)   
             DY = YPL(IST,2)-YPL(IST,1)
             PHPL(IST) = -1.*DATAN2(DX,DZ_PL(IST))
             ALPL(IST) = DATAN2(DY,DSQRT(DX**2+DZ_PL(IST)**2))
          ENDIF   
          ISTOLD = IST
       ENDDO   
**          print *,' vect= ',vect(1),vect(2),vect(3),vect(4),vect(5),
**     &        vect(6),vect(7)

        
       RETURN
       END
*******************************************************************************
       SUBROUTINE FCN(NPAR,GRAD,FVAL,PEST,IFLAG,FUTIL)
*******************************************************************************
*   Calculate FVAL=CHI2 the function minimized by minuit for a given track
*
*******************************************************************************
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)

       PARAMETER(NBSTATION=5)

*       DIMENSION PEST(*),GRAD(*)
       DIMENSION PEST(5),GRAD(5)
       DIMENSION PEEST(NBSTATION)

       COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
       
     
       COMMON/PRECIS/EEXM(NBSTATION),EEYM(NBSTATION),EEPH(NBSTATION),
     &     EEAL(NBSTATION)

       COMMON /MEASUR/XMES(NBSTATION),YMES(NBSTATION),IZMES(NBSTATION),
     &               PHMES(NBSTATION),ALMES(NBSTATION), MPOS(NBSTATION),
     &               MANG(NBSTATION)

       COMMON /PLANE/XPL(NBSTATION,2),YPL(NBSTATION,2),PHPL(NBSTATION),
     &               ALPL(NBSTATION),CHI2PL

       COMMON/VERTEX/ERRV,IVERTEX

       DIMENSION XC(NBSTATION),YC(NBSTATION)

       EXTERNAL RECOCHI2

       PEEST(1)=PEST(1)
       PEEST(2)=PEST(2)
       PEEST(3)=PEST(3)
       IF(IVERTEX.EQ.1) THEN
          PEEST(4)=PEST(4)	! position du vertex
          PEEST(5)=PEST(5)
       ELSE
          PEEST(4)=0.0D0
          PEEST(5)=0.0D0
       ENDIF

       ALAM = DATAN(PEST(3))
       PXZ = DABS(1./PEST(1))
       PTOT = PXZ/DCOS(ALAM)
       

       CALL FOLLOW (0.0D0,PEEST,IFLAG,XPL,YPL,PHPL,ALPL) ! calcul des 
       IF(IFLAG.EQ.1) THEN                       ! points d impacts dans les
          PRINT *,'FCN ',XPL(4,1),XMES(4)        ! plans
          PRINT *,'FCN ',YPL(4,1),YMES(4)
          PRINT *,'FCN ',XPL(5,1),XMES(5)
          PRINT *,'FCN ',YPL(5,1),YMES(5)
       ENDIF

       DO I = 1,NBSTATION
          XC(I) = XPL(I,1)
          YC(I) = YPL(I,1)
       ENDDO

       IF (IVERTEX.EQ.1) THEN 

          FVAL = RECOCHI2(MPOS,MANG,XMES,YMES,ALMES,PHMES,
     &     XC,YC,ALPL,PHPL,PTOT,IZMES,NPLPL)

       ELSE

          FVAL = 0.0D0

          NPLPL = 0 
          DO J = 1,NBSTATION           
             IF (MPOS(J).EQ.1) THEN
                NPLPL = NPLPL+1
                XPLC = XPL(J,IZMES(J)) 
                YPLC = YPL(J,IZMES(J))
                FF = (XMES(J) - XPLC)/EEXM(J)
                FVAL = FVAL + FF**2
                FF = (YMES(J) - YPLC)/EEYM(J)
                FVAL = FVAL + FF**2
             ENDIF
             IF (MANG(J).EQ.1) THEN
                NPLPL = NPLPL+1
                FF = (PHMES(J) - PHPL(J))/EEPH(J)
                FVAL = FVAL + FF**2
                FF = (ALMES(J) - ALPL(J))/EEAL(J)
                FVAL = FVAL + FF**2
             ENDIF
          ENDDO

       ENDIF   

       NPARAM = 3
       IF (IVERTEX.EQ.1) NPARAM = 5
       CHI2PL = FVAL/FLOAT(2*NPLPL-NPARAM)       

       RETURN
       END
********************************************************************************
      SUBROUTINE STOCK_CANDIDAT(ICH1,JHITCH1,ICH2,IFIND,IFIND2,EXM,EYM,
     &     EPH,EAL,NCAN,ICODE)
********************************************************************************
*   Fill common CANDIDAT with track candidates from the search in stations 4&5
*
********************************************************************************
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      PARAMETER (MAXHITTOT=20000,MAXHITCH=10000,NBSTATION=5,MAXCAN=1000)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 
      
      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

            
      COMMON/CANDIDAT/JCAN(NBSTATION,MAXCAN),JCANTYP(MAXCAN),
     &          EEX(MAXCAN),EEY(MAXCAN),EEP(MAXCAN),EEA(MAXCAN)
      DIMENSION IFIND2(10)
      
**      PRINT *,'STOCK st. init.=',ICH1,'id. init.=',ID(IP(ICH1,JHITCH1))
**      PRINT *,'STOCK st. finale=',ICH2,'id. final=',ID(IP(ICH2,IFIND))
**      PRINT *,'STOCK ifind=',IFIND 
**      PRINT *,'STOCK icode=',ICODE 

      DO I = 1,10 
         IF (IFIND2(I).GT.0) THEN
            NCAN = NCAN+1
            JCAN(ICH1,NCAN) = JHITCH1
            JCAN(ICH2,NCAN) = IFIND2(I)
            JCANTYP(NCAN) = ICODE 
            EEX(NCAN) = EXM
            EEY(NCAN) = EYM
            EEP(NCAN) = EPH
            EEA(NCAN) = EAL
         ENDIF   
      ENDDO
   
      RETURN
      END
*******************************************************************************
       SUBROUTINE FCNOLD(NPAR,GRAD,FVAL,PEST,IFLAG,FUTIL)
*******************************************************************************
*   Calculate FVAL=CHI2 the function minimized by minuit for a given track
*
*******************************************************************************
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)

       PARAMETER(NBSTATION=5)

*       DIMENSION PEST(*),GRAD(*)
       DIMENSION PEST(5),GRAD(5)
       DIMENSION PEEST(NBSTATION)

       COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)
       
     
       COMMON/PRECIS/EEXM(NBSTATION),EEYM(NBSTATION),EEPH(NBSTATION),
     &     EEAL(NBSTATION)

       COMMON /MEASUR/XMES(NBSTATION),YMES(NBSTATION),IZMES(NBSTATION),
     &               PHMES(NBSTATION),ALMES(NBSTATION), MPOS(NBSTATION),
     &               MANG(NBSTATION)

       COMMON /PLANE/XPL(NBSTATION,2),YPL(NBSTATION,2),PHPL(NBSTATION),
     &               ALPL(NBSTATION),CHI2PL

       COMMON/VERTEX/ERRV,IVERTEX

       PEEST(1)=PEST(1)
       PEEST(2)=PEST(2)
       PEEST(3)=PEST(3)
       IF(IVERTEX.EQ.1) THEN
          PEEST(4)=PEST(4)	! position du vertex
          PEEST(5)=PEST(5)
       ELSE
          PEEST(4)=0.0D0
          PEEST(5)=0.0D0
       ENDIF

       CALL FOLLOW (0.0D0,PEEST,IFLAG,XPL,YPL,PHPL,ALPL) ! calcul des 
       IF(IFLAG.EQ.1) THEN                       ! points d impacts dans les
          PRINT *,'FCN ',XPL(4,1),XMES(4)        ! plans
          PRINT *,'FCN ',YPL(4,1),YMES(4)
          PRINT *,'FCN ',XPL(5,1),XMES(5)
          PRINT *,'FCN ',YPL(5,1),YMES(5)
       ENDIF

*       IF (IVERTEX.EQ.1) THEN
*          FVAL = (PEST(4)/ERRV)**2 + (PEST(5)/ERRV)**2
*       ELSE
          FVAL = 0.0D0
*       ENDIF
       NPLPL = 0 
       DO J = 1,NBSTATION           
          IF (MPOS(J).EQ.1) THEN
             NPLPL = NPLPL+1
             XPLC = XPL(J,IZMES(J)) 
             YPLC = YPL(J,IZMES(J))
             FF = (XMES(J) - XPLC)/EEXM(J)
             FVAL =FVAL + FF**2
             FF = (YMES(J) - YPLC)/EEYM(J)
             FVAL =FVAL + FF**2
          ENDIF
          IF (MANG(J).EQ.1) THEN
             NPLPL = NPLPL+1
             FF = (PHMES(J) - PHPL(J))/EEPH(J)
             FVAL =FVAL + FF**2
             FF = (ALMES(J) - ALPL(J))/EEAL(J)
             FVAL =FVAL + FF**2
          ENDIF
       ENDDO
**       PRINT *,'ST 1',XPL(1,1),XMES(1),YPL(1,1),YMES(1)          
**       PRINT *,'ST 2',XPL(2,1),XMES(2),YPL(2,1),YMES(2)          
**       PRINT *,'ST 3',XPL(3,1),XMES(3),YPL(3,1),YMES(3)          
**       PRINT *,'ST 4',XPL(4,1),XMES(4),YPL(4,1),YMES(4)          
**       PRINT *,'ST 5',XPL(5,1),XMES(5),YPL(5,1),YMES(5)          
       NPARAM = 3
       IF (IVERTEX.EQ.1) NPARAM = 5
       CHI2PL = FVAL/FLOAT(2*NPLPL-NPARAM)       

       RETURN
       END
****************************************************************************
      SUBROUTINE CHECK_HISTO4(IDHIST,ICH2,IHIT2,ICH1,IHIT1,
     &                        X1,Y1,PHI1,ALAM1,P1,EXM,EYM,EPH,EAL)   
*****************************************************************************
*   Check hit IHIT2 with GEANT informations from hit HIT1
*
*    INPUT : ICH2 : No st. de recherche
*            IDCH1,X1,Y1,PHI1,ALAM1,P1 : Trace de reference
*    OUTPUT : JOK : No hit dans ICH2 appartenant a la meme trace.
* 
*****************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      PARAMETER (MAXHITCH=10000,MAXHITTOT=20000,NBSTATION=5,
     &           MAXTRK=50000)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 

      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT

      COMMON/VERIFGEANT/ITTROUGH(MAXTRK,NBSTATION),
     &                  IT_LIST(MAXTRK),IT_NP(MAXTRK),ITCHECK(MAXTRK),
     &                  ITRACK(MAXHITTOT)
* 
      COMMON/DEBEVT/IDEBUG

      REAL*4 R2
      DATA R2/1./

      JOK = 0
       
      DO I=1,JHIT(ICH2)
         IF (PHM(ICH2,I).LE.6.3) THEN ! vector measurement
            IF (ID(IP(ICH1,IHIT1)).EQ.ID(IP(ICH2,I))) THEN
               JOK = I 
               IF (IDHIST.GT.0) THEN
*                  CALL CHF1(IDHIST,SNGL(P1),SNGL((X1-XM(ICH2,I))**2))
*                  CALL CHF1(IDHIST+1,SNGL(P1),SNGL((Y1-YM(ICH2,I))**2))
*                  CALL CHF1(IDHIST+2,SNGL(P1),
*     &                     SNGL((PHI1-PHM(ICH2,I))**2))
*                  CALL CHF1(IDHIST+3,SNGL(P1),
*     &                     SNGL((ALAM1-ALM(ICH2,I))**2))
*                  CALL CHF1(IDHIST+4,SNGL(P1),R2)
               ENDIF   
            ENDIF
         ENDIF
      ENDDO
      
      IF (JOK.GT.0) THEN
         IF (ITCHECK(ITRACK(IP(ICH1,IHIT1))).EQ.1) THEN
            IF (IDEBUG.GE.2) THEN 
               IF (IHIT2.EQ.0) THEN
                  PRINT *,'CHECK4 histo nb:',IDHIST 
                  PRINT *,'CHECK4 p de st.',ICH1,'=',P1
                  PRINT *,'CHECK4 track not found in st.',ICH2
                  PRINT *,'CHECK4 error X :',(XM(ICH2,JOK)-X1), EXM
                  PRINT *,'CHECK4 error Y :',(YM(ICH2,JOK)-Y1), EYM
                  PRINT *,'CHECK4 error PHI :',(PHM(ICH2,JOK)-PHI1),EPH
                  PRINT *,'CHECK4 error ALAM :',(ALM(ICH2,JOK)-ALAM1),
     &                     EAL
               ELSEIF (IHIT2.NE.JOK) THEN
                  PRINT *,'CHECK4 histo nb:',IDHIST 
                  PRINT *,'CHECK4 p de st.',ICH1,'=',P1
                  PRINT *,'CHECK4 ghost in st.',ICH2
                  PRINT *,'CHECK4 id part. recherchee:',
     &                    ID(IP(ICH1,IHIT1))
                  PRINT *,'CHECK4 id ghost trouve    :',
     &                    ID(IP(ICH2,IHIT2))  
                  PRINT *,'CHECK4 JOK=',JOK,' IHIT2=',IHIT2
               ENDIF
            ENDIF
         ENDIF
      ENDIF
               
      RETURN
      END     
*****************************************************************************
      SUBROUTINE CHECK_HISTO2(IDHIST,ICH2,IHIT2,ICH1,IHIT1,
     &                        X1,Y1,X2,Y2,P1,EX1,EY1,EX2,EY2)   
*****************************************************************************
*   Check hit IHIT2 with GEANT informations from hit HIT1
*
*    INPUT : IDHIST : No histo
*            ICH2 : No st. de recherche
*            IDCH1,X1,Y1,PHI1,ALAM1,P1 : Trace de reference
*    OUTPUT : JOK : No hit dans ICH2 appartenant a la meme trace.
* 
*****************************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      PARAMETER (MAXHITCH=10000,MAXHITTOT=20000,NBSTATION=5,
     &           MAXTRK=50000)
      
      COMMON/CHHIT/XM(NBSTATION,MAXHITCH),YM(NBSTATION,MAXHITCH),
     &             PHM(NBSTATION,MAXHITCH),ALM(NBSTATION,MAXHITCH),
     &             IZM(NBSTATION,MAXHITCH),
     &             IP(NBSTATION,MAXHITCH),JHIT(NBSTATION),
     &             XMR(NBSTATION,MAXHITCH,2),YMR(NBSTATION,MAXHITCH,2) 

      COMMON/RHIT/ITYP(MAXHITTOT),XTR(MAXHITTOT),YTR(MAXHITTOT),
     &            PTOT(MAXHITTOT),ID(MAXHITTOT),IZST(MAXHITTOT),
     &            PVERT1(MAXHITTOT),PVERT2(MAXHITTOT),PVERT3(MAXHITTOT),
     &            ZVERT(MAXHITTOT),NHITTOT
     
      COMMON/VERIFGEANT/ITTROUGH(MAXTRK,NBSTATION),
     &                  IT_LIST(MAXTRK),IT_NP(MAXTRK),ITCHECK(MAXTRK),
     &                  ITRACK(MAXHITTOT)
* 
      COMMON/DEBEVT/IDEBUG
 
      REAL*4 R2
      DATA R2/1./

      JOK = 0
       
      DO I=1,JHIT(ICH2)
         IF (ID(IP(ICH1,IHIT1)).EQ.ID(IP(ICH2,I))) THEN
            JOK = I 
            IF (IDHIST.GT.0) THEN
               IF (IZM(ICH2,I).EQ.1) THEN ! 1st chamber
                  X = X1
                  Y = Y1
               ELSE  ! 2nd chamber
                  X = X2
                  Y = Y2
               ENDIF      
               CALL CHF1(IDHIST,SNGL(P1),SNGL((X-XM(ICH2,I))**2))
               CALL CHF1(IDHIST+1,SNGL(P1),SNGL((Y-YM(ICH2,I))**2))
               CALL CHF1(IDHIST+4,SNGL(P1),R2)
            ENDIF   
         ENDIF
      ENDDO
      
      IF (JOK.GT.0) THEN
         IF (ITCHECK(ITRACK(IP(ICH1,IHIT1))).EQ.1) THEN
            EXM = EX1
            EYM = EY1 
            IF (IZM(ICH2,JOK).EQ.1) THEN
               X = X1
               Y = Y1
               IF (ICH2.EQ.4.OR.ICH2.EQ.5) THEN
                  EXM = EX1
                  EYM = EY1
               ENDIF   
            ELSE
               X = X2
               Y = Y2
               IF (ICH2.EQ.4.OR.ICH2.EQ.5) THEN
                  EXM = EX2
                  EYM = EY2
               ENDIF   
            ENDIF      
            IF (IDEBUG.GE.2) THEN 
               IF (IHIT2.EQ.0) THEN
                  PRINT *,'CHECK2 histo nb:',IDHIST 
                  PRINT *,'CHECK2 p de st.',ICH1,'=',P1
                  PRINT *,'CHECK2 track not found in st.',ICH2
                  PRINT *,'CHECK2 error X :',(XM(ICH2,JOK)-X), EXM
                  PRINT *,'CHECK2 error Y :',(YM(ICH2,JOK)-Y), EYM
               ELSEIF(IHIT2.NE.JOK) THEN
                  PRINT *,'CHECK2 histo nb:',IDHIST 
                  PRINT *,'CHECK2 p de st.',ICH1,'=',P1
                  PRINT *,'CHECK2 ghost in st.',ICH2
                  PRINT *,'CHECK2 id part. recherchee:',
     &                    ID(IP(ICH1,IHIT1))
                  PRINT *,'CHECK2 id ghost trouve    :',
     &                    ID(IP(ICH2,IHIT2))  
                  PRINT *,'CHECK2 JOK=',JOK,' IHIT2=',IHIT2  
               ENDIF
            ENDIF
         ENDIF
      ENDIF
               
      RETURN
      END

      DOUBLE PRECISION FUNCTION DEDX(P,THET,XEA,YEA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REA = DSQRT(XEA**2+YEA**2)
      IF (REA.lT.26.3611) THEN
         if (p .lt. 15.) then
            DP=2.737+0.0494*p-0.001123*p*p
         else 
	    DP=3.0643+0.01346*p
         endif
      ELSE
         if (p .lt. 15.) then
	    DP = 2.1380+0.0351*p-0.000853*p*p
         else 
	    DP = 2.407+0.00702*p
         endif
      ENDIF
      P=P+DP/DCOS(THET)
      DEDX=P
      RETURN
      END

************************************************************************        
      DOUBLE PRECISION FUNCTION DEDX_oldold(P,THET,XEA,YEA)
************************************************************************        
*    DEDX est la nouvelle impulsion au vertex, corrigee de la perte    
*    d'energie dans l'absorbeur
*
************************************************************************        
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PPB(6), PW(6)
*   FIT RESULT FOR PB REGION (5TH ORDER POLY)
*   1  p0           2.24358e+03   2.55765e+00   1.06982e-03   6.30474e-07
*   2  p1           1.16393e+01   4.45081e-02   6.58627e-06  -1.15822e-04
*   3  p2          -1.82314e-01   3.28429e-04   8.69340e-08  -1.27281e-02
*   4  p3           1.60930e-03   2.23812e-06   7.67374e-10  -4.15573e-01
*   5  p4          -6.96885e-06   1.35405e-08   3.32301e-12   1.22136e+02
*   6  p5           1.16339e-08   5.91665e-11   1.06532e-14   3.33965e+04     

      DATA PPB /2.24358d+03, 1.16393d+01, -1.82314d-01, 1.60930d-03
     +         ,-6.96885d-06, 1.16339d-08/
*   FIT RESULT FOR W REGION (5TH ORDER POLY)
*   1  p0           2.90155e+03   3.49066e+00   1.38357e-03  -2.79916e-05
*   2  p1           1.57716e+01   6.09946e-02   9.03687e-06  -6.63098e-03
*   3  p2          -2.48349e-01   4.50365e-04   1.18422e-07  -8.27199e-01
*   4  p3           2.19908e-03   3.07148e-06   1.04860e-09  -1.03290e+02
*   5  p4          -9.54046e-06   1.85908e-08   4.54924e-12  -1.51284e+04
*   6  p5           1.59463e-08   8.11346e-11   1.46446e-14  -2.69491e+06     
      DATA PW /2.90155d+03, 1.57716d+01, -2.48349d-01, 2.19908d-03 
     +         , -9.54046d-06, 1.59463d-08/
* 
      REA = DSQRT(XEA**2+YEA**2)
      IF (REA.GT.26.3611) THEN
         DP=PPB(1)+PPB(2)*P+PPB(3)*P**2
     &    +PPB(4)*P**3+PPB(5)*P**4+PPB(6)*P**5
      ELSE
         DP=PW(1)+PW(2)*P+PW(3)*P**2
     &    +PW(4)*P**3+PW(5)*P**4+PW(6)*P**5
      ENDIF
      P=P+DP/1000.D0/DCOS(THET)
      DEDX=P
      RETURN
      END
************************************************************************        
      DOUBLE PRECISION FUNCTION DEDX_OLD(P,THET,XEA,YEA)
************************************************************************        
*    DEDX est la nouvelle impulsion au vertex, corrigee de la perte    
*    d'energie dans l'absorbeur
*
************************************************************************        
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      REA = DSQRT(XEA**2+YEA**2)
      IF (REA.GT.26.3611) THEN
*   Plomb      
         SPB = 5./DCOS(THET)
         P = P+SPB/1000.*11.35*(16.66D-3 * P+1.33) 
*   Polyethylene      
         SPO = 5./DCOS(THET)
         P = P+SPO/1000.*0.935*(2.22D-3 * P+2.17) 
*   Plomb      
         SPB = 5./DCOS(THET)
         P = P+SPB/1000.*11.35*(16.66D-3 * P+1.33) 
*   Polyethylene      
         SPO = 5./DCOS(THET)
         P = P+SPO/1000.*0.935*(2.22D-3 * P+2.17) 
*   Plomb      
         SPB = 5./DCOS(THET)
         P = P+SPB/1000.*11.35*(16.66D-3 * P+1.33) 
*   Polyethylene      
         SPO = 5./DCOS(THET)
         P = P+SPO/1000.*0.935*(2.22D-3 * P+2.17) 
*   Plomb      
         SPB = 5./DCOS(THET)
         P = P+SPB/1000.*11.35*(16.66D-3 * P+1.33) 

      ELSE   
*   Tungstene     
         SW = (503.-467.)/DCOS(THET) 
         P = P+SW/1000.*19.3*(16.66D-3 * P+1.33)
      ENDIF

*   Concrete      
      SCONC = (467.-315.)/DCOS(THET)
      P = P+SCONC/1000.*2.5*(2.22D-3*P+2.17)
      
*   Carbone      
      SC = (315.-90.)/DCOS(THET)
      P = P+SC/1000.*1.93*(2.22D-3*P+2.17) ! Carbone
               
      DEDX = P  
      
      RETURN
 
      END
*/

************************************************************************        
      SUBROUTINE BRANSON(PXZ,PHI,ALAM,XEA,YEA)
************************************************************************        
*
*   Correction de Branson du multiple scattering dans l'absorbeur
*
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER(NPLANE=10) 
      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
     
      ASIGN = 1.      
      IF (PXZ.LT.0.) ASIGN = -1.       
      PXZ = DABS(PXZ)
      PX = PXZ*DSIN(PHI)
      PY = PXZ*DTAN(ALAM)
      PZ = PXZ*DCOS(PHI)
      
      PTOT = PXZ/DCOS(ALAM)
      
      ZEA = ZABS
      
      REA = DSQRT(XEA**2+YEA**2)
      IF (REA.GT.26.3611) THEN 
         ZBP = ZBP1
      ELSE            ! Abso. W  pour theta < 3 deg
         ZBP = ZBP2
      ENDIF
*      ZBP = ZEA ! Andreas
      XBP = XEA-PX/PZ*(ZEA-ZBP)
      YBP = YEA-PY/PZ*(ZEA-ZBP)
      PZ = PTOT*ZBP/DSQRT(XBP**2+YBP**2+ZBP**2)
      PX = PZ*XBP/ZBP
      PY = PZ*YBP/ZBP
      PXZ = DSQRT(PX**2+PZ**2)
      PHI = DATAN2(PX,PZ)
      ALAM = DATAN2(PY,PXZ)

**      THET = DATAN2(REA,ZEA)

      PT = DSQRT(PX**2+PY**2)      
      THET = DATAN2(PT,PZ) 
      PTOT =  DEDX(PTOT,THET,XEA,YEA)
      
      PXZ = ASIGN*PTOT*DCOS(ALAM)
      
      RETURN
      END

***************************************************************
      SUBROUTINE DSINV(N,A,IDIM,IFAIL)
***************************************************************

      DOUBLE PRECISION    A(IDIM,*),  ZERO,  ONE,  X, Y

      REAL                PIVOTF
      CHARACTER*6         HNAME

      DOUBLE PRECISION    S1, S31, S32, S33,  DOTF

      PIVOTF(X)    =  SNGL(X)
      DOTF(X,Y,S1)  =  X * Y + S1

      DATA      HNAME               /  'DSINV '  /
      DATA      ZERO, ONE           /  0.D0, 1.D0 /

      IF(IDIM .LT. N  .OR.  N .LE. 0)  GOTO 900
*     
* sfact.inc
*
      IFAIL  =  0
      DO 144    J  =  1, N
         IF(PIVOTF(A(J,J)) .LE. 0.)  GOTO 150
         A(J,J)  =  ONE / A(J,J)
         IF(J .EQ. N)  GOTO 199
 140     JP1  =  J+1
         DO 143   L  =  JP1, N
            A(J,L)  =  A(J,J)*A(L,J)
            S1      =  -A(L,J+1)
            DO 141  I  =  1, J
               S1  =  DOTF(A(L,I),A(I,J+1),S1)
 141        CONTINUE
            A(L,J+1)  =  -S1
 143     CONTINUE
 144  CONTINUE
 150  IFAIL  =  -1
      RETURN
 199  CONTINUE
*     
* sfinv.inc
*     
      IF(N .EQ. 1)  GOTO 399
      A(1,2)  =  -A(1,2)
      A(2,1)  =   A(1,2)*A(2,2)
      IF(N .EQ. 2)  GOTO 320
      DO 314    J  =  3, N
         JM2  =  J - 2
         DO 312 K  =  1, JM2
            S31  =  A(K,J)
            DO 311  I  =  K, JM2
               S31  =  DOTF(A(K,I+1),A(I+1,J),S31)
 311        CONTINUE
            A(K,J)  =  -S31
            A(J,K)  =  -S31*A(J,J)
 312     CONTINUE
         A(J-1,J)  =  -A(J-1,J)
         A(J,J-1)  =   A(J-1,J)*A(J,J)
 314  CONTINUE
 320  J  =  1
 323  S33  =  A(J,J)
      IF(J .EQ. N)  GOTO 325
      JP1  =  J + 1
      DO 324 I  =  JP1, N
         S33  =  DOTF(A(J,I),A(I,J),S33)
 324  CONTINUE
 325  A(J,J)  =  S33
      JM1  =  J
      J    =  JP1
      DO 328 K  =  1, JM1
         S32  =  ZERO
         DO 327  I  =  J, N
            S32  =  DOTF(A(K,I),A(I,J),S32)
 327     CONTINUE
         A(K,J)  =  S32
         A(J,K)  =  S32
 328  CONTINUE
      IF(J .LT. N)  GOTO 323
 399  CONTINUE

      RETURN
 900  CALL TMPRNT(HNAME,N,IDIM,0)
      RETURN
      END

*******************************************************
      SUBROUTINE TMPRNT(NAME,N,IDIM,K)
*******************************************************

      CHARACTER*6         NAME
      LOGICAL             MFLAG,    RFLAG

      IF(NAME(2:2) .EQ. 'S') THEN
         CALL KERMTR('F012.1',LGFILE,MFLAG,RFLAG)
      ELSE
         CALL KERMTR('F011.1',LGFILE,MFLAG,RFLAG)
      ENDIF
      IF(NAME(3:6) .EQ. 'FEQN') ASSIGN 1002 TO IFMT
      IF(NAME(3:6) .NE. 'FEQN') ASSIGN 1001 TO IFMT
      IF(MFLAG) THEN
         IF(LGFILE .EQ. 0) THEN
            IF(NAME(3:6) .EQ. 'FEQN') THEN
               WRITE(*,IFMT) NAME, N, IDIM, K
            ELSE
               WRITE(*,IFMT) NAME, N, IDIM
            ENDIF
         ELSE
            IF(NAME(3:6) .EQ. 'FEQN') THEN
               WRITE(LGFILE,IFMT) NAME, N, IDIM, K
            ELSE
               WRITE(LGFILE,IFMT) NAME, N, IDIM
            ENDIF
         ENDIF
      ENDIF
      IF(.NOT. RFLAG) CALL ABEND
      RETURN
 1001 FORMAT(7X, 31H PARAMETER ERROR IN SUBROUTINE , A6,
     +     27H ... (N.LT.1 OR IDIM.LT.N).,
     +             5X, 3HN =, I4, 5X, 6HIDIM =, I4, 1H. )
 1002 FORMAT(7X, 31H PARAMETER ERROR IN SUBROUTINE , A6,
     +     37H ... (N.LT.1 OR IDIM.LT.N OR K.LT.1).,
     +     5X, 3HN =, I4, 5X, 6HIDIM =, I4, 5X, 3HK =, I4,1H.)
      END
*
* $Id$
*
* $Log$
* Revision 1.4.4.2  2000/04/26 15:48:37  morsch
* Some routines from obsolete algo.F are needed by reco_muon.F and have been
* copied there.
*
* Revision 1.4.4.1  2000/01/12 16:00:55  morsch
* New version of MUON code
*
* Revision 1.1.1.1  1996/02/15 17:48:35  mclareni
* Kernlib
*
*

***********************************************************
      SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
***********************************************************

      PARAMETER(KOUNTE  =  27)
      CHARACTER*6         ERCODE,   CODE(KOUNTE)
      LOGICAL             MFLAG,    RFLAG
      INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)

      DATA      LOGF      /  0  /
      DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 255, 255 /
      DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 255, 255 /
      DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 255, 255 /
      DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 255, 255 /
      DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 255, 255 /
      DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 255, 255 /
      DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 255, 255 /
      DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 255, 255 /
      DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 255, 255 /
      DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 255, 255 /
      DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 255, 255 /
      DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 255, 255 /
      DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 255, 255 /
      DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 255, 255 /
      DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 255, 255 /
      DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 255, 255 /
      DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 255, 255 /
      DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 255, 255 /
      DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 255, 255 /
      DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 255, 255 /
      DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 255, 255 /
      DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 255,   0 /
      DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 255,   0 /
      DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 255,   0 /
      DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 255,   0 /
      DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 255, 255 /
      DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 255, 255 /

      LOGF  =  LGFILE
      L  =  0
      IF(ERCODE .NE. ' ')  THEN
         DO 10  L = 1, 6
            IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
 10      CONTINUE
 12      CONTINUE
      ENDIF
      DO 14     I  =  1, KOUNTE
         IF(L .EQ. 0)  GOTO 13
         IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
 13      IF(LIMITM.GE.0) KNTM(I)  =  LIMITM
         IF(LIMITR.GE.0) KNTR(I)  =  LIMITR
 14   CONTINUE
      RETURN
      ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
      LOG  =  LOGF
      DO 20     I  =  1, KOUNTE
         IF(ERCODE .EQ. CODE(I))  GOTO 21
 20   CONTINUE
      WRITE(*,1000)  ERCODE
      CALL ABEND
      RETURN
 21   RFLAG  =  KNTR(I) .GE. 1
      IF(RFLAG  .AND.  (KNTR(I) .LT. 255))  KNTR(I)  =  KNTR(I) - 1
      MFLAG  =  KNTM(I) .GE. 1
      IF(MFLAG  .AND.  (KNTM(I) .LT. 255))  KNTM(I)  =  KNTM(I) - 1
      IF(.NOT. RFLAG)  THEN
         IF(LOGF .LT. 1)  THEN
            WRITE(*,1001)  CODE(I)
         ELSE
            WRITE(LOGF,1001)  CODE(I)
         ENDIF
      ENDIF
      IF(MFLAG .AND. RFLAG)  THEN
         IF(LOGF .LT. 1)  THEN
            WRITE(*,1002)  CODE(I)
         ELSE
            WRITE(LOGF,1002)  CODE(I)
         ENDIF
      ENDIF
      RETURN
 1000 FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +     ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +     ' ERROR MONITOR. RUN ABORTED.')
 1001 FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +     'CONDITION ',A6)
 1002 FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
      END

**************************
      subroutine abend
**************************

      stop 'abend!'
      end

************************************************************************        
      SUBROUTINE FCNFIT(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)
************************************************************************        
*    With magnetic Field Map GRKUTA
*      
*    Calcule FVAL: la fonction minimisee par MINUIT
*    With magnetic field map
*      
************************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(NPLANE=10)
            
      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
     
      COMMON/MEAS/LPLANE(NPLANE),XMP(NPLANE),YMP(NPLANE)
     
      COMMON/FCNOUT/PXZEA,ALAMEA,PHIEA,XEA,YEA,NPLU,CHI2
                 

      DIMENSION GRAD(*),XVAL(*),AMS(500),DISTAZ(500)
  
      DIMENSION XP(NPLANE),YP(NPLANE),
     &          COV(NPLANE,NPLANE),AP(NPLANE),COVY(NPLANE,NPLANE)
      DIMENSION VECT(7),VOUT(7) 
      
      STEP = 2.   ! 1 cm
      NSTEPMAX = 5000
      PITODEG = 57.295
      XV = XVAL(4)   
      YV = XVAL(5)   
              
      ASIGN = 1.
      IF (XVAL(1).LT.0.) ASIGN = -1.
      PHI  = XVAL(2)
      ALAM = XVAL(3)
      PXZ = DABS(1./XVAL(1))
           
      PX = PXZ*DSIN(PHI)
      PY = PXZ*DTAN(ALAM)
      PZ = PXZ*DCOS(PHI)
      PTOT = PXZ/DCOS(ALAM)
      
      A12 = 0.
      NTMAX = 0
      
      ZEA = ZABS
      XEA = XV
      YEA = YV
      PXEA = PX
      PYEA = PY
      PHIEA = PHI
      PXZEA = ASIGN*PXZ
      ALAMEA = ALAM

      VECT(1) = XV
      VECT(2) = YV
      VECT(3) = ZABS
      VECT(4) = PX/PTOT
      VECT(5) = PY/PTOT
      VECT(6) = PZ/PTOT
      VECT(7) = PTOT
      
      R = SQRT(VECT(1)*VECT(1)+VECT(2)*VECT(2))

      Z = VECT(3)
      NSTEP = 0
      IX = 0
      IY = 0 
      IZ = 0
**      PRINT *,' AV GRKUTA ASIGN',ASIGN,' THET',THET
      DO ICH = 1,NPLANE 
         DO WHILE (Z.GE.ZABS.AND.Z.LT.ZPLANEP(ICH)
     &            .AND.NSTEP.LE.NSTEPMAX)
**     &            .AND.(THETA*PITODEG).GT.2.
**     &            .AND. (THETA*PITODEG).LT.9.) 
            NSTEP = NSTEP+1 
**          WRITE(6,*) NSTEP,(VECT(I),I=1,7)
**            CALL RECO_GRKUTA(ASIGN,STEP,VECT,VOUT) ! CCC
            CALL RECO_GHELIX(ASIGN,STEP,VECT,VOUT)
            DO I = 1,7
               VECT(I) = VOUT(I)
            ENDDO   
            Z = VECT(3)
            R = SQRT(VECT(1)*VECT(1)+VECT(2)*VECT(2))
         ENDDO
         IF (NSTEP.EQ.NSTEPMAX) RETURN
         XP(ICH) = VECT(1)-(Z-ZPLANEP(ICH))*VECT(4)/VECT(6)
         YP(ICH) = VECT(2)-(Z-ZPLANEP(ICH))*VECT(5)/VECT(6)
         AL   = THICK/ VECT(6)
         AP(ICH) = (0.0136D0/PTOT)*DSQRT(AL)*(1+0.038D0*DLOG(AL))
      ENDDO
**    PRINT *,' AP GRKUTA ASIGN',ASIGN,' THET',THET


** Matrice de covariance      
      I = 0
      DO II = 1,NPLANE
        IF (LPLANE(II).EQ.1) THEN
           I = I + 1
*     I = II
           J = I - 1
           DO JJ = II, NPLANE
              IF (LPLANE(JJ).EQ.1) THEN
                 J = J + 1
*     J = JJ
                 COV (I,J) = 0.0D0
                 COV (J,I) = A12
                 IF (I .EQ. J) THEN
                    COV(J,I) =COV(J,I) + XPREC**2
                 ENDIF      
                 
*     IF (I .EQ. 10 .AND. J .EQ. 10) PRINT *,'10 10   ',COV(J,I)
**                 DO L = 1,NTMAX
**                    COV(J,I) = COV(J,I)
**     &                   +  (ZPLANEP(II) + DISTAZ(L))*(ZPLANEP(JJ) + 
**     &                   DISTAZ(L))*AMS(L)**2
**                 ENDDO
                 DO K = 1, II-1
                    COV(J,I) = COV(J,I)
     &                   +  (ZPLANEP(II)-ZPLANEP(K))*
     &                   (ZPLANEP(JJ)-ZPLANEP(K))*AP(K)**2
*     IF (I .EQ. 10 .AND. J .EQ. 10) PRINT *,'10 10   ',COV(J,I)
                 ENDDO
                 COVY(I,J) = 0.0D0
                 COVY(J,I) = COV(J,I)
                 IF (I .EQ. J) THEN
                    COVY(J,I) = COVY(J,I) - XPREC**2 + YPREC**2
                 ENDIF
              ENDIF   
           ENDDO
        ENDIF
      ENDDO
 
*  Inversion des matrices de covariance
      NPLU = I
 
      IFAIL = 0
      CALL DSINV(NPLU, COV, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'ERROR'
      IFAIL = 0
      CALL DSINV(NPLU, COVY, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'ERROR'
*      PRINT *,' COVARIANCE MATRIX AFTER'
*      DO I = 1, NPLANE
*         PRINT *,(COV(J,I),J=1,NPLANE)
*      ENDDO
 
** Calcul de FVAL ou CHI2
      FVAL = 0.0D0
      I = 0
      DO II = 1,NPLANE
        IF (LPLANE(II).EQ.1) THEN
        I = I+1
*        I = II
        J = 0
        DO JJ = 1,NPLANE
           IF (LPLANE(JJ).EQ.1) THEN
              J = J+1
*             J = JJ
              FVAL = FVAL + COV(J,I)*(XMP(II)-XP(II))*(XMP(JJ)-XP(JJ))
              FVAL = FVAL + COVY(J,I)*(YMP(II)-YP(II))
     &                               *(YMP(JJ)-YP(JJ))
**             IF (JJ.EQ.II) THEN
**                 FVAL = FVAL + (XM(II)-XP(II))*(XM(JJ)-XP(JJ))/XPREC**2
**                 FVAL = FVAL + (YM(II)-YP(II))
**     &                               *(YM(JJ)-YP(JJ))/YPREC**2
**             ENDIF
           ENDIF
        ENDDO
        ENDIF
      ENDDO
      CHI2 = FVAL

**      IF (CHI2.GT.1.E4) THEN
**         PRINT *,'FCNFIT CHI2= ',CHI2
**         FVAL = 0.
**      ENDIF

      
 1000 FORMAT(I5,7F12.6)
 
      RETURN
      END

************************************************************************        
      SUBROUTINE NEWFCNFIT(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)
************************************************************************        
*    With magnetic Field Map GRKUTA
*      trackfinding
*    Calcule FVAL: la fonction minimisee par MINUIT
*    With magnetic field map
*      
************************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(NPLANE=10)
            
      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
     
      COMMON/MEAS/LPLANE(NPLANE),XMP(NPLANE),YMP(NPLANE)
     
      COMMON/FCNOUT/PXZEA,ALAMEA,PHIEA,XEA,YEA,NPLU,CHI2
                 
      DIMENSION GRAD(*),XVAL(*),AMS(500),DISTAZ(500)
  
      DIMENSION XP(NPLANE),YP(NPLANE),
     &          COV(NPLANE,NPLANE),AP(NPLANE),COVY(NPLANE,NPLANE)
      DIMENSION VECT(7),VOUT(7) 
      
      STEP = 2.   ! 1 cm
      NSTEPMAX = 5000
      PITODEG = 57.295
      XV = XVAL(4)   
      YV = XVAL(5)   
              
      ASIGN = 1.
      IF (XVAL(1).LT.0.) ASIGN = -1.
      PHI  = XVAL(2)
      ALAM = XVAL(3)
      PXZ = DABS(1./XVAL(1))
           
      PX = PXZ*DSIN(PHI)
      PY = PXZ*DTAN(ALAM)
      PZ = PXZ*DCOS(PHI)
      PTOT = PXZ/DCOS(ALAM)
      
      A12 = 0.
      NTMAX = 0
      
      ZEA = ZABS
      XEA = XV
      YEA = YV
      PXEA = PX
      PYEA = PY
      PHIEA = PHI
      PXZEA = ASIGN*PXZ
      ALAMEA = ALAM

      VECT(1) = XV
      VECT(2) = YV
      VECT(3) = ZABS
      VECT(4) = PX/PTOT
      VECT(5) = PY/PTOT
      VECT(6) = PZ/PTOT
      VECT(7) = PTOT
      
      R = SQRT(VECT(1)*VECT(1)+VECT(2)*VECT(2))

      Z = VECT(3)
      NSTEP = 0
      IX = 0
      IY = 0 
      IZ = 0
**      PRINT *,' AV GRKUTA ASIGN',ASIGN,' THET',THET
      DO ICH = 1,NPLANE 
         DO WHILE (Z.GE.ZABS.AND.Z.LT.ZPLANEP(ICH)
     &            .AND.NSTEP.LE.NSTEPMAX)
**     &            .AND.(THETA*PITODEG).GT.2.
**     &            .AND. (THETA*PITODEG).LT.9.) 
            NSTEP = NSTEP+1 
**          WRITE(6,*) NSTEP,(VECT(I),I=1,7)
            CALL RECO_GRKUTA (ASIGN,STEP,VECT,VOUT)
            DO I = 1,7
               VECT(I) = VOUT(I)
            ENDDO   
            Z = VECT(3)
            R = SQRT(VECT(1)*VECT(1)+VECT(2)*VECT(2))
         ENDDO
         IF (NSTEP.EQ.NSTEPMAX) RETURN
         XP(ICH) = VECT(1)-(Z-ZPLANEP(ICH))*VECT(4)/VECT(6)
         YP(ICH) = VECT(2)-(Z-ZPLANEP(ICH))*VECT(5)/VECT(6)
         AL   = THICK/ VECT(6)
         AP(ICH) = (0.0136D0/PTOT)*DSQRT(AL)*(1+0.038D0*DLOG(AL))
      ENDDO
**    PRINT *,' AP GRKUTA ASIGN',ASIGN,' THET',THET


** Matrice de covariance      
      I = 0
      DO II = 1,NPLANE
        IF (LPLANE(II).EQ.1) THEN
           I = I + 1
*     I = II
           J = I - 1
           DO JJ = II, NPLANE
              IF (LPLANE(JJ).EQ.1) THEN
                 J = J + 1
*     J = JJ
                 COV (I,J) = 0.0D0
                 COV (J,I) = A12
                 IF (I .EQ. J) THEN
                    COV(J,I) =COV(J,I) + XPREC**2
                 ENDIF      
                 
*     IF (I .EQ. 10 .AND. J .EQ. 10) PRINT *,'10 10   ',COV(J,I)
                 DO L = 1,NTMAX
                    COV(J,I) = COV(J,I)
     &                   +  (ZPLANEP(II) + DISTAZ(L))*(ZPLANEP(JJ) + 
     &                   DISTAZ(L))*AMS(L)**2
                 ENDDO
                 DO K = 1, II-1
                    COV(J,I) = COV(J,I)
     &                   +  (ZPLANEP(II)-ZPLANEP(K))*
     &                   (ZPLANEP(JJ)-ZPLANEP(K))*AP(K)**2
*     IF (I .EQ. 10 .AND. J .EQ. 10) PRINT *,'10 10   ',COV(J,I)
                 ENDDO
                 COVY(I,J) = 0.0D0
                 COVY(J,I) = COV(J,I)
                 IF (I .EQ. J) THEN
                    COVY(J,I) = COVY(J,I) - XPREC**2 + YPREC**2
                 ENDIF
              ENDIF   
           ENDDO
        ENDIF
      ENDDO
 
*  Inversion des matrices de covariance
      NPLU = I
 
      IFAIL = 0
      CALL DSINV(NPLU, COV, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'ERROR'
      IFAIL = 0
      CALL DSINV(NPLU, COVY, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'ERROR'
*      PRINT *,' COVARIANCE MATRIX AFTER'
*      DO I = 1, NPLANE
*         PRINT *,(COV(J,I),J=1,NPLANE)
*      ENDDO
 
** Calcul de FVAL ou CHI2
      FVAL = 0.0D0
      I = 0
      DO II = 1,NPLANE
        IF (LPLANE(II).EQ.1) THEN
        I = I+1
*        I = II
        J = 0
        DO JJ = 1,NPLANE
           IF (LPLANE(JJ).EQ.1) THEN
              J = J+1
*             J = JJ
              FVAL = FVAL + COV(J,I)*(XMP(II)-XP(II))*(XMP(JJ)-XP(JJ))
              FVAL = FVAL + COVY(J,I)*(YMP(II)-YP(II))
     &                               *(YMP(JJ)-YP(JJ))
**             IF (JJ.EQ.II) THEN
**                 FVAL = FVAL + (XM(II)-XP(II))*(XM(JJ)-XP(JJ))/XPREC**2
**                 FVAL = FVAL + (YM(II)-YP(II))
**     &                               *(YM(JJ)-YP(JJ))/YPREC**2
**             ENDIF
           ENDIF
        ENDDO
        ENDIF
      ENDDO
      CHI2 = FVAL

**      IF (CHI2.GT.1.E4) THEN
**         PRINT *,'FCNFIT CHI2= ',CHI2
**         FVAL = 0.
**      ENDIF

      
 1000 FORMAT(I5,7F12.6)
 
      RETURN
      END

***********************************************************************
      SUBROUTINE INITFIELDOLD
*
*   Galina
***********************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

**      IMPLICIT REAL*8(A-H,O-Z)
**      REAL *4 BX,BY,BZ
      COMMON/DAT1/Z(81),X(81),Y(81,44),DX,DZ,LPX,LPY,LPZ
      COMMON/DAT2/BX(81,81,44),BY(81,81,44),BZ(81,81,44)
      COMMON/REG1/ZMAX,ZMIN,XMAX,XMIN
      COMMON/REG2/AY1,CY1,AY2,CY2
**      REAL *4 BXP,BYP,BZP
      COMMON/SDAT1/ZP(51),RAD(10),FI(33),DZP,DFI,DR,YY0,LPPZ,NR,NFI
      COMMON/SDAT2/BXP(51,10,33),BYP(51,10,33),BZP(51,10,33)
      COMMON/SDAT4/B(2,2,32)
      COMMON/REG3/ZPMAX,ZPMIN,RMAX,RMIN      
cc      COMMON/CONST/PI2,EPS
      REWIND 40
 1000 FORMAT(5(1X,D15.7))
 2000 FORMAT(5(1X,I5))
      READ(40,2000) LPX,LPY,LPZ
      READ(40,1000) (Z(K),K=1,81)
      READ(40,1000) (X(K),K=1,81)
      READ(40,1000) DX,DY,DZ
      READ(40,1000) ZMAX,ZMIN,XMAX,XMIN
c      write(*,*) 'zmin zmax',ZMIN,ZMAX
c      write(*,*) 'xmin xmax',XMIN,XMAX      
      READ(40,1000) AY1,CY1,AY2,CY2
c      write(*,*) 'ay1,cy1,ay2,cy2', AY1,CY1,AY2,CY2 
cc      READ(40,1000) PI2,EPS
      READ(40,1000) (((BX(K,L,M),K=1,81),L=1,81),M=1,44)
      READ(40,1000) (((BY(K,L,M),K=1,81),L=1,81),M=1,44)
      READ(40,1000) (((BZ(K,L,M),K=1,81),L=1,81),M=1,44)
**      RETURN
**      END
c Polar part
      READ(40,2000) LPPZ,NR,NFI
      READ(40,1000) (ZP(K),K=1,51)
      READ(40,1000) (RAD(K),K=1,10)
      READ(40,1000) (FI(L),L=1,33)
      READ(40,1000) DZP,DFI,DR
c      write(*,*) 'dzp dfi dR',DZP,DFI,DR
      READ(40,1000) ZPMAX,ZPMIN,RMAX,RMIN
c      write(*,*) 'zmin zmax',ZPMIN,ZPMAX
c     write(*,*) 'Rmin Rmax',RMIN,RMAX
      READ(40,1000) (((BXP(K,L,M),K=1,51),L=1,10),M=1,33)
      READ(40,1000) (((BYP(K,L,M),K=1,51),L=1,10),M=1,33)
      READ(40,1000) (((BZP(K,L,M),K=1,51),L=1,10),M=1,33)
      READ(40,1000) (((B(K,L,M),K=1,2),L=1,2),M=1,32)
**      RETURN
**      END
      



      RETURN
      END

***********************************************************************
      SUBROUTINE RECO_GUFLDOLD(X,F)
C     ^^^^^^^^^^^^^^^^^^^^^^
C   field map G. Chabratova
C
C  Field map 31/05/99
***********************************************************************


      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON/MAGERR/IMAGERR

      DIMENSION X(7),F(3)

      XT = X(2)
      X(2) = X(1)
      X(1) = XT

      X0 = X(1)/100.
      Y0 = X(2)/100.
      Z0 = X(3)/100.

      CALL FREG1(Z0,X0,Y0,FZ0,FX0,FY0,IND)
**    PRINT 3000,Z0,X0,Y0,FZ0,FX0,FY0,IND

      IF(IND.EQ.0) GOTO 1
      CALL FREG2(Z0,X0,Y0,FZ0,FX0,FY0,IND)
      IMAGERR = 0
      IF(IND.EQ.2) THEN
        IMAGERR = 1 
**        print 1000
**        PRINT 3000,Z0,X0,Y0,FZ0,FX0,FY0,IND        
      ENDIF 
 1000 format(1x,'Attention!!! The point is out of range!!!')

 3000 FORMAT(1X,'Z=',D13.7,1X,'X=',D13.7,1X,'Y=',D13.7,1X,
     & 'BZ=',D13.7,1X,'BX=',D13.7,1X,'BY=',D13.7,1X,'IND=',I3)

  1   F(1) = FX0*10.
      F(2) = FY0*10.
      F(3) = FZ0*10. 

**      X(1) = X0*100.
**      X(2) = Y0*100.
**      X(3) = Z0*100.
      
           
      FT = F(2)
      F(2) = F(1)
      F(1) = FT

      XT = X(2)
      X(2) = X(1)
      X(1) = XT

      RETURN
      END
      
      
**************************************************
      SUBROUTINE FREG1(Z0,X0,Y0,FZ0,FX0,FY0,IND)
**************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
**      REAL *4 BX,BY,BZ
      COMMON/DAT1/Z(81),X(81),Y(81,44),DX,DZ,LPX,LPY,LPZ
      COMMON/DAT2/BX(81,81,44),BY(81,81,44),BZ(81,81,44)
      COMMON/REG1/ZMAX,ZMIN,XMAX,XMIN
      COMMON/REG2/AY1,CY1,AY2,CY2
      KC=1
      LC=1
      MC=1
      IND=0
      IF(Z0.LT.ZMIN.OR.Z0.GT.ZMAX)GO TO 100
      IF(X0.LT.XMIN.OR.X0.GT.XMAX)GO TO 100
      YY1=AY1*Z0+CY1
      YY2=AY2*Z0+CY2
 2000 FORMAT(1X,'YY1=',D15.7,1X,'YY2=',D15.7,1X,'Y0=',D15.7)
c      PRINT 2000,YY1,YY2,Y0
      IF(Y0.LT.YY1)GO TO 100
      IF(Y0.GT.YY2)GO TO 100
      CALL FIZ(Z0,Z,DZ,KC,K0,Z1,Z2,81)
      CALL FIZ(X0,X,DX,LC,L0,X1,X2,81)
      DY=(YY2-YY1)/DFLOAT(LPY-1)
      YY=(Y0-YY1)
      M0=(YY/DY)
      
      IF(Y0.GE.(YY1+DFLOAT(M0)*DY).AND.Y0.LE.(YY1+DFLOAT(M0+1)*DY))
     &GO TO 700
      M0=M0+1
 700  CONTINUE
      Y2=(Y0-(YY1+DFLOAT(M0)*DY))/DY
      Y1=1.-Y2
**    write(*,*) 'm0 Y1 Y2',m0,Y1,Y2
**    print *,' k0=',k0,' l0=',l0,' m0=',m0 
**    print *,' z1=',z1,' z2=',z2
      FX1=Z1*BX(K0,L0,M0)+Z2*BX(K0+1,L0,M0)
      FX2=Z2*BX(K0+1,L0+1,M0)+Z1*BX(K0,L0+1,M0)
      FFX1=X1*FX1+X2*FX2
      GX1=Z1*BX(K0,L0,M0+1)+Z2*BX(K0+1,L0,M0+1)
      GX2=Z2*BX(K0+1,L0+1,M0+1)+Z1*BX(K0,L0+1,M0+1)
      GGX1=X1*GX1+X2*GX2
      FX0=Y1*FFX1+Y2*GGX1
      FX1=Z1*BY(K0,L0,M0)+Z2*BY(K0+1,L0,M0)
      FX2=Z2*BY(K0+1,L0+1,M0)+Z1*BY(K0,L0+1,M0)
      FFX1=X1*FX1+X2*FX2
      GX1=Z1*BY(K0,L0,M0+1)+Z2*BY(K0+1,L0,M0+1)
      GX2=Z2*BY(K0+1,L0+1,M0+1)+Z1*BY(K0,L0+1,M0+1)
      GGX1=X1*GX1+X2*GX2
      FY0=Y1*FFX1+Y2*GGX1
      FX1=Z1*BZ(K0,L0,M0)+Z2*BZ(K0+1,L0,M0)
      FX2=Z2*BZ(K0+1,L0+1,M0)+Z1*BZ(K0,L0+1,M0)
      FFX1=X1*FX1+X2*FX2
      GX1=Z1*BZ(K0,L0,M0+1)+Z2*BZ(K0+1,L0,M0+1)
      GX2=Z2*BZ(K0+1,L0+1,M0+1)+Z1*BZ(K0,L0+1,M0+1)
      GGX1=X1*GX1+X2*GX2
      FZ0=Y1*FFX1+Y2*GGX1
      RETURN
 100  CONTINUE
      IND=1
 1000 format(1x,'Attention!!! The point is out of range!!!')
C      print 1000
      RETURN
      END
      
*************************************************
      SUBROUTINE FIZ(Z0,Z,DEL,KI,K0,Z1,Z2,NDZ)
*************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
** CC      DIMENSION Z(NDZ)
      DIMENSION Z(10000)
      DDEL=Z0-Z(KI)
      KDEL=INT(DDEL/DEL)
      KJ=KI+KDEL
      K0 = NDZ - 1 ! CCCC
*      if (k0.gt.81) print*,'ndz=',ndz
      IF (KJ.GT.NDZ) THEN ! CCC
         K0 = NDZ-1
         GO TO 100
      ENDIF 
      DO 1 K=KJ,NDZ-1 ! CCC
      IF(Z0.LT.Z(K)) THEN
        K0 = K 
        GO TO 100
      ENDIF
 1    CONTINUE
 100  CONTINUE
*      print *,'K0=',K0,' Z0',z0, Z(K0), Z(K0+1),z2
      if (k0.gt.81) print*,'k0=',k0
      Z2=(Z0-Z(K0))/(Z(K0+1)-Z(K0))
      Z1=1.-Z2
**      write(*,*) 'ko z1 z2', K0,Z1,Z2,' ki=',ki,' kj=',kj,' K=',K
**      write(*,*)' NDZ Z(K0) Z(K0+1)',NDZ,Z(K0), Z(K0+1)
      RETURN
      END

***************************************************
      SUBROUTINE FREG2(Z0,X0,Y0,FZ0,FX0,FY0,IND)
**************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
**      REAL *4 BXP,BYP,BZP
      COMMON/SDAT1/ZP(51),RAD(10),FI(33),DZP,DFI,DR,YY0,LPPZ,NR,NFI
      COMMON/SDAT2/BXP(51,10,33),BYP(51,10,33),BZP(51,10,33)
      COMMON/SDAT4/B(2,2,32)
      COMMON/REG3/ZPMAX,ZPMIN,RMAX,RMIN
cc      COMMON/CONST/PI2,EPS
      KP=32+1
      LP=1
      MP=1
      YY0=0.3
      EPS=0.1D-6
      PI2=0.6283185E+01
      R0=DSQRT((X0-YY0)**2+Y0**2)
c       write (*,*)'ro=',R0
      IF(Z0.LT.ZPMIN.OR.Z0.GT.ZPMAX)GO TO 100
      IF(R0.LT.RMIN.OR.R0.GT.RMAX)GO TO 100
      IF(R0.LE.DR)GO TO 3000
      CALL FIZ(Z0,ZP,DZP,KP,K0,Z1,Z2,51)
      CALL FIZ(R0,RAD,DR,LP,L0,X1,X2,10)
**      print *,' r0=',r0,' rad=',rad,' dr=',dr,' lp=',lp,' l0=',l0,
**     &     ' x1=',x1,' x2=',x2
      FI0=DACOS((X0-YY0)/R0)
      IF(Y0.LT.0.D+0)FI0=PI2-FI0
      CALL FIZ(FI0,FI,DFI,MP,M0,Y1,Y2,32)
**      print *,' Apres FIZ',' k0=',k0,' l0=',l0,' m0=',m0
      FX1=Z1*BXP(K0,L0,M0)+Z2*BXP(K0+1,L0,M0)
      FX2=Z2*BXP(K0+1,L0+1,M0)+Z1*BXP(K0,L0+1,M0)
      FFX1=X1*FX1+X2*FX2
      GX1=Z1*BXP(K0,L0,M0+1)+Z2*BXP(K0+1,L0,M0+1)
      GX2=Z2*BXP(K0+1,L0+1,M0+1)+Z1*BXP(K0,L0+1,M0+1)
      GGX1=X1*GX1+X2*GX2
      FX0=Y1*FFX1+Y2*GGX1
      FX1=Z1*BYP(K0,L0,M0)+Z2*BYP(K0+1,L0,M0)
      FX2=Z2*BYP(K0+1,L0+1,M0)+Z1*BYP(K0,L0+1,M0)
      FFX1=X1*FX1+X2*FX2
      GX1=Z1*BYP(K0,L0,M0+1)+Z2*BYP(K0+1,L0,M0+1)
      GX2=Z2*BYP(K0+1,L0+1,M0+1)+Z1*BYP(K0,L0+1,M0+1)
      GGX1=X1*GX1+X2*GX2
      FY0=Y1*FFX1+Y2*GGX1
      FX1=Z1*BZP(K0,L0,M0)+Z2*BZP(K0+1,L0,M0)
** CCC      FX2=Z2*BZ(K0+1,L0+1,M0)+Z1*BZ(K0,L0+1,M0)
      FX2=Z2*BZP(K0+1,L0+1,M0)+Z1*BZP(K0,L0+1,M0)
      FFX1=X1*FX1+X2*FX2
      GX1=Z1*BZP(K0,L0,M0+1)+Z2*BZP(K0+1,L0,M0+1)
      GX2=Z2*BZP(K0+1,L0+1,M0+1)+Z1*BZP(K0,L0+1,M0+1)
      GGX1=X1*GX1+X2*GX2
      FZ0=Y1*FFX1+Y2*GGX1
      IND=0
      RETURN
 100  CONTINUE
      IND=2
 1000 format(1x,'Attention!!! The point is out of range!!!')
C      print 1000
      RETURN
 3000 CONTINUE
      IF(R0.LT.EPS)GO TO 4000
      CALL FIZ(Z0,ZP,DZP,KP,K0,Z1,Z2,51)
      XX=X0-YY0
      FI0=DACOS(XX/R0)
      IF(Y0.LT.0.D+0)FI0=PI2-FI0
      CALL FIZ(FI0,FI,DFI,MP,M0,Y1,Y2,32)
      ALF2=B(1,1,M0)*XX+B(1,2,M0)*Y0
      ALF3=B(2,1,M0)*XX+B(2,2,M0)*Y0
      ALF1=1.-ALF2-ALF3
      FX1=ALF1*BXP(K0,1,1)+ALF2*BXP(K0,1,M0)+ALF3*BXP(K0,1,M0+1)
      FX2=ALF1*BXP(K0+1,1,1)+ALF2*BXP(K0+1,1,M0)+ALF3*BXP(K0+1,1,M0+1)
      FX0=Z1*FX1+Z2*FX2
      FX1=ALF1*BYP(K0,1,1)+ALF2*BYP(K0,1,M0)+ALF3*BYP(K0,1,M0+1)
      FX2=ALF1*BYP(K0+1,1,1)+ALF2*BYP(K0+1,1,M0)+ALF3*BYP(K0+1,1,M0+1)
      FY0=Z1*FX1+Z2*FX2
      FX1=ALF1*BZP(K0,1,1)+ALF2*BZP(K0,1,M0)+ALF3*BZP(K0,1,M0+1)
      FX2=ALF1*BZP(K0+1,1,1)+ALF2*BZP(K0+1,1,M0)+ALF3*BZP(K0+1,1,M0+1)
      FZ0=Z1*FX1+Z2*FX2 
c      write(*,*) 'R<Dr:B(1,1,m0) B(1,2,m0) B(2,1,m0) B(2,2,m0)',
c     +B(1,1,M0),B(1,2,M0),B(2,1,M0),B(2,2,M0)
c      write(*,*)'BX(K0,1,1) BX(K0,1,M0)','BX(K0,1,M0+1) BX(K0+1,1,1)' 
c     + ,'BX(K0+1,1,M0) BX(K0+1,1,M0+1)', BX(K0,1,1),BX(K0,1,M0) ,
c     + BX(K0,1,M0+1),BX(K0+1,1,1),BX(K0+1,1,M0),BX(K0+1,1,M0+1) 
c      write(*,*)'By(K0,1,1) By(K0,1,M0)','By(K0,1,M0+1) By(K0+1,1,1)' 
c     + ,'By(K0+1,1,M0) By(K0+1,1,M0+1)', By(K0,1,1),By(K0,1,M0) ,
c     + By(K0,1,M0+1),By(K0+1,1,1),By(K0+1,1,M0),By(K0+1,1,M0+1) 
ccc      write (*,*)'Bz(K0,1,1) Bz(K0,1,M0) Bz(K0,1,M0+1) Bz(K0+1,1,1) 
  77  FORMAT(5x,E15.7,2x,E15.7)
c      PRINT 70
**   70 FORMAT(22hBz(K0,1,1) Bz(K0,1,M0))
cc      PRINT 77 , BzP(K0,1,1),BzP(K0,1,M0)
cc      PRINT 71
**   71 FORMAT(26hBz(K0,1,M0+1) Bz(K0+1,1,1))
cc     PRINT 77, BzP(K0,1,M0+1),BzP(K0+1,1,1)       
cc     PRINT 72
**   72 FORMAT(29hBz(K0+1,1,M0) Bz(K0+1,1,M0+1))  
cc      PRINT 77 ,BzP(K0+1,1,M0),BzP(K0+1,1,M0+1) 
cc   77 FORMAT(5x,D15.7,5x,D15.7)
          
 
      IND=0
      RETURN
 4000 CONTINUE
      CALL FIZ(Z0,ZP,DZP,KP,K0,Z1,Z2,51)
      FX0=Z1*BXP(K0,1,1)+Z2*BXP(K0+1,1,1)
      FY0=Z1*BYP(K0,1,1)+Z2*BYP(K0+1,1,1)
      FZ0=Z1*BZP(K0,1,1)+Z2*BZP(K0+1,1,1)
c      write(*,*) ' R<eps: Bx(k) Bx(k+1) By(k) By(k+1) Bz(k) Bz(k+1)',
c     + BXP(K0,1,1),BXP(K0+1,1,1),BYP(K0,1,1),BYP(K0+1,1,1),BZP(K0,1,1),
c     + BZP(K0+1,1,1)
      IND=0
      RETURN
      END

***********************************************************************        
      
      SUBROUTINE RECO_GRKUTA (CHARGE,STEP,VECT,VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Runge-Kutta method for tracking a particle through a magnetic *
C.    *  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of     *
C.    *  Standards, procedure 25.5.20)                                 *
C.    *                                                                *
C.    *  Input parameters                                              *
C.    *       CHARGE    Particle charge                                *
C.    *       STEP      Step size                                      *
C.    *       VECT      Initial co-ords,direction cosines,momentum     *
C.    *  Output parameters                                             *
C.    *       VOUT      Output co-ords,direction cosines,momentum      *
C.    *  User routine called                                           *
C.    *       CALL GUFLD(X,F)                                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Authors    R.Brun, M.Hansroul  *********                 *
C.    *                  V.Perevoztchikov (CUT STEP implementation)    *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
**      REAL CHARGE, STEP, VECT(*), VOUT(*), F(4)
**      REAL XYZT(3), XYZ(3), X, Y, Z, XT, YT, ZT
      DIMENSION VECT(*), VOUT(*), F(3)
      DIMENSION XYZT(3), XYZ(3)
      DIMENSION SECXS(4),SECYS(4),SECZS(4),HXP(3)
      EQUIVALENCE (X,XYZ(1)),(Y,XYZ(2)),(Z,XYZ(3)),
     +            (XT,XYZT(1)),(YT,XYZT(2)),(ZT,XYZT(3))
*
      PARAMETER (MAXIT = 1992, MAXCUT = 11)
      PARAMETER (EC=2.9979251D-4,DLT=1D-4,DLT32=DLT/32)
      PARAMETER (ZERO=0, ONE=1, TWO=2, THREE=3)
      PARAMETER (THIRD=ONE/THREE, HALF=ONE/TWO)
      PARAMETER (PISQUA=.986960440109D+01)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6)
*.
*.    ------------------------------------------------------------------
*.
*             This constant is for units CM,GEV/C and KGAUSS
*
      ITER = 0
      NCUT = 0
      DO 10 J=1,7
         VOUT(J)=VECT(J)
   10 CONTINUE
      PINV   = EC * CHARGE / VECT(7)
      TL = 0.
      H      = STEP
*
*
   20 REST  = STEP-TL
      IF (ABS(H).GT.ABS(REST)) H = REST
      CALL RECO_GUFLD(VOUT,F)
*
*             Start of integration
*
      X      = VOUT(1)
      Y      = VOUT(2)
      Z      = VOUT(3)
      A      = VOUT(4)
      B      = VOUT(5)
      C      = VOUT(6)
*
      H2     = HALF * H
      H4     = HALF * H2
      PH     = PINV * H
      PH2    = HALF * PH
      SECXS(1) = (B * F(3) - C * F(2)) * PH2
      SECYS(1) = (C * F(1) - A * F(3)) * PH2
      SECZS(1) = (A * F(2) - B * F(1)) * PH2
      ANG2 = (SECXS(1)**2 + SECYS(1)**2 + SECZS(1)**2)
      IF (ANG2.GT.PISQUA) GO TO 40
      DXT    = H2 * A + H4 * SECXS(1)
      DYT    = H2 * B + H4 * SECYS(1)
      DZT    = H2 * C + H4 * SECZS(1)
      XT     = X + DXT
      YT     = Y + DYT
      ZT     = Z + DZT
*
*              Second intermediate point
*
      EST = ABS(DXT)+ABS(DYT)+ABS(DZT)
      IF (EST.GT.H) GO TO 30
 
      CALL RECO_GUFLD(XYZT,F)
      AT     = A + SECXS(1)
      BT     = B + SECYS(1)
      CT     = C + SECZS(1)
*
      SECXS(2) = (BT * F(3) - CT * F(2)) * PH2
      SECYS(2) = (CT * F(1) - AT * F(3)) * PH2
      SECZS(2) = (AT * F(2) - BT * F(1)) * PH2
      AT     = A + SECXS(2)
      BT     = B + SECYS(2)
      CT     = C + SECZS(2)
      SECXS(3) = (BT * F(3) - CT * F(2)) * PH2
      SECYS(3) = (CT * F(1) - AT * F(3)) * PH2
      SECZS(3) = (AT * F(2) - BT * F(1)) * PH2
      DXT    = H * (A + SECXS(3))
      DYT    = H * (B + SECYS(3))
      DZT    = H * (C + SECZS(3))
      XT     = X + DXT
      YT     = Y + DYT
      ZT     = Z + DZT
      AT     = A + TWO*SECXS(3)
      BT     = B + TWO*SECYS(3)
      CT     = C + TWO*SECZS(3)
*
      EST = ABS(DXT)+ABS(DYT)+ABS(DZT)
      IF (EST.GT.2.*ABS(H)) GO TO 30
 
      CALL RECO_GUFLD(XYZT,F)
*
      Z      = Z + (C + (SECZS(1) + SECZS(2) + SECZS(3)) * THIRD) * H
      Y      = Y + (B + (SECYS(1) + SECYS(2) + SECYS(3)) * THIRD) * H
      X      = X + (A + (SECXS(1) + SECXS(2) + SECXS(3)) * THIRD) * H
*
      SECXS(4) = (BT*F(3) - CT*F(2))* PH2
      SECYS(4) = (CT*F(1) - AT*F(3))* PH2
      SECZS(4) = (AT*F(2) - BT*F(1))* PH2
      A      = A+(SECXS(1)+SECXS(4)+TWO * (SECXS(2)+SECXS(3))) * THIRD
      B      = B+(SECYS(1)+SECYS(4)+TWO * (SECYS(2)+SECYS(3))) * THIRD
      C      = C+(SECZS(1)+SECZS(4)+TWO * (SECZS(2)+SECZS(3))) * THIRD
*
      EST    = ABS(SECXS(1)+SECXS(4) - (SECXS(2)+SECXS(3)))
     ++        ABS(SECYS(1)+SECYS(4) - (SECYS(2)+SECYS(3)))
     ++        ABS(SECZS(1)+SECZS(4) - (SECZS(2)+SECZS(3)))
*
      IF (EST.GT.DLT .AND. ABS(H).GT.1.E-4) GO TO 30
      ITER = ITER + 1
      NCUT = 0
*               If too many iterations, go to HELIX
      IF (ITER.GT.MAXIT) GO TO 40
*
      TL = TL + H
      IF (EST.LT.(DLT32)) THEN
         H = H*TWO
      ENDIF
      CBA    = ONE/ SQRT(A*A + B*B + C*C)
      VOUT(1) = X
      VOUT(2) = Y
      VOUT(3) = Z
      VOUT(4) = CBA*A
      VOUT(5) = CBA*B
      VOUT(6) = CBA*C
      REST = STEP - TL
      IF (STEP.LT.0.) REST = -REST
      IF (REST .GT. 1.E-5*ABS(STEP)) GO TO 20
*
      GO TO 999
*
**              CUT STEP
   30 NCUT = NCUT + 1
*               If too many cuts , go to HELIX
      IF (NCUT.GT.MAXCUT)       GO TO 40
      H = H*HALF
      GO TO 20
*
**              ANGLE TOO BIG, USE HELIX
   40 F1  = F(1)
      F2  = F(2)
      F3  = F(3)
      F4  = SQRT(F1**2+F2**2+F3**2)
      RHO = -F4*PINV
      TET = RHO * STEP
      IF(TET.NE.0.) THEN
         HNORM = ONE/F4
         F1 = F1*HNORM
         F2 = F2*HNORM
         F3 = F3*HNORM
*
         HXP(1) = F2*VECT(IPZ) - F3*VECT(IPY)
         HXP(2) = F3*VECT(IPX) - F1*VECT(IPZ)
         HXP(3) = F1*VECT(IPY) - F2*VECT(IPX)
 
         HP = F1*VECT(IPX) + F2*VECT(IPY) + F3*VECT(IPZ)
*
         RHO1 = ONE/RHO
         SINT = SIN(TET)
         COST = TWO*SIN(HALF*TET)**2
*
         G1 = SINT*RHO1
         G2 = COST*RHO1
         G3 = (TET-SINT) * HP*RHO1
         G4 = -COST
         G5 = SINT
         G6 = COST * HP
 
         VOUT(IX) = VECT(IX) + (G1*VECT(IPX) + G2*HXP(1) + G3*F1)
         VOUT(IY) = VECT(IY) + (G1*VECT(IPY) + G2*HXP(2) + G3*F2)
         VOUT(IZ) = VECT(IZ) + (G1*VECT(IPZ) + G2*HXP(3) + G3*F3)
 
         VOUT(IPX) = VECT(IPX) + (G4*VECT(IPX) + G5*HXP(1) + G6*F1)
         VOUT(IPY) = VECT(IPY) + (G4*VECT(IPY) + G5*HXP(2) + G6*F2)
         VOUT(IPZ) = VECT(IPZ) + (G4*VECT(IPZ) + G5*HXP(3) + G6*F3)
*
      ELSE
         VOUT(IX) = VECT(IX) + STEP*VECT(IPX)
         VOUT(IY) = VECT(IY) + STEP*VECT(IPY)
         VOUT(IZ) = VECT(IZ) + STEP*VECT(IPZ)
*
      ENDIF
*
  999 END

  
*******************************************************************
      SUBROUTINE RECO_GUFLDOLD1(X,B)
C
C    CONSTANT FIELD
C
C *** ROUTINE DESCRIBING THE MAGNETIC FIELD IN THE ALICE SETUP ***
C *** NVE 14-NOV-1990 CERN GENEVA ***
C
C CALLED BY : GUFLD
C ORIGIN    : NICK VAN EIJNDHOVEN
C
C Input :
C -------
C X = (X,Y,Z) coordinates in cm
C
C Output :
C --------
C B    = Magnetic field components (BX,BY,BZ) in KG
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

      PARAMETER(NBSTATION=5)
C
      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)


C
      DIMENSION X(3),B(3)

      XT = X(2)
      X(2) = X(1)
      X(1) = XT

      B(1) = 0.
      B(2) = 0.
      B(3) = 0.
      
      IF (X(3).LT.(-1.*ZCOIL)) THEN
         B(3) = 2.
      ELSEIF ( X(3).LT.(-1.*ZMAGEND)) THEN
         B(1) = -10.
      ENDIF

      BT = B(2)
      B(2) = B(1)
      B(1) = BT

      XT = X(2)
      X(2) = X(1)
      X(1) = XT
**      print *,' x =',X(1),X(2),X(3)
**      print *,' B =',B(1),B(2),B(3)

C


 999  END


*******************************************************************
      SUBROUTINE RECO_GUFLD(X,B)
C
C    Field map Mariana
C
C *** ROUTINE DESCRIBING THE MAGNETIC FIELD IN THE ALICE SETUP ***
C *** NVE 14-NOV-1990 CERN GENEVA ***
C
C CALLED BY : GUFLD
C ORIGIN    : NICK VAN EIJNDHOVEN
C
C Input :
C -------
C X = (X,Y,Z) coordinates in cm
C
C Output :
C --------
C B    = Magnetic field components (BX,BY,BZ) in KG
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
C

C --- Common containing magnetic field map data
      REAL DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV
      INTEGER NX,NY,NZ
 
      PARAMETER(MAXFLD=250000)
      COMMON /SCXMFD/ NX,NY,NZ,DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV(MAXFLD)
C

C
      DIMENSION X(3),B(3)
      DOUBLE PRECISION ONE, RATX, RATY, RATZ, HIX, HIY, HIZ
     $,  RATX1, RATY1, RATZ1
     $,  BHYHZ, BHYLZ, BLYHZ, BLYLZ, BHZ, BLZ
     $,  XL(3)

      PARAMETER (ONE=1)

      EXTERNAL BX,BY,BZ

      XT = X(2)
      X(2) = X(1)
      X(1) = XT


      ISXFMAP = 3 

**      BX(JX,JY,JZ)=BV(3*((JZ-1)*(NX*NY)+(JY-1)*NX+(JX-1))+1)
**      BY(JX,JY,JZ)=BV(3*((JZ-1)*(NX*NY)+(JY-1)*NX+(JX-1))+2)
**      BZ(JX,JY,JZ)=BV(3*((JZ-1)*(NX*NY)+(JY-1)*NX+(JX-1))+3)



*
* --- Act accordingly to ISXFMAP
*
      IF(ISXFMAP.EQ.1) THEN
         IF (ABS(X(3)) .LT. 700.
     +      .AND. X(1)**2+(X(2)+30.)**2 .LT. 560.**2 ) THEN
            B(1)=0.
            B(2)=0.
            B(3)=2.
         ELSE IF (X(3) .GE. 725. .AND. X(3) .LT. 1225.) THEN
            DZ=ABS(975.-X(3))/100.
            B(1)=(1.-0.1*DZ*DZ)*7.0
            B(2)=0.
            B(3)=0.
         ELSE
            B(1)=0.
            B(2)=0.
            B(3)=0.
         ENDIF
      ELSE IF(ISXFMAP.LE.3) THEN
         IF (-700.LT.X(3).AND.X(3).LT.ZMBEG
     +      .AND. X(1)**2+(X(2)+30.)**2 .LT. 560.**2 ) THEN
            B(1)=0.
            B(2)=0.
            B(3)=0.
         ELSE IF ((X(3) .GE. ZMBEG .AND. X(3) .LT. ZMEND) .AND.
     +            (XMBEG.LE.ABS(X(1)).AND.ABS(X(1)).LT.XMEND) .AND.
     +            (YMBEG.LE.ABS(X(2)).AND.ABS(X(2)).LT.YMEND)) THEN



C --- find the position in the grid ---

           XL(1)=ABS(X(1))-XMBEG
           XL(2)=ABS(X(2))-YMBEG
           XL(3)=X(3)-ZMBEG

C --- Start with X

           HIX=XL(1)*UDX
           RATX=HIX-AINT(HIX)
           IX=HIX+1

           HIY=XL(2)*UDY
           RATY=HIY-AINT(HIY)
           IY=HIY+1

           HIZ=XL(3)*UDZ
           RATZ=HIZ-AINT(HIZ)
           IZ=HIZ+1

           IF(ISXFMAP.EQ.2) THEN
* ... Simple interpolation
              
              B(1) = BX(IX,IY,IZ)*(ONE-RATX) + BX(IX+1,IY+1,IZ+1)*RATX
              B(2) = BY(IX,IY,IZ)*(ONE-RATY) + BY(IX+1,IY+1,IZ+1)*RATY
              B(3) = BZ(IX,IY,IZ)*(ONE-RATZ) + BZ(IX+1,IY+1,IZ+1)*RATZ
           ELSE IF(ISXFMAP.EQ.3) THEN
* ... more complicated interpolation
              RATX1=ONE-RATX
              RATY1=ONE-RATY
              RATZ1=ONE-RATZ
**              print *,' bx by bz', BX(IX  ,IY+1,IZ+1),BY(IX  ,IY+1,IZ+1)
**     &        , BZ(IX  ,IY+1,IZ+1)       
              BHYHZ = BX(IX  ,IY+1,IZ+1)*RATX1+BX(IX+1,IY+1,IZ+1)*RATX
              BHYLZ = BX(IX  ,IY+1,IZ  )*RATX1+BX(IX+1,IY+1,IZ  )*RATX
              BLYHZ = BX(IX  ,IY  ,IZ+1)*RATX1+BX(IX+1,IY  ,IZ+1)*RATX
              BLYLZ = BX(IX  ,IY  ,IZ  )*RATX1+BX(IX+1,IY  ,IZ  )*RATX
              BHZ   = BLYHZ             *RATY1+BHYHZ             *RATY
              BLZ   = BLYLZ             *RATY1+BHYLZ             *RATY
              B(1)  = BLZ               *RATZ1+BHZ               *RATZ
*
              BHYHZ = BY(IX  ,IY+1,IZ+1)*RATX1+BY(IX+1,IY+1,IZ+1)*RATX
              BHYLZ = BY(IX  ,IY+1,IZ  )*RATX1+BY(IX+1,IY+1,IZ  )*RATX
              BLYHZ = BY(IX  ,IY  ,IZ+1)*RATX1+BY(IX+1,IY  ,IZ+1)*RATX
              BLYLZ = BY(IX  ,IY  ,IZ  )*RATX1+BY(IX+1,IY  ,IZ  )*RATX
              BHZ   = BLYHZ             *RATY1+BHYHZ             *RATY
              BLZ   = BLYLZ             *RATY1+BHYLZ             *RATY
              B(2)  = BLZ               *RATZ1+BHZ               *RATZ
*
              BHYHZ = BZ(IX  ,IY+1,IZ+1)*RATX1+BZ(IX+1,IY+1,IZ+1)*RATX
              BHYLZ = BZ(IX  ,IY+1,IZ  )*RATX1+BZ(IX+1,IY+1,IZ  )*RATX
              BLYHZ = BZ(IX  ,IY  ,IZ+1)*RATX1+BZ(IX+1,IY  ,IZ+1)*RATX
              BLYLZ = BZ(IX  ,IY  ,IZ  )*RATX1+BZ(IX+1,IY  ,IZ  )*RATX
              BHZ   = BLYHZ             *RATY1+BHYHZ             *RATY
              BLZ   = BLYLZ             *RATY1+BHYLZ             *RATY
              B(3)  = BLZ               *RATZ1+BHZ               *RATZ
*
           ENDIF
* ... use the dipole symmetry

           IF (X(1)*X(2).LT.0) B(2)=-B(2)
           IF (X(1).LT.0) B(3)=-B(3)

*

      ELSE

           B(1)=0.
           B(2)=0.
           B(3)=0.

      ENDIF ! z-coord for m.f.

      ENDIF ! endif ISXFMAP

      SXMAGN = 1.
      SXMGMX = 20.
 12   B(1)=B(1)*SXMAGN
      B(2)=B(2)*SXMAGN
      B(3)=B(3)*SXMAGN
      BTOT=SQRT(B(1)**2+B(2)**2+B(3)**2)
      IF(BTOT.GT.SXMGMX) THEN
         PRINT 10100, BTOT,SXMGMX
10100    FORMAT(' *GUFLD* Field ',G10.4,' larger than max ',G10.4)
      ENDIF


      BT = B(2)
      B(2) = B(1)
      B(1) = BT

      XT = X(2)
      X(2) = X(1)
      X(1) = XT

C

      RETURN
      END 


*******************************************
      DOUBLE PRECISION FUNCTION BX(JX,JY,JZ)
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C --- Common containing magnetic field map data
      REAL DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV
      INTEGER NX,NY,NZ
 
      PARAMETER(MAXFLD=250000)
      COMMON /SCXMFD/ NX,NY,NZ,DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV(MAXFLD)

      BX=BV(3*((JZ-1)*(NX*NY)+(JY-1)*NX+(JX-1))+1)

      END
 
*******************************************
       DOUBLE PRECISION FUNCTION BY(JX,JY,JZ)
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C --- Common containing magnetic field map data
      REAL DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV
      INTEGER NX,NY,NZ
 
      PARAMETER(MAXFLD=250000)
      COMMON /SCXMFD/ NX,NY,NZ,DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV(MAXFLD)

      BY=BV(3*((JZ-1)*(NX*NY)+(JY-1)*NX+(JX-1))+2)

      END

*******************************************
      DOUBLE PRECISION FUNCTION BZ(JX,JY,JZ)
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

C --- Common containing magnetic field map data
      REAL DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV
      INTEGER NX,NY,NZ
 
      PARAMETER(MAXFLD=250000)
      COMMON /SCXMFD/ NX,NY,NZ,DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV(MAXFLD)

      BZ=BV(3*((JZ-1)*(NX*NY)+(JY-1)*NX+(JX-1))+3)

      END

***********************************************************************
      SUBROUTINE INITFIELD

C
C    Marianna
C
C *** INITIALISATION OF THE FIELD MAP ***
C *** FCA 24-AUG-1998 CERN GENEVA ***
C
C CALLED BY : GALICE
C ORIGIN    : NICK VAN EIJNDHOVEN
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)


C --- Common containing magnetic field map data
      REAL DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV
      INTEGER NX,NY,NZ
 
      PARAMETER(MAXFLD=250000)
      COMMON /SCXMFD/ NX,NY,NZ,DZ,DX,DY,UDX,UDY,UDZ
     $,XMBEG,YMBEG,ZMBEG,XMEND,YMEND,ZMEND
     $,BV(MAXFLD)

      ISXFMAP = 3
  
      IF(ISXFMAP.EQ.1) THEN
* ... constant field, nothing to do
      ELSE IF(ISXFMAP.LE.3) THEN
* ... constant mesh field
         PRINT 10000, ISXFMAP
10000    FORMAT(' *SXFMAP* Magnetic field map flag: ',I5
     $        ,'; Reading magnetic field map data ')
*
         OPEN(77,FILE='/home/morsch/AliRoot/V3.02/data/field01.dat',
     $        FORM='FORMATTED',STATUS='OLD')
         READ(77,*) NX,NY,NZ,DX,DY,DZ,XMBEG,YMBEG,ZMBEG
         PRINT*,'NX,NY,NZ,DX,DY,DZ,XMBEG,YMBEG,ZMBEG',
     $   NX,NY,NZ,DX,DY,DZ,XMBEG,YMBEG,ZMBEG
         IF(3*NX*NY*NZ.GT.MAXFLD) THEN
            WRITE(6,10100) 3*NX*NY*NZ,MAXFLD
            STOP 'Increase MAXFLD'
         ENDIF
         UDX=1/DX
         UDY=1/DY
         UDZ=1/DZ
         XMEND=XMBEG+(NX-1)*DX
         YMEND=YMBEG+(NY-1)*DY
         ZMEND=ZMBEG+(NZ-1)*DZ
         DO IZ=1,NZ
            IPZ=3*(IZ-1)*(NX*NY)
            DO IY=1,NY
               IPY=IPZ+3*(IY-1)*NX
               DO IX=1,NX
                  IPX=IPY+3*(IX-1)
                  READ(77,*) BV(IPX+3),BV(IPX+2),BV(IPX+1)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(77)
      ENDIF                     ! endif ISXFMAP
*
10100 FORMAT('*** SXFMAP: Need ',I7,' have ',I7)
      END

****************************************
      SUBROUTINE RANNOR (A,B)
****************************************
C
C CERN PROGLIB# V100    RANNOR          .VERSION KERNFOR  4.18  880425
C ORIG. 18/10/77
C
      Y = RNDM()
      IF (Y.EQ.0.)  Y = RNDM()
      Z = RNDM()

      X = 6.283185*Z
      R = SQRT (-2.0*LOG(Y))
      A = R*SIN (X)
      B = R*COS (X)
      RETURN
      END
************************************************************************  
      SUBROUTINE OLDFCNFIT(NPAR,GRAD,FVAL,XVAL,IFLAG,FUTIL)
************************************************************************ 
*    Calcule FVAL: la fonction minimisee par MINUIT
*    with constant magnetic Field
*      
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      PARAMETER(NPLANE=10)
            
      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
     
      COMMON/MEAS/LPLANE(NPLANE),XMP(NPLANE),YMP(NPLANE)
     
      COMMON/FCNOUT/PXZEA,ALAMEA,PHIEA,XEA,YEA,NPLU,CHI2
           
      DIMENSION GRAD(*),XVAL(*),AMS(500),DISTAZ(500)
  
      DIMENSION XP(NPLANE),YP(NPLANE),
     &          COV(NPLANE,NPLANE),AP(NPLANE),COVY(NPLANE,NPLANE)
      DIMENSION VECT(7),VOUT(7) 

      XV = XVAL(4)   
      YV = XVAL(5)   
              
      ASIGN = 1.
      IF (XVAL(1).LT.0.) ASIGN = -1.
      PHI  = XVAL(2)
      ALAM = XVAL(3)
      PXZ = DABS(1./XVAL(1))
           
      PX = PXZ*DSIN(PHI)
      PY = PXZ*DTAN(ALAM)
      PZ = PXZ*DCOS(PHI)
      PTOT = PXZ/DCOS(ALAM)
      TTHET = DSQRT(DTAN(ALAM)**2+DSIN(PHI)**2)/DCOS(PHI)
      THET = DATAN(TTHET)
      PT = DSQRT(PX**2+PY**2)
      
      RL3 = ASIGN*PT / (0.299792458D-3 * BL3)
      ALPHA = DATAN2(PY,PX)
      XC = XV+RL3*DSIN(ALPHA)
      YC = YV-RL3*DCOS(ALPHA)
      
      A12 = 0.
      NTMAX = 0
      
      ZEA = ZABS
      XEA = XV
      YEA = YV
      PXEA = PX
      PYEA = PY
      PHIEA = PHI
      PXZEA = ASIGN*PXZ
      ALAMEA = ALAM
* 1er plan      
      ANGDEV = (ZPLANEP(1)-ZEA)*TTHET/RL3 
      XP(1) = XC+RL3*DSIN(ANGDEV-ALPHA)
      YP(1) = YC+RL3*DCOS(ANGDEV-ALPHA)
      AL   = THICK/ DCOS(THET)
      AP(1)  = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
* 2eme plan      
      ANGDEV = (ZPLANEP(2)-ZEA)*TTHET/RL3 
      XP(2) = XC+RL3*DSIN(ANGDEV-ALPHA)
      YP(2) = YC+RL3*DCOS(ANGDEV-ALPHA)
      AL   = THICK/ DCOS(THET)
      AP(2) = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
* 3eme plan      
      ANGDEV = (ZPLANEP(3)-ZEA)*TTHET/RL3 
      XP(3) = XC+RL3*DSIN(ANGDEV-ALPHA)
      YP(3) = YC+RL3*DCOS(ANGDEV-ALPHA)
      AL   = THICK/ DCOS(THET)
      AP(3)  = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
* 4eme plan      
      ANGDEV = (ZPLANEP(4)-ZEA)*TTHET/RL3 
      XP(4) = XC+RL3*DSIN(ANGDEV-ALPHA)
      YP(4) = YC+RL3*DCOS(ANGDEV-ALPHA)
      AL   = THICK/ DCOS(THET)
      AP(4) = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
* Fin de L3     
      ANGDEV = (700.D0-ZEA)*TTHET/RL3 
      XPL3 = XC+RL3*DSIN(ANGDEV-ALPHA)
      YPL3 = YC+RL3*DCOS(ANGDEV-ALPHA)
      PX = PT*DCOS(ANGDEV-ALPHA)
      PY = -PT*DSIN(ANGDEV-ALPHA)
      PHIC = DATAN2(PY,PX)
      CX = DSIN(THET)*DCOS(PHIC) 
      CY = DSIN(THET)*DSIN(PHIC) 
      CZ = DCOS(THET)
* Entree du dipole      
      VECT(1) = XPL3+(ZMAGS-700.)*CX/CZ
      VECT(2) = YPL3+(ZMAGS-700.)*CY/CZ
      VECT(3) = ZMAGS
      VECT(4) = CX
      VECT(5) = CY
      VECT(6) = CZ
      VECT(7) = PTOT
      
      PXZ  = PTOT*DSQRT(VECT(4)**2+VECT(6)**2)
      RDIP = ASIGN*PXZ / (0.299792458D-3 * B0)
      PHI1 = DATAN2(VECT(4),VECT(6)) 
      XC   = VECT(1)+RDIP*DCOS(PHI1)
      ZC   = VECT(3)-RDIP*DSIN(PHI1)
      
      IF (DABS((ZMAGE-ZPLANEP(5))/RDIP).GE.1.D0) RETURN ! Particule boucle dans le champ
      XP(5) = XC-ASIGN*DSQRT(RDIP**2-(ZPLANEP(5)-ZC)**2)
      YP(5) = VECT(2)+(ZPLANEP(5)-ZMAGS)*VECT(5)/VECT(6)
      CX = (ZPLANEP(5)-ZC)/RDIP
      CY = VECT(5)
      CZ = DABS((XP(5)-XC)/RDIP)
      THET = DATAN2(DSQRT(CX**2+CY**2),CZ)
      AL   = THICK/ DCOS(THET)
      AP(5) = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
      
      IF (DABS((ZMAGE-ZPLANEP(6))/RDIP).GE.1.D0) RETURN ! Particule boucle dans le champ
      XP(6) = XC-ASIGN*DSQRT(RDIP**2-(ZPLANEP(6)-ZC)**2)
      YP(6) = VECT(2)+(ZPLANEP(6)-ZMAGS)*VECT(5)/VECT(6)
      CX = (ZPLANEP(6)-ZC)/RDIP
      CY = VECT(5)
      CZ = DABS((XP(6)-XC)/RDIP)
      THET = DATAN2(DSQRT(CX**2+CY**2),CZ)
      AL   = THICK/ DCOS(THET)
      AP(6) = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
      
      
      IF (DABS((ZMAGE-ZC)/RDIP).GE.1.D0) RETURN ! Particule boucle dans le champ
      VOUT(1) = XC-ASIGN*DSQRT(RDIP**2-(ZMAGE-ZC)**2)        
      VOUT(2) = VECT(2)+(ZMAGE-ZMAGS)*VECT(5)/VECT(6)
      VOUT(3) = ZMAGE
      VOUT(4) = (ZMAGE-ZC)/RDIP
      VOUT(5) = VECT(5)
      VOUT(6) = DABS((VOUT(1)-XC)/RDIP)
      VOUT(7) = PTOT
      
      DO IV = 1,7
         VECT(IV) = VOUT(IV)
      ENDDO
       
* 7eme plan      
      THET = DATAN2(DSQRT(VECT(4)**2+VECT(5)**2),VECT(6))
      XP(7) = VECT(1)+(ZPLANEP(7)-ZMAGE)*VECT(4)/VECT(6)
      YP(7) = VECT(2)+(ZPLANEP(7)-ZMAGE)*VECT(5)/VECT(6)
      AL   = THICK/ DCOS(THET)
      AP(7)  = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
* 8eme plan      
      XP(8) = XP(7)+(ZPLANEP(8)-ZPLANEP(7))*VECT(4)/VECT(6)
      YP(8) = YP(7)+(ZPLANEP(8)-ZPLANEP(7))*VECT(5)/VECT(6)
      AL   = THICK/ DCOS(THET)
      AP(8)  = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))

* 9eme plan      
      XP(9) = XP(8)+(ZPLANEP(9)-ZPLANEP(8))*VECT(4)/VECT(6)
      YP(9) = YP(8)+(ZPLANEP(9)-ZPLANEP(8))*VECT(5)/VECT(6)
      AL   = THICK/ DCOS(THET)
      AP(9)  = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))
      
* 10eme plan      
      XP(10) = XP(9)+(ZPLANEP(10)-ZPLANEP(9))*VECT(4)/VECT(6)
      YP(10) = YP(9)+(ZPLANEP(10)-ZPLANEP(9))*VECT(5)/VECT(6)
      AL   = THICK/ DCOS(THET)
      AP(10)  = (0.0136D0/PTOT) * DSQRT(AL) * (1 + 0.038D0*DLOG(AL))

** Matrice de covariance      
      I = 0
      DO II = 1,NPLANE
        IF (LPLANE(II).EQ.1) THEN
        I = I + 1
*        I = II
        J = I - 1
        DO JJ = II, NPLANE
           IF (LPLANE(JJ).EQ.1) THEN
           J = J + 1
*           J = JJ
           COV (I,J) = 0.0D0
           COV (J,I) = A12
           IF (I .EQ. J) THEN
                 COV(J,I) =COV(J,I) + XPREC**2
           ENDIF      
                   
*           IF (I .EQ. 10 .AND. J .EQ. 10) PRINT *,'10 10   ',COV(J,I)
           DO L = 1,NTMAX
              COV(J,I) = COV(J,I)
     &  +  (ZPLANEP(II) + DISTAZ(L))*(ZPLANEP(JJ) + DISTAZ(L))*AMS(L)**2
           ENDDO
           DO K = 1, II-1
              COV(J,I) = COV(J,I)
     &  + (ZPLANEP(II)-ZPLANEP(K))*(ZPLANEP(JJ)-ZPLANEP(K))*AP(K)**2
*              IF (I .EQ. 10 .AND. J .EQ. 10) PRINT *,'10 10   ',COV(J,I)
           ENDDO
           COVY(I,J) = 0.0D0
           COVY(J,I) = COV(J,I)
           IF (I .EQ. J) THEN
                 COVY(J,I) = COVY(J,I) - XPREC**2 + YPREC**2
           ENDIF
        ENDIF   
        ENDDO
        ENDIF
      ENDDO
 
*  Inversion des matrices de covariance
      NPLU = I
 
      IFAIL = 0
      CALL DSINV(NPLU, COV, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'ERROR'
      IFAIL = 0
      CALL DSINV(NPLU, COVY, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'ERROR'
*      PRINT *,' COVARIANCE MATRIX AFTER'
*      DO I = 1, NPLANE
*         PRINT *,(COV(J,I),J=1,NPLANE)
*      ENDDO
 
** Calcul de FVAL ou CHI2
      FVAL = 0.0D0
      I = 0
      DO II = 1,NPLANE
        IF (LPLANE(II).EQ.1) THEN
        I = I+1
*        I = II
        J = 0
        DO JJ = 1,NPLANE
           IF (LPLANE(JJ).EQ.1) THEN
              J = J+1
*             J = JJ
              FVAL = FVAL + COV(J,I)*(XMP(II)-XP(II))*(XMP(JJ)-XP(JJ))
              FVAL = FVAL + COVY(J,I)*(YMP(II)-YP(II))
     &                               *(YMP(JJ)-YP(JJ))
**             IF (JJ.EQ.II) THEN
**                 FVAL = FVAL + (XMP(II)-XP(II))*(XMP(JJ)-XP(JJ))/XPREC**2
**                 FVAL = FVAL + (YMP(II)-YP(II))
**     &                               *(YMP(JJ)-YP(JJ))/YPREC**2
**             ENDIF
           ENDIF
        ENDDO
        ENDIF
      ENDDO
      CHI2 = FVAL
         print *,' fcnfit pxz tphi talam xvert yvert chi2',
     &       PXZEA,PHIEA,ALAMEA,
     &       XEA,YEA,CHI2/FLOAT(2*NPLU-5)  
      
 1000 FORMAT(I5,7F12.6)
 
      RETURN
      END

*
* $Id$
*
* $Log$
* Revision 1.4.4.2  2000/04/26 15:48:37  morsch
* Some routines from obsolete algo.F are needed by reco_muon.F and have been
* copied there.
*
* Revision 1.4.4.1  2000/01/12 16:00:55  morsch
* New version of MUON code
*
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE RECO_GHELIX (CHARGE, STEP, VECT, VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Performs the tracking of one step in a magnetic field         *
C.    *  The trajectory is assumed to be a helix in a constant field   *
C.    *  taken at the mid point of the step.                           *
C.    *  Parameters:                                                   *
C.    *   input                                                        *
C.    *     STEP =arc length of the step asked                         *
C.    *     VECT =input vector (position,direction cos and momentum)   *
C.    *     CHARGE=  electric charge of the particle                   *
C.    *   output                                                       *
C.    *     VOUT = same as VECT after completion of the step           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Author    M.Hansroul  *********                          *
C.    *       Modified  S.Egli, S.V.Levonian                           *
C.    *       Modified  V.Perevoztchikov
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION      VECT(7),VOUT(7)
      DIMENSION      XYZ(3),H(4),HXP(3)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6,IPP=7)
      PARAMETER      (SIXTH = 1./6.)
      PARAMETER      (EC=2.9979251E-4)
C.
C.    ------------------------------------------------------------------
C.
C       units are kgauss,centimeters,gev/c
C
      VOUT(IPP) = VECT(IPP)
      IF (CHARGE.EQ.0.)         GO TO 10
      XYZ(1)    = VECT(IX) + 0.5 * STEP * VECT(IPX)
      XYZ(2)    = VECT(IY) + 0.5 * STEP * VECT(IPY)
      XYZ(3)    = VECT(IZ) + 0.5 * STEP * VECT(IPZ)
C
      CALL RECO_GUFLD (XYZ, H)
 
      H2XY = H(1)**2 + H(2)**2
      H(4) = H(3)**2 + H2XY
      IF (H(4).LE.1.E-12)       GO TO 10
      IF (H2XY.LE.1.E-12*H(4))  THEN
         CALL RECO_GHELX3 (CHARGE*H(3), STEP, VECT, VOUT)
         GO TO 999
      ENDIF
      H(4) = SQRT(H(4))
      H(1) = H(1)/H(4)
      H(2) = H(2)/H(4)
      H(3) = H(3)/H(4)
      H(4) = H(4)*EC
*
      HXP(1) = H(2)*VECT(IPZ) - H(3)*VECT(IPY)
      HXP(2) = H(3)*VECT(IPX) - H(1)*VECT(IPZ)
      HXP(3) = H(1)*VECT(IPY) - H(2)*VECT(IPX)
 
      HP = H(1)*VECT(IPX) + H(2)*VECT(IPY) + H(3)*VECT(IPZ)
*
      RHO = -CHARGE*H(4)/VECT(IPP)
      TET = RHO * STEP
      IF (ABS(TET).GT.0.15)     THEN
         SINT = SIN(TET)
         SINTT = (SINT/TET)
         TSINT = (TET-SINT)/TET
         COS1T = 2.*(SIN(0.5*TET))**2/TET
      ELSE
         TSINT = SIXTH*TET**2
         SINTT = (1. - TSINT)
         SINT = TET*SINTT
         COS1T = 0.5*TET
      ENDIF
*
      F1 = STEP * SINTT
      F2 = STEP * COS1T
      F3 = STEP * TSINT * HP
      F4 = -TET*COS1T
      F5 = SINT
      F6 = TET * COS1T * HP
 
      VOUT(IX) = VECT(IX) + (F1*VECT(IPX) + F2*HXP(1) + F3*H(1))
      VOUT(IY) = VECT(IY) + (F1*VECT(IPY) + F2*HXP(2) + F3*H(2))
      VOUT(IZ) = VECT(IZ) + (F1*VECT(IPZ) + F2*HXP(3) + F3*H(3))
 
      VOUT(IPX) = VECT(IPX) + (F4*VECT(IPX) + F5*HXP(1) + F6*H(1))
      VOUT(IPY) = VECT(IPY) + (F4*VECT(IPY) + F5*HXP(2) + F6*H(2))
      VOUT(IPZ) = VECT(IPZ) + (F4*VECT(IPZ) + F5*HXP(3) + F6*H(3))
 
      GO TO 999
 
   10 CONTINUE
      DO 20 I   = 1,3
         VOUT(I) = VECT(I) + STEP * VECT(I+3)
         VOUT(I+3) = VECT(I+3)
   20 CONTINUE
C
  999 END

*
* $Id$
*
* $Log$
* Revision 1.4.4.2  2000/04/26 15:48:37  morsch
* Some routines from obsolete algo.F are needed by reco_muon.F and have been
* copied there.
*
* Revision 1.4.4.1  2000/01/12 16:00:55  morsch
* New version of MUON code
*
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*

*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE RECO_GHELX3 (FIELD, STEP, VECT, VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Tracking routine in a constant field oriented            *
C.    *       along axis 3                                             *
C.    *       Tracking is performed with a conventional                *
C.    *       helix step method                                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Authors    R.Brun, M.Hansroul  *********                 *
C     *       Rewritten  V.Perevoztchikov
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION      VECT(7),VOUT(7),HXP(3)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6,IPP=7)
      PARAMETER      (SIXTH = 1./6.)
      PARAMETER      (EC=2.9979251E-4)
C.
C.    ------------------------------------------------------------------
C.
C       units are kgauss,centimeters,gev/c
C
      VOUT(IPP) = VECT(IPP)
      H4 = FIELD * EC
*
      HXP(1) = - VECT(IPY)
      HXP(2) = + VECT(IPX)
 
      HP = VECT(IPZ)
*
      RHO = -H4/VECT(IPP)
      TET = RHO * STEP
      IF (ABS(TET).GT.0.15)     THEN
         SINT = SIN(TET)
         SINTT = (SINT/TET)
         TSINT = (TET-SINT)/TET
         COS1T = 2.*(SIN(0.5*TET))**2/TET
      ELSE
         TSINT = SIXTH*TET**2
         SINTT = (1. - TSINT)
         SINT = TET*SINTT
         COS1T = 0.5*TET
      ENDIF
*
      F1 = STEP * SINTT
      F2 = STEP * COS1T
      F3 = STEP * TSINT * HP
      F4 = -TET*COS1T
      F5 = SINT
      F6 = TET * COS1T * HP
 
      VOUT(IX) = VECT(IX) + (F1*VECT(IPX) + F2*HXP(1))
      VOUT(IY) = VECT(IY) + (F1*VECT(IPY) + F2*HXP(2))
      VOUT(IZ) = VECT(IZ) + (F1*VECT(IPZ) + F3)
 
      VOUT(IPX) = VECT(IPX) + (F4*VECT(IPX) + F5*HXP(1))
      VOUT(IPY) = VECT(IPY) + (F4*VECT(IPY) + F5*HXP(2))
      VOUT(IPZ) = VECT(IPZ) + (F4*VECT(IPZ) + F6)
 
C
  999 END

************************************************************************
      DOUBLE PRECISION FUNCTION RECOCHI2 (MPOS,MANG,XM,YM,ALAMM,APHIM,
     &     XC,YC,ALAMC,APHIC,PTOT,IZMES,NPLPL)
C.
C.    ******************************************************************
C.    *   Calculate chi2 taking into account MSC                       *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NBSTATION=5,NPLANE=10)

      COMMON/ZDEFIN/ZPLANE(NBSTATION),ZCOIL,ZMAGEND,DZ_PL(NBSTATION)

      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST
C
      DIMENSION MPOS(NBSTATION),MANG(NBSTATION),
     &    XM(NBSTATION),YM(NBSTATION),ALAMM(NBSTATION),APHIM(NBSTATION),
     &    XC(NBSTATION),YC(NBSTATION),ALAMC(NBSTATION),APHIC(NBSTATION),
     &    ALAMP(NBSTATION),APHIP(NBSTATION)

      DIMENSION AP(NPLANE),IZMES(NBSTATION),IPLANE(NPLANE)  
      DIMENSION COVX(NPLANE,NPLANE),COVY(NPLANE,NPLANE)
      DIMENSION XM1(NPLANE),YM1(NPLANE),XC1(NPLANE),YC1(NPLANE)


      ICH = 0 
      NPLPL = 0 
      DO IZ = 1,NBSTATION    
          ICH = ICH + 1 
          IPLANE(ICH) = 0   
          AL   = THICK/(COS(ALAMC(IZ))*COS(APHIC(IZ)))
          AP(ICH) = (0.0136D0/PTOT)*DSQRT(AL)*(1.+0.038D0*DLOG(AL))
          IF (MPOS(IZ).EQ.1.AND.IZMES(IZ).EQ.1) THEN
             IPLANE(ICH) = 1
             XM1(ICH) = XM(IZ)
             YM1(ICH) = YM(IZ)
          ENDIF
          XC1(ICH) = XC(IZ)
          YC1(ICH) = YC(IZ)

          ICH = ICH + 1    
          IPLANE(ICH) = 0   
          AP(ICH) = (0.0136D0/PTOT)*DSQRT(AL)*(1.+0.038D0*DLOG(AL))
          IF (MPOS(IZ).EQ.1.AND.IZMES(IZ).EQ.2) THEN
             IPLANE(ICH) = 1
             XM1(ICH) = XM(IZ)
             YM1(ICH) = YM(IZ)
          ENDIF
          IF (MANG(IZ).EQ.1) THEN
             IPLANE(ICH) = 1
             XM1(ICH) = XM(IZ) - DZ_PL(IZ) * TAN(APHIM(IZ))
             YM1(ICH) = YM(IZ) + DZ_PL(IZ)/COS(APHIM(IZ))*TAN(ALAMM(IZ))
          ENDIF 
          XC1(ICH) = XC(IZ) - DZ_PL(IZ) * TAN(APHIC(IZ))
          YC1(ICH) = YC(IZ) + DZ_PL(IZ)/COS(APHIC(IZ))*TAN(ALAMC(IZ))
          IF (IPLANE(ICH).EQ.1) NPLPL = NPLPL+1          
      ENDDO



** Matrice de covariance X et Y   
      I = 0
      DO II = 1,NPLANE

        IF (IPLANE(II).EQ.1) THEN
           I = I + 1

           J = I - 1
           DO JJ = II, NPLANE
              IF (IPLANE(JJ).EQ.1) THEN
                 J = J + 1

                 COVX (I,J) = 0.
                 COVX (J,I) = 0.
                 IF (I .EQ. J) THEN
                    COVX(J,I) =COVX(J,I) + XPREC**2
                 ENDIF      

                 DO K = 1, II-1
                    COVX(J,I) = COVX(J,I)
     &                   +  (-ZPLANEP(II)+ZPLANEP(K))*
     &                   (-ZPLANEP(JJ)+ZPLANEP(K))*AP(K)**2

                 ENDDO
                 COVY(I,J) = 0.
                 COVY(J,I) = COVX(J,I)
                 IF (I .EQ. J) THEN
                    COVY(J,I) = COVY(J,I) - XPREC**2 + YPREC**2
                 ENDIF
              ENDIF   
           ENDDO
        ENDIF
      ENDDO

      NPLUP = I

 
*  Inversion des matrices de covariance
      IFAIL = 0
      CALL DSINV(NPLUP, COVX, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'RECOCHI2 ERROR COVX'

      IFAIL = 0
      CALL DSINV(NPLUP, COVY, NPLANE, IFAIL)
**      IF (JFAIL.NE.0 .AND. IFAIL .NE. 0) STOP 'ERROR'
      IF (IFAIL .NE. 0) STOP 'RECOCHI2 ERROR COVY'

** Calcul de FVAL ou CHI2
      FVAL = 0.

      I = 0
      DO II = 1,NPLANE
        IF (IPLANE(II).EQ.1) THEN
           I = I+1
           J = 0
**           print *,' II=',ii,' XM1,YM1 =',xm1(ii),ym1(ii),
**     &          ' XC,YC =',xc1(ii),yc1(ii)
           DO JJ = 1,NPLANE
              IF (IPLANE(JJ).EQ.1) THEN
                 J = J+1
                 FVAL = FVAL + COVX(J,I)*(XM1(II)-XC1(II))
     &                *(XM1(JJ)-XC1(JJ))
                 FVAL = FVAL + COVY(J,I)*(YM1(II)-YC1(II))
     &                *(YM1(JJ)-YC1(JJ))
              ENDIF
           ENDDO
        ENDIF
      ENDDO

      RECOCHI2 = FVAL
**      print *,' recochi2 =',recochi2  
        
C
      RETURN
      END

************************************************************************        
      SUBROUTINE RECO_SELECT(ISEL)
************************************************************************       *    ISEL(I) = 1 if track number I is OK, ISEL(I) = 0 otherwise  
************************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NPLANE=10,NBCHAMBER=10,NTRMAX=500)

      COMMON/DEBEVT/IDEBUG

      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &             ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST

      COMMON/TRACKFOUT/IEVOUT,NTREVT,JJOUT(NBCHAMBER,NTRMAX),
     &               ISTAT(NTRMAX),PXZOUT(NTRMAX),TPHIOUT(NTRMAX),
     &               TALAMOUT(NTRMAX),XVERTOUT(NTRMAX),YVERTOUT(NTRMAX),
     &               CHI2OUT(NTRMAX),
     &               XMESOUT(NBCHAMBER,NTRMAX),YMESOUT(NBCHAMBER,NTRMAX)  
     &              ,PXVOUT(NTRMAX),PYVOUT(NTRMAX),PZVOUT(NTRMAX)

      DIMENSION ISEL(NTRMAX)

      DO I = 1,NTREVT
         ISEL(I) = 1
      ENDDO   
 
      DO I = 1,NTREVT
      
         ICH1 = 9
         JJ1 = JJOUT(ICH1,I)
         IF (JJ1.EQ.0) ICH1 = 10
         X1 = XMESOUT(ICH1,I)
         Y1 = YMESOUT(ICH1,I)

         DO J = I+1,NTREVT
            ICH2 = 9
            JJ2 = JJOUT(ICH2,J)
            IF (JJ2.EQ.0) ICH2 = 10
            X2 = XMESOUT(ICH2,J)
            Y2 = YMESOUT(ICH2,J)
            DIST = SQRT(((X2-X1)/(10.*XPREC))**2
     &           +((Y2-Y1)/(10.*YPREC))**2)
            IF (DIST.LT.2.) THEN
               CHI21 =  CHI2OUT(I)
               CHI22 =  CHI2OUT(J)
               IF (CHI22.LT.CHI21) THEN
                  ISEL(I) = 0
                  IF (IDEBUG.EQ.2) THEN 
                     PRINT *,' RECO_SELECT I,ISEL= ',I,ISEL(I)
                  ENDIF
               ELSE
                  ISEL(J) = 0
                  IF (IDEBUG.EQ.2) THEN
                     PRINT *,' RECO_SELECT J,ISEL= ',J,ISEL(J)
                  ENDIF   
               ENDIF   
            ENDIF
         ENDDO
 
      ENDDO   
            

      RETURN
      END

***************************************************
      SUBROUTINE reconstmuon2(IFIT,IDEBUGC,NEV)
***************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      PARAMETER(NPLANE=10,MAXIDG=28,NTRMAX=500)

      REAL*4 PXR,PYR,PZR,ZVR,CHI2R,PXV,PYV,PZV

      COMMON/TRACKFI/EFF,EFF1,EFF2,XPREC,YPREC,PHIPREC,ALAMPREC,
     &                 HCUT,LBKG,SIGCUT,ALPHATOP,HTOP

      COMMON/PAWCR4/IEVR,NTRACKR,ISTATR(NTRMAX),ISIGNR(NTRMAX),
     &              PXR(NTRMAX),PYR(NTRMAX),PZR(NTRMAX),ZVR(NTRMAX),
     &              CHI2R(NTRMAX),PXV(NTRMAX),PYV(NTRMAX),PZV(NTRMAX)

      COMMON/MEAS/LPLANE(NPLANE),XMP(NPLANE),YMP(NPLANE)

      COMMON/FIT/NHITTOT1,izch(maxidg),xgeant(maxidg),
     &     ygeant(nplane)

      COMMON/PRECCUT/PCUT,PTCUT,CHI2CUT

      REAL*4 RN1,RN2

** Read events          
      CALL trackf_read_fit(IEVR,NEV,NHITTOT1,IZCH,XGEANT,YGEANT)
      print*,'nhittot1 ',nhittot1

*      do i=1,NHITTOT1
*         print*,'x=',xgeant(i),' y=',ygeant(i),' ch=',izch(i)
*      enddo

      if (nhittot1.ne.20) goto 55
      nhit1=1
      nhit2=10
      do ntr=1,2                ! loop over tracks
         do nhit=nhit1,nhit2 ! loop over hits
            ich=izch(nhit)
            lplane(ich)=1
            
            CALL RANNOR(RN1,RN2) ! CCC
*            xmp(ich)=xgeant(nhit)
*            ymp(ich)=ygeant(nhit)
            xmp(ich)=xgeant(nhit) + RN1 * XPREC
            ymp(ich)=ygeant(nhit) + RN2 * YPREC
         end do
         call fit_trace(ntr)
         nhit1=11
         nhit2=20
      end do   

      CALL CHFNT(IEVR,NTRACK,ISTATR,ISIGNR,
     &     PXR,PYR,PZR,ZVR,CHI2R,PXV,PYV,PZV)
            
 55   continue

      END

*********************************
      subroutine fit_trace(ntr)
*********************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(NPLANE=10,NTRACK=2,NTRMAX=500)

      COMMON/FCNOUT/PXZEA,ALAMEA,PHIEA,XEA,YEA,NPLU,CHI2

      REAL*4 PXR,PYR,PZR,ZVR,CHI2R,PXV,PYV,PZV
      COMMON/PAWCR4/IEVR,NTRACKR,ISTATR(NTRMAX),ISIGNR(NTRMAX),
     &              PXR(NTRMAX),PYR(NTRMAX),PZR(NTRMAX),ZVR(NTRMAX),
     &              CHI2R(NTRMAX),PXV(NTRMAX),PYV(NTRMAX),PZV(NTRMAX)

      COMMON/MEAS/LPLANE(NPLANE),XMP(NPLANE),YMP(NPLANE)

      COMMON/PARAM/ZPLANEP(NPLANE),THICK,XPREC,YPREC,B0,BL3,ZMAGS,
     &     ZMAGE,ZABS,XMAG,ZBP1,ZBP2,CONST

      COMMON/PRECCUT/PCUT,PTCUT,CHI2CUT 
 
      dimension PXZOUT(NTRMAX),TPHIOUT(NTRMAX),TALAMOUT(NTRMAX),
     &     CHI2OUT(NTRMAX)

      x1=xmp(1)
      y1=ymp(1)
      ipl1=1
      x2=xmp(3)
      y2=ymp(3)
      ipl2=3
      x3=xmp(7)
      ipl3=7
      x4=xmp(9)
      ipl4=9

      PHIAV = DATAN2((X2-X1),(ZPLANEP(IPL2)-ZPLANEP(IPL1)))          
      PHIAP = DATAN2((X4-X3),(ZPLANEP(IPL4)-ZPLANEP(IPL3)))

      DPHI = (PHIAP-PHIAV)
      ASIGN = 1.
      IF (DPHI.LT.0.) ASIGN = -1. ! CCC
      PXZ = CONST/DABS(DPHI)

*     * Cuts PXZ           
      IF (PXZ.LT.PCUT) GO TO 66       
      
      PXZINVI = ASIGN/PXZ       ! CCC
      PHII = PHIAV
      ALAMI = DATAN2((Y2-Y1),DSQRT((X2-X1)**2
     &     +(ZPLANEP(IPL2)-ZPLANEP(IPL1))**2))
      XVR = X1
      YVR = Y1
      
*         print *,' avant prec_fit pxzi phii alami x y',1./ PXZINVI,
*     &         PHII, ALAMI ,XVR,YVR             
*         PRINT *,' X1 X2 X3 X4',X1,X2,X3,X4
*         PRINT *,' Z1 Z2 Z3 Z4',ZPLANEP(IPL1),ZPLANEP(IPL2),
*     &          ZPLANEP(IPL3),ZPLANEP(IPL4)
*         PRINT *,' CONST= ',CONST 

*     * Fit des traces apres l'absorbeur           
      CALL PREC_FIT (PXZINVI,PHII,ALAMI,XVR,YVR,
     &     PXZINVF,PHIF,ALAMF,XVERTF,YVERTF,EPXZINV,EPHI,EALAM,
     &     EXVERT,EYVERT)
      
*     * Correction de Branson       
      CALL BRANSON(PXZEA,PHIEA,ALAMEA,XEA,YEA)

      PXZ1 = DABS(PXZEA)
      PX1 = PXZ1*DSIN(PHIEA)
      PY1 = PXZ1*DTAN(ALAMEA)
      PT1 = DSQRT(PX1**2+PY1**2)
 
*      print*,'pt1=',pt1
*      print*,'ptcut=',ptcut
*      print*,'chi2=',CHI2/FLOAT(2*NPLU-5)
*      print*,'chi2cut=',CHI2CUT

*     * Cuts PT
      IF (PT1.LT.PTCUT) GO TO 66
*     * Cuts CHI2
      IF ((CHI2/FLOAT(2*NPLU-5)).GT.CHI2CUT) GO TO 66

      PXZOUT(NTR) = PXZEA
      TPHIOUT(NTR) = DTAN(PHIEA)
      TALAMOUT(NTR) = DTAN(ALAMEA) 
      CHI2OUT(NTR) =  CHI2/FLOAT(2*NPLU-5) 

      ISIGNR(NTR) = 1 
      IF (PXZOUT(NTR).LT.0.) ISIGNR(NTR) = -1
      PXZ = ABS(PXZOUT(NTR))
      PHI = ATAN(TPHIOUT(NTR)) 
      ALAM = ATAN(TALAMOUT(NTR)) 
      PYR(NTR) = PXZ*SIN(PHI)
      PXR(NTR) = PXZ*TAN(ALAM)
      PZR(NTR) = PXZ*COS(PHI)

 66   CONTINUE 
      
      return
      end

      SUBROUTINE SORTZV (A,INDEX,N1,MODE,NWAY,NSORT)
C
C CERN PROGLIB# M101    SORTZV          .VERSION KERNFOR  3.15  820113
C ORIG. 02/10/75
C
      DIMENSION A(N1),INDEX(N1)
C
C
      N = N1
      IF (N.LE.0)            RETURN
      IF (NSORT.NE.0) GO TO 2
      DO 1 I=1,N
    1 INDEX(I)=I
C
    2 IF (N.EQ.1)            RETURN
      IF (MODE)    10,20,30
   10 CALL SORTTI (A,INDEX,N)
      GO TO 40
C
   20 CALL SORTTC(A,INDEX,N)
      GO TO 40
C
   30 CALL SORTTF (A,INDEX,N)
C
   40 IF (NWAY.EQ.0) GO TO 50
      N2 = N/2
      DO 41 I=1,N2
      ISWAP = INDEX(I)
      K = N+1-I
      INDEX(I) = INDEX(K)
   41 INDEX(K) = ISWAP
   50 RETURN
      END

      SUBROUTINE SORTTI (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END

*     ========================================
      SUBROUTINE SORTTC (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF(ICMPCH(AI,A(I22)))3,3,21
   21 INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (ICMPCH(A(I22),A(I222))) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (ICMPCH(AI,A(I22))) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      FUNCTION ICMPCH(IC1,IC2)
C     FUNCTION TO COMPARE TWO 4 CHARACTER EBCDIC STRINGS - IC1,IC2
C     ICMPCH=-1 IF HEX VALUE OF IC1 IS LESS THAN IC2
C     ICMPCH=0  IF HEX VALUES OF IC1 AND IC2 ARE THE SAME
C     ICMPCH=+1 IF HEX VALUES OF IC1 IS GREATER THAN IC2
      I1=IC1
      I2=IC2
      IF(I1.GE.0.AND.I2.GE.0)GOTO 40
      IF(I1.GE.0)GOTO 60
      IF(I2.GE.0)GOTO 80
      I1=-I1
      I2=-I2
      IF(I1-I2)80,70,60
 40   IF(I1-I2)60,70,80
 60   ICMPCH=-1
      RETURN
 70   ICMPCH=0
      RETURN
 80   ICMPCH=1
      RETURN
      END

      SUBROUTINE SORTTF (A,INDEX,N1)
C
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
