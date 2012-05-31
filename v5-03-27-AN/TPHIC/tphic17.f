CDECK  ID>, GGCDES. 
*=======================================================================
CDECK  ID>, GGINIT. 
      SUBROUTINE gginit
*-----------------------------------------------------------------------
*       GGHIC initialization
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
     &               amin,amax,ymin,ymax,nmas,ny, kferm,
     &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
     &               thetamin, costhv1, kv1,kv2,gvpar(4)
      CHARACTER lumfil*80
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
     &  xstote, ssbr(10)
      COMMON /ggevnt/ nrun,ievent,wsq,ygg,xmg1,xmg2, p2g(5),
     &                ptag1(4),ptag2(4), ngg, kgg(10),pgg(20,5)
      COMMON /ggmssm/ xm1,   xm2,    xmg,    xms,    xmtl,   xmtr,
     &                xmll,  xmlr,   xmnl,   xtanb,  xmha,   xmu,
     &                xmt,   xat,    xmbr,   xab,    u11,    v11
      COMMON/D2LParam/lz,la,Rion,Gamma
      INTEGER lz,la
      REAL    Rion,Gamma
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      DOUBLE PRECISION CKIN
      SAVE  /PYSUBS/
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      DOUBLE PRECISION PARP,PARI
      SAVE /PYPARS/
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      DOUBLE PRECISION P,V
      SAVE  /PYJETS/
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      DOUBLE PRECISION PARU,PARJ
      SAVE  /PYDAT1/
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      DOUBLE PRECISION PMAS,PARF,VCKM
      SAVE  /PYDAT2/
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      SAVE  /PYDAT4/
C          MXSS                 = maximum number of modes
C          NSSMOD               = number of modes
C          ISSMOD               = initial particle
C          JSSMOD               = final particles
C          GSSMOD               = width
C          BSSMOD               = branching ratio
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS)
      INTEGER NSSMOD,ISSMOD,JSSMOD
      REAL GSSMOD,BSSMOD
      SAVE /SSMODE/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMQKSS               = squark mass
C          AMLLSS               = left-slepton mass
C          AMLRSS               = right-slepton mass
C          AMNLSS               = sneutrino mass
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
      COMMON/SSPAR/AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
      REAL xpar(10),ypar(6)
      INTEGER ipar(2)
      INTEGER pycomp
      REAL*8 ggrnd
      EXTERNAL pydata
      la = ia
      lz = iz
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*       Check bounds on 2 photon mass
*
      IF (iproc .EQ. 1) THEN
        IF (amin .LE. 2.31) THEN
          WRITE (6,*)
     &      '  %GGINIT: minimal 2-photon mass is too low, ',
     &      'set to 3m(pho)'
          amin = 2.31
        END IF
*
      ELSE IF (iproc .EQ. 2) THEN
        kc = pycomp (kf_onium)
        xmres = pmas(kc,1)
        IF (kf_onium .EQ. 441) THEN
          xgtres = 0.013
        ELSE
          xgtres = pmas(kc,2)
        END IF
c.        xggres = ggwid(kf_onium)
        IF (xgtres .EQ. 0.) xgtres = 0.1
        IF (xmres.LE.amin .OR. xmres.GE.amax) THEN
          WRITE(6,*)
     &      '  %GGINIT: Xmres is out of the Luminosity mass range'
          STOP
        ENDIF
*
      ELSE IF (iproc .EQ. 3) THEN
        kc = pycomp (kferm)
        amprt(1) = pmas(kc,1)
        IF (amin .LT. 2*amprt(1)) THEN
          WRITE (6,*)
     &      '  %GGINIT: minimal 2-photon mass is too low, ',
     &      'set to 2m(fermion)'
          amin = 2*amprt(1)
        END IF
*
      ELSE IF (iproc .EQ. 4) THEN
        IF (amprt(4) .LE. 0.) amprt(4) = pmas(pycomp(24),1)
        IF (amin .LT. 2*amprt(4)) then
          WRITE (6,*)
     &      '  %GGINIT: minimal 2-photon mass is too low, ',
     &      'set to 2m(W+-)'
          amin = 2*amprt(4)
        END IF
*
      ELSE IF (iproc .EQ. 5) THEN
         WRITE (6,*) 'Process 5 is not implemented yet'
         STOP
C         CALL ssmssm (xm1, xm2, xmg, xms, xmtl, xmtr, xmll, xmlr, xmnl,
C      &               xtanb, xmha, xmu, xmt, xat, xmbr, xab, iallow)
C         IF (iallow .NE. 0) THEN
C           WRITE (6,*) 'Lightest neutralino is not LSP, ',
C      &                'input other MSSM parameters'
C           STOP
C         END IF
C         amprt(1)  = abs (amw1ss)
C         amprt(2)  = abs (amz1ss)
C         WRITE (6,'(''   %GGINIT'','//
C      &            ''', M(chi+) = '',f5.1,'//
C      &            ''', M(chi0) = '',f5.1)') amprt(1), amprt(2)
C *
C         IF (amin .LT. 2*amprt(1)) THEN
C           WRITE (6,*)
C      &      '  %GGINIT: minimal 2-photon mass is too low, ',
C      &      'set to 2M(chi+)'
C           amin = 2*amprt(1)
C         END IF
*
      ELSE IF (iproc .EQ. 6) THEN
        kc1 = pycomp (kv1)
        amprt(1) = pmas(kc1,1) + pmas(kc1,3)
        kc2 = pycomp (kv2)
        amprt(2) = pmas(kc2,1) + pmas(kc2,3)
        IF (amin .LT. amprt(1)+amprt(2)) THEN
          WRITE (6,*)
     &      '  %GGINIT: minimal 2-photon mass is too low, ',
     &      'set to m(kv1)+m(kv2)'
          amin = amprt(1)+amprt(2)
        END IF
*
      END IF
*
      IF (amin .GT. amax) THEN
        WRITE (6,*) '  %GGINIT: AMIN > AMAX, execution stopped'
        STOP
      END IF
*
*  -----  TWOGAM Initialization  -----
*
      gamma = ebmn*ia/amas
      xpar(1)= ebmn*ia  !  Beam energy
      xpar(2)= amin     !  Minimum gg mass
      xpar(3)= amax     !  Maximum gg mass
      xpar(8)= 100.     !  Weight, but not too long
*
      ipar(1)= 1        !  Unweighted events
      ipar(2)= 0        !  Continuum generation
*
      IF (iproc .EQ. 2) THEN
        xpar(9) = xmres     !  Resonance mass
        xpar(10)= xgtres    !  Resonance total width
        ipar(2) = 1         !  Resonance generation
      END IF
*
      ypar(1)= iz       !  Ion charge
      ypar(2)= ia       !  Ion atomic number
      ypar(3)= amas     !  Ion mass
      ypar(4)= gamma    !  Gamma factor of ion
      ypar(5)= ymin     !  Minimum gg-rapidity
      ypar(6)= ymax     !  Maximum gg-rapidity
      eb = xpar(1)
*
      CALL twoini (6,ipar,xpar,ypar)
*
*       The Luminosity function initialization
*
*     ilumf = +1         !   Reading the Luminosity from file
*     ilumf = -1         !   New calculation of the Luminosity function
*
      CALL lumfun (ilumf,lumfil,ymin,ymax,ny,amin,amax,nmas)
*
*       Integral of lum. function over the (M,y) region
      xlumint = 0.0
      hmas = (amax-amin)/nmas
      f1   = dldmx(amin   ,ny)
      DO i=1,nmas
        xm  = amin + i*hmas
        f2  = dldmx(xm-hmas/2.,ny)
        f3  = dldmx(xm        ,ny)
        xlumint = xlumint +(f1 + 4.*f2 + f3)/6.
        f1 = f3
      END DO
      xlumint = xlumint * hmas
*  ----- end of TWOGAM initialization -----
*
      IF (iproc .EQ. 1) THEN
*  ----- PYTHIA initialization -----
        parp(2) = dble(amin)
        mstp(14) = 20
        mstp(82) = 1
        mstp(171) = 1
        mstp(172) = 1
        DO ip = 1,2
          p(1,ip) = 0.D0
          p(2,ip) = 0.D0
        END DO
        p(1,3) = dble( amax/2)
        p(2,3) = dble(-amax/2)
        p(1,4) = dble( amax/2)
        p(2,4) = dble( amax/2)
        p(1,5) =  0.D0
        p(2,5) =  0.D0
        CALL pyinit ('5MOM','gamma','gamma',0.0D0)
*  ----- end of PYTHIA initialization -----
*
      ELSE IF (iproc .EQ. 2) THEN
*   ---  Cross section for narrow Resonance production  ---
        totxsc = dldmx(xmres,ny)                    !  dL/dMx in GeV^-1
        njres = kf_onium / 10
        njres = kf_onium - 10 * njres
        xnorm = 4. * pi**2 * njres * xggres * gev2nb / xmres**2
        xstot = xnorm * totxsc
        xstote = 0.
*
      ELSE IF      (iproc .EQ. 3) THEN
        qf = abs(kchg(kc,1))/3.
        nc = abs(kchg(kc,2))*2. + 1.
      ELSE IF (iproc .EQ. 5) THEN
*           Chargino storing in JETSET
        kc = pycomp (41)
        chaf(kc,1) = 'chi_1'
        kchg(kc,1) = -3
        kchg(kc,2) =  0
        kchg(kc,3) =  1
        pmas(kc,1) = amprt(1)
        qf = 1.
        nc = 1.
*           Neutralino storing in JETSET
        kc = pycomp (42)
        chaf(kc,1) = 'chi_1'
        kchg(kc,1) =  0
        kchg(kc,2) =  0
        kchg(kc,3) =  1
        pmas(kc,1) = amprt(2)
*
*           Gaugino mixing parameters
        u11 = -zmixss(4,1)*sin(gammar)/sqrt(2.) +
     &         zmixss(2,1)*cos(gammar)
        v11 =  zmixss(3,1)*sin(gammal)/sqrt(2.) +
     &         zmixss(2,1)*cos(gammal)
        DO idc = 1,nssmod
          IF (issmod(idc) .EQ. 39) THEN
            IF (moddcy .EQ. 1) THEN
              DO jmod = 1,5
                IF (jssmod(jmod,idc) .EQ. -12) THEN
                  xsbra = bssmod(idc)
                  xsbra = xsbra*xsbra*4
                  GOTO 10
                END IF
              END DO
            ELSE
              DO jmod = 1,5
                IF (jssmod(jmod,idc) .EQ. -12) THEN       ! e nu_e
                  ssbr(1) = bssmod(idc)
                  GOTO 20
                ELSE IF (jssmod(jmod,idc) .EQ. -14) THEN  ! mu nu_mu
                  ssbr(2) = bssmod(idc)
                  GOTO 20
                ELSE IF (jssmod(jmod,idc) .EQ. -16) THEN  ! tau nu_tau
                  ssbr(3) = bssmod(idc)
                  GOTO 20
                ELSE IF (jssmod(jmod,idc) .EQ. -2) THEN   ! u dbar
                  ssbr(4) = bssmod(idc)
                  GOTO 20
                ELSE IF (jssmod(jmod,idc) .EQ.  4) THEN   ! c sbar
                  ssbr(5) = bssmod(idc)
                  GOTO 20
                END IF
              END DO
   20         CONTINUE
            END IF
          END IF
        END DO
*
   10   CONTINUE
      ELSE IF (iproc .EQ. 6) THEN
*         Fix parameters of process gg -> V1V2
        IF (kv1 .GT. kv2) THEN
          kvtmp = kv1
          kv1   = kv2
          kv2   = kvtmp
        END IF
        IF      (kv1.EQ.113 .AND. kv2.EQ.113) THEN
          idx_vv = 1
        ELSE IF (kv1.EQ.113 .AND. kv2.EQ.223) THEN
          idx_vv = 2
        ELSE IF (kv1.EQ.113 .AND. kv2.EQ.333) THEN
          idx_vv = 3
        ELSE IF (kv1.EQ.113 .AND. kv2.EQ.443) THEN
          idx_vv = 4
        ELSE IF (kv1.EQ.223 .AND. kv2.EQ.223) THEN
          idx_vv = 5
        ELSE IF (kv1.EQ.223 .AND. kv2.EQ.333) THEN
          idx_vv = 6
        ELSE IF (kv1.EQ.223 .AND. kv2.EQ.443) THEN
          idx_vv = 7
        ELSE IF (kv1.EQ.333 .AND. kv2.EQ.333) THEN
          idx_vv = 8
        ELSE IF (kv1.EQ.333 .AND. kv2.EQ.443) THEN
          idx_vv = 9
        ELSE IF (kv1.EQ.443 .AND. kv2.EQ.443) THEN
          idx_vv = 10
        ELSE
          WRITE (6,*) '  %GGINIT: unknown (V1 V2) combination'
          WRITE (6,'(3x,''KV1 = '',i3,'', KV2 = '',i3)') kv1, kv2
          STOP
        END IF
        DO 101 igvpar=1,4
  101   gvpar(igvpar) = gvconst(igvpar,idx_vv)
        IF (gvpar(1) .EQ. 0.0) THEN
          WRITE (6,*) '  %GGINIT: (V1 V2) combination not implemented'
          WRITE (6,'(3x,''KV1 = '',i3,'', KV2 = '',i3)') kv1, kv2
          STOP
        END IF
      END IF
*
      IF (iproc.EQ.1 .OR. iproc.EQ.3 .OR.
     &    iproc.EQ.4 .OR. iproc.EQ.5 .OR. iproc.EQ.6) THEN
*         Max cross section estimation
        xsmax0 = 0.
        DO nest = 1,100
          xm  = amin + ggrnd(0) * (amax-amin)
          wsq = xm*xm
          CALL ggxsec
        END DO
        xssum = 0.
        ntry  = 0
      END IF
*
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GGWID.  
      FUNCTION ggwid (kf_res)
*-----------------------------------------------------------------------
*       Routine to define 2-gamma width of a resonance
*       with JETSET code KF_RES
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      REAL wid(5)
      DATA wid /9.E-09, 1.E-06, 5.E-06, 6.3E-06, 0.41E-06/
      IF      (kf_res .EQ. 111) THEN
        ggwid = wid(1)
      ELSE IF (kf_res .EQ. 221) THEN
        ggwid = wid(2)
      ELSE IF (kf_res .EQ. 331) THEN
        ggwid = wid(3)
      ELSE IF (kf_res .EQ. 441
     &    .OR. kf_res .EQ. 10441
     &    .OR. kf_res .EQ. 445) THEN
        ggwid = wid(4)
      ELSE IF (kf_res .EQ. 551
     &    .OR. kf_res .EQ. 10551
     &    .OR. kf_res .EQ. 555) THEN
        ggwid = wid(5)
      END IF
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GGRUN.  
      SUBROUTINE ggrun
*-----------------------------------------------------------------------
*       Single event generation
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
     &               amin,amax,ymin,ymax,nmas,ny, kferm,
     &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
     &               thetamin, costhv1, kv1,kv2,gvpar(4)
      CHARACTER lumfil*80
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      COMMON /ggevnt/ nrun,ievent,wsq,ygg,xmg1,xmg2, p2g(5),
     &                ptag1(4),ptag2(4), ngg, kgg(10),pgg(20,5)
      COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
     &  xstote, ssbr(10)
      COMMON /ggkin/ xkin(10)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      DOUBLE PRECISION P,V
      SAVE  /PYJETS/
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      DOUBLE PRECISION PARU,PARJ
      SAVE  /PYDAT1/
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      DOUBLE PRECISION PMAS,PARF,VCKM
      SAVE  /PYDAT2/
      REAL pgam1(4),pgam2(4)
      REAL p3(4), p4(4)
      REAL*8 dbetaz
      REAL*8 ggrnd
*
    1 CALL gg2gam (amas,wsq,ygg,q1sq,q2sq,
     &             pgam1,pgam2,ptag1,ptag2,weight)
*
*         Cross section calculation
      IF (iproc.EQ.1 .OR. iproc.EQ.3 .OR.
     &    iproc.EQ.4 .OR. iproc.EQ.5 .OR. iproc.EQ.6) THEN
        CALL ggxsec
        xssum = xssum + xscur
        ntry  = ntry  + 1
*         Select event on cross section of process gg -> XX
        IF (xscur0/xsmax0 .LT. sngl (ggrnd(0))) GOTO 1
*         Total cross section calculation due to current statistics
        xstot  = xssum/ntry
        xstote = xstot / sqrt(float(ntry))
      END IF
*
      xmgg   = sqrt (wsq)
      xmg1   = sqrt (q1sq)
      xmg2   = sqrt (q2sq)
      dbetaz = dtanh(1.d0*ygg)
*
*       calculation of the 4-vector of two-photon system
*       and filling arrays pgg(1:2,1:5), kgg(1:2)
      DO ip=1,4
        p2g(ip)   = pgam1(ip) + pgam2(ip)
        pgg(1,ip) = pgam1(ip)
        pgg(2,ip) = pgam2(ip)
      END DO
      p2g(5)   =  xmgg
      pgg(1,5) = -xmg1
      pgg(2,5) = -xmg2
      kgg(1)   = 22
      kgg(2)   = 22
      ngg      = 2
*
* =======> PROCESS # 1
      IF (iproc .EQ. 1) THEN
        DO ip = 1,4
          p(1,ip) = dble (pgam1(ip))
          p(2,ip) = dble (pgam2(ip))
        END DO
        p(1,5) = dble (-xmg1)
        p(2,5) = dble (-xmg2)
        CALL pyevnt
*
* =======> PROCESS # 2
      ELSE IF (iproc .EQ. 2) THEN
        ngg = ngg + 1
        DO ip = 1,5
          pgg(ngg,ip) = p2g(ip)
        END DO
        kgg(ngg) =  kf_onium
*         Fill the JETSET common block
        CALL gg2py (ngg)
        CALL pyexec
*
* =======> PROCESS # 3, # 4, # 5, # 6
      ELSE IF (iproc.EQ.3 .OR. iproc.EQ.4 .OR.
     &         iproc.EQ.5 .OR. iproc.EQ.6) THEN
        phi    = 2*pi*ggrnd(0)
        IF (iproc.EQ.3 .OR. iproc.EQ.4 .OR. iproc.EQ.5) THEN
          e3     = 0.5*sqrt(wsq)
          e4     = e3
          beta   = xkin(1)
          costhe = xkin(2)
          am1    = xkin(3)
          am2    = xkin(3)
          p3abs  = e3 * beta
        ELSE IF (iproc.EQ.6) THEN
          am1    = xkin(1)
          am2    = xkin(2)
          costhe = xkin(3)
          p3abs  = xkin(4)
          e3     = xkin(5)
          e4     = xmgg - e3
        END IF
        sinthe = sqrt ((1.-costhe)*(1.+costhe))
        p3(1) = p3abs * cos(phi) * sinthe
        p3(2) = p3abs * sin(phi) * sinthe
        p3(3) = p3abs * costhe
        p3(4) = e3
        p4(4) = e4
        DO ip = 1,3
          p4(ip) = -p3(ip)
        END DO
        igg1 = ngg + 1
        igg2 = ngg + 2
        DO ip = 1,4
          pgg(igg1,ip) = p3(ip)
          pgg(igg2,ip) = p4(ip)
        END DO
        pgg(igg1,5) = am1
        pgg(igg2,5) = am2
        ngg = ngg + 2
*
        IF      (iproc .EQ. 3) THEN
          kgg(igg1) =  kferm
          kgg(igg2) = -kferm
          CALL gg2py (ngg)
        ELSE IF (iproc .EQ. 4) THEN
          kgg(igg1) =  24
          kgg(igg2) = -24
          CALL gg2py (ngg)
        ELSE IF (iproc .EQ. 5) THEN
          kgg(igg1) =  41
          kgg(igg2) = -41
*           Decay of charginos into W+neutralino
          CALL ggdecy (igg1)
          CALL ggdecy (igg2)
*
          CALL gg2py (ngg)
          k(igg1,1) = 11
          k(igg2,1) = 11
        ELSE IF (iproc.EQ.6) THEN
          kgg(igg1) = kv1
          kgg(igg2) = kv2
          CALL gg2py (ngg)
        END IF
*
*         Boost event into initial frame
        mstu(33) = 1
        CALL pyrobo(3,ngg, 0.d0,0.d0, 0.d0,0.d0, dbetaz)
        CALL pyexec
*
      END IF
*
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GGEXIT. 
      SUBROUTINE ggexit
*-----------------------------------------------------------------------
*       Cross section calculation and end of run and
*       print out cross section and related variables
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
     &               amin,amax,ymin,ymax,nmas,ny, kferm,
     &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
     &               thetamin, costhv1, kv1,kv2,gvpar(4)
      CHARACTER lumfil*80
      COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
     &  xstote, ssbr(10)
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
*
      WRITE (6,'(/1x,79(''=''))')
      WRITE (6,'(1x,''I'',36x,''TPHIC'',36x,      ''I'')')
      WRITE (6,'(1x,''I'',34x,''Process '',i1,34x,''I'')') iproc
      WRITE (6,'(1x,''I'',34x,''---------'',  34x,''I'')')
      WRITE (6,'(1x,''I'',77x,                    ''I'')')
      WRITE (6,'(1x,''I'',20x,''Events thrown   : '',i8,31x,''I'')')
     &       ntry
      WRITE (6,'(1x,''I'',20x,''Events accepted : '',i8,31x,''I'')')
     &       nevent
      WRITE (6,'(1x,''I'',20x,''Cross section   : '',e10.3,'' +- '','//
     &      'e10.3,'' nb'',12x ''I'')') xstot, xstote
      WRITE (6,'(1x,''I'',77x,                    ''I'')')
      WRITE (6,'(1x,79(''=''))')
*
      IF (iproc .EQ. 1) THEN
*   ---  Alternative cross section for minimum bias events ---
        totxsc = 0.0
        hmas   = (amax-amin)/nmas
        f1     = dldmx(amin   ,ny) * dchad(amin)
        DO i=1,nmas
          xm  = amin + i*hmas
          f2  = dldmx(xm-hmas/2.,ny) * dchad(xm-hmas/2.)
          f3  = dldmx(xm        ,ny) * dchad(xm)
          totxsc = totxsc +(f1 + 4.*f2 + f3)/6.
          f1 = f3
        END DO
        sigtt = totxsc * hmas
        sigmx = sigtt/sqrt(float(nevent))
        WRITE (6,*) ' CrosSec =',sigtt,' +/-',sigmx,'  mb'
      END IF
*
      RETURN
      END
*
      REAL FUNCTION dchad(xm)
*-----------------------------------------------------------------------
*       Total cross section for gamma-gamma collisions vs. Xm in mb
*       Parametrization is given by G.A.Schuler, T.Sjostrand,
*       CERN-TH.7193/94, formula (9).
*       Coded by Yu.Kharlov 20.03.1995.
*-----------------------------------------------------------------------
      DATA epsilon, eta /0.0808, 0.4525/
      dchad = 0.211e-3 * (xm**2)**epsilon + 0.297e-3/(xm**2)**eta
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GGXSEC. 
      SUBROUTINE ggxsec
*-----------------------------------------------------------------------
*       Cross section calculation in current event for iproc=1,3,4,5,6
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
     &               amin,amax,ymin,ymax,nmas,ny, kferm,
     &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
     &               thetamin, costhv1, kv1,kv2,gvpar(4)
      CHARACTER lumfil*80
      COMMON /ggevnt/ nrun,ievent,wsq,ygg,xmg1,xmg2, p2g(5),
     &                ptag1(4),ptag2(4), ngg, kgg(10),pgg(20,5)
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
     &  xstote, ssbr(10)
      COMMON /ggkin/ xkin(10)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      DOUBLE PRECISION PMAS,PARF,VCKM
      SAVE  /PYDAT2/
      PARAMETER (epsilon = 0.0808, eta = -0.4525)
      REAL*8 ggrnd
      REAL*8 pymass
*
      IF (iproc .EQ. 1) THEN
*         Parametrization by G.A.Schuler, T.Sjostrand,
*         CERN-TH.7193/94, formula (9).
        xscur0 = 211. * wsq**epsilon + 297. * wsq**eta
        xscur  = xscur0 * xlumint
      ELSE IF (iproc .EQ. 3 .OR. iproc .EQ. 5) THEN
        beta2  = 1.0 - 4.0*amprt(1)**2/wsq
        beta2  = min  (beta2, 1.0)
        beta   = sqrt (beta2)
        IF (beta2 .LT. 0.9) THEN
    1     costhe = 2.*ggrnd(0) - 1.
          costh2 = costhe*costhe
          xnum   = 1. + 2.*beta2 * (1.-beta2) * (1.-costh2) -
     &             (beta2*costh2)**2
          xden   = (1 - beta2*costh2)**2
          ffcur  = xnum / xden
          ffmax  = (1. + beta2) / (1. - beta2)
          rff    = ggrnd(0)
          IF (ffcur .LE. rff*ffmax) GOTO 1
        ELSE
          small  = 1.0 - cos(thetamin)
    2     zz     = exp (ggrnd(0) * alog(2./small))
C          sgn    = sign (1., 2.*ggrnd(0)-1.)
          sgn = 2.*ggrnd(0)-1.
          sgn    = sign (1., sgn)
          costhe = sgn * (zz-1.) / (zz+1.)
          ffcur  = (1.+costhe*costhe) / (1.-beta2*costhe*costhe)
          ffmax  =                 2. / (1.-      costhe*costhe)
          rff    = ggrnd(0)
          IF (ffcur .LE. rff*ffmax) GOTO 2
        END IF
*
        xkin(1) = beta
        xkin(2) = costhe
        xkin(3) = amprt(1)
*
*           Yu.Kharlov, 3.06.1995
        totfac = 4. * pi * alpha**2 * gev2nb * qf**4 * nc / wsq
        IF (beta2 .LT. 0.99) THEN
          xscur0 =((3.-beta2**2)/(2.*beta) * alog ((1.+beta)/(1.-beta))-
     &             2. + beta2) * beta
        ELSE
          xscur0 = 2. * alog (sqrt(wsq)/amprt(1)) - 1.0
        END IF
        xscur0 = totfac * xscur0 * xsbra
        xscur = xscur0 * xlumint
      ELSE IF (iproc .EQ. 4) THEN
        xmw = amprt(4)
        beta2 = 1. - 4.*xmw**2/wsq
        beta = sqrt (beta2)
    3   costhe = 2.*ggrnd(0) - 1.
        costh2 = costhe*costhe
        xnum = 19. - 6.*beta2 * (1.-beta2) +
     &         2.*beta2 * (8. - 3.*beta2) * costh2 +
     &         3.*(beta2*costh2)**2
        xden = (1. - beta2*costh2)**2
        ffcur = xnum / xden
        ffmax = (19. + 10.*beta2 + beta2**2) / (1. - beta2)**2
        rff = ggrnd(0)
        IF (ffcur .LE. rff*ffmax) GOTO 3
*
        xkin(1) = beta
        xkin(2) = costhe
        xkin(3) = xmw
*
*           Yu.Kharlov, 3.06.1995
        totfac = pi * alpha**2 * beta  * gev2nb/ wsq
        xscur0 = 2. * (22. - 9.*beta2 + 3.*beta2**2) / (1. - beta2) -
     &           3. * (1.-beta2**2)/beta * alog ((1.+beta)/(1.-beta))
        xscur0 = totfac * xscur0 * xsbra
        xscur  = xscur0 * xlumint
      ELSE IF (iproc .EQ. 6) THEN
        xmgg = sqrt (wsq)
        bb   = gvpar(3) + gvpar(4)*alog(xmgg/50.0)
        am1  = sngl (pymass(kv1))
        am2  = sngl (pymass(kv2))
        p1   = pfork (xmgg,am1,am2)
        e1   = sqrt (p1*p1 + am1*am1)
        tmin = am1*am1 - xmgg*(e1+p1*costhv1)
        tmax = am1*am1 - xmgg*(e1-p1*costhv1)
        fmin = exp (bb*tmin)
        fmax = exp (bb*tmax)
        IF (fmax .GT. 1.e-10) THEN
          tvar = 1./bb * dlog (fmax - ggrnd(0)*(fmax-fmin))
        ELSE
          tvar = tmax + 1./bb * dlog (1. - ggrnd(0))
        END IF
        costhe = (e1*xmgg - am1*am1 + tvar) / xmgg / p1
        IF (costhe .GT. 1.0) costhe =  1.0
        IF (costhe .LT.-1.0) costhe = -1.0
        IF (ggrnd(0) .LT. 0.5) costhe = -costhe
*
        xkin(1) = am1
        xkin(2) = am2
        xkin(3) = costhe
        xkin(4) = p1
        xkin(5) = e1
*
        IF (kv1 .EQ. 443 .AND. kv2 .EQ. 443) THEN
          factor = 1.0 / alog(xmgg/3.097)
        ELSE
          factor = 1.0
        END IF
        xscur0 = factor * gvpar(1) * xmgg**gvpar(2) / bb * exp(bb*tmax)
        xscur  = xscur0 * xlumint
*
      END IF
*
*       Find max cross section of gg-process
*
      IF (xscur0 .GT. xsmax0) xsmax0 = xscur0
*
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GGDECY. 
      SUBROUTINE ggdecy (igg)
*-----------------------------------------------------------------------
*       Routine to decay a particle number I in GGEVNT common block
*       Input  : IGG - particle number in common GGEVNT arrays
*       Author : Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
     &               amin,amax,ymin,ymax,nmas,ny, kferm,
     &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
     &               thetamin, costhv1, kv1,kv2,gvpar(4)
      CHARACTER lumfil*80
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      COMMON /ggevnt/ nrun,ievent,wsq,ygg,xmg1,xmg2, p2g(5),
     &                ptag1(4),ptag2(4), ngg, kgg(10),pgg(20,5)
      COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
     &  xstote, ssbr(10)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      DOUBLE PRECISION PMAS,PARF,VCKM
      SAVE  /PYDAT2/
      COMMON /ggmssm/ xm1,   xm2,    xmg,    xms,    xmtl,   xmtr,
     &                xmll,  xmlr,   xmnl,   xtanb,  xmha,   xmu,
     &                xmt,   xat,    xmbr,   xab,    u11,    v11
      REAL p1(4), pcm(4), p2(4), p3(4), p4(4), ptmp(3)
      REAL o1(3,3), o2(3,3)
      INTEGER pycomp
      DATA ifl /0/
      REAL*8 ggrnd
*
      DO ip = 1,4
        p1(ip) = pgg(igg,ip)
      END DO
*
      IF (abs(kgg(igg)) .EQ. 41 .AND. iproc .EQ. 3) THEN
*         Decay chi_1^+- --> chi_1^0 W^+- in CMS, uniformly
        eb = (amprt(1)**2 + amprt(2)**2 - amprt(4)**2) / (2 * amprt(1))
        ec = (amprt(1)**2 - amprt(2)**2 + amprt(4)**2) / (2 * amprt(1))
        pb = sqrt ((eb-amprt(2))*(eb+amprt(2)))
        pc = pb
        costhe = 2*ggrnd(0) - 1.
        sinthe = sqrt((1.-costhe)*(1.+costhe))
        phi = 2*pi*ggrnd(0)
        pcm(1) = pb * sinthe * cos(phi)
        pcm(2) = pb * sinthe * sin(phi)
        pcm(3) = pb * costhe
        pcm(4) = eb
        CALL lorenb (amprt(1),p1,pcm,p2)
        DO ip = 1,3
          pcm(ip) = -pcm(ip)
        END DO
        pcm(ip) = ec
        CALL lorenb (amprt(1),p1,pcm,p3)
        DO ip = 1,4
          pgg(ngg+1,ip) = p2(ip)
          pgg(ngg+2,ip) = p3(ip)
        END DO
        pgg(ngg+1,5) = amprt(2)
        pgg(ngg+2,5) = amprt(4)
        kgg(ngg+1) = 42 * sign(1.,1.*kgg(igg))
        kgg(ngg+2) = 24 * sign(1.,1.*kgg(igg))
        ngg = ngg + 2
*
      ELSE IF (abs(kgg(igg)) .EQ. 41 .AND. iproc .EQ. 5) THEN
*         Decay chi_1^+- --> e^+- nu_e chi_1^0 in CMS,
*         according to the exact matrix element
*
        xm23max = amprt(1) - amprt(2)
        xmw     = sngl (pmas(pycomp(24),1))
        gwt     = sngl (pmas(pycomp(24),2))
        atmin   = - atan (xmw/gwt)
        atmax   =   atan ((xm23max**2 - xmw**2) / (xmw*gwt))
        am1     = amprt(1)
        am4     = amprt(2)
*
    1   CONTINUE
*
*       Random mass of (2+3) according to Breit-Wigner
*
        xm23 = xmw**2 + xmw*gwt *
     *       tan (atmin + ggrnd(0)*(atmax - atmin))
        xm23 = max (0., xm23)
        xm23 = sqrt (xm23)
        IF (xm23 .GT. xm23max) GOTO 1       ! to avoid boundary effect
        e23  = (am1**2 + xm23**2 - am4**2) / 2. / am1
*
*       Rest frame of (2+3)
*
        phi    = 2.*pi * ggrnd(0)
        costhe = 2.*ggrnd(0) - 1.
        sinthe = sqrt ((1.-costhe)*(1.+costhe))
*
        pp = xm23 / 2.
        p2(1) = pp * sinthe * cos(phi)
        p2(2) = pp * sinthe * sin(phi)
        p2(3) = pp * costhe
        p2(4) = pp
*
        DO ip = 1,3
          p3(ip) = -p2(ip)
        END DO
        p3(4) = pp
*
*       Transformation from rest frame of (2+3) to rest frame of (1),
*       3-momenta p2, p3 are along z-axis
*
        gamma = e23 / xm23
        gambe = sqrt ((gamma-1.)*(gamma+1.))
*
        p2z = gamma*p2(3) + gambe*p2(4)
        e2  = gamma*p2(4) + gambe*p2(3)
        p3z = gamma*p3(3) + gambe*p3(4)
        e3  = gamma*p3(4) + gambe*p3(3)
*
        p2(3) = p2z
        p2(4) = e2
        p3(3) = p3z
        p3(4) = e3
*
*       Arbitrary rotation of rest frame of (1)
*
        phi    = 2.*pi * ggrnd(0)
        costhe = 2.*ggrnd(0) - 1.
        sinthe = sqrt ((1.-costhe)*(1.+costhe))
        cosphi = cos(phi)
        sinphi = sin(phi)
*
        o1(1,1) =  costhe
        o1(1,2) =  0.
        o1(1,3) = -sinthe
        o1(2,1) =  0.
        o1(2,2) =  1.
        o1(2,3) =  0.
        o1(3,1) =  sinthe
        o1(3,2) =  0.
        o1(3,3) =  costhe
*
        o2(1,1) =  cosphi
        o2(1,2) = -sinphi
        o2(1,3) =  0.
        o2(2,1) =  sinphi
        o2(2,2) =  cosphi
        o2(2,3) =  0.
        o2(3,1) =  0.
        o2(3,2) =  0.
        o2(3,3) =  1.
*
        DO ix = 1,3
          ptmp(ix) = 0.
          DO jx = 1,3
            ptmp(ix) = ptmp(ix) + o1(ix,jx) * p2(jx)
          END DO
        END DO
        DO ix = 1,3
          p2(ix) = ptmp(ix)
        END DO
        DO ix = 1,3
          ptmp(ix) = 0.
          DO jx = 1,3
            ptmp(ix) = ptmp(ix) + o2(ix,jx) * p2(jx)
          END DO
        END DO
        DO ix = 1,3
          p2(ix) = ptmp(ix)
        END DO
*
        DO ix = 1,3
          ptmp(ix) = 0.
          DO jx = 1,3
            ptmp(ix) = ptmp(ix) + o1(ix,jx) * p3(jx)
          END DO
        END DO
        DO ix = 1,3
          p3(ix) = ptmp(ix)
        END DO
        DO ix = 1,3
          ptmp(ix) = 0.
          DO jx = 1,3
            ptmp(ix) = ptmp(ix) + o2(ix,jx) * p3(jx)
          END DO
        END DO
        DO ix = 1,3
          p3(ix) = ptmp(ix)
        END DO
*
        p12 = am1 * p2(4)
*
        e4 = am1 - e23
        pp =sqrt ((e4-am4)*(e4+am4))
        p4(1) =  pp * sinthe * cos(phi)
        p4(2) =  pp * sinthe * sin(phi)
        p4(3) = -pp * costhe
        p4(4) =  e4
*
        width = (- 4.*p12**2*u11**2      - 4.*p12**2*v11**2 +
     &           4.*p12*xm23**2*u11**2 - 2.*p12*am4**2*u11**2 -
     &           2.*p12*am4**2*v11**2  + 2.*p12*am1**2*u11**2 +
     &           2.*p12*am1**2*v11**2  - xm23**4*u11**2 +
     &           xm23**2*am4**2*u11**2 - 2.*xm23**2*am4*am1*u11*v11 -
     &           xm23**2*am1**2*u11**2) /
     &           ((xm23**2-xmw**2)**2 + (xmw*gwt)**2)
        width = width * pp * xm23
*
        widmax = (2.*am1*xm23max *
     &    (2.*(xm23max*u11)**2 + am1**2*(u11**2 + v11**2)) +
     &    (xm23max*am4*u11)**2) * (am1**2 - am4**2) / 2./am1 *
     &    xm23 / ((xm23**2-xmw**2)**2 + (xmw*gwt)**2)
*
*       Select current phase space configuration
*
        IF (width .GT. widmax) WRITE (6,*)
     &    '%GGDECY: Maximum chi+ width violated'
        rwid = widmax * ggrnd(0)
        IF (width .LE. rwid) GOTO 1
*
*         Pick up final lepton flavour
        rflav = ggrnd(0)
        IF (moddcy .EQ. 1) THEN
          IF (rflav .LE. 0.5) THEN
            lflav1 =  11
            lflav2 = -12
          ELSE
            lflav1 =  13
            lflav2 = -14
          END IF
        ELSE
          rr1 = ssbr(1)
          rr2 = ssbr(1)+ssbr(2)
          rr3 = ssbr(1)+ssbr(2)+ssbr(3)
          rr4 = ssbr(1)+ssbr(2)+ssbr(3)+ssbr(4)
          rr5 = ssbr(1)+ssbr(2)+ssbr(3)+ssbr(4)+ssbr(5)
          IF (rflav .LE. rr1) THEN
            lflav1 =  11
            lflav2 = -12
          ELSE IF (rflav .LE. rr2) THEN
            lflav1 =  13
            lflav2 = -14
          ELSE IF (rflav .LE. rr3) THEN
            lflav1 =  15
            lflav2 = -16
          ELSE IF (rflav .LE. rr4) THEN
            lflav1 = -2
            lflav2 =  1
          ELSE IF (rflav .LE. rr5) THEN
            lflav1 = -4
            lflav2 =  3
          ELSE
            GOTO 1
          END IF
        END IF
*
*         Transformation to CMS of 2-gamma
        CALL lorenb (amprt(1),p1,p2,pcm)
        DO ip = 1,4
          pgg(ngg+1,ip) = pcm(ip)
        END DO
        pgg(ngg+1,5) = 0.
        kgg(ngg+1) = lflav1 * sign(1.,1.*kgg(igg))
        ngg = ngg + 1
*
        CALL lorenb (amprt(1),p1,p3,pcm)
        DO ip = 1,4
          pgg(ngg+1,ip) = pcm(ip)
        END DO
        pgg(ngg+1,5) = 0.
        kgg(ngg+1) = lflav2 * sign(1.,1.*kgg(igg))
        ngg = ngg + 1
*
        CALL lorenb (amprt(1),p1,p4,pcm)
        DO ip = 1,4
          pgg(ngg+1,ip) = pcm(ip)
        END DO
        pgg(ngg+1,5) = amprt(2)
        kgg(ngg+1) = 42 * sign(1.,1.*kgg(igg))
        ngg = ngg + 1
      END IF
*
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GG2PY.  
      SUBROUTINE gg2py (nngg)
*-----------------------------------------------------------------------
*       Routine to fill PYTHIA event record by particles from GG
*       Input: NNGG : number of particles
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggevnt/ nrun,ievent,wsq,ygg,xmg1,xmg2, p2g(5),
     &                ptag1(4),ptag2(4), ngg, kgg(10),pgg(20,5)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      DOUBLE PRECISION P,V
      SAVE  /PYJETS/
      DO igg = 1,nngg
        IF (igg .LE. 2) THEN
          k(igg,1) = 21
        ELSE
          k(igg,1)=1
        END IF
        k(igg,2) = kgg(igg)
        DO ip = 1,5
          p(igg,ip) = pgg(igg,ip)
        END DO
      END DO
      n = nngg
      RETURN
      END
*
*=======================================================================
CDECK  ID>, PFORK.  
      FUNCTION pfork (am0,am1,am2)
*-----------------------------------------------------------------------
*       Routine to find a momentum of a secondary particle in the decay
*       of a particle with mass AM0 in the rest to 2 particles with
*       masses AM1 and AM2.
*-----------------------------------------------------------------------
      denom = am0**4 + am1**4 + am2**4 -
     &        2.*(am0*am1)**2 - 2.*(am0*am2)**2 - 2.*(am1*am2)**2
      IF (denom .GT. 0.) THEN
        pfork = sqrt (denom) / 2./am0
      ELSE
        pfork = -1.
      END IF
      RETURN
      END
*=======================================================================
CDECK  ID>, GG2GAM. 
      SUBROUTINE gg2gam (amas,wsq,y3,qs1,qs2,pgam1,pgam2,ptag1,ptag2,wt)
C***********************************************************************
C
C     Generator of the Direct Two Photon Processes in      S.A.Sadovsky
C     coherent heavy ion collisions:    A+B --> 1+2+3       17.01.1995
C     with double differencial gg-luminosity function      Yu.Kharlov
C                                                           20.03.1995
C
C      LUMFUN must be called first to initialise the generator
C
C***********************************************************************
C
C   DETAILED DESCRIPTION
C   ======== ===========
C
C   INPUTS:  only through common /TWOGENC/
C   =======
C
C   OUTPUTS:     AMAS   The mass of ion.
C   ========     WSQ    Square mass of the gamma-gamma system
C                Y3     The rapidity of the gamma-gamma system
C                QS1    Virtual mass squared of photon 1.
C                QS2    Virtual mass squared of photon 2.
C                PGAM1  Four-vector of photon 1
C                PGAM2  Four-vector of photon 2
C                PTAG1  Four-vector of tag 1
C                PTAG2  Four-vector of tag 2
C                WT     The weight of the event (dummy so far)
C
C   The remaining variables are internal to the program.
C
C***********************************************************************
      INTEGER   icnter,ncntrs,nevts,lout,debug
      PARAMETER (ncntrs =7)
      COMMON /twgenc/s,eb,wmin,wmax,th1pmn,th1pmx,th2pmn,th2pmx,
     +   st12mn,st12mx,st12rt,st22mn,st22mx,st22rt,ommin,omrat,eps,
     +   wghtsm,wghtgt,wghtmx,wtmxfd,xmres,xgres,
     +   icnter(ncntrs),nevts,lout,lres,lwate,debug
      SAVE   /twgenc/
      REAL      st12mn,st12mx,st12rt,st22mn,st22mx,st22rt
      REAL      wmin,wmax,ommin,omrat,eps,th1pmn,th1pmx
      REAL      th2pmn,th2pmx,s,eb,wghtsm,wghtmx,wghtgt,wtmxfd
      REAL      xmres,xgres
      LOGICAL   lres,lwate
C
      REAL       alpha,pi
      PARAMETER (alpha=1./137.036,pi=3.14159265)
C
      REAL ptag1(10),ptag2(10),pgam1(10),pgam2(10),wt
      REAL amas,wsq
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL*4  y3,xm,q0,qs,qs1,qs2
      REAL*8  pa(5),pb(5)
C
      REAL*8  pi2
C
      LOGICAL     lstr
      DATA q0,pi2,lstr /  0.060, 6.28318530718 d0, .false. /
      REAL*8 ggrnd
C
  999 CONTINUE
      icnter(1) = icnter(1) + 1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF    (lstr) go to 10
      lstr =.true.
      qs   = q0/sqrt(2.)
C
      pa(1)= 0D0
      pa(2)= 0D0
      pa(3)= dsqrt(1D0*((eb-amas)*(eb+amas)))
      pa(4)= eb*1D0
      pa(5)= amas*1D0
C
      pb(1)= 0D0
      pb(2)= 0D0
      pb(3)=-pa(3)
      pb(4)= pa(4)
      pb(5)= amas*1D0
C
   10 xm   = xmres
*       Generate 2gamma rapidity y3 and mass xm
      CALL ranlum(y3,xm,wt)
C
      wsq  = xm*xm
C
*       Energy and z-momentum of 2gamma system
      egg = xm*cosh(y3)
      pgg = xm*sinh(y3)
*
*       Squared masses of the photons
      qs1 = - qs**2 * dlog(ggrnd(0))    ! is it true that squared mass
      qs2 = - qs**2 * dlog(ggrnd(0))    ! is distributed in exponent?
      xm1 = -qs1
      xm2 = -qs2
*
*       Energies and z-momenta of the photons
      xfun12 = wsq + xm1 - xm2
*
      eg1 = (egg + pgg*sqrt(1 - 4*wsq*xm1/(xfun12*xfun12))) *
     &      xfun12/(2*wsq)
      eg2 = egg - eg1
      pg1 =  sqrt (eg1**2 - xm1)
      pg2 = -sqrt (eg2**2 - xm2)
*
*       4-momenta of the photons
      DO ip=1,2
        pgam1(ip) = 0.
        pgam2(ip) = 0.
      END DO
      pgam1(3) = pg1
      pgam1(4) = eg1
      pgam2(3) = pg2
      pgam2(4) = eg2
*
*       4-momenta of the recoil ions
      DO ip=1,4
        ptag1(ip) = pa(ip) - pgam1(ip)
        ptag2(ip) = pb(ip) - pgam2(ip)
      END DO
*
      RETURN
      END
*
*=======================================================================
CDECK  ID>, RANLUM. 
      SUBROUTINE ranlum (ry,rm,wlum)
C                                                     G.V.Khaoustov
C                                                       12.12.1994
C                                             Corrected: Yu.Kharlov
C                                                       23.05.1995
C                                                       24-JUN-1997
C
C generate two random numbers (ry and rm) with probability proportional
C to luminosity function in ymin-ymax rapidity and xmin-xmax mass range
C
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      PARAMETER (isize=100000)
      COMMON/lumfunct/anorm,cexp,ymin,ymax,ymn,ymx,nny,
     +                           xmin,xmax,amn,amx,nnms,alfun(isize)
      PARAMETER (NCNTRS =7)
      COMMON /TWGENC/S,EB,WMIN,WMAX,TH1PMN,TH1PMX,TH2PMN,TH2PMX,
     +   ST12MN,ST12MX,ST12RT,ST22MN,ST22MX,ST22RT,OMMIN,OMRAT,EPS,
     +   WGHTSM,WGHTGT,WGHTMX,WTMXFD,XMRES,XGRES,
     +   ICNTER(NCNTRS),NEVTS,LOUT,LRES,LWATE,DEBUG
      LOGICAL lres,lwate
      REAL*4  wlum
      DATA    ifl/0/, stpy/0./, stpm/0./
      REAL*8 ggrnd
C
      IF (ifl.EQ.0) THEN
        IF (ymin.LT.ymn .OR. ymin.GT.ymx .OR.
     *      ymax.LT.ymn .OR. ymax.GT.ymx .OR.
     *      xmin.LT.amn .OR. xmin.GT.amx .OR.
     *      xmax.LT.amn .OR. xmax.GT.amx) THEN
          WRITE (6,*) 'Luminosity function is not defined in this range'
          STOP
        ENDIF
        aky = ymax-ymin
        akm = xmax-xmin
        stpy = (ymx-ymn)/nny
        stpm = (amx-amn)/nnms
        rmin = exp(cexp*(xmax))
        rmax = exp(cexp*(xmin))
        dr = rmax-rmin
        ifl = 1
      ENDIF
C
    1 CONTINUE
      IF (.NOT.lres) THEN
        rm  = ggrnd(0)                     ! random mass
        rm  = rmin + rm*dr
        rm  = alog(rm)/cexp                ! exponetial law
      ELSE
        wsq = xmres**2 + xmres*xgres * tan(pi*(ggrnd(0)-0.5))
        IF (wsq .LT. xmin*xmin) GOTO 1
        IF (wsq .GT. xmax*xmax) GOTO 1
        rm  = sqrt (wsq)
      END IF
C
      ry = ggrnd(0)                        ! random Y
*
*       Find y-limits for generated gg-mass: y_max = log(E_max/M_gg)
*
      ym_max  = alog (egg_max/rm)
      ygg_max = min (ymax, ym_max)
      ygg_min = max (ymin,-ym_max)
      ry = ygg_min + ry*(ygg_max-ygg_min)
      iy = (ry-ymn)/stpy+1.
C
      im = (rm-amn)/stpm
      r = anorm*exp(cexp*rm)*ggrnd(0)      ! normalization
      p = (ry-(iy-1)*stpy-ymn)/stpy        ! four point interpolation
      q = (rm-im*stpm-amn)/stpm
      k = iy+im*(nny+1)
      f = (1.-p)*(1.-q)*alfun(k)       + p*(1.-q)*alfun(k+1) +
     *         q*(1.-p)*alfun(k+nny+1) +      q*p*alfun(k+2+nny)
      IF (r .GT. f) GOTO 1
      wlum =  1000.*f                      ! differential luminosity in GeV^-1
*
      RETURN
      END
*
*=======================================================================
CDECK  ID>, GGDATA. 
      BLOCK DATA ggdata
*-----------------------------------------------------------------------
*       Some useful constants and masses of MSSM particles
*       Author: Yu.Kharlov
*-----------------------------------------------------------------------
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      COMMON /ggini/ iproc,nevent,ilumf,lumfil,ebmn,eb,iz,ia,amas,
     &               amin,amax,ymin,ymax,nmas,ny, kferm,
     &               kf_onium,xmres,xgtres,xggres, xlumint, moddcy,
     &               thetamin, costhv1, kv1,kv2,gvpar(4)
      CHARACTER lumfil*80
      COMMON /ggxs/ xsmax0, xscur0, xscur, xsbra, xssum, ntry, xstot,
     &  xstote, ssbr(10)
      COMMON /ggmssm/ xm1,   xm2,    xmg,    xms,    xmtl,   xmtr,
     &                xmll,  xmlr,   xmnl,   xtanb,  xmha,   xmu,
     &                xmt,   xat,    xmbr,   xab,    u11,    v11
      COMMON/D2LParam/lz,la,Rion,Gamma
      INTEGER lz,la
      REAL    Rion,Gamma
      PARAMETER (al = 1./128.)
      DATA pi /3.14159276/, alpha /al/, gev2nb /389379./,
     &     hbarc /0.197327/
      DATA iproc, nevent, ilumf /1, 10000, -1/
      DATA lumfil /'gglum.dat'/
      DATA Rion /0.0/
      DATA amprt /110.,10.,150.,-1.,0.106/
      DATA thetamin /1.e-03/
      DATA costhv1 /1.0/
      DATA xsbra /1./
      DATA ny, nmas /50, 50/
      DATA xm1,   xm2,    xmg,    xms,    xmtl,   xmtr,
     &     xmll,  xmlr,   xmnl,   xtanb,  xmha,   xmu,
     &     xmt,   xat,    xmbr,   xab
     &     /30.,   50.,   250.,   700.,   600.,   600.,
     &     400.,  400.,   400.,   4.,     300.,  -300.,
     &     174.,  300.,   300.,   300./
      DATA gvconst / 63.0, 0.22, 6.0, 1.0,    ! rho rho
     &               14.0, 0.22, 6.0, 1.0,    ! rho omega
     &                7.0, 0.22, 4.0, 1.0,    ! rho phi
     &               0.05, 0.80, 2.5, 0.0,    ! rho J/psi
     &                0.0, 0.00, 0.0, 0.0,    ! omega omega
     &                0.8, 0.22, 4.0, 1.0,    ! omega phi
     &              6.E-3, 0.80, 2.5, 0.0,    ! omega J/psi
     &                0.2, 0.22, 2.0, 1.0,    ! phi phi
     &              3.E-3, 0.80, 1.5, 0.0,    ! phi J/psi
     &              6.E-5, 2.00, 1.5, 0.0/    ! J/psi J/psi
      END
*
*=======================================================================
CDECK  ID>, EVPRINT.
      SUBROUTINE evprint(Lun,Iev,P3,Ygg,ifl,Eb,ia,iz,amas,ptag1,ptag2)
*-----------------------------------------------------------------------
*       Event output for internal use by GGHIC authors
*-----------------------------------------------------------------------
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      DOUBLE PRECISION P,V
      SAVE  /PYJETS/
      REAL ptag1(4),ptag2(4),p3(5)
C
      WRITE(lun,100) iev,n,Eb,p3(5),ygg,ifl,-ifl
      WRITE(lun,101) 90000+ia,iz,ptag1,amas
      WRITE(lun,101) 90000+ia,iz,ptag2,amas
      WRITE(lun,101) 90000+ 3, 0,(p3(j),j=1,5)
      IF (n.NE.0) THEN
        DO i=1,n
          ichg=pychge(k(i,2))/3
          WRITE(lun,101) k(i,2),ichg,(p(i,j),j=1,5)
        ENDDO
      ENDIF
      RETURN
  100 FORMAT(i6,i3,3e16.8,2x,2i8)
  101 FORMAT(i6,i3,3(1x,e14.8),e13.8,e13.8)
      END
*
*=======================================================================
C
CDECK  ID>, TWOINI. 
      SUBROUTINE TWOINI(LOUTX,IPAR,XPAR,YPAR)                           TWOG0965
C***********************************************************************TWOG0966
C                                                                       TWOG0967
C     Initialises the gamma-gamma generator of Ion collisions TWOGAM.   TWOG0968
C                                                                       TWOG0969
C   INPUTS:             LOUTX    Logical unit for print-out             TWOG0970
C   ======= XPAR(1)     EBM      Beam energy                        GeV TWOG0971
C           XPAR(2)     WMN      Minimum two-gamma mass             GeV TWOG0972
C           XPAR(3)     WMX      Maximum two-gamma mass             GeV TWOG0973
C                                                                       TWOG0974
C           XPAR(8)     WGHTMX   Maximum weight                         TWOG0978
C           XPAR(9)     XMRES    Mass of resonance                  GeV TWOG0979
C           XPAR(10)    XGRES    Total width of resonance           GeV TWOG0980
C                                                                       TWOG0981
C           IPAR(1)     0        Unweighted events                      TWOG0982
C                       1        Weighted events                        TWOG0983
C           IPAR(2)     0        Continuum production                   TWOG0984
C                       1        Resonance formation                    TWOG0985
C                                                                       TWOG0986
C   EXTERNALS:   none                                                   TWOG0987
C   ==========                                                          TWOG0988
C                                                                       TWOG0989
C   AUTHOR/DATE:  S.A. Sadovsky,   17 January 1995                      TWOG0990
C   ===========                                                         TWOG0991
C                                                                       TWOG0992
C***********************************************************************TWOG0993
C     IMPLICIT NONE                                                     TWOG0994
      INTEGER   ICNTER,NCNTRS,NEVTS,LOUT,DEBUG                          TWOG0995
      PARAMETER (NCNTRS =7)                                             TWOG0996
      COMMON /TWGENC/S,EB,WMIN,WMAX,TH1PMN,TH1PMX,TH2PMN,TH2PMX,        TWOG0997
     +   ST12MN,ST12MX,ST12RT,ST22MN,ST22MX,ST22RT,OMMIN,OMRAT,EPS,     TWOG0998
     +   WGHTSM,WGHTGT,WGHTMX,WTMXFD,XMRES,XGRES,                       TWOG0999
     +   ICNTER(NCNTRS),NEVTS,LOUT,LRES,LWATE,DEBUG                     TWOG1000
      SAVE   /TWGENC/                                                   TWOG1001
      REAL      ST12MN,ST12MX,ST12RT,ST22MN,ST22MX,ST22RT               TWOG1002
      REAL      WMIN,WMAX,OMMIN,OMRAT,EPS,TH1PMN,TH1PMX                 TWOG1003
      REAL      TH2PMN,TH2PMX,S,EB,WGHTSM,WGHTMX,WGHTGT,WTMXFD          TWOG1004
      REAL      XMRES,XGRES                                             TWOG1005
      LOGICAL   LRES,LWATE                                              TWOG1006
C                                                                       TWOG1007
      REAL      ALPHA,PI,TWOPI,XME,XMU                                  TWOG1008
      PARAMETER (ALPHA=1./137.036,PI=3.14159265)                        TWOG1009
      PARAMETER (TWOPI=2.*PI,XME=0.000511,XMU=0.105658)                 TWOG1010
C                                                                       TWOG1011
      REAL     XPAR(11),YPAR(6)                                         TWOG1012
      INTEGER  LOUTX,IPAR(10),I                                         TWOG1013
C                                                                       TWOG1014
      EB    = XPAR(1)                                                   TWOG1015
      WMIN  = XPAR(2)                                                   TWOG1016
      WMAX  = XPAR(3)                                                   TWOG1017
      WGHTMX= XPAR(8)                                                   TWOG1022
      WTMXFD= 0                                                         TWOG1023
      XMRES = XPAR(9)                                                   TWOG1024
      XGRES = XPAR(10)                                                  TWOG1025
      LOUT  = LOUTX                                                     TWOG1026
      LWATE = IPAR(1) .EQ. 1                                            TWOG1027
      LRES  = IPAR(2) .EQ. 1                                            TWOG1028
      EPS   = 1./64.*(XME*WMIN**2/EB**3)**2/100.                        TWOG1029
C+                                                                      TWOG1030
C EPS is a small constant. We throw 1/2 COS(th/2)/(eps+SIN(th/2))       TWOG1031
C rather than 1/2 COT(th/2) to be able to start theta from zero.        TWOG1032
C-                                                                      TWOG1033
      S     = (2.*EB)**2                                                TWOG1054
      OMMIN = WMIN**2/(4.*EB)                                           TWOG1055
      OMRAT = EB/OMMIN                                                  TWOG1056
C+                                                                      TWOG1057
C Print out initial conditions.                                         TWOG1058
C-                                                                      TWOG1059
      WRITE(LOUT,1000)                                                  TWOG1060
 1000 FORMAT(1H ,10X,27('-'),' TWINIT ',27('-'))                        TWOG1061
      WRITE(LOUT,1001)                                                  TWOG1062
 1001 FORMAT(1H ,10X,'I',60X,'I')                                       TWOG1063
C
      WRITE(LOUT,1002)EB,YPAR(3),YPAR(4),YPAR(1),YPAR(2),WMIN,WMAX,
     +                                   YPAR(5),YPAR(6)
C                                                                       TWOG1065
 1002 FORMAT                                                            TWOG1066
     +  (1H ,10X,'I    Beam energy               ',E15.5,               TWOG1067
     +  ' GeV',11X,'I',/,                                               TWOG1068
     +  1H ,10X,'I    Mass of ion                  ',F12.5,             TWOG1067
     +  ' GeV',11X,'I',/,                                               TWOG1068
     +  1H ,10X,'I    Gamma of ion                 ',F12.5,             TWOG1067
     +  '    ',11X,'I',/,                                               TWOG1068
     +  1H ,10X,'I    Charge of ion                ',F12.5,             TWOG1067
     +  '    ',11X,'I',/,                                               TWOG1068
     +  1H ,10X,'I    Atomic number of ion         ',F12.5,             TWOG1067
     +  '    ',11X,'I',/,                                               TWOG1068
     +  1H ,10X,'I    Minimum gamma-gamma mass     ',F12.5,             TWOG1069
     +  ' GeV',11X,'I',/,                                               TWOG1070
     +  1H ,10X,'I    Maximum gamma-gamma mass     ',F12.5,             TWOG1071
     +  ' GeV',11X,'I',/,                                               TWOG1080
     +  1H ,10X,'I    Minimum rapidity             ',F12.5,             TWOG1069
     +  '    ',11X,'I',/,                                               TWOG1070
     +  1H ,10X,'I    Maximum rapidity             ',F12.5,             TWOG1071
     +  '    ',11X,'I'  )                                               TWOG1080
      IF (LRES) THEN                                                    TWOG1081
      WRITE(LOUT,1001)                                                  TWOG1082
      WRITE(LOUT,1003)XMRES,XGRES                                       TWOG1083
      ENDIF                                                             TWOG1084
 1003 FORMAT                                                            TWOG1085
     +  (1H ,10X,'I    Resonance mass               ',F12.5,            TWOG1086
     +  ' GeV',11X,'I',/,                                               TWOG1087
     +  1H ,10X,'I    Resonance total width        ',F12.5,             TWOG1088
     +  ' GeV',11X,'I')                                                 TWOG1089
C
      IF (LWATE) THEN                                                   TWOG1090
      WRITE(LOUT,1001)                                                  TWOG1091
      WRITE(LOUT,1004)                                                  TWOG1092
      ELSE                                                              TWOG1095
      WRITE(LOUT,1001)                                                  TWOG1096
      WRITE(LOUT,1005)WGHTMX                                            TWOG1097
 1004 FORMAT(1H ,10X,'I    Produce weighted events      ',27X,'I')      TWOG1094
 1005 FORMAT(1H ,10X,'I    Maximum weight               ',F12.5,        TWOG1099
     +                                            '    ',11X,'I')       TWOG1100
      ENDIF                                                             TWOG1101
C                                                                       TWOG1102
      WRITE(LOUT,1001)                                                  TWOG1103
      WRITE(LOUT,1000)                                                  TWOG1104
C                                                                       TWOG1105
      DO 10 I=1,NCNTRS                                                  TWOG1106
   10 ICNTER(I)=0                                                       TWOG1107
      WGHTGT = 0.                                                       TWOG1108
      WGHTSM = 0.                                                       TWOG1109
      RETURN                                                            TWOG1110
      END                                                               TWOG1111
C=========================================================================
CDECK  ID>, DLDMX.  
      FUNCTION DLDMX(Rm,Ny)
      PARAMETER (isize=100000)
C                                                      S.A.Sadovsky
C   Calculate the Luminosity integral over dY            5.02.1995
C             in the rapidity range ymin-ymax
C             using luminosity interpolation table
C      Ny/2 - number of subintervals for Simpson integration
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      COMMON/lumfunct/anorm,cexp,ymin,ymax,ymn,ymx,nny,
     +                           xmin,xmax,amn,amx,nnms,alfun(isize)
      stpy=(ymx-ymn)/nny
      stpm=(amx-amn)/nnms                 ! Interpolation steps
C
      im=(rm-amn)/stpm
      q =(rm-im*stpm-amn)/stpm
C
      My = Ny/2
      Hy =(ymax-ymin)/My                  ! Integration step over Y
      ry = ymin
      iy =(ry-ymn)/stpy+1.
      p  =(ry-(iy-1)*stpy-ymn)/stpy       ! Four point interpolation:  1
      k  = iy+im*(nny+1)
      F1 =(1.-p)*(1.-q)*alfun(k)+p*(1.-q)*alfun(k+1)+
     *        q*(1.-p)*alfun(k+nny+1)+q*p*alfun(k+2+nny)
C
      SIM = 0.0
      DO 10 i=1,My
C
        ry = ymin + Hy*i
        iy =(ry-ymn)/stpy+1.
        p  =(ry-(iy-1)*stpy-ymn)/stpy       ! Four point interpolation:  3
        k  = iy+im*(nny+1)
        F3 =(1.-p)*(1.-q)*alfun(k)+p*(1.-q)*alfun(k+1)+
     *        q*(1.-p)*alfun(k+nny+1)+q*p*alfun(k+2+nny)
C
        ry = ry   - Hy*0.5
        iy =(ry-ymn)/stpy+1.
        p  =(ry-(iy-1)*stpy-ymn)/stpy       ! Four point interpolation:  2
        k  = iy+im*(nny+1)
        F2 =(1.-p)*(1.-q)*alfun(k)+p*(1.-q)*alfun(k+1)+
     *        q*(1.-p)*alfun(k+nny+1)+q*p*alfun(k+2+nny)
C
        SIM = SIM + (F1 + 4.*F2 + F3)/6.
   10 F1  = F3
C
      DLDMX = 1000.*SIM*HY                ! dL/dMx  in GeV^-1
      RETURN
      END
*=======================================================================
CDECK  ID>, LUMFUN. 
      SUBROUTINE LUMFUN(ifl,file,ymin,ymax,ny,xmn,xmx,nms)
C-----------------------------------------------------------------------
C preparation of luminosity function for generator
C
C     ifl=-1 - calculate function + save in file 'file'
C     ifl= 0 - calculate function
C     ifl= 1 - read function from file
C       ymin,ymax,ny - rapidity range & number of points
C       xmn,xmx,nms  - mass range & number of points
C       Attention: dimensions of alfun are alfun(ny+1,nms+1)
C
C     Called by GGINIT
C
C     Author: G.V.Khaoustov 12.12.1994
C     Modified: Yu.Kharlov 18-JUN-1997
C-----------------------------------------------------------------------
      COMMON /ggpar/ pi,hbarc,gev2nb,alpha, amprt(5), qf,nc, egg_max,
     &               gvconst(4,10)
      REAL nc
      COMMON/D2LParam/lz,la,Rion,Gamma
      INTEGER lz,la
      REAL    Rion,Gamma
      PARAMETER (isize=100000)
      CHARACTER*(*) file
      REAL work(1000)
      COMMON/lumfunct/anorm,cexp,umn,umx,ymn,ymx,nny,
     +                           wmn,wmx,amn,amx,nnms,alfun(isize)
C
      umn = ymin                    !  parameters for luminosity generation
      umx = ymax
      wmn = xmn
      wmx = xmx
C
C Define R : radii of the ions (in fm) (if set non-zero, use it)
      IF (Rion .LE. 0) THEN
        rion = 1.2*float(la)**(1.0/3.0)
      ENDIF
*
*       Fing max energy of gg-system (with safity factor 2)
*
      egg_max = 2.*hbarc * gamma / rion
*
      IF(ifl.EQ.1) THEN             !  read lum. function from file
        OPEN(unit=10,err=1,file=file,status='old',
     *                                        form='formatted')
        READ(10, *,err=3) anorm,cexp
        READ(10,99,err=3) ymn,ymx,nny,amn,amx,nnms
        READ(10,100,err=3) (alfun(I),i=1,(nny+1)*(nnms+1))
        CLOSE(unit=10)
        PRINT *,'The Luminosity has been read from File '
        RETURN
      ENDIF
c
      IF (ny*nms .GT. isize) GOTO 2
      PRINT *,'Start the Luminosity calculation'
      nny = ny
      ymn = ymin
      ymx = ymax
      nnms= nms
      amn = xmn
      amx = xmx
c
      stepy=(ymx-ymn)/nny
      stepm=(amx-amn)/nnms
*
*         Lum. function calculation
*
      DO IM = 0,nnms
        am=amn+im*stepm
*         Find Y-limits for mass AM (Yu.Kharlov)
        ym_max  = alog (egg_max/am)
        ym_max  = max(ym_max,5.*stepy)
        ygg_max = min (ymax, ym_max)
        ygg_min = max (ymin,-ym_max)
        DO IY = 1,nny+1
          Y = ymn+(iy-1)*stepy
          IF (Y.LT.ygg_max .AND. Y.GT.ygg_min) THEN
            dl = D2LDMDY(1000.*am,Y)   ! Differential Luminosity in MeV^-1
          ELSE
            dl = 0.0
          END IF
          alfun(iy+(nny+1)*im) = DL
        ENDDO
      ENDDO
c
      s=0.
      DO i=1,nnms+1   ! find max values for exp. parameters calcul.
        am=0.
        DO j=1,nny+1
          IF(am.LT.alfun(j+(nny+1)*(i-1))) am=alfun(j+(nny+1)*(i-1))
        ENDDO
        work(i)=am
        s=s+am
      ENDDO
      s=(s-work(nnms+1))*stepm
c
      smax=0.      ! calculation optimized parameters for exp. generator
      iflag=-1
      DO 101 i=1,nnms
        DO 102 j=i+1,nnms+1
          IF (work(j) .LE. 0.0) GOTO 102
          c=alog(work(i)/work(j))/((i-j)*stepm)
          am=work(i)/exp(c*(i-1)*stepm)
          se=am/c*(exp(c*(amx-amn))-1.)
          IF(s/se.lt.smax) GOTO 102
          DO k=1,nnms+1
            amp=am*exp(c*(k-1)*stepm)
            IF((work(k)-amp)/amp.GT.0.00001) THEN  !  Accuracy problem
              GOTO 102
            ENDIF
          ENDDO
          smax =s/se
          anorm=am
          cexp =c
          iflag=1
  102   CONTINUE
  101 CONTINUE
      IF(iflag.LT.0) THEN
        PRINT *,'Fault to find best exponential approximation'
        STOP
      ENDIF
c
      IF(ifl.EQ.0) RETURN        ! save lum. function
      OPEN(unit=10,err=1,file=file,status='unknown',
     *                                          form='formatted')
      WRITE(10, *,err=3) anorm,cexp
      WRITE(10,99,err=3) ymn,ymx,nny,amn,amx,nnms
      WRITE(10,100,err=3) (alfun(i),i=1,(nny+1)*(nnms+1))
      CLOSE(unit=10)
      PRINT *,'The Luminosity calculation is finished'
      RETURN
    1 PRINT *,' file open error, file:',file
      STOP
    2 PRINT *,' array "ALFUN" to small to contain luminosity function'
      STOP
    3 PRINT *,' file read/writr error, file:',file
      STOP
  100 FORMAT(5e13.5)
   99 FORMAT(2e13.5,i5,2e13.5,i5)
      END
C======================================================================
CDECK  ID>, D2LDMDY.
      REAL FUNCTION  D2LDMDY(M,Y)
C-----------------------------------------------------------------------
C Calculate the Luminosity-function d^2L/dMdY
C where M is the invariant mass and Y the rapidity
C of the photon-photon-system
C
C Calculations are based on the formulae of
C N. Baron, G. Baur, Phys. Rev. C
C C. Cahn, J. D. Jackson, Nucl. Phys.
C
C Kai Hencken, 29. 9. 1994
C
C M : invariant mass (MeV)
C Y : rapidity
C
C D2LDMDY : diff. Lum. (MeV**-1)
C
C     Called by LUMFUN
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/D2LParam/lz,la,Rion,Gamma
      INTEGER lz,la
      REAL    Rion,Gamma
      REAL M,Y
      REAL RM,W1,W2
      COMMON/D2LIParam/RM,W1,W2

      REAL NClass,DeltaL
      EXTERNAL NClass,DeltaL

C
C physical constants hbar*c and 1/alpha
C
      REAL hbarc,falphai
      PARAMETER (hbarc=197.327053,falphai=137.0359895)

      REAL LClass,DL
C
C calculate photon frequencies from M,Y
C
      W1 = M/2.0 * EXP ( Y)
      W2 = M/2.0 * EXP (-Y)
C
C calculate Rho in MeV^-1 instead of fm
C
      RM = Rion / hbarc
C
C calculate the "classical" value first.
C

      LClass = NClass (W1,RM,Gamma) * NClass (W2,RM,Gamma)

C
C substract from this the Delta-value
C
      DL = DeltaL ()

      D2LDMDY = 2.0 / M * float(LZ)**4 / falphai**2 * (LClass - DL)

      RETURN
      END
*=======================================================================
CDECK  ID>, NCLASS. 
      REAL FUNCTION NClass (W, R, Gamma)
C
C the classical photon-number without impact parameter cutoff:
C
      IMPLICIT NONE
      REAL W,R, Gamma

      REAL Pi
      PARAMETER (Pi = 3.141592)

      REAL BesselK0,BesselK1
      EXTERNAL BesselK0,BesselK1

      REAL Xi

      Xi = W*R/Gamma

      NClass = 2.0/PI * (Xi*BesselK0 (Xi)*BesselK1 (Xi) -
     &          Xi**2/2.0 *(BesselK1(Xi)**2 - BesselK0(Xi)**2))

      RETURN
      END
*=======================================================================
CDECK  ID>, DELTAL. 
      REAL FUNCTION DeltaL ()
C
C the difference to the classical value
C     Called by d2LdMdY
C
      IMPLICIT NONE
      COMMON/D2LParam/lz,la,Rion,Gamma
      INTEGER lz,la
      REAL    Rion,Gamma
      REAL RM,W1,W2
      COMMON/D2LIParam/RM,W1,W2

      DOUBLE PRECISION XGAUSS(126),WGAUSS(126)
      COMMON /GAUSSD/XGAUSS,WGAUSS

      INTEGER N,I
      REAL    b1
      REAL*8  t,tmin,tmax
      REAL*8  Sum,Int,Int2

      REAL Itgb1
      EXTERNAL Itgb1

      REAL   Pi
      REAL*8 Accur
      PARAMETER (Pi = 3.141592,Accur = 1D -2)

C
C integrate first over b1
C

C
C Loop incrementing the boundary
C
      tmin = 0.00
      tmax = 0.25
      Sum  = 0.00

  100 CONTINUE
C
C Loop for the Gauss integration
C
      Int=0.0
      DO N=1,6
        Int2 = Int
        Int=0.0
        DO I=2**N-1,2**(N+1)-2
          t   = (tmax-tmin)/2.0*XGAUSS(I) + (tmax+tmin)/2.0
          b1  = RM * DEXP (t)
          Int = Int + WGAUSS(I) * Itgb1(b1) * b1**2
        ENDDO
        Int = (tmax-tmin)/2.0*Int
        IF (Int .NE. 0.) THEN
          IF (DABS ((Int2-Int)/Int) .LT. Accur) GOTO 200
        END IF
      ENDDO
      IF (int .LT. 1.D-20) THEN
        int = 0.0D0
      ELSE
        WRITE(*,*) ' (b1) GAUSS MAY BE INACCURATE'
      END IF

  200 CONTINUE
      Sum = Sum + Int
      IF (sum .NE. 0.0) THEN
        IF (ABS (Int2/Sum) .LT. Accur) GOTO 300
      ELSE
        IF (b1 .GT. 1.E+06)            GOTO 300
      END IF
      tmin = tmax
      tmax = tmax + 0.5
      GOTO 100

  300 CONTINUE

      DeltaL = 4.0*Pi * Sum
      RETURN
      END

CDECK  ID>, ITGB1.  
      REAL FUNCTION Itgb1(b)
C
C integration then over b2
C     Called by DeltaL
C
C                                  Corrected by S.A.Sadovsky
C                                             25.10.1994
      IMPLICIT NONE
c  +seq,d2lpar. Not used Yu.Kh. 18-JUN-1997
      REAL     b
      REAL*8   b1,bmin,bmax
      REAL     RM,W1,W2
      COMMON/D2LIParam/RM,W1,W2
      DOUBLE PRECISION XGAUSS(126),WGAUSS(126)
      COMMON /GAUSSD/XGAUSS,WGAUSS
      INTEGER  N,I
      REAL*8   b2,Int,Int2
      REAL*8   Itgb2
      EXTERNAL Itgb2
      REAL*8   Accur
      PARAMETER (Accur = 1.D-2)

      b1 = b
      bmin = b1 - 2.0*RM
      IF (RM .GT. bmin) THEN
        bmin = RM
      ENDIF
      bmax = b1 + 2.0 * RM
      Int = 0.0
      DO N=1,6
        Int2 = Int
        Int = 0.0
        DO I=2**N-1,2**(N+1)-2
          b2 = (bmax-bmin)/2.0*XGAUSS(I) + (bmax+bmin)/2.0
          Int = Int + WGAUSS(I) * Itgb2(b1,b2) * b2
        END DO
        Int = (bmax-bmin)/2.0*Int
        IF (Int.NE.0.) THEN
          IF (DABS((Int2 - Int)/Int) .LT. Accur) GOTO 100
        ENDIF
      ENDDO
      IF (int .LT. 1.D-20) THEN
        int = 0.0D0
      ELSE
        WRITE(*,*) ' (b2) GAUSS MAY BE INACCURATE, Itgb1=',Int
      END IF
  100 CONTINUE
      Itgb1 = Int
      RETURN
      END

CDECK  ID>, ITGB2.  
      REAL*8 FUNCTION Itgb2 (b1,b2)
C
C       The function to be integrated over
C       Called by Itgb1
C
C                                  Corrected by S.A.Sadovsky
C                                             25.10.1994
      IMPLICIT NONE
      COMMON/D2LParam/lz,la,Rion,Gamma
      INTEGER lz,la
      REAL    Rion,Gamma
      REAL*8 b1,b2,arg
      REAL RM,W1,W2
      COMMON/D2LIParam/RM,W1,W2

      REAL Nphoton
      EXTERNAL Nphoton

      arg = (b1**2 + b2**2 - 4.0 * RM**2) / (2.0*b1*b2)
      IF (arg .GT. 1.D0) arg = 1.D0
      Itgb2 = Nphoton(W1,sngl(b1),Gamma) *
     &        Nphoton(W2,sngl(b2),Gamma) *
     &        DACOS (arg)

      RETURN
      END

CDECK  ID>, NPHOTON.
      REAL FUNCTION Nphoton (W,Rho,Gamma)
C
C       The differential photonnumber for a nucleus
C       (without form factor)
C
      IMPLICIT NONE
      REAL W,Rho,Gamma

      REAL BESSELK1
      EXTERNAL BESSELK1

      REAL Wphib,WGamma

      REAL PI
      PARAMETER (PI=3.141592)

      WGamma = W/Gamma
C
C without form factor
C
      Wphib = WGamma*BESSELK1 (WGamma*Rho)

      Nphoton = 1.0/PI**2 * Wphib**2

      RETURN
      END

CDECK  ID>, BESSELK0.   
      REAL FUNCTION BESSELK0(X)
C
C  Special functions I0,K0,I1,K1
C  see Abramowitz&Stegun
C
      IMPLICIT NONE
      REAL X,Y

      REAL BESSELI0
      EXTERNAL BESSELI0

      IF (X .LT. 2.0) THEN

        Y = X**2/4.0

        BESSELK0 =
     &    (-LOG(X/2.0)*BESSELI0(X))+(-.57721566+Y*(0.42278420
     &    +Y*(0.23069756+Y*(0.3488590E-1+Y*(0.262698E-2
     &    +Y*(0.10750E-3+Y*0.740E-5))))))

      ELSE

        Y = 2.0/X

        BESSELK0 =
     &    (EXP(-X)/SQRT(X))*(1.25331414+Y*(-0.7832358E-1
     &    +Y*(0.2189568E-1+Y*(-0.1062446E-1+Y*(0.587872E-2
     &    +Y*(-0.251540E-2+Y*0.53208E-3))))))

      ENDIF

      RETURN
      END


CDECK  ID>, BESSELI0.   
      REAL FUNCTION BESSELI0(X)
      IMPLICIT NONE

      REAL X,Y,AX

      AX = ABS(X)

      IF (AX .LT. 3.75) THEN

        Y = (X/3.75)**2

        BESSELI0 =
     &    1.0+Y*(3.5156229+Y*(3.0899424+Y*(1.2067492
     &    +Y*(0.2659732+Y*(0.360768E-1+Y*0.45813E-2)))))

      ELSE

        Y = 3.75/AX

        BESSELI0 =
     &    (EXP(AX)/SQRT(AX))*(0.39894228+Y*(0.1328592E-1
     &    +Y*(0.225319E-2+Y*(-0.157565E-2+Y*(0.916281E-2
     &    +Y*(-0.2057706E-1+Y*(0.2635537E-1+Y*(-0.1647633E-1
     &    +Y*0.392377E-2))))))))

      ENDIF

      RETURN
      END


CDECK  ID>, BESSELK1.   
      REAL FUNCTION BESSELK1(X)
      IMPLICIT NONE

      REAL X,Y

      REAL BESSELI1
      EXTERNAL BESSELI1

      IF (X .LT. 2.0) THEN

        Y = X**2/4.0
        BESSELK1 =
     &    (LOG(X/2.0)*BESSELI1(X))+(1.0/X)*(1.0+Y*(0.15443144
     &    +Y*(-0.67278579+Y*(-0.18156897+Y*(-0.1919402E-1
     &    +Y*(-0.110404E-2+Y*(-0.4686E-4)))))))

      ELSE

        Y=2.0/X
        BESSELK1 =
     &    (EXP(-X)/SQRT(X))*(1.25331414+Y*(0.23498619
     &    +Y*(-0.3655620E-1+Y*(0.1504268E-1+Y*(-0.780353E-2
     &    +Y*(0.325614E-2+Y*(-0.68245E-3)))))))

      ENDIF

      RETURN
      END

CDECK  ID>, BESSELI1.   
      REAL FUNCTION BESSELI1(X)
      IMPLICIT NONE

      REAL X,Y,AX

      AX = ABS(X)

      IF (AX .LT. 3.75) THEN

        Y = (X/3.75)**2

        BESSELI1 =
     &    AX*(0.5+Y*(0.87890594+Y*(0.51498869+Y*(0.15084934
     &    +Y*(0.2658733E-1+Y*(0.301532E-2+Y*0.32411E-3))))))


      ELSE

        Y = 3.75/AX

        BESSELI1 =
     &    0.2282967E-1+Y*(-0.2895312E-1+Y*(0.1787654E-1
     &    -Y*0.420059E-2))

        BESSELI1 =
     &    0.39894228+Y*(-0.3988024E-1+Y*(-0.362018E-2
     &    +Y*(0.163801E-2+Y*(-0.1031555E-1+Y*BESSELI1))))

        BESSELI1 = BESSELI1 *
     &    EXP(AX)/SQRT(AX)

      ENDIF

      IF (X .LT. 0) THEN
        BESSELI1 = -BESSELI1
      ENDIF

      RETURN
      END

CDECK  ID>, GAUSSDAT.   
C DATA for the Gauss integration, x-values and weight for a
C N=2,4,8,16,32,64 point Gauss integration.
C
C based on program GAUSS64 of D. Trautmann
C data from Abramowitz & Stegun
C
C KAI HENCKEN, FEBRUAR 1993
C
      BLOCKDATA GAUSSDATA
      IMPLICIT NONE
      DOUBLE PRECISION XGAUSS(126),WGAUSS(126)
      COMMON /GAUSSD/XGAUSS,WGAUSS

      DATA XGAUSS(1)/ .57735026918962576D0/
      DATA XGAUSS(2)/-.57735026918962576D0/
      DATA WGAUSS(1)/ 1.00000000000000000D0/
      DATA WGAUSS(2)/ 1.00000000000000000D0/

      DATA XGAUSS(3)/ .33998104358485627D0/
      DATA XGAUSS(4)/ .86113631159405258D0/
      DATA XGAUSS(5)/-.33998104358485627D0/
      DATA XGAUSS(6)/-.86113631159405258D0/
      DATA WGAUSS(3)/ .65214515486254613D0/
      DATA WGAUSS(4)/ .34785484513745385D0/
      DATA WGAUSS(5)/ .65214515486254613D0/
      DATA WGAUSS(6)/ .34785484513745385D0/

      DATA XGAUSS(7)/ .18343464249564981D0/
      DATA XGAUSS(8)/ .52553240991632899D0/
      DATA XGAUSS(9)/ .79666647741362674D0/
      DATA XGAUSS(10)/ .96028985649753623D0/
      DATA XGAUSS(11)/-.18343464249564981D0/
      DATA XGAUSS(12)/-.52553240991632899D0/
      DATA XGAUSS(13)/-.79666647741362674D0/
      DATA XGAUSS(14)/-.96028985649753623D0/
      DATA WGAUSS(7)/ .36268378337836198D0/
      DATA WGAUSS(8)/ .31370664587788727D0/
      DATA WGAUSS(9)/ .22238103445337448D0/
      DATA WGAUSS(10)/ .10122853629037627D0/
      DATA WGAUSS(11)/ .36268378337836198D0/
      DATA WGAUSS(12)/ .31370664587788727D0/
      DATA WGAUSS(13)/ .22238103445337448D0/
      DATA WGAUSS(14)/ .10122853629037627D0/

      DATA XGAUSS(15)/ .0950125098376374402D0/
      DATA XGAUSS(16)/ .281603550779258913D0/
      DATA XGAUSS(17)/ .458016777657227386D0/
      DATA XGAUSS(18)/ .617876244402643748D0/
      DATA XGAUSS(19)/ .755404408355003034D0/
      DATA XGAUSS(20)/ .865631202387831744D0/
      DATA XGAUSS(21)/ .944575023073232576D0/
      DATA XGAUSS(22)/ .989400934991649933D0/
      DATA XGAUSS(23)/-.0950125098376374402D0/
      DATA XGAUSS(24)/-.281603550779258913D0/
      DATA XGAUSS(25)/-.458016777657227386D0/
      DATA XGAUSS(26)/-.617876244402643748D0/
      DATA XGAUSS(27)/-.755404408355003034D0/
      DATA XGAUSS(28)/-.865631202387831744D0/
      DATA XGAUSS(29)/-.944575023073232576D0/
      DATA XGAUSS(30)/-.989400934991649933D0/
      DATA WGAUSS(15)/ .189450610455068496D0/
      DATA WGAUSS(16)/ .182603415044923589D0/
      DATA WGAUSS(17)/ .169156519395002538D0/
      DATA WGAUSS(18)/ .149595988816576732D0/
      DATA WGAUSS(19)/ .124628971255533872D0/
      DATA WGAUSS(20)/ .0951585116824927848D0/
      DATA WGAUSS(21)/ .0622535239386478929D0/
      DATA WGAUSS(22)/ .0271524594117540949D0/
      DATA WGAUSS(23)/ .189450610455068496D0/
      DATA WGAUSS(24)/ .182603415044923589D0/
      DATA WGAUSS(25)/ .169156519395002538D0/
      DATA WGAUSS(26)/ .149595988816576732D0/
      DATA WGAUSS(27)/ .124628971255533872D0/
      DATA WGAUSS(28)/ .0951585116824927848D0/
      DATA WGAUSS(29)/ .0622535239386478929D0/
      DATA WGAUSS(30)/ .0271524594117540949D0/

      DATA XGAUSS(31)/ .0483076656877383162D0/
      DATA XGAUSS(32)/ .144471961582796493D0/
      DATA XGAUSS(33)/ .239287362252137075D0/
      DATA XGAUSS(34)/ .331868602282127650D0/
      DATA XGAUSS(35)/ .421351276130635345D0/
      DATA XGAUSS(36)/ .506899908932229390D0/
      DATA XGAUSS(37)/ .587715757240762329D0/
      DATA XGAUSS(38)/ .663044266930215201D0/
      DATA XGAUSS(39)/ .732182118740289680D0/
      DATA XGAUSS(40)/ .794483795967942407D0/
      DATA XGAUSS(41)/ .849367613732569970D0/
      DATA XGAUSS(42)/ .896321155766052124D0/
      DATA XGAUSS(43)/ .934906075937739689D0/
      DATA XGAUSS(44)/ .964762255587506430D0/
      DATA XGAUSS(45)/ .985611511545268335D0/
      DATA XGAUSS(46)/ .997263861849481564D0/
      DATA XGAUSS(47)/-.0483076656877383162D0/
      DATA XGAUSS(48)/-.144471961582796493D0/
      DATA XGAUSS(49)/-.239287362252137075D0/
      DATA XGAUSS(50)/-.331868602282127650D0/
      DATA XGAUSS(51)/-.421351276130635345D0/
      DATA XGAUSS(52)/-.506899908932229390D0/
      DATA XGAUSS(53)/-.587715757240762329D0/
      DATA XGAUSS(54)/-.663044266930215201D0/
      DATA XGAUSS(55)/-.732182118740289680D0/
      DATA XGAUSS(56)/-.794483795967942407D0/
      DATA XGAUSS(57)/-.849367613732569970D0/
      DATA XGAUSS(58)/-.896321155766052124D0/
      DATA XGAUSS(59)/-.934906075937739689D0/
      DATA XGAUSS(60)/-.964762255587506430D0/
      DATA XGAUSS(61)/-.985611511545268335D0/
      DATA XGAUSS(62)/-.997263861849481564D0/
      DATA WGAUSS(31)/ .0965400885147278006D0/
      DATA WGAUSS(32)/ .0956387200792748594D0/
      DATA WGAUSS(33)/ .0938443990808045654D0/
      DATA WGAUSS(34)/ .0911738786957638847D0/
      DATA WGAUSS(35)/ .0876520930044038111D0/
      DATA WGAUSS(36)/ .0833119242269467552D0/
      DATA WGAUSS(37)/ .0781938957870703065D0/
      DATA WGAUSS(38)/ .0723457941088485062D0/
      DATA WGAUSS(39)/ .0658222227763618468D0/
      DATA WGAUSS(40)/ .0586840934785355471D0/
      DATA WGAUSS(41)/ .0509980592623761762D0/
      DATA WGAUSS(42)/ .0428358980222266807D0/
      DATA WGAUSS(43)/ .0342738629130214331D0/
      DATA WGAUSS(44)/ .0253920653092620595D0/
      DATA WGAUSS(45)/ .0162743947309056706D0/
      DATA WGAUSS(46)/ .00701861000947009660D0/
      DATA WGAUSS(47)/ .0965400885147278006D0/
      DATA WGAUSS(48)/ .0956387200792748594D0/
      DATA WGAUSS(49)/ .0938443990808045654D0/
      DATA WGAUSS(50)/ .0911738786957638847D0/
      DATA WGAUSS(51)/ .0876520930044038111D0/
      DATA WGAUSS(52)/ .0833119242269467552D0/
      DATA WGAUSS(53)/ .0781938957870703065D0/
      DATA WGAUSS(54)/ .0723457941088485062D0/
      DATA WGAUSS(55)/ .0658222227763618468D0/
      DATA WGAUSS(56)/ .0586840934785355471D0/
      DATA WGAUSS(57)/ .0509980592623761762D0/
      DATA WGAUSS(58)/ .0428358980222266807D0/
      DATA WGAUSS(59)/ .0342738629130214331D0/
      DATA WGAUSS(60)/ .0253920653092620595D0/
      DATA WGAUSS(61)/ .0162743947309056706D0/
      DATA WGAUSS(62)/ .00701861000947009660D0/

      DATA XGAUSS(63)/ .02435029266342443250D0/
      DATA XGAUSS(64)/ .0729931217877990394D0/
      DATA XGAUSS(65)/ .121462819296120554D0/
      DATA XGAUSS(66)/ .169644420423992818D0/
      DATA XGAUSS(67)/ .217423643740007084D0/
      DATA XGAUSS(68)/ .264687162208767416D0/
      DATA XGAUSS(69)/ .311322871990210956D0/
      DATA XGAUSS(70)/ .357220158337668116D0/
      DATA XGAUSS(71)/ .402270157963991604D0/
      DATA XGAUSS(72)/ .446366017253464088D0/
      DATA XGAUSS(73)/ .489403145707052957D0/
      DATA XGAUSS(74)/ .531279464019894546D0/
      DATA XGAUSS(75)/ .571895646202634034D0/
      DATA XGAUSS(76)/ .611155355172393250D0/
      DATA XGAUSS(77)/ .648965471254657340D0/
      DATA XGAUSS(78)/ .685236313054233243D0/
      DATA XGAUSS(79)/ .719881850171610827D0/
      DATA XGAUSS(80)/ .752819907260531897D0/
      DATA XGAUSS(81)/ .783972358943341408D0/
      DATA XGAUSS(82)/ .813265315122797560D0/
      DATA XGAUSS(83)/ .840629296252580363D0/
      DATA XGAUSS(84)/ .865999398154092820D0/
      DATA XGAUSS(85)/ .889315445995114106D0/
      DATA XGAUSS(86)/ .910522137078502806D0/
      DATA XGAUSS(87)/ .929569172131939576D0/
      DATA XGAUSS(88)/ .946411374858402816D0/
      DATA XGAUSS(89)/ .961008799652053719D0/
      DATA XGAUSS(90)/ .973326827789910964D0/
      DATA XGAUSS(91)/ .983336253884625957D0/
      DATA XGAUSS(92)/ .991013371476744321D0/
      DATA XGAUSS(93)/ .996340116771955279D0/
      DATA XGAUSS(94)/ .999305041735772139D0/
      DATA XGAUSS(95)/-.02435029266342443250D0/
      DATA XGAUSS(96)/-.0729931217877990394D0/
      DATA XGAUSS(97)/-.121462819296120554D0/
      DATA XGAUSS(98)/-.169644420423992818D0/
      DATA XGAUSS(99)/-.217423643740007084D0/
      DATA XGAUSS(100)/-.264687162208767416D0/
      DATA XGAUSS(101)/-.311322871990210956D0/
      DATA XGAUSS(102)/-.357220158337668116D0/
      DATA XGAUSS(103)/-.402270157963991604D0/
      DATA XGAUSS(104)/-.446366017253464088D0/
      DATA XGAUSS(105)/-.489403145707052957D0/
      DATA XGAUSS(106)/-.531279464019894546D0/
      DATA XGAUSS(107)/-.571895646202634034D0/
      DATA XGAUSS(108)/-.611155355172393250D0/
      DATA XGAUSS(109)/-.648965471254657340D0/
      DATA XGAUSS(110)/-.685236313054233243D0/
      DATA XGAUSS(111)/-.719881850171610827D0/
      DATA XGAUSS(112)/-.752819907260531897D0/
      DATA XGAUSS(113)/-.783972358943341408D0/
      DATA XGAUSS(114)/-.813265315122797560D0/
      DATA XGAUSS(115)/-.840629296252580363D0/
      DATA XGAUSS(116)/-.865999398154092820D0/
      DATA XGAUSS(117)/-.889315445995114106D0/
      DATA XGAUSS(118)/-.910522137078502806D0/
      DATA XGAUSS(119)/-.929569172131939576D0/
      DATA XGAUSS(120)/-.946411374858402816D0/
      DATA XGAUSS(121)/-.961008799652053719D0/
      DATA XGAUSS(122)/-.973326827789910964D0/
      DATA XGAUSS(123)/-.983336253884625957D0/
      DATA XGAUSS(124)/-.991013371476744321D0/
      DATA XGAUSS(125)/-.996340116771955279D0/
      DATA XGAUSS(126)/-.999305041735772139D0/
      DATA WGAUSS(63)/ .0486909570091397204D0/
      DATA WGAUSS(64)/ .0485754674415034269D0/
      DATA WGAUSS(65)/ .0483447622348029572D0/
      DATA WGAUSS(66)/ .0479993885964583077D0/
      DATA WGAUSS(67)/ .0475401657148303087D0/
      DATA WGAUSS(68)/ .0469681828162100173D0/
      DATA WGAUSS(69)/ .0462847965813144172D0/
      DATA WGAUSS(70)/ .0454916279274181445D0/
      DATA WGAUSS(71)/ .0445905581637565631D0/
      DATA WGAUSS(72)/ .0435837245293234534D0/
      DATA WGAUSS(73)/ .0424735151236535890D0/
      DATA WGAUSS(74)/ .0412625632426235286D0/
      DATA WGAUSS(75)/ .0399537411327203414D0/
      DATA WGAUSS(76)/ .0385501531786156291D0/
      DATA WGAUSS(77)/ .0370551285402400460D0/
      DATA WGAUSS(78)/ .0354722132568823838D0/
      DATA WGAUSS(79)/ .0338051618371416094D0/
      DATA WGAUSS(80)/ .0320579283548515535D0/
      DATA WGAUSS(81)/ .0302346570724024789D0/
      DATA WGAUSS(82)/ .0283396726142594832D0/
      DATA WGAUSS(83)/ .0263774697150546587D0/
      DATA WGAUSS(84)/ .0243527025687108733D0/
      DATA WGAUSS(85)/ .0222701738083832542D0/
      DATA WGAUSS(86)/ .0201348231535302094D0/
      DATA WGAUSS(87)/ .0179517157756973431D0/
      DATA WGAUSS(88)/ .0157260304760247193D0/
      DATA WGAUSS(89)/ .0134630478967186426D0/
      DATA WGAUSS(90)/ .0111681394601311288D0/
      DATA WGAUSS(91)/ .00884675982636394772D0/
      DATA WGAUSS(92)/ .00650445796897836286D0/
      DATA WGAUSS(93)/ .00414703326056246764D0/
      DATA WGAUSS(94)/ .00178328072169643295D0/
      DATA WGAUSS(95)/ .0486909570091397204D0/
      DATA WGAUSS(96)/ .0485754674415034269D0/
      DATA WGAUSS(97)/ .0483447622348029572D0/
      DATA WGAUSS(98)/ .0479993885964583077D0/
      DATA WGAUSS(99)/ .0475401657148303087D0/
      DATA WGAUSS(100)/ .0469681828162100173D0/
      DATA WGAUSS(101)/ .0462847965813144172D0/
      DATA WGAUSS(102)/ .0454916279274181445D0/
      DATA WGAUSS(103)/ .0445905581637565631D0/
      DATA WGAUSS(104)/ .0435837245293234534D0/
      DATA WGAUSS(105)/ .0424735151236535890D0/
      DATA WGAUSS(106)/ .0412625632426235286D0/
      DATA WGAUSS(107)/ .0399537411327203414D0/
      DATA WGAUSS(108)/ .0385501531786156291D0/
      DATA WGAUSS(109)/ .0370551285402400460D0/
      DATA WGAUSS(110)/ .0354722132568823838D0/
      DATA WGAUSS(111)/ .0338051618371416094D0/
      DATA WGAUSS(112)/ .0320579283548515535D0/
      DATA WGAUSS(113)/ .0302346570724024789D0/
      DATA WGAUSS(114)/ .0283396726142594832D0/
      DATA WGAUSS(115)/ .0263774697150546587D0/
      DATA WGAUSS(116)/ .0243527025687108733D0/
      DATA WGAUSS(117)/ .0222701738083832542D0/
      DATA WGAUSS(118)/ .0201348231535302094D0/
      DATA WGAUSS(119)/ .0179517157756973431D0/
      DATA WGAUSS(120)/ .0157260304760247193D0/
      DATA WGAUSS(121)/ .0134630478967186426D0/
      DATA WGAUSS(122)/ .0111681394601311288D0/
      DATA WGAUSS(123)/ .00884675982636394772D0/
      DATA WGAUSS(124)/ .00650445796897836286D0/
      DATA WGAUSS(125)/ .00414703326056246764D0/
      DATA WGAUSS(126)/ .00178328072169643295D0/

      END
*=======================================================================
      SUBROUTINE LORENB (U,PS,PI,PF)
C
C CERN PROGLIB# U102    LORENB          .VERSION KERNFOR  4.04  821124
C ORIG. 20/08/75 L.PAPE
C
      DOUBLE PRECISION PF4, FN
      DIMENSION      PS(4),PI(4),PF(4)

      IF (PS(4).EQ.U) GO TO 17
      PF4  = (PI(4)*PS(4)+PI(3)*PS(3)+PI(2)*PS(2)+PI(1)*PS(1)) / U
      FN   = (PF4+PI(4)) / (PS(4)+U)
      PF(1)= PI(1) + FN*PS(1)
      PF(2)= PI(2) + FN*PS(2)
      PF(3)= PI(3) + FN*PS(3)
      PF(4)= PF4
      GO TO 18
C
   17 PF(1)= PI(1)
      PF(2)= PI(2)
      PF(3)= PI(3)
      PF(4)= PI(4)
C
   18 CONTINUE
C
      RETURN
C
      END
