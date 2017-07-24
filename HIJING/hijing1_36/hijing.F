* $Id$
C     Version 1.36
C     Nothing important has been changed here. A few 'garbage' has been
C     cleaned up here, like common block HIJJET3 for the sea quark strings
C     which were originally created to implement the DPM scheme which
C     later was abadoned in the final version. The lines which operate
C     on these data are also deleted in the program.
C
C
C     Version 1.35
C     There are some changes in the program: subroutine HARDJET is now
C     consolidated with HIJHRD. HARDJET is used to re-initiate PYTHIA
C     for the triggered hard processes. Now that is done  altogether
C     with other normal hard processes in modified JETINI. In the new
C     version one calls JETINI every time one calls HIJHRD. In the new
C     version the effect of the isospin of the nucleon on hard processes,
C     especially direct photons is correctly considered.
C     For A+A collisions, one has to initilize pythia
C     separately for each type of collisions, pp, pn,np and nn,
C     or hp and hn for hA collisions. In JETINI we use the following
C     catalogue for different types of collisions:
C     h+h: h+h (I_TYPE=1)
C     h+A: h+p (I_TYPE=1), h+n (I_TYPE=2)
C     A+h: p+h (I_TYPE=1), n+h (I_TYPE=2)
C     A+A: p+p (I_TYPE=1), p+n (I_TYPE=2), n+p (I_TYPE=3), n+n (I_TYPE=4)
C*****************************************************************
c
C
C     Version 1.34
C     Last modification on January 5, 1998. Two misstakes are corrected in
C     function G. A Misstake in the subroutine Parton is also corrected.
C     (These are pointed out by Ysushi Nara).
C
C
C       Last modifcation on April 10, 1996. To conduct final
C       state radiation, PYTHIA reorganize the two scattered
C       partons and their final momenta will be a little
C       different. The summed total momenta of the partons
C       from the final state radiation are stored in HINT1(26-29)
C       and HINT1(36-39) which are little different from 
C       HINT1(21-24) and HINT1(41-44).
C
C       Version 1.33
C
C       Last modfication  on September 11, 1995. When HIJING and
C       PYTHIA are initialized, the shadowing is evaluated at
C       b=0 which is the maximum. This will cause overestimate
C       of shadowing for peripheral interactions. To correct this
C       problem, shadowing is set to zero when initializing. Then
C       use these maximum  cross section without shadowing as a
C       normalization of the Monte Carlo. This however increase
C       the computing time. IHNT2(16) is used to indicate whether
C       the sturcture function is called for (IHNT2(16)=1) initialization
C       or for (IHNT2(16)=0)normal collisions simulation
C
C       Last modification on Aagust 28, 1994. Two bugs associate
C       with the impact parameter dependence of the shadowing is
C       corrected.
C
C
c       Last modification on October 14, 1994. One bug is corrected
c       in the direct photon production option in subroutine
C       HIJHRD.( this problem was reported by Jim Carroll and Mike Beddo).
C       Another bug associated with keeping the decay history
C       in the particle information is also corrected.(this problem
C       was reported by Matt Bloomer)
C
C
C       Last modification on July 15, 1994. The option to trig on
C       heavy quark production (charm IHPR2(18)=0 or beauty IHPR2(18)=1) 
C       is added. To do this, set IHPR2(3)=3. For inclusive production,
C       one should reset HIPR1(10)=0.0. One can also trig larger pt
C       QQbar production by giving HIPR1(10) a nonvanishing value.
C       The mass of the heavy quark in the calculation of the cross
C       section (HINT1(59)--HINT1(65)) is given by HIPR1(7) (the
C       default is the charm mass D=1.5). We also include a separate
C       K-factor for heavy quark and direct photon production by
C       HIPR1(23)(D=2.0).
C
C       Last modification on May 24, 1994.  The option to
C       retain the information of all particles including those
C       who have decayed is IHPR(21)=1 (default=0). KATT(I,3) is 
C       added to contain the line number of the parent particle 
C       of the current line which is produced via a decay. 
C       KATT(I,4) is the status number of the particle: 11=particle
C       which has decayed; 1=finally produced particle.
C
C
C       Last modification on May 24, 1994( in HIJSFT when valence quark
C       is quenched, the following error is corrected. 1.2*IHNT2(1) --> 
C       1.2*IHNT2(1)**0.333333, 1.2*IHNT2(3) -->1.2*IHNT(3)**0.333333)
C
C
C       Last modification on March 16, 1994 (heavy flavor production
C       processes MSUB(81)=1 MSUB(82)=1 have been switched on,
C       charm production is the default, B-quark option is
C       IHPR2(18), when it is switched on, charm quark is 
C       automatically off)
C
C
C       Last modification on March 23, 1994 (an error is corrected
C       in the impact parameter dependence of the jet cross section)
C
C       Last modification Oct. 1993 to comply with non-vax
C       machines' compiler 
C
C*********************************************
C	LAST MODIFICATION April 5, 1991
CQUARK DISTRIBUTIOIN (1-X)**A/(X**2+C**2/S)**B 
C(A=HIPR1(44),B=HIPR1(46),C=HIPR1(45))
C STRING FLIP, VENUS OPTION IHPR2(15)=1,IN WHICH ONE CAN HAVE ONE AND
C TWO COLOR CHANGES, (1-W)**2,W*(1-W),W*(1-W),AND W*2, W=HIPR1(18), 
C AMONG PT DISTRIBUTION OF SEA QUARKS IS CONTROLLED BY HIPR1(42)
C
C	gluon jets can form a single string system
C
C	initial state radiation is included
C	
C	all QCD subprocesses are included
c
c	direct particles production is included(currently only direct
C		photon)
c
C	Effect of high P_T trigger bias on multiple jets distribution
c
C******************************************************************
C	                        HIJING.10                         *
C	          Heavy Ion Jet INteraction Generator        	  *
C	                           by                       	  *
C		   X. N. Wang      and   M. Gyulassy           	  *
C	 	      Lawrence Berkeley Laboratory		  *
C								  *
C******************************************************************
C
C******************************************************************
C NFP(K,1),NFP(K,2)=flavor of q and di-q, NFP(K,3)=present ID of  *
C proj, NFP(K,4) original ID of proj.  NFP(K,5)=colli status(0=no,*
C 1=elastic,2=the diffrac one in single-diffrac,3= excited string.*
C |NFP(K,6)| is the total # of jet production, if NFP(K,6)<0 it   *
C can not produce jet anymore. NFP(K,10)=valence quarks scattering*
C (0=has not been,1=is going to be, -1=has already been scattered *
C NFP(k,11) total number of interactions this proj has suffered   *
C PP(K,1)=PX,PP(K,2)=PY,PP(K,3)=PZ,PP(K,4)=E,PP(K,5)=M(invariant  *
C mass), PP(K,6,7),PP(K,8,9)=transverse momentum of quark and     *
C diquark,PP(K,10)=PT of the hard scattering between the valence  *
C quarks; PP(K,14,15)=the mass of quark,diquark.       		  * 
C******************************************************************
C
C****************************************************************
C
C	SUBROUTINE HIJING
C
C****************************************************************
	SUBROUTINE HIJING(FRAME,BMIN0,BMAX0)
	CHARACTER FRAME*8
	DIMENSION SCIP(300,300),RNIP(300,300),SJIP(300,300),JTP(3),
     &			IPCOL(90000),ITCOL(90000),SCIP2(300,300)
#define BLANKET_SAVE
#include "hiparnt.inc"
C
#include "hijcrdn.inc"
#include "hijglbr.inc"
#include "himain1.inc"
#include "himain2.inc"
#include "histrng.inc"
#include "hijjet1.inc"
#include "hijjet2.inc"
#include "hijjet4.inc"
C
#include "lujets_hijing.inc"
#include "ludat1_hijing.inc"
        SAVE

	BMAX=MIN(BMAX0,HIPR1(34)+HIPR1(35))
	BMIN=MIN(BMIN0,BMAX)
	IF(IHNT2(1).LE.1 .AND. IHNT2(3).LE.1) THEN
		BMIN=0.0
		BMAX=2.5*SQRT(HIPR1(31)*0.1/HIPR1(40))
	ENDIF
C			********HIPR1(31) is in mb =0.1fm**2
C*******THE FOLLOWING IS TO SELECT THE COORDINATIONS OF NUCLEONS 
C       BOTH IN PROJECTILE AND TARGET NUCLEAR( in fm)
C
	YP(1,1)=0.0
	YP(2,1)=0.0
	YP(3,1)=0.0
	IF(IHNT2(1).LE.1) GO TO 14
	DO 10 KP=1,IHNT2(1)
5	R=HIRND(1)
c
        if(IHNT2(1).EQ.2) then
           rnd1=max(RLU_HIJING(NSEED),1.0e-20)
           rnd2=max(RLU_HIJING(NSEED),1.0e-20)
           rnd3=max(RLU_HIJING(NSEED),1.0e-20)
           R=-0.5*(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
        endif
c
	X=RLU_HIJING(0)
	CX=2.0*X-1.0
	SX=SQRT(1.0-CX*CX)
C		********choose theta from uniform cos(theta) distr
	PHI=RLU_HIJING(0)*2.0*HIPR1(40)
C		********choose phi form uniform phi distr 0 to 2*pi
	YP(1,KP)=R*SX*COS(PHI)
	YP(2,KP)=R*SX*SIN(PHI)
	YP(3,KP)=R*CX
	IF(HIPR1(29).EQ.0.0) GO TO 10
	DO 8  KP2=1,KP-1
		DNBP1=(YP(1,KP)-YP(1,KP2))**2
		DNBP2=(YP(2,KP)-YP(2,KP2))**2
		DNBP3=(YP(3,KP)-YP(3,KP2))**2
		DNBP=DNBP1+DNBP2+DNBP3
		IF(DNBP.LT.HIPR1(29)*HIPR1(29)) GO TO 5
C			********two neighbors cannot be closer than 
C				HIPR1(29)
8	CONTINUE
10	CONTINUE
c*******************************
        if(IHNT2(1).EQ.2) then
           YP(1,2)=-YP(1,1)
           YP(2,2)=-YP(2,1)
           YP(3,2)=-YP(3,1)
        endif
c********************************
	DO 12 I=1,IHNT2(1)-1
	DO 12 J=I+1,IHNT2(1)
	IF(YP(3,I).GT.YP(3,J)) GO TO 12
	Y1=YP(1,I)
	Y2=YP(2,I)
	Y3=YP(3,I)
	YP(1,I)=YP(1,J)
	YP(2,I)=YP(2,J)
	YP(3,I)=YP(3,J)
	YP(1,J)=Y1
	YP(2,J)=Y2
	YP(3,J)=Y3
12	CONTINUE
C
C******************************
14	YT(1,1)=0.0
	YT(2,1)=0.0
	YT(3,1)=0.0
	IF(IHNT2(3).LE.1) GO TO 24
	DO 20 KT=1,IHNT2(3)
15	R=HIRND(2)
c
         if(IHNT2(3).EQ.2) then
            rnd1=max(RLU_HIJING(NSEED),1.0e-20)
            rnd2=max(RLU_HIJING(NSEED),1.0e-20)
            rnd3=max(RLU_HIJING(NSEED),1.0e-20)
            R=-0.5*(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &           +4.38*0.85*log(rnd3)/(4.38+0.85))
         endif
c
	X=RLU_HIJING(0)
	CX=2.0*X-1.0
	SX=SQRT(1.0-CX*CX)
C		********choose theta from uniform cos(theta) distr
	PHI=RLU_HIJING(0)*2.0*HIPR1(40)
C		********chose phi form uniform phi distr 0 to 2*pi
	YT(1,KT)=R*SX*COS(PHI)
	YT(2,KT)=R*SX*SIN(PHI)
	YT(3,KT)=R*CX
	IF(HIPR1(29).EQ.0.0) GO TO 20
	DO 18  KT2=1,KT-1
		DNBT1=(YT(1,KT)-YT(1,KT2))**2
		DNBT2=(YT(2,KT)-YT(2,KT2))**2
		DNBT3=(YT(3,KT)-YT(3,KT2))**2
		DNBT=DNBT1+DNBT2+DNBT3
		IF(DNBT.LT.HIPR1(29)*HIPR1(29)) GO TO 15
C			********two neighbors cannot be closer than 
C				HIPR1(29)
18	CONTINUE
20	CONTINUE
c**********************************
         if(IHNT2(3).EQ.2) then
            YT(1,2)=-YT(1,1)
            YT(2,2)=-YT(2,1)
            YT(3,2)=-YT(3,1)
         endif
c*********************************
	DO 22 I=1,IHNT2(3)-1
	DO 22 J=I+1,IHNT2(3)
	IF(YT(3,I).LT.YT(3,J)) GO TO 22
	Y1=YT(1,I)
	Y2=YT(2,I)
	Y3=YT(3,I)
	YT(1,I)=YT(1,J)
	YT(2,I)=YT(2,J)
	YT(3,I)=YT(3,J)
	YT(1,J)=Y1
	YT(2,J)=Y2
	YT(3,J)=Y3
22	CONTINUE

C********************
24	MISS=-1
C****	BB IS THE ABSOLUTE VALUE OF IMPACT PARAMETER,BB**2 IS 
C       RANDOMLY GENERATED AND ITS ORIENTATION IS RANDOMLY SET 
C       BY THE ANGLE PHI  FOR EACH COLLISION.******************
C
	BB=SQRT(BMIN**2+RLU_HIJING(0)*(BMAX**2-BMIN**2))
	PHI=2.0*HIPR1(40)*RLU_HIJING(0)

50	MISS=MISS+1
	IF(MISS.GT.50) THEN
	   WRITE(6,*) 'infinite loop happened in  HIJING'
	   STOP
	ENDIF

	NATT=0
	JATT=0
	EATT=0.0
	CALL HIJINI
        NLOP=0
C			********Initialize for a new event
60	NT=0
	NP=0
	N0=0
	N01=0
	N10=0
	N11=0
        NELT=0
        NINT=0
        NELP=0
        NINP=0
	NSG=0
	NCOLT=0
        NPSPECP=0
        NNSPECP=0
        NPSPECT=0
        NNSPECT=0 

	BBX=BB*COS(PHI)
	BBY=BB*SIN(PHI)
	HINT1(19)=BB
	HINT1(20)=PHI
C
	DO 70 JP=1,IHNT2(1)
	DO 70 JT=1,IHNT2(3)
	   SCIP(JP,JT)=-1.0
	   B2=(YP(1,JP)+BBX-YT(1,JT))**2+(YP(2,JP)+BBY-YT(2,JT))**2
	   R2=B2*HIPR1(40)/HIPR1(31)/0.1
	   SCIP2(JP,JT)=R2
C		********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
	   RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)
     &          /1.2**2/REAL(IHNT2(1))**0.6666667,1.0)
	   RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)
     &          /1.2**2/REAL(IHNT2(3))**0.6666667,1.0)
	   APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
	   APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
	   HINT1(18)=HINT1(14)-APHX1*HINT1(15)
     &			-APHX2*HINT1(16)+APHX1*APHX2*HINT1(17)
	   IF(IHPR2(14).EQ.0.OR.
     &          (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) THEN
	      GS=1.0-EXP(-(HIPR1(30)+HINT1(18))*ROMG(R2)/HIPR1(31))
	      RANTOT=RLU_HIJING(0)
	      IF(RANTOT.GT.GS) GO TO 70
	      GO TO 65
	   ENDIF
	   GSTOT_0=2.0*(1.0-EXP(-(HIPR1(30)+HINT1(18))
     &             /HIPR1(31)/2.0*ROMG(0.0)))
	   R2=R2/GSTOT_0
	   GS=1.0-EXP(-(HIPR1(30)+HINT1(18))/HIPR1(31)*ROMG(R2))
	   GSTOT=2.0*(1.0-SQRT(1.0-GS))
	   RANTOT=RLU_HIJING(0)*GSTOT_0
	   IF(RANTOT.GT.GSTOT) GO TO 70
	   IF(RANTOT.GT.GS) THEN
	      CALL HIJCSC(JP,JT)
	      GO TO 70
C			********perform elastic collisions
	   ENDIF
 65	   SCIP(JP,JT)=R2
	   RNIP(JP,JT)=RANTOT
	   SJIP(JP,JT)=HINT1(18)
	   NCOLT=NCOLT+1
c           if (R2.GT.2.D0) THEN
c              write (8,*) R2
c           ENDIF
	   IPCOL(NCOLT)=JP
	   ITCOL(NCOLT)=JT
70	CONTINUE

c *** cl glauber ***
        npart=0
        xmeana=0D0
        ymeana=0D0
        xmeanb=0D0
        ymeanb=0D0
        xmeanp=0D0
        ymeanp=0D0
        xm2=0D0
        ym2=0D0
        xym=0D0
c
	IF(NCOLT>0) THEN 
           DO 1110 JP=1,IHNT2(1)
              xmeana=xmeana+YP(1,JP)+BBX
              ymeana=ymeana+YP(2,JP)+BBY
              DO 1120 JT=1,IHNT2(3)
                 IF(SCIP2(JP,JT).LT.2.0D0) THEN
                    npart=npart+1
                    xmeanp=xmeanp+YP(1,JP)+BBX
                    ymeanp=ymeanp+YP(2,JP)+BBY
                    xm2=xm2+(YP(1,JP)+BBX)*(YP(1,JP)+BBX)
                    ym2=ym2+(YP(2,JP)+BBY)*(YP(2,JP)+BBY)
                    xym=xym+(YP(1,JP)+BBX)*(YP(2,JP)+BBY)
                    goto 1110
                 end if
 1120         continue
 1110      continue
c
           DO 1130 JT=1,IHNT2(3)
              xmeanb=xmeanb+YT(1,JT)
              ymeanb=ymeanb+YT(2,JT)
              DO 1140 JP=1,IHNT2(1)
                 IF(SCIP2(JP,JT).LT.2.0D0) THEN
                    npart=npart+1
                    xmeanp=xmeanp+YT(1,JT)
                    ymeanp=ymeanp+YT(2,JT)
                    xm2=xm2+YT(1,JT)*YT(1,JT)
                    ym2=ym2+YT(2,JT)*YT(2,JT)
                    xym=xym+YT(1,JT)*YT(2,JT)
                    goto 1130
                 end if
 1140         continue
 1130      continue

           IF (npart.GT.0) THEN
              xmeana=xmeana/IHNT2(1)
              ymeana=ymeana/IHNT2(1)
              xmeanb=xmeanb/IHNT2(3)
              ymeanb=ymeanb/IHNT2(3)
              xmeanp=xmeanp/npart
              ymeanp=ymeanp/npart
              xm2=xm2/npart
              ym2=ym2/npart
              xym=xym/npart
c
              sx2=xm2-xmeanp*xmeanp
              sy2=ym2-ymeanp*ymeanp
              sxy=xym-xmeanp*ymeanp
c           
              delx=xmeanb-xmeana
              dely=ymeanb-ymeana
              dtmp=delx**2+dely**2
              bbtrue=sqrt(dtmp)
              dnumt=(sy2-sx2)*(delx**2-dely**2)-4D0*sxy*delx*dely
              ddent=(sy2+sx2)*bbtrue**2
              eccrp=dnumt/ddent
              dtmp=(sy2-sx2)*(sy2-sx2)+4D0*sxy*sxy
              eccpart=sqrt(dtmp)/(sx2+sy2)
              eccmc=(sy2-sx2)/(sy2+sx2)
              write(*,*),'HOUT: ',bb,' ',bbtrue,' ',ncolt,' ',npart,
     1             ' ',eccrp,' ',eccpart, bbx, bby
           end if
        end if
c
C		********total number interactions proj and targ has
C				suffered

	IF(NCOLT.EQ.0) THEN
	   NLOP=NLOP+1
           IF(NLOP.LE.20.OR.
     &          (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) THEN 
              BB=SQRT(BMIN**2+RLU_HIJING(0)*(BMAX**2-BMIN**2))
              GO TO 60
           ENDIF

           RETURN
	ENDIF
C               ********At large impact parameter, there maybe no
C                       interaction at all. For NN collision
C                       repeat the event until interaction happens
C
	IF(IHPR2(3).NE.0) THEN
	   NHARD=1+INT(RLU_HIJING(0)*(NCOLT-1)+0.5)
	   NHARD=MIN(NHARD,NCOLT)
	   JPHARD=IPCOL(NHARD)
	   JTHARD=ITCOL(NHARD)
	ENDIF
C
	IF(IHPR2(9).EQ.1) THEN
		NMINI=1+INT(RLU_HIJING(0)*(NCOLT-1)+0.5)
		NMINI=MIN(NMINI,NCOLT)
		JPMINI=IPCOL(NMINI)
		JTMINI=ITCOL(NMINI)
	ENDIF
C		********Specifying the location of the hard and
C			minijet if they are enforced by user
C
	DO 200 JP=1,IHNT2(1)
	DO 200 JT=1,IHNT2(3)
	IF(SCIP(JP,JT).EQ.-1.0) GO TO 200
		NFP(JP,11)=NFP(JP,11)+1
		NFT(JT,11)=NFT(JT,11)+1
	IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).GT.1) THEN
		NP=NP+1
		N01=N01+1
	ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).LE.1) THEN
		NT=NT+1
		N10=N10+1
	ELSE IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).LE.1) THEN
		NP=NP+1
		NT=NT+1
		N0=N0+1
	ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).GT.1) THEN
		N11=N11+1
	ENDIF
	JOUT=0
	NFP(JP,10)=0
	NFT(JT,10)=0
C*****************************************************************
	IF(IHPR2(8).EQ.0 .AND. IHPR2(3).EQ.0) GO TO 160
C		********When IHPR2(8)=0 no jets are produced
	IF(NFP(JP,6).LT.0 .OR. NFT(JT,6).LT.0) GO TO 160
C		********jets can not be produced for (JP,JT)
C			because not enough energy avaible for 
C				JP or JT 
	R2=SCIP(JP,JT)
	HINT1(18)=SJIP(JP,JT)
	TT=ROMG(R2)*HINT1(18)/HIPR1(31)
	TTS=HIPR1(30)*ROMG(R2)/HIPR1(31)
	NJET=0
	IF(IHPR2(3).NE.0 .AND. JP.EQ.JPHARD .AND. JT.EQ.JTHARD) THEN
           CALL JETINI(JP,JT,1)
           CALL HIJHRD(JP,JT,0,JFLG,0)
           HINT1(26)=HINT1(47)
           HINT1(27)=HINT1(48)
           HINT1(28)=HINT1(49)
           HINT1(29)=HINT1(50)
           HINT1(36)=HINT1(67)
           HINT1(37)=HINT1(68)
           HINT1(38)=HINT1(69)
           HINT1(39)=HINT1(70)
C
	   IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
	   IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
	   IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &				JFLG.GE.3) IASG(NSG,3)=1
	   IHNT2(9)=IHNT2(14)
	   IHNT2(10)=IHNT2(15)
	   DO 105 I05=1,5
	      HINT1(20+I05)=HINT1(40+I05)
	      HINT1(30+I05)=HINT1(50+I05)
 105	   CONTINUE
	   JOUT=1
	   IF(IHPR2(8).EQ.0) GO TO 160
	   RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)/1.2**2
     &		/REAL(IHNT2(1))**0.6666667,1.0)
	   RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)/1.2**2
     &		/REAL(IHNT2(3))**0.6666667,1.0)
	   APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
	   APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
	   HINT1(65)=HINT1(61)-APHX1*HINT1(62)
     &			-APHX2*HINT1(63)+APHX1*APHX2*HINT1(64)
	   TTRIG=ROMG(R2)*HINT1(65)/HIPR1(31)
	   NJET=-1
C		********subtract the trigger jet from total number
C			of jet production  to be done since it has
C				already been produced here
	   XR1=-ALOG(EXP(-TTRIG)+RLU_HIJING(0)*(1.0-EXP(-TTRIG)))
 106	   NJET=NJET+1
	   XR1=XR1-ALOG(RLU_HIJING(0))
	   IF(XR1.LT.TTRIG) GO TO 106
	   XR=0.0
 107	   NJET=NJET+1
	   XR=XR-ALOG(RLU_HIJING(0))
	   IF(XR.LT.TT-TTRIG) GO TO 107
	   NJET=NJET-1
	   GO TO 112
	ENDIF
C		********create a hard interaction with specified P_T
c				 when IHPR2(3)>0
	IF(IHPR2(9).EQ.1.AND.JP.EQ.JPMINI.AND.JT.EQ.JTMINI) GO TO 110
C		********create at least one pair of mini jets 
C			when IHPR2(9)=1
C
	IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LT.EXP(-TT)*
     &		(1.0-EXP(-TTS))) GO TO 160
C		********this is the probability for no jet production
110	XR=-ALOG(EXP(-TT)+RLU_HIJING(0)*(1.0-EXP(-TT)))
111	NJET=NJET+1
	XR=XR-ALOG(RLU_HIJING(0))
	IF(XR.LT.TT) GO TO 111
112	NJET=MIN(NJET,IHPR2(8))
	IF(IHPR2(8).LT.0)  NJET=ABS(IHPR2(8))
C		******** Determine number of mini jet production
C
	DO 150 I_JET=1,NJET
           CALL JETINI(JP,JT,0)
	   CALL HIJHRD(JP,JT,JOUT,JFLG,1)
C		********JFLG=1 jets valence quarks, JFLG=2 with 
C			gluon jet, JFLG=3 with q-qbar prod for
C			(JP,JT). If JFLG=0 jets can not be produced 
C			this time. If JFLG=-1, error occured abandon
C			this event. JOUT is the total hard scat for
C			(JP,JT) up to now.
	   IF(JFLG.EQ.0) GO TO 160
	   IF(JFLG.LT.0) THEN
	      IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJHRD'
	      GO TO 50
	   ENDIF
	   JOUT=JOUT+1
	   IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
	   IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
	   IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &			JFLG.GE.3) IASG(NSG,3)=1
C		******** jet with PT>HIPR1(11) will be quenched
 150	CONTINUE
 160	CONTINUE
	CALL HIJSFT(JP,JT,JOUT,IERROR)
	IF(IERROR.NE.0) THEN
	   IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJSFT'
	   GO TO 50
	ENDIF
C
C		********conduct soft scattering between JP and JT
	JATT=JATT+JOUT

200	CONTINUE

c
c**************************
c
	DO 201 JP=1,IHNT2(1)
c           write(6,*) JP, NFP(JP,3), NFP(JP,4), NFP(JP,5)
           IF(NFP(JP,5).GT.2) THEN
              NINP=NINP+1
           ELSE IF(NFP(JP,5).EQ.2.OR.NFP(JP,5).EQ.1) THEN
              NELP=NELP+1
           ENDIF

           IF(NFP(JP,5).LE.2) THEN
              IF (NFP(JP,3) .EQ. 2212) THEN
                 NPSPECP = NPSPECP + 1
              ELSE IF (NFP(JP,3) .EQ. 2112) THEN
                 NNSPECP = NNSPECP + 1
              ENDIF
           ENDIF
 201    continue
	DO 202 JT=1,IHNT2(3)
           IF(NFT(JT,5).GT.2) THEN
              NINT=NINT+1
           ELSE IF(NFT(JT,5).EQ.2.OR.NFT(JT,5).EQ.1) THEN
              NELT=NELT+1
           ENDIF

           IF(NFT(JT,5).LE.2) THEN
              IF (NFT(JT,3) .EQ. 2212) THEN
                 NPSPECT = NPSPECT + 1
              ELSE IF (NFT(JT,3) .EQ. 2112) THEN
                 NNSPECT = NNSPECT + 1
              ENDIF
           ENDIF
 202    continue
c     
c*******************************

C********perform jet quenching for jets with PT>HIPR1(11)**********

	IF((IHPR2(8).NE.0.OR.IHPR2(3).NE.0).AND.IHPR2(4).GT.0.AND.
     &			IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
		DO 271 I=1,IHNT2(1)
			IF(NFP(I,7).EQ.1) CALL QUENCH(I,1)
271		CONTINUE
		DO 272 I=1,IHNT2(3)
			IF(NFT(I,7).EQ.1) CALL QUENCH(I,2)
272		CONTINUE
		DO 273 ISG=1,NSG
			IF(IASG(ISG,3).EQ.1) CALL QUENCH(ISG,3)
273		CONTINUE
	ENDIF
C
C**************fragment all the string systems in the following*****
C
C********N_ST is where particle information starts
C********N_STR+1 is the number of strings in fragmentation
C********the number of strings before a line is stored in K(I,4)
C********IDSTR is id number of the string system (91,92 or 93)
C
c
        IF(IHPR2(20).NE.0) THEN
	   DO 360 ISG=1,NSG
		CALL HIJFRG(ISG,3,IERROR)
        	IF(MSTU(24).NE.0 .OR.IERROR.GT.0) THEN
		   MSTU(24)=0
		   MSTU(28)=0
		   IF(IHPR2(10).NE.0) THEN
		      call LULIST_HIJING(1)
		      WRITE(6,*) 'error occured, repeat the event'
		   ENDIF
		   GO TO 50
		ENDIF
C			********Check errors
C
		N_ST=1
		IDSTR=92
		IF(IHPR2(21).EQ.0) THEN
		   CALL LUEDIT_HIJING(2)
		ELSE
351                N_ST=N_ST+1
		   IF((N_ST .LE. 9000)
     &            .AND. (K(N_ST,2).LT.91.OR.K(N_ST,2).GT.93)) GO TO  351
                   IDSTR=K(N_ST,2)
                   N_ST=N_ST+1
		ENDIF
C
		IF(FRAME.EQ.'LAB') THEN
			CALL HIBOOST
		ENDIF
C		******** boost back to lab frame(if it was in)
C
		N_STR=0
		DO 360 I=N_ST,N
		   IF(K(I,2).EQ.IDSTR) THEN
		      N_STR=N_STR+1
		      GO TO 360
		   ENDIF
		   K(I,4)=N_STR
		   NATT=NATT+1
		   KATT(NATT,1)=K(I,2)
		   KATT(NATT,2)=20
		   KATT(NATT,4)=K(I,1)
                   IF(K(I,3).EQ.0 ) THEN
                      KATT(NATT,3)=0
                   ELSE
                      IF(K(K(I,3),2).EQ.IDSTR) THEN
                         KATT(NATT,3)=0
                      ELSE
                         KATT(NATT,3)=NATT-I+K(I,3)+N_STR-K(K(I,3),4)
                      ENDIF
                   ENDIF
C       ****** identify the mother particle
		   PATT(NATT,1)=P(I,1)
		   PATT(NATT,2)=P(I,2)
		   PATT(NATT,3)=P(I,3)
		   PATT(NATT,4)=P(I,4)
                   VATT(NATT,1)=V(I,1)
                   VATT(NATT,2)=V(I,2)
                   VATT(NATT,3)=V(I,3)
                   VATT(NATT,4)=V(I,4)

		   EATT=EATT+P(I,4)
360	   CONTINUE
C		********Fragment the q-qbar jets systems *****
C
	   JTP(1)=IHNT2(1)
	   JTP(2)=IHNT2(3)
	   DO 400 NTP=1,2
	   DO 400 J_JTP=1,JTP(NTP)
		CALL HIJFRG(J_JTP,NTP,IERROR)
        	IF(MSTU(24).NE.0 .OR. IERROR.GT.0) THEN
		   MSTU(24)=0
		   MSTU(28)=0
		   IF(IHPR2(10).NE.0) THEN
		      call LULIST_HIJING(1)
		      WRITE(6,*) 'error occured, repeat the event'
		   ENDIF
		   GO TO 50
		ENDIF
C			********check errors
C
		N_ST=1
		IDSTR=92

		NFTP=NFP(J_JTP,5)
		IF(NTP.EQ.2) NFTP=10+NFT(J_JTP,5)

		IF(IHPR2(21).EQ.0) THEN
		   CALL LUEDIT_HIJING(2)
		ELSE IF (NFTP.EQ. 3 .OR. NFTP .EQ. 13) THEN
381		   N_ST=N_ST+1
		   IF(K(N_ST,2).LT.91.OR.K(N_ST,2).GT.93) GO TO  381
		   IDSTR=K(N_ST,2)
		   N_ST=N_ST+1
                ENDIF
		IF(FRAME.EQ.'LAB') THEN
			CALL HIBOOST
		ENDIF
C		******** boost back to lab frame(if it was in)
C

		N_STR=0
		DO 390 I=N_ST,N
		   IF(K(I,2).EQ.IDSTR) THEN
		      N_STR=N_STR+1
		      GO TO 390
		   ENDIF
		   K(I,4)=N_STR
		   NATT=NATT+1
		   KATT(NATT,1)=K(I,2)
		   KATT(NATT,2)=NFTP
		   KATT(NATT,4)=K(I,1)
                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
		      KATT(NATT,3)=0
		   ELSE
		      KATT(NATT,3)=NATT-I+K(I,3)+N_STR-K(K(I,3),4)
		   ENDIF
C       ****** identify the mother particle
		   PATT(NATT,1)=P(I,1)
		   PATT(NATT,2)=P(I,2)
		   PATT(NATT,3)=P(I,3)
		   PATT(NATT,4)=P(I,4)
		   EATT=EATT+P(I,4)
                   VATT(NATT,1)=V(I,1)  
                   VATT(NATT,2)=V(I,2)  
                   VATT(NATT,3)=V(I,3)  
                   VATT(NATT,4)=V(I,4)
390		CONTINUE 
400	   CONTINUE
C		********Fragment the q-qq related string systems
	ENDIF

	DO 450 I=1,NDR
		NATT=NATT+1
		KATT(NATT,1)=KFDR(I)
		KATT(NATT,2)=40
		KATT(NATT,3)=0
		PATT(NATT,1)=PDR(I,1)
		PATT(NATT,2)=PDR(I,2)
		PATT(NATT,3)=PDR(I,3)
		PATT(NATT,4)=PDR(I,4)
                VATT(NATT,1)=V(I,1)  
                VATT(NATT,2)=V(I,2)  
                VATT(NATT,3)=V(I,3)  
                VATT(NATT,4)=V(I,4)
		EATT=EATT+PDR(I,4)
450	CONTINUE
C			********store the direct-produced particles
C
	DENGY=EATT/(IHNT2(1)*HINT1(6)+IHNT2(3)*HINT1(7))-1.0
	IF(ABS(DENGY).GT.HIPR1(43).AND.IHPR2(20).NE.0
     &     .AND.IHPR2(21).EQ.0) THEN
	IF(IHPR2(10).NE.0) WRITE(6,*) 'Energy not conserved, repeat the
     &	   event'
C		call LULIST_HIJING(1)
		GO TO 50
	ENDIF

	RETURN
	END
