!********************************************************************** 
!                                                                       
!     TABLE OF CONTENTS                                                 
!                                                                       
!     1. evolinit -> Initializes the evolution code by filling the COMMO
!                    BLOCKS which contain the constants                 
!     2. alphasNNPDF,ffn,vfn,funcs -> Routines related to the computatio
!     3. evolfactx -> Computes the evolution kernels in x-space.        
!     4. xgrid  -> Gives the x pts for the convolution                  
!     5. gammax -> FT algorithm for Mellin inversion                    
!     6. gamma -> evaluates gamma coefficients                          
!     7. evolfactn -> N space quantities evaluation                     
!     8. evolint -> convolution routines                                
!     9. andim_IPT -> Evaluates an dimension                            
!     10.matutils.f, psi.f -> Matrix operation, special functions and dg
!                                                                       
!     NB: pdfin & NN are in wrapXNN.f                                   
!                                                                       
!********************************************************************** 
                                                                        
!****                                                                   
! 1 *-------------------------------------------------------------------
!****                                                                   
      SUBROUTINE LH_EVOLINIT 
      IMPLICIT none 
!                                                                       
      INTEGER ipt,imodev,ivfn,itmc 
      COMMON/nnpdf10EVFLAGS/ipt,imodev,ivfn,itmc 
!                                                                       
      REAL*8 ca,cf,tr 
      COMMON/NNPDF10COLFACT/ca,cf,tr 
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON /NNPDF10BETA/ beta0,beta1,beta2,b1,b2 
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON/nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      REAL*8 q20,q2 
      COMMON/NNPDF10EVSCALE/q20,q2 
!                                                                       
      INTEGER i, nff 
      REAL*8 alphas,LH_AlphaExactBeta,LH_AlphaMZ 
      EXTERNAL LH_AlphaExactBeta,LH_AlphaMZ 
!                                                                       
!     Color Factors                                                     
!                                                                       
      CA= 3D0 
      CF= 4D0/3D0 
      TR= 0.5D0 
!                                                                       
!     Constants                                                         
!                                                                       
      EMC   = 0.5772156649D0 
      ZETA2 = 1.644934067D0 
      ZETA3 = 1.2020569031D0 
      ZETA4 = 1.0823232337D0 
!                                                                       
!     Beta function coefficients                                        
!                                                                       
      DO I=3,6 
         BETA0(I) = (33d0-2d0*I)/3d0 
         BETA1(I) = 102d0-38d0/3d0*I 
         BETA2(I) = 2857d0/2d0-5033d0/18d0*I+325d0/54d0*I**2d0 
         B1(I) = BETA1(I)/BETA0(I) 
         B2(I) = BETA2(I)/BETA0(I) 
      ENDDO 
!                                                                       
      IF(IVFN.EQ.0) THEN 
         NFF = 3
         IF(IMODEV.EQ.1)THEN 
            CALL LH_FFN(q20,alphas,ipt,LH_AlphaExactBeta,NFF) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_FFN(q20,alphas,ipt,LH_AlphaMZ,NFF) 
         ENDIF 
         AS0 = ALPHAS/4d0/PI 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_FFN(q2TH(4),alphas,ipt,LH_AlphaExactBeta,NFF) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_FFN(q2TH(4),alphas,ipt,LH_AlphaMZ,NFF) 
         ENDIF 
         ASC = ALPHAS/4d0/PI 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_FFN(q2TH(5),alphas,ipt,LH_AlphaExactBeta,NFF) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_FFN(q2TH(5),alphas,ipt,LH_AlphaMZ,NFF) 
         ENDIF 
         ASB = ALPHAS/4d0/PI 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_FFN(q2TH(6),alphas,ipt,LH_AlphaExactBeta,NFF) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_FFN(q2TH(6),alphas,ipt,LH_AlphaMZ,NFF) 
         ENDIF 
         AST = ALPHAS/4d0/PI 
      ELSE 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_VFN(q20,alphas,ipt,LH_AlphaExactBeta) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_VFN(q20,alphas,ipt,LH_AlphaMZ) 
         ENDIF 
         AS0 = ALPHAS/4d0/PI 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_VFN(q2TH(4),alphas,ipt,LH_AlphaExactBeta) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_VFN(q2TH(4),alphas,ipt,LH_AlphaMZ) 
         ENDIF 
         ASC = ALPHAS/4d0/PI 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_VFN(q2TH(5),alphas,ipt,LH_AlphaExactBeta) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_VFN(q2TH(5),alphas,ipt,LH_AlphaMZ) 
         ENDIF 
         ASB = ALPHAS/4d0/PI 
         IF(IMODEV.EQ.1)THEN 
            CALL LH_VFN(q2TH(6),alphas,ipt,LH_AlphaExactBeta) 
         ELSEIF(IMODEV.EQ.0)THEN 
            CALL LH_VFN(q2TH(6),alphas,ipt,LH_AlphaMZ) 
         ENDIF 
         AST = ALPHAS/4d0/PI 
      ENDIF 
!                                                                       
      ASQ = 0.d0 
!                                                                       
      RETURN 
      END                                           
                                                                        
!****                                                                   
! 2 *-------------------------------------------------------------------
!****                                                                   
      FUNCTION alphaNNPDF(QQ) 
      IMPLICIT none 
!                                                                       
      INTEGER nff 
      INTEGER ipt,imodev,ivfn,itmc 
      COMMON/NNPDF10EVFLAGS/ipt,imodev,ivfn,itmc 
!                                                                       
      REAL*8 alphaNNPDF 
      REAL*8 qq,qq2 
      REAL*8 alphas,LH_AlphaExactBeta,LH_AlphaMZ 
      EXTERNAL LH_AlphaExactBeta,LH_AlphaMZ 
!                                                                       
      qq2 = qq**2. 
      CALL lh_evolinit 
!                                                                       
      IF(ivfn.eq.0)THEN
         WRITE(*,*)"FFN not available!"
         NFF = 3
         IF(imodev.eq.1)THEN 
            CALL LH_FFN(qq2,alphas,ipt,LH_AlphaExactBeta,NFF) 
         ELSEIF(imodev.eq.0)THEN 
            CALL LH_FFN(qq2,alphas,ipt,LH_AlphaMZ,NFF) 
         ENDIF 
      ELSEIF(ivfn.EQ.1)THEN 
         IF (imodev.EQ.1)THEN 
            CALL LH_VFN(qq2,alphas,ipt,LH_AlphaExactBeta) 
         ELSEIF(imodev.EQ.0)THEN 
            CALL LH_VFN(qq2,alphas,ipt,LH_AlphaMZ) 
         ENDIF 
      ENDIF 
!                                                                       
      alphaNNPDF = alphas 
!                                                                       
      RETURN 
      END                                           
                                                                        
!                                                                       
!     FFN: returns the value of alpha_s computed for Fixed              
!          Flavour Number.                                              
!     VFN: returns the value of alpha_s computed for Variable           
!          Flavour Number.                                              
!                                                                       
      SUBROUTINE LH_FFN(q2,alphas,ipt,FUNC,nf) 
      IMPLICIT none 
!                                                                       
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON/nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      INTEGER nf 
      INTEGER ipt 
      REAL*8  q2,qq2 
      REAL*8  alphas,func,alphasref,qq2ref 
!                                                                       
      alphasref = asref 
      qq2ref    = q2ref 
      qq2       = q2 
                                                                        
      alphas = FUNC(nf,qq2ref,alphasref,qq2,ipt) 
                                                                        
      RETURN 
      END                                           
!                                                                       
!                                                                       
      SUBROUTINE LH_VFN(qq2,alphas,ipt,FUNC) 
      IMPLICIT none 
!
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 

      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON/nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!
      INTEGER nfi,nff,dnf,snf,ipt
      double precision q2,qq2,FUNC,alphasref,qq2ref
      double precision alphas,c2,asi
!     c2 is the same coefficient used in eq. 2.42
!     of hep-ph/0408244 and obtained in eq. 10 
!     of hep-ph/9706430. In the following it is divided
!     by (4*pi)^2 to match the notations
      parameter(c2=14d0/3d0)
      external FUNC

      q2=qq2

      if(q2.ge.q2th(6))then
         nff=6
      elseif(q2.ge.q2th(5))then
         nff=5
      elseif(q2.ge.q2th(4))then
         nff=4
      else
         nff=3
      endif
      if(q2ref.gt.q2th(6))then
         nfi=6
      elseif(q2ref.gt.q2th(5))then
         nfi=5
      elseif(q2ref.gt.q2th(4))then
         nfi=4
      else
         nfi=3
      endif
!
      alphasref = asref
      qq2ref = q2ref

 10   if(nff.eq.nfi) then
         alphas = FUNC(nfi,qq2ref,alphasref,q2,ipt)
         return
      else
         if(nff.gt.nfi)then
            dnf=1
            snf=1
         else
            dnf=-1
            snf=0
         endif
         asi = FUNC(nfi,qq2ref,alphasref,q2th(nfi+snf),ipt)
         if(ipt.ge.2)then
            if(nff.gt.nfi) asi=asi+(c2/(4d0*pi)**2d0)*asi**3d0
            if(nff.lt.nfi) asi=asi-(c2/(4d0*pi)**2d0)*asi**3d0
         endif
         alphasref = asi
         qq2ref = q2th(nfi+snf)
         nfi = nfi+dnf
         goto 10
      endif   
      end
!                                                                       
!     FUNCS for the computation of alpha_s with diff methods            
!     - AlphaMZ: computes alpha_s as function of alpha_s at a given     
!                refernce scale.                                        
!     - AlphaExactBeta: Exact solution of the QCD beta function         
!                       equation using fourth order Runge-Kutta         
!                       algorithm.                                      
!                                                                       
      FUNCTION LH_AlphaMZ(nfi,mz2,asz,q2,ipt) 
      IMPLICIT none 
!                                                                       
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON/nnpdf10BETA/beta0,beta1,beta2,b1,b2 
!                                                                       
      INTEGER nfi,ipt 
      REAL*8 q2ref,q2,asi,asref 
      REAL*8 alo,t,as,den,mz2,asz,LH_AlphaMZ 
                                                                        
      q2ref = mz2 
      asref = asz/4./pi 
      asi = asref 
      t   = log(q2/q2ref) 
      den = 1.+beta0(nfi)*asi*t 
      alo = asi/den 
                                                                        
!     LO                                                                
      as = alo 
                                                                        
!     NLO                                                               
      IF(ipt.GE.1) as = alo*(1-b1(nfi)*alo*log(den)) 
                                                                        
!     NNLO                                                              
      IF(ipt.GE.2)THEN 
         as = alo*(1.+(alo*(alo-asi)*(b2(nfi)-b1(nfi)**2.)              &
     &        +as*b1(nfi)*dlog(as/asi)))                                
      ENDIF 
!                                                                       
      LH_AlphaMZ = 4d0*pi*as 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
      FUNCTION LH_AlphaExactBeta(nf,r20,as0,r2,ipt) 
!                                                                       
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON/nnpdf10BETA/beta0,beta1,beta2,b1,b2 
!                                                                       
      INTEGER NFMIN, NFMAX, NF, NSTEP, K1,ipt,nnf 
      REAL*8 LH_fbeta,LH_AlphaExactBeta 
      REAL*8 as,as0,r20,r2 
      REAL*8 dlr,lrrat,sxth 
      REAL*8 xk0,xk1,xk2,xk3 
                                                                        
      PARAMETER (NFMIN = 3, NFMAX = 6) 
      PARAMETER (NSTEP=50) 
      PARAMETER (SXTH = 0.166666666666666D0 ) 
!                                                                       
      NNF = NF 
      AS0 = AS0/4./pi 
      AS  = AS0 
      LRRAT = dLOG (R2/R20) 
      DLR   = LRRAT / NSTEP 
!                                                                       
!     ..Solution of the evolution equation depending on  NAORD          
!   (fourth-order Runge-Kutta beyond the leading order)                 
      IF (IPT.EQ.0) THEN 
         AS = AS0 / (1d0+ BETA0(NF) * AS0 * LRRAT) 
      ELSE IF (IPT.EQ.1) THEN 
         DO 2 K1 = 1, NSTEP 
            XK0 = DLR * LH_FBETA (AS,NNF,IPT) 
            XK1 = DLR * LH_FBETA (AS + 0.5d0 * XK0,NNF,IPT) 
            XK2 = DLR * LH_FBETA (AS + 0.5d0 * XK1,NNF,IPT) 
            XK3 = DLR * LH_FBETA (AS + XK2,NNF,IPT) 
            AS = AS + SXTH * (XK0 + 2d0* XK1 + 2d0* XK2 + XK3) 
    2     CONTINUE 
      ELSE IF (IPT.EQ.2) THEN 
         DO 3 K1 = 1, NSTEP 
            XK0 = DLR * LH_FBETA (AS,nnf,ipt) 
            XK1 = DLR * LH_FBETA (AS + 0.5d0 * XK0,nnf,ipt) 
            XK2 = DLR * LH_FBETA (AS + 0.5d0 * XK1,nnf,ipt) 
            XK3 = DLR * LH_FBETA (AS + XK2,nnf,ipt) 
            AS = AS + SXTH * (XK0 + 2d0* XK1 + 2d0* XK2 + XK3) 
    3    CONTINUE 
      END IF 
                                                                        
      LH_alphaexactbeta = as*(4d0*pi) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_fbeta(a,nf,ipt) 
      IMPLICIT none 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON/nnpdf10BETA/beta0,beta1,beta2,b1,b2 
!                                                                       
      REAL*8 lh_fbeta,a 
      INTEGER nf,ipt 
                                                                        
      IF(ipt.EQ.1)THEN 
         LH_FBETA = - A**2d0 * ( BETA0(NF) + A * BETA1(NF) ) 
      ELSEIF(ipt.EQ.2)THEN 
         LH_FBETA = - A**2d0 * ( BETA0(NF) + A * ( BETA1(NF)            &
     &     + A * BETA2(NF) ) )                                          
      ENDIF 
                                                                        
      RETURN 
      END                                           
!                                                                       
!     Find \LambdaQCD(nf) from the reference scale                      
!     Needed for LHAPDF interface                                       
!                                                                       
      function LH_LAMBDAZ(NF) 
      implicit none 
!                                                                       
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON/nnpdf10BETA/beta0,beta1,beta2,b1,b2 
!                                                                       
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON/nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      INTEGER nf,i 
      REAL*8 LH_lambdaz 
      REAL*8 aref,lambda5sq,lambdaz5,qth(4:6) 
!                                                                       
      AREF = ASREF /(4.d0*PI) 
      DO I= 3,6 
         QTH(i) = dsqrt(Q2TH(I)) 
      ENDDO 
!                                                                       
      LAMBDA5SQ = Q2REF * DEXP(-1d0/(BETA0(5)*AREF)                     &
     &     - (BETA1(5)/BETA0(5)**2d0)                                   &
     &     * DLOG(BETA1(5)*AREF/(BETA0(5)+BETA1(5)*AREF)))              
      LAMBDAZ5 = DSQRT(LAMBDA5SQ) 
      IF(NF.EQ.5)THEN 
         LH_LAMBDAZ = LAMBDAZ5 
      ELSEIF(NF.EQ.6)THEN 
                               !!!!!! add formula!!!!                   
         LH_LAMBDAZ = LAMBDAZ5 
      ELSEIF(NF.LE.4)THEN 
         LH_LAMBDAZ = LAMBDAZ5 * ((QTH(5)/LAMBDAZ5)**(2d0/25d0))        &
     &        * dLOG(Q2TH(5)/LAMBDA5SQ)**(963d0/14375d0)                
      ENDIF 
                                                                        
!                                                                       
      RETURN 
      END                                           
!****                                                                   
! 3 *-------------------------------------------------------------------
!****                                                                   
!                                                                       
!     evolfactx: Computes the ev. kernels in x-space.(DIFF FROM CODE!! N
!                                                                       
      SUBROUTINE LH_EVOLFACTX(obs,x) 
      IMPLICIT none 
!                                                                       
      INTEGER npoints 
      PARAMETER(npoints=2**8) 
      INTEGER IPT,IMODEV,IVFN,ITMC 
      COMMON /NNPDF10EVFLAGS/ IPT,IMODEV,IVFN,ITMC 
      INTEGER ieval,niter,nmax 
      COMMON /nnpdf10NGRID/ ieval,niter,nmax 
      REAL*8 xmin,xm1,xm2,xmax 
      COMMON/nnpdf10GRID/xmin,xm1,xm2,xmax 
!                                                                       
                                             ! NO MORE FNCTS OF NDATA   
      REAL*8 gm2ns(4),gm2ns24(2),gm2sg(2,2) 
      COMMON /nnpdf10GM2/gm2ns,gm2ns24,gm2sg 
      REAL*8 xx(NPOINTS),wn(NPOINTS),evkns(4,NPOINTS),                  &
     &     evksg(2,2,NPOINTS),evkns24(2,NPOINTS)                        
      COMMON /nnpdf10XEVK/ xx,wn,evkns,evksg,evkns24 
!                                                                       
      REAL*8 q20,q2 
      COMMON/nnpdf10EVSCALE/q20,q2 
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON/nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      INTEGER j,obs 
      REAL*8 x 
      REAL*8 flx 
      REAL*8 LH_gamm2ns,LH_gamm2sg,LH_gamm2ns24 
      REAL*8 xxtmp(NPOINTS),wtmp(NPOINTS) 
      REAL*8 LH_gammax_nsp,LH_gammax_nsm,LH_gammax_nsvu,LH_gammax_nsvd 
      REAL*8 LH_gammax_ns24q,LH_gammax_ns24g 
      REAL*8 LH_gammax_qq, LH_gammax_qg, LH_gammax_gq, LH_gammax_gg 
!                                                                       
      flx=0.2d0 
      if(ipt.eq.2)flx=0.3d0 
!                                                                       
      do j=1,NPOINTS 
         xxtmp(j)=0.d0 
         wtmp(j)=0.d0 
      enddo 
!                                                                       
      call LH_xgrid(x,xm1,xm2,xmax,flx,niter,xxtmp,wtmp,nmax) 
!                                                                       
      DO j=1,NMAX 
         xx(j)          = xxtmp(j) 
         wn(j)          = wtmp(j) 
         evkns(1,j)     = LH_gammax_nsp(xx(j),Q20,Q2) 
         evkns(2,j)     = LH_gammax_nsm(xx(j),Q20,Q2) 
         evkns(3,j)     = LH_gammax_nsvu(xx(j),Q20,Q2) 
         evkns(4,j)     = LH_gammax_nsvd(xx(j),Q20,Q2) 
         evkns24(1,j)   = LH_gammax_ns24q(xx(j),Q20,Q2) 
         evkns24(2,j)   = LH_gammax_ns24g(xx(j),Q20,Q2) 
         evksg(1,1,j)   = LH_gammax_qq(xx(j),Q20,Q2) 
         evksg(1,2,j)   = LH_gammax_qg(xx(j),Q20,Q2) 
         evksg(2,1,j)   = LH_gammax_gq(xx(j),Q20,Q2) 
         evksg(2,2,j)   = LH_gammax_gg(xx(j),Q20,Q2) 
      ENDDO 
!                                                                       
      gm2ns(1)     = LH_gamm2ns(1,x,Q20,Q2) 
      gm2ns(2)     = LH_gamm2ns(2,x,Q20,Q2) 
      gm2ns(3)     = LH_gamm2ns(3,x,Q20,Q2) 
      gm2ns(4)     = LH_gamm2ns(4,x,Q20,Q2) 
      gm2ns24(1)   = LH_gamm2ns24(1,x,Q20,Q2) 
      gm2ns24(2)   = LH_gamm2ns24(2,x,Q20,Q2) 
      gm2sg(1,1)   = LH_gamm2sg(1,1,x,Q20,Q2) 
      gm2sg(1,2)   = LH_gamm2sg(1,2,x,Q20,Q2) 
      gm2sg(2,1)   = LH_gamm2sg(2,1,x,Q20,Q2) 
      gm2sg(2,2)   = LH_gamm2sg(2,2,x,Q20,Q2) 
!                                                                       
      RETURN 
      END                                           
!****                                                                   
! 4 *------------------------------------------------------------------ 
!****                                                                   
!                                                                       
!     Gives the pts in x for the convolution grid                       
!                                                                       
                                                                        
      SUBROUTINE LH_xgrid(xmin,xmm1,xmm2,xmax,flx,NITER,xx,wn,NMAX) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER npoints 
      PARAMETER(npoints=2**8) 
      INTEGER NINT,NINT1,NINT2,NINT3 
      INTEGER ntmp,NITER,NGAUSS,NMAX 
      PARAMETER(NGAUSS=4) 
      INTEGER i,k,ii,k2,k3 
      REAL*8 xmin,xmax,xm1,rtmp,xdum,xm2,xmm2,flx,xmm1 
      REAL*8 sum,diff,z 
      REAL*8 xx(NPOINTS),wn(NPOINTS) 
      REAL*8 down(NPOINTS),up(NPOINTS) 
      REAL*8 ZS(NGAUSS),WZ(NGAUSS) 
!                                                                       
      DATA WZ                                                           &
     &       / 0.347854845137454D0,				      &
     &         0.652145154862546D0,				      &
     &         0.652145154862546D0,				      &
     &         0.347854845137454D0/				      
                                                                        
      DATA ZS                                                           &
     &       / -0.861136311594053D0,				      &
     &         -0.339981043584856D0,				      &
     &          0.339981043584856D0,				      &
     &          0.861136311594053D0 /				      
                                                                        
      k=0 
                                                                        
      xm1=xmm1 
                                                                        
      DO ii=0,NITER 
                                                                        
!     X < 0.99: 3 kind of division                                      
                                                                        
         IF(xmin.lt.xmm2)THEN 
                                                                        
            NINT=2**ii 
            NINT3=int(flx*NINT) 
            xm2=xmm2 
            if(NINT3.lt.1)xm2=xmax 
                                                                        
            if(xmin.lt.xm1)then 
               rtmp=abs((dlog(xm1-xmin))/(xm2-xmin)) 
               ntmp=int(NINT/rtmp) 
               NINT1=ntmp 
            else 
               NINT1=0 
            endif 
                                                                        
            NINT2=NINT-NINT1-NINT3 
                                                                        
            if(xmin.le.xm1)then 
               do i=1,NINT1 
                  xdum=xmin*(xm1/xmin)**(dble(i-1)/dble(NINT1)) 
                  down(i)=xdum 
               enddo 
            endif 
                                                                        
            if(xmin.gt.xm1)xm1=xmin 
            if(xmin.le.xm1)then 
               do i=1,NINT2 
                  xdum=xm1+(xm2-xm1)*(dble(i-1)/dble(NINT2)) 
                  down(i+NINT1)=xdum 
               enddo 
            endif 
                                                                        
            do i=1,NINT3 
               xdum=1.d0-(1.d0-xm2)*((1.d0-xmax)/                       &
     &              (1.d0-xm2))**(dble(i-1)/dble(NINT3))                
               down(i+NINT1+NINT2)=xdum 
            enddo 
                                                                        
!     X>= 0.99: only intervals logarithmic in log(1-x)                  
         ELSE 
!                                                                       
            NINT = 2**ii 
                                                                        
            DO i=1,NINT 
               xdum=1.d0-(1.d0-xmin)*((1.d0-xmax)/                      &
     &              (1.d0-xmin))**(dble(i-1)/dble(NINT))                
               down(i)=xdum 
            ENDDO 
                                                                        
         ENDIF 
                                                                        
         do i = 2, NINT 
            up(i-1) = down(i) 
         enddo 
         up(NINT)  = xmax 
                                                                        
         DO K2 = 1, NINT 
            SUM  = UP(K2) + DOWN(K2) 
            DIFF = UP(K2) - DOWN(K2) 
            DO K3 = 1, NGAUSS 
               K = K + 1 
               Z = (SUM + DIFF * ZS(K3)) * 0.5d0 
               WN(K) = DIFF * 0.5d0 * WZ(K3) 
               xx(K) = Z 
               IF(Z.LT.XMIN)WRITE(*,*)"WARNING",K,XX(K),k2,DOWN(K2),UP(K2)                                 
            ENDDO 
         ENDDO 
                                                                        
!                                                                       
      ENDDO 
!                                                                       
      NMAX=k 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!****                                                                   
! 5 *-------------------------------------------------------------------
!****                                                                   
!                                                                       
!     File: gammax.f                                                    
!     Contains the functions to compute x*Gamma(x,Q2i,Q2f) as           
!     Mellin inversions of the evolution operators computed             
!     in N-space.                                                       
!                                                                       
      FUNCTION LH_gammax_nsp(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_nsp 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_nsp(evpdf,x,Q2i,Q2f) 
      LH_gammax_nsp = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_nsm(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_nsm 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_nsm(evpdf,x,Q2i,Q2f) 
      LH_gammax_nsm = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_nsvu(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_nsvu 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_nsvu(evpdf,x,Q2i,Q2f) 
      LH_gammax_nsvu = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_nsvd(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_nsvd 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_nsvd(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_nsvd = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_ns24q(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_ns24q 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_ns24q(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_ns24q = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_ns24g(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_ns24g 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_ns24g(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_ns24g = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_qq(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_qq 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_qq(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_qq = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_qg(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_qg 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_qg(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_qg = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_gq(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_gq 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_gq(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_gq = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_gg(x,Q2i,Q2f) 
      IMPLICIT none 
      REAL*8 LH_gammax_gg 
      REAL*8 x,Q2i,Q2f 
      REAL*8 evpdf 
                                                                        
      call LH_fixedtalbot_gg(evpdf,x,Q2i,Q2f) 
                                                                        
      LH_gammax_gg = x*evpdf 
                                                                        
      RETURN 
      END                                           
!                                                                       
!     Fixed Talbot algorithm                                            
                                                                        
      SUBROUTINE LH_fixedtalbot_nsp(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m = 16 
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_ZFUNC(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(1) 
         tmp=tmp+dreal(tmp2) 
      enddo 
                                                                        
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(1))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_nsm(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m = 16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(2) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(2))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_nsvu(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m = 16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(3) 
         tmp=tmp+dreal(tmp2) 
      enddo 
                                                                        
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(3))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_nsvd(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m = 16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(4) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(4))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_ns24q(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m=16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_ZFUNC(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns24(1) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns24(1))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_ns24g(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m=16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_ZFUNC(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns24(2) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns24(2))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_qq(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m=16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(1,1) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(1,1))) 
                                                                        
      RETURN 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_qg(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m=16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(1,2) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(1,2))) 
                                                                        
      return 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_gq(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m=16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(2,1) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,zpdfns,ZPDFNS24,zpdfsg) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(2,1))) 
                                                                        
      return 
      END                                           
!                                                                       
      SUBROUTINE LH_fixedtalbot_gg(xfunc,x,Q2i,Q2f) 
      IMPLICIT none 
                                                                        
      REAL*8 x,xfunc 
      REAL*8 Q2i,Q2f 
                                                                        
      integer m,j 
      REAL*8 theta,pi,sigma,t,tmp,r 
      COMPLEX*16 s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2) 
                                                                        
      pi=acos(-1d0) 
      m=16 
                                                                        
      t=-dlog(x) 
                                                                        
      tmp  = 0d0 
      tmp2 = 0d0 
                                                                        
      r=2d0*m/5d0/t 
                                                                        
      do j=1,m-1 
         theta = dble(j)*pi/dble(m) 
         sigma=theta+(theta/tan(theta)-1d0)/tan(theta) 
         s=r*theta*DCMPLX(1d0/tan(theta),1d0) 
         s1 = s + (1d0,0d0) 
         call LH_zfunc(s1,x,Q2i,Q2f,zpdfns,ZPDFNS24,zpdfsg) 
         tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(2,2) 
         tmp=tmp+dreal(tmp2) 
      enddo 
      call LH_zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,zpdfns,ZPDFNS24,zpdfsg) 
                                                                        
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(2,2))) 
                                                                        
      return 
      END                                           
                                                                        
!****                                                                   
! 6 *-------------------------------------------------------------------
!****                                                                   
!     File: gamma.f                                                     
!     Function for computing Gamma(2), second moment                    
!     of the evolution factor, for the Singlet-Gluon                    
!     coupled evolution                                                 
!                                                                       
                                                                        
      FUNCTION LH_gamm2ns(nsflav,x,Q2i,Q2f) 
      IMPLICIT none 
!                                                                       
      REAL*8 Q2i,Q2f 
!                                                                       
      REAL*8 Q2it,Q2ft 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
!                                                                       
      integer nsflav 
      COMPLEX*16 zn 
      REAL*8 ain,afin,eps,x 
      REAL*8 LH_gamm2ns,LH_DGAUSS_INTERN 
      REAL*8 LH_gammax_nsp_wrp,LH_gammax_nsm_wrp,LH_gammax_nsvu_wrp,    &
     &     LH_gammax_nsvd_wrp,LH_gammax_nsvct_wrp,LH_gammax_nsvsb_wrp   
      external LH_gammax_nsp_wrp,LH_gammax_nsm_wrp,LH_gammax_nsvu_wrp,  &
     &     LH_gammax_nsvd_wrp,LH_gammax_nsvct_wrp,LH_gammax_nsvsb_wrp   
                                                                        
      parameter(eps=1d-8) 
      COMPLEX*16 efnns(4),efnns24(2),efnsg(2,2) 
                                                                        
!     NB: 2 moment                                                      
      zn = (2d0,0d0) 
!                                                                       
      call LH_ZFUNC(zn,x,Q2i,Q2f,EFNNS,EFNNS24,EFNSG) 
!                                                                       
!     Fill auxiliary variables needed by wrapper functions              
!                                                                       
      Q2it = Q2i 
      Q2ft = Q2f 
!                                                                       
      ain=0d0 
      afin=x 
!                                                                       
      IF(NSFLAV.EQ.1) THEN 
         LH_gamm2ns = DREAL(efnns(1))-                                  &
     &        LH_DGAUSS_INTERN(LH_gammax_nsp_wrp,ain,afin,eps)          
       ELSEIF(NSFLAV.EQ.2) THEN 
         LH_gamm2ns = DREAL(efnns(2))                                   &
     &         -LH_DGAUSS_INTERN(LH_gammax_nsm_wrp,ain,afin,eps)        
      ELSEIF(NSFLAV.EQ.3) THEN 
         LH_gamm2ns = DREAL(efnns(3))                                   &
     &        -LH_DGAUSS_INTERN(LH_gammax_nsvu_wrp,ain,afin,eps)        
      ELSEIF(NSFLAV.EQ.4) THEN 
         LH_gamm2ns = DREAL(efnns(4))                                   &
     &        -LH_DGAUSS_INTERN(LH_gammax_nsvd_wrp,ain,afin,eps)        
      ENDIF 
!                                                                       
      return 
      END                                           
!                                                                       
      FUNCTION LH_gamm2ns24(ns24flav,x,Q2i,Q2f) 
      IMPLICIT none 
!                                                                       
      REAL*8 Q2i,Q2f 
!                                                                       
      REAL*8 Q2it,Q2ft 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
!                                                                       
      integer ns24flav 
      COMPLEX*16 zn 
      REAL*8 ain,afin,eps,x 
      REAL*8 LH_gamm2ns24,LH_DGAUSS_INTERN 
      REAL*8 LH_gammax_ns24q_wrp,LH_gammax_ns24g_wrp 
      external LH_gammax_ns24q_wrp,LH_gammax_ns24g_wrp 
                                                                        
      parameter(eps=1d-8) 
      COMPLEX*16 efnns(4),efnns24(2),efnsg(2,2) 
                                                                        
!     NB: 2 moment                                                      
      zn = (2d0,0d0) 
!                                                                       
      call LH_ZFUNC(zn,x,Q2i,Q2f,EFNNS,EFNNS24,EFNSG) 
!                                                                       
!     Fill auxiliary variables needed by wrapper functions              
!                                                                       
      Q2it = Q2i 
      Q2ft = Q2f 
!                                                                       
      ain=0d0 
      afin=x 
!                                                                       
      IF(ns24flav.EQ.1) THEN 
         LH_gamm2ns24 = DREAL(efnns24(1))                               &
     &               -LH_DGAUSS_INTERN(LH_gammax_ns24q_wrp,ain,afin,eps)
      ELSEIF(ns24flav.EQ.2) THEN 
         LH_gamm2ns24 = DREAL(efnns24(2))                               &
     &               -LH_DGAUSS_INTERN(LH_gammax_ns24g_wrp,ain,afin,eps)
                                                                        
      ENDIF 
!                                                                       
      RETURN 
      END                                           
                                                                        
!                                                                       
      FUNCTION LH_gamm2sg(l,c,x,Q2i,Q2f) 
      IMPLICIT none 
!                                                                       
      REAL*8 Q2i,Q2f 
!                                                                       
      REAL*8 Q2it,Q2ft 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
!                                                                       
      integer l,c 
      COMPLEX*16 zn 
      REAL*8 ain,afin,eps,x 
      REAL*8 LH_gamm2sg,LH_DGAUSS_INTERN 
      REAL*8 LH_gammax_qq_wrp,LH_gammax_qg_wrp 
      REAL*8 LH_gammax_gq_wrp,LH_gammax_gg_wrp 
      external LH_gammax_qq_wrp,LH_gammax_qg_wrp 
      external LH_gammax_gq_wrp,LH_gammax_gg_wrp 
                                                                        
      parameter(eps=1d-8) 
      COMPLEX*16 efnns(4),efnns24(2),efnsg(2,2) 
                                                                        
!     NB: 2 moment                                                      
      zn=(2d0,0d0) 
!                                                                       
      call LH_ZFUNC(zn,x,Q2i,Q2f,EFNNS,EFNNS24,EFNSG) 
!                                                                       
!     Fill auxiliary variables needed by wrapper functions              
!                                                                       
      Q2it = Q2i 
      Q2ft = Q2f 
!                                                                       
      ain=0d0 
      afin=x 
!                                                                       
      IF(l.EQ.1.AND.c.EQ.1) THEN 
         LH_gamm2sg=DREAL(efnsg(1,1))-                                  &
     &        LH_DGAUSS_INTERN(LH_gammax_qq_wrp,ain,afin,eps)           
      ELSEIF(l.EQ.1.AND.c.EQ.2) THEN 
         LH_gamm2sg=DREAL(efnsg(1,2))                                   &
     &        -LH_DGAUSS_INTERN(LH_gammax_qg_wrp,ain,afin,eps)          
      ELSEIF(l.EQ.2.AND.c.EQ.1) THEN 
         LH_gamm2sg=DREAL(efnsg(2,1))                                   &
     &        -LH_DGAUSS_INTERN(LH_gammax_gq_wrp,ain,afin,eps)          
      ELSEIF(l.EQ.2.AND.c.EQ.2) THEN 
         LH_gamm2sg=DREAL(efnsg(2,2))                                   &
     &        -LH_DGAUSS_INTERN(LH_gammax_gg_wrp,ain,afin,eps)          
      ENDIF 
                                                                        
      RETURN 
      END                                           
!                                                                       
!     Internal Wrapper functions needed to perform x-integration        
!     of the Gamma functions using Dgauss_Intern.                       
!                                                                       
      FUNCTION LH_gammax_nsp_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_nsp_wrp 
      REAL*8 y 
      REAL*8 QQ2,QQ0 
      REAL*8 Q2it,Q2ft 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_nsp 
      EXTERNAL LH_gammax_nsp 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_nsp_wrp = LH_gammax_nsp(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
                                                                        
      FUNCTION LH_gammax_nsm_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_nsm_wrp 
      REAL*8 y 
      REAL*8 QQ2,QQ0 
      REAL*8 Q2it,Q2ft 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_nsm 
      EXTERNAL LH_gammax_nsm 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_nsm_wrp = LH_gammax_nsm(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_nsvu_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_nsvu_wrp 
      REAL*8 y 
      REAL*8 QQ2,QQ0 
      REAL*8 Q2it,Q2ft 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_nsvu 
      EXTERNAL LH_gammax_nsvu 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_nsvu_wrp = LH_gammax_nsvu(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_nsvd_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_nsvd_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_nsvd 
      EXTERNAL LH_gammax_nsvd 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_nsvd_wrp = LH_gammax_nsvd(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_ns24q_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_ns24q_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_ns24q 
      EXTERNAL LH_gammax_ns24q 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_ns24q_wrp = LH_gammax_ns24q(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_ns24g_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_ns24g_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_ns24g 
      EXTERNAL LH_gammax_ns24g 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_ns24g_wrp = LH_gammax_ns24g(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_qq_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_qq_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_qq 
      EXTERNAL LH_gammax_qq 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_qq_wrp = LH_gammax_qq(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_qg_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_qg_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_qg 
      EXTERNAL LH_gammax_qg 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_qg_wrp = LH_gammax_qg(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_gq_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_gq_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_gq 
      EXTERNAL LH_gammax_gq 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_gq_wrp = LH_gammax_gq(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_gammax_gg_wrp(y) 
      IMPLICIT none 
!                                                                       
      REAL*8 LH_gammax_gg_wrp 
      REAL*8 y 
      REAL*8 Q2it,Q2ft 
      REAL*8 QQ2,QQ0 
      COMMON / NNPDF10QTRANSFER / Q2it,Q2ft 
      REAL*8 LH_gammax_gg 
      EXTERNAL LH_gammax_qg 
!                                                                       
      QQ0 = Q2IT 
      QQ2 = Q2FT 
      LH_gammax_gg_wrp = LH_gammax_gg(y,QQ0,QQ2) 
!                                                                       
      RETURN 
      END                                           
                                                                        
!****                                                                   
! 7 *-------------------------------------------------------------------
!****                                                                   
!     File: evolfactn.f                                                 
!     Returns the evolution factors in N space                          
!      - Fixed Flavour Number Scheme (IVFN=0)                           
!      - Zero Mass-Variable Flavour Number Scheme (IVFN=1)              
!                                                                       
      SUBROUTINE LH_ZFUNC(ZN,X,Q20,Q2,ZFUNCNS,ZFUNCNS24,ZFUNCSG) 
      IMPLICIT none 
!                                                                       
      INTEGER IPT,IMODEV,IVFN,ITMC 
      COMMON /NNPDF10EVFLAGS/ IPT,IMODEV,IVFN,ITMC 
!                                                                       
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON /nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      INTEGER I,J 
      INTEGER NFI,NFF,STEP 
!                                                                       
      REAL*8 Q20,Q2 
      REAL*8 X 
!                                                                       
      COMPLEX*16 ZN 
      COMPLEX*16 EFNNS(4),EFNSG(2,2) 
      COMPLEX*16 EFNNS0(4),EFNSG0(2,2),EFNSG1(2,2) 
      COMPLEX*16 ZFUNCNS(4),ZFUNCNS24(2),ZFUNCSG(2,2) 
!                                                                       
      COMPLEX*16 EFNNSBQ(4),EFNNSCB(4),EFNNSBT(4),EFNNSTQ(4) 
      COMPLEX*16 EFNSGBQ(2,2),EFNSGCB(2,2),EFNSGBT(2,2),EFNSGTQ(2,2) 
!                                                                       
!     Allow backward evolution                                          
!                                                                       
      IF(Q20.LE.Q2) STEP =  1 
      IF(Q20.GT.Q2) STEP = -1 
!                                                                       
!     Determine # of active flavours at the initial (NFI) and           
!     final (NFF) scale.                                                
!                                                                       
      IF(Q20.GE.Q2TH(6)) THEN 
         NFI=6 
      ELSEIF(Q20.GE.Q2TH(5)) THEN 
         NFI=5 
      ELSEIF(Q20.GE.Q2TH(4)) THEN 
         NFI=4 
      ELSE 
         NFI=3 
      ENDIF 
!                                                                       
      IF(Q2.GT.Q2TH(6)) THEN 
         NFF=6 
      ELSEIF(Q2.GT.Q2TH(5)) THEN 
         NFF=5 
      ELSEIF(Q2.GT.Q2TH(4)) THEN 
         NFF=4 
      ELSE 
         NFF=3 
      ENDIF 
      IF(IVFN.EQ.0) NFF=4 
!                                                                       
!     No call to Coefficient Functions!                                 
!     Evolution Kernel                                                  
      DO I=1,4 
         EFNNS0(I) = (0d0,0d0) 
      ENDDO 
      DO I=1,2 
         DO J=1,2 
            EFNSG0(I,J) = (0.D0,0.D0) 
            EFNSG1(I,J) = (0.D0,0.D0) 
            EFNSG(I,J) = (0.D0,0.D0) 
         END DO 
      END DO 
!                                                                       
!     Fixed Flavour Number Scheme (NF=4)                                
      IF (IVFN.EQ.0) THEN 
         NFI = 4 
         CALL LH_EVOLFACTN(ZN,Q20,Q2,NFI,EFNNS,EFNSG1) 
         ZFUNCNS24(1) = EFNSG1(1,1) 
         ZFUNCNS24(2) = EFNSG1(1,2) 
!                                                                       
!     (Zero Mass-) Variable Flavour Number Scheme                       
      ELSEIF(IVFN.EQ.1) THEN 
         IF(NFF.EQ.NFI) THEN 
            CALL LH_EVOLFACTN(ZN,Q20,Q2,NFI,EFNNS,EFNSG1) 
         ELSE 
            CALL LH_EVOLFACTN(ZN,Q20,Q2TH(NFI+STEP),NFI,EFNNS0,EFNSG0) 
            DO I=1,4 
               EFNNS(I) = EFNNS0(I) 
            ENDDO 
            CALL LH_MEQUAL_C(EFNSG,EFNSG0,2,2) 
   10       NFI = NFI+STEP 
            IF(NFI.NE.NFF) THEN 
               CALL LH_EVOLFACTN(ZN,Q2TH(NFI),Q2TH(NFI+STEP),NFI,       &
     &              EFNNS0,EFNSG0)                                      
            ELSE 
               CALL LH_EVOLFACTN(ZN,Q2TH(NFI),Q2,NFI,EFNNS0,EFNSG0) 
            ENDIF 
            DO I=1,4 
               EFNNS(I) = EFNNS(I)*EFNNS0(I) 
            ENDDO 
            CALL LH_MMULT_C(EFNSG0,2,2,EFNSG,2,2,EFNSG1) 
                                                                        
            IF(NFI.NE.NFF) THEN 
               CALL LH_MEQUAL_C(EFNSG,EFNSG1,2,2) 
               GOTO 10 
            ENDIF 
         ENDIF 
!                                                                       
!     Evolution of T_24                                                 
         IF(Q2.LE.Q2TH(5))THEN 
            ZFUNCNS24(1) = EFNSG1(1,1) 
            ZFUNCNS24(2) = EFNSG1(1,2) 
         ELSEIF(Q2.LE.Q2TH(6)) THEN 
            CALL LH_EVOLFACTN(ZN,Q20,Q2TH(5),4,EFNNSCB,EFNSGCB) 
            CALL LH_EVOLFACTN(ZN,Q2TH(5),Q2 ,5,EFNNSBQ,EFNSGBQ) 
            ZFUNCNS24(1) = EFNNSBQ(1)*EFNSGCB(1,1) 
            ZFUNCNS24(2) = EFNNSBQ(1)*EFNSGCB(1,2) 
         ELSE 
            CALL LH_EVOLFACTN(ZN,Q20,Q2TH(5),4,EFNNSCB,EFNSGCB) 
            CALL LH_EVOLFACTN(ZN,Q2TH(5),Q2TH(6),5,EFNNSBT,EFNSGBT) 
            CALL LH_EVOLFACTN(ZN,Q2TH(6),Q2,6,EFNNSTQ,EFNSGTQ) 
            ZFUNCNS24(1) = EFNNSTQ(1)*EFNNSBT(1)*EFNSGCB(1,1) 
            ZFUNCNS24(2) = EFNNSTQ(1)*EFNNSBT(1)*EFNSGCB(1,2) 
         ENDIF 
      ENDIF 
!                                                                       
!     Output                                                            
      DO I=1,4 
         ZFUNCNS(I) =  EFNNS(I) 
      ENDDO 
                                                                        
      DO I=1,2 
         ZFUNCNS24(I)= ZFUNCNS24(I) 
      ENDDO 
                                                                        
      DO I=1,2 
         DO J=1,2 
            ZFUNCSG(I,J)= EFNSG1(I,J) 
         ENDDO 
      ENDDO 
                                                                        
      RETURN 
      END                                           
!                                                                       
!     Returns the evolution factors, in N space, from Q20 to Q2 with NF 
!     active flavours, for Non Singlet and Singlet-Gluon.               
!                                                                       
      SUBROUTINE LH_EVOLFACTN(ZN,Q2I,Q2F,NF,EFNNS,EFNSG) 
      IMPLICIT none 
!                                                                       
      INTEGER IPT,IMODEV,IVFN,ITMC 
      COMMON /NNPDF10EVFLAGS/ IPT,IMODEV,IVFN,ITMC 
!                                                                       
      REAL*8 q20,q2 
      COMMON/nnpdf10EVSCALE/q20,q2 
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON /nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON/nnpdf10BETA/beta0,beta1,beta2,b1,b2 
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      INTEGER NF 
      INTEGER I,J,K 
      REAL*8 Q2I,Q2F,ASI,ASF,T 
      REAL*8 TMP 
      COMPLEX*16 ZN 
      COMPLEX*16 EFNNS(4),EFNSG(2,2),EFNSGTMP(2,2) 
      COMPLEX*16 EXPNS,EXPM,EXPP 
      COMPLEX*16 UNS0,UNS1(3),LP,LM 
      COMPLEX*16 U(20,2,2) 
      COMPLEX*16 EM(2,2),EP(2,2) 
!                                                                       
      COMPLEX*16 U1(2,2) 
      COMPLEX*16 L(2,2),LU1(2,2),U1L(2,2) 
      COMPLEX*16 USUM(2,2),USUMTMP(2,2),UINV(2,2),USUML(2,2) 
      COMPLEX*16 DETINV 
!                                                                       
!     Compute alphas                                                    
      IF(Q2I.EQ.Q20) THEN 
         ASI = AS0 
      ELSEIF(Q2I.EQ.Q2TH(4)) THEN 
         ASI = ASC 
      ELSEIF(Q2I.EQ.Q2TH(5)) THEN 
         ASI = ASB 
      ELSEIF(Q2I.EQ.Q2TH(6)) THEN 
         ASI = AST 
      ENDIF 
                                                                        
      IF(Q2F.EQ.Q20) THEN 
         ASF = AS0 
      ELSEIF(Q2F.EQ.Q2TH(4)) THEN 
         ASF = ASC 
      ELSEIF(Q2F.EQ.Q2TH(5)) THEN 
         ASF = ASB 
      ELSEIF(Q2F.EQ.Q2TH(6)) THEN 
         ASF = AST 
      ELSE 
         ASF = ASQ 
      ENDIF 
!                                                                       
!     Evaluation of evolution factors                                   
      T = DLOG(ASF/ASI) 
      DO I=1,4 
         EFNNS(I) = (0d0,0d0) 
      ENDDO 
                                                                        
      DETINV=0.D0 
      DO I=1,2 
         DO J=1,2 
            UINV(I,J)=(0d0,0d0) 
            EFNSG(I,J) = (0d0,0d0) 
            EFNSGTMP(I,J) = (0d0,0d0) 
         ENDDO 
      ENDDO 
!                                                                       
!     U matrices                                                        
      CALL LH_UMATRIX(ZN,NF,LP,LM,EP,EM,U,UNS0,UNS1) 
!                                                                       
!     LO evolution factor                                               
!                                                                       
!     Non singlet                                                       
      EXPNS  = EXP ( - UNS0 * T ) 
!                                                                       
!     Singlet                                                           
      EXPM   = EXP ( - LM * T ) 
      EXPP   = EXP ( - LP * T ) 
                                                                        
      DO I=1,2 
        DO J=1,2 
                                                    ! eq.(22)           
           L(I,J) = EXPM * EM(I,J) + EXPP * EP(I,J) 
        ENDDO 
      ENDDO 
                                                                        
      IF(IPT.EQ.0)THEN 
!                                                                       
!     LO solution                                                       
         DO I= 1,4 
            EFNNS(I) = EXPNS 
         ENDDO 
         DO I=1,2 
            DO J=1,2 
               EFNSG(I,J) = L(I,J) 
            ENDDO 
         ENDDO 
                                                                        
      ELSEIF(IPT.EQ.1)THEN 
!                                                                       
!     NLO solution                                                      
         IF(IMODEV.EQ.0)THEN 
!                                                                       
!     Truncated solution IMODEV=0                                       
            DO I=1,3 
               EFNNS(I)= EXPNS * ((1D0,0D0) + UNS1(I) * (ASF-ASI)) 
            ENDDO 
            EFNNS(4)=EFNNS(3) 
                                                                        
            DO I=1,2 
               DO J=1,2 
                  U1(I,J)=U(1,I,J) 
               ENDDO 
            ENDDO 
                                                                        
            CALL LH_MMULT_C(U1,2,2,L,2,2,U1L) 
            CALL LH_MMULT_C(L,2,2,U1,2,2,LU1) 
            DO I=1,2 
               DO J=1,2 
                  EFNSG(I,J) = L(I,J) + ASF * U1L(I,J) - ASI * LU1(I,J) 
               ENDDO 
            ENDDO 
                                                                        
         ELSEIF(IMODEV.EQ.1)THEN 
!                                                                       
!     Iterated solution IMODEV=1                                        
            TMP = 0D0 
            TMP = DLOG(( 1D0 + B1(NF) * ASF )/(1D0 + B1(NF)*ASI))/B1(NF) 
            DO I=1,3 
               EFNNS(I) = EXPNS * EXP( TMP * UNS1(I)) 
            ENDDO 
            EFNNS(4)=EFNNS(3) 
!                                                                       
            USUM(1,1) = (1d0,0d0) 
            USUM(1,2) = (0d0,0d0) 
            USUM(2,1) = (0d0,0d0) 
            USUM(2,2) = (1d0,0d0) 
            USUMTMP(1,1) = (1d0,0d0) 
            USUMTMP(1,2) = (0d0,0d0) 
            USUMTMP(2,1) = (0d0,0d0) 
            USUMTMP(2,2) = (1d0,0d0) 
!                                                                       
            DO K=1,20 
               DO I=1,2 
                  DO J=1,2 
                     USUM(I,J) = USUM(I,J)+                             &
     &                    (ASF**(DBLE(K)) * U(K,I,J))                   
                     USUMTMP(I,J) = USUMTMP(I,J) +                      &
     &                    (ASI**(DBLE(K)) * U(K,I,J))                   
                   ENDDO 
               ENDDO 
            ENDDO 
!                                                                       
            DETINV = 1.D0 / ( USUMTMP(1,1) * USUMTMP(2,2)-              &
     &           USUMTMP(1,2) * USUMTMP(2,1) )                          
            UINV(1,1) = DETINV * USUMTMP(2,2) 
            UINV(1,2) = - DETINV * USUMTMP(1,2) 
            UINV(2,1) = - DETINV * USUMTMP(2,1) 
            UINV(2,2) = DETINV * USUMTMP(1,1) 
!                                                                       
            CALL LH_MMULT_C(USUM,2,2,L,2,2,USUML) 
            CALL LH_MMULT_C(USUML,2,2,UINV,2,2,EFNSGTMP) 
                                                                        
            DO I=1,2 
               DO J=1,2 
                  EFNSG(I,J)= EFNSGTMP(I,J) 
               ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END                                           
                                                                        
!                                                                       
!     File: umatrix.f                                                   
!                                                                       
!     Computes the matrices used for the singlet evolution.             
!     Equation numbers refer to the "Notes on PT evolution"             
!                                                                       
      SUBROUTINE LH_UMATRIX(ZN,NF,LP,LM,EP,EM,U,UNS0,UNS1) 
      IMPLICIT none 
!                                                                       
      INTEGER IPT,IMODEV,IVFN,ITMC 
      COMMON /NNPDF10EVFLAGS/ IPT,IMODEV,IVFN,ITMC 
!                                                                       
      REAL*8 beta0(3:6),beta1(3:6),beta2(3:6) 
      REAL*8 b1(3:6),b2(3:6) 
      COMMON/nnpdf10BETA/beta0,beta1,beta2,b1,b2 
!                                                                       
      INTEGER NF 
      COMPLEX*16 ZN 
      COMPLEX*16 UNS0,UNS1(3) 
      COMPLEX*16 U(20,2,2) 
      COMPLEX*16 LP,LM 
      COMPLEX*16 EP(2,2),EM(2,2) 
!                                                                       
      INTEGER I,J,K,JJ,KK,II 
      COMPLEX*16 QQ0,QG0,GQ0,GG0 
      COMPLEX*16 SQ,LDIFF 
      COMPLEX*16 DP(2,2),DM(2,2) 
      COMPLEX*16 G0(2,2),G1(2,2) 
      COMPLEX*16 G0NS(3),G1NS(3) 
      COMPLEX*16 R(20,2,2),RT(20,2,2) 
      COMPLEX*16 R0(2,2),R1(2,2) 
      COMPLEX*16 R1PP(2,2),R1PM(2,2),R1MP(2,2),R1MM(2,2) 
      COMPLEX*16 U1PP(2,2),U1PM(2,2),U1MP(2,2),U1MM(2,2) 
!                                                                       
      CALL LH_ANDIM_LO(ZN,NF,G0NS,G0) 
      UNS0 = G0NS(1) / BETA0(NF) 
!                                                                       
      DO K=1,20 
         DO I=1,2 
            DO J=1,2 
               R(K,I,J)=(0.D0,0.D0) 
               RT(K,I,J)=(0.D0,0.D0) 
               U(K,I,J)=(0.D0,0.D0) 
            ENDDO 
         ENDDO 
      ENDDO 
!                                                                       
      DO I=1,2 
         DO J=1,2 
                                          ! eq. (18)                    
            R0(I,J) = G0(I,J) / BETA0(NF) 
         END DO 
      END DO 
!                                                                       
      QQ0 = R0(1,1) 
      QG0 = R0(1,2) 
      GQ0 = R0(2,1) 
      GG0 = R0(2,2) 
!                                                                       
                                                 ! eq.(20)              
      SQ = SQRT((QQ0-GG0)**2d0 + 4.d0*QG0*GQ0) 
      LP = .5d0*(QQ0 + GG0 + SQ) 
      LM = .5d0*(QQ0 + GG0 - SQ) 
!                                                                       
                                               ! eq.(21)                
      EM(1,1) = (QQ0 - LP) / (LM - LP) 
      EM(1,2)= QG0 / (LM - LP) 
      EM(2,1) = GQ0 / (LM - LP) 
      EM(2,2) = (GG0 - LP) / (LM - LP) 
!                                                                       
      EP(1,1) = (QQ0 - LM) / (LP - LM) 
      EP(1,2) = QG0 / (LP - LM) 
      EP(2,1) = GQ0 / (LP - LM) 
      EP(2,2) = (GG0 - LM) / (LP - LM) 
!                                                                       
!     NLO coefficients                                                  
      IF (IPT.GE.1) THEN 
         CALL LH_ANDIM_NLO(ZN,NF,G1NS,G1) 
!                                                                       
!     NON SINGLET                                                       
         DO I=1,3 
            UNS1(I)= -(G1NS(I)/BETA0(NF)) + B1(NF)*UNS0 
         ENDDO 
!                                                                       
!     SINGLET                                                           
         DO I=1,2 
            DO J=1,2 
                                                              ! eq. (18)
               R(1,I,J)  = G1(I,J)/BETA0(NF) - B1(NF)*R0(I,J) 
                                                               ! R1_T=R1
               RT(1,I,J)  = G1(I,J)/BETA0(NF) - B1(NF)*R0(I,J) 
               R1(I,J) = R(1,I,J) 
           ENDDO 
         ENDDO 
!                                                                       
!     Computation of U1 according to eq.(25)                            
         CALL LH_MMULT_C(EP,2,2,R1,2,2,DP) 
         CALL LH_MMULT_C(EM,2,2,R1,2,2,DM) 
!                                                                       
         CALL LH_MMULT_C(DP,2,2,EP,2,2,R1PP) 
         CALL LH_MMULT_C(DP,2,2,EM,2,2,R1PM) 
         CALL LH_MMULT_C(DM,2,2,EP,2,2,R1MP) 
         CALL LH_MMULT_C(DM,2,2,EM,2,2,R1MM) 
!                                                                       
         DO I = 1,2 
            DO J = 1,2 
               U1PP(I,J) = R1PP(I,J) 
               U1PM(I,J) = R1PM(I,J)/(LM - LP - 1d0) 
               U1MP(I,J) = R1MP(I,J)/(LP - LM - 1d0) 
               U1MM(I,J) = R1MM(I,J) 
               U(1,I,J) = -(U1MM(I,J)+U1PP(I,J)) + U1PM(I,J) + U1MP(I,J) 
            ENDDO 
         ENDDO 
!                                                                       
         IF(IMODEV.EQ.0)THEN 
            DO K=2,20 
               DO I=1,2 
                  DO J=1,2 
                     U(K,I,J)=0.D0 
                  ENDDO 
               ENDDO 
            ENDDO 
                                                                        
         ELSEIF(IMODEV.EQ.1)THEN 
!                                                                       
!     Computation of U_k (k>=2) at NLO according to eqs (23,26,27)      
            LDIFF = LM - LP 
                                                                        
            DO K=2,20 
               DO I=1,2 
                  DO J=1,2 
                                                      ! eq (27)         
                     R(K,I,J) = - B1(NF) * R(K-1,I,J) 
                  ENDDO 
               ENDDO 
            ENDDO 
                                                                        
            DO K=2,20 
               DO I=1,2 
                  DO J=1,2 
                     RT(K,I,J) = R(K,I,J) 
                     DO JJ=1,K-1 
                        DO KK=1,2 
                           RT(K,I,J) = RT(K,I,J) +                      &
     &                          ( R(JJ,I,KK) * U(K-JJ,KK,J) )           
                                                              ! eq (26) 
                        ENDDO 
                     ENDDO 
                  ENDDO 
               ENDDO 
!                                                                       
               DO I=1,2 
                  DO J=1,2 
                     DO II=1,2 
                        DO JJ=1,2 
                                               ! eq(23)                 
                           U(K,I,J) = U(K,I,J)                          &
     &                          -(EM(I,II)*RT(K,II,JJ)*EM(JJ,J)/DBLE(K))&
     &                          -(EP(I,II)*RT(K,II,JJ)*EP(JJ,J)/DBLE(K))&
     &                          -(EM(I,II)*RT(K,II,JJ)* EP(JJ,J)/       &
     &                          ( DBLE(K) + LDIFF ))                    &
     &                          -(EP(I,II)*RT(K,II,JJ)                  &
     &                          * EM(JJ,J) / ( DBLE(K) - LDIFF ))       
                        ENDDO 
                     ENDDO 
                  ENDDO 
               ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END                                           
                                                                        
!****                                                                   
! 8 *-------------------------------------------------------------------
!****                                                                   
!                                                                       
!     EVOLINT.F Evaluate convolution in x space (DIFF FROM SRC CODE: NO 
!                                                                       
      SUBROUTINE LH_PDFEVOLX(x,pdftmp) 
      IMPLICIT none 
!                                                                       
      INTEGER IPT,IMODEV,IVFN,ITMC 
      COMMON /NNPDF10EVFLAGS/ IPT,IMODEV,IVFN,ITMC 
      REAL*8 EPSEVOL 
      COMMON /NNPDF10EVOLACC/ EPSEVOL 
!                                                                       
      REAL*8 gm2ns(4),gm2ns24(2),gm2sg(2,2) 
      COMMON /nnpdf10GM2/gm2ns,gm2ns24,gm2sg 
!                                                                       
      REAL*8 q20,q2 
      COMMON/nnpdf10EVSCALE/q20,q2 
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON /nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      INTEGER i,j,k 
!                                                                       
      REAL*8 x 
      REAL*8 pdfevolx1(13),pdf0(13),pdftmp(13) 
      REAL*8 pdfevolx24(2) 
      REAL*8 LH_pdfevolxint_ns,LH_pdfevolxint_sg,LH_pdfevolxint_ns24 
      EXTERNAL LH_pdfevolxint_ns,LH_pdfevolxint_ns24,LH_pdfevolxint_sg 
!                                                                       
      call LH_PDFIN(x,pdf0) 
!                                                                       
      do i=1,13 
         pdfevolx1(i)=0d0 
         pdftmp(i)=0d0 
      enddo 
!                                                                       
      DO k=1,2 
         DO j=1,2 
            pdfevolx1(j) = LH_pdfevolxint_sg(k,j,x) 
         ENDDO 
         pdftmp(k) = pdfevolx1(1)+gm2sg(k,1)*pdf0(1)+                   &
     &                    pdfevolx1(2)+gm2sg(k,2)*pdf0(2)               
      ENDDO 
!                                                                       
      pdfevolx1(3) = LH_pdfevolxint_ns(3,3,x) 
      pdftmp(3)    = pdfevolx1(3) + gm2ns(3)*pdf0(3) 
!                                                                       
      pdfevolx1(4) = LH_pdfevolxint_ns(4,4,x) 
      pdftmp(4)    = pdfevolx1(4) + gm2ns(4)*pdf0(4) 
!                                                                       
      pdfevolx1(5) = LH_pdfevolxint_ns(4,5,x) 
      pdftmp(5)    = pdfevolx1(5) + gm2ns(4)*pdf0(5) 
!                                                                       
      pdfevolx1(6) = LH_pdfevolxint_ns(3,6,x) 
      pdftmp(6)    = pdfevolx1(6) + gm2ns(3)*pdf0(6) 
!                                                                       
      pdfevolx1(7) = LH_pdfevolxint_ns(4,7,x) 
      pdftmp(7)    = pdfevolx1(7) + gm2ns(4)*pdf0(7) 
!                                                                       
      pdfevolx1(8) = LH_pdfevolxint_ns(3,8,x) 
      pdftmp(8)    = pdfevolx1(8) + gm2ns(3)*pdf0(8) 
!                                                                       
      DO k=9,11 
         pdfevolx1(k) = LH_pdfevolxint_ns(1,k,x) 
         pdftmp(k)    = pdfevolx1(k) + gm2ns(1)*pdf0(k) 
      ENDDO 
!                                                                       
      IF(IVFN.EQ.0) THEN 
         DO k=12,13 
            pdftmp(k)    = pdftmp(1) 
         ENDDO 
      ELSE 
         DO j=1,2 
            pdfevolx24(j) = LH_pdfevolxint_ns24(j,x) 
         ENDDO 
         pdftmp(12) = pdfevolx24(1) + gm2ns24(1)*pdf0(1)                &
     &        + pdfevolx24(2) + gm2ns24(2)*pdf0(2)                      
         pdftmp(13) = pdftmp(1) 
      ENDIF 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_pdfevolxint_sg(ll,cc,xc) 
      IMPLICIT none 
!                                                                       
      INTEGER npoints 
      PARAMETER(npoints=2**8) 
      INTEGER ieval,niter,nmax 
      COMMON /nnpdf10NGRID/ ieval,niter,nmax 
      REAL*8 EPSEVOL 
      COMMON /NNPDF10EVOLACC/ EPSEVOL 
!                                                                       
!  qnt ridefinite gm2ns(4,ndata)-> gm2ns(4)                             
      REAL*8 xx(NPOINTS),wn(NPOINTS),evkns(4,NPOINTS),                  &
     &     evksg(2,2,NPOINTS),evkns24(2,NPOINTS)                        
      COMMON /nnpdf10XEVK/ xx,wn,evkns,evksg,evkns24 
!                                                                       
      REAL*8 q20,q2 
      COMMON/nnpdf10EVSCALE/q20,q2 
      REAL*8 q2th(4:6),asref,q2ref 
      REAL*8 as0,asc,asb,ast,asq 
      COMMON/nnpdf10VFNS/ q2th,asref,q2ref 
      COMMON /NNPDF10AS/ as0,asc,asb,ast,asq 
!                                                                       
      REAL*8 LH_pdfevolxint_sg 
      INTEGER NGAUSS 
      PARAMETER(NGAUSS=4) 
      INTEGER ll,cc 
      INTEGER j,k,idum,kk 
      INTEGER nint 
      REAL*8 xc,EPSTMP 
      REAL*8 trap1,trap2 
      REAL*8 y,pdf0(13),pdf1(13) 
!                                                                       
      CALL LH_PDFIN(xc,pdf0) 
!                                                                       
      trap2=0.d0 
      idum=0 
      k=0 
      kk=0 
!                                                                       
  187 NINT=2**kk 
      kk=kk+1 
      trap1=trap2 
      trap2=0.d0 
!                                                                       
      DO j=1+k,NGAUSS*NINT+k 
         idum=idum+1 
         y=xx(j) 
         CALL LH_PDFIN(xc/y,pdf1) 
         trap2 = trap2                                                  &
     &          + wn(j)*evksg(ll,cc,j)*(pdf1(cc)/y-y*pdf0(cc))/y        
      ENDDO 
                                                                        
      k=idum 
!                                                                       
      EPSTMP=DABS((trap2-trap1)/trap2) 
      IF(EPSTMP.GT.EPSEVOL.and.idum.lt.NMAX)then 
         GOTO 187 
      ELSE 
         LH_pdfevolxint_sg=trap2 
      ENDIF 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_pdfevolxint_ns(nsflag,iflav,xc) 
      IMPLICIT none 
!                                                                       
      INTEGER npoints 
      PARAMETER(npoints=2**8) 
      INTEGER ieval,niter,nmax 
      COMMON /nnpdf10NGRID/ ieval,niter,nmax 
      REAL*8 EPSEVOL 
      COMMON /NNPDF10EVOLACC/ EPSEVOL 
!                                                                       
!  qnt ridefinite gm2ns(4,ndata)-> gm2ns(4)                             
      REAL*8 xx(NPOINTS),wn(NPOINTS),evkns(4,NPOINTS),                  &
     &     evksg(2,2,NPOINTS),evkns24(2,NPOINTS)                        
      COMMON /nnpdf10XEVK/ xx,wn,evkns,evksg,evkns24 
!                                                                       
      REAL*8 LH_PDFEVOLXINT_NS 
!                                                                       
      INTEGER nsflag,iflav 
      INTEGER j,k,idum,kk 
      INTEGER nint 
      REAL*8 xc,EPSTMP 
      REAL*8 trap1,trap2 
      REAL*8 y,pdf0(13),pdf1(13) 
      INTEGER NGAUSS 
      PARAMETER(NGAUSS=4) 
!                                                                       
      call LH_PDFIN(xc,pdf0) 
!                                                                       
      trap2=0.d0 
      idum=0 
      k=0 
      kk=0 
!                                                                       
  187 NINT=2**kk 
      kk=kk+1 
      trap1=trap2 
      trap2=0.d0 
!                                                                       
      do j=1+k,NGAUSS*NINT+k 
         idum=idum+1 
         y=xx(j) 
         call LH_PDFIN(xc/y,pdf1) 
                                                                        
         trap2 = trap2                                                  &
     &     + wn(j)*evkns(nsflag,j)*(pdf1(iflav)/y-y*pdf0(iflav))/y      
      enddo 
      k=idum 
!                                                                       
      EPSTMP=DABS((trap2-trap1)/trap2) 
      IF(EPSTMP.GT.EPSEVOL.and.idum.lt.NMAX)then 
         goto 187 
      else 
         LH_pdfevolxint_ns=trap2 
      endif 
!                                                                       
      RETURN 
      END                                           
!                                                                       
      FUNCTION LH_pdfevolxint_ns24(ns24flag,xc) 
      IMPLICIT none 
!                                                                       
      INTEGER npoints 
      PARAMETER(npoints=2**8) 
      INTEGER ieval,niter,nmax 
      COMMON /nnpdf10NGRID/ ieval,niter,nmax 
      REAL*8 EPSEVOL 
      COMMON /NNPDF10EVOLACC/ EPSEVOL 
!                                                                       
!  qnt ridefinite gm2ns(4,ndata)-> gm2ns(4)                             
      REAL*8 xx(NPOINTS),wn(NPOINTS),evkns(4,NPOINTS),                  &
     &     evksg(2,2,NPOINTS),evkns24(2,NPOINTS)                        
      COMMON /nnpdf10XEVK/ xx,wn,evkns,evksg,evkns24 
!                                                                       
      REAL*8 LH_pdfevolxint_ns24 
!                                                                       
      INTEGER ns24flag 
      INTEGER j,k,idum,kk 
      INTEGER nint 
      REAL*8 xc,epstmp 
      REAL*8 trap1,trap2 
      REAL*8 y,pdf0(13),pdf1(13) 
      INTEGER NGAUSS 
      PARAMETER(NGAUSS=4) 
!                                                                       
      call LH_PDFIN(xc,pdf0) 
!                                                                       
      trap2=0.d0 
      idum=0 
      k=0 
      kk=0 
!                                                                       
  187 NINT=2**kk 
      kk=kk+1 
      trap1=trap2 
      trap2=0.d0 
!                                                                       
      do j=1+k,NGAUSS*NINT+k 
         idum=idum+1 
         y=xx(j) 
         call LH_PDFIN(xc/y,pdf1) 
         trap2 = trap2                                                  &
     &        + wn(j)*evkns24(ns24flag,j)                               &
     &        *(pdf1(ns24flag)/y-y*pdf0(ns24flag))/y                    
      enddo 
      k=idum 
!                                                                       
      EPSTMP=DABS((trap2-trap1)/trap2) 
      IF(EPSTMP.GT.EPSEVOL.and.idum.lt.NMAX)then 
         goto 187 
      else 
         LH_pdfevolxint_ns24=trap2 
      endif 
!                                                                       
      RETURN 
      END                                           
                                                                        
!****                                                                   
! 9 *-------------------------------------------------------------------
!****                                                                   
!                                                                       
!     Change of basis from NN basis to pdfparam                         
!     based on Flavour Assumptions see eq. 86 of notes PTevol.pdf       
!                                                                       
!     Convert from the INPUTPARAMETRIZATION convention for PDF ordering:
!     to the  EVOLUTION convention:                                     
                                                                        
!       1    2   3   4   5   6   7   8    9    10   11    12    13      
!     Sigma  g  u_v d_v s_v c_v b_v t_v  T_3  T_8  T_15  T_24  T_35     
                                                                        
!     Convention for pdfs parametrized with nets                        
!     PDFIN(1) -> Singlet                                               
!     PDFIN(2) -> Gluon                                                 
!     PDFIN(3) -> Triplet T_3                                           
!     PDFIN(4) -> TotalValence V                                        
!     PDFIN(5) -> SeaAsymmetry Delta_S                                  
                                                                        
!     Convention for input PDFs to evolution equations                  
!     PDFOUT(1) = Singlet                                               
!     PDFOUT(2) = Gluon                                                 
!     PDFOUT(3) = u_v = ( V + T_3 + 2 Delta_S )/2                       
!     PDFOUT(4) = d_v = ( V - T_3 - 2 Delta_S )/2                       
!     PDFOUT(5) = s_v = 0                                               
!     PDFOUT(6) = c_v = 0                                               
!     PDFOUT(7) = b_v = 0                                               
!     PDFOUT(8) = t_v = 0                                               
!     PDFOUT(9) = T_3                                                   
!     PDFOUT(10)= T_8 = ( 1 / ( 2 + C ) ) * ( 3*C*V + 2(1-C)*Singlet )  
!     PDFOUT(11) = T_15 = Singlet                                       
!     PDFOUT(12) = T_24 = Singlet                                       
!     PDFOUT(13) = T_35 = Singlet                                       
                                                                        
!***************************                                            
!                                                                       
      SUBROUTINE LH_pdfinpar2evln(PDFIN,PDFOUT) 
      implicit none 
!                                                                       
      integer I,J 
      integer MXPDF 
      parameter(MXPDF=13) 
      integer NTOTPDF 
      parameter(NTOTPDF=5) 
!                                                                       
      REAL*8 PDFIN(NTOTPDF),PDFOUT(MXPDF) 
      REAL*8 par2evln(mxpdf,ntotpdf) 
      COMMON/nnpdf10PPAR2EVLN/PAR2EVLN 
!                                                                       
      CALL LH_initpar2evln 
!                                                                       
      DO I=1,MXPDF 
         PDFOUT(I) = 0D0 
         DO J=1,NTOTPDF 
            PDFOUT(I) = PDFOUT(I) + PAR2EVLN(I,J) * PDFIN(J) 
         ENDDO 
      ENDDO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
      SUBROUTINE LH_initpar2evln 
      IMPLICIT none 
!                                                                       
      integer I,J 
      integer MXPDF 
      parameter(MXPDF=13) 
      integer NTOTPDF 
      parameter(NTOTPDF=5) 
      integer MXPAR 
      parameter(MXPAR=2e2) 
      real*8 CS 
      parameter(CS = 0.5) 
!                                                                       
      REAL*8 par2evln(mxpdf,ntotpdf) 
      COMMON/nnpdf10PPAR2EVLN/PAR2EVLN 
!                                                                       
      do i = 1,mxpdf 
         do j = 1, ntotpdf 
            par2evln(i,j) = 0d0 
         enddo 
      enddo 
      PAR2EVLN(1,1) = 1d0 
      PAR2EVLN(2,2) = 1d0 
      PAR2EVLN(3,3) = 0.5d0 
      PAR2EVLN(3,4) = 0.5d0 
      PAR2EVLN(3,5) = 1d0 
      PAR2EVLN(4,3) = - 0.5d0 
      PAR2EVLN(4,4) = 0.5d0 
      PAR2EVLN(4,5) = - 1d0 
      PAR2EVLN(9,3) = 1d0 
      PAR2EVLN(10,1) = 2.d0*(1d0 - CS)/(2d0 + CS) 
      PAR2EVLN(10,4) = 3d0*CS/(2d0 + CS) 
      PAR2EVLN(11,1) = 1d0 
      PAR2EVLN(12,1) = 1d0 
      PAR2EVLN(13,1) = 1d0 
!                                                                       
      return 
      END                                           
!                                                                       
!     Change of basis from our basis to LHAPDF one:                     
!     Computes the PDFs of the single partons starting from the T_i and 
!     combinations which are used in the evolution.                     
!                                                                       
!     The combinations used in the evolution are numbered               
!     according to the following table:                                 
!                                                                       
!       1    2   3   4   5   6   7   8    9    10   11    12    13      
!     Sigma  g  u_v d_v s_v c_v b_v t_v  T_3  T_8  T_15  T_24  T_35     
!                                                                       
!     The output PDFs are numbered according to the Les Houches         
!     convention:                                                       
!                                                                       
!      -6   -5   -4   -3   -2   -1   0   1   2   3   4   5   6          
!     tbar bbar cbar sbar ubar dbar  g   d   u   s   c   b   t          
!                                                                       
!                                                                       
      SUBROUTINE LH_PDFEVLN2LHA(pdfin,pdfout) 
      IMPLICIT none 
!                                                                       
      INTEGER I,J 
      integer MXPDF 
      parameter(MXPDF=13) 
!                                                                       
      REAL*8 pdfin(MXPDF),pdfout(-6:6) 
      REAL*8 evln2lha(mxpdf,mxpdf) 
      COMMON/nnpdf10EEVLN2LHA/EVLN2LHA 
!                                                                       
      CALL LH_INITEVLN2LHA 
!                                                                       
      DO I = 1,MXPDF 
         PDFOUT(I-7)=0D0 
         DO J = 1,MXPDF 
            PDFOUT(I-7) = PDFOUT(I-7) + EVLN2LHA(I,J)*PDFIN(J) 
         ENDDO 
      ENDDO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!                                                                       
      SUBROUTINE LH_INITEVLN2LHA 
      IMPLICIT none 
!                                                                       
      INTEGER I,J 
!                                                                       
      INTEGER IPT,IMODEV,IVFN,ITMC 
      COMMON /NNPDF10EVFLAGS/ IPT,IMODEV,IVFN,ITMC 
      integer MXPDF 
      parameter(MXPDF=13) 
!                                                                       
      REAL*8 evln2lha(mxpdf,mxpdf) 
      COMMON/nnpdf10EEVLN2LHA/EVLN2LHA 
!                                                                       
!     Write it as a matrix for LHAPDF interface                         
!     PDFOUT(I-7) = EVLN2LHA(I,J) * PDFIN(J)                            
!                                                                       
      IF(IVFN.EQ.1) THEN 
!                                                                       
         do i = 1,mxpdf 
            do j = 1,mxpdf 
               evln2lha(i,j) = 0d0 
            enddo 
         enddo 
!                                                                       
         DO i=1,mxpdf 
            evln2lha(i,1) = 10d0 
         ENDDO 
         evln2lha(7,1)  = 0d0 
!                                                                       
         evln2lha(7,2)  = 120d0 
!                                                                       
         evln2lha(5,3)  = - 60d0 
         evln2lha(9,3)  = + 60d0 
!                                                                       
         evln2lha(6,4)  = - 60d0 
         evln2lha(8,4)  = + 60d0 
!                                                                       
         evln2lha(4,5)  = - 60d0 
         evln2lha(10,5) = + 60d0 
!                                                                       
         evln2lha(3,6)  = - 60d0 
         evln2lha(11,6) = + 60d0 
!                                                                       
         evln2lha(2,7)  = - 60d0 
         evln2lha(12,7) = + 60d0 
!                                                                       
         evln2lha(1,8)  = - 60d0 
         evln2lha(13,8) = + 60d0 
!                                                                       
         evln2lha(5,9)  = + 30d0 
         evln2lha(6,9)  = - 30d0 
         evln2lha(8,9)  = - 30d0 
         evln2lha(9,9)  = + 30d0 
!                                                                       
         evln2lha(4,10) = - 20d0 
         evln2lha(5,10) = + 10d0 
         evln2lha(6,10) = + 10d0 
         evln2lha(8,10) = + 10d0 
         evln2lha(9,10) = + 10d0 
         evln2lha(10,10)= - 20d0 
!                                                                       
         evln2lha(3,11) = - 15d0 
         evln2lha(4,11) = + 5d0 
         evln2lha(5,11) = + 5d0 
         evln2lha(6,11) = + 5d0 
         evln2lha(8,11) = + 5d0 
         evln2lha(9,11) = + 5d0 
         evln2lha(10,11)= + 5d0 
         evln2lha(11,11)= - 15d0 
!                                                                       
         evln2lha(2,12) = - 12d0 
         evln2lha(3,12) = + 3d0 
         evln2lha(4,12) = + 3d0 
         evln2lha(5,12) = + 3d0 
         evln2lha(6,12) = + 3d0 
         evln2lha(8,12) = + 3d0 
         evln2lha(9,12) = + 3d0 
         evln2lha(10,12)= + 3d0 
         evln2lha(11,12)= + 3d0 
         evln2lha(12,12)= - 12d0 
!                                                                       
         evln2lha(1,13) = - 10d0 
         evln2lha(2,13) = + 2d0 
         evln2lha(3,13) = + 2d0 
         evln2lha(4,13) = + 2d0 
         evln2lha(5,13) = + 2d0 
         evln2lha(6,13) = + 2d0 
         evln2lha(8,13) = + 2d0 
         evln2lha(9,13) = + 2d0 
         evln2lha(10,13)= + 2d0 
         evln2lha(11,13)= + 2d0 
         evln2lha(12,13)= + 2d0 
         evln2lha(13,13)= - 10d0 
!                                                                       
         do i = 1,mxpdf 
            do j = 1,mxpdf 
               evln2lha(i,j) = evln2lha(i,j) / 120d0 
            enddo 
         enddo 
                                                                        
      ELSE 
!                                                                       
         do i = 1,mxpdf 
            do j = 1,mxpdf 
               evln2lha(i,j) = 0d0 
            enddo 
         enddo 
!                                                                       
         DO i=1,mxpdf 
            evln2lha(i,1) = 3d0 
         ENDDO 
         evln2lha(7,1)  = 0d0 
!                                                                       
         evln2lha(7,2)  = 24d0 
!                                                                       
         evln2lha(5,3)  = - 12d0 
         evln2lha(9,3)  = + 12d0 
!                                                                       
         evln2lha(6,4)  = - 12d0 
         evln2lha(8,4)  = + 12d0 
!                                                                       
         evln2lha(4,5)  = - 12d0 
         evln2lha(10,5) = + 12d0 
!                                                                       
         evln2lha(3,6)  = - 12d0 
         evln2lha(11,6) = + 12d0 
!                                                                       
         evln2lha(5,9)  = + 6d0 
         evln2lha(6,9)  = - 6d0 
         evln2lha(8,9)  = - 6d0 
         evln2lha(9,9)  = + 6d0 
!                                                                       
         evln2lha(4,10) = - 4d0 
         evln2lha(5,10) = + 2d0 
         evln2lha(6,10) = + 2d0 
         evln2lha(8,10) = + 2d0 
         evln2lha(9,10) = + 2d0 
         evln2lha(10,10)= - 4d0 
!                                                                       
         evln2lha(3,11) = - 3d0 
         evln2lha(4,11) = + 1d0 
         evln2lha(5,11) = + 1d0 
         evln2lha(6,11) = + 1d0 
         evln2lha(8,11) = + 1d0 
         evln2lha(9,11) = + 1d0 
         evln2lha(10,11)= + 1d0 
         evln2lha(11,11)= - 3d0 
!                                                                       
         do i = 1,mxpdf 
            do j = 1,mxpdf 
               evln2lha(i,j) = evln2lha(i,j) / 24d0 
            enddo 
         enddo 
!                                                                       
      ENDIF 
!                                                                       
      RETURN 
      END                                           
                                                                        
!                                                                       
!     pdfevln2inpar.f                                                   
!                                                                       
!     Convert from the  EVOLUTION convention:                           
!       1    2   3   4   5   6   7   8    9    10   11    12    13      
!     Sigma  g  u_v d_v s_v c_v b_v t_v  T_3  T_8  T_15  T_24  T_35     
!                                                                       
!     to the INPUTPARAMETRIZATION convention for PDF ordering:          
!     based on Flavour Assumptions                                      
!                                                                       
                                                                        
      subroutine LH_PDFEVLN2INPAR(PDFIN,PDFOUT) 
      implicit none 
!                                                                       
      integer MXPDF 
      parameter(MXPDF=13) 
      integer NTOTPDF 
      parameter(NTOTPDF=5) 
      real*8 PDFIN(MXPDF),PDFOUT(NTOTPDF) 
                                                                        
!     JPDF=1 -> Singlet                                                 
      PDFOUT(1) = PDFIN(1) 
                                                                        
!     JPDF=2 -> Gluon                                                   
      PDFOUT(2) = PDFIN(2) 
                                                                        
!     JPDF=3 -> Triplet                                                 
      PDFOUT(3) = PDFIN(9) 
                                                                        
!     JPDF=4 -> Total Valence -> V_T = u_V + d_V                        
      PDFOUT(4) = PDFIN(3) + PDFIN(4) 
                                                                        
!     JPDF=5 -> Sea Asymmetry -> Delta_S = ( u_V - d_V - T_3 )/2        
      PDFOUT(5) = ( PDFIN(3) - PDFIN(4) - PDFIN(9) )/2d0 
                                                                        
      return 
      END                                           
                                                                        
!*****                                                                  
! 10 *------------------------------------------------------------------
!*****                                                                  
!                                                                       
!     File: andim_lo.f                                                  
!     Returns the LO anomalous dimensions in the                        
!                                                                       
       SUBROUTINE LH_ANDIM_LO(N,NF,P0NS,P0SG) 
       IMPLICIT none 
!                                                                       
       REAL*8 ca,cf,tr 
       COMMON/NNPDF10COLFACT/ca,cf,tr 
       REAL*8 emc,zeta2,zeta3,zeta4,pi 
       PARAMETER (pi    = 3.1415926535897932385) 
       COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
       COMPLEX*16 LH_PSI, S1 
!                                                                       
       COMPLEX*16 NS,N1,N2,NM 
       COMPLEX*16 PQQA,PQGA,PGQA,PGGA,PGGB 
       COMPLEX*16 N 
       INTEGER NF 
       COMPLEX*16 P0NS 
       COMPLEX*16 P0SG (2,2) 
!                                                                       
       NS = N * N 
       N1 = N + cmplx(1d0,0d0) 
       N2 = N + (2d0,0d0) 
       NM = N - (1d0,0d0) 
!                                                                       
       S1 = cmplx(EMC,0d0) + LH_PSI(N1) 
       PQQA = (3d0,0d0) - 4d0* S1 + 2d0/(N * N1) 
       PQGA = 4d0* (NS + N + 2d0) / (N * N1 * N2) 
       PGQA = 2d0 * (NS + N + 2d0) / (N * N1 * NM) 
       PGGA = 11d0/3D0 - 4d0* S1 + 4d0/(N * NM) + 4d0/(N1 * N2) 
       PGGB = - 4d0/3D0 
!                                                                       
!     Output to the array                                               
       P0NS      = CF * PQQA 
!                                                                       
       P0SG(1,1) = CF * PQQA 
       P0SG(1,2) = TR * dble(NF) * PQGA 
       P0SG(2,1) = CF * PGQA 
       P0SG(2,2) = CA * PGGA + TR *dble(NF) * PGGB 
!                                                                       
       RETURN 
      END                                           
!                                                                       
!     File: andim_nlo.f                                                 
!     Returns the NLO anomalous dimensions in the                       
!                                                                       
      SUBROUTINE LH_ANDIM_NLO(N,NF,P1NS,P1SG) 
      IMPLICIT none 
!                                                                       
      REAL*8 ca,cf,tr 
      COMMON/NNPDF10COLFACT/ca,cf,tr 
      REAL*8 emc,zeta2,zeta3,zeta4,pi 
      PARAMETER (pi    = 3.1415926535897932385) 
      COMMON /NNPDF10CONSTS/ emc,zeta2,zeta3,zeta4 
!                                                                       
      COMPLEX*16 LH_DPSI,LH_PSI, S1, S2 
      INTEGER I,J 
      COMPLEX*16 NS,NT,NFO,NFI,NSI,NSE,NE,NN 
      COMPLEX*16 N1,N2,NM,NMS,N1S,N1T,N2S,N2T 
      COMPLEX*16 N3,N4,N5,N6 
      COMPLEX*16 S11,S12,S13,S14,S15,S16 
      COMPLEX*16 SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP 
      COMPLEX*16 SSTR2P,SSTR3P 
      COMPLEX*16 PPSA,PQGA,PGQA,PGGA,PQGB,PGQB,PGGB,PGQC,PGGC 
      COMPLEX*16 PNPA,PNMA,PNSB,PNSC 
      COMPLEX*16 N 
      INTEGER NF 
      COMPLEX*16 P1NS(3) 
      COMPLEX*16 P1SG(2,2) 
!                                                                       
      S1 = EMC + LH_PSI(N+1d0) 
      S2 = ZETA2 - LH_DPSI(N+1d0,1) 
!                                                                       
      NS = N * N 
      NT = NS * N 
      NFO = NT * N 
      NFI = NFO * N 
      NSI = NFI * N 
      NSE = NSI * N 
      NE = NSE * N 
      NN = NE * N 
!                                                                       
      NM = N - 1d0 
      N1 = N + 1d0 
      N2 = N + 2d0 
      NMS = NM * NM 
      N1S = N1 * N1 
      N1T = N1S * N1 
      N2S = N2 * N2 
      N2T = N2S * N2 
                                                                        
! ..Analytic continuations of the occuring sums as given in GRV (1990)  
!   (with an improved parametrization of the moments of  Sp(x)/(1+x).)  
!                                                                       
       N3 = N + 3D0 
       N4 = N + 4D0 
       N5 = N + 5D0 
       N6 = N + 6D0 
       S11 = S1  + 1D0/N1 
       S12 = S11 + 1D0/N2 
       S13 = S12 + 1D0/N3 
       S14 = S13 + 1D0/N4 
       S15 = S14 + 1D0/N5 
       S16 = S15 + 1D0/N6 
       SPMOM = 1.0000D0 * (ZETA2 - S1 / N ) / N  -                      &
     &         0.9992D0 * (ZETA2 - S11/ N1) / N1 +                      &
     &         0.9851D0 * (ZETA2 - S12/ N2) / N2 -                      &
     &         0.9005D0 * (ZETA2 - S13/ N3) / N3 +                      &
     &         0.6621D0 * (ZETA2 - S14/ N4) / N4 -                      &
     &         0.3174D0 * (ZETA2 - S15/ N5) / N5 +                      &
     &         0.0699D0 * (ZETA2 - S16/ N6) / N6                        
!                                                                       
       SLC = - 5D0/8D0 * ZETA3 
       SLV = - ZETA2/2D0* (LH_PSI(N1/2D0) - LH_PSI(N/2D0))              &
     & + S1/NS + SPMOM                                                  
       SSCHLM = SLC - SLV 
       SSTR2M = ZETA2 - LH_DPSI (N1/2D0,1) 
       SSTR3M = 0.5D0 * LH_DPSI (N1/2D0,2) + ZETA3 
                                                                        
       SSCHLP = SLC + SLV 
       SSTR2P = ZETA2 - LH_DPSI(N2/2D0,1) 
       SSTR3P = 0.5D0 * LH_DPSI(N2/2D0,2) + ZETA3 
                                                                        
!                                                                       
! ..The contributions to P1NS as given in Gonzalez-Arroyo et al. (1979) 
!   (Note that the anomalous dimensions in the literature often differ  
!    from these moments of the splitting functions by factors -1 or -2, 
!    in addition to possible different normalizations of the coupling)  
!                                                                       
                                                                        
       PNMA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +                &
     &      16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2M ) +          &
     &      64D0* SSCHLM + 24D0* S2 - 3D0 - 8D0* SSTR3M -               &
     &      8D0* (3D0* NT + NS -1D0) / (NT * N1T) +                     &
     &      16D0* (2D0* NS + 2D0* N +1D0) / (NT * N1T) ) * (-0.5D0)     
       PNPA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +                &
     &      16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2P ) +          &
     &      64D0* SSCHLP + 24D0* S2 - 3D0 - 8D0* SSTR3P -               &
     &      8D0* (3D0* NT + NS -1D0) / (NT * N1T) -                     &
     &      16D0* (2D0* NS + 2D0* N +1D0)/(NT * N1T) ) * (-0.5D0)       
                                                                        
       PNSB = ( S1 * (536D0/9D0 + 8D0* (2D0* N + 1D0) / (NS * N1S)) -   &
     &      (16D0* S1 + 52D0/3D0- 8D0/(N * N1)) * S2 - 43D0/6D0 -       &
     &      (151D0* NFO + 263D0* NT + 97D0* NS + 3D0* N + 9D0) *        &
     &      4D0/ (9D0* NT * N1T) ) * (-0.5D0)                           
       PNSC = ( -160D0/9D0* S1 + 32D0/3.* S2 + 4D0/3D0 +                &
     &      16D0*(11D0*NS+5D0*N-3D0)/(9D0* NS * N1S))*(-0.5D0)          
!                                                                       
!     ..The contributions to P1SG as given in Floratos et al. (1981)    
!     ..Pure singlet (PS) and QG                                        
!                                                                       
       PPSA = (5d0* NFI + 32d0* NFO + 49d0* NT+38d0* NS + 28d0* N + 8d0)&
     &      / (NM * NT * N1T * N2S) * 2d0                               
!                                                                       
       PQGA = (-2d0* S1 * S1 + 2d0* S2 - 2d0* SSTR2P)                   &
     &      * (NS + N + 2d0) / (N * N1 * N2)                            &
     &      + (8d0* S1 * (2d0* N + 3d0)) / (N1S * N2S)                  &
     &      + 2d0* (NN + 6d0* NE + 15d0* NSE + 25d0* NSI + 36d0* NFI    &
     &      + 85d0* NFO + 128d0* NT + 104d0* NS + 64d0* N + 16d0)       &
     &      / (NM * NT * N1T * N2T)                                     
       PQGB = (2d0* S1 * S1 - 2d0* S2 + 5d0) * (NS + N + 2d0)           &
     &      / (N * N1 * N2) - 4d0* S1 / NS                              &
     &      + (11d0* NFO + 26d0* NT + 15d0* NS + 8d0* N + 4d0)          &
     &      / (NT * N1T * N2)                                           
                                                                        
!                                                                       
!     ..GQ and GG                                                       
       PGQA = (- S1 * S1 + 5d0* S1 - S2) * (NS + N + 2d0)               &
     &      / (NM * N * N1)  -  2d0* S1 / N1S                           &
     &      - (12d0* NSI + 30d0* NFI + 43d0* NFO + 28d0* NT - NS        &
     &      - 12d0* N - 4d0) / (2d0* NM * NT * N1T)                     
       PGQB = (S1*S1 + S2 - SSTR2P) * (NS + N + 2d0) / (NM * N * N1)    &
     &      - S1 * (17d0* NFO + 41d0* NS - 22d0* N - 12d0)              &
     &      / (3d0* NMS * NS * N1)                                      &
     &      + (109d0* NN + 621d0* NE + 1400d0* NSE + 1678d0* NSI        &
     &      + 695d0* NFI - 1031d0* NFO - 1304d0* NT - 152d0* NS         &
     &      + 432d0* N + 144d0) / (9d0* NMS * NT * N1T * N2S)           
       PGQC = (S1 - 8d0/3d0) * (NS + N + 2d0) / (NM * N * N1) + 1d0/ N1S 
       PGQC = 4d0/3d0* PGQC 
!                                                                       
       PGGA = - (2d0* NFI + 5d0* NFO + 8d0* NT + 7d0* NS- 2d0* N - 2d0) &
     &      * 8d0* S1 / (NMS * NS * N1S * N2S) -  67d0/9d0* S1 + 8d0/3d0&
     &      - 4d0* SSTR2P * (NS + N + 1d0) / (NM * N * N1 * N2)         &
     &      + 2d0* S1 * SSTR2P - 4d0* SSCHLP + 0.5d0 * SSTR3P           &
     &      + (457d0* NN + 2742d0* NE + 6040d0* NSE + 6098d0* NSI       &
     &      + 1567d0* NFI - 2344d0* NFO - 1632d0* NT + 560d0* NS        &
     &      + 1488d0* N + 576d0) / (18d0* NMS * NT * N1T * N2T)         
       PGGB = (38d0* NFO + 76d0* NT + 94d0* NS + 56d0* N + 12d0) *(-2d0)&
     &      / (9d0* NM * NS * N1S * N2)  +  20d0/9d0* S1  -  4d0/3d0    
       PGGC = (2d0* NSI + 4d0* NFI + NFO - 10d0* NT - 5d0* NS - 4d0* N  &
     &      - 4d0) * (-2d0) / (NM * NT * N1T * N2)  -  1d0              
!                                                                       
!     Output to the array                                               
       DO I=1,2 
          DO J=1,2 
             P1SG(I,J)=0d0 
          ENDDO 
       ENDDO 
                                                                        
       do I=1,3 
          P1NS(I)=0d0 
       enddo 
!                                                                       
!     NON SINGLET                                                       
!                                                                       
!     Plus                                                              
       P1NS(1) = CF *((CF-CA/2d0)* PNPA + CA* PNSB + TR*dble(NF)* PNSC) 
                                                                        
!     Minus=Valence                                                     
       P1NS(2) = CF *((CF-CA/2d0)* PNMA + CA* PNSB + TR*dble(NF)* PNSC) 
       P1NS(3) = P1NS(2) 
!                                                                       
!     SINGLET                                                           
       P1SG(1,1) = P1NS(1) + TR*dble(NF)*CF*PPSA*4d0 
       P1SG(1,2) = TR*dble(NF) * (CA * PQGA + CF * PQGB)*4d0 
       P1SG(2,1) = (CF*CF*PGQA + CF*CA*PGQB+TR*dble(NF)*CF*PGQC)*4d0 
       P1SG(2,2) = (CA*CA*PGGA + TR*dble(NF)*(CA*PGGB+CF*PGGC))*4d0 
!                                                                       
       RETURN 
      END                                           
!*****                                                                  
! 10 *------------------------------------------------------------------
!*****                                                                  
!                                                                       
!      matutils.f                                                       
!      A collection of Fortran utilities to deal with matrices:         
!        - multiplication (mmult_c, mmult_r).                           
!        - equating (mequal_c,mequal_r)                                 
!                                                                       
      SUBROUTINE LH_MMULT_C(A,ROWSA,COLSA,B,ROWSB,COLSB,C) 
      IMPLICIT none 
!                                                                       
      INTEGER ROWSA,COLSA,ROWSB,COLSB 
      COMPLEX*16 A(ROWSA,COLSA), B(ROWSB,COLSB), C(ROWSA,COLSB) 
      INTEGER I,J,K 
!                                                                       
!     Check that the matrices we would like to multiply have the        
!     correct dimensions                                                
!                                                                       
      IF (COLSA .NE. ROWSB) THEN 
          WRITE(*,*) 'Coglione... dimensioni delle matrici sbagliate!' 
          RETURN 
      ELSE 
!                                                                       
!     Initialize the output matrix to zero                              
!                                                                       
        DO I = 1,ROWSA 
          DO J = 1,COLSB 
            C(I,J) = (0.0d0,0.0d0) 
          ENDDO 
        ENDDO 
!                                                                       
!     Perform the multiplication according to:                          
!                          C[i,j] = Sum_k(A[i,k]*B[k,i])                
!                                                                       
        DO I = 1,ROWSA 
          DO J = 1,COLSB 
            DO K = 1,COLSA 
              C(I,J) = C(I,J) + A(I,K) * B(K,J) 
            ENDDO 
          ENDDO 
        ENDDO 
!                                                                       
      ENDIF 
      RETURN 
      END                                           
!                                                                       
!                                                                       
      SUBROUTINE LH_MMULT_R(A,ROWSA,COLSA,B,ROWSB,COLSB,C) 
      IMPLICIT none 
!                                                                       
      INTEGER ROWSA,COLSA,ROWSB,COLSB 
!                                                                       
      REAL*8 A(ROWSA,COLSA), B(ROWSB,COLSB), C(ROWSA,COLSB) 
!                                                                       
      INTEGER I,J,K 
!                                                                       
!     Check that the matrices we would like to multiply have the        
!     correct dimensions                                                
      IF (COLSA .NE. ROWSB) THEN 
          WRITE(*,*) 'wrong dimensions in matrix' 
          RETURN 
      ELSE 
!                                                                       
!     Initialize the output matrix to zero                              
        DO I = 1,ROWSA 
          DO J = 1,COLSB 
            C(I,J) = 0.0d0 
          ENDDO 
        ENDDO 
                                                                        
!                                                                       
!     Perform the multiplication according to:                          
!                          C[i,j] = Sum_k(A[i,k]*B[k,i])                
        DO I = 1,ROWSA 
          DO J = 1,COLSB 
            DO K = 1,COLSA 
              C(I,J) = C(I,J) + A(I,K) * B(K,J) 
            ENDDO 
          ENDDO 
        ENDDO 
!                                                                       
      ENDIF 
      RETURN 
      END                                           
!                                                                       
!                                                                       
      SUBROUTINE LH_MEQUAL_C(A,B,I,J) 
      IMPLICIT none 
!                                                                       
      INTEGER I,J 
      INTEGER K,L 
      COMPLEX*16 A(I,J), B(I,J) 
!                                                                       
      DO K=1,I 
         DO L=1,J 
            A(K,L) = B(K,L) 
         ENDDO 
      ENDDO 
      RETURN 
      END                                           
!                                                                       
!                                                                       
      SUBROUTINE LH_MEQUAL_R(A,B,I,J) 
      IMPLICIT none 
!                                                                       
      INTEGER I,J 
      INTEGER K,L 
      REAL*8 A(I,J), B(I,J) 
!                                                                       
      DO K=1,I 
         DO L=1,J 
            A(K,L) = B(K,L) 
         ENDDO 
      ENDDO 
      RETURN 
      END                                           
!                                                                       
!     psi.f                                                             
!     CalculATE complex psi function,  PZI(Z),  and its m'th derivatives
!     DPZI(Z,M) from the asymtotic expansions. The                      
!     functional equations are used for  |Im(Z)| < 10  to shift the     
!     argument to  Re(Z) >= 10  before applying the expansions.         
!                                                                       
       FUNCTION LH_PSI (Z) 
       IMPLICIT COMPLEX*16 (A-Z) 
       SUB = 0D0 
       ZZ = Z 
!                                                                       
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN 
!                                                                       
    1    CONTINUE 
         IF (DBLE(ZZ) .LT. 10.D0) THEN 
           SUB = SUB - 1d0/ ZZ 
           ZZ = ZZ + 1d0 
           GOTO 1 
         END IF 
!                                                                       
       END IF 
!                                                                       
! ..Use of the asymtotic expansion (at the shifted argument)            
!     Abramowitz-Stegun eq. 6.3.18 plus one term                        
!                                                                       
       RZ = 1d0/ ZZ 
       DZ = RZ * RZ 
       LH_PSI = SUB + LOG(ZZ) - 0.5d0 * RZ - DZ/5040D0 * ( 420d0+ DZ *  &
     &       ( - 42d0 + DZ * (20d0 - 21d0 * DZ) ) )                     
!                                                                       
       RETURN 
      END                                           
!                                                                       
!                                                                       
       FUNCTION LH_DPSI (Z,M) 
       IMPLICIT COMPLEX*16 (A-Z) 
!                                                                       
       INTEGER M, K1, K2 
       SUB = 0D0 
       ZZ = Z 
!                                                                       
! ..Shift of the argument using the functional equations                
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN 
!                                                                       
    1    CONTINUE 
         SUBM = -1d0/ZZ 
         DO 10 K1 = 1, M 
           SUBM = - SUBM * K1 / ZZ 
   10    CONTINUE 
!                                                                       
         IF (DBLE(ZZ) .LT. 10.D0) THEN 
           SUB = SUB + SUBM 
           ZZ = ZZ + 1d0 
           GOTO 1 
         END IF 
!                                                                       
       END IF 
!                                                                       
! ..Expansion (Bernoulli) coefficients for the first derivative         
       A1 =  1D0 
       A2 =  1d0/2D0 
       A3 =  1d0/6D0 
       A4 = -1d0/30D0 
       A5 =  1d0/42D0 
       A6 = -1d0/30D0 
       A7 =  5d0/66D0 
!                                                                       
! ..Expansion coefficients for the higher derivatives                   
       IF (M .EQ. 1) GOTO 2 
       DO 11 K2 = 2, M 
         A1 = A1 * (dble(K2)-1d0) 
         A2 = A2 *  dble(K2) 
         A3 = A3 * (dble(K2)+1d0) 
         A4 = A4 * (dble(K2)+3d0) 
         A5 = A5 * (dble(K2)+5d0) 
         A6 = A6 * (dble(K2)+7d0) 
         A7 = A7 * (dble(K2)+9d0) 
   11  CONTINUE 
    2  CONTINUE 
!                                                                       
! ..Use of the asymtotic expansion (at the shifted argument)            
!                                                                       
       RZ = 1d0/ ZZ 
       DZ = RZ * RZ 
       LH_DPSI = SUB + (-1)**(M+1d0) * RZ**M * ( A1 + RZ * (A2 + RZ *   &
     &        (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) )    
!                                                                       
       RETURN 
      END                                           
!                                                                       
!     dgauss_intern.f                                                   
!                                                                       
      DOUBLE PRECISION FUNCTION LH_DGAUSS_INTERN(F,A,B,EPS) 
      DOUBLE PRECISION W(12),X(12),A,B,EPS,DELTA,CONST,AA,BB,Y,C1,C2,S8,&
     &                 S16,U,F                                          
      DATA CONST /1.0D-25/ 
      DATA W                                                            &
     &       / 0.1012285362903762591525313543D0,		   &
     &         0.2223810344533744705443559944D0,		   &
     &         0.3137066458778872873379622020D0,		   &
     &         0.3626837833783619829651504493D0,		   &
     &         0.0271524594117540948517805725D0,		   &
     &         0.0622535239386478928628438370D0,		   &
     &         0.0951585116824927848099251076D0,		   &
     &         0.1246289712555338720524762822D0,		   &
     &         0.1495959888165767320815017305D0,		   &
     &         0.1691565193950025381893120790D0,		   &
     &         0.1826034150449235888667636680D0,		   &
     &         0.1894506104550684962853967232D0 /		   
      DATA X                         	        		   &
     &       / 0.9602898564975362316835608686D0,		   &
     &         0.7966664774136267395915539365D0,		   &
     &         0.5255324099163289858177390492D0,		   &
     &         0.1834346424956498049394761424D0,		   &
     &         0.9894009349916499325961541735D0,		   &
     &         0.9445750230732325760779884155D0,		   &
     &         0.8656312023878317438804678977D0,		   &
     &         0.7554044083550030338951011948D0,		   &
     &         0.6178762444026437484466717640D0,		   &
     &         0.4580167776572273863424194430D0,		   &
     &         0.2816035507792589132304605015D0,		   &
     &         0.0950125098376374401853193354D0 /		   
      DELTA=CONST*DABS(A-B) 
      LH_DGAUSS_INTERN=0. 
      AA=A 
    5 Y=B-AA 
      IF(DABS(Y) .LE. DELTA) RETURN 
    2 BB=AA+Y 
      C1=0.5D0*(AA+BB) 
      C2=C1-AA 
      S8=0. 
      S16=0. 
      DO 1 I = 1,4 
      U=X(I)*C2 
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U)) 
      DO 3 I = 5,12 
      U=X(I)*C2 
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U)) 
      S8=S8*C2 
      S16=S16*C2 
      IF(DABS(S16-S8) .GT. EPS*(1.0D0+DABS(S16))) GOTO 4 
      LH_DGAUSS_INTERN=LH_DGAUSS_INTERN+S16 
      AA=BB 
      GOTO 5 
    4 Y=0.5D0*Y 
      IF(DABS(Y) .GT. DELTA) GOTO 2 
      PRINT 7 
      LH_DGAUSS_INTERN=0. 
      RETURN 
    7 FORMAT(1X,'DGAUSS_INTERN ... TOO HIGH ACCURACY REQUIRED') 
      END                                           
