      SUBROUTINE IFGFACT(IKEY,ITYPE,IINIT)
C stores and manages options for   wid_a1_fit
C ITYPE=0,1,2 means zero value, from full table, from GFACT plus table for KKpi
C IINIT is for GFACT, whether it  has to reinitialize itself or just be used
C IKEY 0,1,2 means  initialization  use by GFUNCT (reini) other use
C IKEY -1 mieans initialization with recalculation og GFUNCT constats

      DATA ITYPE0,IINIT0 /1,0/
      SAVE ITYPE0,IINIT0

      IF    (IKEY.EQ.-1) THEN
        ITYPE0=ITYPE
        IINIT0=IINIT
      ELSEIF(IKEY.EQ.0) THEN
        ITYPE0=ITYPE
      ELSEIF(IKEY.EQ.1) THEN
        IINIT=IINIT0
        IINIT0=1
      ELSEIF(IKEY.EQ.2) THEN
        ITYPE= ITYPE0
      ENDIF
      END

      real function wid_a1_fit(qq)
      implicit none
      real                     qq
C..............................................................
C.    Output:    a1 width  as function of qq (formula (29) of ref [1]) 
C.               but to speed up execution linear interpolation 
C.               from numerical table is used. From read in table 
C.               linear interpolation (extrapolation) is used.
C.    Input:     qq  [GeV**2]
C.               table wida1_qq_tot_table.txt prepared in 
C.               ./a1-tabler directory
C.    References:[1] arXiv:0911.4436 (hep-ph) D.Gomez Dumm et al.(tau ->3pi nu)
C.    Called    : files f3pi_rcht.f  fkkpi.f fkk0pi0.f
C                       ./demo-standalone/tautestroman.f 
C..............................................................

      integer Nq
      PARAMETER (Nq=1001)   
      integer ik
      common /a1_width/ qmax,qmin,qk_min,    w_qmin,    del_qq
      real*8            qmax,qmin,qk_min(Nq),w_qmin(Nq),del_qq
      integer kq
      character path*120 

      real*8 aq,bq,qk_max,w_qmax
      real*8 xqq,GFACT
      real wid_a1_fitKKpi


      INTEGER IFIRST,ITYPE
      DATA    IFIRST,ITYPE /0,1/
      INTEGER  IMODE,IINIT
!      DATA IMODE /0/


C. READING TABLE AND INITIALIZATION
C. ================================

      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C Initialization from automatically created routine
        call InitA1TAB(qk_min    ,w_qmin    )
        qmax = qk_min(Nq) 
        qmin = qk_min(1)  
        del_qq = (qmax-qmin)/float(Nq-1)
      ENDIF

      CALL IFGFACT(2,ITYPE,IINIT)
C ALTERNATIVE CALCULATIONS FOR SPECIAL PURPOSES
C. ============================================
      IF (ITYPE.EQ.0) THEN
       wid_a1_fit=0
       RETURN
      ELSEIF (ITYPE.EQ.2) THEN
       xqq=qq
       wid_a1_fit=wid_a1_fitKKpi(qq)+GFACT(xqq)
!       write(*,*) 'uwaga', wid_a1_fit, wid_a1_fitKKpi(qq)
       RETURN
      ENDIF


C. INTEPOLATION, for extrapolation values at ends are used
C. =======================================================
      if(qq.gt.qmin.and.qq.le.qmax) then
        kq = (qq-qmin)/del_qq
        kq = kq +1

        qk_max = qk_min(kq+1)
        w_qmax = w_qmin(kq+1)

        aq = (w_qmax-w_qmin(kq))/(qk_max-qk_min(kq))
        bq = (w_qmax*qk_min(kq) -w_qmin(kq)*qk_max)
     $      /(qk_min(kq)-qk_max)
        wid_a1_fit = aq*qq+bq

      elseif(qq.ge.qmax) then               ! above maximun
        wid_a1_fit = w_qmin(Nq) 

      else                                  ! below minimum     
        wid_a1_fit = w_qmin(1)

      endif

      return 
      end


