      real function wid_a1_fitKKpi(qq)
      implicit none
      real                     qq
C..............................................................
C.    Output:    KKpi part ofa1 width  as function of qq 
C.               but to speed up execution linear interpolation 
C.               from numerical table is used. From read in table 
C.               linear interpolation (extrapolation) is used.
C.    Input:     qq  [GeV**2]
C..............................................................

      integer Nq
      PARAMETER (Nq=1001)   
      integer ik
      real*8 qmax,qmin,qk_min(Nq),w_qmin(Nq),del_qq
      integer kq
      character path*120 

      real*8 aq,bq,qk_max,w_qmax

      INTEGER IFIRST
      DATA IFIRST/0/
      INTEGER  IMODE
!      DATA IMODE /0/
      COMMON   /IMODE/  IMODE      ! 0 is for calculation of G , 1 is for spectra, 2 is for use of G

      IF (IMODE.EQ.0) THEN
       wid_a1_fitKKpi=0
      ENDIF

C. READING TABLE AND INITIALIZATION
C. ================================

      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C Initialization from automatically created routine
        call InitA1TABKKpi(qk_min    ,w_qmin    )
        qmax = qk_min(Nq) 
        qmin = qk_min(1)  
        del_qq = (qmax-qmin)/float(Nq-1)
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
        wid_a1_fitKKpi = aq*qq+bq

      elseif(qq.ge.qmax) then               ! above maximun
        wid_a1_fitKKpi = w_qmin(Nq) 

      else                                  ! below minimum     
        wid_a1_fitKKpi = w_qmin(1)

      endif

      return 
      end


