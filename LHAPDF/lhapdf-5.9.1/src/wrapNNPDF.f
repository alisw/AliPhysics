!********************************************************
!
!     wrapNNPDF.f: 
!     Routine called by LHAPDF package for the evolved PDF
!     in a (x,Q) point as called by NNPDF.LHpdf file.
!         
!     In 'wrapevolve.f' the package calls:
!     IF(NAME(NSET).EQ.'NNPDF') call NNPDFevolve(x,Q,f)
!     IF(NAME(NSET).EQ.'NNPDF') call NNPDFread(nset)
!     IF(NAME(NSET).EQ.'NNPDF') call NNPDFalfa(alfas,Q)
!     IF(NAME(NSET).EQ.'NNPDF') call NNPDFinit(nset,Eorder,Q2fit)
!
!
!********************************************************

      SUBROUTINE NNPDFevolve(x,Q,f)
      IMPLICIT none
!
      INCLUDE 'parmsetup.inc'
      CHARACTER*16 name(nmxset)
      INTEGER nmem(nmxset),ndef(nmxset),mmem
      COMMON/NAME/name,nmem,ndef,mmem
!
      INTEGER nset,pdfmem
      REAL*8 f(-6:6)
      REAL*8 x,Q
!
      INTEGER order
      REAL*8 alfas,alphaNNPDF,alfas0,q0
      REAL*8 Eorder,Q2fit,mass
      REAL*8 pdfout(13)
      REAL*8 check
      REAL*8 pi
      PARAMETER (pi = 3.1415926535897932385)
!
      INTEGER i,nff 
      REAL*8 xpt
!
      INTEGER ipt,imodev,ivfn,itmc
      COMMON/nnpdf10EVFLAGS/ipt,imodev,ivfn,itmc 
      INTEGER ieval,niter,nmax
      COMMON/nnpdf10ngrid/ ieval,niter,nmax
      REAL*8 xmin,xm1,xm2,xmax
      COMMON/nnpdf10GRID/xmin,xm1,xm2,xmax  
      REAL*8 epsevol
      COMMON/nnpdf10evolacc/epsevol
!
      REAL*8 qq,qq2,qth(4:6)
      REAL*8 q20,q2
      COMMON/nnpdf10EVSCALE/q20,q2
      REAL*8 q2th(4:6),asref,q2ref
      COMMON/nnpdf10vfns/q2th,asref,q2ref
      REAL*8 as0,asc,asb,ast,asq
      COMMON/nnpdf10AS/ as0,asc,asb,ast,asq
!
      q2    = q**2.                ! for the commons
      qq2   = q**2.
      qq = q
      asq   = alphaNNPDF(qq)   
      asq   = asq / (4.*pi) 
!       
      xpt = x
      CALL lh_evolfactx(0,xpt)
      CALL lh_pdfevolx(xpt,pdfout)
      CALL lh_pdfevln2lha(pdfout,f)
      DO i = -6,6,1
         f(i) = x * f(i)
      ENDDO
!
      RETURN

!********************************************************

      ENTRY NNPDFread(nset)
!
!     Read evol params from NNPDF08.LHpdf
      READ(1,*) ivfn,imodev,itmc,niter,xmin,xm1,xm2,xmax,epsevol
!
      RETURN

!********************************************************

      ENTRY NNPDFalfa(alfas,Q)
!      
      QQ = Q             
      alfas =  alphaNNPDF(QQ)
!
      RETURN

!********************************************************

      ENTRY NNPDFinit(nset,Eorder,Q2fit)
!      
      CALL GetOrderPDFM(nset,order)
      IPT = order 
!
      CALL GetQ2fitM(nset,Q2)
      Q2fit = Q2
      Q20   = Q2
!
      call GetQmassM(nset,4,mass)
      QTH(4) = mass
      call GetQmassM(nset,5,mass)
      QTH(5) = mass
      call GetQmassM(nset,6,mass)
      QTH(6) = mass
      DO i = 4,6
         q2th(i) = qth(i)**2d0
      ENDDO
!
      CALL GetAlfas(nset,alfas0,Q0)
      asref = 0.119 
      q2ref = q0**2d0
!
      RETURN

!********************************************************

      ENTRY NNPDFpdf(nset)
!
      pdfmem = mmem
      IF (pdfmem.LT.0) THEN
         WRITE(*,*) 'NNPDF set:'
         WRITE(*,*) 'PDF member out of range:'
         WRITE(*,*) 'member = ',pdfmem
         STOP
      ENDIF
      RETURN
!
      END
