      subroutine GetOrderPDF(order)
      implicit none
      integer nset,order
      real*8 Q,Q2
      nset = 1
      call GetOrderPDFM(nset,order)
      return
c
      entry GetRenFac(Q)
      nset = 1
      call GetRenFacM(nset,Q)
      return
c
      entry GetQ2fit(Q2)
      nset = 1
      call GetQ2fitM(nset,Q2)
      return
      end
c
      subroutine InitEvolve(nset)
      implicit none
      character*16 s2
      real*8 x,Q,Q2,Q2fit,muR
      integer nset,order,EvlOrd
c
      save EvlOrd,Q2fit,muR
*
      EvlOrd=-1
      read(1,*) s2,Q2fit,muR
      if (index(s2,'lo').eq.1) EvlOrd=0
      if (index(s2,'nlo').eq.1) EvlOrd=1
      if (index(s2,'nnlo').eq.1) EvlOrd=2
      if (EvlOrd.lt.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Unknown PDF evolution order ',s2
         stop
      endif
      if (muR.ne.1.0d0) then
         write(*,*) '***********************************************'
         write(*,*) '* Note than the renormalization scale is      *'
         write(*,*) '* unequal to the factorization scale for this *'
         write(*,*) '* particular PDF set.                         *'
         write(*,*) '* See manual for proper use.                  *'
         write(*,*) '***********************************************'
      endif
c      print *,'calling readevolve', nset
      call readevolve(nset)
      return
*
      entry InitEvolveCode(nset)
c      print *,'calling initevolution', nset,Evlord,Q2fit
      call initevolution(nset,EvlOrd,Q2fit)
      return
*
      entry GetOrderPDFM(nset,order)
      order=EvlOrd
      return
*
      entry GetRenFacM(nset,Q)
      Q=muR
      return
*
      entry GetQ2fitM(nset,Q2)
      Q2=Q2fit
      return
*
      end
