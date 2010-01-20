C--------------------------------------------------------------------------
C
C Environment:
C      This software is part of the EvtGen package developed jointly
C      for the BaBar and CLEO collaborations.  If you use all or part
C      of it, please give an appropriate acknowledgement.
C
C Copyright Information: See EvtGen/COPYRIGHT
C      Copyright (C) 1998      Caltech, UCSB
C
C Module: jetsetcont.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
      subroutine jetsetcont(energy,flavor,ndaug,kf,px,py,pz,e)



      implicit none

      common/lujets/n,k(4000,5),p(4000,5),v(4000,5)

      integer n,k
      real p,v

      real*8 energy

      real senergy

      integer lucomp
      external lucomp

      integer flavor,kf(100),i,ndaug
      real*8  px(100),py(100),pz(100),e(100)

      senergy=energy


      call lu2ent(0,flavor,-flavor,senergy)


      ndaug=0

      do 10,i=4,n

c     only wants the particles that are created from the
c     fragmentation, not the secondary daugthers.

        if (k(i,3).eq.3) then
          kf(i-3)=k(i,2)
          px(i-3)=p(i,1)
          py(i-3)=p(i,2)
          pz(i-3)=p(i,3)
          e(i-3)=p(i,4)
          ndaug=ndaug+1
         endif

10    continue


      end



