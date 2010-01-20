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
C Module: continuum.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
      subroutine continuumcdf(energy,flavor,ndaug,kf,px,py,pz,e)

      implicit none

      common/lyjets/n,k(4000,5),p(4000,5),v(4000,5)

      integer n,k
      real p,v

      real*8 pxsum,pysum,pzsum

      real*8 energy

      real senergy

      integer lycomp
      external lycomp

      integer flavor,kf(100),i,ndaug,istart
      real*8  px(100),py(100),pz(100),e(100)

      senergy=energy

C      if(flavor.ne.0) kf_user=flavor
C      flavor = kf_user


C     call lyeevt(kf_user,senergy)  ! e+e- --> qq shower parton model
 4    call lyeevt(flavor,senergy) ! e+e- --> qq shower parton model

C      call lylist(1)

      ndaug=0

C--- Save the primary quarks 

      do 1,i=1,n
        if(abs(k(i,2)).gt.6) goto 1
        if(k(i,2).eq.0)      goto 1
        if(k(i,3).ne.0)      goto 1
        ndaug=ndaug+1
        kf(ndaug)=k(i,2)
        px(ndaug)=p(i,1)
        py(ndaug)=p(i,2)
        pz(ndaug)=p(i,3)
        e(ndaug) =p(i,4)
1     continue

C     sum to check that we preserve momentum
      pxsum=0.0
      pysum=0.0
      pzsum=0.0

      istart=0

      do 2, i=1,n
         if(k(i,3).le.0.and.k(i,2).ne.22) goto 2
c        if particle is jet or cluster

         if (k(i,2).eq.22.and.istart.eq.0) goto 22
         if ( k(i,3).le.0) goto 2

         if(k(k(i,3),2).eq.92.or.k(k(i,3),2).eq.91) goto 22

         goto 2
         
 22      continue
C     Copy particle. If parent is cluster or string or photon if
C     before first cluster or string. Uggly...
         if (k(i,2).ne.22) istart=1
         ndaug = ndaug + 1
         kf(ndaug)=k(i,2)
         px(ndaug)=p(i,1)
         py(ndaug)=p(i,2)
         pz(ndaug)=p(i,3)
         pxsum=pxsum+px(ndaug)
         pysum=pysum+py(ndaug)
         pzsum=pzsum+pz(ndaug)
         e(ndaug) =p(i,4)
 2    continue


      if (abs(pxsum).gt.0.001.or.
     +       abs(pysum).gt.0.001.or.
     +       abs(pzsum).gt.0.001) then


         print *, 'Momentum not conserved in jetset fragmentation:'
         print *,'dPx:',pxsum,' dPy:',pysum,' dPz:',pzsum

         call lylist(1)

         goto 4

      endif

      return
      end






