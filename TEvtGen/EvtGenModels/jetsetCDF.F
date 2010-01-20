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
C Module: jetset1.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
      subroutine jetsetcdf(ip,m,ndaug,kf,km,px,py,pz,e)

C
C interface to JETSET 7.4 to have one particle decayed
C including possibly fragmentation, if the decay products include
C partons.
C

      implicit none

      INTEGER MSTU,MSTJ
      REAL PARU,PARJ
      COMMON/LYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LYDAT1/

      common/lyjets/n,k(4000,5),p(4000,5),v(4000,5)

      integer n,k
      real p,v

      integer ndmax
      parameter (ndmax=100)

      integer lycomp
      external lycomp

      integer ip,kf(ndmax),i,ndaug,km(ndmax)
      integer kp,kid,ipart1,ipart
      real*8  m,px(ndmax),py(ndmax),pz(ndmax),e(ndmax)

      integer ip1
      real qmax

c      used to use ly1ent both since it does not set the mass
c      of the daughter we have to manipulate the common blocks our
c      self. ryd April 25-1999

c      call ly1ent(1,ip,0.0,0.0,0.0)

      K(1,1)=1
      K(1,2)=ip
      P(1,5)=m
      P(1,4)=m
      P(1,1)=0.0
      P(1,2)=0.0
      P(1,3)=0.0
      n=1


c     now we can decay this particle
      call lydecy(1)

C code copied from LUEXEC to avoid error with shower
C switched on
C...Decay products may develop a shower.
        IF(MSTJ(92).GT.0) THEN
          IP1=MSTJ(92)
          QMAX=SQRT(MAX(0.,(P(IP1,4)+P(IP1+1,4))**2-(P(IP1,1)+P(IP1+1,
     &    1))**2-(P(IP1,2)+P(IP1+1,2))**2-(P(IP1,3)+P(IP1+1,3))**2))
          CALL LYSHOW(IP1,IP1+1,QMAX)
          CALL LYPREP(IP1)
          MSTJ(92)=0
        ELSEIF(MSTJ(92).LT.0) THEN
          IP1=-MSTJ(92)
          CALL LYSHOW(IP1,-3,P(IP1,5))
          CALL LYPREP(IP1)
          MSTJ(92)=0
        ENDIF

c
c for debugging:
c      call lylist(1)
c
      mstj(21)=0
      call lyexec
      mstj(21)=2
 
c find partons, delete secondary partons, set mother pointers

      ndaug = 0
      ipart1 = 1
      ipart = 1

      do 10 i=2,n
        kp = k(i,3)
        kid = k(i,2)
        if (abs(kid) .ge. 1 .and. abs(kid) .le. 8
     1  .or. kid .eq. 21
     2  .or. kid .ge. 91 .and. kid .le. 94)  then
           if (ipart1 .eq. 1) ipart1 = i
           ipart = i
           if (kp .ne. 1)  goto 10
           kp = 0
        else
           if (kp .gt. ipart)  then
              goto 10
           elseif (kp .ge. ipart1)  then
              kp = ipart1-1
           else
              kp = 0
              endif
           endif
        ndaug = ndaug + 1
        km(ndaug)=kp
        kf(ndaug)=kid
        px(ndaug)=p(i,1)
        py(ndaug)=p(i,2)
        pz(ndaug)=p(i,3)
        e(ndaug)=p(i,4)
c
c for debugging:
c        print '( 2I5,I12,4F12.4 )',ndaug,km(ndaug),kf(ndaug),
c     1    px(ndaug),py(ndaug),pz(ndaug),e(ndaug)
c
   10 continue

      end


