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
C Module: begevtgenstore.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
      subroutine begevtgenstore(nument,eventnum,numparticle,istat,
     +                          partnum,motherfirst,motherlast,
     +                          daugfirst,dauglast,
     +                          px,py,pz,e,m,x,y,z,t)
*
* routine to fill the stdhep common blocks from
* evtgen (C++). This routine allows the C++ program not to
* have to mess with common blocks.
*
* Anders Ryd,  Dec 96   Created.
*
* 
      implicit none
#include "EvtGenModels/stdhep.inc"

        integer nument,eventnum,numparticle,istat,partnum
        integer motherfirst,motherlast,daugfirst,dauglast

        integer stdhepnum

        double precision px,py,pz,e,m,x,y,z,t

        stdhepnum=partnum

        if (istat.eq.1) then

c       make the quarks unstable so that geant will not try
c       to put them through the detector simulation!

         if (
     +        stdhepnum.eq.-1.or.
     +        stdhepnum.eq. 1.or.
     +        stdhepnum.eq.-2.or.
     +        stdhepnum.eq. 2.or.
     +        stdhepnum.eq.-3.or.
     +        stdhepnum.eq. 3.or.
     +        stdhepnum.eq.-4.or.
     +        stdhepnum.eq. 4.or.
     +        stdhepnum.eq.-5.or.
     +        stdhepnum.eq. 5.or.
     +        stdhepnum.eq.-6.or.
     +        stdhepnum.eq. 6) then
           istat=2
         endif




c         if (
c     +        stdhepnum.ne.-1.and.
c     +        stdhepnum.ne. 1.and.
c     +        stdhepnum.ne.-2.and.
c     +        stdhepnum.ne. 2.and.
c     +        stdhepnum.ne.-3.and.
c     +        stdhepnum.ne. 3.and.
c     +        stdhepnum.ne.-4.and.
c     +        stdhepnum.ne. 4.and.
c     +        stdhepnum.ne.-5.and.
c     +        stdhepnum.ne. 5.and.
c     +        stdhepnum.ne.-6.and.
c     +        stdhepnum.ne. 6.and.
c     +        stdhepnum.ne.-11.and.
c     +        stdhepnum.ne.11.and.
c     +        stdhepnum.ne.-13.and.
c     +        stdhepnum.ne.13.and.
c     +        stdhepnum.ne.12.and.
c     +        stdhepnum.ne.-12.and.
c     +        stdhepnum.ne.14.and.
c     +        stdhepnum.ne.-14.and.
c     +        stdhepnum.ne.16.and.
c     +        stdhepnum.ne.-16.and.
c     +        stdhepnum.ne.2212.and.
c     +        stdhepnum.ne.-2212.and.
c     +        stdhepnum.ne.2112.and.
c     +        stdhepnum.ne.-2112.and.
c     +        stdhepnum.ne.211.and.
c     +        stdhepnum.ne.-211.and.
c     +        stdhepnum.ne.321.and.
c     +        stdhepnum.ne.-321.and.
c     +        stdhepnum.ne.130.and.
c     +        stdhepnum.ne.310.and.
c     +        stdhepnum.ne.22) then
c    
c            print *,'particle stdhep:',stdhepnum,' is stable'
c
c         endif

        endif

        nevhep=eventnum
        nhep=numparticle
        isthep(nument)=istat
        idhep(nument)=stdhepnum
        jmohep(1,nument)=motherfirst
        jmohep(2,nument)=motherlast
        jdahep(1,nument)=daugfirst
        jdahep(2,nument)=dauglast
        phep(1,nument)=px
        phep(2,nument)=py
        phep(3,nument)=pz
        phep(4,nument)=e
        phep(5,nument)=m
        vhep(1,nument)=x
        vhep(2,nument)=y
        vhep(3,nument)=z
        vhep(4,nument)=t

 	 return

      end





