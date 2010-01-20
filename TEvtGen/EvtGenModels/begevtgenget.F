      subroutine begevtgenget(entry,eventnum,numparticle,istat,
     +                          partnum,mother,daugfirst,dauglast,
     +                          px,py,pz,e,m,x,y,z,t)
*
* Copyright Information: See EvtGen/COPYRIGHT
*
* routine to read from the stdhep common block from
* evtgen (C++). This routine allows the C++ program not to
* have to mess with common blocks.
*
* Anders Ryd,  Oct 1   Created.
*
* 
      implicit none
#include "EvtGenModels/stdhep.inc"

        integer entry,eventnum,numparticle,istat,partnum
        integer mother,daugfirst,dauglast

        integer stdhepnum

        double precision px,py,pz,e,m,x,y,z,t

        stdhepnum=partnum



        eventnum=nevhep
        numparticle=nhep
        istat=isthep(entry)
        stdhepnum=idhep(entry)
        mother=jmohep(1,entry)
        daugfirst=jdahep(1,entry)
        dauglast=jdahep(2,entry)
        px=phep(1,entry)
        py=phep(2,entry)
        pz=phep(3,entry)
        e=phep(4,entry)
        m=phep(5,entry)
        x=vhep(1,entry)
        y=vhep(2,entry)
        z=vhep(3,entry)
        t=vhep(4,entry)

 	 return

      end








