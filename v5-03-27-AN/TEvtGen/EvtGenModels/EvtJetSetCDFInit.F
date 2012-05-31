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
C Module: EvtJetSetInit.F
C
C Description:
C
C Modification history:
C
C    DJL/RYD     August 11, 1998         Module created
C
C------------------------------------------------------------------------
      subroutine evtjetsetcdfinit(fname)

      implicit none

      external lydata

      character*(*) fname


      open(54,STATUS='OLD',FILE=fname)

      call LYUPDA(2,54)

      close(54)      


C      call LYLIST(12)
      

      end






