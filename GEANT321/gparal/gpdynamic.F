*CMZ :          11/11/94  14.40.52  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   13/07/94
C------------------------------------------------------------------------
C     gpdynamic
C           Function: informs parallel Geant that events will be 
C                     dynamically distributed between nodes, ie
C                     each "worker" node will request (the data of) an event
C                     from a "master" or "leader". 
C
C     called by:  <USER>, (if CERNLIB_PARA switch is used)
C------------------------------------------------------------------------
#if defined(CERNLIB_PARA)
      subroutine gpdynamic
      implicit none
#include "geant321/multiprox.inc"

c
c     This routine overrides the default (set in gpinit) 
c       of a static decomposition of events , with which 
c       each processor would get 1/(number of processors) of the share 
c       of events. 
c     It should be used if an input file exists or the time per event 
c       has a large standard deviation, the events can be dynamically assigned.
c       You will need to look at the routines gukine and muxread
c       in gexam3 for how to distribute events.

c
c     The following initialises the "dynamic" configuration,
c       in which there is no "master" or "leader".
c
c         -- Id of task expected to be a "master"
      npleader = 0
c         -- Id of task expected to be 1st worker
c       		          (representative for printout ...)
      nfirstworker= 1

      return
      end
#else
      SUBROUTINE GPDYNAMIC_DUMMY
      END
#endif
