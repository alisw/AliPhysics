*CMZ :          11/11/94  15.39.18  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   13/07/94
**DOC
C  Parallel Geant - main routines, to be used in parallelising GEANT programs.
C  Target architectures: scalable parallel processors and (for appropriate
C  problems) workstation clusters
C  Implemenation relies on availability of:
C  An implementation of the MPI (Message Passing Interface Standard) for the
C   relevant ensemble of machines or parallel processor.
C
C  A common file base for input and output files, such as NFS or AFS.
C
C  [ NOTE: MPI Implementations exist for most Unix machines, including
C    the portable, public domain, implementation. For more information see
C    the MPI Web page at URL http://www.mcs.anl.gov/mpi/index.html  ]
C  ---
C  First released version  March 1996, John Apostolakis japost@dxcern.cern.ch
C  --------------------------------------------------------------------------
C     subroutine gprun
C
C           Function: in some cases of static event partitioning, it
C             makes sure each node only does its portion of the events to
C             to be simulated.
C
C     called by:  grun (if CERNLIB_PARA switch is used)
C------------------------------------------------------------------------
**ENDDOC
#if defined(CERNLIB_PARA)
      subroutine gprun
c     Routine is called at the start of grun, to normalize Number of Events
c       in the case of a "static" decomposition of events
      implicit none
*
#include "geant321/gcflag.inc"
#include "geant321/multiprox.inc"
*
*
      integer modsize
*
c     The default method of partitioning events between processors is
c     the Static method.

c     To choose the dynamic method, in which events are shared out 
c     to different nodes by a master, the user must call gpdynamic
c     and write a multiplexing routine (modifying the routine muxread 
c     in the example gexam3).


      if ( npleader .eq. -1) then
c
c     Static distribution of events, the new default
c     All nodes share the events equally. It will work well if events take
c     close to the same CPU time to simulate, ie the standard deviation
c     of the event time is small compared to the average event time.
c
c     In this case do a static distribution of events ...
c     (assumes 1<=nprank<npsize, so it is OK for MPI, where 0<=nprank<npsize )
c
c
	  nevtot  = nevent
	  nevent  = nevtot / npsize
	  modsize = mod(nevtot,npsize)
	  if( nprank .lt. modsize ) then
	      nevent = nevent + 1
	  endif
      else
c
c     Dynamic distribution of events, the new default
c
c     This version does dynamic distribution of events, and terminates at
c     the end of the input file. The number of events each node will
c     simulate is not determined in advance. Each node requests work
c     from the "master" node, which gives it out (the "farming" model).
c     The number of events that each node will simulate can vary, so it
c     is best to make the maximum number of events
c     per node equal be the original nevent [ ... or else one node could
c     terminate unexpectedly ]
c
	  nevtot  = nevent
c
c     This code must be used in conjunction with additional code in
c     the routine "gukine" to handle sharing the work load. An example
c     of such a routine can be found in "muxread".

      endif

      return
      end
#endif
