*CMZ :          11/11/94  14.40.52  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   13/07/94
C------------------------------------------------------------------------
C     gpinit
C           Function: initialises the parallel environment (MPI) and
C           parallel's Geant's variables that track it.
C
C     called by:  <USER>, uginit, (if CERNLIB_PARA switch is used)
C------------------------------------------------------------------------
#if defined(CERNLIB_PARA)
      subroutine gpinit
      implicit none
#include "geant321/mpifinc.inc"
#include "geant321/multiprox.inc"
#include "geant321/multseeds.inc"
#include "geant321/gcflag.inc"
      integer             msgerrpas, ierr
      common  /gpinitpas/ msgerrpas, ierr

      integer ierr1, ierr2, ierr3
      integer ibase, lunerr
      integer msgerror, rngerror
      data    ibase / 1 /
      data    msgerror, rngerror / 0, 0 /
      parameter (lunerr= 0)

c     Initialise Event Level Parallel Geant (message passing)
c
c     Tasks:
c       Initialise Message Passing and my relevant variables (npsize etc)

      nprank= -1

c     Initialise Message Passing

      call MPI_init( ierr1 )
      call MPI_Comm_rank( MPI_COMM_WORLD, nprank, ierr2 )
      call MPI_Comm_size( MPI_COMM_WORLD, npsize, ierr3 )
      if ( (ierr1 .gt. 0)  .or. (ierr2 .gt. 0)  .or. (ierr3 .gt. 0) )
     +    msgerror=1

      write(*,*)  ' Node ', nprank , ' initialised - one of ', npsize

c     Set Error flag ...

      if( msgerror .ne. 0) then
          write (lunerr, *)  ' Message Passing Library error ',
     &				        ' in initialisation. '
          if ( (ierr1.ne.0) .or. (nprank .ne. -1)) then
             write (lunerr, *)  ' Occured in node ', nprank
          else
             write (lunerr, *)  ' MASSIVE ERROR> in an unknown node. '
             write (lunerr, *)   
     &       ' This is non-recoverable, MPI did not initialise at all'
c
c            This is non-recoverable, MPI did not initialise at all
             call gpabort()
c            ieorun= 1
             stop
          endif
      endif
c
c     By default a static decomposition of events is done, with
c       each processor geting 1/number of processors of the share 
c       of events. 
c     If an input file exists or the time per event has a large standard
c       deviation, the events can be dynamically assigned. You will then
c       need to look at the routines gukine and muxread in gexam3.
c
c     The following just initialises the default "static" configuration,
c       in which there is no "master" or "leader".
c
c         -- Id of task expected to be a "master"
      npleader = -1
c         -- Id of task expected to be 1st worker
c       		          (representative for printout ...)
      nfirstworker= 0

      msgerrpas= msgerror

      return
      end
#else
      SUBROUTINE GPINIT_DUMMY
      END
#endif
