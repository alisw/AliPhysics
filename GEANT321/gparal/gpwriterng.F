*CMZ :          26/02/96  11.38.47  by  S.Ravndal
*-- Author :
C------------------------------------------------------------------------
C     gpwriterng( lunwrite, filename, ierr )
C
C     Function: write out the state of the pseudo random number
C               generator on all the nodes into a file.
C
c     A convention: if lunwrite is >= 0, use that unit.
c       If not, create a new file called <filename>
C
C     Possible error codes:
c	  ierr= -102         The number of seeds required for a subsequence of
c                              the RNG (numperseq) is more than allocated
c			       in the array iseedrng (SEEDS_PER_SEQ)
c
C     called by:  <USER> (if PARALLEL switch is used)
C------------------------------------------------------------------------
#if defined(CERNLIB_PARA)
      subroutine gpwriterng( lunwrite, filename, ierr )
      implicit none
      integer lunwrite, ierr
      character filename*(*)

#include "geant321/mpifinc.inc"
#include "geant321/multseeds.inc"
      integer    numprocs, nrank, nleader
      integer    MAXNODES, SEEDS_PER_SEQ
      parameter (MAXNODES=4096)
      character*12 chNameRng

c     In order to accomodate RANMAR, which has the largest table of seeds.
      parameter (SEEDS_PER_SEQ=103)

c     For  RANECU it would be enough: (SEEDS_PER_SEQ=2).
c     For ranlux            parameter (SEEDS_PER_SEQ= 25)

      integer   myseeds(SEEDS_PER_SEQ)
      integer   iseedrng( SEEDS_PER_SEQ, MAXNODES )
c     version of the seed file format
      real   version
      integer nsubseq, numperseq, nperline
      integer lunwr, lunerr
      integer isubseq, inumb, imax, ino
      data  version / 1.0 /

      if (lunwrite .lt. 0) then
c     Default number
          lunwr= 74   
          open (unit=lunwr, file=filename )
      else
          lunwr= lunwrite
      endif

c     Make sure to go to the start of the file
c
      rewind( lunwr )
#if defined(CERNLIB_PARA_RANECU)
c     -----------------------------------------------------------------
c     The RANECU pseudo random number generator, used by GEANT's GRNDM.
c     -----------------------------------------------------------------
#if !defined(CERNLIB_PARA_TASKRNG)
c     Currently the number of sequences is the number of processors!
      call gprocs( numprocs, nrank, nleader )
      nsubseq=   numprocs

      if( numprocs .gt. MAXNODES ) then
         write(*, *) 'GPwriteRng: ERROR Too many nodes ', numprocs
	 ierr= -101
	 goto 9999
      endif

#endif
      chnamerng= 'RANECU'
      numperseq= 2              !  2 "seed" integers per RANECU sequence
      nperline=  2              !  2 numbers per line
#endif
#if defined(CERNLIB_PARA_RANMAR)
      chnamerng= 'RANMAR'
      numperseq= 103            !  103 "seed" integers per RANMAR sequence
      nperline=  5              !    5 numbers per line
#endif
#if defined(CERNLIB_PARA_RANLUX)
      chnamerng= 'RANLUX'
      numperseq= 25             !   25 "seed" integers per RANMAR sequence
      nperline=  5              !    5 numbers per line
#endif

c     "Gather" all seeds to node 0
c
#if defined(CERNLIB_PARA_RANECU)
#if !defined(CERNLIB_PARA_TASKRNG)
c     In this version there is one sequence per node
      call grndmq( myseeds(1), myseeds(2), iseqnc, 'G' )

      if (numperseq .gt. SEEDS_PER_SEQ ) then
	  ierr= -102
	  goto 9999
      endif

c     The following has all nodes send their seed, and node 0 collect them.
c      It should work if  numperseq <= SEEDS_PER_SEQ
c
      call MPI_Gather( myseeds, numperseq, MPI_INTEGER, iseedrng,
     $     SEEDS_PER_SEQ, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      if ( ierr .gt. 0 ) then
         write( lunerr, '(a,i,a)' ) 'GpWr(iteRng)Seeds: Error ', ierr,
     $        ' in Mpi_gather '
         goto 9999
      endif
#endif
#endif
c     Write out the header and then all the seeds, on node 0 only !!
c
      if ( nrank .eq. 0	  ) then
	  write(lunwr, 900 ) 'Version ', version,
     $         ' of pRNG seed file for parallel Geant '

	  write(lunwr, 905)  'RNG    ',   chNameRng,
     $         ' Identity of pseudorandom Number generator '
	  write(lunwr, 910)  'SUBSEQ   ',   nsubseq,
     $         ' Number  of  subsequences '
	  write(lunwr, 910)  'NOPERSEQ ', numperseq,
     $         ' Numbers per subsequence '
	  write(lunwr, 910)  'NPERLINE ', nperline,
     $         ' Numbers per line '
	  write(lunwr, 920)  'ENDPROLOG'


	  DO isubseq= 1, nsubseq
	      write(lunwr, 930 )  ' Subsequence ' , isubseq
	      DO inumb= 1, numperseq, nperline
		 imax=max(inumb+nperline-1,numperseq)
		 write( lunwr, 940 ) ( iseedrng( ino, isubseq ),
     $ 	             ino=inumb, imax)
	      ENDDO
	  ENDDO

	  if (lunwrite .lt. 0) then
	      close (unit=lunwr)
	  endif

      endif

 900  format( a8, f5.2, a )
 905  format( a7, a12, a )
 910  format( a9, i6, a )
 920  format( a )
 930  format( a14, i6 )
 940  format( 8i15 )


9999  continue

      return
      end
#endif
