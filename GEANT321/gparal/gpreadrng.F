*CMZ :          26/02/96  11.38.47  by  S.Ravndal
*-- Author :
C------------------------------------------------------------------------
C     SUBROUTINE gpreadrng( lunread, filename, ierr )
C
C     Function: read a file containing a state of the pseudo random number
C               generator on all the nodes and set the RNG seeds using it.
c
c     A convention: if lunread is >= 0, use that unit.
c     If not, read from the file called <filename>
c
C     Possible error codes:
c	  ierr= -102         The number of seeds required for a subsequence of
c                            the RNG (numperseq) is more than allocated
c                	     in the array iseedrng (SEEDS_PER_SEQ)
c
c     called by:  <USER> (if PARALLEL switch is used)
c
c     Author:     John Apostolakis,  January 1996.
c
c     Last Modified: February 22, 1996   J.A.
c------------------------------------------------------------------------
#if defined(CERNLIB_PARA)
      subroutine gpreadrng( lunread, filename, ierr )
      implicit none
#if defined(CERNLIB_PARA_COMM)
#include "geant321/mpifinc.inc"
#endif
      integer lunread, ierr
      character filename*(*)
c     Containes iseeda, iseedb, iseqnc
#include "geant321/multseeds.inc" 
      integer    numprocs, nrank, nleader
      integer    MAXNODES, SEEDS_PER_SEQ
      parameter (MAXNODES=4096)
c     RANMAR
      parameter (SEEDS_PER_SEQ=103)

      integer   myseeds(SEEDS_PER_SEQ)
      integer   iseedrng( SEEDS_PER_SEQ, MAXNODES )

#include "geant321/gcunit.inc" 
c     version of the seed file format
      real   version, versionrd
      integer nsubseq, numperseq, nperline, nsubseqexp, nsubseqrd
      integer lunrd, lunerr
      integer isubseq, inumb, imax, ino, ichkerr(4), iia, iseqrd
      character*12 chSUBSEQ, chNOPERSEQ, chNPERLINE, chENDPROLOG,
     $     chsubsequence , chVersion
      character*12 chNameRng, chNameRead, chIdRng
      data  version / 1.0 /

      lunerr= lout

      if (lunread .lt. 0) then
          lunrd= 74   !! Some number ...
          open (unit=lunrd, file=filename, status='old')
      else
          lunrd= lunread
      endif

c #if defined(CERNLIB_PARA_RANECU)
c     if( ParaRng .eq. 'RANECU' )
c     --------------------------------------------------------------------
c     The RANECU pseudo random number generator, used by GEANT's GRNDM(Q).
c     --------------------------------------------------------------------
#if  !defined(CERNLIB_PARA_TASKRNG)
c     Currently the number of sequences must be the number of processors!
      call gprocs( numprocs, nrank, nleader )
      nsubseqexp=  numprocs

      if( numprocs .gt. MAXNODES )then
c        write(lunerr, *) 'GPreadRng: ERROR Too many nodes ', numprocs
         write(*, *) 'GPreadRng: ERROR Too many nodes ', numprocs
      endif
#endif
      chnamerng= 'RANECU'
      numperseq= 2              !  2 "seed" integers per RANECU sequence
      nperline=  2              !  2 numbers per line
c #endif 

#if defined(CERNLIB_PARA_RANMAR)
#endif

      read ( lunrd, 905 )  chVersion, versionrd
 905  format( a8, f5.2, a )
      write ( *, * ) ' Reading seedfile version ', version

      read(lunrd, 910) chIdRng,    chNameRead
 910  format( a7, a12 )
      read(lunrd, 920) chSUBSEQ,   nsubseqrd
 920  format( a9, i6 )
      read(lunrd, 920) chNOPERSEQ, numperseq
      read(lunrd, 920) chNPERLINE, nperline
      read(lunrd, 930) chENDPROLOG
 930  format( a9 )

c     Check that the strings before the values were correctly read in.
c
      call checkstr( chSUBSEQ,    'SUBSEQ ',   7, 
     $     'Number  of  subsequences ', ichkerr(1) )
      call checkstr( chNOPERSEQ,  'NOPERSEQ ', 9,
     $     ' Numbers per subsequence ',  ichkerr(2) )
      call checkstr( chNPERLINE,  'NPERLINE ', 9,
     $     ' Numbers per line ',         ichkerr(3) )
      call checkstr( chENDPROLOG, 'ENDPROLOG', 9,
     $     ' No variable to be read ',   ichkerr(4) )
      DO iia=1,4
          IF( ichkerr(iia) .ne. 0 ) then
             write(*,*) ' GpReadRng : STOPping because of error ',
     $            '  in reading overall header, part ', iia
             stop
          ENDIF
      ENDDO

C     CHECK the VALUES read in
C
      IF ( nsubseqexp .ne. nsubseqrd ) THEN
         write(*,*)
     $        ' GpReadRng: Error: Mismatch in number of subsequences, ',
     $        ' expected = ' , nsubseqexp, ' file = ', nsubseqrd
         IF ( nsubseqexp .lt. nsubseqrd ) THEN
            write(*,*) ' Enough seeds to work with, continuing '
            nsubseq= nsubseqexp
         ELSE
            write(*,*) ' Not enough subsequences, aborting '
            call gpabort()
         ENDIF
      ELSE
         nsubseq= nsubseqexp
      ENDIF

c     Now read in the values from each subsequence
c
      DO isubseq= 1, nsubseq
          read(lunrd, 970 )  chSubsequence, iseqrd
          call checkstr( chSubsequence, ' Subsequence', 11,
     $         ' Which of the Subsequences ', ichkerr(1) )
          IF( (ichkerr(1).ne.0) .or. (iseqrd .ne. isubseq) ) then
             write(*,*) ' GpReadRng : STOPping because of error ',
     $                    ' in reading subsequence header '
             write(*,*) ' chSubsequence=/', chSubsequence, '/, iseqrd=',
     $                    iseqrd, ' iseq-expected= ', isubseq
             stop
          ENDIF
          DO inumb= 1, numperseq, nperline
             imax=max(inumb+nperline-1,numperseq)
             read( lunrd, 980 ) ( iseedrng( ino, isubseq ),
     $            ino=inumb, imax)
          ENDDO
      ENDDO
 970  format( a14, i6 )
 980  format( 8i15 )
#if defined(CERNLIB_PARA_RANECU)
#if !defined(CERNLIB_PARA_TASKRNG)

#if defined(CERNLIB_PARA_COMM)
c     Now "scatter" all seeds from node 0 to all the nodes
c
c     if numperseq < SEEDS_PER_SEQ the following should still work.
c
      call MPI_Scatter( iseedrng, SEEDS_PER_SEQ, MPI_INTEGER,
     $                  myseeds,  numperseq,     MPI_INTEGER,
     $                  0, MPI_COMM_WORLD, ierr)
c
      if ( ierr .gt. 0 ) then
         write( lunerr, '(a,i,a)' ) 'GpDefRng: Error ', ierr,
     $        ' in Mpi_gather '
      endif
#else
c     Else each node selects its seeds ...
c
      DO inumb= 1, numperseq
	  myseeds( inumb )= iseedrng( inumb, nrank+1 )
      ENDDO
#endif

c     Set the pRNG seeds (In this version there is one sequence per node)
      call grndmq( myseeds(1), myseeds(2), iseqnc, 'S' )
#endif
#endif
      if (lunread .lt. 0) then
          close (unit=lunrd)
      endif

      return
      end
#endif
