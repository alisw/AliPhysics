*CMZ :          02/11/94  15.23.33  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   20/10/94
C   Variables to store information on the parallel configuration
      integer            nprank, npsize, nevtot, nfirstworker, npleader
      common /multiprox/ nprank, npsize, nevtot, nfirstworker, npleader
C     Descriptions:
C         npsize= number of processors being used
C         nprank= the rank of this processor ( 0 to npsize - 1 )
C         nevtot= total number of events requested (for all procs)
C         the last two are only used in the case that one node reads an event
C         input file, while the rest do the real simulation work
C         nfirstworker= rank of first process to simulate
C         npleader= rank of the process that reads the input
