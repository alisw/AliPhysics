CCC   common/coulmb/
C
C     This common block includes the following data:
C
C          o input q mesh (GeV/c) 
C          o interpolated Coulomb correction for a specified
C            Coulomb source size (input as parameter Q0 when
C            switch_coulomb = 3) for both like and unlike pairs.
C        
C          The correction corresponds to John Cramer's tabulation
C          of Scott Pratt's Coulomb wave function integration code
C          results for finite size spherical sources.  The present
C          code reads these tables and interpolates to the requested
C          source radius.
C
C          The input source radii, the number of input files of data
C          and the q-space grid are assumed static and are therefore
C          hard-coded as follows:

      integer*4 max_c2_coul         ! # q-mesh pts for input Coul. correc.
      integer*4 ncoulradsteps       ! # input radii values (= # files)
      parameter (max_c2_coul = 288)
      parameter (ncoulradsteps = 7)

      real*4 coulradmin     ! min source radius (fm) for Coul Correc. tables
      real*4 coulradmax     ! max source radius (fm) for Coul Correc. tables
      real*4 coulradstep    ! Source radius step size for Coul Correc. tables
      parameter (coulradmin  =  6.0)
      parameter (coulradmax  = 18.0)
      parameter (coulradstep =  2.0)

      common/coulmb/
     1   c2_coul_like(max_c2_coul),c2_coul_unlike(max_c2_coul),
     2   q_coul(max_c2_coul)

      real*4 c2_coul_like   ! Finite source size Coul. correc for like pairs
      real*4 c2_coul_unlike ! Finite source size Coul. correc for unlike pairs
      real*4 q_coul         ! Momentum transfer grid (GeV/c).

CCC   END common/coulmb/
