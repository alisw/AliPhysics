CCC   common/correlations/
C
C     This common block includes the following data:
C
C          o Two-Body Correlation Arrays for the following cases:
C            - c2mod*  =  Model
C            - c2fit*  =  Fitted
C            - c2err*  =  Errors for fitted
C
C            Where (*) denotes the following cases:
C                      Like Pairs,   1D dist
C                      Unlike Pairs, 1D dist
C                      Like Pairs,   3D dist, Fine Mesh
C                      Unlike Pairs, 3D dist, Fine Mesh
C                      Like Pairs,   3D dist, Coarse Mesh
C                      Unlike Pairs, 3D dist, Coarse Mesh

      integer*4 max_c2_1d, max_c2_3d
      parameter (max_c2_1d = 100)
      parameter (max_c2_3d =  10)

      common/correlations/
     1  c2mod_like_1d(max_c2_1d), c2mod_unlike_1d(max_c2_1d),
     2  c2fit_like_1d(max_c2_1d), c2fit_unlike_1d(max_c2_1d),
     3  c2err_like_1d(max_c2_1d), c2err_unlike_1d(max_c2_1d),
     4  c2mod_like_3d_fine(max_c2_3d,max_c2_3d,max_c2_3d),
     5  c2mod_unlike_3d_fine(max_c2_3d,max_c2_3d,max_c2_3d),
     6  c2mod_like_3d_coarse(max_c2_3d,max_c2_3d,max_c2_3d),
     7  c2mod_unlike_3d_coarse(max_c2_3d,max_c2_3d,max_c2_3d),
     4  c2fit_like_3d_fine(max_c2_3d,max_c2_3d,max_c2_3d),
     5  c2fit_unlike_3d_fine(max_c2_3d,max_c2_3d,max_c2_3d),
     6  c2fit_like_3d_coarse(max_c2_3d,max_c2_3d,max_c2_3d),
     7  c2fit_unlike_3d_coarse(max_c2_3d,max_c2_3d,max_c2_3d),
     4  c2err_like_3d_fine(max_c2_3d,max_c2_3d,max_c2_3d),
     5  c2err_unlike_3d_fine(max_c2_3d,max_c2_3d,max_c2_3d),
     6  c2err_like_3d_coarse(max_c2_3d,max_c2_3d,max_c2_3d),
     7  c2err_unlike_3d_coarse(max_c2_3d,max_c2_3d,max_c2_3d)

CCC   Variable Type Declarations:

CCC   One-dimensional Functions:

      real*4 c2mod_like_1d, c2mod_unlike_1d
      real*4 c2fit_like_1d, c2fit_unlike_1d
      real*4 c2err_like_1d, c2err_unlike_1d

CCC   Three-dimensional Functions:

      real*4 c2mod_like_3d_fine,   c2mod_unlike_3d_fine
      real*4 c2mod_like_3d_coarse, c2mod_unlike_3d_coarse
      real*4 c2fit_like_3d_fine,   c2fit_unlike_3d_fine
      real*4 c2fit_like_3d_coarse, c2fit_unlike_3d_coarse
      real*4 c2err_like_3d_fine,   c2err_unlike_3d_fine
      real*4 c2err_like_3d_coarse, c2err_unlike_3d_coarse

CCC   END common/correlations/
