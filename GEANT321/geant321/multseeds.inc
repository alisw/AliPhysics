*CMZ :          13/07/94  11.16.35  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   13/07/94
C     Auxiliary variable for calculating pRNG seeds on different nodes.
      integer            iseeda, iseedb, iseqnc
      common /multseeds/ iseeda, iseedb, iseqnc
