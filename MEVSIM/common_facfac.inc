CCC   Common FACFAC for Factorials

      integer factorial_max
      parameter (factorial_max = 10000) ! max # multiplicity per event;
CCC                                     ! for any specific particle ID;
CCC                                     ! also used for log(n!).
      Common/FACFAC/ FACLOG(factorial_max)
      real*4 FACLOG

