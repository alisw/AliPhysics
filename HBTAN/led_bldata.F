
      SUBROUTINE led_bldata
c--mlv      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_AAPI/AAPI(20,2)/FSI_AAND/AAND(20,4)
      COMMON/FSI_AADD/AADD(20,3)/FSI_AAPIN/AAPIN(20,2)
      COMMON/FSI_AAKK/AAKK(9)
      
C---Parameters for GPIPI (I,J), J=1(2) corresp. to isospin=0(2)
      DATA AAPI/.2600D00, .2500D00, .3480D00,-.0637D00, 16*0.D0,
     1         -.0280D00,-.0820D00, .2795D00,-.0086D00, 16*0.D0/
C---Parameters for GPIN (I,J), J=1(2) corresp. to piN isospin=1/2(3/2)
      DATA AAPIN/ .12265D1, .1563D2,     18*0.D0,
     1            -.750D0,  .7688D2,     18*0.D0/
C---Parameters for GND (I,J), J=1(2) corresp. to nd(pd), S=1/2,
C---                          J=3(4) corresp. to nd(pd), S=3/2
      DATA AAND/-.3295D1,-.8819D3, .4537D4, .28645D5, 16*0.D0,
     1          -.13837D2, .11505D2, .0D0, .10416D2,  16*0.D0,
     2          -.32180D2, .1014D2,  .0D0, .0D0,      16*0.D0,
     3          -.60213D2, .1333D2,  .0D0,-.70309D2,  16*0.D0/
      DATA AADD/ .10617D4, .3194D-2, .56849D3, 17*0.D0,
     1           20*0.D0,
     2          -.1085D4, .1987D5, 18*0.D0/
C--- AAKK= m_K^2, m_pi^2, m_eta^2, m_S*^2, m_delta^2,
C---       gam(S*-->KK-b), gam(S*-->pipi), gam(delta-->KK-b),
C---       gam(delta-->pi eta)       
      DATA AAKK/.247677D00,.01947977D00,.2997015D00,.9604D00,
     1          .96511D00,
cc     2          .792D00, .199D00, .333D00, .222D00/ ! Martin (77)
     2          .094D00, .110D00, .333D00, .222D00/ ! Morgan (93)      
  
c-       write(*,*)'aand',(aand(i,1),i=1,20)

      END
