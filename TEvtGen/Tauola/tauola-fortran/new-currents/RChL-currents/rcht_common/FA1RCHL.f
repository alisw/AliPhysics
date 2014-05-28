      COMPLEX FUNCTION FA1RCHL(QQ)
      IMPLICIT NONE
      REAL                     QQ
      DOUBLE PRECISION M1,M2,M3
      REAL GGMA1
      REAL wid_a1_fit 
C.......................................................................
C.
C.    FA1CHL - RchT version of  the A1 propagator
C.
C.    Inputs    : QQ - invariant masses**2  [GeV**2]
C.    Outputs   : FA1RCHL formfactor value at QQ
C.
C.    COMMON    : RCHT_3PI content is defined in this routine
C.
C.    Calls     : functions from file ./wid_a1_fit.f
C.    Called    : from file f3pi_rcht.f, fkkpi.f, fkk0pi0.f
C************************************************************************
       include '../parameter.inc'
       include '../funct_declar.inc'
C******************************************
C    Initilisation of the mass of the particles
C*****************************************
        call rchl_parameters(5)

c$$$C we impose isospin symmetry requesting that charged and neutral pion mass
c$$$C are equal. This may need to be changed
c$$$        MMPI_AV = (2.*MPIC+MPIZ)/3.

	M1 = MMPI_AV      
        M2 = MMPI_AV     
        M3 = MMPI_AV  
c
C   Function wid_a1_fit.f calculates the energy dependence of 
C   the a1 meson width
C   
      IF(QQ.GE.(M1+M2+M3)**2) THEN
        GGMA1 = wid_a1_fit(QQ)
      ELSE 
        GGMA1 = 0.
      ENDIF


      FA1RCHL = 1./(QQ-MMA1*MMA1+i*MMA1*GGMA1)

      RETURN
      END

C to switch on/off remotely option for scalar contr. 
C to FORM1 FORM2 of 3 pi mode.     
      subroutine getFF3PISCAL(INUM)
      INTEGER INUM
      include '../parameter.inc'
      INUM=FF3PISCAL
      return
      end

      subroutine setFF3PISCAL(INUM)
      INTEGER INUM
      include '../parameter.inc'
      FF3PISCAL=INUM
      return
      end
