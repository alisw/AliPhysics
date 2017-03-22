      SUBROUTINE SETPOWWGHT(POWER)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      COMMON/PYWGHT/WGHTPOW
      SAVE /PYWGHT/
C.. Divide power by 2 because we get p_T^2 from the Pythia
      WGHTPOW = POWER/2.D0
      RETURN
      END
      
C...PYEVWT
C...Implementation routine for event weighting
C...Uses a pt-hard based weighting with w ~ (p_T,hard/p0)^n
C...Where n can be set using a common block variable
C...p0 = 5 GeV/c
C...Multiplies the
C...standard PYTHIA differential cross-section by a process- and
C...kinematics-dependent factor WTXS. For MSTP(142)=1 this corresponds
C...to generation of weighted events, with weight 1/WTXS, while for
C...MSTP(142)=2 it corresponds to a modification of the underlying
C...physics.
 
      SUBROUTINE PYEVWT(WTXS)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYWGHT/WGHTPOW
      SAVE /PYDAT1/,/PYINT1/,/PYINT2/,/PYWGHT/
      DATA P02 /25.D0/

      PT2=VINT(48)
 
      WTXS=(PT2/P02)**WGHTPOW
 
      RETURN
      END
 
C*********************************************************************
