*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUPTDI(KFL,PX,PY)

C...Purpose: to generate transverse momentum according to a Gaussian.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.

C...Generate p_T and azimuthal angle, gives p_x and p_y.
      KFLA=IABS(KFL)
      PT=PARJ(21)*SQRT(-LOG(MAX(1E-10,RLU(0))))
      IF(MSTJ(91).EQ.1) PT=PARJ(22)*PT
      IF(KFLA.EQ.0.AND.MSTJ(13).LE.0) PT=0.
      PHI=PARU(2)*RLU(0)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)

      RETURN
      END
