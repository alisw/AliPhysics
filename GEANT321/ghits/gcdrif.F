*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:08  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE  GCDRIF (RADD,ZMIN,ZMAX,DETREP,HITREP,IOUT)
C.
C.    **************************************************************************
C.    *                                                                        *
C.    *              Digitisation of Drift-  Chambers                          *
C.    *              --------------------------------                          *
C.    *                                                                        *
C.    *   Digitisation routine for a cylindrical drift chamber.                *
C.    * RADD      radius of cylinder in cm                                     *
C.    * ZMIN      z of lower end of cylinder                                   *
C.    * ZMAX      z of upper end of "                                          *
C.    * DETREP(1) number of wires                                              *
C.    * DETREP(2) wire spacing in PHI (radians)                                *
C.    * DETREP(3) cosine of wire angle                                         *
C.    * DETREP(4) sine of wire angle (signed like dphi/dz)                     *
C.    * DETREP(5) dphi/dz along wire                                           *
C.    * DETREP(6)  phi of point with z=0 on wire 1                             *
C.    * DETREP(7) drift velocity (cm/nsec)                                     *
C.    * DETREP(8) quantity describing the drift angle                          *
C.    *           if.ne.0 ==> user routine GUDTIM                              *
C.    * HITREP(1) phi coordinate of intersection                               *
C.    * HITREP(2) z coordinate                                                 *
C.    * HITREP(3) dphi/dr                                                      *
C.    * HITREP(4) dz/dr                                                        *
C.    * IOUT(1)   wire number (1..NWI with  increasing phi)  (-1 for           *
C.    *           bad DETREP parameters)                                       *
C.    * IOUT(2)   drift time (nsec) (+/- for phi(hit)>/< phi(wire)             *
C.    * IOUT(3)   digitised current division information (rel.  pos.           *
C.    *           along wire of charge) (per mille)                            *
C.    * IOUT(4)   amount of charge deposited to wire                           *
C.    * Coordinate systems along wire                                          *
C.    *     I.             Charge                   I.                         *
C.    *     .              |                        .                          *
C.    *     |              .                        |                          *
C.    *     =========================================  SENSE WIRE              *
C.    * ...................................................> Z (cm)            *
C.    *     Z              Z.                       Z.                         *
C.    *      L                                                                 *
C.    * ...............................................> X (arbitrary scale)   *
C.    *     0              X.                       L                          *
C.    *           X.                  (L-X.)                                   *
C.    * The scaling used is such that L . 1000.                                *
C.    * Knowing the position Z. of the deposit of charge,                      *
C.    *            Z.-ZL                                                       *
C.    *    X.                          =                           L           *
C.    * .          .....                                                       *
C.    *            Z.-ZL                                                       *
C.    * This information is stored into IOUT(3).                               *
C.    *   Routine to  calculate the  error on  the current  division           *
C.    * information as obtained by "GCDRIFT".                                  *
C.    * ICD       digitized     current     division     information           *
C.    *           (0 ... 1000)                                                 *
C.    * ERP       variance of Gaussian  distributed pedestral errors           *
C.    *           on the measured pulse heights  relative to the sum           *
C.    *           of the pulse heights                                         *
C.    * ERS       variance of  Gaussian distributed slope  errors on           *
C.    *           the measured  pulse heights  relative to  the each           *
C.    *           pulse heights                                                *
C.    *   Here we  assume that X.  has been determined  by measuring           *
C.    * the pulse heights I., I. with some statistical errors.                 *
C.    * X.   is then given by the formula                                      *
C.    *     X. = L . I./I.   with  I. . I.+I.                                  *
C.    * and its error is determined by                                         *
C.    *     .X. = -(X./I.) .I. + (L-X./I.) .I.                                 *
C.    * with the errors on measuring the pulse heights                         *
C.    *      .I. = .. + ...I.                                                  *
C.    *      .I. = .. + ...I.                                                  *
C.    * ..,  .. are of dimension  (I)  and represent the "pedestral"           *
C.    * errors;                                                                *
C.    * .., .. are the "slope" errors.                                         *
C.    * All  are   assumed  to  be  distributed   independently  (no           *
C.    * correlations), randomly and Gaussian around zero. This gives           *
C.    * the final result                                                       *
C.    *           ..      ..                  X.(L-X.)                         *
C.    *   .X. = - .. X. + .. (L-X.) + (..-..) ........                         *
C.    *           I.      I.                     L                             *
C.    *         ..................   .................                         *
C.    *               "pedestal"          "slope"                              *
C.    * In  GCDERR,  the X. derived from GCDRIF is set to                      *
C.    *      X. = X. + .X.    (but 0 . X. . L)                                 *
C.    * using  ERP .....  variance for ./I.  ,  ../I.  distributions           *
C.    *        ERS .....  variance for .., ..  distributions.                  *
C.    *                                                                        *
C.    *    ==>Called by : <USER>, GUDIGI                                       *
C.    *       Author    D.Mitaroff *********                                   *
C.    *                                                                        *
C.    **************************************************************************
C.
      DIMENSION  DETREP(8), HITREP(4), IOUT(4)
#include "geant321/gconsp.inc"
C.
C.    -----------------------------------------------------------------
C.
      ZREL    = 1000.
      IOUT(1) = -1
      NWI = DETREP(1)
      WSP = DETREP(2)
      DVL = DETREP(7)
      IF (WSP .EQ. 0.)   GOTO 99
      IF (DVL .EQ. 0.)   GOTO 99
C
C---- CALCULATE WIRE NUMBER.
      FI = HITREP(1)
      ZZ = HITREP(2)
      FI0 = DETREP(6) + ZZ * DETREP(5)
      DFI = FI - FI0
   10 IF (DFI .GE. 0.)    GOTO 11
      DFI = DFI + TWOPI
      GOTO 10
   11 IF (DFI .LT. TWOPI)   GOTO 12
      DFI = DFI - TWOPI
      GOTO 11
   12 IW = DFI / WSP + 0.5
      DIS = DFI - IW * WSP
      IF (IW .EQ. NWI)   IW = 0
      IOUT(1) = IW + 1
C
C---- CALCULATE DRIFT TIME.
      DIS = DIS * RADD * DETREP(3)
      IF ( DETREP(8) .NE. 0. )   GOTO 2
      IOUT(2) = DIS / DVL
      GOTO 3
C
C---- DRIFT TIME BY USER ROUTINE.
    2 IOUT(2) = GUDTIM (DETREP,HITREP,IW+1,DIS)
C
C---- CALCULATE CURRENT DIVISION INFORMATION.
    3 Z0 = ZZ + DIS * DETREP(4)
      IF (Z0 .LT. ZMIN)   Z0 = ZMIN
      IF (Z0 .GT. ZMAX)   Z0 = ZMAX
      IOUT(3) = ZREL * (Z0 - ZMIN) / (ZMAX - ZMIN)
      IOUT(4) = 0
C.
  99  RETURN
      END
