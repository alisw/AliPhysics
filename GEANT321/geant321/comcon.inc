*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_COMCON_INC
#define CERNLIB_GEANT321_COMCON_INC
*
*
* comcon.inc
*
*
*=== comcon ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     include file: comcon copy                 created 26/11/86 by pa *
*                                                                      *
*     Last change  on  30-oct-92   by   Alfredo Ferrari, INFN - Milan  *
*     included in the following subroutines or functions: not updated  *
*                                                                      *
*     description of the common block(s) and variable(s)               *
*                                                                      *
*                                                                      *
*     /comcon/ contains information about compounds                    *
*         icompl(imat) = number of materials in a compound             *
*         matnum(i)    = material numbers (compounds consist of the    *
*                        materials)                                    *
*         aocmbc(i)    = atomic densities in barn**-1 cm**-1           *
*                        (Atoms Over Cm times Barn for Compounds)      *
*         cona13(i)    = A^1/3 of the constituents                     *
*         sumzmf(i)    = cumulative z^2 for multiple scattering        *
*         contnt(i)    = partial densities                             *
*         sumacn(i)    = cumulative atom contents of the compounds     *
*         cabinx(i)    = cumulative at.cont.*inelastic cross section   *
*         cabelx(i)    = cumulative at.cont.*  elastic cross section   *
*         anxnor       = total inelastic cross section                 *
*         elxnor       = total   elastic cross section                 *
*----------------------------------------------------------------------*
*
      COMMON /FKCOMC/ AOCMBC (ICOMAX), CONA13 (ICOMAX), SUMZMF (ICOMAX),
     &                CONTNT (ICOMAX), SUMACN (ICOMAX), CABINX (ICOMAX),
     &                CABELX (ICOMAX), ANXNOR, ELXNOR,  ICOMPL (MXXMDF),
     &                MATNUM (ICOMAX)
 

#endif
