*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:40  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_MAPA_INC
#define CERNLIB_GEANT321_MAPA_INC
*
*
* mapa.inc
*
*----------------------------------------------------------------------*
*     include file: mapa copy                    created 26/11/86 by pa*
*     note: change also mapa2 copy if you touch this file              *
*     changes: none                                                    *
*     included in the following subroutines or functions: not updated  *
*                                                                      *
*     description of the common block(s) and variable(s)               *
*                                                                      *
*                                                                      *
*     /mapa/ contains the parameters describing the target materials   *
*        matnam = matnam(i)=name of the material number i              *
*        medium = medium(nreg)=material number of the region nreg      *
*        amss   = the atomic weight of the materials                   *
*        aocmbm = atomic densities in barn**-1 cm**-1                  *
*                 (Atoms Over Cm times Barn for Materials)             *
*        ztar   = atomic numbers of the materials                      *
*        rho    = densities of the materials                           *
*        zlin   = inelastic scattering lengths of the materials        *
*                 for beam particles at the average beam energy in cm  *
*        zlel   = elastic scattering lengths of the materials for      *
*                 beam particles at average beam energy in cm          *
*        zlrad  = radiation lengths of the materials in cm             *
*        zul    = inelastic scattering lengths of the materials        *
*                 for neutrons at threshold energy in cm               *
*        mulflg = flags for multiple scattering options                *
*        icomp  = 0 if the material is not a compound                  *
*        mssnum = mass number of the target nucleus, if =< 0 it means  *
*                 that it is in the natural isotopic composition       *
*        nregs  = total number of the regions                          *
*        nmat   = total number of the materials                        *
*        mtbsnm = medium for which inelastic interaction biasing must  *
*                 be done                                              *
*                                                                      *
*                        Mxxrgn = maximum number of regions            *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON / FKMAPA / AOCMBM (MXXMDF), AMSS (MXXMDF), ZTAR   (MXXMDF),
     &                  RHO    (MXXMDF), ZLIN (MXXMDF), ZLEL   (MXXMDF),
     &                  ZLRAD  (MXXMDF), ZUL  (MXXMDF), MEDIUM (MXXRGN),
     &                  MULFLG (MXXMDF), ICOMP(MXXMDF), MSSNUM (MXXMDF),
     &                  NREGS, NMAT, MTBSNM
 

#endif
