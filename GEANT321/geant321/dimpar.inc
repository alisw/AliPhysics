*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:38  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_DIMPAR_INC
#define CERNLIB_GEANT321_DIMPAR_INC
*
*
* dimpar.inc
*
*                                                                     *
*=== dimpar ==========================================================*
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*      DIMPAR: included in any routine                                *
*                                                                     *
*          Mxxrgn = maximum number of regions                         *
*          Mxxmdf = maximum number of media in Fluka                  *
*          Mxxmde = maximum number of media in Emf                    *
*          Mfstck = stack dimension in Fluka                          *
*          Mestck = stack dimension in Emf                            *
*          Nallwp = number of allowed particles                       *
*          Mpdpdx = number of particle types for which EM dE/dx pro-  *
*                   cesses (ion,pair,bremss) have to be computed      *
*          Icomax = maximum number of materials for compounds (equal  *
*                   to the sum of the number of materials for every   *
*                   compound )                                        *
*          Nstbis = number of stable isotopes recorded in common iso- *
*                   top                                               *
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
      PARAMETER ( MXXRGN = 500  )
      PARAMETER ( MXXMDF = 50   )
      PARAMETER ( MXXMDE = 18   )
      PARAMETER ( MFSTCK = 1000 )
      PARAMETER ( MESTCK = 100  )
      PARAMETER ( NALLWP = 39   )
      PARAMETER ( MPDPDX = 8    )
      PARAMETER ( ICOMAX = 120  )
      PARAMETER ( NSTBIS = 304  )
 

#endif
