*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXSCAH(ID,NAME,CHOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates and plot an histogram of physics quantities such as *
C.    *  the total number of radiation lengths or interaction lengths  *
C.    *  in function of the SCAN parameter TETA for a value of PHI.    *
C.    *    CHOPT='O' histogram is generated at Exit  of VOLUME.        *
C.    *    CHOPT='I' histogram is generated at Entry of VOLUME.        *
C.    *    CHOPT='X' radiation lengths                                 *
C.    *    CHOPT='L' Interaction lengths                               *
C.    *    CHOPT='P' Plot the histogram                                *
C.    *  If VOLUME='XXXX' Mother volume is used.                       *
C.    *  The histogram identifier IDPHI is used to also identify which *
C.    *  PHI division to plot. IPHI=MOD(IDPHI,1000).                   *
C.    *  If IPHI=0, then all PHI divisions are generated (not plotted) *
C.    *  with histogram identifiers IDPHI+PHI division number.         *
C.    *                                                                *
C.    *       Author:    R.Brun      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcscan.inc"
      CHARACTER*(*) NAME,CHOPT
C.
C.    ------------------------------------------------------------------
C.
  99  END
 
