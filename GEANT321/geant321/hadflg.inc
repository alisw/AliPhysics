*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:39  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_HADFLG_INC
#define CERNLIB_GEANT321_HADFLG_INC
*
*
* hadflg.inc
*
*
*=== hadflg ===========================================================*
*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 13 september 1991 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on  1-dec-91     by    Alfredo Ferrari               *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                                         Calumo                       *
*                                         Datesh                       *
*                                         Eventv                       *
*                                         Ferhav                       *
*                                         Hadriv                       *
*                                         Rchanv                       *
*                                                                      *
*                    Ielflg = Flag to control the elastic channel      *
*                             in Hadriv:                               *
*                             -1: elastic channel reduced according    *
*                                 to a rough Pauli blocking factor     *
*                              0: elastic channel fully accounted for  *
*                             +1: elastic channel suppressed           *
*                    Icxflg = Flag to control the charge exchange      *
*                             channel in Hadriv:                       *
*                             -1: ch. ex. channel reduced according    *
*                                 to a rough Pauli blocking factor     *
*                              0: ch. ex. channel fully accounted for  *
*                             +1: ch. ex. channel suppressed           *
*             Innure(1,1,k) = Projectile index for the Hadriv reaction *
*                             k (first entrance channel)               *
*             Innure(2,1,k) = Target index for the Hadriv reaction     *
*                             k (first entrance channel)               *
*             Innure(1,2,k) = Projectile index for the Hadriv reaction *
*                             k (second entrance channel, if any)      *
*             Innure(2,2,k) = Target index for the Hadriv reaction     *
*                             k (second entrance channel, if any)      *
*                             It holds:                                *
*                        Nure (Innure(1,j,k),(Innure(2,j,k)/8+1)) = k  *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER ( PPAMXB = 0.6       D+00 )
      PARAMETER ( PAUMXB = 2.D+00 * 0.8 D+00 )
      PARAMETER ( PPAMXM = 0.6       D+00 )
      PARAMETER ( PAUMXM = 2.D+00 * 0.4 D+00 )
*
      COMMON / FKHADF / IELFLG, ICXFLG, IKCHXG (16), INNURE (2,2,16),
     &                  IOLDHD
 

#endif
