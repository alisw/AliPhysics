*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_AADAT_INC
#define CERNLIB_GEANT321_AADAT_INC
*
*
* aadat.inc
*
*                                                                      *
*=== aadat ============================================================*
*                                                                      *
      PARAMETER (IPROM = 100)
      PARAMETER (ITARM = 100)
      COMMON/FKAADA/ENPERN,PLPERN,SIGAA(IPROM,ITARM),SIGNN,TMASS,AAEVNO,
     +             SELAA(IPROM,ITARM),RLASTP(IPROM,ITARM),
     +             MATPRO,MATTAR,LASTM,IPROA,IPROZ,ITARA,ITARZ
 

#endif
