*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:39  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_HIGFIS_INC
#define CERNLIB_GEANT321_HIGFIS_INC
*
*
* higfis.inc
*
*
*=== higfis ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     HIGh energy FISsion common:                                      *
*                                                                      *
*     Created   on   9 april 1993   by   Alfredo Ferrari & Paola Sala  *
*                                             INFN - Milan             *
*     Last change on   11-apr-93    by   Alfredo Ferrari, INFN - Milan *
*                                                                      *
*          Fisinh = logical flag for HE fission                        *
*          Afis(i)= ith fragment mass number just after fission        *
*          Zfis(i)= ith fragment atomic number just after fission      *
*          Ufis(i)= ith fragment excitation energy (MeV) just after    *
*                   fission                                            *
*         Ekfis(i)= ith fragment kinetic(MeV) energy just after fission*
*         Amfis(i)= ith fragment atomic mass (MeV) just after fission  *
*         Ppfis(i)= ith fragment momentum (MeV/c) just after fission   *
*         Atfis(i)= ith fragment mass number after evaporation         *
*         Ztfis(i)= ith fragment atomic number after evaporation       *
*         Utfis(i)= ith fragment excitation energy (MeV) after evapora-*
*                   tion                                               *
*        Recfis(i)= ith fragment kinetic(MeV) energy after evaporation *
*        Atmfis(i)= ith fragment atomic mass (MeV) after evaporation   *
*        Pptfis(i)= ith fragment momentum (MeV/c) after evaporation    *
*          Ebfiss = fission barrier (MeV)                              *
*          Amdiff = Fission fragment kinetic energy (MeV)              *
*            Apr0 = Mass number of the (first) fissioning nucleus      *
*            Zpr0 = Atomic number of the (first) fissioning nucleus    *
*             Uu0 = Excitation energy (MeV) of the (first) fissioning  *
*                   nucleus                                            *
*           Erec0 = Recoil kinetic energy (MeV) of the (first) fissio- *
*                   ning nucleus                                       *
*      Npartf(j,i)= Number of evaporated particle of type j after the  *
*                   evaporation of the ith fragment                    *
*        Hevfis(i)= Evaporated "heavies" kinetic energy (MeV) after the*
*                   evaporation of the ith fragment                    *
*           Nfiss = Number of fission fragments in the stack           *
*                                                                      *
*----------------------------------------------------------------------*
*
      LOGICAL FISINH
      COMMON / FKHIGF / AFIS  (10), ZFIS  (10), UFIS  (10), EKFIS  (10),
     &                  AMFIS (10), PPFIS (10), COSLFF(3,0:10),
     &                  ATFIS (10), ZTFIS (10), UTFIS (10), RECFIS (10),
     &                  AMTFIS(10), PPTFIS(10), EBFISS, AMDIFF,
     &                  APR0, ZPR0, EREC0, UU0, HEVFIS (0:10),
     &                  NPARTF (6,0:10), NFISS, FISINH
      DIMENSION COSLF0 (3)
      EQUIVALENCE ( COSLF0 (1), COSLFF (1,0) )

#endif
