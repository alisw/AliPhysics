#ifndef CERNLIB_ISAJET_ISABNK_INC
#define CERNLIB_ISAJET_ISABNK_INC
*
*
* isabnk.inc
*
C
      COMMON/ISABNK/BANK,FILISA,FILIS2
      SAVE /ISABNK/
      CHARACTER*12 BANK
      CHARACTER*80 FILISA,FILIS2
C
C  If BANK='ISAP' Zebra bank ISAP (particles) will be written out
C  if BANK='ISAC'     "      ISAC (pseudo calorimeter) will be written out
C  If BANK='ISAL'     "      ISAL (leptons) will be written out
C  if BANK='ISAPISAC' both groups will be written out
C  if BANK='ISAPISACISAL' all groups will be written out
C
C  FILISA= name of ISAJET events file (ZEBRA)
C  FILIS2= name of second file if needed for output
C

#endif
