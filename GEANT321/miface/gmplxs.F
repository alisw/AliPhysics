*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.47  by  S.Giani
*-- Author :
      SUBROUTINE GMPLXS(D,LD,IGAMS,LGAM,INABS,LNAB,ITHRMS,LTHRM,
     + IDICTS,LDICT,NTX,NTS,IGCBS,LGCB,AWR,Q,NSEI,NAEI,NMT2,NMT4,
     + NMT16,NMT17,NMT18,NMT22,NMT23,NMT24,NMT28,NMT51,NMT91,
     + NMT102,NMT103,NMT104,NMT105,NMT106,NMT107,NMT108,NMT109,
     + NMT111,NMT112,NMT113,NMT114,IGCBS2,LGCB2,KZ,LR,QLR,
     + IIN,IIM,ID,SIGNN)
C   This routine is a copy of COLISN, used to compute
C   the cross-section of low-energy neutrons processes.
C   The calling sequence is as the one is COLISN plus
C   the ID number of the process and the output value SIGNN
C                            Kati Lassila-Perini 3.12.94
C
#include "geant321/minput.inc"
#include "geant321/mconst.inc"
#include "geant321/mnutrn.inc"
#include "geant321/mapoll.inc"
#include "geant321/mcross.inc"
#include "geant321/mmass.inc"
#include "geant321/mupsca.inc"
#include "geant321/mpstor.inc"
#include "geant321/mmicab.inc"
      DIMENSION D(*),LD(*),IGAMS(*),LGAM(*),INABS(*),LNAB(*),
     + ITHRMS(*),LTHRM(*),IDICTS(NNR,NNUC),LDICT(NNR,NNUC),NTX(*),
     + NTS(*),IGCBS(NGR,NNUC),LGCB(NGR,NNUC),AWR(*),Q(NQ,NNUC),
     + NSEI(*),NAEI(*),NMT2(*),NMT4(*),NMT16(1),NMT17(*),NMT18(*),
     + NMT22(*),NMT23(*),NMT24(*),NMT28(*),NMT51(*),NMT91(*),
     + NMT102(*),NMT103(*),NMT104(*),NMT105(*),NMT106(*),NMT107(*),
     + NMT108(*),NMT109(*),NMT111(*),NMT112(*),NMT113(*),NMT114(*),
     + IGCBS2(NGR,NNUC),LGCB2(NGR,NNUC),KZ(*),LR(NQ,NNUC),QLR(NQ,NNUC),
     + FM(MAXNEU)
C
C
      CALL GTMED(NMED,MED)
      SIGNN = 0.
C
      IF (ID.EQ.99) THEN
C
C       DETERMINE THE TOTAL NEUTRON DISAPPEARANCE (MT-102 TO MT-114
C       AND MT-18)
         L1=LNAB(IIN)
         IF(L1.EQ.0) THEN
            SIGNN = 0.
         ELSE
            LS1=INABS(IIN)+LMOX2
            LEN=L1/2
            CALL TBSPLT(D(LS1),E,LEN,SIGNN)
         ENDIF
      ELSE
C
         L1=LDICT(ID,IIN)
         IF(L1.EQ.0) THEN
            SIGNN = 0.
         ELSE
            LS1=IDICTS(ID,IIN)+LMOX2
            LEN=L1/2
            CALL TBSPLT(D(LS1),E,LEN,SIGNN)
C
            IF(ID.EQ.10) THEN
C       THE TREATMENT OF THE FISSION REACTION ASSUMES THE FISSION
C       CROSS SECTION IS STORED AS NUBAR*SIGF
               L1=LDICT(134,IIN)
               IF(L1.EQ.0)THEN
                  SIGNN = 0.0
               ELSE
                  LS1=IDICTS(134,IIN)+LMOX2
                  LEN=L1
                  CALL GETNU(D(LS1),LD(LS1),EOLD,LEN,XNU)
                  SIGNN=SIGNN/XNU
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
      END
