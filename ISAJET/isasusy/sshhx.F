#include "isajet/pilot.h"
      SUBROUTINE SSHHX
C-----------------------------------------------------------------------
C     Calculates the decays Hi -> Hj + X.
C
C     Includes vertex corrections for triple Higgs couplings due
C     to top and stop quarks effects.
C     See Kunszt and Zwirner CERN-TH.6150/91 for all but hh-hc-hc
C     correction which is in our Higgs-->SUSY paper:
C     Baer et. al. FSU-HEP-920630 or UH-511-749-92.
C 
C     The hh-hl-hl vertex correction now includes both 
C        top & bottom and stop and sbottom squark
C        (non-degenerate with mixing) effects.  
C        A-terms and mu=-2m1 are also included.
C
C
C     Bisset's HIGPRO
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/sspar.inc"
#include "isajet/sssm.inc"
#include "isajet/sstype.inc"
C
      DOUBLE PRECISION PI,SR2,G2,GP2,BETA,ALPHA,SW2,CW2,LGTST,CBMA
     $,SBMA,LAMB1,DWID,DELLPP,MH,M1,M2,LAMB,TEMP,DTEMPL,DTEMPR
     $,DELHLL,DELHPP,DELHCC,CAB2,SAB2
      DOUBLE PRECISION SSDLAM
      REAL WID,ASMT,MTMT,MTQ,SUALFS,HIGFRZ
      DOUBLE PRECISION SSMQCD
C
      PI=4.*ATAN(1.D0)
      SR2=SQRT(2.D0)
      G2=4.0*PI*ALFAEM/SN2THW
      GP2=4*PI*ALFAEM/(1-SN2THW)
      HIGFRZ=SQRT(AMTLSS*AMTRSS)
      ASMT=SUALFS(AMTP**2,.36,AMTP,3)
      MTMT=AMTP/(1.+4*ASMT/3./PI+(16.11-1.04*(5.-6.63/AMTP))*
     $(ASMT/PI)**2)
      MTQ=SSMQCD(DBLE(MTMT),DBLE(HIGFRZ))

      BETA=ATAN(1.0/RV2V1)
      ALPHA=ALFAH
      SW2=SN2THW
      CW2=1.-SN2THW
C
      LGTST=(1+(AMTLSS/MTQ)**2)*(1+(AMTRSS/MTQ)**2)
      LGTST=LOG(LGTST)
      SBMA=SIN(BETA-ALPHA)
      CBMA=COS(BETA-ALPHA)
      CAB2=(DCOS(ALPHA+BETA))**2             
      SAB2=1.0-CAB2
C
C          hl0 -> ha0 + ha0
C
      IF(AMHL.GT.2*AMHA) THEN
        LAMB1=AMHL**2-4.0*AMHA**2
        DWID=SBMA*COS(2.0*BETA) 
C          Now add hl-hp-hp vertex correction
        DELLPP=3.0*G2*CW2*(MTQ**4)*COS(ALPHA)
        DELLPP=DELLPP*(COS(BETA)**2)/(16.0*(PI**2))
        DELLPP=DELLPP/((AMW**4)*(SIN(BETA))**3)
        DELLPP=DELLPP*LGTST
        DWID=(DWID+DELLPP)**2
        DWID=DWID*G2*(AMZ**2)/(128.0*PI*CW2*(AMHL**2))
        DWID=DWID*SQRT(LAMB1)
        WID=DWID
        CALL SSSAVE(ISHL,WID,ISHA,ISHA,0,0,0)
      ENDIF
C
C          hh -> ha + z
C
      IF(AMHH.GT.AMHA+AMZ) THEN
        MH=AMHH
        M1=AMHA
        M2=AMZ
        LAMB=SSDLAM(MH**2,M1**2,M2**2)
        DWID=SQRT(G2*CW2)+SQRT(GP2*SW2)
        DWID=DWID**2*SAB2*SQRT(LAMB)
        DWID=DWID/(64.0*PI*(AMZ**2)*(AMHH**3))
        DWID=DWID*LAMB
        WID=DWID
        CALL SSSAVE(ISHH,WID,ISHA,IDZ,0,0,0)
      ENDIF
C
C          hh -> hl + hl
C
      IF(AMHH.GT.2*AMHL) THEN
        LAMB1=AMHH**2-4.0*AMHL**2
        TEMP=CBMA*COS(2.0*ALPHA)
        TEMP=TEMP+2.0*SBMA*SIN(2.0*ALPHA)
C
C          Now add hh-hl-hl vertex correction
C
C        The following 8 lines calculate the radiative
C        hh-hl-hl vertex correction including only
C        effects from tops and stop squarks.
C
C        DTEMPL=3.0*LOG(1.0+(AMTLSS/MTQ)**2)
C        DTEMPL=DTEMPL-2.0*AMTLSS**2/(AMTLSS**2+MTQ**2)
C        DTEMPR=3.0*LOG(1.0+(AMTRSS/MTQ)**2)
C        DTEMPR=DTEMPR-2.0*AMTRSS**2/(AMTRSS**2+MTQ**2)
C        DELHLL=3.0*G2*CW2*(MTQ**4)*SIN(ALPHA)
C        DELHLL=DELHLL*(COS(ALPHA)**2)/(PI**2)
C        DELHLL=DELHLL/(16.0*(AMW**4)*(SIN(BETA))**3)
C        DELHLL=DELHLL*(DTEMPL+DTEMPR)                  
C
C        The subroutine SSHL calculates the radiative
C        hh-hl-hl vertex correction including both 
C        top & bottom and stop and sbottom squark
C        (non-degenerate with mixing) effects.  
C        A-terms and mu=-2m1 are also included.
C
        CALL SSDHLL(DELHLL)
C
C        Note:  the variable TEMP in the line below 
C        this is the Lagrangian term (as noted on 
C        page 27 of Prof. Tata's personal Lagrangian
C        term notes.  Thus DELHLL must also be the 
C        Lagrangian entry - not the potential entry.
C        The subroutine SSHLL IS set up to yield the
C        the Lagrangian entry. (We must be very careful
C        about the relative sign between TEMP and DELHLL.)
C 
        DWID=G2*(AMZ**2)*(TEMP+DELHLL)**2
        DWID=DWID/(128.0*PI*CW2*(AMHH**2))
        DWID=DWID*SQRT(LAMB1)
        WID=DWID
        CALL SSSAVE(ISHH,WID,ISHL,ISHL,0,0,0)
      ENDIF
C
C          hh -> ha + ha
C
      IF(AMHH.GT.2*AMHA) THEN
        LAMB1=AMHH**2-4.0*AMHA**2
        DWID=CBMA*COS(2*BETA)
C          Now add hh-hp-hp vertex correction
        DELHPP=3.0*G2*CW2*(MTQ**4)*SIN(ALPHA)
        DELHPP=DELHPP*(COS(BETA)**2)/(16.0*(PI**2))
        DELHPP=DELHPP/((AMW**4)*(SIN(BETA))**3)
        DELHPP=DELHPP*LGTST
        DWID=G2*(AMZ**2)*(DWID+DELHPP)**2
        DWID=DWID/(128.0*PI*CW2*(AMHH**2))
        DWID=DWID*SQRT(LAMB1)
        WID=DWID
        CALL SSSAVE(ISHH,WID,ISHA,ISHA,0,0,0)
      ENDIF
C
C          hh -> hc+ + hc-
C
      IF(AMHH.GT.2*AMHC) THEN
        LAMB1=1.0-4.0*(AMHC**2)/(AMHH**2)
        DWID=CBMA*COS(2.0*BETA)/(2.0*CW2)
        DWID=COS(BETA+ALPHA)-DWID                   
C          Now add hh-hc-hc vertex correction
        DELHCC=3.0*G2*MTQ**4*SIN(ALPHA)
        DELHCC=DELHCC/( SIN(BETA)*(DTAN(BETA))**2 )
        DELHCC=DELHCC/(32.0*PI**2*AMW**4)
        DELHCC=DELHCC*LGTST
        DWID=G2*AMW**2*(-DWID+DELHCC)**2
        DWID=DWID*SQRT(LAMB1)/(16.0*PI*AMHH)
        WID=DWID
        CALL SSSAVE(ISHH,WID,ISHC,-ISHC,0,0,0)
      ENDIF
C
C          ha -> hl + z
C
      IF(AMHA.GT.AMHL+AMZ) THEN
         MH=AMHA
         M1=AMHL
         M2=AMZ
         LAMB=SSDLAM(MH**2,M1**2,M2**2)
         DWID=SQRT(G2*CW2)+SQRT(GP2*SW2)
         DWID=DWID**2*CAB2*SQRT(LAMB)
         DWID=DWID/(64.0*PI*(AMZ**2)*(AMHA**3))
         DWID=DWID*LAMB
         WID=DWID
         CALL SSSAVE(ISHA,WID,ISHL,IDZ,0,0,0)
      ENDIF
C
C          hc+ -> w+ + hl
C
      IF(AMHC.GT.AMW+AMHL) THEN
        MH=AMHC
        M1=AMW
        M2=AMHL
        LAMB=SSDLAM(MH**2,M1**2,M2**2)
        DWID=G2*CAB2*SQRT(LAMB)
        DWID=DWID/( 64.0*PI*(AMW**2)*(AMHC**3) )
        DWID=DWID*LAMB
        WID=DWID
        CALL SSSAVE(ISHC,WID,ISHL,IDW,0,0,0)
      ENDIF
C
      RETURN
      END
