*CMZ :          17/07/98  15.54.14  by  Federico Carminati
*-- Author :
      SUBROUTINE SHSTAT
c	=================

c	Print Run Statistics

*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,SHRUNP.
      COMMON /SHRUNP/ VMAJ,IMIN,NRUN,NEVTOT
*KEEP,SHPHYP.
      COMMON /SHPHYP/ JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,
     +                  JKAC,JKA0,JRHO,JOME,JPHI,JPSI,JDRY
*KEEP,SHRATS.
      COMMON /SHRATS/ RETAPI,RPROPI,RKACPI,RRHOPI,ROMEPI,RPHIPI
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEND.

      WRITE (MSTU(11),1001) VMAJ,IMIN,NRUN,NEVTOT,IEVT
     +        ,JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,JKAC,JKA0
     +        ,JRHO,JOME,JPHI,JPSI
      WRITE (MSTU(11),1002)        JDRY
     +        ,RETAPI,RPROPI,RKACPI
     +        ,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,NPHI,NPSI,NDRY
      RETURN

1001	FORMAT(1H1,///,10X,'SHAKER Major Version ',f2.1,
     +  ' Minor Version ',i2,/,10x,45(1h-),
     +  //,5x,'Run Number ...................................... ',i8,
     +   /,5x,'Total Number of Requested Events ................ ',i8,
     +   /,5x,'Total Number of Generated Events ................ ',i8,
     +   /,5x,'Weighted Generation (pi0 and eta only) .......... ',i8,
     +   /,5x,'Rapidity Density ................................ ',i8,
     +   /,5x,'Rapidity Limit .................................. ',f6.2,
     +   /,5x,'pT Limit ........................................ ',f6.2,
     +   /,5x,'Weak Decays Enable .............................. ',i8,
     +   /,5x,'pi0    Generation ............................... ',i8,
     +   /,5x,'eta    Generation ............................... ',i8,
     +   /,5x,'pi+/-  Generation ............................... ',i8,
     +   /,5x,'proton Generation ............................... ',i8,
     +   /,5x,'K+/-   Generation ............................... ',i8,
     +   /,5x,'K0     Generation ............................... ',i8,
     +   /,5x,'rho    Generation ............................... ',i8,
     +   /,5x,'omega  Generation ............................... ',i8,
     +   /,5x,'phi    Generation ............................... ',i8,
     +   /,5x,'J/psi  Generation ............................... ',i8)
 1002   FORMAT(
     +     5x,'Dr-Yan Generation ............................... ',i8,
     +   /,5x,'eta/pi0 Ratio ................................... ',f6.2,
     +   /,5x,'p/pi    Ratio ................................... ',f6.2,
     +   /,5x,'K/pi    Ratio ................................... ',f6.2,
     +   /,5x,'pi0    / Event .................................. ',i8,
     +   /,5x,'eta    / Event .................................. ',i8,
     +   /,5x,'pi+/-  / Event .................................. ',i8,
     +   /,5x,'proton / Event .................................. ',i8,
     +   /,5x,'K+/-   / Event .................................. ',i8,
     +   /,5x,'K0     / Event .................................. ',i8,
     +   /,5x,'rho    / Event .................................. ',i8,
     +   /,5x,'omega  / Event .................................. ',i8,
     +   /,5x,'phi    / Event .................................. ',i8,
     +   /,5x,'J/psi  / Event .................................. ',i8,
     +   /,5x,'Dr-Yan / Event .................................. ',i8)

      END
