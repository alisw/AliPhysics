*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHLIST
C	=================

c	List data output

*KEEP,SHRUNP.
      COMMON /SHRUNP/ VMAJ,IMIN,NRUN,NEVTOT
*KEEP,SHPHYP.
      COMMON /SHPHYP/ JWEI,NDNDY,YLIM,PTLIM,JWEAK,JPI0,JETA,JPIC,JPRO,
     +                  JKAC,JKA0,JRHO,JOME,JPHI,JPSI,JDRY
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.

      WRITE(MSTU(11),*) '------------------------------------------'
      WRITE(MSTU(11),*) 'RUN: ',NRUN
      WRITE(MSTU(11),*) 'EVT: ',IEVT
      CALL LULIST(2)

      RETURN
      END
