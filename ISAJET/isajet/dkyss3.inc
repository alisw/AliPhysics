#ifndef CERNLIB_ISAJET_DKYSS3_INC
#define CERNLIB_ISAJET_DKYSS3_INC
*
*
* dkyss3.inc
*
C
C          Data for SUSY 3-body matrix elements. There is a double 
C          pointer structure, first to modes, and then to poles that
C          make up the matrix element for that mode:
C          MELEM=-I in /DKYTAB/ points to the mode information:
C            J1SS3(I) = start of pole list for this mode
C            J2SS3(I) = end of pole list for this mode
C            WTSS3(I) = maximum weight for this mode
C          J1SS3<J<J2SS3 points to the corresponding poles:
C            KSS3(J)    = pole type
C            AMSS3(J)   = pole mass
C            ZISS3(2,J) = initial couplings
C            ZFSS3(2,J) = final couplings
C          For gaugino -> gaugino f fbar, the pole types are
C            KSS3=1: spin-1 pole in f-fbar channel
C            KSS3=2: spin-0 pole in gaugino-f channel
C            KSS3=3: spin-0 pole in gaugino-fbar channel
C            KSS3=4: spin-0 pole in f-fbar channel
C          The two couplings are the coefficients of 1,gamma_5 or of
C          gamma_mu,gamma_mu*gamma_5. 
C
      INTEGER MXMSS3,MXPSS3
      PARAMETER (MXMSS3=1000)
      PARAMETER (MXPSS3=2000)
      COMMON/DKYSS3/NMSS3,NPSS3,
     $J1SS3(MXMSS3),J2SS3(MXMSS3),WTSS3(MXMSS3),
     $KSS3(MXPSS3),AMSS3(MXPSS3),ZISS3(2,MXPSS3),ZFSS3(2,MXPSS3)
      INTEGER NMSS3,NPSS3,KSS3,J1SS3,J2SS3
      REAL WTSS3,AMSS3
      COMPLEX ZISS3,ZFSS3

#endif
