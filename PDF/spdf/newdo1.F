*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:48  cernlib
* Version 7.04
*
* Revision 1.2  1996/10/28 15:03:18  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:25  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
              SUBROUTINE NEWDO1(R8X,R8QSQ,R8U,R8D,R8SS,R8STR,R8C,R8G)
C--------------------------------------------------------------
C       THIS SUBROUTINE EVALUATES THE UPDATED DUKE-OWENS QUARK
C       DISTRIBUTIONS , SET 1--------------U. BAUR (JUNE 9, 1991)
C       adapted and corrected by---H.Plothow-Besch (Aug 14, 1991)
C       BASED ON FSU-HEP-910606 PREPRINT
C--------------------------------------------------------------
C
#include "pdf/expdp.inc"
     +       R8X,R8QSQ,R8U,R8D,R8SS,R8STR,R8C,R8G
        REAL ND,NUD,X,QSQ
        DATA Q0/2./, ALAM/0.177/
        X = R8X
        QSQ = R8QSQ
        S= LOG( LOG( MAX(Q0,QSQ)/ALAM)/ LOG(Q0/ALAM))
        S2=S*S
        S3=S2*S
        ONEMX = 1. - X
C-------SEA DIST-----------------------------------------------
        AS = 0.9090 - 0.4023*S + 0.006305*S2
        AS1 = -0.3823*S + 0.02766*S2
        BS1 = 7.278 - 0.7904*S + 0.8108*S2
        ALS = -1.6629*S + 0.5719*S2
        BES = -0.01333*S + 0.5299*S2
        GAS = 0.1211*S - 0.1739*S2
C       SS=AS*X**(AS1-1.)*ONEMX**BS1*(1.+ALS*X+BES*X*X+GAS*X**3)/6.
        SS = AS * X**AS1 * ONEMX**BS1 * (1.+ALS*X+BES*X*X+GAS*X*X*X)/ 6.
C-------D-DIST-----------------------------------------
        H3 = 0.8388 - 0.2092*S + 0.02657*S2
C       H4 = 4.6670 + 0.7951*S - 0.1081*S2
        H4 = 4.6670 + 0.7951*S + 0.1081*S2
        GD1 = -1.0232*S + 0.05799*S2
        GD2 = 0.8616*S + 0.1530*S2
        CALL GAMMADO(H3,G3,IER)
        CALL GAMMADO(1.+H3,G31,IER)
        CALL GAMMADO(2.+H3,G32,IER)
        CALL GAMMADO(1.+H4,G4,IER)
        CALL GAMMADO(1.+H3+H4,G34,IER)
        CALL GAMMADO(2.+H3+H4,G341,IER)
        CALL GAMMADO(3.+H3+H4,G342,IER)
        BD1=G3*G4/G34
        BD2=G31*G4/G341
        BD3=G32*G4/G342
        ND=BD1+GD1*BD2+GD2*BD3
C       D=X**(H3-1.)*(1.-X)**H4*(1.+GD1*X+GD2*X**2)/ND
        D=X**H3*(1.-X)**H4*(1.+GD1*X+GD2*X**2)/ND
C-------VALENCE U-DIST-----------------------------------------
        H1 = 0.6650 - 0.1097*S - 0.002442*S2
        H2 = 3.6140 + 0.8395*S - 0.02186*S2
        GUD1 = 0.8673 - 1.6637*S + 0.3420*S2
C       GUD2 = 1.1049*S - 0.22369*S2
        GUD2 = 1.1049*S - 0.2369*S2
        CALL GAMMADO(H1,G1,IER)
        CALL GAMMADO(1.+H1,G11,IER)
        CALL GAMMADO(2.+H1,G21,IER)
        CALL GAMMADO(1.+H2,G2,IER)
        CALL GAMMADO(1.+H1+H2,G12,IER)
        CALL GAMMADO(2.+H1+H2,G121,IER)
        CALL GAMMADO(3.+H1+H2,G122,IER)
        BUD1=G1*G2/G12
        BUD2=G11*G2/G121
        BUD3=G21*G2/G122
        NUD=BUD1+GUD1*BUD2+GUD2*BUD3
C       UPD=3*X**(H1-1.)*(1.-X)**H2*(1.+GUD1*X+GUD2*X**2)/NUD
        UPD=3. * X**H1 * ONEMX**H2 * (1.+GUD1*X+GUD2*X*X)/ NUD
C       U=UPD-D+SS
C       D=D+SS
        U=UPD-D
C-------CHARM DIST---------------------------------------------
        AC = 0.09469*S - 0.07066*S2 + 0.01236*S3
        AC1 = -0.1447 - 0.4020*S + 0.1533*S2 - 0.06479*S3
        BC1 = 6.7599 + 1.6596*S + 0.6798*S2 - 0.8525*S3
        ALC = -4.4559*S + 3.3756*S2 - 0.9468*S3
        BEC = 7.862*S - 3.6591*S2 + 0.03672*S3
        GAC = -0.2472*S - 0.7510*S2 + 0.04870*S3
C       C=AC*X**(AC1-1.)*(1.-X)**BC1*(1.+ALC*X+BEC*X*X+GAC*X**3)
        C=AC * X**AC1 * ONEMX**BC1 * (1.+ALC*X+BEC*X*X+GAC*X*X*X)
C-------GLUE DIST----------------------------------------------
        AG = 3.0170 - 4.7347*S + 3.3594*S2 - 0.9443*S3
        AG1 = -0.9342*S + 0.5454*S2 - 0.1668*S3
        BG1 = 5.3040 + 1.4654*S - 1.4292*S2 + 0.7569*S3
        ALG = -3.9141*S + 2.8445*S2 - 0.8411*S3
        BEG = 9.0176*S - 10.426*S2 + 4.0983*S3
        GAG = -5.9602*S + 7.5150*S2 - 2.7329*S3
C       G=AG*X**(AG1-1.)*(1.-X)**BG1*(1.+ALG*X+BEG*X*X+GAG*X**3)
        G=AG * X**AG1 * ONEMX**BG1 * (1.+ALG*X+BEG*X*X+GAG*X*X*X)
        R8U = U
        R8D = D
        R8SS = SS
        R8STR = SS
        R8C = C
        R8G = G
        RETURN
        END
