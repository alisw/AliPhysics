#ifndef FSTACK_H
#define FSTACK_H 1

#include "cfortran.h"
#include "Rtypes.h"

#include "Fdimpar.h"

extern "C" {
//*$ create stack.add
//*copy stack
//*
//*=== stack ============================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     include file: stack copy                    created 26/11/86 by p*
//*                                                                      *
//*     changes: last change on 16-sep-1999   by    alfredo ferrari      *
//*                                                  infn, milan         *
//*                                                                      *
//*     included in the following subroutines or functions: not updated  *
//*                                                                      *
//*            aemshh                                                    *
//*            beamdv                                                    *
//*            beamem                                                    *
//*            beamso                                                    *
//*            delthr                                                    *
//*            dplstk                                                    *
//*            epilog                                                    *
//*            feeder                                                    *
//*            flnwst                                                    *
//*            flukam                                                    *
//*            geofar                                                    *
//*            geomtr                                                    *
//*            kashea                                                    *
//*            kaskad                                                    *
//*            kasneu                                                    *
//*            kasray                                                    *
//*            mgdraw                                                    *
//*            pphnev                                                    *
//*            soevsv                                                    *
//*            source                                                    *
//*            stckad                                                    *
//*            stuprf                                                    *
//*            zeroin                                                    *
//*                                                                      *
//*     description of the common block(s) and variable(s)               *
//*                                                                      *
//*     /stack/ stack for the primaries                                  *
//*        wt     = weight of the particle                               *
//*        pmom   = laboratory momentum of the particle in gev/c         *
//*        tke    = laboratory kinetic energy of the particle in gev     *
//*        xa     = x-coordinate of the particle                         *
//*        ya     = y-coordinate of the particle                         *
//*        za     = z-coordinate of the particle                         *
//*        tx     = direction cosine of the particle                     *
//*                 with respect to x-axis                               *
//*        ty     = direction cosine of the particle                     *
//*                 with respect to y-axis                               *
//*        tz     = direction cosine of the particle                     *
//*                 with respect to z-axis                               *
//*        txpol  = direction cosine of the particle polarization        *
//*        typol  = direction cosine of the particle polarization        *
//*        tzpol  = direction cosine of the particle polarization        *
//*        txnor  = direction cosine of a (possible) surface normal      *
//*        tynor  = direction cosine of a (possible) surface normal      *
//*        tznor  = direction cosine of a (possible) surface normal      *
//*        dfnear = distance to the nearest boundary                     *
//*        agestk = age of the particle (seconds)                        *
//*        aknshr = kshort component of k0/k0bar                         *
//*        raddly = delay (s) in production wrt the nominal primary "0"  *
//*                 time for particle produced in radioactive decays     *
//*                (i.e. those coming from decays of daughter isotopes)  *
//*        cmpath = cumulative path travelled by the particle since it   *
//*                 was produced (cm)
//*        sparek = spare real variables available for k.w.burn          *
//*        ispark = spare integer variables available for k.w.burn       *
//*        ilo    = type of the particle (see btype in /paprop/)         *
//*        igroup = energy group for low energy neutrons                 *
//*        lo     = generation of the particle                           *
//*        louse  = user flag                                            *
//*        nreg   = number of the region of the particle                 *
//*        nlattc = number of the lattice cell of the particle           *
//*        nhspnt = pointer to the history object (geant4 geometry)      *
//*        nevent = number of the event which created the particle       *
//*        numpar = particle number                                      *
//*        lraddc = flag for particles generated in radioactive decyas   *
//*        nparma = biggest particle number encountered                  *
//*        mstack = size of the stack                                    *
//*        lstmax = highest value of the stack pointer encountered       *
//*                 in the run                                           *
//*        lstack = stack pointer                                        *
//*        lstaol = stack pointer of the last processed particle         *
//*        igroun = energy group number of the last processed particle   *
//*                 if it is a low energy neutron                        *
//*                                                                      *
//*----------------------------------------------------------------------*
//*

typedef struct {
   Double_t xa[mfstck+1];             //(0:MFSTCK)
   Double_t ya[mfstck+1];             //(0:MFSTCK)
   Double_t za[mfstck+1];             //(0:MFSTCK)
   Double_t tx[mfstck+1];             //(0:MFSTCK)
   Double_t ty[mfstck+1];             //(0:MFSTCK)
   Double_t tz[mfstck+1];             //(0:MFSTCK)
   Double_t txpol[mfstck+1];          //(0:MFSTCK)
   Double_t typol[mfstck+1];          //(0:MFSTCK)
   Double_t tzpol[mfstck+1];          //(0:MFSTCK)
   Double_t txnor[mfstck+1];          //(0:MFSTCK)
   Double_t tynor[mfstck+1];          //(0:MFSTCK)
   Double_t tznor[mfstck+1];          //(0:MFSTCK)
   Double_t wt[mfstck+1];             //(0:MFSTCK)
   Double_t pmom[mfstck+1];           //(0:MFSTCK)
   Double_t tke[mfstck+1];            //(0:MFSTCK)
   Double_t dfnear[mfstck+1];         //(0:MFSTCK)
   Double_t agestk[mfstck+1];         //(0:MFSTCK)
   Double_t aknshr[mfstck+1];         //(0:MFSTCK)
   Double_t raddly[mfstck+1];         //(0:MFSTCK)
   Double_t cmpath[mfstck+1];         //(0:MFSTCK)
   Double_t sparek[mfstck+1][mkbmx1]; //(MKBMX1,0:MFSTCK)
   Int_t    ispark[mfstck+1][mkbmx2]; //(MKBMX2,0:MFSTCK)
   Int_t    ilo[mfstck+1];            //(0:MFSTCK)
   Int_t    igroup[mfstck+1];         //(0:MFSTCK)
   Int_t    lo[mfstck+1];             //(0:MFSTCK)
   Int_t    louse[mfstck+1];          //(0:MFSTCK)
   Int_t    nreg[mfstck+1];           //(0:MFSTCK)
   Int_t    nlattc[mfstck+1];         //(0:MFSTCK)
   Int_t    nhspnt[mfstck+1];         //(0:MFSTCK)
   Int_t    nevent[mfstck+1];         //(0:MFSTCK)
   Int_t    numpar[mfstck+1];         //(0:MFSTCK)
   Int_t    lraddc[mfstck+1];         //(0:MFSTCK)
   Int_t    nparma;
   Int_t    mstack;
   Int_t    lstmax;
   Int_t    lstack;
   Int_t    lstaol;
   Int_t    igroun;
} stackCommon;
#define STACK COMMON_BLOCK(STACK,stack)
COMMON_BLOCK_DEF(stackCommon,STACK);
}

#endif
