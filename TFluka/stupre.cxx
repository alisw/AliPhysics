#include <Riostream.h>
#ifndef WIN32
# define stupre stupre_
#else
# define stupre STUPRE
#endif
//
// Fluka include
#include "Fdimpar.h"  //(DIMPAR) fluka include

// Fluka commons
#include "Fdblprc.h"  //(DBLPRC) fluka common
#include "Femfstk.h"  //(EMFSTK) fluka common
#include "Fevtflg.h"  //(EVTFLG) fluka common
#include "Fpaprop.h"  //(PAPROP) fluka common
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Femfrgn.h"  //(EFMRGN) fluka common

//Virtual MC
#include "TFluka.h"
#include "TFlukaCodes.h"
#include "TVirtualMCStack.h"
#include "TVirtualMCApplication.h"
#include "TParticle.h"
#include "TVector3.h"

extern "C" {
void stupre()
{
//*----------------------------------------------------------------------*
//*                                                                      *
//*  SeT User PRoperties for Emf particles                               *
//*                                                                      *
//*----------------------------------------------------------------------*
// Get the pointer to the VMC
  TFluka* fluka =  (TFluka*) gMC;
  static Double_t emassmev = PAPROP.am[kFLUKAelectron + 6] * 1000.;
    
  Int_t lbhabh = 0;
  if (EVTFLG.ldltry == 1) {
    if (EMFSTK.ichemf[EMFSTK.npemf-1] * EMFSTK.ichemf[EMFSTK.npemf-2] < 0) lbhabh = 1;
  }

// mkbmx1         = dimension for kwb real spare array in fluka stack in DIMPAR
// mkbmx2         = dimension for kwb int. spare array in fluka stack in DIMPAR
// EMFSTK.espark  = spare real variables available for 
// EMFSTK.iespak  = spare integer variables available for
// TRACKR.spausr  = user defined spare variables for the current particle
// TRACKR.ispusr  = user defined spare flags for the current particle
// EMFSTK.louemf  = user flag
// TRACKR.llouse  = user defined flag for the current particle

  Int_t npnw, ispr;
  for (npnw = EMFSTK.npstrt-1; npnw <= EMFSTK.npemf-1; npnw++) {

      for (ispr = 0; ispr <= mkbmx1-1; ispr++) 
          EMFSTK.espark[npnw][ispr] = TRACKR.spausr[ispr];

      for (ispr = 0; ispr <= mkbmx2-1; ispr++) 
          EMFSTK.iespak[npnw][ispr] = TRACKR.ispusr[ispr];

      EMFSTK.louemf[npnw] = TRACKR.llouse;
  }


  Int_t verbosityLevel = fluka->GetVerbosityLevel();
  Bool_t debug = (verbosityLevel>=3)?kTRUE:kFALSE;
  fluka->SetTrackIsNew(kTRUE);

// Get the stack produced from the generator
  TVirtualMCStack* cppstack = fluka->GetStack();
  
// EVTFLG.ntrcks = track number
// Increment the track number and put it into the last flag

  Int_t kp;
  
  for (kp = EMFSTK.npstrt - 1; kp <= EMFSTK.npemf - 1; kp++) {
    
// Check transport cut first
    Int_t    ireg   = EMFSTK.iremf[kp];
    Double_t cut    = (TMath::Abs(EMFSTK.ichemf[kp]) == 1) ? EMFRGN.elethr[ireg-1] :  EMFRGN.phothr[ireg-1];
    Double_t e      = EMFSTK.etemf[kp];

    if ((e < cut) 
        && (
            (EMFSTK.ichemf[kp] ==  0) ||
            (EMFSTK.ichemf[kp] == -1) ||
            (EMFSTK.ichemf[kp] ==  1 &&  EMFRGN.phothr[ireg-1] > emassmev)
            )
        )
    {
        EMFSTK.iespak[kp][mkbmx2-1] = TRACKR.ispusr[mkbmx2-1];
        EMFSTK.iespak[kp][mkbmx2-2] =  0;
        EMFSTK.iespak[kp][mkbmx2-3] = TRACKR.jtrack;
        EMFSTK.irlatt[kp] = TRACKR.lt1trk;
        continue;
    }

// Save the parent track number and reset it at each loop
    Int_t done = 0;
    Int_t parent =  TRACKR.ispusr[mkbmx2-1];
    Int_t flukaid = 0;

// Identify particle type
    if      (EMFSTK.ichemf[kp] == -1)  flukaid = kFLUKAelectron;
    else if (EMFSTK.ichemf[kp] ==  0)  flukaid = kFLUKAphoton;
    else if (EMFSTK.ichemf[kp] ==  1)  flukaid = kFLUKApositron;


    e *= emvgev;
    Int_t    pdg    = fluka->PDGFromId(flukaid);
    Double_t p      = sqrt(e * e - PAPROP.am[flukaid+6] * PAPROP.am[flukaid+6]);
    Double_t px     = p * EMFSTK.uemf[kp];
    Double_t py     = p * EMFSTK.vemf[kp];
    Double_t pz     = p * EMFSTK.wemf[kp];
    Double_t tof    = EMFSTK.agemf[kp];
    Double_t polx   = EMFSTK.upol[kp];
    Double_t poly   = EMFSTK.vpol[kp];
    Double_t polz   = EMFSTK.wpol[kp];
    Double_t vx     = EMFSTK.xemf[kp];
    Double_t vy     = EMFSTK.yemf[kp];
    Double_t vz     = EMFSTK.zemf[kp];
    Double_t weight = EMFSTK.wtemf[kp];

    Int_t ntr;
    TMCProcess mech;
    Int_t is = 0;

//* case of no parent left (pair, photoelectric, annihilation):
//* all secondaries are true
    if ((EVTFLG.lpairp == 1) || (EVTFLG.lphoel == 1) ||
        (EVTFLG.lannfl == 1) || (EVTFLG.lannrs == 1)) {

        if (EVTFLG.lpairp == 1) mech = kPPair;
        else if (EVTFLG.lphoel == 1) mech = kPPhotoelectric;
        else mech = kPAnnihilation;
        cppstack->PushTrack(done, parent, pdg,
                           px, py, pz, e, vx, vy, vz, tof,
                           polx, poly, polz, mech, ntr, weight, is);
        if (debug) cout << endl << " !!! stupre (PAIR, ..) : ntr=" << ntr << " pdg " << pdg << " parent=" << parent << " energy " << e-PAPROP.am[flukaid+6] << endl;

        EMFSTK.iespak[kp][mkbmx2-1] = ntr;
        EMFSTK.iespak[kp][mkbmx2-2] = 0;
    } // end of lpairp, lphoel, lannfl, lannrs

//* Compton: secondary is true only if charged (e+, e-)
    else if ((EVTFLG.lcmptn == 1)) {

        if (EMFSTK.ichemf[kp] != 0) {
            mech = kPCompton;
            cppstack->PushTrack(done, parent, pdg,
                               px, py, pz, e, vx, vy, vz, tof,
                               polx, poly, polz, mech, ntr, weight, is);
            if (debug) cout << endl << " !!! stupre (COMPTON) : ntr=" << ntr << " pdg " << pdg << " parent=" << parent << endl;
            EMFSTK.iespak[kp][mkbmx2-1] = ntr;
            EMFSTK.iespak[kp][mkbmx2-2] = 0;
        }
    } // end of lcmptn

//* Bremsstrahlung: true secondary only if charge = 0 (photon)
    else if ((EVTFLG.lbrmsp == 1)) {
        if (EMFSTK.ichemf[kp] == 0) {
            mech = kPBrem;
            cppstack->PushTrack(done, parent, pdg,
                               px, py, pz, e, vx, vy, vz, tof,
                               polx, poly, polz, mech, ntr, weight, is);
            if (debug) cout << endl << " !!! stupre (BREMS) : ntr=" << ntr << " pdg " << pdg << " parent=" << parent << endl;
            EMFSTK.iespak[kp][mkbmx2-1] = ntr;
            EMFSTK.iespak[kp][mkbmx2-2] = 0;
        }
    } // end of lbrmsp

//* Delta ray: If Bhabha, true secondary only if negative (electron)
    else if ((EVTFLG.ldltry == 1)) {
        if (lbhabh == 1) {
            if (EMFSTK.ichemf[kp] == -1) {
                mech = kPDeltaRay;
                cppstack->PushTrack(done, parent, pdg,
                                   px, py, pz, e, vx, vy, vz, tof,
                                   polx, poly, polz, mech, ntr, weight, is);
                EMFSTK.iespak[kp][mkbmx2-1] = ntr;
                EMFSTK.iespak[kp][mkbmx2-2] = 0;
           if (debug) cout << endl << " !!! stupre (BHABA) : ntr=" << ntr << " pdg " << pdg << " parent=" << parent << endl;
            } // end of Bhabha
        } // lbhabh == 1

//* Delta ray: Otherwise Moller: true secondary is the electron with
//*            lower energy, which has been put higher in the stack
        else if (kp == EMFSTK.npemf-1) {
            mech = kPDeltaRay;
            cppstack->PushTrack(done, parent, pdg,
                               px, py, pz, e, vx, vy, vz, tof,
                               polx, poly, polz, mech, ntr, weight, is);
            if (debug) cout << endl << " !!! stupre (Moller) : ntr=" << ntr << " pdg " << pdg << " parent=" << parent << endl;
            EMFSTK.iespak[kp][mkbmx2-1] = ntr;
            EMFSTK.iespak[kp][mkbmx2-2] = 0;
        } // end of Delta ray
    } // end of ldltry

  } // end of loop
} // end of stupre
} // end of extern "C"

