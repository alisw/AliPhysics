#include <Riostream.h>
#include "AliRun.h"
#include "TFluka.h"
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

//Virtual MC
#include "TFluka.h"
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

  Int_t lbhabh = 0;
  if (EVTFLG.ldltry == 1) {
    if (EMFSTK.iq[EMFSTK.np-1] * EMFSTK.iq[EMFSTK.np-2] < 0) lbhabh = 1;
  }

// mkbmx1 = dimension for kwb real spare array in fluka stack in DIMPAR
// mkbmx2 = dimension for kwb int. spare array in fluka stack in DIMPAR
// EMFSTK.espark  = spare real variables available for 
// EMFSTK.iespak  = spare integer variables available for
// TRACKR.spausr = user defined spare variables for the current particle
// TRACKR.ispusr = user defined spare flags for the current particle
// EMFSTK.louemf = user flag
// TRACKR.llouse = user defined flag for the current particle

  Int_t npnw, ispr;
  for (npnw=EMFSTK.npstrt-1; npnw<=EMFSTK.np-1; npnw++) {

    for (ispr=0; ispr<=mkbmx1-1; ispr++) 
      EMFSTK.espark[npnw][ispr] = TRACKR.spausr[ispr];

    for (ispr=0; ispr<=mkbmx2-1; ispr++) 
      EMFSTK.iespak[npnw][ispr] = TRACKR.ispusr[ispr];

    EMFSTK.louemf[npnw] = TRACKR.llouse;
  }

// Get the pointer to the VMC
  TVirtualMC* fluka = TFluka::GetMC();
// Get the stack produced from the generator
  TVirtualMCStack* cppstack = fluka->GetStack();
  
// EVTFLG.ntrcks = track number
// Increment the track number and put it into the last flag

  Int_t kp;
  for (kp=EMFSTK.npstrt-1; npnw<=EMFSTK.np-1; kp++) {

//* save the parent track number and reset it at each loop
    Int_t ncrtrk = EVTFLG.ntrcks;
    
    Int_t done = 0;
    Int_t parent = ncrtrk;
    Int_t flukaid = 0;
    if (EMFSTK.iq[kp] == -1) flukaid = 3;
    else if (EMFSTK.iq[kp] == 0)  flukaid = 7;
    else if (EMFSTK.iq[kp] == 0)  flukaid = 4;
    Int_t pdg = fluka->PDGFromId(flukaid);
    Double_t e  = EMFSTK.e[kp]*emvgev;
    Double_t p = sqrt(e*e - PAPROP.am[flukaid+6] * PAPROP.am[flukaid+6]);
    Double_t px = p * EMFSTK.u[kp];
    Double_t pz = p * EMFSTK.v[kp];
    Double_t py = p * EMFSTK.w[kp];
    Double_t tof = EMFSTK.agemf[kp];
    Double_t polx = EMFSTK.upol[kp];
    Double_t poly = EMFSTK.vpol[kp];
    Double_t polz = EMFSTK.wpol[kp];
    Double_t vx = EMFSTK.x[kp];
    Double_t vy = EMFSTK.y[kp];
    Double_t vz = EMFSTK.z[kp];
    Double_t weight = EMFSTK.wt[kp];
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
        cppstack->SetTrack(done, parent, pdg,
           px, py, pz, e, vx, vy, vz, tof,
           polx, poly, polz, mech, ntr, weight, is);

cout << endl << " !!! stupre: ntr=" << ntr << endl;
      EMFSTK.iespak[kp][mkbmx2-1] = ntr;
    } // end of lpairp, lphoel, lannfl, lannrs
    
//* Compton: secondary is true only if charged (e+, e-)
    else if ((EVTFLG.lcmptn == 1)) {
      if (EMFSTK.iq[kp] != 0) {
        mech = kPCompton;
        cppstack->SetTrack(done, parent, pdg,
           px, py, pz, e, vx, vy, vz, tof,
           polx, poly, polz, mech, ntr, weight, is);
cout << endl << " !!! stupre: ntr=" << ntr << endl;
        EMFSTK.iespak[kp][mkbmx2-1] = ntr;
      }
    } // end of lcmptn

//* Bremsstrahlung: true secondary only if charge = 0 (photon)
    else if ((EVTFLG.lbrmsp == 1)) {
      if (EMFSTK.iq[kp] == 0) {
        mech = kPBrem;
        cppstack->SetTrack(done, parent, pdg,
           px, py, pz, e, vx, vy, vz, tof,
           polx, poly, polz, mech, ntr, weight, is);
cout << endl << " !!! stupre: ntr=" << ntr << endl;
        EMFSTK.iespak[kp][mkbmx2-1] = ntr;
      }
    } // end of lbrmsp

//* Delta ray: If Bhabha, true secondary only if negative (electron)
    else if ((EVTFLG.ldltry == 1)) {
      if (lbhabh == 1) {
        if (EMFSTK.iq[kp] == 0) {
          mech = kPDeltaRay;
          cppstack->SetTrack(done, parent, pdg,
            px, py, pz, e, vx, vy, vz, tof,
            polx, poly, polz, mech, ntr, weight, is);
          EMFSTK.iespak[kp][mkbmx2-1] = ntr;
        } // end of Bhabha
      }

//* Delta ray: Otherwise Moller: true secondary is the electron with
//*            lower energy, which has been put higher in the stack
      else if (kp == EMFSTK.np) {
        mech = kPDeltaRay;
        cppstack->SetTrack(done, parent, pdg,
          px, py, pz, e, vx, vy, vz, tof,
          polx, poly, polz, mech, ntr, weight, is);
cout << endl << " !!! stupre: ntr=" << ntr << endl;
        EMFSTK.iespak[kp][mkbmx2-1] = ntr;
      } // end of Delta ray
    } // end of ldltry

  } // end of loop

// !!! TO BE CONFIRMED !!!
  EVTFLG.ntrcks = EMFSTK.iespak[EMFSTK.np-1][mkbmx2-1];

} // end of stupre
} // end of extern "C"

