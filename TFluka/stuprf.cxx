#include <Riostream.h>
#include "AliRun.h"
#include "TFluka.h"
#ifndef WIN32
# define stuprf stuprf_
#else
# define stuprf STUPRF
#endif
//
// Fluka include
#include "Fdimpar.h"  //(DIMPAR) fluka include
// Fluka commons
#include "Fdblprc.h"  //(DBLPRC) fluka common
#include "Fevtflg.h"  //(EVTFLG) fluka common
#include "Fpaprop.h"  //(PAPROP) fluka common
#include "Fstack.h"   //(STACK)  fluka common
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Ffinuc.h"   //(FINUC)  fluka common

//Virtual MC
#include "TFluka.h"
#include "TVirtualMCStack.h"
#include "TVirtualMCApplication.h"
#include "TParticle.h"
#include "TVector3.h"

extern "C" {
void stuprf(Int_t& ij, Int_t& mreg,
            Double_t& xx, Double_t& yy, Double_t& zz,
	    Int_t& numsec, Int_t& npprmr)
{
//*----------------------------------------------------------------------*
//*                                                                      *
//*  SeT User PRoperties for Fluka particles                             *
//*                                                                      *
//*----------------------------------------------------------------------*

// STACK.lstack  = stack pointer
// STACK.louse   = user flag
// TRACKR.llouse = user defined flag for the current particle
  STACK.louse[STACK.lstack] = TRACKR.llouse;

// mkbmx1 = dimension for kwb real spare array in fluka stack in DIMPAR
// mkbmx2 = dimension for kwb int. spare array in fluka stack in DIMPAR
// STACK.sparek  = spare real variables available for k.w.burn
// STACK.ispark  = spare integer variables available for k.w.burn
// TRACKR.spausr = user defined spare variables for the current particle
// TRACKR.ispusr = user defined spare flags for the current particle
  Int_t ispr;
  for (ispr=0; ispr<=mkbmx1-1; ispr++) {
    STACK.sparek[STACK.lstack][ispr] = TRACKR.spausr[ispr];
  }  
  for (ispr=0; ispr<=mkbmx2-1; ispr++) {
    STACK.ispark[STACK.lstack][ispr] = TRACKR.ispusr[ispr];
  }  
 
// Get the pointer to the VMC
  TVirtualMC* fluka = TFluka::GetMC();
// Get the stack produced from the generator
  TVirtualMCStack* cppstack = fluka->GetStack();
  
// EVTFLG.ntrcks = track number
// Increment the track number and put it into the last flag
  if (numsec-1 > npprmr) {
// Now call the SetTrack(...)
    Int_t done = 0;
    Int_t parent = TRACKR.ispusr[mkbmx2-1];
    Int_t pdg = fluka->PDGFromId(ij);
    
    Double_t px = FINUC.plr[numsec-1]*FINUC.cxr[numsec-1];
    Double_t pz = FINUC.plr[numsec-1]*FINUC.cyr[numsec-1];
    Double_t py = FINUC.plr[numsec-1]*FINUC.czr[numsec-1];
    Double_t e  = FINUC.tki[numsec-1] + PAPROP.am[FINUC.kpart[numsec-1]+6];
    Double_t vx = xx;
    Double_t vy = yy;
    Double_t vz = zz;
    Double_t tof = TRACKR.atrack;
    Double_t polx = FINUC.cxrpol[numsec-1];
    Double_t poly = FINUC.cyrpol[numsec-1];
    Double_t polz = FINUC.czrpol[numsec-1];

    TMCProcess mech = kPHadronic;
    if (EVTFLG.ldecay == 1) mech = kPDecay;
    else if (EVTFLG.ldltry == 1) mech = kPDeltaRay;
    else if (EVTFLG.lpairp == 1) mech = kPPair;
    else if (EVTFLG.lbrmsp == 1) mech = kPBrem;
    Double_t weight = FINUC.wei[numsec-1];
    Int_t is = 0;
    Int_t ntr;  

//virtual void SetTrack(Int_t done, Int_t parent, Int_t pdg,
//Double_t px, Double_t py, Double_t pz, Double_t e,
//Double_t vx, Double_t vy, Double_t vz, Double_t tof,
//Double_t polx, Double_t poly, Double_t polz,
//TMCProcess mech, Int_t& ntr, Double_t weight,
//Int_t is) = 0;

    
    cppstack->SetTrack(done, parent, pdg,
		    px, py, pz, e,
		    vx, vy, vz, tof,
		    polx, poly, polz,
		    mech, ntr, weight, is);

cout << endl << " !!! stuprf: ntr=" << ntr << endl;
    EVTFLG.ntrcks = ntr;
    STACK.ispark[STACK.lstack][mkbmx2-1] = EVTFLG.ntrcks;
  } // end of if (numsec-1 > npprmr)

} // end of stuprf
} // end of extern "C"

