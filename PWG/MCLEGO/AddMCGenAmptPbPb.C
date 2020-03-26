R__LOAD_LIBRARY(libAMPT)
R__LOAD_LIBRARY(libTAmpt)
R__LOAD_LIBRARY(libEGPythia6)
R__LOAD_LIBRARY(libpythia6)
R__LOAD_LIBRARY(libAliPythia6)

#include "AliGenerator.h"
#include "AliGenAmpt.h"

AliGenerator *AddMCGenAmptPbPb(
                          Double_t Energy      = 5020.,   // CM energy
                          Double_t bmin        = 0.0,     // minimum impact parameter
                          Double_t bmax        = 20.0,    // maximum impact parameter
                          Double_t parA        = 0.3,     // parameter a in Lund symmetric splitting function
                          Double_t parB        = 0.4,     // parameter b in Lund symmetric splitting function
                          Double_t ptHardMin   = 3.0,     // minimum pt hard (was 3.0 in previous AMPT productions)
                          Double_t mu          = 2.2814,  // parton screening mass in fm^(-1)
                          Double_t alpha_s     = 1./3.,   // alpa in parton cascade
                          Int_t Flag_SM        = 4,       // flag for string melting
                          Int_t NTmax          = 3,     // NTMAX: number of timesteps
                          Bool_t stringMelting = kTRUE   // string melting option
                           )
{
  AliGenAmpt *genAMPT = new AliGenAmpt(-1);
  //=========================================================================

  // Decayer
  AliDecayer *decayer = new AliDecayerPythia();
  genAMPT->SetForceDecay( kHadronicD );
  genAMPT->SetDecayer( decayer );
  //=========================================================================

  // Collision system
  genAMPT->SetEnergyCMS(Energy);
  genAMPT->SetReferenceFrame("CMS");
  genAMPT->SetProjectile("A", 208, 82);
  genAMPT->SetTarget("A", 208, 82);
  genAMPT->SetPtHardMin (ptHardMin);
  genAMPT->SetImpactParameterRange(bmin,bmax);
  //=========================================================================

  // options
  genAMPT->SetJetQuenching(0);     // enable jet quenching
  genAMPT->SetShadowing(1);        // enable shadowing
  genAMPT->SetDecaysOff(1);        // neutral pion and heavy particle decays switched off
  genAMPT->SetSpectators(0);       // track spectators
  genAMPT->SetIsoft(Flag_SM);      // 4=string melting, 1=standard AMPT
  genAMPT->SetXmu(mu);             // parton xsection
  genAMPT->SetNtMax(NTmax);        // time bins
  genAMPT->SetAlpha(alpha_s);      // alpha =0.333
  genAMPT->SetStringFrag(parA,parB); // string fragmentation parameters
  genAMPT->SetIpop(1);             // enable popcorn mechanism (net-baryon stopping)
  //=========================================================================

  // Boost into LHC lab frame
  genAMPT->SetBoostLHC(1);

  // randomize reaction plane
  genAMPT->SetRandomReactionPlane(kTRUE);

 return genAMPT;

}
