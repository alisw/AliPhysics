AliGenerator *AddMCGenAmpt(
                          Double_t Energy      = 2760.,   // CM energy
                          Double_t bmin        = 0.0,     // minimum impact parameter
                          Double_t bmax        = 20.0,    // maximum impact parameter
                          Double_t ptHardMin   = 3.0,     // minimum pt hard (was 3.0 in previous AMPT productions)
                          Bool_t stringMelting = kTRUE,   // string melting option
                          Bool_t useART        = kTRUE,   // use hadronic rescattering phase (ART)
                          Bool_t pAcollisions  = kFALSE,  // pA instead of AA collisions
                          Bool_t ppcollisions  = kFALSE   // pp instead of AA collisions
                           )
{
  // User defined generator

  gSystem->Load("libAMPT");
  gSystem->Load("libTAmpt");
  gSystem->Load("libEGPythia6"); 
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");


  AliGenAmpt *genAMPT = new AliGenAmpt(-1);
  //=========================================================================


  // User settings
  Int_t Flag_SM    = 4;       // flag for string melting: 1 = default, 4 = String Melting
  Int_t NTmax      = 150;     // NTMAX: number of timesteps (Default = 150), to turn off ART set it to 3
  Double_t mu      = 3.2264;  // parton screening mass in fm^(-1) (Default = 3.2264)
  Double_t alpha_s = 1./3.;   // change mu and alpha_s (Default = 1./3.) to vary scattering cross-section
                              // mu = 3.2 fm^-1 and alpha_s = 0.33 ==> sigma_{partonic} = 1.5mb
  if(!stringMelting)
    Flag_SM = 1;

  if(!useART)
    NTmax = 3;

  if(pAcollisions)
    NTmax = 1500; // this was used in earlier productions for p-Pb
  //=========================================================================


  // Decayer
  AliDecayer *decayer = new AliDecayerPythia();
  genAMPT->SetForceDecay( kHadronicD );
  genAMPT->SetDecayer( decayer );
  //=========================================================================

  // Collision system
  genAMPT->SetEnergyCMS(Energy);
  genAMPT->SetReferenceFrame("CMS");
  if(ppcollisions) {
    genAMPT->SetProjectile("P", 1, 1);
    genAMPT->SetTarget("P", 1, 1);
  }
  else if (pAcollisions){
    genAMPT->SetProjectile("A", 208, 82);
    genAMPT->SetTarget("P", 1, 1);
  } else {
    genAMPT->SetProjectile("A", 208, 82);
    genAMPT->SetTarget("A", 208, 82);
  }
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
  genAMPT->SetStringFrag(0.5,0.9); // string fragmentation parameters
  genAMPT->SetIpop(1);             // enable popcorn mechanism (net-baryon stopping)
  //=========================================================================

  // Boost into LHC lab frame
  genAMPT->SetBoostLHC(1);

  // randomize reaction plane
  genAMPT->SetRandomReactionPlane(kTRUE);

 return genAMPT;

}
