AliGenerator* AddMCEMCocktailV2(  Int_t collisionsSystem      = 200,
                                  Int_t centrality            = 0,
                                  Int_t decayMode             = 1,
                                  Int_t selectedMothers       = 62591,
                                  TString paramFile           = "",
                                  TString paramFileDir        = "",
                                  Int_t numberOfParticles     = 1000,
                                  Double_t minPt              = 0.,
                                  Double_t maxPt              = 20,
                                  Int_t pythiaErrorTolerance  = 2000,
                                  Bool_t externalDecayer      = 0,
                                  Bool_t decayLongLived       = 0,
                                  Bool_t dynamicalPtRange     = 0,
                                  Bool_t useYWeights          = 0
                                )
{
  // collisions systems defined:
  // 0   : pp 900 GeV
  // 100 : pp 2.76 TeV
  // 200 : pp 7 TeV
  // 300 : pPb 5.023 TeV
  // 400 : PbPb 2.76 TeV
  
  // load libraries
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");
  
  // Create and Initialize Generator
  AliGenEMCocktailV2 *gener     = new AliGenEMCocktailV2();
  
  //=======================================================================
  // Set External decayer
  TVirtualMCDecayer *decayer  = new AliDecayerPythia();
  if (externalDecayer) decayer->AliDecayerPythia::SetDecayerExodus();
  if (decayLongLived) decayer->AliDecayerPythia::DecayLongLivedParticles();
  
  gener->SetParametrizationFile(paramFile);
  gener->SetParametrizationFileDirectory(paramFileDir);
  gener->SetNPart(numberOfParticles);                         // source multiplicity per event
  gener->SetPtRange(minPt,maxPt);
  gener->SetDynamicalPtRange(dynamicalPtRange);
  gener->SetUseYWeighting(useYWeights);
  gener->SetYRange(-1.,1.);
  gener->SetPhiRange(0., 360.);
  gener->SetOrigin(0.,0.,0.); 
  gener->SetSigma(0.,0.,0.);
  gener->SetVertexSmear(kPerEvent);
  gener->SetTrackingFlag(0);
  gener->SelectMotherParticles(selectedMothers);
  gener->SetCollisionSystem(collisionsSystem);                //pp 7 TeV
  gener->SetCentrality(centrality);                           // kpp
  (AliPythia::Instance())->SetMSTU(22, pythiaErrorTolerance);   // tolerance for error due to rhos
  
  if (decayMode == 1){
    gener->SetDecayMode(kGammaEM);    	// kGammaEM      => single photon
  } else if (decayMode == 2){
    gener->SetDecayMode(kElectronEM);	 // kElectronEM   => single electron
  } else if (decayMode == 3){
    gener->SetDecayMode(kDiElectronEM); // kDiElectronEM => electron-positron
  }	
  gener->SetDecayer(decayer);
  gener->SetWeightingMode(kNonAnalog); 	// select weighting:
                      // kNonAnalog => weight ~ dN/dp_T
                      // kAnalog    => weight ~ 1
  gener->CreateCocktail();
  gener->Init();

  return gener;
}
