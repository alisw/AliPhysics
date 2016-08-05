AliGenerator* AddMCEMCocktail(Int_t collisionsSystem  = 200,
                              Int_t centrality        = 0,
                              Int_t decayMode         = 1,
                              Int_t selectedMothers   = 62591,
                              Int_t paramPi0          = 0,
                              Int_t paramEta          = 3,
                              Int_t paramOmega        = 3,
                              Int_t paramPhi          = 0,
                              Int_t numberOfParticles	= 1000,
                              Double_t minPt          = 0.,
                              Double_t maxPt          = 20
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
  AliGenEMCocktail *gener     = new AliGenEMCocktail();

  //=======================================================================
  // Set External decayer
  TVirtualMCDecayer *decayer  = new AliDecayerPythia();

  gener->SetNPart(numberOfParticles);               // source multiplicity per event
  gener->SetPtRange(minPt,maxPt);
  gener->SetYRange(-1.,1.);
  gener->SetPhiRange(0., 360.);
  gener->SetOrigin(0.,0.,0.); 
  gener->SetSigma(0.,0.,0.);
  gener->SetVertexSmear(kPerEvent);
  gener->SetTrackingFlag(0);
  gener->SelectMotherParticles(selectedMothers);
  gener->SetPtParamPi0(paramPi0);
  gener->SetPtParamEta(paramEta);
  gener->SetPtParamOmega(paramOmega);
  gener->SetPtParamPhi(paramPhi);
  gener->SetCollisionSystem(collisionsSystem); 			//pp 7 TeV
  gener->SetCentrality(centrality);				// kpp

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
