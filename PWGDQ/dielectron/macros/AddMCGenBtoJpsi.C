AliGenerator* AddMCGenBtoJpsi(Bool_t useEvtGenForB=kFALSE, Double_t energyConfig=13000.){
  gSystem->Load("liblhapdf"); 
  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");  
  // load libraries to use Evtgen
  gSystem->Load("libPhotos");
  gSystem->Load("libEvtGen");
  gSystem->Load("libEvtGenExternal");
  gSystem->Load("libTEvtGen");  

  //Generating a cocktail
  AliGenCocktail *gener = new AliGenCocktail();
  gener->UsePerEventRates();
   // Jpsi from B
  AliGenPythia *pythia = new AliGenPythia(-1);
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-2., 2.);
  pythia->SetPtRange(0, 1000.);
  pythia->SetProcess(kPyBeautyppMNRwmi);
  pythia->SetEnergyCMS(energyConfig);
  pythia->SetTune(320); // Perugia0 tune
  pythia->UseNewMultipleInteractionsScenario();
  if(useEvtGenForB) pythia->SetForceDecay(kNoDecayBeauty);
  else {
  pythia->SetCutOnChild(1);
  pythia->SetPdgCodeParticleforAcceptanceCut(443);
  pythia->SetChildYRange(-2, 2);
  pythia->SetChildPtRange(0, 10000.);
  pythia->SetForceDecay(kBJpsiUndecayed);
  }
  pythia->SetStackFillOpt(AliGenPythia::kHeavyFlavor);
  //
  // 
  AliGenEvtGen *gene = new AliGenEvtGen();
  gene->SetForceDecay(kBJpsiDiElectron);
  if(useEvtGenForB) gene->SetParticleSwitchedOff(AliGenEvtGen::kHFPart);
  else gene->SetParticleSwitchedOff(AliGenEvtGen::kCharmPart);

  gener->AddGenerator(pythia,         "Pythia",         1.);
  gener->AddGenerator(gene, "EvtGen", 1.);
  //
  return gener;


return;
}
