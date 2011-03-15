
void ConfigGammaJet (TString inputfile = "files.txt" ) {

  Int_t debugLevel = 0;
  
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN.so");
  gSystem->Load("libPWG4GammaConv.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");


  TChain * chain = createAODChain (inputfile);
  if(!chain) {
    cout << "Errror in chain cration"<<endl;
    return -1;
  }
  
  AliAODInputHandler * inpHandler = new AliAODInputHandler();
  inpHandler->AddFriend("AliAODGammaConversion.root");
  
  AliAnalysisManager *mgr  = new AliAnalysisManager("GammaJet Manager", "GammaJet Manager");
  mgr->SetInputEventHandler  (inpHandler);
  mgr->SetDebugLevel(debugLevel);

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(), AliAnalysisManager::kOutputContainer, "histos.root");

  // AliAnaConvIsolation * iso3 = new AliAnaConvIsolation(0.7, 0.5, 0.5, 0.5, 0.5);
  // gammaJetAna->AddIsolationAna(dynamic_cast<TObject*>(iso3));


  AliAnalysisTaskGammaJet * gammaJetAna = new AliAnalysisTaskGammaJet("gamma jet analysis");
  gammaJetAna->SetDebugLevel(0);
  gammaJetAna->SetConversionCutId("90022670901120321036000000090");
  gammaJetAna->SetMinNTracks(0);
  mgr->AddTask(gammaJetAna);

  AliAnaConvIsolation * isolation = gammaJetAna->GetIsolation();
  isolation->SetMaxPtThreshold(0.7);
  isolation->SetSumPtThreshold(2.5);
  isolation->SetSumPtFraction(0.1);
  isolation->SetConeSize(0.4);
 
  AliAnaConvCorrPhoton * ghAna = new AliAnaConvCorrPhoton("photon_wDec");
  ghAna->DoDecayParticles();
  gammaJetAna->AddPhotonHadronAna(dynamic_cast<TObject*>(ghAna));

  AliAnaConvCorrPhoton * ghAna2 = new AliAnaConvCorrPhoton("photon_noDec");
  ghAna2->SkipDecayParticles();
  gammaJetAna->AddPhotonHadronAna(dynamic_cast<TObject*>(ghAna2));

  AliAnaConvCorrPion * gPiAna = new AliAnaConvCorrPion("pionHadron");
  gammaJetAna->AddPionHadronAna(dynamic_cast<TObject*>(gPiAna));

  inpHandler->AddFriend("AliAODs.pwg4jets.root");
  AliAnaConvCorrPhotonJet* gJetAna = new AliAnaConvCorrPhotonJet("photonJet");
  gJetAna->SetTriggerPt(0.0);
  gJetAna->SetCorrelatedPt(0.0);
  gammaJetAna->AddPhotonJetAna(gJetAna);

  // gROOT->LoadMacro("AddTaskGammaJet.C");
  // AliAnalysisTaskGammaJet * gj2 =  AddTaskGammaJet("test", 3.0, 40, 0.4);
  // isolation = gj2->GetIsolation();
  // isolation->SetMinPt(0.3);


  mgr->ConnectInput  (gammaJetAna,  0, cinput1  );
  mgr->ConnectOutput (gammaJetAna,  1, coutput2 );
  
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
  //mgr->StartAnalysis("local",chain, 10000);

}


///_________________________________________________________________________________________
TChain * createAODChain (TString inputfile = "allfiles.txt") {

  chain = new TChain("aodTree");

  TString line;
  ifstream in;
  in.open(inputfile.Data());
  while (in.good()) {
    in >> line;
    if (line.Length() == 0) continue;
    chain->Add(line.Data());
  }

  return chain;
  
}
