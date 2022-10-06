//based on 888
TString names("anaFilter_1;anaFilter_2");
TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
Int_t GetN(){return nDie;}
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition);
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition);

//________________________________________________________________
AliAnalysisFilter* Config_hmurakam_ElectronEfficiencyV2(Int_t cutDefinition)
{

  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;

  AliAnalysisFilter *anaFilter = new AliAnalysisFilter(Form("anaFilter_%d",cutDefinition),Form("anaFilter_%d",cutDefinition)); // named constructor seems mandatory!
  // do not change these initial values!
  //Int_t selectedPID=-1;
  //Bool_t isPrefilterCutset=kFALSE;

  AliDielectronV0Cuts *noconv = new AliDielectronV0Cuts("IsGamma","IsGamma");
  // which V0 finder you want to use
  noconv->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
  // add some pdg codes (they are used then by the KF package and important for gamma conversions)
  noconv->SetPdgCodes(22,11,11); // mother, daughter1 and 2
  // add default PID cuts (defined in AliDielectronPID)
  // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
  //noconv->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
  // add the pair cuts for V0 candidates
  noconv->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kR,                             3.0,  90.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  // selection or rejection of V0 tracks
  noconv->SetExcludeTracks(kTRUE);
  
  anaFilter->AddCuts( SetupTrackCuts(cutDefinition) );
  anaFilter->AddCuts( SetupPIDcuts(cutDefinition) );
  if(cutDefinition==1){
	cout << "Conversion !!!!!!!!!!!!!!!" <<endl;
	anaFilter->AddCuts( noconv );
  }
  std::cout << "...cuts added!" <<std::endl;

  anaFilter->Print();

  return anaFilter;
}
//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;

  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  // pT and eta
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2,   1e30);
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8,   0.8);

  //TPC
  trackCutsDiel->SetRequireTPCRefit(kTRUE);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);//(1)
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,  0.8,   1.5);//(2)
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,       0.0,   4.0);//(3)
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);//(4)
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracTPC,    0.0,   0.4);//(5)

  //ITS
  trackCutsDiel->SetRequireITSRefit(kTRUE);
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,         3.0, 100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,       0.0,   4.5);

  //primary selection
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,   3.0);

  printf("Add shared cluster cut\n");
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS, 1.0, 6.0, kTRUE);//accept no shared cls hit (default)

  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("Trackcuts","Trackcuts",AliDielectronCutGroup::kCompAND);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(trackCutsDiel);

  trackCuts->Print();

  return trackCuts;

}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition)
{

  std::cout << "SetupPIDcuts()" <<std::endl;

  //PID1
  AliDielectronPID *pidTPCTOFreq = new AliDielectronPID("pidTPCTOFreq","pidTPCTOFreq");//<--> same as recoverTOF
  pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,     -3. , 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. , 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCTOFreq->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,      -3., 3., 0.4, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kP);// TOF required
  //PID2
  AliDielectronPID *pidTPCHadRej = new AliDielectronPID("pidTPCHadRej","pidTPCHadRej");//pure TPC pid
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,  -3., 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kPion,    -100., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kKaon,      -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kProton,    -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC

  AliAnalysisCuts* fancyCut=0x0;
  AliDielectronCutGroup* combinedPIDcuts = new AliDielectronCutGroup("combinedPIDcuts","combinedPIDcuts",AliDielectronCutGroup::kCompOR);
  combinedPIDcuts->AddCut(pidTPCTOFreq);//PID1
  combinedPIDcuts->AddCut(pidTPCHadRej);//PID2
  fancyCut = combinedPIDcuts;

  return fancyCut;

}

//_____________________________________________________________________________________________
const AliDielectronEventCuts *GetEventCuts(){
  
  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Any && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex(kTRUE);
  //eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD);
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;

}
