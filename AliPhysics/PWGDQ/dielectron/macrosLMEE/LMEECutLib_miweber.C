class LMEECutLib {
  
public:

  static  enum LMMECutSet {
    kpp2010All_pidITSTPCTOFif_trkSPDfirst_1,      // Cut set for testing (pp 2010, all multiplicties, PID like QM2014 preliminaries)
    kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1, // Cut set for testing (pp 2010, all multiplicties, PID like QM2014 preliminaries, tight DCA)
    kpp2010All,                              
    kpp2010Low,                              
    kpp2010Mid,                              
    kpp2010High,                              
    kCUTSETMAX
  };
  
  
  LMEECutLib() {}
  
  AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
  AliAnalysisCuts*            GetCentralityCuts(Int_t centSel);
  AliDielectronTrackRotator*  GetTrackRotator(Int_t cutSet);
  AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);
  
  AliAnalysisCuts* GetPairCutsAna(Int_t cutSet, Int_t togglePC=0); //Bool_t togglePC=kFALSE
  AliAnalysisCuts* GetPairCutsPre(Int_t cutSet);
  
  AliAnalysisCuts* GetPIDCutsAna(Int_t cutSet);
  AliAnalysisCuts* GetPIDCutsPre(Int_t cutSet);
  
  AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet);
  AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet);

  void SetEtaCorrection(AliDielectron *die, Int_t selPID, Int_t selCent, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/

  
};


void LMEECutLib::SetEtaCorrection(AliDielectron *die, Int_t selPID, Int_t selCent, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  // MW: IS THIS STILL NEEDED???
  //
  printf("starting LMEECutLib::SetEtaCorrection()\n");
  printf(" corrZdim = %f\n", corrZdim);
  printf(" corrYdim = %f\n", corrYdim);
  //Bool_t hasMC=die->GetHasMC();
  //Bool_t hasTuneOnData=kFALSE; //((AliAnalysisTaskPIDResponse*)AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0))->GetTuneOnData();
  //printf("tune on data switched: %d \n",hasTuneOnData);
  // printf("name task at 0: %s \n",AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0)->GetName());
  // printf("input event %p \n", AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  // printf("pid response %p \n",((AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->GetPIDResponse());
  // printf("pid response task %p \n",AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0));
  // AliAnalysisManager::GetAnalysisManager()->GetTasks()->At(0)->Dump();;
  
  // AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  // AliPIDResponse* pidResponse = inputHandler->GetPIDResponse();
  // if(pidResponse) hasTuneOnData = pidResponse->IsTunedOnData();
  // printf("man %p inp %p pid %p ====> %d \n",man,inputHandler,pidResponse,hasTuneOnData);
  
  TF2 *fCntrdCorr=0x0;
  TF2 *fWdthCorr=0x0;
  Double_t fitMinZdim, fitMaxZdim;
  Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrZdim==AliDielectronVarManager::kRefMultTPConly) 
  {
    if (selCent == kpp2010All) { fitMinZdim=400.; fitMaxZdim=1400; }
    
    fCntrdCorr = new TF2("fCntrdCorr", "[9]*([0]+[1]*x+[2]*x*x)+[10]*([3]+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y+[8]*y*y*y*y*y)+[11]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    fWdthCorr  = new TF2("fWdthCorr",  "[9]*([0]+[1]*x+[2]*x*x)+[10]*([3]+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y+[8]*y*y*y*y*y)+[11]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    Double_t parCntrd[]={1.984011e-01, -2.538876e-04, 7.636316e-09, 2.731183e-01, 2.190149e-01, -3.437725e+00, -6.438784e-01, 5.317945e+00, 6.089976e-01, 1.018070e+00, 1.002105e+00, 2.428792e-02};
    Double_t parWdth[] ={1.116969e+00, 5.063439e-05, -5.526647e-09, 1.207224e+00, -3.633499e-02, -5.340433e-01, 1.557123e-01, 7.865346e-01, -1.619051e-01, 1.019816e+00, 1.004669e+00, -1.185922e+00};
    fCntrdCorr->SetParameters(parCntrd);
    fWdthCorr ->SetParameters(parWdth);
  }
  else {
    printf(" no eta correction applied!\n");
    return;
  }
  
  die->SetCentroidCorrFunction(fCntrdCorr, corrZdim, corrYdim);
  die->SetWidthCorrFunction(fWdthCorr, corrZdim, corrYdim);
  
  printf(" TPC PID eta correction loaded!!!\n");
}

// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask, currently to 'kpp2010'
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
  case kpp2010All:
    
    AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
    
    //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
    eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
    eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
    eventCuts->SetRequireVertex();
    eventCuts->SetMinVtxContributors(1);
    eventCuts->SetVertexZ(-10.,10.);
    break;
    
  default: cout << "No Event Cut defined" << endl;
  }
  return eventCuts;
}


//Selection of relatively 'flat' centralities
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(Int_t centSel) {
  AliDielectronVarCuts* centCuts = 0x0;
  switch (centSel) {
  case kpp2010All:
    centCuts = new AliDielectronVarCuts("centCuts","Multiplicitypp2010All");
    centCuts->AddCut(AliDielectronVarManager::kNacc,0.,500.);
    break;
  case kpp2010Low:
    centCuts = new AliDielectronVarCuts("centCuts","Multiplicitypp2010Low");
    centCuts->AddCut(AliDielectronVarManager::kNacc,0.,19.);
    break;
  case kpp2010Mid:
    centCuts = new AliDielectronVarCuts("centCuts","Multiplicitypp2010Mid");
    centCuts->AddCut(AliDielectronVarManager::kNacc,20.,39.);
    break;
  case kpp2010High:
    centCuts = new AliDielectronVarCuts("centCuts","Multiplicitypp2010High");
    centCuts->AddCut(AliDielectronVarManager::kNacc,40.,500.);
    break;
  default: cout << "No Centrality selected" << endl;
  }
  return centCuts;
}


//Basic track rotator settings from J/Psi, more investigation needed
AliDielectronTrackRotator* LMEECutLib::GetTrackRotator(Int_t cutSet) {
  AliDielectronTrackRotator* trackRotator = 0x0;
  switch (cutSet) {
    default: cout << "No Rotator defined" << endl;
      //default:
      //  trackRotator = new AliDielectronTrackRotator();
      //  trackRotator->SetIterations(20);
      //  trackRotator->SetConeAnglePhi(TMath::Pi()/180*165);
      //  trackRotator->SetStartAnglePhi(TMath::Pi());
      //  break;
  }
  return trackRotator;
}


AliDielectronMixingHandler* LMEECutLib::GetMixingHandler(Int_t cutSet) {
  AliDielectronMixingHandler* mixingHandler = 0x0;
  switch (cutSet) {
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:
    mixingHandler = new AliDielectronMixingHandler;
    mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
    mixingHandler->AddVariable(AliDielectronVarManager::kNacc,"0,500");
    // for using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
    // mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
    mixingHandler->SetDepth(10);
    mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
    break;
    //[...]
  default: cout << "No Mixer defined" << endl;
  }
  return mixingHandler;
}



//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsAna(Int_t cutSet, Int_t togglePC)  {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  switch (cutSet) {
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:
    cout << "No Pair Cuts used - ok " << endl; 
    break;
  default: cout << "No Pair Cuts defined " << endl;
  }
  return pairCuts;
}


//Pair Cuts for PREFILTER step
// cuts = REJECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsPre(Int_t cutSet)  {  
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  switch (cutSet) {
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:
    AliDielectronVarCuts *pairCutsPhiV = new AliDielectronVarCuts("pairCutsPhiV","pairCutsPhiV");//mass and Phiv together
    pairCutsPhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.05);
    pairCutsPhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.5 , 3.2 );
    pairCuts = pairCutsPhiV;
    break;
    
  default: cout << "No Prefilter Pair Cuts defined " << endl;
  } 
  return pairCuts;
}



AliAnalysisCuts* LMEECutLib::GetPIDCutsAna(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;
  
  //-----------------------------------------------
  // Define different PID Cuts, that are used later
  //-----------------------------------------------
  // PID cuts depend on TPC_inner_p, if not specified
  // PID cut ranges correspond to global momentum P
  // check it again!!!
  //-----------------------------------------------
  
  //
  //
  //TPC: electron inclusion asymmetric
  //     pion     exclusion 3sigma
  //TOF: electron inclusion 3sigma in region where p,K cross electrons in TPC
  AliDielectronPID *pidTPCTOF_Semi1 = new AliDielectronPID("pidTPCTOF_Semi1","pidTPCTOF_Semi1");
  pidTPCTOF_Semi1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
  pidTPCTOF_Semi1->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
  pidTPCTOF_Semi1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,1.7 , kFALSE);
  //
  //
  // LOOSE PID TPC+TOF
  AliDielectronPID *pidTPCTOF_Semi_LOOSE = new AliDielectronPID("pidTPCTOF_Semi_LOOSE","pidTPCTOF_Semi_LOOSE");
  pidTPCTOF_Semi_LOOSE->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
  pidTPCTOF_Semi_LOOSE->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,1.7 , kFALSE);
  //
  //
  // PID TPC only
  AliDielectronPID *pidTPC_3 = new AliDielectronPID("pidTPC_3","pidTPC_3");
  pidTPC_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
  pidTPC_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
  
  
  //TPC: electron inclusion asymmetric
  //     pion     exclusion 3sigma
  //ITS: electron inclusion asymmetric in region where p,K cross electrons in TPC
  //TOF: electron inclusion 3sigma in similar region - BUT ONLY IF AVAILABLE
  AliDielectronPID *pidTPCITS_TOFif1 = new AliDielectronPID("pidTPCITS_TOFif1","pidTPCITS_TOFif1");
  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,1.5 , kFALSE);
  pidTPCITS_TOFif1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,1.7 , kFALSE, AliDielectronPID::kIfAvailable);
  //
  //
  //TPC: electron inclusion asymmetric
  //     pion     exclusion 3sigma
  //ITS: electron inclusion asymmetric OVER FULL MOMENTUM RANGE
  //TOF: electron inclusion 3sigma - BUT ONLY IF AVAILABLE
  AliDielectronPID *pidTPCITS_TOFif2 = new AliDielectronPID("pidTPCITS_TOFif2","pidTPCITS_TOFif2");
  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
  pidTPCITS_TOFif2->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  //
  //
  // LOOSE PID ITS+TPC+TOFif
  AliDielectronPID *pidTPCITS_TOFif_LOOSE = new AliDielectronPID("pidTPCITS_TOFif_LOOSE","pidTPCITS_TOFif_LOOSE");
  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-10. ,20. , 0. ,100., kFALSE);
  pidTPCITS_TOFif_LOOSE->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  //
  //
  // PID ITS+TPC
  AliDielectronPID *pidTPCITS_3 = new AliDielectronPID("pidTPCITS_3","pidTPCITS_3");
  pidTPCITS_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
  pidTPCITS_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
  pidTPCITS_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
  //
  //
  // tighter PID ITS+TPC+TOFif
  // ITS only up to momentum where proton contamination is seen in TPC signal
  AliDielectronPID *pidTPCITS_TOFif56 = new AliDielectronPID("pidTPCITS_TOFif56","pidTPCITS_TOFif56");
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 2.5, 0. ,100., kFALSE);
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 0.5, 0. ,  2., kFALSE);
  pidTPCITS_TOFif56->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  
  
  // PID for V0 task
  AliDielectronPID *pid_V0select_1 = new AliDielectronPID("pid_V0select_1","pid_V0select_1");
  pid_V0select_1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
  pid_V0select_1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -1.5, 1.5, 0. ,100., kFALSE);
  
  
  // eta range:
  AliDielectronVarCuts *etaRange090 = new AliDielectronVarCuts("etaRange090","etaRange090");
  etaRange090->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
  AliDielectronVarCuts *etaRange084 = new AliDielectronVarCuts("etaRange084","etaRange084");
  etaRange084->AddCut(AliDielectronVarManager::kEta, -0.84, 0.84);
  AliDielectronVarCuts *etaRange080 = new AliDielectronVarCuts("etaRange080","etaRange080"); //changed from 0.76 to 0.8 on 2014-08-12!
  etaRange080->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
  // pt range:
  AliDielectronVarCuts *ptRange400to3500 = new AliDielectronVarCuts("ptRange400to3500","ptRange400to3500");
  ptRange400to3500->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
  
  
  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------
  switch (cutSet) {
      
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:
    AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
    cgPIDCutsAna->AddCut(etaRange080);
    cgPIDCutsAna->AddCut(ptRange400to3500);
    cgPIDCutsAna->AddCut(pidTPCITS_TOFif2);
    cgPIDCutsAna->AddCut(GetTrackCutsAna(cutSet));
    pidCuts = cgPIDCutsAna;
    break;
    default: cout << "No Analysis PID Cut defined " << endl;
  }
  return pidCuts;
}


//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

AliAnalysisCuts* LMEECutLib::GetTrackCutsAna(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (cutSet) {
      
    //----------
    // these MAIN settings just load the main track selection directly below:
    //----------
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
    AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
    trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
    AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
    trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2010(kFALSE); loose DCA, 2D cut
    trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    
    cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
    cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
    cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
    trackCuts = cgTrackCutsAnaSPDfirst;
    break;

  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:
    AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
    trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
    AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
    trackCutsDiel->SetAODFilterBit(1<<5); // (=32) filterbit 5! //GetStandardITSTPCTrackCuts2010(kFALSE); 
    trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    
    cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
    cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
    cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
    trackCuts = cgTrackCutsAnaSPDfirst;
    break;
    
  default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
} 



//Relaxed PID cuts for additional rejectin step, do not use blindly
AliAnalysisCuts* LMEECutLib::GetPIDCutsPre(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;
  switch (cutSet) {
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:

    // eta range:
    AliDielectronVarCuts *etaRangePre1 = new AliDielectronVarCuts("etaRangePre1","etaRangePre1");
    etaRangePre1->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
    // pt range:
    AliDielectronVarCuts *ptRangePre1 = new AliDielectronVarCuts("ptRangePre1","ptRangePre1");
    ptRangePre1->AddCut(AliDielectronVarManager::kPt, .2, 3.5); // 0.2 is realistic. turnon at ~180MeV
    
    AliDielectronCutGroup* cgITSTPCTOFpre = new AliDielectronCutGroup("cgITSTPCTOFpre","cgITSTPCTOFpre",AliDielectronCutGroup::kCompAND);
    AliDielectronPID *pidITSTPCTOFpre = new AliDielectronPID("pidITSTPCTOFpre","pidITSTPCTOFpre");
    pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3., 0. ,100., kFALSE);
    pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
    // ITS will be used:
    pidITSTPCTOFpre->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 2., 0. , 1.8, kFALSE);
    // TOF will be used if available, and with pt instead of p:
    pidITSTPCTOFpre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0.4,100., kFALSE, 
			    AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
    cgITSTPCTOFpre->AddCut(pidITSTPCTOFpre);
    cgITSTPCTOFpre->AddCut(etaRangePre1);
    cgITSTPCTOFpre->AddCut(ptRangePre1);
    cgITSTPCTOFpre->AddCut(GetTrackCutsAna(cutSet));
    
    AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
    cgInitialTrackFilter->AddCut(GetPIDCutsAna(cutSet)); // in case the prefilter cuts do not include all needed global tracks.
    cgInitialTrackFilter->AddCut(cgITSTPCTOFpre);
    
    pidCuts = cgInitialTrackFilter;   
    break;
    
  default: cout << "No Prefilter PID Cut defined " << endl;
  }
  return pidCuts;
}


//Possibly different cut sets for Prefilter step
//Not used at the moment
AliAnalysisCuts* LMEECutLib::GetTrackCutsPre(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (cutSet) {
  case kpp2010All_pidITSTPCTOFif_trkSPDfirst_1:
  case kpp2010All_pidITSTPCTOFif_trkSPDfirsttight_1:
    AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
    AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
    trackCutsDiel->SetAODFilterBit(1); 
    
    cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
    cgTrackCutsPre->AddCut(trackCutsDiel);
    cgTrackCutsPre->AddCut(trackCutsAOD);
    trackCuts = cgTrackCutsPre;
    break;
    
  default: cout << "No Prefilter Track Cut defined " << endl;
  }
  return trackCuts;
}



