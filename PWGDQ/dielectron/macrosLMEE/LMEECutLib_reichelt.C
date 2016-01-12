class LMEECutLib {
  
public:
  
  static enum enCentSel {
    kPbPb2011Central=0,
    kPbPb2011_00to10,
    kPbPb2011MidCentral,
    kPbPb2011_10to20,
    kPbPb2011SemiCentral,
    kPbPb2011_20to50,
    kPbPb2011_10to50,
    kPbPb2011_00to50,
    kCENTSELMAX
  };
  
  static enum enPairCut {
    kPairCut_OFF=0,
    kPairCut_mee10_theta30,
    kPairCut_mee20_theta20, // not yet done (maybe for prefilter with SPD+SDD or as default PairCut for final pairs)
    kPairCut_mee20_theta50,
    kPairCut_mee30_theta60,
    kPairCut_mee40_theta80,
    kPairCut_mee60_theta100,
    kPairCut_mee200_theta300, // for testing
    kPairCut_phiv157_mee40,
    kPairCut_phiv157_mee60,
    kPairCut_phiv157_mee80,
    kPairCut_phiv157_mee100,
    kPairCut_phiv236_mee40,
    kPairCut_phiv236_mee60,
    kPairCut_phiv236_mee80,
    kPairCut_phiv236_mee100
  };
  
  static enum enKineCut {
    kKineCut_pt50_eta080=0,
    kKineCut_pt50_eta090,
    kKineCut_pt200_eta080,
    kKineCut_pt200_eta090,
    kKineCut_pt300_eta080,
    kKineCut_pt300_eta090,
    kKineCut_pt400_eta080,
    kKineCut_pt400_eta090
  };
  
  static enum LMMECutSet {
    //kCutSetMIN=0,
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl,          // syst 8
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl,          // syst 7
    kPbPb2011_V0select_2_looseNoTOF,
    kPbPb2011_V0select_1,
    kPbPb2011_V0select_1_Arm,
    kPbPb2011V0_2_loose,  // (NO FULL CUTSET)
    kPbPb2011V0_1,        // (NO FULL CUTSET)
    kPbPb2011MC_pi0Dal_1,
    kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight,  // syst 6 (no train run yet)
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight,   // syst 5 (no train run yet)
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight,       // syst 4
    kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight,       // syst 2
    kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4,            // syst 3
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4,             // syst 1
    kPbPb2011_pidITSTPC_trkSPDfirst_3,            // (cutSet w/o pairing)
    kPbPb2011_pidTPC_trkSPDfirst_3,               // (cutSet w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose, // (cutSet w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose, // (cutSet w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose,      // (cutSet w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose,      // (cutSet w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1,
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1,       // cutSet for Technical Preliminaries for QM2014 (no prefilter used!)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_1,
    kPbPb2011_pidTPCTOF_trkSPDfirst_1,
    kCutSetMAX,
    //
    // PID
    kPbPb2011PID_ITSTPCTOFif_1=kCutSetMAX, // ITS+TPC+TOFf
    kPbPb2011PID_ITSTPCTOFif_2,
    kPbPb2011PID_ITSTPCTOFif_3,
    kPbPb2011PID_ITSTPCTOFif_LOOSE,
    kPbPb2011PID_ITSTPC_1,      // ITS+TPC
    kPbPb2011PID_TPCTOF_1,      // TPC+TOF
    kPbPb2011PID_TPCTOF_LOOSE,
    kPbPb2011PID_TPC_1,         // TPC
    kPbPb2011PID_TPC_2,
    kPbPb2011PID_TPC_pre,
    kPbPb2011PID_ITS_1,         // ITS
    kPIDMAX,
    //
    // Quality
    kPbPb2011TRK_SPDfirst_1=kPIDMAX,    // main track selection. with SPD first
    kPbPb2011TRK_SPDfirst_2,            // main track selection. >5 ITS cls
    kPbPb2011TRK_SDDfirstSPDnone_1,     // complimentary tracks. strictly without SPD, to be combined with others!
    kPbPb2011TRK_SDDfirstSPDnone_2,     // complimentary tracks. >4 ITS cls
    kPbPb2011TRK_Minimal_1,
    kQualityMAX
  };
  
  //char* LMEECutNames[kQualityMAX] = { "PbPb2011TPCandTOF","PbPb2011TPCorTOF"};
  
  
  LMEECutLib();
  
  AliDielectronEventCuts*     GetEventCuts(Int_t cutSet, Bool_t hasMC=kFALSE);
  AliAnalysisCuts*            GetCentralityCuts();
  AliDielectronTrackRotator*  GetTrackRotator();
  AliDielectronMixingHandler* GetMixingHandler();
  
  AliAnalysisCuts* GetPairCutsAna();
  AliAnalysisCuts* GetPairCutsPre(Int_t cutSet=-1);
  // main track cut functions:
  AliAnalysisCuts* GetTrackCutsAna();
  AliAnalysisCuts* GetTrackCutsPre();
  AliAnalysisCuts* GetESDTrackCutsAna();
  
  void      SetIsQATask(Bool_t b=kTRUE)         { fIsQATask=b; }
  void      SetIsRandomRejTask(Bool_t b=kTRUE)  { fIsRandomRejTask=b; }
  void      SetDoRejectionStep(Bool_t b=kTRUE)  { fDoRejectionStep=b; }
  void      SetFillPureMC(Bool_t b=kTRUE)       { fFillPureMC=b; }

  Bool_t    GetDoRejectionStep() { return fDoRejectionStep; }
  
  void      SetITSSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void      SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim);
  void      SetITSSigmaEleCorrectionMC(AliAnalysisTaskElectronEfficiency *task, Int_t corrZdim, Int_t corrYdim);
  void      SetTPCSigmaEleCorrectionMC(AliAnalysisTaskElectronEfficiency *task, Int_t corrZdim, Int_t corrYdim);
  
  void      AddMCSignals(AliDielectron *die, Int_t cutDefinition);
  void      InitHistograms(AliDielectron *die, Int_t cutDefinition);
  void      InitCF(AliDielectron* die, Int_t cutDefinition);
  
  // kept public to avoid writing setters...
  Int_t     selectedCentrality;
  Int_t     selectedPIDAna;
  Int_t     selectedPIDPre;
  Int_t     selectedQualityAna;
  Int_t     selectedQualityPre;
  Int_t     selectedKineCutsAna;
  Int_t     selectedKineCutsPre;
  Int_t     selectedPairCutsAna;
  Int_t     selectedPairCutsPre;
  
  
private:
  // internal track cut functions (called by GetTrackCuts):
  AliAnalysisCuts* GetKineCutsAna();
  AliAnalysisCuts* GetKineCutsPre(Int_t cutSet=-1);
  AliAnalysisCuts* GetPIDCuts(Int_t cutSet=-1);
  AliAnalysisCuts* GetQualityCuts(Int_t cutSet=-1, Int_t doExclusion=0);
  AliAnalysisCuts* GetMCTrackCuts();
  // helper functions
  TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
  TVectorD* GetVector(Int_t var);
  
  Bool_t    fIsQATask;
  Bool_t    fIsRandomRejTask;
  Bool_t    fDoRejectionStep;
  Bool_t    fFillPureMC;
  
  static enum enCutType {
    kInclude = 0,
    kExclude = 1,
    kCUTTYPEMAX
  };
  
  static enum enHistVars {
    kMee=0, kMee500,
    kPtee, kP2D, kRuns,
    kPhiV, kOpAng, kOpAng2,
    kEta2D, kEta3D, kPhi2D, kY3D,
    kSigmaEle, kSigmaOther, kTPCdEdx,
    kPairDCAsigXY, kPairDCAabsXY, kPairLinDCAsigXY, kPairLinDCAabsXY
  };
  
};


//_______________________________________________________________________________________________
LMEECutLib::LMEECutLib() :
selectedCentrality(-1),
selectedPIDAna(-1),
selectedPIDPre(-1),
selectedQualityAna(-1),
selectedQualityPre(-1),
selectedKineCutsAna(LMEECutLib::kKineCut_pt200_eta080),
selectedKineCutsPre(LMEECutLib::kKineCut_pt50_eta090),
selectedPairCutsAna(LMEECutLib::kPairCut_OFF),
selectedPairCutsPre(LMEECutLib::kPairCut_OFF),
fIsQATask(kFALSE),
fIsRandomRejTask(kFALSE),
fDoRejectionStep(kFALSE),
fFillPureMC(kFALSE)
{
  // Constructor
}


//_______________________________________________________________________________________________
void LMEECutLib::SetITSSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the ITS, can be one/two/three-dimensional
  //
  printf("starting LMEECutLib::SetITSSigmaEleCorrection()\n --> for MC task, use ::SetITSSigmaEleCorrectionMC()!\n");
  printf(" correction Zdim = %s\n", AliDielectronVarManager::GetValueName(corrZdim));
  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
  
  TF2 *fCntrdCorr=0x0;
  TF2 *fWdthCorr=0x0;
  Double_t fitMinZdim=   0;
  Double_t fitMaxZdim=4000;
  Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  if (corrZdim==AliDielectronVarManager::kRefMultTPConly)
  {
    fCntrdCorr = new TF2("fCntrdCorr", "[0]+0.*[1]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    fWdthCorr  = new TF2("fWdthCorr",  "[0]+0.*[1]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    Double_t parCntrd[]={0.,0.}; //-0.5 for tests
    Double_t parWdth[] ={1.,0.};
    fCntrdCorr->SetParameters(parCntrd);
    fWdthCorr ->SetParameters(parWdth);
  }
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
//  die->SetCentroidCorrFunctionITS(fCntrdCorr, corrZdim, corrYdim);
//  die->SetWidthCorrFunctionITS(fWdthCorr, corrZdim, corrYdim);
//  printf(" ITS PID eta correction loaded!\n");
}

//_______________________________________________________________________________________________
void LMEECutLib::SetTPCSigmaEleCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  printf("starting LMEECutLib::SetTPCSigmaEleCorrection()\n --> for MC task, use ::SetTPCSigmaEleCorrectionMC()!\n");
  printf(" correction Zdim = %s\n", AliDielectronVarManager::GetValueName(corrZdim));
  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));
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
  Double_t fitMinZdim=   0;
  Double_t fitMaxZdim=4000;
  Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  if (corrZdim==AliDielectronVarManager::kRefMultTPConly)
  {
    fCntrdCorr = new TF2("fCntrdCorr", "[9]*([0]+[1]*x+[2]*x*x)+[10]*([3]+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y+[8]*y*y*y*y*y)+[11]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    fWdthCorr  = new TF2("fWdthCorr",  "[9]*([0]+[1]*x+[2]*x*x)+[10]*([3]+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y+[8]*y*y*y*y*y)+[11]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    Double_t parCntrd[]={2.003967e-01, -2.584087e-04, 1.123938e-08, 1.559636e-01, 2.542528e-01, -3.392378e+00, -7.109151e-01, 5.374647e+00, 6.646497e-01, 1.015413e+00, 1.007024e+00, 1.289175e-01}; // fit of 00-50% (fitrange 400-2300, |eta| < 0.8)
    //  Double_t parCntrd[]={1.984011e-01, -2.538876e-04, 7.636316e-09, 2.731183e-01, 2.190149e-01, -3.437725e+00, -6.438784e-01, 5.317945e+00, 6.089976e-01, 1.018070e+00, 1.002105e+00, 2.428792e-02}; // fit of 20-50% (fitrange 400-1400, |eta| < 0.8)
    Double_t parWdth[] ={1.122035e+00, 3.897153e-05, -2.394857e-10, 1.219613e+00, -4.266204e-02, -5.047624e-01, 1.818052e-01, 7.607361e-01, -2.020611e-01, 1.005953e+00, 1.024095e+00, -1.209148e+00}; // fit of 00-50% (fitrange 400-2300, |eta| < 0.8)
    //  Double_t parWdth[] ={1.116969e+00, 5.063439e-05, -5.526647e-09, 1.207224e+00, -3.633499e-02, -5.340433e-01, 1.557123e-01, 7.865346e-01, -1.619051e-01, 1.019816e+00, 1.004669e+00, -1.185922e+00}; // fit of 20-50% (fitrange 400-1400, |eta| < 0.8)
    fCntrdCorr->SetParameters(parCntrd);
    fWdthCorr ->SetParameters(parWdth);
  }
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
  die->SetCentroidCorrFunction(fCntrdCorr, corrZdim, corrYdim);
  die->SetWidthCorrFunction(fWdthCorr, corrZdim, corrYdim);
  printf(" TPC PID eta correction loaded!\n");
}


//_______________________________________________________________________________________________
void LMEECutLib::SetITSSigmaEleCorrectionMC(AliAnalysisTaskElectronEfficiency *task, Int_t corrZdim, Int_t corrYdim) {
  //
  // MC post-correction for the centroid and width of electron sigmas in the ITS, can be one/two/three-dimensional
  //
  TF2 *fCntrdCorr=0x0;
  TF2 *fWdthCorr=0x0;
  Double_t fitMinZdim=   0;
  Double_t fitMaxZdim=4000;
  Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  if (corrZdim==AliDielectronVarManager::kNacc)
  {
    fCntrdCorr = new TF2("fCntrdCorr", "[0]+0.*[1]*x*y", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    fWdthCorr  = new TF2("fWdthCorr",  "[0]+0.*[1]*x*y", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    Double_t parCntrd[]={0.,0.}; //-0.5 for tests
    Double_t parWdth[] ={1.,0.};
    fCntrdCorr->SetParameters(parCntrd);
    fWdthCorr ->SetParameters(parWdth);
  }
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
//  task->SetCentroidCorrFunctionITS(fCntrdCorr, corrZdim, corrYdim);
//  task->SetWidthCorrFunctionITS(fWdthCorr, corrZdim, corrYdim);
//  printf(" MC ITS PID eta correction loaded!\n");
}
  
//_______________________________________________________________________________________________
void LMEECutLib::SetTPCSigmaEleCorrectionMC(AliAnalysisTaskElectronEfficiency *task, Int_t corrZdim, Int_t corrYdim) {
  //
  // MC post-correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  TF2 *fCntrdCorr=0x0;
  TF2 *fWdthCorr=0x0;
  Double_t fitMinZdim=   0;
  Double_t fitMaxZdim=4000;
  Double_t fitMinEta=-0.9, fitMaxEta=+0.9; // Ydim, usually eta
  
  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }
  
  if (corrZdim==AliDielectronVarManager::kNacc)
  {
    fCntrdCorr = new TF2("fCntrdCorr", "[0]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    fWdthCorr  = new TF2("fWdthCorr",  "[0]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    Double_t parCntrd[]={0.}; //0.2 for tests
    Double_t parWdth[] ={1.0};
    fCntrdCorr->SetParameters(parCntrd);
    fWdthCorr ->SetParameters(parWdth);
  }
  else {
    printf(" no correction available for Zdim = %s!\n", AliDielectronVarManager::GetValueName(corrZdim));
    printf(" no correction applied!\n");
    return;
  }
  
//  task->SetCentroidCorrFunction(fCntrdCorr, corrZdim, corrYdim);
//  task->SetWidthCorrFunction(fWdthCorr, corrZdim, corrYdim);
//  printf(" MC TPC PID eta correction loaded!\n");
}



// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask, currently to 'kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1'
//_______________________________________________________________________________________________
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet, Bool_t hasMC) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
      eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
      eventCuts->SetRequireVertex();
      eventCuts->SetMinVtxContributors(1);
      eventCuts->SetVertexZ(-10.,10.);
      if (hasMC) {
        //eventCuts->SetVertexType(.....);
        eventCuts->SetCentralityRange(0.0,50.0);
      } else {
        eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
        //eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTPC); // kVtxAny // AOD
      }
      break;
    default: cout << "No Event Cut defined" << endl;
  }
  return eventCuts;
}


//Selection of relatively 'flat' centralities is a bit difficult...
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetCentralityCuts() {
  AliAnalysisCuts* centCuts = 0x0;
  
  // Online Trigger MB+Semi
  AliDielectronVarCuts *trgMBorSEMICNT = new AliDielectronVarCuts("trgMBorSEMICNT","trgMBorSEMICNT");
  trgMBorSEMICNT->SetCutType(AliDielectronVarCuts::kAny);
  trgMBorSEMICNT->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  1, kFALSE); // MB
  trgMBorSEMICNT->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  7, kFALSE); // SEMICNT
  //trgMBorSEMICNT->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  4, kFALSE); //CNT
  
  switch (selectedCentrality) {
    case kPbPb2011Central:
    case kPbPb2011_00to10:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_00to10");
      centRange->AddCut(AliDielectronVarManager::kCentrality,0.,10.);
      centCuts = centRange;
      break;
    case kPbPb2011MidCentral:
    case kPbPb2011_10to20:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_10to20");
      centRange->AddCut(AliDielectronVarManager::kCentrality,10.,20.);
      AliDielectronCutGroup* centCutsCG =new AliDielectronCutGroup("centCutsCG","centCutsCG",AliDielectronCutGroup::kCompAND);
      centCutsCG->AddCut(centRange);
      centCutsCG->AddCut(trgMBorSEMICNT);
      centCuts = centCutsCG;
      break;
    case kPbPb2011SemiCentral:
    case kPbPb2011_20to50:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_20to50");
      centRange->AddCut(AliDielectronVarManager::kCentrality,20.,50.);
      centCuts = centRange;
      break;
    case kPbPb2011_00to50:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_00to50");
      centRange->AddCut(AliDielectronVarManager::kCentrality,0.,50.);
      ////centRange->AddBitCut(AliDielectronVarManager::kTriggerExclOFF,  4, kTRUE); //this is wrong...
      ////centRange->AddBitCut(AliDielectronVarManager::kTriggerInclONL,  4, kTRUE); //also wrong...
      AliDielectronCutGroup* centCutsCG =new AliDielectronCutGroup("centCutsCG","centCutsCG",AliDielectronCutGroup::kCompAND);
      centCutsCG->AddCut(centRange);
      centCutsCG->AddCut(trgMBorSEMICNT);
      centCuts = centCutsCG;
      break;
    case kPbPb2011_10to50:
      AliDielectronVarCuts* centRange = new AliDielectronVarCuts("centCuts","CentralityPbPb2011_10to50");
      centRange->AddCut(AliDielectronVarManager::kCentrality,10.,50.);
      AliDielectronCutGroup* centCutsCG =new AliDielectronCutGroup("centCutsCG","centCutsCG",AliDielectronCutGroup::kCompAND);
      centCutsCG->AddCut(centRange);
      centCutsCG->AddCut(trgMBorSEMICNT);
      centCuts = centCutsCG;
      break;
    default: cout << "No Centrality selected" << endl;
  }
  return centCuts;
}


//Basic track rotator settings from J/Psi, more investigation needed
//_______________________________________________________________________________________________
AliDielectronTrackRotator* LMEECutLib::GetTrackRotator() {
  AliDielectronTrackRotator* trackRotator = 0x0;
  Int_t cutSet=-1;
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


//_______________________________________________________________________________________________
AliDielectronMixingHandler* LMEECutLib::GetMixingHandler() {
  AliDielectronMixingHandler* mixingHandler = 0x0;
  Int_t cutSet=1;
  switch (cutSet) {
    case 1:
      mixingHandler = new AliDielectronMixingHandler();
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,30,50,80");
      // now using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
      mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(15);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
      //[...]
    default: cout << "No Mixer defined" << endl;
  }
  return mixingHandler;
}



//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetPairCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronVarCuts* pairVarCuts = (AliDielectronVarCuts*) GetPairCutsPre(selectedPairCutsAna);
  if (!pairVarCuts) return 0x0;
  pairVarCuts->InvertCuts();
  return pairVarCuts;
}

//Pair Cuts for PREFILTER step
// cuts = REJECTION!!!
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetPairCutsPre(Int_t cutSet) {  
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  // for the default function call, pick the cutSet according to 'selectedPairCutsPre', which is set via the Config file.
  if (cutSet<0) { cutSet = selectedPairCutsPre; }
  
  if (cutSet==kPairCut_OFF) {
    cout << "Pair Cuts disabled " << endl;
    return 0x0;
  }
  AliDielectronVarCuts* pairVarCuts = new AliDielectronVarCuts("pairVarCuts","pairVarCuts");
  switch (cutSet) {
      // Mee and opening angle
    case kPairCut_mee10_theta30:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.01);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.03);
      break;
    case kPairCut_mee20_theta20:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.02);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.02);
      break;
    case kPairCut_mee20_theta50:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.02);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05);
      break;
    case kPairCut_mee30_theta60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.03);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.06);
      break;
    case kPairCut_mee40_theta80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.08);
      break;
    case kPairCut_mee60_theta100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.10);
      break;
    case kPairCut_mee200_theta300: // just for testing
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.20);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.30);
      break;
      
      // Mee and phiv
    case kPairCut_phiv157_mee40:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv157_mee60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv157_mee80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.08);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv157_mee100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.10);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee40:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.08);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
    case kPairCut_phiv236_mee100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.10);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      break;
      
    default: cout << "No (Prefilter) Pair Cuts defined " << endl;
  } 
  return pairVarCuts;
}



//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetTrackCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
  cgTrackCutsAna->AddCut( GetPIDCuts(selectedPIDAna) );
  cgTrackCutsAna->AddCut( GetQualityCuts(selectedQualityAna) );
  cgTrackCutsAna->AddCut( GetKineCutsAna() );
  
  return cgTrackCutsAna;
}

//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetTrackCutsPre() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
  cgTrackCutsPre->AddCut( GetPIDCuts(selectedPIDPre) );
  cgTrackCutsPre->AddCut( GetQualityCuts(selectedQualityPre) );
  cgTrackCutsPre->AddCut( GetKineCutsPre() );
  
  // in case the prefilter cuts do not include all needed global tracks, we create an "OR" cutgroup:
  AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
  cgInitialTrackFilter->AddCut( cgTrackCutsPre );
  cgInitialTrackFilter->AddCut( GetTrackCutsAna() );
  
  return cgInitialTrackFilter;
}



//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetKineCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetKineCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronVarCuts* kineCuts = (AliDielectronVarCuts*) GetKineCutsPre(selectedKineCutsAna);
  if (!kineCuts) return 0x0;
  return kineCuts;
}

//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetKineCutsPre(Int_t cutSet) {  
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetKineCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  // for the default function call, pick the cutSet according to 'selectedKineCutsPre', which is set via the Config file.
  if (cutSet<0) { cutSet = selectedKineCutsPre; }
  
  AliDielectronVarCuts* kineCuts = new AliDielectronVarCuts("kineCuts","kineCuts");
  switch (cutSet) {
    case kKineCut_pt50_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .05, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt50_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .05, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
    case kKineCut_pt200_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .2, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt200_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .2, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
    case kKineCut_pt300_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .3, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt300_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .3, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
    case kKineCut_pt400_eta080:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      break;
    case kKineCut_pt400_eta090:
      kineCuts->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
      kineCuts->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
      break;
      
    default: cout << "No (Prefilter) Kine Cuts defined " << endl;
  } 
  return kineCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetPIDCuts(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //-----------------------------------------------
  // PID cuts depend on TPC_inner_p, if not specified
  // PID cut ranges correspond to global momentum P
  // check it again!!!
  //-----------------------------------------------
  
  if (cutSet<0) {
    cout << "Invalid cutSet selected! " << endl;
    return 0x0;
  }
  
  AliDielectronPID *pidCuts = new AliDielectronPID("pidCuts","pidCuts");
  switch (cutSet) {
      ///
      /// ITS + TPC + TOFif
      ///
      // ===== Default PID ITS+TPC+TOFif =====
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl:
    case kPbPb2011MC_pi0Dal_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011PID_ITSTPCTOFif_1:
      //TPC: electron inclusion asymmetric
      //     pion     exclusion 3sigma
      //ITS: electron inclusion asymmetric OVER FULL MOMENTUM RANGE
      //TOF: electron inclusion 3sigma - BUT ONLY IF AVAILABLE
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      // ===== Tighter PID ITS+TPC+TOFif =====
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011PID_ITSTPCTOFif_2:
      //ITS: electron inclusion asymmetric in region where p,K cross electrons in TPC
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 2.5, 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 0.5, 0. ,  2., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2. , 2. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      // ===== Relaxed PID ITS+TPC+TOFif =====
    case kPbPb2011PID_ITSTPCTOFif_3:
      //ITS: electron inclusion asymmetric in region where p,K cross electrons in TPC
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 2., 0. ,  2., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      // ===== LOOSE PID ITS+TPC+TOFif =====
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose: // loose "ITSTPCTOFif" PID - for 2D contamination study
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
    case kPbPb2011PID_ITSTPCTOFif_LOOSE:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-10. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      
      ///
      /// ITS+TPC
      ///
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
    case kPbPb2011PID_ITSTPC_1:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3. , 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 1. , 0. ,100., kFALSE);
      return pidCuts;
      break;
      
      
      ///
      /// TPC + TOF
      ///
      // ===== Default PID TPC+TOF =====
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
    case kPbPb2011PID_TPCTOF_1:
      //TPC: electron inclusion asymmetric
      //     pion     exclusion 3sigma
      //TOF: electron inclusion 3sigma in region where p,K cross electrons in TPC
      //     (may be useful to use TOF with pt instead of p, because pt determines its acceptance.)
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,1.7 , kFALSE);
      return pidCuts;
      break;
      // ===== LOOSE PID TPC+TOF =====
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose: // loose "ITSTPCTOFif" PID - for 1D contamination study in TPC
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
    case kPbPb2011PID_TPCTOF_LOOSE:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,1.7 , kFALSE);
      return pidCuts;
      break;
      
      
      ///
      /// TPC (+TOFif)
      ///
    case kPbPb2011_pidTPC_trkSPDfirst_3:
    case kPbPb2011PID_TPC_1:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      return pidCuts;
      break;
    case kPbPb2011PID_TPC_pre:
    case kPbPb2011PID_TPC_2:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      return pidCuts;
      break;
      
      ///
      /// ITS
      ///
    case kPbPb2011PID_ITS_1:
      pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 3.);
      return pidCuts;
      break;
      
      ///
      /// Other
      ///
    case kPbPb2011_V0select_1_Arm:
    case kPbPb2011_V0select_1:
      // PID for V0 task
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -1.5, 1.5, 0. ,100., kFALSE);
      return pidCuts;
      break;
    case kPbPb2011_V0select_2_looseNoTOF:
      pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
      // no TOF to really see what gets rejected...
      return pidCuts;
      break;
      
    default: cout << "No PID Cut defined " << endl;
      return 0x0;
  }
  return 0x0;
  
}


//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetQualityCuts(Int_t cutSet, Int_t doExclusion) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetQualityCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  
  if (cutSet<0) {
    // for the default function call, pick the cutSet according to 'selectedQualityAna', which is set via the Config file.
    //cutSet = selectedQualityAna;
    cout << "Invalid cutSet selected! " << endl;
    return 0x0;
  }
  
  switch (cutSet) {
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl:
      // combine 
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1));
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011V0_2_loose, kExclude));
      trackCuts = cgTrackCutsAna;
      break;
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl:
      // combine 
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1));
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011V0_2_loose, kExclude));
      //// instead, only FOR DOING A COMPARISON: ////
      ////      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1));
      trackCuts = cgTrackCutsAna;
      break;
      
    case kPbPb2011_V0select_2_looseNoTOF: doExclusion = kInclude; // and go on with kPbPb2011V0_2_loose... (no break)
    case kPbPb2011V0_2_loose:
      // primarily used for exclusion of V0s from track samples.
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      // which V0 finder you want to use
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly // kOnTheFly better for small radii (quote Julian)
      // add some pdg codes (they are used then by the KF package and important for gamma conversions)
      gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      // add default PID cuts (defined in AliDielectronPID)
      // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
      //gammaV0Cuts->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
      // add the pair cuts for V0 candidates
      // variations from Julian in ( ... )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.05),   1.0,  kFALSE); // ( 0.02 -- 0.05 )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.20, kFALSE); // ( 0.05 -- 0.2 )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE); // ( 0.05 -- 0.1 )
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
      //      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // ( |0.35|   ------- 0.3 )
      // selection or rejection of V0 tracks
      if (doExclusion==kExclude)  gammaV0Cuts->SetExcludeTracks(kTRUE);
      else                        gammaV0Cuts->SetExcludeTracks(kFALSE);
      // add the V0cuts directly to the track filter or to some cut group of it
      
      // some more default cuts automatically set by AliDielectronV0Cuts::InitEvent() !!!
      
      cgTrackCutsV0excl = new AliDielectronCutGroup("cgTrackCutsV0excl","cgTrackCutsV0excl",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0excl->AddCut(gammaV0Cuts);
      trackCuts = cgTrackCutsV0excl;
      break;
      
    case kPbPb2011_V0select_1_Arm: doExclusion = kInclude;
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
      gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // should increase purity...
      if (doExclusion==kExclude)  gammaV0Cuts->SetExcludeTracks(kTRUE);
      else                        gammaV0Cuts->SetExcludeTracks(kFALSE);
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      cgTrackCutsV0select->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsV0select;
      break;
      
    case kPbPb2011_V0select_1: doExclusion = kInclude; // and go on with kPbPb2011V0_1... (no break)
    case kPbPb2011V0_1:
      // primarily used for selection of V0s for TPC eta correction. with additional track cuts.
      // inspired by "AliDielectronV0Cuts.cxx"
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      // which V0 finder you want to use
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
      // add some pdg codes (they are used then by the KF package and important for gamma conversions)
      gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      // add default PID cuts (defined in AliDielectronPID)
      // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
      //gammaV0Cuts->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
      // add the pair cuts for V0 candidates
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
      //gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // should increase purity...
      // selection or rejection of V0 tracks
      if (doExclusion==kExclude)  gammaV0Cuts->SetExcludeTracks(kTRUE);
      else                        gammaV0Cuts->SetExcludeTracks(kFALSE);
      // add the V0cuts directly to the track filter or to some cut group of it
      
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      // some more default cuts automatically set by AliDielectronV0Cuts::InitEvent() !!!
      
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      cgTrackCutsV0select->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsV0select;
      break;
      
      //----------
      // MC
      //----------
    case kPbPb2011MC_pi0Dal_1:
      AliDielectronVarCuts* trackCutsMC =new AliDielectronVarCuts("trackCutsMC","trackCutsMC");
      trackCutsMC->SetCutOnMCtruth(kTRUE);
      trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -11.01,  11.01);
      //trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -10.01,  10.01, kTRUE); //excludeRange
      trackCutsMC->AddCut(AliDielectronVarManager::kPdgCodeMother, 111);
      
      cgTrackCutsMC = new AliDielectronCutGroup("cgTrackCutsMC","cgTrackCutsMC",AliDielectronCutGroup::kCompAND);
      cgTrackCutsMC->AddCut(GetQualityCuts(kPbPb2011TRK_SPDfirst_1));
      cgTrackCutsMC->AddCut(trackCutsMC);
      trackCuts = cgTrackCutsMC;
      break;
      
      //----------
      // these MAIN settings have to combine different track selections:
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
      // combine typical and new trackcuts with "kCompOR" condition:
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAnaSPDorSDD","cgTrackCutsAnaSPDorSDD",AliDielectronCutGroup::kCompOR);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SPDfirst_1));         // typical trackcuts with requirement of SPD
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SDDfirstSPDnone_1)); // new additional trackcuts with SDD instead of SPD
      trackCuts = cgTrackCutsAna;
      break;
      
      //----------
      // these MAIN settings just load the main track selection directly below:
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
    case kPbPb2011_pidTPC_trkSPDfirst_3:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      //----------
    case kPbPb2011TRK_SPDfirst_1: // main track selection, now closer to what Hongyan does...
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
      
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      
      
    case kPbPb2011TRK_SDDfirstSPDnone_1: // complimentary tracks, strictly without SPD, to be combined with others!
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0); // means at least 3 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
      
      cgTrackCutsAnaSDDnoSPD = new AliDielectronCutGroup("cgTrackCutsAnaSDDnoSPD","cgTrackCutsAnaSDDnoSPD",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsDiel);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSDDnoSPD;
      break;
      
      //----------
      // MAIN settings - combined trackset - variation 1
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
      // combine typical and new trackcuts with "kCompOR" condition:
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAnaSPDorSDD","cgTrackCutsAnaSPDorSDD",AliDielectronCutGroup::kCompOR);
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SPDfirst_2));         // typical trackcuts with requirement of SPD
      cgTrackCutsAna->AddCut(GetQualityCuts(kPbPb2011TRK_SDDfirstSPDnone_2)); // new additional trackcuts with SDD instead of SPD
      trackCuts = cgTrackCutsAna;
      break;
      
      //----------
      // MAIN settings - single trackset - variation 1
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      //----------
    case kPbPb2011TRK_SPDfirst_2: // main track selection, 5+ ITS clusters
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     5.0, 100.0); // means at least 3 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
      
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      
    case kPbPb2011TRK_SDDfirstSPDnone_2: // complimentary tracks, 4+ ITS clusters, strictly without SPD, to be combined with others!
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 4 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<6); // GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
      
      cgTrackCutsAnaSDDnoSPD = new AliDielectronCutGroup("cgTrackCutsAnaSDDnoSPD","cgTrackCutsAnaSDDnoSPD",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsDiel);
      cgTrackCutsAnaSDDnoSPD->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSDDnoSPD;
      break;
      
      //----------
      // for Prefilter
      //----------
    case kPbPb2011TRK_Minimal_1:
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1); //does nothing for ESDs, ITSSA(???) // maybe use FilterBit(2) instead!
      //        trackCutsDiel->SetRequireITSRefit(kTRUE); //function in AliDielectronTrackCuts
      //        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //function in AliDielectronTrackCuts
      
      cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
      cgTrackCutsPre->AddCut(trackCutsDiel);
      cgTrackCutsPre->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsPre;
      break;
      //----------
      
      //[...]
    default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
} 


//*******************************************************************************
//*******************************************************************************
//** ESD TRACK CUTS TUNED FOR AGREEMENT BETWEEN AODS AND ESDS  ******************
//** NOT NECESSARILY 100% OPTIMIZED FOR DIEL-ANALYSIS          ******************
//*******************************************************************************
//*******************************************************************************

//WHEN RUNNING ON ESDs: LOAD Default Cuts for AODs
//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetESDTrackCutsAna() {
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetESDTrackCutsAna()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> Setting ESD Track Cuts >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> ( do we run on ESD?! ) >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliESDtrackCuts* esdTrackCutsH = 0x0;
  switch (selectedPIDAna) {
    default:
      // standard cuts with very loose DCA: Bit4 (Int: 16), AOD095&115
      
      esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
      esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
      esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
      esdTrackCutsH->SetDCAToVertex2D(kTRUE);
      
      //The cuts below should be the onyl ones that are missing
      //explicitely in the TrackCutsAna method
      //To be sure, StandardITSTPCTrackCuts is loaded however
      /* 
       esdTrackCutsH = new AliESDtrackCuts();
       esdTrackCutsH->SetAcceptKinkDaughters(kFALSE);
       //Not done so far via dielectron cuts:
       */
      /*
       esdTrackCuts->SetDCAToVertex2D(kFALSE);
       esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
       esdTrackCuts->SetMaxChi2PerClusterITS(36);
       */
      
      break;
      //default: cout << "No ESD Track Cut defined " << endl;
  }
  return esdTrackCutsH;
}



//_______________________________________________________________________________________________
AliAnalysisCuts* LMEECutLib::GetMCTrackCuts() {
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetMCTrackCuts()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* trackCuts=0x0;
  switch (selectedPIDAna) {
    default:
      AliDielectronVarCuts* trackCutsMC =new AliDielectronVarCuts("trackCutsMC","trackCutsMC");
      trackCutsMC->SetCutOnMCtruth(kTRUE);
      trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -11.01,  11.01);
      //trackCutsMC->AddCut(AliDielectronVarManager::kPdgCode      , -10.01,  10.01, kTRUE); //excludeRange
      //trackCutsMC->AddCut(AliDielectronVarManager::kPdgCodeMother, 111);
      trackCuts = trackCutsMC;
      break;
  }
  return trackCuts;
}



//_______________________________________________________________________________________________
void LMEECutLib::AddMCSignals(TNamed* task, Int_t cutDefinition){
  //Do we have an MC handler?
  //if (!die->GetHasMC()) return;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  AddMCSignals()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  
  AliDielectronSignalMC* eleFinalState = new AliDielectronSignalMC("eleFinalState","eleFinalState");
  eleFinalState->SetFillPureMCStep(fFillPureMC);
  eleFinalState->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFinalState->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); // for leg electrons: kFinalState = kPrimary
  //mother
  eleFinalState->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary); // has no effect for this leg setting (kFinalState). Still includes injected J/psi.
  
  AliDielectronSignalMC* eleDirect = new AliDielectronSignalMC("eleDirect","eleDirect");
  eleDirect->SetFillPureMCStep(fFillPureMC);
  eleDirect->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleDirect->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleDirect->SetLegSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect); // for leg electrons: kDirect is empty
  
  AliDielectronSignalMC* eleNoCocktail = new AliDielectronSignalMC("eleNoCocktail","eleNoCocktail");
  eleNoCocktail->SetFillPureMCStep(fFillPureMC);
  eleNoCocktail->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleNoCocktail->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleNoCocktail->SetLegSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail); // for leg electrons: kNoCocktail = kFinalState + kSecondary
  
  AliDielectronSignalMC* eleSecondary = new AliDielectronSignalMC("eleSecondary","eleSecondary");
  eleSecondary->SetFillPureMCStep(fFillPureMC);
  eleSecondary->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleSecondary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondary->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC* eleSecondaryWeak = new AliDielectronSignalMC("eleSecondaryWeak","eleSecondaryWeak");
  eleSecondaryWeak->SetFillPureMCStep(fFillPureMC);
  eleSecondaryWeak->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleSecondaryWeak->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondaryWeak->SetLegSources(AliDielectronSignalMC::kSecondaryFromWeakDecay, AliDielectronSignalMC::kSecondaryFromWeakDecay); // for leg electrons: kSecondaryFromWeakDecay is empty
  
  
  AliDielectronSignalMC* eleFromBGEvent = new AliDielectronSignalMC("eleFromBGEvent","eleFromBGEvent");
  eleFromBGEvent->SetFillPureMCStep(fFillPureMC);
  eleFromBGEvent->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFromBGEvent->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFromBGEvent->SetLegSources(AliDielectronSignalMC::kFromBGEvent, AliDielectronSignalMC::kFromBGEvent);

  AliDielectronSignalMC* eleFinalStateFromBGEvent = new AliDielectronSignalMC("eleFinalStateFromBGEvent","eleFinalStateFromBGEvent");
  eleFinalStateFromBGEvent->SetFillPureMCStep(fFillPureMC);
  eleFinalStateFromBGEvent->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFinalStateFromBGEvent->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromBGEvent->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  
  AliDielectronSignalMC* eleFinalStateFromBGEvent_mDir = new AliDielectronSignalMC("eleFinalStateFromBGEvent_mDir","eleFinalStateFromBGEvent_mDir");
  eleFinalStateFromBGEvent_mDir->SetFillPureMCStep(fFillPureMC);
  eleFinalStateFromBGEvent_mDir->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleFinalStateFromBGEvent_mDir->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromBGEvent_mDir->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  //mother
  eleFinalStateFromBGEvent_mDir->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
  
  
  // In Pythia the mother label of primary hadrons gives the (positive) label of the mother quark,
  //  so mother source kNoCocktail should correspond to non-injected signals.
  // In Hijing however the mother label of primary hadrons is negative (they are not traced back to the quarks),
  //  so mother source kNoCocktail gives hadrons which have a mother, i.e. electrons whith have a grandmother.
  AliDielectronSignalMC* eleWithGrandmother = new AliDielectronSignalMC("eleWithGrandmother","eleWithGrandmother");
  eleWithGrandmother->SetFillPureMCStep(fFillPureMC);
  eleWithGrandmother->SetLegPDGs(11,1);  //dummy second leg (never MCtrue)
  eleWithGrandmother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleWithGrandmother->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState); // for leg electrons: kFinalState = kPrimary, also in case of 'SetMotherSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail)'.
  //mother
  eleWithGrandmother->SetMotherSources(AliDielectronSignalMC::kNoCocktail, AliDielectronSignalMC::kNoCocktail); // excludes JPsi, includes conversions.
  eleWithGrandmother->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons
  
  
  // stuff:
  //mother
  //  ele->SetCheckBothChargesMothers(kTRUE,kTRUE);
  //  ele->SetMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
  //  ele->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  //  ele->SetGrandMotherPDGs(902,902,kTRUE,kTRUE); // exclude open charm,beauty hadrons
  //  ele->SetGrandMotherPDGs(500,500,kTRUE,kTRUE); // exclude non-prompt jpsi eletrons
  
  
  AliDielectronSignalMC* gammaConversion = new AliDielectronSignalMC("gammaConversion","gammaConversion");
  gammaConversion->SetFillPureMCStep(fFillPureMC);
  gammaConversion->SetLegPDGs(11,-11);
  gammaConversion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  gammaConversion->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  //mother
  gammaConversion->SetMotherPDGs(22,22);
  gammaConversion->SetMothersRelation(AliDielectronSignalMC::kSame);
  
  // add direct di lepton resonances
  /*
   AliDielectronSignalMC* directP[7];
   TParticlePDG *ap;
   Int_t pdg[] = {111, 113, 221, 223, 331, 333, 443};
   for(Int_t i=0; i<7; i++) {
   ap = TDatabasePDG::Instance()->GetParticle(pdg[i]);
   directP[i] = new AliDielectronSignalMC(Form("direct%s",ap->GetName()),Form("direct%s",ap->GetName()));
   directP[i]->SetLegPDGs(11,-11);
   directP[i]->SetMotherPDGs(pdg[i],pdg[i]);
   directP[i]->SetMothersRelation(AliDielectronSignalMC::kSame);
   directP[i]->SetFillPureMCStep(kTRUE);
   directP[i]->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
   directP[i]->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
   // directP[i]->SetCheckBothChargesLegs(kTRUE,kTRUE);
   // directP[i]->SetCheckBothChargesMothers(kTRUE,kTRUE);
   }
   */
  
  
  if (task->IsA()==AliAnalysisTaskElectronEfficiency::Class()) {
    // selection. only the first signal will be used.
//    (static_cast<AliAnalysisTaskElectronEfficiency*>task)->AddSignalMC(eleFinalStateFromBGEvent);
    (static_cast<AliAnalysisTaskElectronEfficiency*>task)->AddSignalMC(eleFinalStateFromBGEvent_mDir);
  }
  else if (task->IsA()==AliDielectron::Class()) {
    if (! (static_cast<AliDielectron*>task)->GetHasMC()) return;
    // selection. multiple signals are possible.
    (static_cast<AliDielectron*>task)->AddSignalMC(eleFromBGEvent);
    (static_cast<AliDielectron*>task)->AddSignalMC(eleFinalStateFromBGEvent);
    (static_cast<AliDielectron*>task)->AddSignalMC(eleFinalStateFromBGEvent_mDir);
    (static_cast<AliDielectron*>task)->AddSignalMC(eleWithGrandmother);
    //die->AddSignalMC(elePrimary);
    //die->AddSignalMC(eleFinalState);    // for leg electrons: kFinalState = kPrimary
    //die->AddSignalMC(eleDirect);        // for leg electrons: kDirect is empty
    //die->AddSignalMC(eleNoCocktail);    // for leg electrons: kNoCocktail = kFinalState + kSecondary
    //die->AddSignalMC(eleSecondary);     // for leg electrons: kSecondary = number of 'gammaConversion'
    //die->AddSignalMC(eleSecondaryWeak); // for leg electrons: kSecondaryFromWeakDecay is empty
    //die->AddSignalMC(gammaConversion); //
    //for(Int_t i=0; i<7; i++) die->AddSignalMC(directP[i]);
    
    for (Int_t i=0; i<(static_cast<AliDielectron*>task)->GetMCSignals()->GetEntriesFast(); ++i) {
      cout << "added MCsignal: " << (static_cast<AliDielectron*>task)->GetMCSignals()->At(i)->GetName() << endl;
    }
  }
}


//_______________________________________________________________________________________________
void LMEECutLib::InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack;Random");
  
  //Event class
  histos->AddClass("Event"); // all classes will be stored in 'THashList fHistoList'
  
  //Track classes
  //to fill also track info from 2nd event loop until 3
  // in AliDielectron.cxx: fgkTrackClassNames[4] = {"ev1+","ev1-","ev2+","ev2-"};
  if (die->DoEventProcess()) {
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  // fgkPairClassNames[11] = {
  //  "ev1+_ev1+",  "ev1+_ev1-",  "ev1-_ev1-",  // 0-2 (same event)
  //  "ev1+_ev2+",  "ev1-_ev2+",  "ev2+_ev2+",  // 3-4 (+5)
  //  "ev1+_ev2-",  "ev1-_ev2-",                // 6-7
  //  "ev2+_ev2-",  "ev2-_ev2-",  "ev1+_ev1-_TR"
  // };
  
  if (!fIsQATask && !fIsRandomRejTask) // -> analysis with pairing
  {
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      // Legs of final Pairs. Both charges together. No duplicate entries.
      if (!fIsQATask) {
        histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
      }
    }
    
    //Mixed event and track rotation
    if (die->GetMixingHandler()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }
    if (die->GetTrackRotator()) {
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
    }
    
    if (fDoRejectionStep) 
    {
      //PreFilter Classes
      //to fill also track info from 2nd event loop until 2
//      for (Int_t i=0; i<2; ++i){
//        histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
//        // class 'Pre_%s': "Fill Histogram information for tracks after prefilter"(AliDielectron::FillHistogramsTracks(...)), so it is identical to Track_%s if no paircut is used.
//      }
      
      //Create Classes for Rejected Tracks/Pairs:
      for (Int_t i=0; i<3; ++i){
        histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
        // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
        histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
      }
      
      /*
       //track rotation
       histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
       histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
       */
    }// end: (fDoRejectionStep)
    
  }// end: (!fIsQATask)
  
  
  if (fIsRandomRejTask)
  {
    //
    // _____ histograms for AliAnalysisTaskRandomRejection _____
    //
    histos->AddClass("Rand_Pair");
    histos->AddClass("Rand_RejPair");
    const char* cRandomPairClassNames[4] = { "Testpart", "DataEle", "RejTestpart", "RejDataEle" };
    for (Int_t i=0; i<4; ++i){
      histos->AddClass(Form("Random_%s",cRandomPairClassNames[i]));
    }
    histos->UserHistogram("Random","Pt","",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Random","Eta","",200,-2,2,AliDielectronVarManager::kEta);
    histos->UserHistogram("Random","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","Px","",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Random","Py","",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Random","Pz","",200,0,10.,AliDielectronVarManager::kPz);
    histos->UserHistogram("Random","Pt_Eta_Phi","",GetVector(kP2D),GetVector(kEta2D),GetVector(kPhi2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    histos->UserHistogram("Random","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Random","TPC_dEdx_Pt","",
                          100,0.,2.,100,0.,100., AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCsignal);
  }
  
	//add histograms to event class
	histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","RunNumber","",AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600"),AliDielectronVarManager::kRunNumber);
	histos->UserHistogram("Event","Centrality","","-1,0,10,20,30,40,50,60,70,80,90,100,101;#events",AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","centrality","",100,0,100,AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","nESDTracks","",8000,0,80000,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","Nacc","",8000,0,8000,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","RefMultTPConly","",8000,0,8000,AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","epTPC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kTPCrpH2uc);
  histos->UserHistogram("Event","epV0AC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kv0ACrpH2);
  histos->UserHistogram("Event","epV0AC_epTPC","",120,-TMath::PiOver2(),TMath::PiOver2(),120,-TMath::PiOver2(),TMath::PiOver2(),AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);
  
  
  //add histograms to Track classes
  // axis labels are set to the values in 'AliDielectronVarManager.cxx' if the histogram title is empty or starts with ';'! [see AliDielectronHistos::AdaptNameTitle(...)]
  histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",
                        GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  
  // ITS
  histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                        GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","ITSnSigmaEle_Eta",";Eta;n#sigma_{ele}^{ITS}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  // 3D
  if (die->DoEventProcess()) {
    histos->UserHistogram("Track","ITSnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{ITS};p (GeV/c)",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kP);
    histos->UserHistogram("Track","ITSnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{ITS};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
    //histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    //histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
  }
  // TPC
  histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                        GetVector(kP2D), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  if (die->DoEventProcess()) {
    // some histograms with SigmaEleRaw to benchmark the dEdx-eta-correction:
    histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta",";Eta;n#sigma_{ele,Raw}^{TPC}",
                          GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw);
    histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPCnSigmaEleRaw_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele,Raw}^{TPC}",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEleRaw);
    histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPCnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TPC}",
                          GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
    // 3D
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  }
  
  if (fIsQATask) {
    histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
                          GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPC_dEdx_Eta_P",";Eta;TPC signal (arb units);p_{in} (GeV/c)",
                          GetVector(kEta3D), GetVector(kTPCdEdx), GetVector(kP2D),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
    histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TPC signal (arb units)",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kTPCdEdx),
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal);
    
    //histos->UserHistogram("Track","TPCnSigmaEle_Eta_P"," ... // see above
    //histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly", ... // see above
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_RunNumber",";Eta;n#sigma_{ele}^{TPC};run",
                          GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kRuns),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_RunNumber",";N_{TPC ref};n#sigma_{ele}^{TPC};run",
                          BinsToVector(100,0.,5000.), GetVector(kSigmaEle), GetVector(kRuns),
                          AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_P_RunNumber",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};run",
                          GetVector(kP2D), GetVector(kSigmaEle), GetVector(kRuns), 
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Track","RefMultTPConly_Nacc",";N_{TPC ref};N_{acc}",
                          BinsToVector(100,0.,5000.), BinsToVector(100,0.,5000.), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Track","TPCnSigmaEleRaw_Eta_RefMultTPConly",";Eta;n#sigma_{ele,Raw}^{TPC};N_{TPC ref}",
                          GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw,AliDielectronVarManager::kRefMultTPConly);
    
    histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    histos->UserHistogram("Track","TPCnSigmaKao_Eta",";Eta;n#sigma_{kaon}^{TPC}",
                          GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  }
  
  // TRD
  // TRD variables need lot of computing time. since a performance update by Julian, they will not be computed if not needed! (status 7.3.14)
  //  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;p_{in} (GeV/c);TRD prob Electrons",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle);
  //  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;p_{in} (GeV/c);TRD prob Pions",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio);
  
  // TOF
  histos->UserHistogram("Track","TOFbeta_P",";p_{in} (GeV/c);TOF beta",
                        GetVector(kP2D), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  if (fIsQATask) {
    histos->UserHistogram("Track","TOFnSigmaPio_P",";p_{in} (GeV/c);TOF number of sigmas Pions",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    histos->UserHistogram("Track","TOFnSigmaKao_P",";p_{in} (GeV/c);TOF number of sigmas Kaons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
    histos->UserHistogram("Track","TOFnSigmaPro_P",";p_{in} (GeV/c);TOF number of sigmas Protons",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  }    
  // 2D-PID
  if (fIsQATask) {
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
                          50,0.,2.5, 160,-12.,20., 150,-10.,20.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
                          50,0.,2.5, 160,-12.,20., 50,-5.,5.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);
  }
  
  // Eta, Phi, Pt
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Phi","",GetVector(kEta2D),GetVector(kPhi2D), AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Pt","",GetVector(kP2D),GetVector(kEta2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi_Pt","",GetVector(kP2D),GetVector(kPhi2D), AliDielectronVarManager::kPt,AliDielectronVarManager::kPhi);
  
  // DCA
  histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  
  // Quality
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  
  if (fIsQATask) {
    histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
                          160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRows_Pt",";Pt [GeV];TPC crossed rows",
                          GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt",";Pt [GeV];TPC crossed rows over findable",
                          GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
                          100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
                          120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
    //V0
    // these are pair variables, but only internal in the V0-finder...
    //histos->UserHistogram("Track","ArmenterosAlpha_ArmenterosPt",";kArmPt;kArmAlpha",
    //                      100,0.,0.5, 120,-0.6,0.6, AliDielectronVarManager::kArmPt,AliDielectronVarManager::kArmAlpha);
    //histos->UserHistogram("Track","TPCnSigmaEle_ArmenterosPt",";kArmPt;TPCnSigmaEle",
    //                      100,0.,0.5, 160,-12.,20., AliDielectronVarManager::kArmPt,AliDielectronVarManager::kTPCnSigmaEle);
  }
  
  if (!fIsQATask) 
  {
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass","",500,0.,5.,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","PairPt","",160,0.,8., AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Rapidity","",200,-2.,2.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","",240,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
    
    //2D and 3D histograms
    histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                          GetVector(kMee), GetVector(kPtee),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
                          100,-1.,1., 120,0.,TMath::TwoPi(),
                          AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                          GetVector(kMee), GetVector(kPtee), GetVector(kOpAng), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                          GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PairPt_Rapidity",";Inv. Mass [GeV];Pair Pt [GeV];Rapidity",
                          GetVector(kMee), GetVector(kPtee), GetVector(kY3D), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kY);
    
    if (die->DoEventProcess()) {
      //opening angle and PhiV
      histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                            GetVector(kMee), GetVector(kOpAng), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
      histos->UserHistogram("Pair","InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                            GetVector(kMee), GetVector(kPhiV), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
      histos->UserHistogram("Pair","PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                            GetVector(kPtee), GetVector(kOpAng), 
                            AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
      histos->UserHistogram("Pair","PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                            GetVector(kPtee), GetVector(kPhiV), 
                            AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
      histos->UserHistogram("Pair","OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                            GetVector(kOpAng), GetVector(kPhiV), 
                            AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
      
      // pair DCA
      histos->UserHistogram("Pair","InvMass_PairDCAsigXY","",
                            GetVector(kMee), GetVector(kPairDCAsigXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY);
      histos->UserHistogram("Pair","InvMass_PairDCAabsXY","",
                            GetVector(kMee), GetVector(kPairDCAabsXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAabsXY);
      histos->UserHistogram("Pair","InvMass_PairLinDCAsigXY","",
                            GetVector(kMee), GetVector(kPairLinDCAsigXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAsigXY);
      histos->UserHistogram("Pair","InvMass_PairLinDCAabsXY","",
                            GetVector(kMee), GetVector(kPairLinDCAabsXY),
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAabsXY);
      // 3D
      histos->UserHistogram("Pair","InvMass_PairDCAsigXY_PairLinDCAsigXY","",
                            GetVector(kMee), GetVector(kPairDCAsigXY), GetVector(kPairLinDCAsigXY), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPairLinDCAsigXY);
      histos->UserHistogram("Pair","InvMass_PairDCAsigXY_OpeningAngle","",
                            GetVector(kMee), GetVector(kPairDCAsigXY), GetVector(kOpAng), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kOpeningAngle);
      histos->UserHistogram("Pair","InvMass_PairLinDCAsigXY_OpeningAngle","",
                            GetVector(kMee), GetVector(kPairLinDCAsigXY), GetVector(kOpAng), 
                            AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAsigXY, AliDielectronVarManager::kOpeningAngle);
    }
    histos->UserHistogram("Pair","InvMass_PairDCAsigXY_PairPt","",
                          GetVector(kMee), GetVector(kPairDCAsigXY), GetVector(kPtee), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPairDCAsigXY, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","InvMass_PairLinDCAsigXY_PairPt","",
                          GetVector(kMee), GetVector(kPairLinDCAsigXY), GetVector(kPtee), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPairLinDCAsigXY, AliDielectronVarManager::kPt);
    
    //centrality
    histos->UserHistogram("Pair","InvMass_Centrality",";Inv. Mass [GeV];Centrality;#pairs",
                          GetVector(kMee), BinsToVector(102,-1,101), 
                          AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);
    histos->UserHistogram("Pair","PairPt_Centrality",";Pair Pt [GeV];Centrality;#pairs",
                          GetVector(kPtee), BinsToVector(102,-1,101), 
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kCentrality);
    
    // the histograms "Pre" may produce a memory leak:
    //    W-TROOT::Append: Replacing existing TH1: Pt (Potential memory leak).
    //    <message repeated 4 times>
    
    //    //add histograms to Track classes
    //    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
    //    histos->UserHistogram("Pre","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
    //    histos->UserHistogram("Pre","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
    //    histos->UserHistogram("Pre","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
    //    
    //    histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
    //                          GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    //    histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
    //                          GetVector(kP2D), BinsToVector(120,0.,120.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    //    
    //    histos->UserHistogram("Pre","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
    //                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    //    histos->UserHistogram("Pre","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    //    histos->UserHistogram("Pre","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //    //histos->UserHistogram("Pre","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
    //    
    //    histos->UserHistogram("Pre","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    //    histos->UserHistogram("Pre","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    //    histos->UserHistogram("Pre","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    //    histos->UserHistogram("Pre","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
    //                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    
  }// end: (!fIsQATask)
  
  die->SetHistogramManager(histos);
}


//_______________________________________________________________________________________________
TVectorD* LMEECutLib::GetVector(Int_t var) 
{
  switch (var) 
  {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
    case kPhi2D:  return AliDielectronHelper::MakeLinBinning(120, 0., TMath::TwoPi());
    case kY3D:    return AliDielectronHelper::MakeLinBinning( 50,-2.5,2.5);
      
    case kSigmaEle:
      if (fIsQATask) return AliDielectronHelper::MakeLinBinning(100,-10.,10.);
      else          return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      if (fIsQATask) return AliDielectronHelper::MakeLinBinning(100,-20.,20.);
      else          return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);
    case kTPCdEdx:
      if (fIsQATask) return AliDielectronHelper::MakeLinBinning(120,  0.,120.);
      else          return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);
      
    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86, 
                                                                   0.90, 0.94, 0.98, 1.02, 1.06, 
                                                                   1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 
                                                                   3.10, 3.30, 3.50, 4.00, 4.50, 5.00 
                                                                   ");
    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50 
                                                                   ");
    case kPtee:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 
                                                                   3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 
                                                                   5.00, 6.00, 7.00, 8.00 
                                                                   ");
    case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 
                                                                   1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 
                                                                   2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 
                                                                   4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00 
                                                                   ");
      //2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 
      //2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 
      //3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 
      
    case kPairDCAsigXY:
    case kPairLinDCAsigXY: return AliDielectronHelper::MakeLinBinning(50, 0., 10.); // in sigma
    case kPairDCAabsXY:
    case kPairLinDCAabsXY: return AliDielectronHelper::MakeLogBinning(50, 0.0001, 1.); // in cm
      
    case kRuns:   return AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600");
      
    default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
      break;
  } 
  //if ( var.EqualTo("p_2D"      , kIgnoreCase) ) return AliDielectronHelper::MakeLinBinning(160,0.,8.);
}


//_______________________________________________________________________________________________
TVectorD* LMEECutLib::BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //  
  //  return vec;
}


//_______________________________________________________________________________________________
void LMEECutLib::InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,100,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kM,200,-0.01,3.99); //20Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,30.,50.,80.,100.");
  
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,160,0.,8.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,350,0.,700.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,60,0.,120.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);
  
  //only in this case write MC truth info
  if (die->GetHasMC()) {
    cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }
  
  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}
