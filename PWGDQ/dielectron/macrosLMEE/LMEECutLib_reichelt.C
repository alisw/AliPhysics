class LMEECutLib {
  
public:

	static  enum enCentSel {
	  kPbPb2011Central,
    kPbPb2011_00to10,
    kPbPb2011MidCentral,
    kPbPb2011_10to20,
		kPbPb2011SemiCentral,
    kPbPb2011_20to50,
    kPbPb2011_10to50,
    kPbPb2011_00to50,
		kCENTSELMAX
	};
  
	static  enum enPairCut {
    kPairCut_mee10_theta30,
    kPairCut_mee20_theta20, // not yet done (maybe for prefilter with SPD+SDD or as default PairCut for final pairs)
	  kPairCut_mee20_theta50,
    kPairCut_mee30_theta60,
    kPairCut_mee40_theta80,
    kPairCut_mee60_theta100,
    kPairCut_phiv157_mee40,
    kPairCut_phiv157_mee60,
    kPairCut_phiv157_mee80,
    kPairCut_phiv157_mee100,
    kPairCut_phiv236_mee40,
    kPairCut_phiv236_mee60,
    kPairCut_phiv236_mee80,
    kPairCut_phiv236_mee100,
		kPairCut_OFF
	};
  
//	PIDCutsAna
//	PIDCutsPre
//	TrackCutsAna
//	TrackCutsPre
  
  
	static  enum LMMECutSet {
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_TESTING,
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
    kPbPb2011_pidITSTPC_trkSPDfirst_3,            // (cutset w/o pairing)
    kPbPb2011_pidTPC_trkSPDfirst_3,               // (cutset w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose, // (cutset w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose, // (cutset w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose,      // (cutset w/o pairing)
    kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose,      // (cutset w/o pairing)
    kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1,
    kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1,       // Cutset for Technical Preliminaries for QM2014 (no prefilter used!)
    kPbPb2011_pidTPCTOF_trkSPDorSDD_1,
    kPbPb2011_pidTPCTOF_trkSPDfirst_1,
    kPbPb2011PID_ITSTPCTOFif2,        // (NO FULL CUTSET)
    kPbPb2011PID_TPCTOF3,             // (NO FULL CUTSET)
    kPbPb2011TRK_SDDfirstSPDnone,     // (NO FULL CUTSET) complimentary tracks, strictly without SPD, to be combined with others!
    kPbPb2011TRK_SPDfirst,            // (NO FULL CUTSET) main track selection, with SPD first
    kPbPb2011TRK_SDDfirstSPDnone4cls, // (NO FULL CUTSET) complimentary tracks, strictly without SPD, to be combined with others!
    kPbPb2011TRK_SPDfirst5cls,        // (NO FULL CUTSET) main track selection, with SPD first
    kPbPb2011_TPCITS_TOFif1,
    kPbPb2011_TPCTOF_Semi2, // changed PairCutsAna from PhiV to OpeningAngle. prefilter cuts renewed (if applicable)
    // following cutsets are not complete anymore!
    kPbPb2011_TPCTOF_Semi1, // old prefilter cuts (leg & pair), some are confusing
    kPbPb2011NoPID, // pairing disabled in config
    kPbPb2011TPCandTOF, // this was the final one activated by Christoph!
    kPbPb2011TPCandTOFHPT,
		kPbPb2011TPC, //TOF required, more relaxed cut on TPC
		kPbPb2011TPCandTOFwide, //TOF required, more relaxed cut on TPC
		kPbPb2011TPCorTOF,
		kpp2010TPCandTOF,
		kpp2010TPCorTOF,
		kCUTSETMAX
	};
  
	static  enum enCutType {
	  kInclude = 0,
    kExclude = 1,
		kCUTTYPEMAX
	};
  
	//char* LMEECutNames[kCUTSETMAX] = { "PbPb2011TPCandTOF","PbPb2011TPCorTOF"};
  
  
	LMEECutLib();
  
	AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
	AliAnalysisCuts*            GetCentralityCuts();
	AliDielectronTrackRotator*  GetTrackRotator();
	AliDielectronMixingHandler* GetMixingHandler();
  
	AliAnalysisCuts* GetPairCutsAna(Int_t togglePC=0); //Bool_t togglePC=kFALSE
	AliAnalysisCuts* GetPairCutsPre();
  
	AliAnalysisCuts* GetPIDCutsAna();
	AliAnalysisCuts* GetPIDCutsPre();
  
	AliAnalysisCuts* GetTrackCutsAna(Int_t cutSet=-1, Int_t doExclusion=0);
	AliAnalysisCuts* GetTrackCutsPre(Int_t cutSet=-1);
	AliAnalysisCuts* GetESDTrackCutsAna();
  
  void SetEtaCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  
  
  Int_t selectedCentrality;
  Int_t selectedPIDAna;
  Int_t selectedPIDPre;
  Int_t selectedTrackAna;
  Int_t selectedTrackPre;
  Int_t selectedPairCutsPre;
  Int_t selectedPairCutsAna;
  
};


LMEECutLib::LMEECutLib() :
selectedCentrality(-1),
selectedPIDAna(-1),
selectedPIDPre(-1),
selectedTrackAna(-1),
selectedTrackPre(-1),
selectedPairCutsPre(LMEECutLib::kPairCut_OFF),
selectedPairCutsAna(LMEECutLib::kPairCut_OFF)
{
  // Constructor
}


void LMEECutLib::SetEtaCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  printf("starting LMEECutLib::SetEtaCorrection()\n");
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


// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask, currently to 'kPbPb2011_TPCTOF_Semi1'
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
    case kPbPb2011_TPCTOF_Semi1:
      //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
      eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
      eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
      //eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTPC); // AOD
      //           eventCuts->SetCentralityRange(0.0,80.0);
      eventCuts->SetRequireVertex();
      eventCuts->SetMinVtxContributors(1);
      eventCuts->SetVertexZ(-10.,10.);
      break;
    default: cout << "No Event Cut defined" << endl;
  }
  return eventCuts;
}


//Selection of relatively 'flat' centralities is a bit difficult...
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
AliAnalysisCuts* LMEECutLib::GetPairCutsAna(Int_t togglePC)  {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  switch (selectedPairCutsAna) {
    case kPbPb2011MC_pi0Dal_1:
      if (togglePC==1) {
        AliDielectronVarCuts* pairCutsMC =new AliDielectronVarCuts("pairCutsMC","pairCutsMC");
        pairCutsMC->AddCut(AliDielectronVarManager::kHaveSameMother, kTRUE);
        pairCuts = pairCutsMC;
        break;
      }
      //else use what comes below...
      
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      cout << "No Pair Cuts used - ok " << endl; // since 18.02.2014
      break;
      
    case kPbPb2011_TPCITS_TOFif1:
    case kPbPb2011_TPCTOF_Semi2:
      //        AliDielectronVarCuts* pairCutsPhivGood =new AliDielectronVarCuts("pairCutsPhivGood","pairCutsPhivGood");
      //        pairCutsPhivGood->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.0); 
      AliDielectronVarCuts* pairCutsOpAngGood =new AliDielectronVarCuts("pairCutsOpAngGood","pairCutsOpAngGood");
      pairCutsOpAngGood->AddCut(AliDielectronVarManager::kOpeningAngle, 0.05, 999.); // in upgrade: 0.05
      AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
      pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.02); // in upgrade: 0.01
      AliDielectronVarCuts* pairCutsInvMgood =new AliDielectronVarCuts("pairCutsInvMgood","pairCutsInvMgood");
      pairCutsInvMgood->AddCut(AliDielectronVarManager::kM, 0.02, 99999.);
      
      AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
      pairCutsCG->AddCut(pairCutsInvM);
      pairCutsCG->AddCut(pairCutsOpAngGood);
      //        pairCutsCG->AddCut(pairCutsPhivGood);
      
      AliDielectronCutGroup* pairCutsCG2 =new AliDielectronCutGroup("pairCutsCG2","pairCutsCG2",AliDielectronCutGroup::kCompOR);
      pairCutsCG2->AddCut(pairCutsInvMgood);
      pairCutsCG2->AddCut(pairCutsCG);
      pairCuts = pairCutsCG2;
      break;
      
    case kPbPb2011_TPCTOF_Semi1:
      //[...] // PhiV and InvMass
    default: cout << "No Pair Cuts defined " << endl;
  }
  return pairCuts;
}


//Pair Cuts for PREFILTER step
// cuts = REJECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsPre()  {  
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  AliDielectronVarCuts* pairVarCuts =new AliDielectronVarCuts("pairVarCuts","pairVarCuts");
  
  switch (selectedPairCutsPre) {
      
      // Mee and opening angle
    case kPairCut_mee10_theta30:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.01);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.03);
      pairCuts = pairVarCuts;
      break;
    case kPairCut_mee20_theta20:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.02);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.02);
      pairCuts = pairVarCuts;
      break;
    case kPairCut_mee20_theta50:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.02);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.05);
      pairCuts = pairVarCuts;
      break;
    case kPairCut_mee30_theta60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.03);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.06);
      pairCuts = pairVarCuts;
      break;
    case kPairCut_mee40_theta80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.08);
      pairCuts = pairVarCuts;
      break;
    case kPairCut_mee60_theta100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.10);
      pairCuts = pairVarCuts;
      break;
      
      // Mee and phiv
    case kPairCut_phiv157_mee40:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv157_mee60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv157_mee80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.08);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv157_mee100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.10);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.5*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv236_mee40:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.04);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv236_mee60:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.06);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv236_mee80:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.08);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
    case kPairCut_phiv236_mee100:
      pairVarCuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.10);
      pairVarCuts->AddCut(AliDielectronVarManager::kPhivPair, 0.75*TMath::Pi(), 3.2); 
      pairCuts = pairVarCuts;
      break;
      
    case kPairCut_OFF:
    default: cout << "No Prefilter Pair Cuts defined " << endl;
  } 
  return pairCuts;
}



AliAnalysisCuts* LMEECutLib::GetPIDCutsAna() {
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
  
  AliDielectronPID *pid_V0select_2 = new AliDielectronPID("pid_V0select_2","pid_V0select_2");
  pid_V0select_2->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-12. ,20. , 0. ,100., kFALSE);
  // no TOF to really see what gets rejected...
  
  // eta range:
  AliDielectronVarCuts *etaRange090 = new AliDielectronVarCuts("etaRange090","etaRange090");
  etaRange090->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
  AliDielectronVarCuts *etaRange084 = new AliDielectronVarCuts("etaRange084","etaRange084");
  etaRange084->AddCut(AliDielectronVarManager::kEta, -0.84, 0.84);
  AliDielectronVarCuts *etaRange080 = new AliDielectronVarCuts("etaRange080","etaRange080"); //changed from 0.76 to 0.8 on 2014-08-12!
  etaRange080->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
  // pt range:
  AliDielectronVarCuts *ptRange200to3500 = new AliDielectronVarCuts("ptRange200to3500","ptRange200to3500");
  ptRange200to3500->AddCut(AliDielectronVarManager::kPt, .2, 3.5);
  AliDielectronVarCuts *ptRange400to3500 = new AliDielectronVarCuts("ptRange400to3500","ptRange400to3500");
  ptRange400to3500->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
  
  
  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------
  switch (selectedPIDAna) {
      
    case kPbPb2011_V0select_2_looseNoTOF:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange090);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pid_V0select_2);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      
    case kPbPb2011_V0select_1_Arm:
    case kPbPb2011_V0select_1:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange090);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pid_V0select_1);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight: // tighter "ITSTPCTOFif" PID
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCITS_TOFif56);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl:
    case kPbPb2011MC_pi0Dal_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4: // regular "ITSTPCTOFif" PID
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCITS_TOFif2);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_TESTING:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCITS_TOFif2);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange090);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCTOF_Semi1);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      
      /* PID for cutsets without pairing! */
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCITS_3);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2011_pidTPC_trkSPDfirst_3:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange090);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPC_3);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose: // loose "ITSTPCTOFif" PID - for 2D contamination study
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCITS_TOFif_LOOSE);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange090);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCTOF_Semi_LOOSE);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
      /* end - cutsets without pairing */
      
    case kPbPb2011_TPCITS_TOFif1:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange084); // was 0.84 -> not ideal
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCITS_TOFif1);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() );
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2011_TPCTOF_Semi2:
    case kPbPb2011_TPCTOF_Semi1:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange084);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pidTPCTOF_Semi1);
      cgPIDCutsAna->AddCut( GetTrackCutsAna() ); // for 'kPbPb2011_TPCTOF_Semi1', this was called in the Config
      pidCuts = cgPIDCutsAna;
      break;
      //[...]
    default: cout << "No Analysis PID Cut defined " << endl;
  }
  return pidCuts;
}


//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

AliAnalysisCuts* LMEECutLib::GetTrackCutsAna(Int_t cutSet, Int_t doExclusion) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  
  // for the default function call, pick the cutset according to 'selectedTrackAna', which is set via the Config file.
  if (cutSet<0) { cutSet = selectedTrackAna; }
  
  switch (cutSet) {
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_7_V0excl:
      // combine 
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1));
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011V0_2_loose, kExclude));
      trackCuts = cgTrackCutsAna;
      break;
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_7_V0excl:
      // combine 
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1));
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011V0_2_loose, kExclude));
      //// instead, only FOR DOING A COMPARISON: ////
      ////      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1));
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
      cgTrackCutsMC->AddCut(GetTrackCutsAna(kPbPb2011TRK_SPDfirst));
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
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011TRK_SPDfirst));         // typical trackcuts with requirement of SPD
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011TRK_SDDfirstSPDnone)); // new additional trackcuts with SDD instead of SPD
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
    case kPbPb2011TRK_SPDfirst: // main track selection, now closer to what Hongyan does...
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
      
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_TESTING:
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
//      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
//      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
      
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      
      
      
      
    case kPbPb2011TRK_SDDfirstSPDnone: // complimentary tracks, strictly without SPD, to be combined with others!
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
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011TRK_SPDfirst5cls));         // typical trackcuts with requirement of SPD
      cgTrackCutsAna->AddCut(GetTrackCutsAna(kPbPb2011TRK_SDDfirstSPDnone4cls)); // new additional trackcuts with SDD instead of SPD
      trackCuts = cgTrackCutsAna;
      break;
      
      //----------
      // MAIN settings - single trackset - variation 1
      //----------
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
      //----------
    case kPbPb2011TRK_SPDfirst5cls: // main track selection, 5+ ITS clusters
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
      
    case kPbPb2011TRK_SDDfirstSPDnone4cls: // complimentary tracks, 4+ ITS clusters, strictly without SPD, to be combined with others!
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
      
      // ==========
    case kPbPb2011_TPCITS_TOFif1:
    case kPbPb2011_TPCTOF_Semi2: // no pt and eta ranges in the trackcuts anymore!
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      // trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.84, 0.84); // eta commented out later. (05.02.2014)
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     110.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.0);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(16); //does nothing for ESDs // 16=2^4 -> filter bit 4!
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //function in AliDielectronTrackCuts
      
      cgTrackCutsAna = new AliDielectronCutGroup("cgTrackCutsAna","cgTrackCutsAna",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAna->AddCut(trackCutsDiel);
      cgTrackCutsAna->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAna;
      break;
      
    case kPbPb2011_TPCTOF_Semi1:
      //[...]
    default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
} 



//Relaxed PID cuts for additional rejectin step, do not use blindly
AliAnalysisCuts* LMEECutLib::GetPIDCutsPre() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;
  switch (selectedPIDPre) {
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
    case kPbPb2011_pidTPC_trkSPDfirst_3:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
    case kPbPb2011_TPCITS_TOFif1:
    case kPbPb2011_TPCTOF_Semi2:
      
      // eta range:
      AliDielectronVarCuts *etaRangePre1 = new AliDielectronVarCuts("etaRangePre1","etaRangePre1");
      etaRangePre1->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
      // pt range:
      AliDielectronVarCuts *ptRangePre1 = new AliDielectronVarCuts("ptRangePre1","ptRangePre1");
      ptRangePre1->AddCut(AliDielectronVarManager::kPt, .2, 3.5); // 0.2 is realistic. turnon at ~180MeV
      //AliDielectronVarCuts *ptRangePre2 = new AliDielectronVarCuts("ptRangePre2","ptRangePre2");
      //ptRangePre2->AddCut(AliDielectronVarManager::kPt, .4, 3.5);
      //AliDielectronVarCuts *ptRangePre3 = new AliDielectronVarCuts("ptRangePre3","ptRangePre3");
      //ptRangePre3->AddCut(AliDielectronVarManager::kPt, 0.05, 1.5);
      
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
      cgITSTPCTOFpre->AddCut( GetTrackCutsAna() );
      
      //        AliDielectronCutGroup* cgTPCpre = new AliDielectronCutGroup("cgTPCpre","cgTPCpre",AliDielectronCutGroup::kCompAND);
      //        AliDielectronPID *pidTPCpre = new AliDielectronPID("pidTPCpre","pidTPCpre");
      //        pidTPCpre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE);
      //        pidTPCpre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
      //        // TOF will be used if available, and with pt instead of p:
      //        pidTPCpre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0.4,5.  , kFALSE, 
      //                          AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
      //        cgTPCpre->AddCut(pidTPCpre);
      //        cgTPCpre->AddCut(etaRangePre1);
      //        cgTPCpre->AddCut(ptRangePre2);
      //        cgTPCpre->AddCut( GetTrackCutsAna() );
      
      //        AliDielectronCutGroup* cgITSSA = new AliDielectronCutGroup("cgITSSA","cgITSSA",AliDielectronCutGroup::kCompAND);
      //        AliDielectronPID *pidITSSA = new AliDielectronPID("pidITSSA","pidITSSA");
      //        pidITSSA->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 3.);
      //        // this means that very many pions will be used for rejection!
      //        cgITSSA->AddCut(pidITSSA);
      //        cgITSSA->AddCut(etaRangePre1);
      //        cgITSSA->AddCut(ptRangePre3);
      //        cgITSSA->AddCut( GetTrackCutsPre() );
      
      AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
      cgInitialTrackFilter->AddCut( GetPIDCutsAna() ); // in case the prefilter cuts do not include all needed global tracks.
      cgInitialTrackFilter->AddCut(cgITSTPCTOFpre);
      //cgInitialTrackFilter->AddCut(cgTPCpre);
      //cgInitialTrackFilter->AddCut(cgITSSA);
      pidCuts = cgInitialTrackFilter;   // kCompOR works!!!
      //cout << " ========== pidCuts prefilter: ========== " << endl;
      //pidCuts->Print();
      break;
      
    case kPbPb2011_TPCTOF_Semi1:
      //[...]
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
    case kPbPb2011_pidITSTPC_trkSPDfirst_3:
    case kPbPb2011_pidTPC_trkSPDfirst_3:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_2_loose:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_2_loose:
    case kPbPb2011_pidTPCTOF_trkSPDfirst_2_loose:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPD5orSDD4cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst5cls_6_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDorSDD_5_tight:
    case kPbPb2011_pidITS2gevTPCTOFif_trkSPDfirst_5_tight:
    case kPbPb2011_pidITSTPCTOFif_trkSPD5orSDD4cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst5cls_4:
    case kPbPb2011_pidITSTPCTOFif_trkSPDorSDD_1:
    case kPbPb2011_pidITSTPCTOFif_trkSPDfirst_1:
    case kPbPb2011_pidTPCTOF_trkSPDorSDD_1:
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
      
    case kPbPb2011_pidTPCTOF_trkSPDfirst_1:
      //[...]
    case kPbPb2011_TPCITS_TOFif1:
    case kPbPb2011_TPCTOF_Semi2:
      //[...]
    case kPbPb2011_TPCTOF_Semi1:
      //[...]
    default: cout << "No Prefilter Track Cut defined " << endl;
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
