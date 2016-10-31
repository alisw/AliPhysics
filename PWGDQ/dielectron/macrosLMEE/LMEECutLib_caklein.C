#ifndef LMEECutLib_caklein
#define LMEECutLib_caklein

class AnalysisCut{
public:
  AnalysisCut(){};

  //Setter
  void SetPIDAna(Int_t sPIDAna){PIDAna = sPIDAna;};
  void SetPIDPre(Int_t sPIDPre){PIDPre = sPIDPre;}
  void SetTrackSelectionAna(Int_t sTrackSelectionAna){TrackSelectionAna = sTrackSelectionAna;}
  void SetTrackSelectionPre(Int_t sTrackSelectionPre){TrackSelectionPre = sTrackSelectionPre;}
  void SetPairCutsAna(Int_t sPairCutsAna){PairCutsAna = sPairCutsAna;}
  void SetPairCutsPre(Int_t sPairCutsPre){PairCutsPre = sPairCutsPre;}

  void SetCentrality(Int_t sCentrality){Centrality = sCentrality;}
  void SetPreFilterType(Int_t sPreFilterType){PreFilterType = sPreFilterType;}
  void SetMixing(Int_t sMixing){Mixing = sMixing;}
  void SetESDTrackSelection(Int_t sESDTrackSelection){ESDTrackSelection = sESDTrackSelection;}


  //Getter
  Int_t GetPIDAna(){return PIDAna;}
  Int_t GetPIDPre(){return PIDPre;}
  Int_t GetTrackSelectionAna(){return TrackSelectionAna;}
  Int_t GetTrackSelectionPre(){return TrackSelectionPre;}
  Int_t GetPairCutsAna(){return PairCutsAna;}
  Int_t GetPairCutsPre(){return PairCutsPre;}

  Int_t GetCentrality(){return Centrality;}
  Int_t GetPreFilterType(){return PreFilterType;}
  Int_t GetMixing(){return Mixing;}
  Int_t GetESDTrackSelection(){return ESDTrackSelection;}
private:
  Int_t PIDAna;
  Int_t PIDPre;
  Int_t TrackSelectionAna;
  Int_t TrackSelectionPre;
  Int_t PairCutsAna;
  Int_t PairCutsPre;
  Int_t PreFilterType;
  Int_t Centrality;
  Int_t Mixing;
  Int_t ESDTrackSelection;
};

class LMEECutLib {
public:
  // Possible PID Settings
  enum LMEEPIDAna{
    kPbPb2015_pidV0_pt400,
    kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton,
    kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton2,
    kITSTPCTOFif_trkSPDfirst_kINT7_pt100_woPID,
    kPbPb2015_Pt400_TPCele_symITS_tightTOFif,
    kPbPb2015_Pt400_TPCele_symITS_tightTOFreq,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq,
    kPbPb2015_Pt400_ITSSA
  };
  enum LMEEPIDPre{
    kStandardPre
  };
  // Possible Track Selections
  enum LMEETrackSelectionAna{
    kV0,
    kSPDfirst,
    kSPDorSDD_1,
    kITSSA
  };
  enum LMEETrackSelectionPre{
    kPrefilter_cut1
  };
  enum LMEETrackCuts{
    kPbPb2015_V0_tight,
    kSPD_bit4,
    kSDD_bit6,
    kITSSA_bit1
  };
  enum LMEEPairCutsAna{
    kMC_pi0Dal, // havesamemother = true
    kPairCutsAna, // Cut off (theta < 0.05) && (Minv < 0.02)
    kNoPairCutsAna // No Cuts applied, since 18.02.2014
  };
  enum LMEEPairCutsPre{
    kInvM0to030MeV_OpAng0to060mrad,
    kNoPairCutsPre
  };
  enum LMEEPreFilterType{
    kPreFilterAllSigns,
    kPreFilterUnlikeOnly,
    kNoPreFilter
  };
  // Possible Centralty Selections
  enum LMEECentSel{
    kPbPbCentral,     // 0%-10%
    kPbPbMidCentral,  //10%-20%
    kPbPbSemiCentral, //20%-50%
    kPbPbPeripheral,  //50%-90%
    kPbPb_00to50,     //0%-50%
    kPbPb_00to90,     //0%-90%
    kPbPb_00to100
  };
  enum LMEEEventMixing{
    kEventMixing_1  // 5 classes Z-Vertex, 7 classes centrality
  };
  enum LMEEESDTrackSelection{
    kStandardESD
  };
  enum LMEEEventCut{
    kStandard
  };
  LMEECutLib() {}

  AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
  AliAnalysisCuts*            GetCentralityCuts(AnalysisCut AnaCut);
  AliDielectronTrackRotator*  GetTrackRotator(AnalysisCut AnaCut);
  AliDielectronMixingHandler* GetMixingHandler(AnalysisCut AnaCut);

  AliAnalysisCuts* GetPairCutsAna(AnalysisCut AnaCut, Int_t togglePC=0); //Bool_t togglePC=kFALSE
  AliAnalysisCuts* GetPairCutsPre(AnalysisCut AnaCut);

  AliAnalysisCuts* GetPIDCutsAna(AnalysisCut AnaCut);
  AliAnalysisCuts* GetPIDCutsPre(AnalysisCut AnaCut);

  AliAnalysisCuts* GetTrackSelectionAna(AnalysisCut AnaCut);
  AliAnalysisCuts* GetTrackSelectionPre(AnalysisCut AnaCut);
  AliDielectronCutGroup* GetTrackCuts(Int_t cutSet);
  AliAnalysisCuts* GetESDTrackCutsAna(AnalysisCut AnaCut);

  void SetEtaCorrection(AliDielectron *die, Int_t selPID, Int_t selCent, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/

};

void LMEECutLib::SetEtaCorrection(AliDielectron *die, Int_t selPID, Int_t selCent, Int_t corrZdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrection()\n";
  std::cout << " corrZdim = " << corrZdim << "\n";
  std::cout << " corrYdim = " << corrYdim << std::endl;
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

  if (corrZdim==AliDielectronVarManager::kRefMultTPConly) {
    if (selCent == kPbPb2011SemiCentral){
      fitMinZdim=400.; fitMaxZdim=1400;
    }
    fCntrdCorr = new TF2("fCntrdCorr", "[9]*([0]+[1]*x+[2]*x*x)+[10]*([3]+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y+[8]*y*y*y*y*y)+[11]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    fWdthCorr  = new TF2("fWdthCorr",  "[9]*([0]+[1]*x+[2]*x*x)+[10]*([3]+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y+[8]*y*y*y*y*y)+[11]", fitMinZdim, fitMaxZdim, fitMinEta, fitMaxEta);
    Double_t parCntrd[]={1.984011e-01, -2.538876e-04, 7.636316e-09, 2.731183e-01, 2.190149e-01, -3.437725e+00, -6.438784e-01, 5.317945e+00, 6.089976e-01, 1.018070e+00, 1.002105e+00, 2.428792e-02};
    Double_t parWdth[] ={1.116969e+00, 5.063439e-05, -5.526647e-09, 1.207224e+00, -3.633499e-02, -5.340433e-01, 1.557123e-01, 7.865346e-01, -1.619051e-01, 1.019816e+00, 1.004669e+00, -1.185922e+00};
    fCntrdCorr->SetParameters(parCntrd);
    fWdthCorr ->SetParameters(parWdth);
  }
  else {
    std::cout << " no eta correction applied!" << std::endl;
    return;
  }

  die->SetCentroidCorrFunction(fCntrdCorr, corrZdim, corrYdim);
  die->SetWidthCorrFunction(fWdthCorr, corrZdim, corrYdim);
  std::cout << " TPC PID eta correction loaded!!!" << std::endl;
  return;
}


// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask, currently to 'kPbPb2011_TPCTOF_Semi1'
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
    case kStandard:
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


//Selection of relatively 'flat' centralities
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(AnalysisCut AnaCut) {
  AliDielectronVarCuts* centCuts = 0x0;
  switch (AnaCut.GetCentrality()) {
    case kPbPbCentral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbCentral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.,10.);
      break;
    case kPbPbMidCentral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbMidCentral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,10.,20.);
      break;
    case kPbPbSemiCentral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbSemiCentral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,10.,50.);
      break;
    case kPbPbPeripheral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbPeripheral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,50.,90.);
      break;
    case kPbPb_00to50:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb_00to50");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.,50.);
      break;
    case kPbPb_00to90:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb_00to90");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.01,90.);
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


AliDielectronMixingHandler* LMEECutLib::GetMixingHandler(AnalysisCut AnaCut) {
  AliDielectronMixingHandler* mixingHandler = 0x0;
  switch (AnaCut.GetMixing()) {
    case kEventMixing_1:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,20,30,50,80");
      // now using TPC event plane, uncorrected. (also, the old phi range was wrong, now same effective binning.)
      mixingHandler->AddVariable(AliDielectronVarManager::kTPCrpH2uc, 6, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(15);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    default: cout << "No Mixer defined" << endl;
  }
  return mixingHandler;
}


//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsAna(AnalysisCut AnaCut, Int_t togglePC)  {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  switch (AnaCut.GetPairCutsAna()) {
    case kMC_pi0Dal /*kPbPb2011MC_pi0Dal_1*/:
      if (togglePC==1) {
        AliDielectronVarCuts* pairCutsMC =new AliDielectronVarCuts("pairCutsMC","pairCutsMC");
        pairCutsMC->AddCut(AliDielectronVarManager::kHaveSameMother, kTRUE);
        pairCuts = pairCutsMC;
        break;
      }
    case kPairCutsAna:
      //  AliDielectronVarCuts* pairCutsPhivGood =new AliDielectronVarCuts("pairCutsPhivGood","pairCutsPhivGood");
      //  pairCutsPhivGood->AddCut(AliDielectronVarManager::kPhivPair, 0.0, 2.0);
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
    case kNoPairCutsAna:
      cout << "No Pair Cuts used - ok " << endl; // since 18.02.2014
      break;

    default: cout << "No Pair Cuts defined " << endl;
  }
  return pairCuts;
}

//Pair Cuts for PREFILTER step
// cuts = REJECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsPre(AnalysisCut AnaCut)  {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  switch (AnaCut.GetPairCutsPre()) {
    case kInvM0to030MeV_OpAng0to060mrad:
      AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
      pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.03); // in upgrade: 0.01
      AliDielectronVarCuts* pairCutsOpAng =new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");
      pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.06); // in upgrade: 0.05

      AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
      pairCutsCG->AddCut(pairCutsInvM);
      pairCutsCG->AddCut(pairCutsOpAng);
      //pairCutsCG->AddCut(pairCutsPhiv);
      pairCuts = pairCutsCG;
      break;
    case kNoPairCutsPre:
      //[...] // PhiV and InvMass
    default: cout << "No Prefilter Pair Cuts defined " << endl;
  }
  return pairCuts;
}

AliAnalysisCuts* LMEECutLib::GetPIDCutsAna(AnalysisCut AnaCut) {
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
  //TPC: electron inclusion asymmetric
  //     pion     exclusion 3sigma
  //ITS: electron inclusion asymmetric OVER FULL MOMENTUM RANGE
  //TOF: electron inclusion 3sigma - BUT ONLY IF AVAILABLE

  // PID for V0 task
  AliDielectronPID *pid_TOFonly = new AliDielectronPID("pid_TOFonly","pid_TOFonly");
  pid_TOFonly->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE);
  pid_TOFonly->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 2. , 0. ,100., kTRUE);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif","pid_TPCele_AsymITS_tightTOFif");
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFreq = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFreq","pid_TPCele_AsymITS_tightTOFreq");
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);


  // eta range:
  AliDielectronVarCuts *etaRange090 = new AliDielectronVarCuts("etaRange090","etaRange090");
  etaRange090->AddCut(AliDielectronVarManager::kEta, -0.90, 0.90);
  AliDielectronVarCuts *etaRange084 = new AliDielectronVarCuts("etaRange084","etaRange084");
  etaRange084->AddCut(AliDielectronVarManager::kEta, -0.84, 0.84);
  AliDielectronVarCuts *etaRange080 = new AliDielectronVarCuts("etaRange080","etaRange080");
  etaRange080->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
  AliDielectronVarCuts *etaRange076 = new AliDielectronVarCuts("etaRange076","etaRange076");
  etaRange076->AddCut(AliDielectronVarManager::kEta, -0.76, 0.76);
  // pt range:
  AliDielectronVarCuts *ptRange500to3500 = new AliDielectronVarCuts("ptRange500to3500","ptRange500to3500");
  ptRange500to3500->AddCut(AliDielectronVarManager::kPt, 0.5, 3.5);
  AliDielectronVarCuts *ptRange400to3500 = new AliDielectronVarCuts("ptRange400to3500","ptRange400to3500");
  ptRange400to3500->AddCut(AliDielectronVarManager::kPt, 0.4, 3.5);
  AliDielectronVarCuts *ptRange300to3500 = new AliDielectronVarCuts("ptRange300to3500","ptRange300to3500");
  ptRange300to3500->AddCut(AliDielectronVarManager::kPt, 0.3, 3.5);
  AliDielectronVarCuts *ptRange200to3500 = new AliDielectronVarCuts("ptRange200to3500","ptRange200to3500");
  ptRange200to3500->AddCut(AliDielectronVarManager::kPt, 0.2, 3.5);
  AliDielectronVarCuts *ptRange100to3500 = new AliDielectronVarCuts("ptRange100to3500","ptRange100to3500");
  ptRange100to3500->AddCut(AliDielectronVarManager::kPt, 0.1, 3.5);



  AliDielectronPID *PID_cutoff_pion_kaon_proton = new AliDielectronPID("PID_cutoff_pion_kaon_proton","PID_cutoff_pion_kaon_proton");
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -2.0, 2. , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -2.0, 2. , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PID_cutoff_pion_kaon_proton_EleIncl = new AliDielectronPID("PID_cutoff_pion_kaon_proton_EleIncl","PID_cutoff_pion_kaon_proton_EleIncl");
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.5, 3. , 0. ,100., kFALSE);
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronCutGroup* PID_cutoff_pion_kaon_proton_cg = new AliDielectronCutGroup("PID_cutoff_pion_kaon_proton_cg","PID_cutoff_pion_kaon_proton_cg",AliDielectronCutGroup::kCompOR);
  PID_cutoff_pion_kaon_proton_cg->AddCut(PID_cutoff_pion_kaon_proton);
  PID_cutoff_pion_kaon_proton_cg->AddCut(PID_cutoff_pion_kaon_proton_EleIncl);


  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------
  switch (AnaCut.GetPIDAna()) {
    case kITSTPCTOFif_trkSPDfirst_kINT7_pt100_woPID:
    // NO PID
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange100to3500);
      // cgPIDCutsAna->AddCut(pidTPCITS_TOFif2);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;

    case kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton:
    // Cut out pion/kaon/proton band LHC15o but refilled when particle in TOF electron band
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(PID_cutoff_pion_kaon_proton_cg);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;

    case kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton2:
    // Cut out pion/kaon/proton band LHC15o
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(PID_cutoff_pion_kaon_proton);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;

    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif:
    // ITS & TOFifavailable & TPC sigma cut for LHC15o
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pid_TPCele_AsymITS_tightTOFif);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pid_TPCele_AsymITS_tightTOFreq);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;

    case kPbPb2015_Pt400_ITSSA:
    // only ITSSA tracks
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2015_pidV0_pt400:
    // Only TOF PID for V0 selection
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to3500);
      cgPIDCutsAna->AddCut(pid_TOFonly);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;
    default: cout << "No Analysis PID Cut defined " << endl;
  }
  return pidCuts;
}

AliAnalysisCuts* LMEECutLib::GetTrackSelectionAna(AnalysisCut AnaCut) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (AnaCut.GetTrackSelectionAna()) {
    case kV0:
      trackCuts = GetTrackCuts(kPbPb2015_V0_tight);
      break;
    case kSPDfirst:
      trackCuts = GetTrackCuts(kSPD_bit4);
      break;
    case kITSSA:
      trackCuts = GetTrackCuts(kITSSA_bit1);
      break;
    case kSPDorSDD_1:
      AliDielectronCutGroup* cgTrackSelAna = new AliDielectronCutGroup("cgTrackSelAna","cgTrackSelAna",AliDielectronCutGroup::kCompOR);
      cgTrackSelAna->AddCut(GetTrackCuts(kSPD_bit4));
      cgTrackSelAna->AddCut(GetTrackCuts(kSDD_bit6));
      trackCuts = cgTrackSelAna;
      break;

    default: cout << "No Analysis Track Selection defined " << endl;
  }
  return trackCuts;
}

AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (cutSet) {
    case kPbPb2015_V0_tight:
      // primarily meant for inclusion, for quite pure sample...
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
      // gammaV0Cuts->SetExcludeTracks(kTRUE);
      gammaV0Cuts->SetExcludeTracks(kFALSE);
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      cgTrackCutsV0select->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsV0select;
      break;

    case kITSSA_bit1:
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     4.0, 100.0); // means at least 2 with PID
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<1); // ITSSA
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
    case kSPD_bit4:
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
    case kSDD_bit6:  // CHECK IN PATRICKS LMEELIB
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0); // means at least 3 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<6); //GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first

      cgTrackCutsAnaSDDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSDDfirst","cgTrackCutsAnaSDDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSDDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSDDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSDDfirst;
      break;

    default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
}


//Relaxed PID cuts for additional rejectin step, do not use blindly
AliAnalysisCuts* LMEECutLib::GetPIDCutsPre(AnalysisCut AnaCut) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;
  switch (AnaCut.GetPIDPre()) {
    case kStandardPre:

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
      pidITSTPCTOFpre->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 2., 0. ,  2., kFALSE);
      // TOF will be used if available, and with pt instead of p:
      //  pidITSTPCTOFpre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0.4,100., kFALSE,
      //                          AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
      cgITSTPCTOFpre->AddCut(pidITSTPCTOFpre);
      cgITSTPCTOFpre->AddCut(etaRangePre1);
      cgITSTPCTOFpre->AddCut(ptRangePre1);
      cgITSTPCTOFpre->AddCut(GetTrackSelectionPre(AnaCut));

      AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
      // in case the prefilter cuts do not include all needed global tracks.
      cgInitialTrackFilter->AddCut(GetPIDCutsAna(AnaCut));
      cgInitialTrackFilter->AddCut(cgITSTPCTOFpre);

      pidCuts = cgInitialTrackFilter;   // kCompOR works!!! <- checked with 'SetNoPairing()' and commented out 'GetPIDCutsAna(selectedPID)'
      //cout << " ========== pidCuts prefilter: ========== " << endl;
      //pidCuts->Print();
      break;

    case default:
    default: cout << "No Prefilter PID Cut defined " << endl;
  }
  return pidCuts;
}


//Possibly different cut sets for Prefilter step
//Not used at the moment
AliAnalysisCuts* LMEECutLib::GetTrackSelectionPre(AnalysisCut AnaCut) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (AnaCut.GetTrackSelectionPre()) {
    case kPrefilter_cut1:
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
AliAnalysisCuts* LMEECutLib::GetESDTrackCutsAna(AnalysisCut AnaCut) {
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetESDTrackCutsAna()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> Setting ESD Track Cuts >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> ( do we run on ESD?! ) >>>>>>>>>>>>>>>>>>>>>> " << endl;
  //cout << " >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliESDtrackCuts* esdTrackCutsH = 0x0;
  switch (AnaCut.GetESDTrackSelection()) {
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

#endif
