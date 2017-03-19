#ifndef LMEECutLib_amichalik
#define LMEECutLib_amichalik

class AnalysisCut{
public:
  AnalysisCut(){};

  //Setter
  void SetPIDAna(Int_t sPIDAna){PIDAna = sPIDAna;};
  void SetPIDPre(Int_t sPIDPre){PIDPre = sPIDPre;}
  void SetTrackSelectionAna(Int_t sTrackSelectionAna){TrackSelectionAna = sTrackSelectionAna;}
  void SetAdditionalTrackCuts(Int_t sAdditionalTrackCuts){AdditionalTrackCuts = sAdditionalTrackCuts;}
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
  Int_t GetAdditionalTrackCuts(){return AdditionalTrackCuts;}
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
  Int_t AdditionalTrackCuts;
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
    kPDGLepton,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif
  };
  enum LMEEPIDPre{
    kStandardPre
  };
  // Possible Track Selections
  enum LMEETrackSelectionAna{
    kSPDfirst
  };
  enum LMEETrackSelectionPre{
    kPrefilter_cut1
  };
  enum LMEETrackCuts{
    kSPD_bit4
  };
  enum LMEEAdditionalCuts{
    // kOmega,
    // kEta,
    // kEtaprime,
    // kPhi,
    // kRho,
    // kJPsi,
    // kPion,
    // kAll
  };
  enum LMEEPairCutsAna{
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
    kPbPbSemiCentral, //20%-50%
    kPbPbSemiCentralRun1,
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
  AliDielectronCutGroup* GetAdditionalTrackCuts(AnalysisCut AnaCut);
  AliAnalysisCuts* GetESDTrackCutsAna(AnalysisCut AnaCut);


};





// #######################################################################################################################
AliAnalysisCuts* LMEECutLib::GetPIDCutsAna(AnalysisCut AnaCut) {
  AliAnalysisCuts* pidCuts=0x0;


  // eta range:
  AliDielectronVarCuts *etaRange076 = new AliDielectronVarCuts("etaRange076","etaRange076");
  etaRange076->AddCut(AliDielectronVarManager::kEta, -0.76, 0.76);
  AliDielectronVarCuts *etaRange080 = new AliDielectronVarCuts("etaRange076","etaRange076");
  etaRange080->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
  // pt range:
  AliDielectronVarCuts *ptRange400to3500 = new AliDielectronVarCuts("ptRange400to3500","ptRange400to3500");
  ptRange400to3500->AddCut(AliDielectronVarManager::kPt, 0.4, 3.5);
  AliDielectronVarCuts *ptRange200to3500 = new AliDielectronVarCuts("ptRange200to3500","ptRange200to3500");
  ptRange200to3500->AddCut(AliDielectronVarManager::kPt, 0.2, 3.5);


  // PID
  AliDielectronVarCuts *electrons = new AliDielectronVarCuts("electrons","electrons");
  electrons->AddCut(AliDielectronVarManager::kPdgCode, 11.);
  AliDielectronVarCuts *positrons = new AliDielectronVarCuts("positrons","positrons");
  positrons->AddCut(AliDielectronVarManager::kPdgCode, -11.);

  AliDielectronVarCuts *pion = new AliDielectronVarCuts("pion","pion");
  pion->AddCut(AliDielectronVarManager::kPdgCode, 211.);
  AliDielectronVarCuts *anti_pion = new AliDielectronVarCuts("anti_pion","anti_pion");
  anti_pion->AddCut(AliDielectronVarManager::kPdgCode, -211.);

  AliDielectronVarCuts *kaon = new AliDielectronVarCuts("kaon","kaon");
  kaon->AddCut(AliDielectronVarManager::kPdgCode, 321.);
  AliDielectronVarCuts *anti_kaon = new AliDielectronVarCuts("anti_kaon","anti_kaon");
  anti_kaon->AddCut(AliDielectronVarManager::kPdgCode, -321.);

  AliDielectronVarCuts *proton = new AliDielectronVarCuts("proton","proton");
  proton->AddCut(AliDielectronVarManager::kPdgCode, 2212.);
  AliDielectronVarCuts *anti_proton = new AliDielectronVarCuts("anti_proton","anti_proton");
  anti_proton->AddCut(AliDielectronVarManager::kPdgCode, -2212.);

  AliDielectronVarCuts *gamma = new AliDielectronVarCuts("gamma","gamma");
  gamma->AddCut(AliDielectronVarManager::kPdgCode, 22.);


  AliDielectronCutGroup* leptons = new AliDielectronCutGroup("leptons","leptons",AliDielectronCutGroup::kCompOR);
  leptons->AddCut(electrons);
  leptons->AddCut(positrons);
  // leptons->AddCut(pion);
  // leptons->AddCut(anti_pion);
  // leptons->AddCut(kaon);
  // leptons->AddCut(anti_kaon);
  // leptons->AddCut(proton);
  // leptons->AddCut(anti_proton);
  // leptons->AddCut(gamma);



  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif","pid_TPCele_AsymITS_tightTOFif");
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  switch (AnaCut.GetPIDAna()){
    case kPDGLepton:
    // APPLY ALL THE CUTS HERE
    AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
    cgPIDCutsAna->AddCut(etaRange080);
    cgPIDCutsAna->AddCut(ptRange400to3500);
    cgPIDCutsAna->AddCut(leptons);
    cgPIDCutsAna->AddCut(GetTrackCuts(kSPD_bit4));
    pidCuts = cgPIDCutsAna;
    break;

    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif:
    // APPLY ALL THE CUTS HERE
    AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
    cgPIDCutsAna->AddCut(etaRange080);
    cgPIDCutsAna->AddCut(ptRange400to3500);
    cgPIDCutsAna->AddCut(pid_TPCele_AsymITS_tightTOFif);
    cgPIDCutsAna->AddCut(GetTrackCuts(kSPD_bit4));
    pidCuts = cgPIDCutsAna;
    break;
  }
  return pidCuts;
}


// #######################################################################################################################
// AliDielectronCutGroup* LMEECutLib::GetAdditionalTrackCuts(AnalysisCut AnaCut) {
//   cout << " >>>>>>>>>>>>>>>>>>>>>> GetAdditionalTrackCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
//   // AliDielectronCutGroup* trackCuts=0x0;
//
//   AliDielectronCutGroup* cgTrackCuts = new AliDielectronCutGroup("cgTrackCuts","cgTrackCuts",AliDielectronCutGroup::kCompOR);
//
//   switch (AnaCut.GetAdditionalTrackCuts()) {
//     case kPion:
//     std::cout << "Pions selected" << std::endl;
//     AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
//     trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  111.);
//     AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
//     trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -111.);
//     cgTrackCuts->AddCut(trackCutsMother1);
//     cgTrackCuts->AddCut(trackCutsMother2);
//     break;
//     case kOmega:
//       std::cout << "Omegas selected" << std::endl;
//       AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
//       trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  223.);
//       AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
//       trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -223.);
//       cgTrackCuts->AddCut(trackCutsMother1);
//       cgTrackCuts->AddCut(trackCutsMother2);
//     break;
    // case kEta:
    //   std::cout << "Etas selected" << std::endl;
    //   AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
    //   trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  221.);
    //   AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
    //   trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -221.);
    //   cgTrackCuts->AddCut(trackCutsMother1);
    //   cgTrackCuts->AddCut(trackCutsMother2);
    // break;
    // case kEtaprime:
    //   std::cout << "Etas selected" << std::endl;
    //   AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
    //   trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  331.);
    //   AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
    //   trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -331.);
    //   cgTrackCuts->AddCut(trackCutsMother1);
    //   cgTrackCuts->AddCut(trackCutsMother2);
    // break;
    // case kPhi:
    //   std::cout << "Etas selected" << std::endl;
    //   AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
    //   trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  333.);
    //   AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
    //   trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -333.);
    //   cgTrackCuts->AddCut(trackCutsMother1);
    //   cgTrackCuts->AddCut(trackCutsMother2);
    // break;
    // case kRho:
    //   std::cout << "Etas selected" << std::endl;
    //   AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
    //   trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  113.);
    //   AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
    //   trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -113.);
    //   cgTrackCuts->AddCut(trackCutsMother1);
    //   cgTrackCuts->AddCut(trackCutsMother2);
    // break;
    // case kJPsi:
    //   std::cout << "Etas selected" << std::endl;
    //   AliDielectronVarCuts* trackCutsMother1 = new AliDielectronVarCuts("trackCutsMother1","trackCutsMother1");
    //   trackCutsMother1->AddCut(AliDielectronVarManager::kPdgCodeMother,  443.);
    //   AliDielectronVarCuts* trackCutsMother2 = new AliDielectronVarCuts("trackCutsMother2","trackCutsMother2");
    //   trackCutsMother2->AddCut(AliDielectronVarManager::kPdgCodeMother, -443.);
    //   cgTrackCuts->AddCut(trackCutsMother1);
    //   cgTrackCuts->AddCut(trackCutsMother2);
    // break;
//     case kAll:
//       std::cout << "All selected" << std::endl;
//       AliDielectronVarCuts* trackCutsMother = new AliDielectronVarCuts("trackCutsMother","trackCutsMother");
//       // trackCutsMother->AddCut(AliDielectronVarManager::kPdgCodeMother, -10000., 10000.);
//       cgTrackCuts->AddCut(trackCutsMother);
//     break;
//     // default: cout << "No Analysis Track Cut defined " << endl;
//   }
//
//   return cgTrackCuts;
// }



// #######################################################################################################################
AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (cutSet) {
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

    default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
}


// #######################################################################################################################
//Selection of relatively 'flat' centralities
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(AnalysisCut AnaCut) {
  AliDielectronVarCuts* centCuts = 0x0;
  switch (AnaCut.GetCentrality()) {
    case kPbPbCentral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbCentral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.,10.);
      break;
    case kPbPbSemiCentral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbSemiCentral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,10.,50.);
      break;
    case kPbPbSemiCentralRun1:
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
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.001,90.);
      break;
    default: cout << "No Centrality selected" << endl;
  }
  return centCuts;
}

// #######################################################################################################################
// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask, currently to 'kPbPb2011_TPCTOF_Semi1'
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
    case kStandard:
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

// #######################################################################################################################
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



// #######################################################################################################################
//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
AliAnalysisCuts* LMEECutLib::GetPairCutsAna(AnalysisCut AnaCut, Int_t togglePC)  {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pairCuts=0x0;
  switch (AnaCut.GetPairCutsAna()) {
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



// //Pair Cuts for PREFILTER step
// // cuts = REJECTION!!!
// AliAnalysisCuts* LMEECutLib::GetPairCutsPre(AnalysisCut AnaCut)  {
//   cout << " >>>>>>>>>>>>>>>>>>>>>> GetPairCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
//   AliAnalysisCuts* pairCuts=0x0;
//   switch (AnaCut.GetPairCutsPre()) {
//     case kInvM0to030MeV_OpAng0to060mrad:
//       AliDielectronVarCuts* pairCutsInvM =new AliDielectronVarCuts("pairCutsInvM","pairCutsInvM");
//       pairCutsInvM->AddCut(AliDielectronVarManager::kM, 0.0, 0.03); // in upgrade: 0.01
//       AliDielectronVarCuts* pairCutsOpAng =new AliDielectronVarCuts("pairCutsOpAng","pairCutsOpAng");
//       pairCutsOpAng->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0, 0.06); // in upgrade: 0.05
//
//       AliDielectronCutGroup* pairCutsCG =new AliDielectronCutGroup("pairCutsCG","pairCutsCG",AliDielectronCutGroup::kCompAND);
//       pairCutsCG->AddCut(pairCutsInvM);
//       pairCutsCG->AddCut(pairCutsOpAng);
//       //pairCutsCG->AddCut(pairCutsPhiv);
//       pairCuts = pairCutsCG;
//       break;
//     case kNoPairCutsPre:
//       //[...] // PhiV and InvMass
//     default: cout << "No Prefilter Pair Cuts defined " << endl;
//   }
//   return pairCuts;
// }

// //Relaxed PID cuts for additional rejectin step, do not use blindly
// AliAnalysisCuts* LMEECutLib::GetPIDCutsPre(AnalysisCut AnaCut) {
//   cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
//   AliAnalysisCuts* pidCuts=0x0;
//   switch (AnaCut.GetPIDPre()) {
//     case kStandardPre:
//
//       // eta range:
//       AliDielectronVarCuts *etaRangePre1 = new AliDielectronVarCuts("etaRangePre1","etaRangePre1");
//       etaRangePre1->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
//       // pt range:
//       AliDielectronVarCuts *ptRangePre1 = new AliDielectronVarCuts("ptRangePre1","ptRangePre1");
//       ptRangePre1->AddCut(AliDielectronVarManager::kPt, .2, 3.5); // 0.2 is realistic. turnon at ~180MeV
//
//       AliDielectronCutGroup* cgITSTPCTOFpre = new AliDielectronCutGroup("cgITSTPCTOFpre","cgITSTPCTOFpre",AliDielectronCutGroup::kCompAND);
//       AliDielectronPID *pidITSTPCTOFpre = new AliDielectronPID("pidITSTPCTOFpre","pidITSTPCTOFpre");
//       pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3., 0. ,100., kFALSE);
//       pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -3. , 3., 0. ,100., kTRUE);
//       // ITS will be used:
//       pidITSTPCTOFpre->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -4. , 2., 0. ,  2., kFALSE);
//       // TOF will be used if available, and with pt instead of p:
//       //  pidITSTPCTOFpre->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3., 0.4,100., kFALSE,
//       //                          AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
//       cgITSTPCTOFpre->AddCut(pidITSTPCTOFpre);
//       cgITSTPCTOFpre->AddCut(etaRangePre1);
//       cgITSTPCTOFpre->AddCut(ptRangePre1);
//       cgITSTPCTOFpre->AddCut(GetTrackSelectionPre(AnaCut));
//
//       AliDielectronCutGroup* cgInitialTrackFilter = new AliDielectronCutGroup("cgInitialTrackFilter","cgInitialTrackFilter",AliDielectronCutGroup::kCompOR);
//       // in case the prefilter cuts do not include all needed global tracks.
//       cgInitialTrackFilter->AddCut(GetPIDCutsAna(AnaCut));
//       cgInitialTrackFilter->AddCut(cgITSTPCTOFpre);
//
//       pidCuts = cgInitialTrackFilter;   // kCompOR works!!! <- checked with 'SetNoPairing()' and commented out 'GetPIDCutsAna(selectedPID)'
//       //cout << " ========== pidCuts prefilter: ========== " << endl;
//       //pidCuts->Print();
//       break;
//
//     case default:
//     default: cout << "No Prefilter PID Cut defined " << endl;
//   }
//   return pidCuts;
// }


// //Possibly different cut sets for Prefilter step
// //Not used at the moment
// AliAnalysisCuts* LMEECutLib::GetTrackSelectionPre(AnalysisCut AnaCut) {
//   cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionPre() >>>>>>>>>>>>>>>>>>>>>> " << endl;
//   AliDielectronCutGroup* trackCuts=0x0;
//   switch (AnaCut.GetTrackSelectionPre()) {
//     case kPrefilter_cut1:
//       AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
//       trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
//       trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
//       trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,     3.0, 100.0);
//       AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
//       trackCutsDiel->SetAODFilterBit(1); //does nothing for ESDs, ITSSA(???) // maybe use FilterBit(2) instead!
//       //        trackCutsDiel->SetRequireITSRefit(kTRUE); //function in AliDielectronTrackCuts
//       //        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst); //function in AliDielectronTrackCuts
//
//       cgTrackCutsPre = new AliDielectronCutGroup("cgTrackCutsPre","cgTrackCutsPre",AliDielectronCutGroup::kCompAND);
//       cgTrackCutsPre->AddCut(trackCutsDiel);
//       cgTrackCutsPre->AddCut(trackCutsAOD);
//       trackCuts = cgTrackCutsPre;
//       break;
//
//     default: cout << "No Prefilter Track Cut defined " << endl;
//   }
//   return trackCuts;
// }


// //Basic track rotator settings from J/Psi, more investigation needed
// AliDielectronTrackRotator* LMEECutLib::GetTrackRotator(Int_t cutSet) {
//   AliDielectronTrackRotator* trackRotator = 0x0;
//   switch (cutSet) {
//     default: cout << "No Rotator defined" << endl;
//       //default:
//       //  trackRotator = new AliDielectronTrackRotator();
//       //  trackRotator->SetIterations(20);
//       //  trackRotator->SetConeAnglePhi(TMath::Pi()/180*165);
//       //  trackRotator->SetStartAnglePhi(TMath::Pi());
//       //  break;
//   }
//   return trackRotator;
// }

//
AliAnalysisCuts* LMEECutLib::GetESDTrackCutsAna(AnalysisCut AnaCut) {
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetESDTrackCutsAna()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliESDtrackCuts* esdTrackCutsH = 0x0;
  switch (AnaCut.GetESDTrackSelection()) {
    default:
      esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
      esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
      esdTrackCutsH->SetDCAToVertex2D(kTRUE);
      break;
  }
  return esdTrackCutsH;
}

#endif
