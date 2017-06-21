#ifndef LMEECutLib_caklein
#define LMEECutLib_caklein

// #include "AliDielectron.h"

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
  void SetStandardCut(){
    SetPairCutsAna(LMEECutLib::kNoPairCutsAna);

    SetPreFilterType(LMEECutLib::kNoPreFilter);
    SetPIDPre(LMEECutLib::kStandardPre);
    SetTrackSelectionPre(LMEECutLib::kPrefilter_cut1);
    SetPairCutsPre(LMEECutLib::kNoPairCutsPre);

    SetMixing(LMEECutLib::kEventMixing_1);
    SetESDTrackSelection(LMEECutLib::kStandardESD);
  };

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
    // Clean Samples
    kPbPb2015_pidV0_electron_pt400_woPionRej,
    kPbPb2015_pidV0_electron_pt400_wITScut,
    kPbPb2015_pure_pion_pt400_wITScut,
    kPbPb2015_pure_kaon_pt400_wITScut,
    kPbPb2015_pure_proton_pt400_wITScut,
    kPbPb2015_pidV0_electron_pt400_wTPCcut,
    kPbPb2015_pure_pion_pt400_wTPCcut,
    kPbPb2015_pure_kaon_pt400_wTPCcut,
    kPbPb2015_pure_proton_pt400_wTPCcut,
    kPbPb2015_pure_electron_pt400,
    kPbPb2015_pure_electron_pt400_woTOF,
    kPbPb2015_pure_pion_pt400,
    kPbPb2015_pure_kaon_pt400,
    kPbPb2015_pure_proton_pt400,
    // no PID
    kITSTPCTOFif_trkSPDfirst_kINT7_pt400_woPID,
    // Analysis cuts
    kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton,
    kPbPb2015_Pt400_tightTOFreq,
    kPbPb2015_Pt400_symTOFreq,
    kPbPb2015_Pt400_tightTOFif,
    kPbPb2015_Pt400_looseTOFif,
    kPIDcut_0,
    kPIDcut_1,
    kPIDcut_2,
    kPIDcut_3,
    kPIDcut_4,
    kPIDcut_5,
    kPIDcut_6,
    kPIDcut_7,
    kPIDcut_8,
    kPIDcut_9,
    kPIDcut_10,
    kPIDcut_11,
    kPIDcut_12,
    kPIDcut_13,
    kPIDcut_14,
    kPIDcut_15,
    kPIDcut_16,
    kPIDcut_17,
    kPIDcut_18,
    kPIDcut_19,
    // Contamination Study
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noTPCcut,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noITScut,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noPionRej,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noTPCcut,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noITScut,
    kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noPionRej,
    // MISC
    knanoAODTOFeffCut,
    knanoAODTOFeffCut_wTOF,
    kPionTOFeff_wTOF,
    kPionTOFeff_woTOF,
    kPbPb2015_Pt100_ResolutionCuts
  };
  enum LMEEPIDPre{
    kStandardPre
  };
  // Possible Track Selections
  enum LMEETrackSelectionAna{
    kV0,
    kSPDfirst,
    kITSSA,
    kNone,
    kTRACKcut_0,
    kTRACKcut_1,
    kTRACKcut_2,
    kTRACKcut_3,
    kTRACKcut_4,
    kTRACKcut_5,
    kTRACKcut_6,
    kTRACKcut_7,
    kTRACKcut_8,
    kTRACKcut_9,
    kTRACKcut_10,
    kTRACKcut_11,
    kTRACKcut_12,
    kTRACKcut_13,
    kTRACKcut_14,
    kTRACKcut_15,
    kTRACKcut_16,
    kTRACKcut_17,
    kTRACKcut_18,
    kTRACKcut_19
  };
  enum LMEETrackSelectionPre{
    kPrefilter_cut1
  };
  enum LMEETrackCuts{
    kPbPb2015_V0_tight,
    kSPD_bit4,
    kITSSA_bit1,
    kNoTrackCuts
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
    kPbPbSemiCentral, //10%-50%
    kPbPbPeripheral,  //50%-90%
    kPbPb_00to80,     //0%-80%
  };
  enum LMEEEventMixing{
    kEventMixing_1  // 5 classes Z-Vertex, 7 classes centrality
  };
  enum LMEEESDTrackSelection{
    kStandardESD
  };
  enum LMEEEventCut{
    kStandard,
    kStandard_run2cuts
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
  AliDielectronCutGroup* LMEECutLib::SetKinematics(AliDielectronVarCuts *etaRange, AliDielectronVarCuts *ptRange, AliDielectronPID* PID, AnalysisCut AnaCut) {
    AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
    cgPIDCutsAna->AddCut(etaRange);
    cgPIDCutsAna->AddCut(ptRange);
    cgPIDCutsAna->AddCut(PID);
    cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
    return cgPIDCutsAna;
  }


  void SetEtaCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/

};

void LMEECutLib::SetEtaCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrection()\n";
  std::string file_name = "/home/cklein/LMeeAnaFW/005_TPC_Recalibration/summary/output.root";

  TFile* _file = TFile::Open(file_name.c_str());
  std::cout << _file << std::endl;
  if (_file == 0x0){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/data/output.root .");
    std::cout << "Copy TPC correction from Alien" << std::endl;
    _file = TFile::Open("output.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }
  if (runwise){
    TObjArray* arr_mean = dynamic_cast<TObjArray*>(_file->Get("mean_correction_arr"));
    TObjArray* arr_width =dynamic_cast<TObjArray*>( _file->Get("width_correction_arr"));
    std::cout << arr_mean << " " << arr_width << std::endl;

    die->SetWidthCorrArr(arr_width, kTRUE, corrXdim, corrYdim, corrZdim);
    die->SetCentroidCorrArr(arr_mean, kTRUE, corrXdim, corrYdim, corrZdim);
  }
  else{
    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));
    die->SetCentroidCorrFunction(mean, corrXdim, corrYdim, corrZdim);
    die->SetWidthCorrFunction(width, corrXdim, corrYdim, corrZdim);
  }

}


// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {
  AliDielectronEventCuts* eventCuts = 0x0;
  switch (cutSet) {
    case 0:
      std::cout << "Event Cuts 0" << std::endl;

      //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
      eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
      eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
      eventCuts->SetRequireVertex();
      eventCuts->SetMinVtxContributors(1);
      eventCuts->SetVertexZ(-10.,10.);
      break;
    case 1:
      std::cout << "Event Cuts 1" << std::endl; // Cuts really bad in Centrality distribution. Do not use at the moment 17.05.2017
      //Basic Event Cuts for pp and Pb-Pb, additional cuts may be in the AddTask
      eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
      eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
      eventCuts->SetRequireVertex();
      eventCuts->SetMinVtxContributors(1);
      eventCuts->SetVertexZ(-10.,10.);

      // Event cut copied from Pascal and adapted for 10-50% centrality
      TF1* fRefMultVZEROmultUp = new TF1("cutRefMultVZEROmult", "[0]+[1]*x+[2]*x*x", 0., 3.e+5);
      fRefMultVZEROmultUp->SetParameters(600., 0.35, 3e-6);

      TF1* fRefMultVZEROmultLow = new TF1("cutRefMultVZEROmult", "[0]+[1]*x+[2]*x*x", 0., 3.e+5);
      fRefMultVZEROmultLow->SetParameters(-700., 0.3, 3e-6);

    	eventCuts->SetMinCorrCutFunction(fRefMultVZEROmultLow, AliDielectronVarManager::kMultV0, AliDielectronVarManager::kRefMult);
    	eventCuts->SetMaxCorrCutFunction(fRefMultVZEROmultUp, AliDielectronVarManager::kMultV0, AliDielectronVarManager::kRefMult);
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
    case kPbPbSemiCentral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbSemiCentral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,10.,50.);
      break;
    case kPbPbPeripheral:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPbPeripheral");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,50.,80.);
      break;
    case kPbPb_00to80:
      centCuts = new AliDielectronVarCuts("centCuts","CentralityPbPb_00to80");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.01,80.);
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
      mixingHandler->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 6, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(40);
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

AliAnalysisCuts* LMEECutLib::GetPIDCutsAna(AnalysisCut AnaCut) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;

  //-----------------------------------------------
  // Define different PID Cuts, that are used later
  //-----------------------------------------------


  // eta range:
  AliDielectronVarCuts *etaRange080 = new AliDielectronVarCuts("etaRange080","etaRange080");
  etaRange080->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
  // pt range:
  AliDielectronVarCuts *ptRange400to8000 = new AliDielectronVarCuts("ptRange400to8000","ptRange400to8000");
  ptRange400to8000->AddCut(AliDielectronVarManager::kPt, 0.4, 8.0);
  AliDielectronVarCuts *ptRange100to8000 = new AliDielectronVarCuts("ptRange100to8000","ptRange100to8000");
  ptRange100to8000->AddCut(AliDielectronVarManager::kPt, 0.1, 8.0);


  // Normal Analysis Cuts
  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif","pid_TPCele_AsymITS_tightTOFif");
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_looseTOFif = new AliDielectronPID("pid_TPCele_AsymITS_looseTOFif","pid_TPCele_AsymITS_looseTOFif");
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFreq = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFreq","pid_TPCele_AsymITS_tightTOFreq");
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *pid_TPCele_AsymITS_symTOFreq = new AliDielectronPID("pid_TPCele_AsymITS_symTOFreq","pid_TPCele_AsymITS_symTOFreq");
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 4. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

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

  AliDielectronPID *PIDcut_0 = new AliDielectronPID("PIDcut_0","PIDcut_0");
  PIDcut_0->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.5 , 0. ,100., kFALSE);
  PIDcut_0->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_0->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_0->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0. , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_1 = new AliDielectronPID("PIDcut_1","PIDcut_1");
  PIDcut_1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_1->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_1->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_2 = new AliDielectronPID("PIDcut_2","PIDcut_2");
  PIDcut_2->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_2->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_2->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_2->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_3 = new AliDielectronPID("PIDcut_3","PIDcut_3");
  PIDcut_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_3->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_4 = new AliDielectronPID("PIDcut_4","PIDcut_4");
  PIDcut_4->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_4->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_4->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_4->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_5 = new AliDielectronPID("PIDcut_5","PIDcut_5");
  PIDcut_5->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_5->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_5->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_5->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_6 = new AliDielectronPID("PIDcut_6","PIDcut_6");
  PIDcut_6->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_6->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_6->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  PIDcut_6->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_7 = new AliDielectronPID("PIDcut_7","PIDcut_7");
  PIDcut_7->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_7->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_7->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.5 , 0. ,100., kFALSE);
  PIDcut_7->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_8 = new AliDielectronPID("PIDcut_8","PIDcut_8");
  PIDcut_8->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_8->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_8->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_8->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_9 = new AliDielectronPID("PIDcut_9","PIDcut_9");
  PIDcut_9->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_9->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_9->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_9->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_10 = new AliDielectronPID("PIDcut_10","PIDcut_10");
  PIDcut_10->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_10->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_10->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.5 , 0. ,100., kFALSE);
  PIDcut_10->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_11 = new AliDielectronPID("PIDcut_11","PIDcut_11");
  PIDcut_11->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_11->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_11->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_11->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.5. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_12 = new AliDielectronPID("PIDcut_12","PIDcut_12");
  PIDcut_12->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_12->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_12->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_12->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_13 = new AliDielectronPID("PIDcut_13","PIDcut_13");
  PIDcut_13->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.5 , 0. ,100., kFALSE);
  PIDcut_13->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_13->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_13->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_14 = new AliDielectronPID("PIDcut_14","PIDcut_14");
  PIDcut_14->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_14->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_14->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 0.5 , 0. ,100., kFALSE);
  PIDcut_14->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_15 = new AliDielectronPID("PIDcut_15","PIDcut_15");
  PIDcut_15->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_15->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_15->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_15->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_16 = new AliDielectronPID("PIDcut_16","PIDcut_16");
  PIDcut_16->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  PIDcut_16->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_16->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_16->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_17 = new AliDielectronPID("PIDcut_17","PIDcut_17");
  PIDcut_17->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  PIDcut_17->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_17->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  PIDcut_17->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_18 = new AliDielectronPID("PIDcut_18","PIDcut_18");
  PIDcut_18->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_18->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 5.0 , 0. ,100., kTRUE);
  PIDcut_18->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_18->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_19 = new AliDielectronPID("PIDcut_19","PIDcut_19");
  PIDcut_19->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_19->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.0 , 0. ,100., kTRUE);
  PIDcut_19->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_19->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5. , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);



  // ########### Cuts for PID Contamination Analysis
  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif_noTPCcut = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif_noTPCcut","pid_TPCele_AsymITS_tightTOFif_noTPCcut");
  pid_TPCele_AsymITS_tightTOFif_noTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFif_noTPCcut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif_noTPCcut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif_noITScut = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif_noITScut","pid_TPCele_AsymITS_tightTOFif_noITScut");
  pid_TPCele_AsymITS_tightTOFif_noITScut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFif_noITScut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif_noITScut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif_noPionRej = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif_noPionRej","pid_TPCele_AsymITS_tightTOFif_noPionRej");
  pid_TPCele_AsymITS_tightTOFif_noPionRej->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif_noPionRej->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif_noPionRej->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFreq_noTPCcut = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFreq_noTPCcut","pid_TPCele_AsymITS_tightTOFreq_noTPCcut");
  pid_TPCele_AsymITS_tightTOFreq_noTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFreq_noTPCcut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq_noTPCcut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFreq_noITScut = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFreq_noITScut","pid_TPCele_AsymITS_tightTOFreq_noITScut");
  pid_TPCele_AsymITS_tightTOFreq_noITScut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFreq_noITScut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2. , 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq_noITScut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFreq_noPionRej = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFreq_noPionRej","pid_TPCele_AsymITS_tightTOFreq_noPionRej");
  pid_TPCele_AsymITS_tightTOFreq_noPionRej->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq_noPionRej->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq_noPionRej->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  // ################## PID for Clean Samples
  AliDielectronPID *pid_electron_woPionRej = new AliDielectronPID("pid_electron_woPionRej","pid_electron_woPionRej");
  pid_electron_woPionRej->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron_woPionRej->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1., 0. ,100., kFALSE); // ITS Cut is asymmetric to compare to TOFif Cuts for Contamination studies
  pid_electron_woPionRej->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);

  AliDielectronPID *pid_electron_wITScut = new AliDielectronPID("pid_electron_wITScut","pid_electron_wITScut");
  pid_electron_wITScut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron_wITScut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1., 0. ,100., kFALSE); // ITS Cut is asymmetric to compare to TOFif Cuts for Contamination studies
  pid_electron_wITScut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_pion_wITScut = new AliDielectronPID("pid_pion_wITScut","pid_pion_wITScut");
  pid_pion_wITScut->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_pion_wITScut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1., 0. ,100., kFALSE); // ITS Cut is asymmetric to compare to TOFif Cuts for Contamination studies

  AliDielectronPID *pid_kaon_wITScut = new AliDielectronPID("pid_kaon_wITScut","pid_kaon_wITScut");
  pid_kaon_wITScut->AddCut(AliDielectronPID::kTOF,AliPID::kKaon, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_kaon_wITScut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1., 0. ,100., kFALSE); // ITS Cut is asymmetric to compare to TOFif Cuts for Contamination studies
  pid_kaon_wITScut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_proton_wITScut = new AliDielectronPID("pid_proton_wITScut","pid_proton_wITScut");
  pid_proton_wITScut->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_proton_wITScut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1., 0. ,100., kFALSE); // ITS Cut is asymmetric to compare to TOFif Cuts for Contamination studies
  pid_proton_wITScut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);


  AliDielectronPID *pid_electron_wTPCcut = new AliDielectronPID("pid_electron_wTPCcut","pid_electron_wTPCcut");
  pid_electron_wTPCcut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_electron_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_pion_wTPCcut = new AliDielectronPID("pid_pion_wTPCcut","pid_pion_wTPCcut");
  pid_pion_wTPCcut->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_pion_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);

  AliDielectronPID *pid_kaon_wTPCcut = new AliDielectronPID("pid_kaon_wTPCcut","pid_kaon_wTPCcut");
  pid_kaon_wTPCcut->AddCut(AliDielectronPID::kTOF,AliPID::kKaon, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_kaon_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_kaon_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_proton_wTPCcut = new AliDielectronPID("pid_proton_wTPCcut","pid_proton_wTPCcut");
  pid_proton_wTPCcut->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_proton_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_proton_wTPCcut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_electron = new AliDielectronPID("pid_electron","pid_electron");
  pid_electron->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  // pid_electron->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 4. , 0. ,100., kTRUE);

  AliDielectronPID *pid_electron_woTOF = new AliDielectronPID("pid_electron_woTOF","pid_electron_woTOF");
  pid_electron_woTOF->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron_woTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  // pid_electron_woTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 4. , 0. ,100., kTRUE);

  AliDielectronPID *pid_pion = new AliDielectronPID("pid_pion","pid_pion");
  pid_pion->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_pion->AddCut(AliDielectronPID::kITS,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);

  AliDielectronPID *pid_kaon = new AliDielectronPID("pid_kaon","pid_kaon");
  pid_kaon->AddCut(AliDielectronPID::kTOF,AliPID::kKaon, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_kaon->AddCut(AliDielectronPID::kITS,AliPID::kKaon, -3.0, 3. , 0. ,100., kFALSE);
  pid_kaon->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_proton = new AliDielectronPID("pid_proton","pid_proton");
  pid_proton->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_proton->AddCut(AliDielectronPID::kITS,AliPID::kProton, -3.0, 3. , 0. ,100., kFALSE);
  pid_proton->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  // Cuts to analyse TOF Efficiency
  AliDielectronPID *nanoAODTOFeffCut = new AliDielectronPID("nanoAODTOFeffCut","nanoAODTOFeffCut");
  nanoAODTOFeffCut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  nanoAODTOFeffCut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  AliDielectronPID *nanoAODTOFeffCut_wTOF = new AliDielectronPID("nanoAODTOFeffCut_wTOF","nanoAODTOFeffCut_wTOF");
  nanoAODTOFeffCut_wTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  nanoAODTOFeffCut_wTOF->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  nanoAODTOFeffCut_wTOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -99.0, 99. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  AliDielectronPID *pid_pion_wTOF = new AliDielectronPID("pid_pion_wTOF","pid_pion_wTOF");
  pid_pion_wTOF->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -99., 99., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_pion_wTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);
  AliDielectronPID *pid_pion_woTOF = new AliDielectronPID("pid_pion_woTOF","pid_pion_woTOF");
  pid_pion_woTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);



  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------
  switch (AnaCut.GetPIDAna()) {
    case kITSTPCTOFif_trkSPDfirst_kINT7_pt400_woPID:
    // NO PID
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to8000);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;

    case kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton:
    // Cut out pion/kaon/proton band LHC15o but refilled when particle in TOF electron band
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PID_cutoff_pion_kaon_proton_cg, AnaCut);;
      break;

    case kPbPb2015_Pt400_tightTOFif:
    // ITS & TOFifavailable & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFif, AnaCut);
      break;

    case kPbPb2015_Pt400_looseTOFif:
    // ITS & TOFifavailable & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_looseTOFif, AnaCut);;
      break;
    case kPbPb2015_Pt400_tightTOFreq:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFreq, AnaCut);;
      break;

    case kPbPb2015_Pt400_symTOFreq:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_symTOFreq, AnaCut);;
      break;
    case kPIDcut_0:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_0, AnaCut);;
      break;
    case kPIDcut_1:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_1, AnaCut);;
      break;
    case kPIDcut_2:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_2, AnaCut);;
      break;
    case kPIDcut_3:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_3, AnaCut);;
      break;
    case kPIDcut_4:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_4, AnaCut);;
      break;
    case kPIDcut_5:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_5, AnaCut);;
      break;
    case kPIDcut_6:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_6, AnaCut);;
      break;
    case kPIDcut_7:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_7, AnaCut);;
      break;
    case kPIDcut_8:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_8, AnaCut);;
      break;
    case kPIDcut_9:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_9, AnaCut);;
      break;
    case kPIDcut_10:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_10, AnaCut);;
      break;
    case kPIDcut_11:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_11, AnaCut);;
      break;
    case kPIDcut_12:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_12, AnaCut);;
      break;
    case kPIDcut_13:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_13, AnaCut);;
      break;
    case kPIDcut_14:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_14, AnaCut);;
      break;
    case kPIDcut_15:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_15, AnaCut);;
      break;
    case kPIDcut_16:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_16, AnaCut);;
      break;
    case kPIDcut_17:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_17, AnaCut);;
      break;
    case kPIDcut_18:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_18, AnaCut);;
      break;
    case kPIDcut_19:
    // ITS & TOFrequired & TPC sigma cut for LHC15o
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_19, AnaCut);;
      break;

    // ################## CONTAMINATION STUDY STUFF

    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noTPCcut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFif_noTPCcut, AnaCut);;
      break;

    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noITScut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFif_noITScut, AnaCut);;
      break;

      case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noPionRej:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFif_noPionRej, AnaCut);;
      break;
    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noTPCcut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFreq_noTPCcut, AnaCut);;
      break;

    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noITScut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFreq_noITScut, AnaCut);;
      break;

    case kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noPionRej:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_TPCele_AsymITS_tightTOFreq_noPionRej, AnaCut);;
      break;
      // ################# CLEAN SAMPLES ##################
    case kPbPb2015_pidV0_electron_pt400_woPionRej:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_electron_woPionRej, AnaCut);;
      break;

    case kPbPb2015_pure_electron_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_electron, AnaCut);;
      break;
    case kPbPb2015_pure_electron_pt400_woTOF:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_electron_woTOF, AnaCut);;
      break;
    case kPbPb2015_pure_kaon_pt400_wITScut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_kaon_wITScut, AnaCut);;
      break;
    case kPbPb2015_pure_proton_pt400_wITScut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_proton_wITScut, AnaCut);;
      break;
    case kPbPb2015_pidV0_electron_pt400_wTPCcut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_electron_wTPCcut, AnaCut);;
      break;
    case kPbPb2015_pure_pion_pt400_wTPCcut:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange080);
      cgPIDCutsAna->AddCut(ptRange400to8000);
      // cgPIDCutsAna->AddCut(pid_pion_wTPCcut);
      cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
      pidCuts = cgPIDCutsAna;
      break;
    case kPbPb2015_pure_kaon_pt400_wTPCcut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_kaon_wTPCcut, AnaCut);;
      break;
    case kPbPb2015_pure_proton_pt400_wTPCcut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_proton_wTPCcut, AnaCut);;
      break;
    case kPbPb2015_pure_pion_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_pion, AnaCut);;
      break;
    case kPbPb2015_pure_kaon_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_kaon, AnaCut);;
      break;
    case kPbPb2015_pure_proton_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_proton, AnaCut);;
      break;

    // MISC Cut settings
    case knanoAODTOFeffCut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, nanoAODTOFeffCut, AnaCut);;
      break;
    case knanoAODTOFeffCut_wTOF:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, nanoAODTOFeffCut_wTOF, AnaCut);;
      break;
    case kPionTOFeff_wTOF:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_pion_wTOF, AnaCut);;
      break;
    case kPionTOFeff_woTOF:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, pid_pion_woTOF, AnaCut);;
      break;
    case kPbPb2015_Pt100_ResolutionCuts:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange120);
      cgPIDCutsAna->AddCut(ptRange100to8000);
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
      std::cout << "SPDfirst" << std::endl;
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);

      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0); // means at least 2 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1); // means 0 and 1 shared Cluster

      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014

      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      // trackCuts = GetTrackCuts(kSPD_bit4);
      // break;
    case kITSSA:
      trackCuts = GetTrackCuts(kITSSA_bit1);
      break;
      case kTRACKcut_0:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_1:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_2:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_3:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_4:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   0.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_5:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_6:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   2.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_7:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   2.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_8:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_9:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_10:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_11:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_12:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   2.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   0.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   2.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_13:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_14:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_15:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_16:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_17:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   2.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_18:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   0.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_19:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   1.1);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   2.5);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;

    case kNone:
      trackCuts = GetTrackCuts(kNoTrackCuts);
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

      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0); // means at least 2 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1); // means 0 and 1 shared Cluster

      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    130.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014

      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
    case kNoTrackCuts:
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); //GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

      cgTrackCutsAnaSDDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSDDfirst","cgTrackCutsAnaSDDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSDDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSDDfirst->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsAnaSDDfirst;
      break;

    default: cout << "No Analysis Track Cut defined " << endl;
  }
  return trackCuts;
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
      ptRangePre1->AddCut(AliDielectronVarManager::kPt, 0.2, 8.0); // 0.2 is realistic. turnon at ~180MeV

      AliDielectronCutGroup* cgITSTPCTOFpre = new AliDielectronCutGroup("cgITSTPCTOFpre","cgITSTPCTOFpre",AliDielectronCutGroup::kCompAND);
      AliDielectronPID *pidITSTPCTOFpre = new AliDielectronPID("pidITSTPCTOFpre","pidITSTPCTOFpre");
      pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3., 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      pidITSTPCTOFpre->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 3., 0. ,100., kTRUE,  AliDielectronPID::kIfAvailable);
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
      trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.9,   0.9);
      trackCutsAOD->AddCut(AliDielectronVarManager::kPt, -0.2,   8.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

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
  cout << " >>>>>>>>>>>>>>>>>>>>>>  GetESDTrackCutsAna()  >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliESDtrackCuts* esdTrackCutsH = 0x0;
  switch (AnaCut.GetESDTrackSelection()) {
    default:
      // standard cuts with very loose DCA: Bit4 (Int: 16), AOD095&115
      esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
      esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
      esdTrackCutsH->SetDCAToVertex2D(kTRUE);

      break;
      //default: cout << "No ESD Track Cut defined " << endl;
  }
  return esdTrackCutsH;
}
