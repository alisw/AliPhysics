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
    kPbPb2015_pure_electron_pt200,
    kPbPb2015_pure_pion_pt200,

    kPbPb2015_Pt200_cut5_woTPCelecut,
    kPbPb2015_Pt200_cut5_woPionRej,
    kPbPb2015_pure_electron_pt200_woTPCelecut,
    kPbPb2015_pure_kaon_pt200_woTPCelecut,
    kPbPb2015_pure_proton_pt200_woTPCelecut,
    kPbPb2015_pure_pion_pt200_wTPCelecut,

    // Analysis cuts
    kPbPb2015_Pt200_PID_cutoff_pion_kaon_proton,
    kPbPb2015_Pt200_noPID,
    kPDGelectron,
    kPDGelectronMotherPhoton,
    kPIDcut_0_onlyLooseTPC,
    kPID_Jeromian_01,
    kPIDcut_0,
    kPIDcut_5_pt400,
    kPIDcut_11_pt400,
    kPIDcut_23_pt400,
    kPID_Jeromian_00,
    kPID_Jeromian_00_TPConly,
    kPID_Jeromian_00_TOFonly,

    kPIDcut_1_pt200,
    kPIDcut_2_pt200,
    kPIDcut_3_pt200,
    kPIDcut_4_pt200,
    kPIDcut_5_pt200,
    kPIDcut_6_pt200,
    kPIDcut_7_pt200,
    kPIDcut_8_pt200,
    kPIDcut_9_pt200,
    kPIDcut_10_pt200,
    kPIDcut_11_pt200,
    kPIDcut_12_pt200,
    kPIDcut_13_pt200,
    kPIDcut_14_pt200,
    kPIDcut_15_pt200,
    kPIDcut_16_pt200,
    kPIDcut_17_pt200,
    kPIDcut_18_pt200,
    kPIDcut_19_pt200,
    kPIDcut_20_pt200,
    kPIDcut_21_pt200,
    kPIDcut_22_pt200,
    kPIDcut_23_pt200,
    kPIDcut_24_pt200,
    kPIDcut_25_pt200,
    kPIDcut_26_pt200,
    kPIDcut_27_pt200,
    kPIDcut_28_pt200,
    kPIDcut_29_pt200,
    kPIDcut_30_pt200,
    kPIDcut_23_pt200_TRDincl,
    kPIDcut_5_pt200_looserPionRejection,
    // Contamination Study
    // MISC
    kPbPb2015_pure_electron_pt200_woTOF,
    kPbPb2015_pure_electron_pt200_wTOFhit,
    kPbPb2015_pure_electron_pt200_wTOFif,
    kPbPb2015_pure_electron_pt200_wTOFreq,
    kPbPb2015_PDG_pure_electron_pt200_woTOF,
    kPbPb2015_PDG_pure_electron_pt200_wTOFhit,
    kPbPb2015_PDG_pure_electron_pt200_wTOFif,
    kPbPb2015_PDG_pure_electron_pt200_wTOFreq,

    knanoAODTOFeffCut,
    kPbPb2015_Pt100_ResolutionCuts
  };
  enum LMEEPIDPre{
    kStandardPre
  };
  // Possible Track Selections
  enum LMEETrackSelectionAna{
    kV0,
    kNoTrackCuts,
    kResolutionTrackCuts,
    kZ0cuts,
    kSPDfirst,
    kSPDfirst_0SITSCl,
    kSPDfirst_ITSMAP,
    kSPDfirst_InverseSharedCluster,
    kITSdEdxConvRej,
    kSPDfirst_PDGCodePion,
    kSPDfirst_PDGCodeElectron,
    kSPDfirst_noConversion,
    kSPDfirst_Charm,
    kITSSA,
    kNone,
    kTRACKcut_0,
    kTRACKcut_5_without_CrossedOverFindableCut,
    kTRACKcut_5_woSharedCluster,
    kTRACKcut_5_0SharedCluster,
    kTRACKcut_5_woSharedCluster_noSPDfirst,
    kTRACKcut_5_0SharedCluster_noSPDfirst,
    kTRACKcut_5_noSPDfirst,
    kTRACKcut_5_noConvRejection,
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
    kTRACKcut_19,
    kTRACKcut_20,
    kTRACKcut_21,
    kTRACKcut_22,
    kTRACKcut_23,
    kTRACKcut_24,
    kTRACKcut_25,
    kTRACKcut_26,
    kTRACKcut_27,
    kTRACKcut_28,
    kTRACKcut_29,
    kTRACKcut_30,
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
    kPbPb0090,     // 0%-90%
    kPbPb0010,     // 0%-10%
    kPbPb1050,     // 10%-50%

  };
  enum LMEEEventMixing{
    kEventMixing_1,  // 5 classes Z-Vertex, 7 classes centrality
    kEventMixing_2,  // 5 classes Z-Vertex, 7 classes centrality
    kEventMixing_3,  // 5 classes Z-Vertex, 7 classes centrality
    kEventMixing_4,  // 5 classes Z-Vertex, 7 classes centrality
    kEventMixing_5  // 5 classes Z-Vertex, 7 classes centrality
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
  AliDielectronCutGroup* LMEECutLib::SetKinematics(AliDielectronVarCuts *etaRange, AliDielectronVarCuts *ptRange, AliAnalysisCuts* PID, AnalysisCut AnaCut) {
    AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
    cgPIDCutsAna->AddCut(etaRange);
    cgPIDCutsAna->AddCut(ptRange);
    cgPIDCutsAna->AddCut(PID);
    cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
    return cgPIDCutsAna;
  }


  void SetEtaCorrectionTPC(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void SetEtaCorrectionITS(AliDielectron *die, Int_t corrZdim, Int_t corrYdim, Bool_t data); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void SetEtaCorrectionTOF(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/

};

void LMEECutLib::SetEtaCorrectionTPC(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  std::string file_name = "/home/cklein/LMeeAnaFW/005_TPC_Recalibration/summary/output.root";

  TFile* _file = TFile::Open(file_name.c_str());
  std::cout << _file << std::endl;
  if (_file == 0x0){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/data/recalibration/recalib_data_tpc_nsigmaele.root .");
    std::cout << "Copy TPC correction from Alien" << std::endl;
    _file = TFile::Open("recalib_data_tpc_nsigmaele.root");
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
void LMEECutLib::SetEtaCorrectionITS(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionITS()\n";
  std::string file_name = "/home/cklein/LMeeAnaFW/005_TPC_Recalibration/summary/output_ITS.root";

  TFile* _file = TFile::Open(file_name.c_str());
  std::cout << _file << std::endl;
  if (_file == 0x0){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/data/recalibration/recalib_data_its_nsigmaele.root .");
    std::cout << "Copy ITS correction from Alien" << std::endl;
    _file = TFile::Open("recalib_data_its_nsigmaele.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }

  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));
  die->SetCentroidCorrFunctionITS(mean, corrXdim, corrYdim, corrZdim);
  die->SetWidthCorrFunctionITS(width, corrXdim, corrYdim, corrZdim);

}
void LMEECutLib::SetEtaCorrectionTOF(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionTOF()\n";
  std::string file_name = "/home/cklein/LMeeAnaFW/005_TPC_Recalibration/summary/output_TOF.root";

  TFile* _file = TFile::Open(file_name.c_str());
  std::cout << _file << std::endl;
  if (_file == 0x0){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/data/recalibration/recalib_data_tof_nsigmaele.root .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    _file = TFile::Open("recalib_data_tof_nsigmaele.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }

  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));
  die->SetCentroidCorrFunctionTOF(mean, corrXdim, corrYdim, corrZdim);
  die->SetWidthCorrFunctionTOF(width, corrXdim, corrYdim, corrZdim);

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
    case kPbPb0090:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0090");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 0., 90.);
      break;
    case kPbPb0010:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0010");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 0., 10.);
      break;
    case kPbPb1050:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb1050");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 10., 50.);
      break;
    case kPbPb0080:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0080");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 0., 80.);
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
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kEventMixing_2:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,0,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,20,30,50,80");
      mixingHandler->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 6, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kEventMixing_3:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,20,30,50,80");
      mixingHandler->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 1, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kEventMixing_4:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0,10,20,30,50,80");
      mixingHandler->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 6, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kEventMixing_5:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,20,30,50,80");
      mixingHandler->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 12, TMath::Pi()/-2., TMath::Pi()/2.);
      mixingHandler->SetDepth(60);
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
  AliDielectronVarCuts *etaRange120 = new AliDielectronVarCuts("etaRange120","etaRange120");
  etaRange120->AddCut(AliDielectronVarManager::kEta, -1.20, 1.20);
  // pt range:
  AliDielectronVarCuts *ptRange400to8000 = new AliDielectronVarCuts("ptRange400to8000","ptRange400to8000");
  ptRange400to8000->AddCut(AliDielectronVarManager::kPt, 0.4, 8.0);
  AliDielectronVarCuts *ptRange20000to100000 = new AliDielectronVarCuts("ptRange20000to100000","ptRange20000to100000");
  ptRange20000to100000->AddCut(AliDielectronVarManager::kPt, 4., 100.0);
  AliDielectronVarCuts *ptRange100to8000 = new AliDielectronVarCuts("ptRange100to8000","ptRange100to8000");
  ptRange100to8000->AddCut(AliDielectronVarManager::kPt, 0.1, 8.0);
  AliDielectronVarCuts *ptRange200to8000 = new AliDielectronVarCuts("ptRange200to8000","ptRange200to8000");
  AliDielectronVarCuts *ptRange100toINF = new AliDielectronVarCuts("ptRange100toINF","ptRange100toINF");
  ptRange100toINF->AddCut(AliDielectronVarManager::kPt, 0.1, 100.0);
  AliDielectronVarCuts *ptRange200to8000 = new AliDielectronVarCuts("ptRange200to8000","ptRange200to8000");
  ptRange200to8000->AddCut(AliDielectronVarManager::kPt, 0.2, 8.0);

  // PDG Code PID
  AliDielectronVarCuts *PDGelectron = new AliDielectronVarCuts("PDGelectron","PDGelectron");
  PDGelectron->AddCut(AliDielectronVarManager::kPdgCode, 11.);
  AliDielectronVarCuts *PDGpositron = new AliDielectronVarCuts("PDGelectron","PDGelectron");
  PDGpositron->AddCut(AliDielectronVarManager::kPdgCode, -11.);
  AliDielectronCutGroup* PDGlepton = new AliDielectronCutGroup("PDGlepton","PDGlepton",AliDielectronCutGroup::kCompOR);
  PDGlepton->AddCut(PDGelectron);
  PDGlepton->AddCut(PDGpositron);

  // PDG Code PID
  AliDielectronVarCuts *MotherIsPhoton = new AliDielectronVarCuts("MotherIsPhoton","MotherIsPhoton");
  MotherIsPhoton->AddCut(AliDielectronVarManager::kPdgCodeMother, 22.);

  AliDielectronCutGroup* PDGleptonMotherPhoton = new AliDielectronCutGroup("PDGleptonMotherPhoton","PDGleptonMotherPhoton",AliDielectronCutGroup::kCompOR);
  PDGleptonMotherPhoton->AddCut(PDGelectron);
  PDGleptonMotherPhoton->AddCut(PDGpositron);
  PDGleptonMotherPhoton->AddCut(MotherIsPhoton);

  // Normal Analysis Cuts
  AliDielectronPID *pid_TPCele_AsymITS_tightTOFif = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFif","pid_TPCele_AsymITS_tightTOFif");
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_looseTOFif = new AliDielectronPID("pid_TPCele_AsymITS_looseTOFif","pid_TPCele_AsymITS_looseTOFif");
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 4. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_looseTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_TPCele_AsymITS_noTOF = new AliDielectronPID("pid_TPCele_AsymITS_noTOF","pid_TPCele_AsymITS_noTOF");
  pid_TPCele_AsymITS_noTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_noTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 5. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_noTOF->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);

  AliDielectronPID *pid_TPCele_AsymITS_tightTOFreq = new AliDielectronPID("pid_TPCele_AsymITS_tightTOFreq","pid_TPCele_AsymITS_tightTOFreq");
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 4. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 1. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_tightTOFreq->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *pid_TPCele_AsymITS_symTOFreq = new AliDielectronPID("pid_TPCele_AsymITS_symTOFreq","pid_TPCele_AsymITS_symTOFreq");
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99., 4. , 0. ,100., kTRUE);
  pid_TPCele_AsymITS_symTOFreq->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PID_cutoff_pion_kaon_proton = new AliDielectronPID("PID_cutoff_pion_kaon_proton","PID_cutoff_pion_kaon_proton");
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -2.5, 2.5 , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -2.5, 2.5 , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  // PID_cutoff_pion_kaon_proton->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PID_cutoff_pion_kaon_proton_EleIncl = new AliDielectronPID("PID_cutoff_pion_kaon_proton_EleIncl","PID_cutoff_pion_kaon_proton_EleIncl");
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE);
  PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  // PID_cutoff_pion_kaon_proton_EleIncl->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronCutGroup* PID_cutoff_pion_kaon_proton_cg = new AliDielectronCutGroup("PID_cutoff_pion_kaon_proton_cg","PID_cutoff_pion_kaon_proton_cg",AliDielectronCutGroup::kCompOR);
  PID_cutoff_pion_kaon_proton_cg->AddCut(PID_cutoff_pion_kaon_proton);
  PID_cutoff_pion_kaon_proton_cg->AddCut(PID_cutoff_pion_kaon_proton_EleIncl);


  AliDielectronPID *Jeromian_01_hadron_cut = new AliDielectronPID("Jeromian_01_hadron_cut","Jeromian_01_hadron_cut");
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -2.5, 2.5 , 0. ,100., kTRUE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -2.5, 2.5 , 0. ,100., kTRUE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *Jeromian_01_ele_incl = new AliDielectronPID("Jeromian_01_ele_incl","Jeromian_01_ele_incl");
  Jeromian_01_ele_incl->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  Jeromian_01_ele_incl->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE);
  Jeromian_01_ele_incl->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronCutGroup* Jeromian_01_sum = new AliDielectronCutGroup("Jeromian_01_sum","Jeromian_01_sum",AliDielectronCutGroup::kCompOR);
  Jeromian_01_sum->AddCut(Jeromian_01_hadron_cut);
  // Jeromian_01_sum->AddCut(Jeromian_01_ele_incl);



  AliDielectronPID *PIDcut_1 = new AliDielectronPID("PIDcut_1","PIDcut_1");
  PIDcut_1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_1->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_1->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_2 = new AliDielectronPID("PIDcut_2","PIDcut_2");
  PIDcut_2->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.5 , 0. ,100., kFALSE);
  PIDcut_2->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_2->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_2->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_3 = new AliDielectronPID("PIDcut_3","PIDcut_3");
  PIDcut_3->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_3->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_3->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  PIDcut_3->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_4 = new AliDielectronPID("PIDcut_4","PIDcut_4");
  PIDcut_4->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  PIDcut_4->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_4->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  PIDcut_4->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_5 = new AliDielectronPID("PIDcut_5","PIDcut_5");
  PIDcut_5->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_5->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_5->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_5->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_6 = new AliDielectronPID("PIDcut_6","PIDcut_6");
  PIDcut_6->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_6->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_6->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_6->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_7 = new AliDielectronPID("PIDcut_7","PIDcut_7");
  PIDcut_7->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_7->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_7->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_7->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_8 = new AliDielectronPID("PIDcut_8","PIDcut_8");
  PIDcut_8->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_8->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_8->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_8->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_9 = new AliDielectronPID("PIDcut_9","PIDcut_9");
  PIDcut_9->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_9->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_9->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_9->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_10 = new AliDielectronPID("PIDcut_10","PIDcut_10");
  PIDcut_10->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_10->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_10->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  PIDcut_10->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_11 = new AliDielectronPID("PIDcut_11","PIDcut_11");
  PIDcut_11->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_11->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_11->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_11->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 2.5 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_12 = new AliDielectronPID("PIDcut_12","PIDcut_12");
  PIDcut_12->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  PIDcut_12->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_12->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_12->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 2.5 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_13 = new AliDielectronPID("PIDcut_13","PIDcut_13");
  PIDcut_13->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_13->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_13->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_13->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_14 = new AliDielectronPID("PIDcut_14","PIDcut_14");
  PIDcut_14->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_14->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_14->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 0.5 , 0. ,100., kFALSE);
  PIDcut_14->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_15 = new AliDielectronPID("PIDcut_15","PIDcut_15");
  PIDcut_15->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_15->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_15->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  PIDcut_15->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 2.5 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_16 = new AliDielectronPID("PIDcut_16","PIDcut_16");
  PIDcut_16->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.5 , 0. ,100., kFALSE);
  PIDcut_16->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_16->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 0.5 , 0. ,100., kFALSE);
  PIDcut_16->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 2.5 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_17 = new AliDielectronPID("PIDcut_17","PIDcut_17");
  PIDcut_17->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_17->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_17->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_17->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_18 = new AliDielectronPID("PIDcut_18","PIDcut_18");
  PIDcut_18->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  PIDcut_18->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 4.5 , 0. ,100., kTRUE);
  PIDcut_18->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_18->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 2.5 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_19 = new AliDielectronPID("PIDcut_19","PIDcut_19");
  PIDcut_19->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  PIDcut_19->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_19->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 0.5 , 0. ,100., kFALSE);
  PIDcut_19->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *PIDcut_20 = new AliDielectronPID("PIDcut_20","PIDcut_20");
  PIDcut_20->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_20->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_20->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.0 , 0. ,100., kFALSE);
  PIDcut_20->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 2.5 , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);

  double kaon_offset = 0.; double proton_offset = 0; double nSigmaPion = 4.;

                         kaon_offset = 0.5; proton_offset = 0.; nSigmaPion = 3.5;
  AliDielectronPID*      Jeromian_00_hadron_cut = new AliDielectronPID     ("Jeromian_00_hadron_cut","Jeromian_00_hadron_cut");
  AliDielectronPID*      Jeromian_00_ele_incl   = new AliDielectronPID     ("Jeromian_00_ele_incl",  "Jeromian_00_ele_incl");
  AliDielectronCutGroup* Jeromian_00_sum        = new AliDielectronCutGroup("Jeromian_00_sum",       "Jeromian_00_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_00_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_00_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_00_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0, 3.0+kaon_offset , 0. ,100., kTRUE);
                         Jeromian_00_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_00_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_00_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_00_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_00_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_00_sum->AddCut(Jeromian_00_hadron_cut);
                         Jeromian_00_sum->AddCut(Jeromian_00_ele_incl);




                         kaon_offset = 0.5; proton_offset = 0.; nSigmaPion = 3.5;
  AliDielectronPID*      Jeromian_01_hadron_cut = new AliDielectronPID     ("Jeromian_01_hadron_cut","Jeromian_01_hadron_cut");
  AliDielectronPID*      Jeromian_01_ele_incl   = new AliDielectronPID     ("Jeromian_01_ele_incl",  "Jeromian_01_ele_incl");
  AliDielectronCutGroup* Jeromian_01_sum        = new AliDielectronCutGroup("Jeromian_01_sum",       "Jeromian_01_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0, 3.0+kaon_offset , 0. ,100., kTRUE);
                         Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_01_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_01_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_01_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_01_sum->AddCut(Jeromian_01_hadron_cut);
                         Jeromian_01_sum->AddCut(Jeromian_01_ele_incl);

                         kaon_offset = 0.5; proton_offset = 0.; nSigmaPion = 3.5;
  AliDielectronPID*      Jeromian_02_hadron_cut = new AliDielectronPID     ("Jeromian_02_hadron_cut","Jeromian_02_hadron_cut");
  AliDielectronPID*      Jeromian_02_ele_incl   = new AliDielectronPID     ("Jeromian_02_ele_incl",  "Jeromian_02_ele_incl");
  AliDielectronCutGroup* Jeromian_02_sum        = new AliDielectronCutGroup("Jeromian_02_sum",       "Jeromian_02_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_02_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_02_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_02_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_02_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_02_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_02_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_02_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_02_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_02_sum->AddCut(Jeromian_02_hadron_cut);
                         Jeromian_02_sum->AddCut(Jeromian_02_ele_incl);

                         kaon_offset = -0.5; proton_offset = -0.5; nSigmaPion = 4;
  AliDielectronPID*      Jeromian_03_hadron_cut = new AliDielectronPID     ("Jeromian_03_hadron_cut","Jeromian_03_hadron_cut");
  AliDielectronPID*      Jeromian_03_ele_incl   = new AliDielectronPID     ("Jeromian_03_ele_incl",  "Jeromian_03_ele_incl");
  AliDielectronCutGroup* Jeromian_03_sum        = new AliDielectronCutGroup("Jeromian_03_sum",       "Jeromian_03_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_03_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_03_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_03_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_03_sum->AddCut(Jeromian_03_hadron_cut);
                         Jeromian_03_sum->AddCut(Jeromian_03_ele_incl);

                         kaon_offset = -0.5; proton_offset = -0.5; nSigmaPion = 4;
  // AliDielectronPID*      Jeromian_03_hadron_cut = new AliDielectronPID     ("Jeromian_03_hadron_cut","Jeromian_03_hadron_cut");
  AliDielectronPID*      Jeromian_03_ele_incl_TOF   = new AliDielectronPID     ("Jeromian_03_ele_incl_TOF",  "Jeromian_03_ele_incl_TOF");
  AliDielectronPID*      Jeromian_03_ele_incl_TRD   = new AliDielectronPID     ("Jeromian_03_ele_incl_TRD",  "Jeromian_03_ele_incl_TRD");
  AliDielectronCutGroup* Jeromian_03_sum_TRDincl= new AliDielectronCutGroup("Jeromian_03_sum_TRDincl",       "Jeromian_03_sum_TRDincl",AliDielectronCutGroup::kCompOR);
                         // Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         // Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         // Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         // Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         // Jeromian_03_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_03_ele_incl_TOF  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_03_ele_incl_TOF  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_03_ele_incl_TOF  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);


                        //  Double_t trdEleProb = 0.99;
                        // AliDielectronVarCuts *TRDtrkCuts = new AliDielectronVarCuts("TRDtrkCuts","TRDtrkCuts");
                       	// TRDtrkCuts->AddCut(ADVM::kTRDpidQuality, 2.5, 7.0, kTRUE); // tracklets for PID
                       	// TRDtrkCuts->AddCut(ADVM::kTRDchi2Trklt, .0, 5.0, kTRUE);
                       	// TRDtrkCuts->AddCut(ADVM::kTRDprob2DEle, trdEleProb, 1., kTRUE);


  AliDielectronCutGroup* CutGroupOnlyTPCPlusTRD= new AliDielectronCutGroup("CutGroupOnlyTPCPlusTRD",       "CutGroupOnlyTPCPlusTRD",AliDielectronCutGroup::kCompAND);
  AliDielectronPID*      Jeromian_03_onlyTPCnSigma   = new AliDielectronPID     ("Jeromian_03_onlyTPCnSigma",  "Jeromian_03_onlyTPCnSigma");
                          Jeromian_03_onlyTPCnSigma  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.,  0., 100., kFALSE);


  AliDielectronVarCuts *TRD_PIDcuts = new AliDielectronVarCuts("TRD_PIDcuts","TRD_PIDcuts");
                         TRD_PIDcuts->AddCut(AliDielectronVarManager::kTRDpidQuality, 3.5, 7.0, kFALSE); // tracklets for PID
                         TRD_PIDcuts->AddCut(AliDielectronVarManager::kTRDchi2Trklt, .0, 5.0, kFALSE);
                         TRD_PIDcuts->AddCut(AliDielectronVarManager::kTRDprob2DEle, 0.99, 1., kFALSE); // Double_t trdEleProb = 0.99;

                         CutGroupOnlyTPCPlusTRD->AddCut(Jeromian_03_onlyTPCnSigma);
                         CutGroupOnlyTPCPlusTRD->AddCut(TRD_PIDcuts);


                         Jeromian_03_sum_TRDincl->AddCut(Jeromian_03_hadron_cut);
                         Jeromian_03_sum_TRDincl->AddCut(Jeromian_03_ele_incl_TOF);
                         Jeromian_03_sum_TRDincl->AddCut(CutGroupOnlyTPCPlusTRD);

                         kaon_offset = +0.5; proton_offset = +0.5; nSigmaPion = 3.5;
  AliDielectronPID*      Jeromian_04_hadron_cut = new AliDielectronPID     ("Jeromian_04_hadron_cut","Jeromian_04_hadron_cut");
  AliDielectronPID*      Jeromian_04_ele_incl   = new AliDielectronPID     ("Jeromian_04_ele_incl",  "Jeromian_04_ele_incl");
  AliDielectronCutGroup* Jeromian_04_sum        = new AliDielectronCutGroup("Jeromian_04_sum",       "Jeromian_04_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_04_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_04_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_04_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_04_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_04_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_04_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_04_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_04_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_04_sum->AddCut(Jeromian_04_hadron_cut);
                         Jeromian_04_sum->AddCut(Jeromian_04_ele_incl);

                         kaon_offset = +0.0; proton_offset = +0.5; nSigmaPion = 4.;
  AliDielectronPID*      Jeromian_05_hadron_cut = new AliDielectronPID     ("Jeromian_05_hadron_cut","Jeromian_05_hadron_cut");
  AliDielectronPID*      Jeromian_05_ele_incl   = new AliDielectronPID     ("Jeromian_05_ele_incl",  "Jeromian_05_ele_incl");
  AliDielectronCutGroup* Jeromian_05_sum        = new AliDielectronCutGroup("Jeromian_05_sum",       "Jeromian_05_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_05_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_05_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_05_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_05_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_05_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_05_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_05_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_05_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_05_sum->AddCut(Jeromian_05_hadron_cut);
                         Jeromian_05_sum->AddCut(Jeromian_05_ele_incl);

                         kaon_offset = -0.5; proton_offset = +0.5; nSigmaPion = 4.;
  AliDielectronPID*      Jeromian_06_hadron_cut = new AliDielectronPID     ("Jeromian_06_hadron_cut","Jeromian_06_hadron_cut");
  AliDielectronPID*      Jeromian_06_ele_incl   = new AliDielectronPID     ("Jeromian_06_ele_incl",  "Jeromian_06_ele_incl");
  AliDielectronCutGroup* Jeromian_06_sum        = new AliDielectronCutGroup("Jeromian_06_sum",       "Jeromian_06_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_06_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_06_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_06_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_06_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_06_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_06_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_06_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_06_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_06_sum->AddCut(Jeromian_06_hadron_cut);
                         Jeromian_06_sum->AddCut(Jeromian_06_ele_incl);

                         kaon_offset = -0.5; proton_offset = -0.5; nSigmaPion = 3.5;
  AliDielectronPID*      Jeromian_07_hadron_cut = new AliDielectronPID     ("Jeromian_07_hadron_cut","Jeromian_07_hadron_cut");
  AliDielectronPID*      Jeromian_07_ele_incl   = new AliDielectronPID     ("Jeromian_07_ele_incl",  "Jeromian_07_ele_incl");
  AliDielectronCutGroup* Jeromian_07_sum        = new AliDielectronCutGroup("Jeromian_07_sum",       "Jeromian_07_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_07_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_07_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_07_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_07_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_07_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_07_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_07_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_07_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_07_sum->AddCut(Jeromian_07_hadron_cut);
                         Jeromian_07_sum->AddCut(Jeromian_07_ele_incl);

                         kaon_offset = 0; proton_offset = 0; nSigmaPion = 3.5;
  AliDielectronPID*      Jeromian_08_hadron_cut = new AliDielectronPID     ("Jeromian_08_hadron_cut","Jeromian_08_hadron_cut");
  AliDielectronPID*      Jeromian_08_ele_incl   = new AliDielectronPID     ("Jeromian_08_ele_incl",  "Jeromian_08_ele_incl");
  AliDielectronCutGroup* Jeromian_08_sum        = new AliDielectronCutGroup("Jeromian_08_sum",       "Jeromian_08_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_08_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_08_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_08_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_08_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_08_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_08_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_08_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_08_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_08_sum->AddCut(Jeromian_08_hadron_cut);
                         Jeromian_08_sum->AddCut(Jeromian_08_ele_incl);

                         kaon_offset = 0.5; proton_offset = -0.5; nSigmaPion = 4.5;
  AliDielectronPID*      Jeromian_09_hadron_cut = new AliDielectronPID     ("Jeromian_09_hadron_cut","Jeromian_09_hadron_cut");
  AliDielectronPID*      Jeromian_09_ele_incl   = new AliDielectronPID     ("Jeromian_09_ele_incl",  "Jeromian_09_ele_incl");
  AliDielectronCutGroup* Jeromian_09_sum        = new AliDielectronCutGroup("Jeromian_09_sum",       "Jeromian_09_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_09_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_09_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_09_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_09_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_09_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_09_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_09_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_09_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_09_sum->AddCut(Jeromian_09_hadron_cut);
                         Jeromian_09_sum->AddCut(Jeromian_09_ele_incl);

                         kaon_offset = 0.5; proton_offset = 0.; nSigmaPion = 4.5;
  AliDielectronPID*      Jeromian_10_hadron_cut = new AliDielectronPID     ("Jeromian_10_hadron_cut","Jeromian_10_hadron_cut");
  AliDielectronPID*      Jeromian_10_ele_incl   = new AliDielectronPID     ("Jeromian_10_ele_incl",  "Jeromian_10_ele_incl");
  AliDielectronCutGroup* Jeromian_10_sum        = new AliDielectronCutGroup("Jeromian_10_sum",       "Jeromian_10_sum",AliDielectronCutGroup::kCompOR);
                         Jeromian_10_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
                         Jeromian_10_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_10_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0+kaon_offset,   3.0+kaon_offset   , 0. ,100., kTRUE);
                         Jeromian_10_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0+proton_offset, 3.0+proton_offset , 0. ,100., kTRUE);
                         Jeromian_10_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
                         Jeromian_10_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
                         Jeromian_10_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
                         Jeromian_10_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
                         Jeromian_10_sum->AddCut(Jeromian_10_hadron_cut);
                         Jeromian_10_sum->AddCut(Jeromian_10_ele_incl);




  AliDielectronPID *PIDcut_5_looserPionRejection = new AliDielectronPID("PIDcut_5_looserPionRejection","PIDcut_5_looserPionRejection");
  PIDcut_5_looserPionRejection->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.0 , 0. ,100., kFALSE);
  PIDcut_5_looserPionRejection->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_5_looserPionRejection->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5 , 0. ,100., kFALSE);
  PIDcut_5_looserPionRejection->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDnoneExisting = new AliDielectronPID("PIDnoneExisting","PIDnoneExisting");
  PIDnoneExisting->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1000.0, 1000.0 , 0. ,100., kFALSE);


  // ################## PID for Clean Samples
  AliDielectronPID *pid_Pt200_cut5_woTPCelecut = new AliDielectronPID("pid_Pt200_cut5_woTPCelecut","pid_Pt200_cut5_woTPCelecut");
  pid_Pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5, 0. ,100., kFALSE);
  pid_Pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 4.5 , 0. ,100., kTRUE);
  pid_Pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_Pt200_cut5_woPionRej = new AliDielectronPID("pid_Pt200_cut5_woPionRej","pid_Pt200_cut5_woPionRej");
  pid_Pt200_cut5_woPionRej->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5, 0. ,100., kFALSE);
  pid_Pt200_cut5_woPionRej->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0, 0. ,100., kFALSE);
  pid_Pt200_cut5_woPionRej->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *electron_pt200_cut5_woTPCelecut = new AliDielectronPID("electron_pt200_cut5_woTPCelecut","electron_pt200_cut5_woTPCelecut");
  electron_pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  electron_pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 4.5 , 0. ,100., kTRUE);

  AliDielectronPID *kaon_pt200_cut5_woTPCelecut = new AliDielectronPID("kaon_pt200_cut5_woTPCelecut","kaon_pt200_cut5_woTPCelecut");
  kaon_pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kKaon, -3, 3, 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  kaon_pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99. , 4.5 , 0. ,100., kTRUE);

  AliDielectronPID *proton_pt200_cut5_woTPCelecut = new AliDielectronPID("proton_pt200_cut5_woTPCelecut","proton_pt200_cut5_woTPCelecut");
  proton_pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -3, 3, 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  proton_pt200_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99. , 4.5 , 0. ,100., kTRUE);

  AliDielectronPID *pion_pt200_cut5_wTPCelecut = new AliDielectronPID("pion_pt200_cut5_wTPCelecut","pion_pt200_cut5_wTPCelecut");
  pion_pt200_cut5_wTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3, 3., 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  pion_pt200_cut5_wTPCelecut->AddCut(AliDielectronPID::kITS,AliPID::kPion, -3., 3., 0. ,100., kFALSE);
  pion_pt200_cut5_wTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3., 0. ,100., kFALSE);




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
  pid_electron->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0, 3. , 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE);
  pid_electron->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99.0,-5., 0. ,100., kTRUE);

  AliDielectronPID *pid_electron_woTOF = new AliDielectronPID("pid_electron_woTOF","pid_electron_woTOF");
  pid_electron_woTOF->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron_woTOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 2. , 0. ,100., kFALSE);
  pid_electron_woTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_electron_wTOFreq = new AliDielectronPID("pid_electron_wTOFreq","pid_electron_wTOFreq");
  pid_electron_wTOFreq->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron_wTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 2. , 0. ,100., kFALSE);
  pid_electron_wTOFreq->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron_wTOFreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_electron_wTOFif = new AliDielectronPID("pid_electron_wTOFif","pid_electron_wTOFif");
  pid_electron_wTOFif->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron_wTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 2. , 0. ,100., kFALSE);
  pid_electron_wTOFif->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0, 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron_wTOFif->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  AliDielectronPID *pid_electron_wTOFhit = new AliDielectronPID("pid_electron_wTOFhit","pid_electron_wTOFhit");
  pid_electron_wTOFhit->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  pid_electron_wTOFhit->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 2. , 0. ,100., kFALSE);
  pid_electron_wTOFhit->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -99.0, 99. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  pid_electron_wTOFhit->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 5. , 0. ,100., kTRUE);

  // AliDielectronCutGroup* PDGlepton = new AliDielectronCutGroup("PDGlepton","PDGlepton",AliDielectronCutGroup::kCompOR);
  AliDielectronCutGroup* pid_PDG_electron_woTOF = new AliDielectronCutGroup("pid_PDG_electron_woTOF","pid_PDG_electron_woTOF",AliDielectronCutGroup::kCompAND);
  pid_PDG_electron_woTOF->AddCut(PDGlepton);
  pid_PDG_electron_woTOF->AddCut(pid_electron_woTOF);
  AliDielectronCutGroup* pid_PDG_electron_wTOFreq = new AliDielectronCutGroup("pid_PDG_electron_wTOFreq","pid_PDG_electron_wTOFreq",AliDielectronCutGroup::kCompAND);
  pid_PDG_electron_wTOFreq->AddCut(PDGlepton);
  pid_PDG_electron_wTOFreq->AddCut(pid_electron_wTOFreq);
  AliDielectronCutGroup* pid_PDG_electron_wTOFif = new AliDielectronCutGroup("pid_PDG_electron_wTOFif","pid_PDG_electron_wTOFif",AliDielectronCutGroup::kCompAND);
  pid_PDG_electron_wTOFif->AddCut(PDGlepton);
  pid_PDG_electron_wTOFif->AddCut(pid_electron_wTOFif);
  AliDielectronCutGroup* pid_PDG_electron_wTOFhit = new AliDielectronCutGroup("pid_PDG_electron_wTOFhit","pid_PDG_electron_wTOFhit",AliDielectronCutGroup::kCompAND);
  pid_PDG_electron_wTOFhit->AddCut(PDGlepton);
  pid_PDG_electron_wTOFhit->AddCut(pid_electron_wTOFhit);


  AliDielectronPID *pid_pion_looseelctron = new AliDielectronPID("pid_pion_looseelctron","pid_pion_looseelctron");
  pid_pion_looseelctron->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -100., 100., 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronPID *pid_pion = new AliDielectronPID("pid_pion","pid_pion");
  pid_pion->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3., 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
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
  // pid_pion_wTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);

  AliDielectronPID *pid_pion_woTOF = new AliDielectronPID("pid_pion_woTOF","pid_pion_woTOF");
  // pid_pion_woTOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);

  AliDielectronPID *pid_pion_wTOF3Sreq = new AliDielectronPID("pid_pion_wTOF3Sreq","pid_pion_wTOF3Sreq");
  pid_pion_wTOF3Sreq->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kRequire);
  // pid_pion_wTOF3Sreq->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);

  AliDielectronPID *pid_pion_wTOF3Sif = new AliDielectronPID("pid_pion_wTOF3Sif","pid_pion_wTOF3Sif");
  pid_pion_wTOF3Sif->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3., 3., 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  // pid_pion_wTOF3Sif->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -3.0, 3. , 0. ,100., kFALSE);



  //-----------------------------------------------
  // Now see what Config actually loads and assemble final cuts
  //-----------------------------------------------
  switch (AnaCut.GetPIDAna()) {

    case kPbPb2015_Pt200_PID_cutoff_pion_kaon_proton:
    // Cut out pion/kaon/proton band LHC15o but refilled when particle in TOF electron band
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PID_cutoff_pion_kaon_proton_cg, AnaCut);;
      break;
    case kPbPb2015_Pt200_noPID:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDnoneExisting, AnaCut);;
      break;

    case kPID_Jeromian_00:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_00_sum, AnaCut);;
      break;
    case kPID_Jeromian_00_TPConly:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_00_hadron_cut, AnaCut);;
      break;
    case kPID_Jeromian_00_TOFonly:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_00_ele_incl, AnaCut);;
      break;

    case kPID_Jeromian_01:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_01_sum, AnaCut);;
      break;

    case kPIDcut_0_onlyLooseTPC:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange20000to100000, pid_pion_looseelctron, AnaCut);;
      break;

    case kPDGelectron:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PDGlepton, AnaCut);;
      break;

    case kPIDcut_5_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_5, AnaCut);;
      break;
    case kPIDcut_11_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, PIDcut_11, AnaCut);;
      break;
    case kPIDcut_23_pt400:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange400to8000, Jeromian_03_sum, AnaCut);;
      break;


    case kPIDcut_1_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_1, AnaCut);;
      break;
    case kPIDcut_2_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_2, AnaCut);;
      break;
    case kPIDcut_3_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_3, AnaCut);;
      break;
    case kPIDcut_4_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_4, AnaCut);;
      break;
    case kPIDcut_5_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_5, AnaCut);;
      break;
    case kPIDcut_5_pt200_looserPionRejection:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_5_looserPionRejection, AnaCut);;
      break;
    case kPIDcut_6_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_6, AnaCut);;
      break;
    case kPIDcut_7_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_7, AnaCut);;
      break;
    case kPIDcut_8_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_8, AnaCut);;
      break;
    case kPIDcut_9_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_9, AnaCut);;
      break;
    case kPIDcut_10_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_10, AnaCut);;
      break;
    case kPIDcut_11_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_11, AnaCut);;
      break;
    case kPIDcut_12_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_12, AnaCut);;
      break;
    case kPIDcut_13_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_13, AnaCut);;
      break;
    case kPIDcut_14_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_14, AnaCut);;
      break;
    case kPIDcut_15_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_15, AnaCut);;
      break;
    case kPIDcut_16_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_16, AnaCut);;
      break;
    case kPIDcut_17_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_17, AnaCut);;
      break;
    case kPIDcut_18_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_18, AnaCut);;
      break;
    case kPIDcut_19_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_19, AnaCut);;
      break;
    case kPIDcut_20_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDcut_20, AnaCut);;
      break;
    case kPIDcut_21_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_01_sum, AnaCut);;
      break;
    case kPIDcut_22_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_02_sum, AnaCut);;
      break;
    case kPIDcut_23_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_03_sum, AnaCut);;
      break;
    case kPIDcut_24_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_04_sum, AnaCut);;
      break;
    case kPIDcut_25_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_05_sum, AnaCut);;
      break;
    case kPIDcut_26_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_06_sum, AnaCut);;
      break;
    case kPIDcut_27_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_07_sum, AnaCut);;
      break;
    case kPIDcut_28_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_08_sum, AnaCut);;
      break;
    case kPIDcut_29_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_09_sum, AnaCut);;
      break;
    case kPIDcut_30_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_10_sum, AnaCut);;
      break;
    case kPIDcut_23_pt200_TRDincl:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_03_sum_TRDincl, AnaCut);;
      break;


      // ###################### CONTAMINATION FOR CUT 5 ##########################################

    case kPbPb2015_Pt200_cut5_woTPCelecut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_Pt200_cut5_woTPCelecut, AnaCut);;
      break;
    case kPbPb2015_Pt200_cut5_woPionRej:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_Pt200_cut5_woPionRej, AnaCut);;
      break;
    case kPbPb2015_pure_electron_pt200_woTPCelecut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, electron_pt200_cut5_woTPCelecut, AnaCut);;
      break;
    case kPbPb2015_pure_kaon_pt200_woTPCelecut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, kaon_pt200_cut5_woTPCelecut, AnaCut);;
      break;
    case kPbPb2015_pure_proton_pt200_woTPCelecut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, proton_pt200_cut5_woTPCelecut, AnaCut);;
      break;
    case kPbPb2015_pure_pion_pt200_wTPCelecut:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pion_pt200_cut5_wTPCelecut, AnaCut);;
      break;


      // ################# CLEAN SAMPLES ##################

      case kPbPb2015_pure_electron_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_electron, AnaCut);;
      break;
    case kPbPb2015_pure_pion_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_pion, AnaCut);;
      break;

      // TOF Efficiency Study with V0 Electrons:
      case kPbPb2015_pure_electron_pt200_woTOF:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_electron_woTOF, AnaCut);;
        break;
      case kPbPb2015_pure_electron_pt200_wTOFhit:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_electron_wTOFhit, AnaCut);;
        break;
      case kPbPb2015_pure_electron_pt200_wTOFreq:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_electron_wTOFreq, AnaCut);;
        break;
      case kPbPb2015_pure_electron_pt200_wTOFif:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_electron_wTOFif, AnaCut);;
        break;
      case kPbPb2015_PDG_pure_electron_pt200_woTOF:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_PDG_electron_woTOF, AnaCut);;
        break;
      case kPbPb2015_PDG_pure_electron_pt200_wTOFhit:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_PDG_electron_wTOFhit, AnaCut);;
        break;
      case kPbPb2015_PDG_pure_electron_pt200_wTOFreq:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_PDG_electron_wTOFreq, AnaCut);;
        break;
      case kPbPb2015_PDG_pure_electron_pt200_wTOFif:
        pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, pid_PDG_electron_wTOFif, AnaCut);;
        break;

    case kPbPb2015_Pt100_ResolutionCuts:
      AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
      cgPIDCutsAna->AddCut(etaRange120);
      cgPIDCutsAna->AddCut(ptRange100toINF);
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
    case kResolutionTrackCuts:
      std::cout << "kResolutionTrackCuts" << std::endl;
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);

      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0); // means at least 2 with PID
      trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,  15.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   2.1); // means 0 and 1 shared Cluster

      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
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
    case kNoTrackCuts:
      std::cout << "kNoTrackCuts" << std::endl;
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -100.0,   100.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);

      // trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0); // means at least 2 with PID
      // trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster
      // trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014

      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      // trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      // cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      // cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
    case kZ0cuts:
      std::cout << "kZ0cuts" << std::endl;
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);

      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      2.0, 100.0); // means at least 2 with PID
      // trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1); // lower limit 0.8 in most filterbits! // 1.1 since 26.02.2014

      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      trackCutsDiel->SetAODFilterBit(1<<4); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      // trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      // cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
    case kSPDfirst_PDGCodePion:
      std::cout << "SPDfirstPion" << std::endl;
      AliDielectronVarCuts* pdgCodepionP = new AliDielectronVarCuts("PDGPionP","PDGPionP");
      pdgCodepionP->AddCut(AliDielectronVarManager::kPdgCode,  210.5,   211.5);

      AliDielectronVarCuts* pdgCodepionN = new AliDielectronVarCuts("PDGPionN","PDGPionN");
      pdgCodepionN->AddCut(AliDielectronVarManager::kPdgCode, -210.5,  -211.5);
      AliDielectronCutGroup* pdgCodepion = new AliDielectronCutGroup("pdgCodepion","pdgCodepion",AliDielectronCutGroup::kCompOR);
      pdgCodepion->AddCut(pdgCodepionP);
      pdgCodepion->AddCut(pdgCodepionN);

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
      cgTrackCutsAnaSPDfirst->AddCut(pdgCodepion);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      // trackCuts = GetTrackCuts(kSPD_bit4);
      // break;
    case kITSSA:
      trackCuts = GetTrackCuts(kITSSA_bit1);
      break;
      case kTRACKcut_5_woSharedCluster:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_5_0SharedCluster:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;

      case kTRACKcut_5_woSharedCluster_noSPDfirst:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_5_0SharedCluster_noSPDfirst:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_5_noSPDfirst:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_5_noConvRejection:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;


      case kTRACKcut_1:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_2:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);

        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_3:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_4:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_5:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_6:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_7:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_8:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_9:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_10:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_11:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_12:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_13:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_14:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_15:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_16:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_17:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_18:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_19:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_20:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_21:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_22:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_23:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_24:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_25:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_26:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_27:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_28:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_29:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     -0.1,   0.1); // means 0 shared Cluster

        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        trackCuts = cgTrackCutsAnaSPDfirst;
      break;
      case kTRACKcut_30:
        AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   3.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
        trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
        AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
        double delta = 0.00001;
        AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
        trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
        AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
        trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
        AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
        trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
        AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
        trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
        AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
        trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
        AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
        trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
        SharedClusterCut->AddCut(trackCutsSharedCluster0);
        SharedClusterCut->AddCut(trackCutsSharedCluster2);
        SharedClusterCut->AddCut(trackCutsSharedCluster4);
        SharedClusterCut->AddCut(trackCutsSharedCluster8);
        SharedClusterCut->AddCut(trackCutsSharedCluster16);
        SharedClusterCut->AddCut(trackCutsSharedCluster32);
        AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
        trackCutsDiel->SetAODFilterBit(1<<4);
        trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
        cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
        cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
        cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
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
      // gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
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
      cout << "Bit 4" << endl;
      esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
      esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
      esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
      esdTrackCutsH->SetDCAToVertex2D(kTRUE);

      break;
      //default: cout << "No ESD Track Cut defined " << endl;
  }
  return esdTrackCutsH;
}
