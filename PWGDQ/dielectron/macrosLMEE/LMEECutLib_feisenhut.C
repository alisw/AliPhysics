#ifndef LMEECutLib_feisenhut
#define LMEECutLib_feisenhut

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
    // SetPairCutsAna(LMEECutLib::kNoPairCutsAna);

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
    // Analysis cuts
    kNoPID_Pt20,
    kNoPID_Pt75,
    kNoPID_Pt200,
    kPID_Jeromian_01,
    kPID_Jeromian_01_PreFilter,
    kPID_Jeromian_01_pt200,
    kPIDcut_1_pt75,
    kPIDcut_TEST,
    kNoPID_noKinCuts
    };
  enum LMEEPIDPre{
    kStandardPre
  };
  // Possible Track Selections
  enum LMEETrackSelectionAna{
    kNoTrackCuts,
    kTRACKcut_1,
    kTRACKcut_2,
    kTRACKcut_2_PreFilter,
    kTRACKcut_1_secondary,
    kDefaultNoTrackCuts,
    kV0track
  };
  enum LMEETrackSelectionPre{
    kPrefilter_cut1
  };
  // enum LMEETrackCuts{
  //   kPbPb2015_V0_tight,
  //   kSPD_bit4,
  //   kITSSA_bit1,
  //   kNoTrackCuts
  // };
  enum LMEEPairCutsAna{
    kPairCutsAna, // Cut off (theta < 0.05) && (Minv < 0.02)
    kNoPairCutsAna, // No Cuts applied, since 18.02.2014
    kV0pair,
    kV0pair_PreFilter,
    kV0_onlyCos,
    kV0_onlyChi2NDF,
    kV0_onlyLegDist,
    kV0_onlyR,
    kV0_onlyPsiPair,
    kV0_onlyM,
    kV0_onlyArmPt,
    kV0_onlyArmAlpha,
    kV0_wCosChi2,
    kV0_wCosChi2LegDist,
    kV0_wCosChi2LegDistR,
    kV0_wCosChi2LegDistRPsiPair,
    kV0_wCosChi2LegDistRPsiPairM,
    kV0_wCosChi2LegDistRPsiPairMArmPt,
    kV0_wCosChi2LegDistRPsiPairMArmPtAlpha,
    kV0_wPairM,
    kV0_wPairMArmPt,
    kV0_wPairMArmPtAlpha
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
    kPbPb0010,     // 0%-10%
    kPbPb1020,
    kPbPb2030,
    kPbPb3040,
    kPbPb4050,
    kPbPb5060,
    kPbPb6070,
    kPbPb7080,
    kPbPb8090,
    kPbPb1050, //10%-50%
    kPbPb5080,  //50%-80%
    kPbPb0080,     //0%-80%
    kPbPb0020,
    kPbPb2080,
    kPbPb2040,
    kPbPb6080,
    kPbPb0005,
    kPbPb0510,
    kPbPb1015,
    kPbPb1520
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
                                                                                // Printf("%d DEBUG CutLib before PID " , __LINE__);
    cgPIDCutsAna->AddCut(PID);
                                                                                // Printf("%d DEBUG CutLib after PID " , __LINE__);

    cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
                                                                                // Printf("%d DEBUG CutLib after GetTrackSelectionAna " , __LINE__);
    cgPIDCutsAna->AddCut(GetPairCutsAna(AnaCut));
                                                                                // Printf("%d DEBUG CutLib after GetPairCutsAna " , __LINE__);
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
    case kPbPb0010:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0010");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 0.001, 10.);
      break;
    case kPbPb1020:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb1020");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 10., 20.);
      break;
    case kPbPb2030:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb2030");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 20., 30.);
      break;
    case kPbPb3040:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb3040");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 30., 40.);
      break;
    case kPbPb4050:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb4050");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 40., 50.);
      break;
    case kPbPb5060:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb5060");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 50., 60.);
      break;
    case kPbPb6070:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb6070");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 60., 70.);
      break;
    case kPbPb7080:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb7080");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 70., 80.);
      break;
    case kPbPb8090:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb8090");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 80., 90.);
      break;
    case kPbPb1050:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb1050");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 10., 50.);
      break;
    case kPbPb5080:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb5080");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 50., 80.);
      break;
    case kPbPb0080:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0080");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 0.001, 80.);
      break;
    case kPbPb2040:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0020");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 20., 40.);
      break;
    case kPbPb6080:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb6080");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 60., 80.);
      break;
    case kPbPb0005:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0005");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 0., 5.);
      break;
    case kPbPb0510:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb0510");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 5., 10.);
      break;
    case kPbPb1015:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb1015");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 10., 15.);
      break;
    case kPbPb1520:
      centCuts = new AliDielectronVarCuts("centCuts","kPbPb1520");
      centCuts->AddCut(AliDielectronVarManager::kCentralityNew, 15., 20.);
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

    case kV0pair:
      // primarily meant for inclusion, for quite pure sample...
      std::cout << "Using kV0 Cutsetting" << std::endl;
      // AliDielectronV0Cuts *gammaV0Finder = new AliDielectronV0Cuts("gammaV0Finder","gammaV0Finder");
      // gammaV0Finder->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
      AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      // gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.06, kFALSE);
      gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.9,  0.9, kFALSE); // should increase purity...
      // gammaV0Cuts->SetExcludeTracks(kTRUE);
      // gammaV0Cuts->SetExcludeTracks(kFALSE); (standard?)
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      // trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      // cgTrackCutsV0select->AddCut(gammaV0Finder);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      // cgTrackCutsV0select->AddCut(trackCutsAOD);
      pairCuts = cgTrackCutsV0select;
      break;

      case kV0pair_PreFilter:
        // primarily meant for inclusion, for quite pure sample...
        std::cout << "Using kV0_PreFilter Cutsetting" << std::endl;
        // AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
        AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
        gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.8,  1.0,  kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.3, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.2, kFALSE);
        gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -1.0,   1.0,  kFALSE); // should increase purity...
        cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
        cgTrackCutsV0select->AddCut(gammaV0Cuts);
        pairCuts = cgTrackCutsV0select;
        break;


      // case kV0_onlyCos:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyCos Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyChi2NDF:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyChi2NDF Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyLegDist:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyLegDist Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyR:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyR Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyPsiPair:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyPsiPair Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyM:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyM Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyArmPt:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyArmPt Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.06, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_onlyArmAlpha:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_onlyArmAlpha Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.9,  0.9, kFALSE); // should increase purity...
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2 Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2LegDist:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2LegDist Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2LegDistR:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2LegDistR Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2LegDistRPsiPair:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2LegDistRPsiPair Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2LegDistRPsiPairM:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2LegDistRPsiPairM Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2LegDistRPsiPairMArmPt:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2LegDistRPsiPairMArmPt Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.06, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wCosChi2LegDistRPsiPairMArmPtAlpha:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wCosChi2LegDistRPsiPairMArmPtAlpha Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle,              0.95,  1.0,  kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.06, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.9,  0.9, kFALSE); // should increase purity...
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wPairM:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wM Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wPairMArmPt:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wMArmPt Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.06, kFALSE);
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;
      //
      // case kV0_wPairMArmPtAlpha:
      //   // primarily meant for inclusion, for quite pure sample...
      //   std::cout << "Using kV0_wMArmPtAlpha Cutsetting" << std::endl;
      //   AliDielectronVarCuts* gammaV0Cuts =new AliDielectronVarCuts("gammaV0Cuts","gammaV0Cuts");
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.16, kFALSE);
      //   // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.06, kFALSE);
      //   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.9,  0.9, kFALSE); // should increase purity...
      //   cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      //   cgTrackCutsV0select->AddCut(gammaV0Cuts);
      //   pairCuts = cgTrackCutsV0select;
      //   break;

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
  AliDielectronVarCuts *etaRange110 = new AliDielectronVarCuts("etaRange110","etaRange110");
  etaRange110->AddCut(AliDielectronVarManager::kEta, -1.10, 1.10);
  AliDielectronVarCuts *etaRange120 = new AliDielectronVarCuts("etaRange120","etaRange120");
  etaRange120->AddCut(AliDielectronVarManager::kEta, -1.20, 1.20);
  AliDielectronVarCuts *etaRange150 = new AliDielectronVarCuts("etaRange150","etaRange150");
  etaRange150->AddCut(AliDielectronVarManager::kEta, -1.50, 1.50);
  // pt range:
  AliDielectronVarCuts *ptRange400to8000 = new AliDielectronVarCuts("ptRange400to8000","ptRange400to8000");
  ptRange400to8000->AddCut(AliDielectronVarManager::kPt, 0.4, 8.0);
  AliDielectronVarCuts *ptRange20000to100000 = new AliDielectronVarCuts("ptRange20000to100000","ptRange20000to100000");
  ptRange20000to100000->AddCut(AliDielectronVarManager::kPt, 4., 100.0);
  AliDielectronVarCuts *ptRange100to8000 = new AliDielectronVarCuts("ptRange100to8000","ptRange100to8000");
  ptRange100to8000->AddCut(AliDielectronVarManager::kPt, 0.1, 8.0);
  AliDielectronVarCuts *ptRange100toINF = new AliDielectronVarCuts("ptRange100toINF","ptRange100toINF");
  ptRange100toINF->AddCut(AliDielectronVarManager::kPt, 0.1, 100.0);
  AliDielectronVarCuts *ptRange75to8000 = new AliDielectronVarCuts("ptRange75to8000","ptRange75to8000");
  ptRange75to8000->AddCut(AliDielectronVarManager::kPt, 0.075, 8.0);
  AliDielectronVarCuts *ptRange50to8000 = new AliDielectronVarCuts("ptRange50to8000","ptRange50to8000");
  ptRange50to8000->AddCut(AliDielectronVarManager::kPt, 0.050, 8.0);
  AliDielectronVarCuts *ptRange20to8000 = new AliDielectronVarCuts("ptRange20to8000","ptRange20to8000");
  ptRange20to8000->AddCut(AliDielectronVarManager::kPt, 0.02, 8.0);
  AliDielectronVarCuts *ptRange200to8000 = new AliDielectronVarCuts("ptRange200to8000","ptRange200to8000");
  ptRange200to8000->AddCut(AliDielectronVarManager::kPt, 0.2, 8.0);
  AliDielectronVarCuts *ptRange20to100000 = new AliDielectronVarCuts("ptRange20to100000","ptRange20to100000");
  ptRange20to100000->AddCut(AliDielectronVarManager::kPt, 0.02, 100.0);


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
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE);
  // Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, 4. , 0. ,100., kTRUE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -3.5, 3.5. , 0. ,100., kTRUE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -2.5, 2.5 , 0. ,100., kTRUE);
  Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -2.5, 2.5 , 0. ,100., kTRUE);
  // Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *Jeromian_01_ele_incl = new AliDielectronPID("Jeromian_01_ele_incl","Jeromian_01_ele_incl");
  Jeromian_01_ele_incl->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -3.5, 3.5 , 0.3 ,100., kTRUE);
  // Jeromian_01_ele_incl->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE);
  Jeromian_01_ele_incl->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);

  AliDielectronCutGroup* Jeromian_01_sum = new AliDielectronCutGroup("Jeromian_01_sum","Jeromian_01_sum",AliDielectronCutGroup::kCompOR);
  Jeromian_01_sum->AddCut(Jeromian_01_hadron_cut);
  Jeromian_01_sum->AddCut(Jeromian_01_ele_incl);


  AliDielectronPID *Jeromian_01_PreFilter_hadron_cut = new AliDielectronPID("Jeromian_01_PreFilter_hadron_cut","Jeromian_01_PreFilter_hadron_cut");
  Jeromian_01_PreFilter_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -3. , 3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable    ,AliDielectronVarManager::kP);

  AliDielectronCutGroup* Jeromian_01_sum_PreFilter = new AliDielectronCutGroup("Jeromian_01_sum_PreFilter","Jeromian_01_sum_PreFilter",AliDielectronCutGroup::kCompOR);
  Jeromian_01_sum_PreFilter->AddCut(Jeromian_01_PreFilter_hadron_cut);


  AliDielectronPID *PIDcut_1 = new AliDielectronPID("PIDcut_1","PIDcut_1");
  PIDcut_1->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  PIDcut_1->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  PIDcut_1->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  PIDcut_1->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *PIDcut_TEST = new AliDielectronPID("PIDcut_TEST","PIDcut_TEST");
  PIDcut_TEST->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1000.0, 1000.0 , 0. ,100., kFALSE);
  // PIDcut_TEST->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2.0, 3.5 , 0. ,100., kFALSE);
  // PIDcut_TEST->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99, 3.5 , 0. ,100., kTRUE);
  // PIDcut_TEST->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  // PIDcut_TEST->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);




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




  //                        kaon_offset = 0.5; proton_offset = 0.; nSigmaPion = 3.5;
  // AliDielectronPID*      Jeromian_01_hadron_cut = new AliDielectronPID     ("Jeromian_01_hadron_cut","Jeromian_01_hadron_cut");
  // AliDielectronPID*      Jeromian_01_ele_incl   = new AliDielectronPID     ("Jeromian_01_ele_incl",  "Jeromian_01_ele_incl");
  // AliDielectronCutGroup* Jeromian_01_sum        = new AliDielectronCutGroup("Jeromian_01_sum",       "Jeromian_01_sum",AliDielectronCutGroup::kCompOR);
  //                        Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3.0, 3.  , 0. ,100., kFALSE);
  //                        Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
  //                        Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,     -3.0, 3.0+kaon_offset , 0. ,100., kTRUE);
  //                        Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTPC,AliPID::kProton,   -3.0, 3.0+proton_offset , 0. ,100., kTRUE);
  //                        Jeromian_01_hadron_cut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  //                        Jeromian_01_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -99.0, nSigmaPion  , 0. ,100., kTRUE);
  //                        Jeromian_01_ele_incl  ->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE);
  //                        Jeromian_01_ele_incl  ->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3.  , 0. ,100., kFALSE, AliDielectronPID::kRequire);
  //                        Jeromian_01_sum->AddCut(Jeromian_01_hadron_cut);
  //                        Jeromian_01_sum->AddCut(Jeromian_01_ele_incl);

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
  // PIDnoneExisting->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1000.0, 1000.0 , 0. ,100., kFALSE);


  // ################## PID for Clean Samples
  AliDielectronPID *pid_Pt75_cut5_woTPCelecut = new AliDielectronPID("pid_Pt75_cut5_woTPCelecut","pid_Pt75_cut5_woTPCelecut");
  pid_Pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5, 0. ,100., kFALSE);
  pid_Pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 4.5 , 0. ,100., kTRUE);
  pid_Pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *pid_Pt75_cut5_woPionRej = new AliDielectronPID("pid_Pt75_cut5_woPionRej","pid_Pt75_cut5_woPionRej");
  pid_Pt75_cut5_woPionRej->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5, 0. ,100., kFALSE);
  pid_Pt75_cut5_woPionRej->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0, 0. ,100., kFALSE);
  pid_Pt75_cut5_woPionRej->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);

  AliDielectronPID *electron_pt75_cut5_woTPCelecut = new AliDielectronPID("electron_pt75_cut5_woTPCelecut","electron_pt75_cut5_woTPCelecut");
  electron_pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3, 3, 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  electron_pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion,     -99. , 4.5 , 0. ,100., kTRUE);

  AliDielectronPID *kaon_pt75_cut5_woTPCelecut = new AliDielectronPID("kaon_pt75_cut5_woTPCelecut","kaon_pt75_cut5_woTPCelecut");
  kaon_pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kKaon, -3, 3, 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  kaon_pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99. , 4.5 , 0. ,100., kTRUE);

  AliDielectronPID *proton_pt75_cut5_woTPCelecut = new AliDielectronPID("proton_pt75_cut5_woTPCelecut","proton_pt75_cut5_woTPCelecut");
  proton_pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -3, 3, 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  proton_pt75_cut5_woTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -99. , 4.5 , 0. ,100., kTRUE);

  AliDielectronPID *pion_pt75_cut5_wTPCelecut = new AliDielectronPID("pion_pt75_cut5_wTPCelecut","pion_pt75_cut5_wTPCelecut");
  pion_pt75_cut5_wTPCelecut->AddCut(AliDielectronPID::kTOF,AliPID::kPion, -3, 3., 0.4 ,100., kFALSE, AliDielectronPID::kRequire);
  pion_pt75_cut5_wTPCelecut->AddCut(AliDielectronPID::kITS,AliPID::kPion, -3., 3., 0. ,100., kFALSE);
  pion_pt75_cut5_wTPCelecut->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3., 0. ,100., kFALSE);




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

    case kNoPID_Pt75:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange75to8000, PIDnoneExisting, AnaCut);;
      break;
    case kNoPID_Pt20:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange20to8000, PIDnoneExisting, AnaCut);;
      break;
    case kNoPID_Pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, PIDnoneExisting, AnaCut);;
      break;
    case kNoPID_noKinCuts:
      pidCuts = LMEECutLib::SetKinematics(etaRange150, ptRange20to100000, PIDnoneExisting, AnaCut);;
      break;
    case kPID_Jeromian_01:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange75to8000, Jeromian_01_sum, AnaCut);;
      break;
    case kPID_Jeromian_01_PreFilter:
      pidCuts = LMEECutLib::SetKinematics(etaRange110, ptRange50to8000, Jeromian_01_sum_PreFilter, AnaCut);;
      break;
    case kPID_Jeromian_01_pt200:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange200to8000, Jeromian_01_sum, AnaCut);;
      break;
    case kPIDcut_1_pt75:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange75to8000, PIDcut_1, AnaCut);;
      break;
    case kPIDcut_TEST:
      pidCuts = LMEECutLib::SetKinematics(etaRange080, ptRange75to8000, PIDcut_TEST, AnaCut);;
      break;

    default: cout << "No Analysis PID Cut defined " << endl;
  }
  return pidCuts;
}

AliAnalysisCuts* LMEECutLib::GetTrackSelectionAna(AnalysisCut AnaCut) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
  switch (AnaCut.GetTrackSelectionAna()) {
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
      trackCutsDiel->SetAODFilterBit(1<<0); // (=16) filterbit 4! //GetStandardITSTPCTrackCuts2011(kFALSE); loose DCA, 2D cut
      // trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      // cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
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
      // trackCutsDiel->SetAODFilterBit(1<<4);
      trackCutsDiel->SetAODFilterBit(1<<0);
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
      // trackCutsDiel->SetAODFilterBit(1<<4);
      trackCutsDiel->SetAODFilterBit(1<<0);
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;

    case kTRACKcut_2_PreFilter:
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kKinkIndex0,          0. );
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
      // AliDielectronCutGroup* SharedClusterCut = new AliDielectronCutGroup("SharedClusterCut","SharedClusterCut",AliDielectronCutGroup::kCompOR);
      // double delta = 0.00001;
      // AliDielectronVarCuts* trackCutsSharedCluster0 = new AliDielectronVarCuts("trackCutsSharedCluster0", "trackCutsSharedCluster0");
      // trackCutsSharedCluster0->AddCut(AliDielectronVarManager::kNclsSMapITS, 0-delta, 0+delta);
      // AliDielectronVarCuts* trackCutsSharedCluster2 = new AliDielectronVarCuts("trackCutsSharedCluster2", "trackCutsSharedCluster2");
      // trackCutsSharedCluster2->AddCut(AliDielectronVarManager::kNclsSMapITS, 2-delta, 2+delta);
      // AliDielectronVarCuts* trackCutsSharedCluster4 = new AliDielectronVarCuts("trackCutsSharedCluster4", "trackCutsSharedCluster4");
      // trackCutsSharedCluster4->AddCut(AliDielectronVarManager::kNclsSMapITS, 4-delta, 4+delta);
      // AliDielectronVarCuts* trackCutsSharedCluster8 = new AliDielectronVarCuts("trackCutsSharedCluster8", "trackCutsSharedCluster8");
      // trackCutsSharedCluster8->AddCut(AliDielectronVarManager::kNclsSMapITS, 8-delta, 8+delta);
      // AliDielectronVarCuts* trackCutsSharedCluster16 = new AliDielectronVarCuts("trackCutsSharedCluster16", "trackCutsSharedCluster16");
      // trackCutsSharedCluster16->AddCut(AliDielectronVarManager::kNclsSMapITS, 16-delta, 16+delta);
      // AliDielectronVarCuts* trackCutsSharedCluster32 = new AliDielectronVarCuts("trackCutsSharedCluster32", "trackCutsSharedCluster32");
      // trackCutsSharedCluster32->AddCut(AliDielectronVarManager::kNclsSMapITS, 32-delta, 32+delta);
      // SharedClusterCut->AddCut(trackCutsSharedCluster0);
      // SharedClusterCut->AddCut(trackCutsSharedCluster2);
      // SharedClusterCut->AddCut(trackCutsSharedCluster4);
      // SharedClusterCut->AddCut(trackCutsSharedCluster8);
      // SharedClusterCut->AddCut(trackCutsSharedCluster16);
      // SharedClusterCut->AddCut(trackCutsSharedCluster32);
      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      // trackCutsDiel->SetAODFilterBit(1<<4);
      trackCutsDiel->SetAODFilterBit(1<<0);
      // trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
      trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      // cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;

    case kTRACKcut_1_secondary:
      std::cout << "kTRACKcut_1_secondary" << std::endl;
      AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -0.1,   0.1  , kTRUE);  // kTrue in order to exclude selection
      trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      0.0, 100.0);     // offen
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
      trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 161.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
      trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS,     0.9,   10.); //different SharedClusterCut: means 1 and 10 shared Cluster
                                                                              // std::cout << "Number of AODCuts: " <<  trackCutsAOD->GetNCuts() << std::endl;
                                                                              // for (size_t i = 0; i < trackCutsAOD->GetNCuts(); i++) {
                                                                              //   std::cout << i+1 <<"-th cut name: " <<  trackCutsAOD->GetCutName(i) << std::endl;
                                                                              // }

      AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
      // trackCutsDiel->SetAODFilterBit(1<<4);
      trackCutsDiel->SetAODFilterBit(1<<0);
      // trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
      cgTrackCutsAnaSPDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSPDfirst","cgTrackCutsAnaSPDfirst",AliDielectronCutGroup::kCompAND);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsDiel);
      cgTrackCutsAnaSPDfirst->AddCut(trackCutsAOD);
      // cgTrackCutsAnaSPDfirst->AddCut(SharedClusterCut);
      trackCuts = cgTrackCutsAnaSPDfirst;
      break;

    case kV0track:
      // primarily meant for inclusion, for quite pure sample...
      AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("gammaV0Cuts","gammaV0Cuts");
      gammaV0Cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);  // kAll(default), kOffline or kOnTheFly
      // gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0,  kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0,  kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0,  kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.02, kFALSE);
      // gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // should increase purity...
      // gammaV0Cuts->SetExcludeTracks(kTRUE);
      // // gammaV0Cuts->SetExcludeTracks(kFALSE); (standard?)
      // AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
      // trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
      // trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);
      cgTrackCutsV0select = new AliDielectronCutGroup("cgTrackCutsV0select","cgTrackCutsV0select",AliDielectronCutGroup::kCompAND);
      cgTrackCutsV0select->AddCut(gammaV0Cuts);
      // cgTrackCutsV0select->AddCut(trackCutsAOD);
      trackCuts = cgTrackCutsV0select;
      break;

    // case kNone:
      // trackCuts = GetTrackCuts(kNoTrackCuts);
      // break;

    case kDefaultNoTrackCuts:
      cout << "No Track Cuts used - ok " << endl; // since 30.01.2020
      break;

    default: cout << "No Analysis Track Selection defined " << endl;
  }
  return trackCuts;
}

// AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet) {
//   cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
//   AliDielectronCutGroup* trackCuts=0x0;
//   switch (cutSet) {
//     case kNoTrackCuts:
//       AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
//       AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
//       trackCutsDiel->SetAODFilterBit(1<<4); //GetStandardITSTPCTrackCuts2011(kTRUE), SPD none, SDD first
//       trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
//
//       cgTrackCutsAnaSDDfirst = new AliDielectronCutGroup("cgTrackCutsAnaSDDfirst","cgTrackCutsAnaSDDfirst",AliDielectronCutGroup::kCompAND);
//       cgTrackCutsAnaSDDfirst->AddCut(trackCutsDiel);
//       cgTrackCutsAnaSDDfirst->AddCut(trackCutsAOD);
//       trackCuts = cgTrackCutsAnaSDDfirst;
//       break;
//
//     default: cout << "No Analysis Track Cut defined " << endl;
//   }
//   return trackCuts;
// }


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
