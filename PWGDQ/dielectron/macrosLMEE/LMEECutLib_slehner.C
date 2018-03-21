#ifndef LMEECutLib_slehner
#define LMEECutLib_slehner

//#include "AliDielectron.h"

class AnalysisCut{
public:
  AnalysisCut(){};

  //Setter
  void SetPIDAna(Int_t sPIDAna){PIDAna = sPIDAna;};
  void SetPIDPre(Int_t sPIDPre){PIDPre = sPIDPre;}
  void SetTrackSelectionAna(Int_t sTrackSelectionAna){TrackSelectionAna = sTrackSelectionAna;}
  void SetTrackSelectionPre(Int_t sTrackSelectionPre){TrackSelectionPre = sTrackSelectionPre;}


  void SetCentrality(Int_t sCentrality){Centrality = sCentrality;}
  void SetPreFilterType(Int_t sPreFilterType){PreFilterType = sPreFilterType;}
  void SetMixing(Int_t sMixing){Mixing = sMixing;}
  void SetESDTrackSelection(Int_t sESDTrackSelection){ESDTrackSelection = sESDTrackSelection;}
  void SetStandardCut(){


    SetMixing(LMEECutLib::kEventMixing_1);
    SetESDTrackSelection(LMEECutLib::kStandardESD);
  };

  //Getter
  Int_t GetPIDAna(){return PIDAna;}
  Int_t GetPIDPre(){return PIDPre;}
  Int_t GetTrackSelectionAna(){return TrackSelectionAna;}
  Int_t GetTrackSelectionPre(){return TrackSelectionPre;}


  Int_t GetCentrality(){return Centrality;}
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

  LMEECutLib() {}

  AliDielectronCutGroup*     GetEventCuts(Int_t cutSet);
  AliAnalysisCuts*            GetCentralityCuts(AnalysisCut AnaCut);
  AliDielectronTrackRotator*  GetTrackRotator(AnalysisCut AnaCut);


  AliDielectronPID* GetPIDCutsAna();

  AliAnalysisCuts* GetTrackSelectionAna(AnalysisCut AnaCut);
  AliDielectronCutGroup* GetTrackCuts(Int_t cutSet=0);
  AliDielectronCutGroup* LMEECutLib::SetKinematics(AliDielectronVarCuts *etaRange, AliDielectronVarCuts *ptRange, AliDielectronPID* PID, AnalysisCut AnaCut) {
    AliDielectronCutGroup* cgPIDCutsAna = new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
    cgPIDCutsAna->AddCut(etaRange);
    cgPIDCutsAna->AddCut(ptRange);
    cgPIDCutsAna->AddCut(PID);
    cgPIDCutsAna->AddCut(GetTrackSelectionAna(AnaCut));
    return cgPIDCutsAna;
  }


//  void SetEtaCorrection(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/

};

//void LMEECutLib::SetEtaCorrection(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
//  //
//  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
//  //
//  std::cout << "starting LMEECutLib::SetEtaCorrection()\n";
//  std::string file_name = "/home/cklein/LMeeAnaFW/005_TPC_Recalibration/summary/output.root";
//
//  TFile* _file = TFile::Open(file_name.c_str());
//  std::cout << _file << std::endl;
//  if (_file == 0x0){
//    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/c/cklein/data/output.root .");
//    std::cout << "Copy TPC correction from Alien" << std::endl;
//    _file = TFile::Open("output.root");
//  }
//  else {
//    std::cout << "Correction loaded" << std::endl;
//  }
//  if (runwise){
//    TObjArray* arr_mean = dynamic_cast<TObjArray*>(_file->Get("mean_correction_arr"));
//    TObjArray* arr_width =dynamic_cast<TObjArray*>( _file->Get("width_correction_arr"));
//    std::cout << arr_mean << " " << arr_width << std::endl;
//
//    die->SetWidthCorrArr(arr_width, kTRUE, corrXdim, corrYdim, corrZdim);
//    die->SetCentroidCorrArr(arr_mean, kTRUE, corrXdim, corrYdim, corrZdim);
//  }
//  else{
//    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
//    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));
//    die->SetCentroidCorrFunction(mean, corrXdim, corrYdim, corrZdim);
//    die->SetWidthCorrFunction(width, corrXdim, corrYdim, corrZdim);
//  }
//
//}


// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
AliDielectronCutGroup* LMEECutLib::GetEventCuts() {
  
  cuts     = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);    
  
  AliDielectronEventCuts* eventCuts = 0x0;
  
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
    
//  cuts->AddCut(eventCuts);
//  cuts->AddCut(AliDielectronVarManager::kQnTPCrpH2,-999.,kTRUE); // makes sure that the event has an eventplane
  return cuts;
}


//Selection of relatively 'flat' centralities
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(AnalysisCut AnaCut) {
  AliDielectronVarCuts* centCuts = 0x0;

  return centCuts;
}



//AliDielectronMixingHandler* LMEECutLib::GetMixingHandler(AnalysisCut AnaCut) {
//  AliDielectronMixingHandler* mixingHandler = 0x0;
//  switch (AnaCut.GetMixing()) {
//    case kEventMixing_1:
//      mixingHandler = new AliDielectronMixingHandler;
//      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10,-5,0,5,10");
//      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,20,30,50,80");
//      mixingHandler->AddVariable(AliDielectronVarManager::kQnTPCrpH2, 6, TMath::Pi()/-2., TMath::Pi()/2.);
//      mixingHandler->SetDepth(40);
//      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
//      break;
//    default: cout << "No Mixer defined" << endl;
//  }
//  return mixingHandler;
//}


  AliDielectronPID* LMEECutLib::GetPIDCutsAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  //nanoAOD Prefilter cuts
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);
  
  //cut set 3 cuts
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.0 , 0. ,100., kTRUE);
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1.0 , 0. ,100., kFALSE);
  pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
  
  
  return   pidCuts;
}

AliAnalysisCuts* LMEECutLib::GetTrackSelectionAna() {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* trackCuts=0x0;
  AliDielectronCutGroup* trackCutsAna= new AliDielectronCutGroup("cgPIDCutsAna","cgPIDCutsAna",AliDielectronCutGroup::kCompAND);
  trackCutsAna->AddCut(GetTrackCuts());
  trackCutsAna->AddCut(GetPIDCutsAna()); 
  trackCuts=trackCutsAna; 
  return trackCuts;
}

AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet) {
  //TRACK CUTS AS USED IN CongifLMEE_nano_PbPb2015.C and in LMEECutLib_caklein.C (cut set 3) (20.03.2018)
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=new AliDielectronCutGroup("trackCuts","trackCuts",AliDielectronCutGroup::kCompAND);
  
  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");

  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
//  varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);

  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);
  
trackCutsAOD->AddCut(AliDielectronVarManager::kPt,           0.2, 8.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
//  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      5.0, 100.0);
//  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   5.0);
//  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    120.0, 160.0);
trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
  
trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);//(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst) -> error
trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//(1<<4) -> error


trackCuts->AddCut(varCuts);
trackCuts->AddCut(trkCuts);

trackCuts->AddCut(trackCutsDiel);
trackCuts->AddCut(trackCutsAOD);

trackCuts->Print();

    
return trackCuts;
}




