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

};



// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
AliDielectronCutGroup* LMEECutLib::GetEventCuts() {
  
  cuts     = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);    
  
  AliDielectronEventCuts* eventCuts = 0x0;
  
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
    
  return cuts;
}


AliDielectronPID* LMEECutLib::GetPIDCutsAna(int sel) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetPIDCutsAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliAnalysisCuts* pidCuts=0x0;
  
  AliDielectronVarCuts *etaRange080 = new AliDielectronVarCuts("etaRange080","etaRange080");
  etaRange080->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
  AliDielectronVarCuts *ptRange200to8000 = new AliDielectronVarCuts("ptRange200to8000","ptRange200to8000");
  ptRange200to8000->AddCut(AliDielectronVarManager::kPt, 0.2, 8.0);
  

  AliDielectronPID *pidFilterCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
 
  //nanoAOD Prefilter cuts - should always be applied  if working on nanoAODs in real data, otherwise MC and real data might not use same cuts
  pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
  pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);
  

  
  switch (sel) {
      case 1:
      // additional PID cuts: carsten analysis PID cut (Physics Forum 12.04.18)
      pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
      pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
      pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
      pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      break;
  }
  
  return(pidFilterCuts);       //Add nanoAODfilter PID cuts

}

AliAnalysisCuts* LMEECutLib::GetTrackSelectionAna(int selTrm, int selPID) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts=0x0;
    
  //Add nanoAOD filter cuts
  AliDielectronVarCuts *varCutsFilter   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCutsFilter = new AliDielectronTrackCuts("TrkCuts","TrkCuts");

  trkCutsFilter->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCutsFilter->SetRequireITSRefit(kTRUE);
  trkCutsFilter->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  varCutsFilter->AddCut(AliDielectronVarManager::kPt,           0.2, 8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
//  varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
  varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);

  varCutsFilter->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCutsFilter->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);
  
  trackCuts->AddCut(varCutsFilter);
  trackCuts->AddCut(trkCutsFilter);
  
  trackCuts->AddCut(GetPIDCutsAna(selPID));
  
  switch (sel) {

    cout<<"chose track cut "<<sel<<endl;      
      
      case 1:

            AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");     
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
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
            trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);   //(1<<4) -> error
            trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::Detector(0),AliDielectronTrackCuts::ITSClusterRequirement(3));  //(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst) -> error

            cgTrackCarstCutsAna->AddCut(trackCutsDiel);
            cgTrackCarstCutsAna->AddCut(trackCutsAOD);
            cgTrackCarstCutsAna->AddCut(SharedClusterCut);
            
            trackCuts = cgTrackCarstCutsAna; 
            break;
    }
  
  return trackCuts;
}






