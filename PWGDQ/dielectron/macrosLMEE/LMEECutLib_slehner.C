#ifndef LMEECutLib_slehner
#define LMEECutLib_slehner

class LMEECutLib {
public:

  LMEECutLib() {}


  AliDielectronPID* GetPIDCutsAna();

  AliDielectronCutGroup* GetTrackCuts(int trsel=0, int pidsel=0);
  AliDielectronEventCuts* GetEventCuts(int sel);
};



// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
AliDielectronEventCuts* LMEECutLib::GetEventCuts(int sel) {
  cout<<"setting event cuts"<<endl;
  
  AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCuts","evcuts");
  
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
    
  return eventCuts;
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
  

  
//  switch (sel) {
      cout<<"chose PID cut "<<sel<<endl; 
//      case 0:
      // additional PID cuts: carsten analysis PID cut (Physics Forum 12.04.18)
      pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
      pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
      pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
      pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
      break;
//  }
  
  return(pidFilterCuts);       //Add nanoAODfilter PID cuts

}

AliDielectronCutGroup* LMEECutLib::GetTrackCuts(int selTr, int selPID) {
  cout << " >>>>>>>>>>>>>>>>>>>>>> GetTrackSelectionAna() >>>>>>>>>>>>>>>>>>>>>> " << endl;
  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("CutsAna","CutsAna",AliDielectronCutGroup::kCompAND);
    
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

//  switch (selTr) {

    std::cout<<"chose track cut "<<selTr<<endl;      
      
//      case 0:

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
//            trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::Detector(0),AliDielectronTrackCuts::ITSClusterRequirement(3));  //(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst) -> error
            trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);

            trackCuts->AddCut(trackCutsDiel);
            trackCuts->AddCut(trackCutsAOD);
            trackCuts->AddCut(SharedClusterCut);
            

            break;
//    }
  
  return trackCuts;
}






