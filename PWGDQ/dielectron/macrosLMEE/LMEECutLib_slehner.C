#ifndef LMEECutLib_slehner
#define LMEECutLib_slehner

class LMEECutLib {
    
  
private:
      
    
public:
  
  LMEECutLib() {  
  ::Info("LMEECutLib","CREATE NEW LMEECUTLIB slehner");
    pidFilterCuts = new AliDielectronPID("PIDCuts1","PIDCuts1");
    fUsedVars= new TBits(AliDielectronVarManager::kNMaxValues);
    
  }
  static AliDielectronPID* GetPIDCutsAna();
  AliDielectronCutGroup* GetTrackCuts(int trsel=0, int pidsel=0, Int_t MVACut=0, Bool_t useAODFilterCuts, TString TMVAweight="TMVAClassification_BDTG.weights_094.xml");
  AliDielectronEventCuts* GetEventCuts(int sel);
  static TH3D SetEtaCorrectionTPC( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise);
  static AliDielectronPID* pidFilterCuts;
  static TBits *fUsedVars;               // used variables
  TH1 *fPostPIDCntrdCorr;   // post pid correction object for electron sigma centroids in TPC
};

static TH3D LMEECutLib::SetEtaCorrectionTPC( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise, int sel) {
//For usage with TreeMaker
//For efficiency task postcalibration is set in AddTask_slehner_ElectronEfficiency.C
    
  ::Info("LMEECutLib::SetEtaCorrectionTPC","starting LMEECutLib::SetEtaCorrectionTPC()\n");
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/recalib_data_tpc_nsigmaele.root";
  gSystem->Exec(TString::Format("alien_cp %s .",path.Data()));
  ::Info("LMEECutLib::SetEtaCorrectionTPC","Copy TPC correction from Alien: %s",path.Data());
  _file = TFile::Open("recalib_data_tpc_nsigmaele.root");
  
  TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));
  
  // AliDielectron::SetCentroidCorrFunction
  UInt_t valType[20] = {0};
  valType[0]=corrXdim;     valType[1]=corrYdim;     valType[2]=corrZdim;
  AliDielectronHistos::StoreVariables(mean, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",corrXdim,corrYdim,corrZdim);
  fPostPIDCntrdCorr = (TH1*)mean->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorr)  {
    printf("POST TPC PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDCntrdCorr->GetZaxis()->GetName());
    case 2: printf(" %s, ",fPostPIDCntrdCorr->GetYaxis()->GetName());
    case 1: printf(" %s ",fPostPIDCntrdCorr->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }
  
  if(fPostPIDCntrdCorr)     AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorr);
  
  
  
  // AliDielectron::SetWidthCorrFunction
  {
  UInt_t valType[20] = {0};
  valType[0]=corrXdim;     valType[1]=corrYdim;     valType[2]=corrZdim;
  AliDielectronHistos::StoreVariables(width, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",corrXdim,corrYdim,corrZdim);
  fPostPIDWdthCorr = (TH1*)width->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorr)  {
    printf("POST TPC PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDWdthCorr->GetZaxis()->GetName());
    case 2: printf(" %s, ",fPostPIDWdthCorr->GetYaxis()->GetName());
    case 1: printf(" %s ",fPostPIDWdthCorr->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
        }
    }
  
  if(fPostPIDWdthCorr)      AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorr);

  if(sel==1){
        if(mean)   ::Info("LMEECutLib::SetEtaCorrectionTPC","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
        else {
        ::Info("LMEECutLib::SetEtaCorrectionTPC","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
        return 0;
        }
      return *mean;
  }
  else{
        if(width)   ::Info("LMEECutLib::SetEtaCorrectionTPC","Width Correction Histo loaded, entries: %f",width->GetEntries());
        else {
        ::Info("LMEECutLib::SetEtaCorrectionTPC","Width Correction Histo not loaded! entries: %f",width->GetEntries());
        return 0;
        }
      return *width;
  }
}  

// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Double_t centMin, Double_t centMax) {
  ::Info("LMEE_CutLib_slehner","setting event cuts");
  
  AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCutsSlehner","evcuts");
  
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  eventCuts->SetCentralityRange(centMin,centMax,kTRUE);    //isRun2 = true 
  return eventCuts;
}


AliDielectronPID* LMEECutLib::GetPIDCutsAna(int sel, Bool_t useAODFilterCuts) {
    
  ::Info("LMEE_CutLib_slehner","setting PID cuts");
//  pidFilterCuts = new AliDielectronPID("PIDCuts1","PIDCuts1");
  //nanoAOD Prefilter cuts - should always be applied  if working on nanoAODs in real data, otherwise MC and real data might not use same cuts
  if(useAODFilterCuts){
    pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
    pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
    pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);
  }
  
  switch (sel) {

      case 0:
        // additional PID cuts: carsten analysis PID cut (Physics Forum 12.04.18)
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
        pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
        break;
    
      case 1:
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1, 2.0 , 0. ,100., kFALSE); //<-
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
        pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
        break;
    
      case 2:
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);   //<-
        pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
        break;
    
      case 3:
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
        pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 0.5 , 0. ,100., kFALSE);  //<-
        pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
        break;
    
      case 4:
        // additional PID cuts: carsten analysis PID cut (Physics Forum 12.04.18)
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -4, 4.0 , 0. ,100., kFALSE);   //<-
        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
        pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
        pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
        break;      

            
  }
  
  return(pidFilterCuts);       //Add nanoAODfilter PID cuts

}

AliDielectronCutGroup* LMEECutLib::GetTrackCuts(int selTr, int selPID,  Int_t MVACut, Bool_t useAODFilterCuts, TString TMVAweight) {
  
  ::Info("LMEE_CutLib_slehner","setting Track cuts");
  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup(TString::Format("CutTr%d_PID%d_MVA%f",selTr, selPID,MVACut),TString::Format("CutTr%d_PID%d_MVA%f",selTr, selPID,MVACut),AliDielectronCutGroup::kCompAND);
    
  ////Add nanoAOD filter cuts
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
//varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
  varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCutsFilter->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");   
  switch (selTr) {     
    case 0:     // Carsten's cuts from Physics Forum 12.04.18, except for conversion rejection cuts (commented out)
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
            break;
     case 1:
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -0.5,   0.5);                     //<-
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -1.0,   1.0);                     //<-
//          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
            break;
      case 2:
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);                    
//          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);            
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.97, 1.03);          //<- 
            break;
      case 3:
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);                    
//          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      130.0, 160.0);     //<-
//          trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);            
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);                   
            break;
      case 4:
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
//          trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);   //<-
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);   //<-
            break;
      case 5:     // Carsten's cuts from Physics Forum 12.04.18, with conversion rejection cuts
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
            trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
            break;
  }
  
  // /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TMVA vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

//  TString weightFile="alien:///alice/cern.ch/user/s/slehner/TMVAClassification_BDTG.weights_094.xml";
  AliDielectronTMVACuts *TMVACuts=0;
  
  if(MVACut!=0){
    TString weightFile="alien:///alice/cern.ch/user/s/selehner/TMVAweights/"+TMVAweight;

    Printf("Use TMVA weight input file: %s",weightFile.Data());

    TMVACuts = new AliDielectronTMVACuts(TString::Format("TMVA%d",MVACut),TString::Format("TMVA%d",MVACut));
    TMVACuts->AddTMVAInput("nITS", (Float_t) AliDielectronVarManager::kNclsITS);
    TMVACuts->AddTMVAInput("ITS1Shared", (Float_t) AliDielectronVarManager::kClsS1ITS);
    TMVACuts->AddTMVAInput("ITS2Shared", (Float_t) AliDielectronVarManager::kClsS2ITS);
    TMVACuts->AddTMVAInput("ITS3Shared", (Float_t) AliDielectronVarManager::kClsS3ITS);
    TMVACuts->AddTMVAInput("ITS4Shared", (Float_t) AliDielectronVarManager::kClsS4ITS);
    TMVACuts->AddTMVAInput("ITS5Shared", (Float_t) AliDielectronVarManager::kClsS5ITS);
    TMVACuts->AddTMVAInput("ITS6Shared", (Float_t) AliDielectronVarManager::kClsS6ITS);
    TMVACuts->AddTMVAInput("nITSshared_frac",(Float_t) AliDielectronVarManager::kNclsSFracITS);
    TMVACuts->AddTMVAInput("NCrossedRowsTPC",(Float_t) AliDielectronVarManager::kNclsCrTPC);
    TMVACuts->AddTMVAInput("NClustersTPC",(Float_t) AliDielectronVarManager::kNclsTPC);
    TMVACuts->AddTMVAInput("NTPCSignal",(Float_t) AliDielectronVarManager::kTPCsignalN);
    TMVACuts->AddTMVAInput("log(abs(DCAxy))",(Float_t) AliDielectronVarManager::kLogDCAXY);
    TMVACuts->AddTMVAInput("log(abs(DCAz))",(Float_t) AliDielectronVarManager::kLogDCAZ);  
    TMVACuts->AddTMVAInput("chi2GlobalPerNDF",(Float_t) AliDielectronVarManager::kChi2GlobalNDF);
    TMVACuts->AddTMVAInput("chi2ITS",(Float_t) AliDielectronVarManager::kITSchi2);
    TMVACuts->AddTMVAInput("eta",(Float_t) AliDielectronVarManager::kEta);
    TMVACuts->AddTMVAInput("phi",(Float_t) AliDielectronVarManager::kPhi);
    TMVACuts->AddTMVAInput("pt",(Float_t) AliDielectronVarManager::kPt);  
    TMVACuts->AddTMVAInput("centrality",(Float_t) AliDielectronVarManager::kCentrality);

    TMVACuts->SetTMVAWeights("BDTG method", weightFile.Data());
    Printf("Use TMVA cut value = %f _________",-1. + 0.2*MVACut);
    TMVACuts->SetTMVACutValue(-1. + 0.2*MVACut);
  }
  
  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);   //(1<<4) -> error
  trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
  
  if(useAODFilterCuts){
    trackCuts->AddCut(varCutsFilter);
    trackCuts->AddCut(trkCutsFilter);
  }
  trackCuts->AddCut(trackCutsDiel);
  trackCuts->AddCut(trackCutsAOD);
  if(MVACut!=0) trackCuts->AddCut(TMVACuts);

  trackCuts->AddCut(GetPIDCutsAna(selPID,useAODFilterCuts));
  return trackCuts;
}