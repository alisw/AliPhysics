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
  AliDielectronCutGroup* GetTrackCuts(int trsel=0, int pidsel=0, Int_t MVACut=0, Bool_t useAODFilterCuts, TString TMVAweight="NONE", Bool_t isUPC=kFALSE);
  AliDielectronCutGroup* GetPairCuts(int pairsel=0);
  AliDielectronEventCuts* GetEventCuts(Double_t centMin, Double_t centMax, Bool_t usePileUpRej=kFALSE, Bool_t isUPC=kFALSE);
  void SetEtaCorrection( Int_t det, Bool_t isMC, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, int sel);
  static AliDielectronPID* pidFilterCuts;
  static TBits *fUsedVars;               // used variables
  TH1 *fPostPIDCntrdCorr;   // post pid correction object for electron sigma centroids in TPC
};

TH3D LMEECutLib::SetEtaCorrection(Int_t det, Bool_t isMC, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, int sel) {
//For usage with TreeMaker
//For efficiency task postcalibration is set in AddTask_slehner_ElectronEfficiency.C
  TString detstr;
  TString type;
  switch(det) {
    case 1: detstr="its";break;
    case 2: detstr="tpc";break;
    case 3: detstr="tof";break;
    } 
  switch(isMC) {
    case kTRUE: type="mc";break;
    case kFALSE: type="data";break;
    }
    
  ::Info("LMEECutLib::SetEtaCorrection",(TString("Starting Correction for ")+detstr).Data());
  TString path="alien:///alice/cern.ch/user/s/selehner/recal/";
  TString fName= "recalib_"+type+"_"+detstr+"_nsigmaele.root";
  
  TFile* corrfile;
  corrfile= TFile::Open(fName.Data());
  if(!corrfile){
    ::Info("LMEECutLib::SetEtaCorrection",(TString("Couldn't find correctiion for ")+detstr).Data()+TString(" -> get it from grid "));
    gSystem->Exec(TString::Format("alien_cp %s .",(path+fName).Data()));
    corrfile = TFile::Open(fName.Data());
    if(!corrfile) ::Error("LMEECutLib::SetEtaCorrection",(TString("Cannot get correction from Alien: ")+path+fName).Data());
    else  ::Info("LMEECutLib::SetEtaCorrection",(TString("Copied correction from Alien: ")+path+fName).Data());
  }    
    
  TH3D* mean = dynamic_cast<TH3D*>(corrfile->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(corrfile->Get("sum_width_correction"));
  
  // AliDielectron::SetCentroidCorrFunction
  UInt_t valType[20] = {0};
  valType[0]=corrXdim;     valType[1]=corrYdim;     valType[2]=corrZdim;
  AliDielectronHistos::StoreVariables(mean, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",corrXdim,corrYdim,corrZdim);
  fPostPIDCntrdCorr = (TH1*)mean->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorr)  {
    printf("POST %s on %s PID CORRECTION added for centroids:  ",detstr.Data(),type.Data());
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
  
  if(fPostPIDCntrdCorr)  {   
    switch(det) {
      case 1: AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorr); break;
      case 2: AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorr); break;
      case 3: AliDielectronPID::SetCentroidCorrFunctionTOF(fPostPIDCntrdCorr); break;
    } 
  }
  
  
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
    printf("POST %s on %s PID CORRECTION added for widths:  ",detstr.Data(),type.Data());
    switch(fPostPIDWdthCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDWdthCorr->GetZaxis()->GetName());break;
    case 2: printf(" %s, ",fPostPIDWdthCorr->GetYaxis()->GetName());break;
    case 1: printf(" %s ",fPostPIDWdthCorr->GetXaxis()->GetName());break;
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
        }
    }
  
  if(fPostPIDWdthCorr){
    switch(det) {
      case 1: AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDCntrdCorr); break;
      case 2: AliDielectronPID::SetWidthCorrFunction(fPostPIDCntrdCorr); break;
      case 3: AliDielectronPID::SetWidthCorrFunctionTOF(fPostPIDCntrdCorr);  break;
  }
  if(sel==1){
        if(mean)   ::Info("LMEECutLib::SetEtaCorrection","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
        else {
        ::Info("LMEECutLib::SetEtaCorrection","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
        return 0;
        }
      return *mean;
  }
  else{
        if(width)   ::Info("LMEECutLib::SetEtaCorrection","Width Correction Histo loaded, entries: %f",width->GetEntries());
        else {
        ::Info("LMEECutLib::SetEtaCorrection","Width Correction Histo not loaded! entries: %f",width->GetEntries());
        return 0;
        }
      return *width;
  }
}
}  

// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Double_t centMin, Double_t centMax, Bool_t usePileUpRej, Bool_t isUPC) {
  ::Info("LMEE_CutLib_slehner","setting event cuts");
  
  AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCuts","evcuts");
  if(!isUPC){
    eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
    eventCuts->SetRequireVertex();
    eventCuts->SetMinVtxContributors(1);
    eventCuts->SetVertexZ(-10.,+10.);
    if(centMax!=0) eventCuts->SetCentralityRange(centMin,centMax,kTRUE);    //isRun2 = true 
  }
  
  if(usePileUpRej){
    Double_t pileUpCutsTPCClustersMin = 2.0e+06;
    Double_t pileUpCutsTPCClustersMax = 3.1e+06;        

    printf("Adding TPC-Cluster based Pile-up rejection! \n");
    TF1* fFitMin = new TF1("fFit","pol6",0,90);
    fFitMin->SetParameters(pileUpCutsTPCClustersMin,-95678.946999,2152.010478,-50.119000,0.780528,-0.006150,0.000019);
    TF1* fFitMax = new TF1("fFit","pol6",0,90);
    fFitMax->SetParameters(pileUpCutsTPCClustersMax,-95678.946999,2152.010478,-50.119000,0.780528,-0.006150,0.000019);
        
    eventCuts->SetMinCorrCutFunction(fFitMin, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);
    eventCuts->SetMaxCorrCutFunction(fFitMax, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);
}   
  
  return eventCuts;
}


AliDielectronPID* LMEECutLib::GetPIDCutsAna(int sel, Bool_t useAODFilterCuts) {

  cout<<"!!!!LMEE_CutLib_slehner: setting PID cuts "<<sel<<endl;
  ::Info("LMEE_CutLib_slehner","setting PID cuts %d",sel);  
  
  //ITS tracklets
  if(sel==-999){
    ::Info("LMEE_CutLib_slehner","setting PID cuts -999 for ITS tracklets");
  }
  else if(sel==-998){
    ::Info("LMEE_CutLib_slehner","setting PID cuts -998 for ITS tracklets");
    pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -1.0, 1. , 0. ,100., kFALSE);
  }
  else{  
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
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;

        case 1:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1., 2.0 , 0. ,100., kFALSE); //<-
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;

        case 2:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.0 , 0. ,100., kTRUE);   //<-
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;

        case 3:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -2.5, 0.5 , 0. ,100., kFALSE);  //<-
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;

        case 4:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -4., 4.0 , 0. ,100., kFALSE);   //<-
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;      

        case 5:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2, 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 5.5 , 0. ,100., kTRUE);   //<-
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;

        case 6:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire); //<-
          break;    

        case 7:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100., 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
  //        pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable); //<-
          break;

        case 8:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
  //        pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          break;

        case 9:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100., 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.0 , 2.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire); //<-
          break;         

        case 10:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100., 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.5 , 3.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire); //<-
          break;  

        //cuts 11-20: Adopt Carsten's Cut Variation (AN 8.2.19) for the variables that are not in the MVA ()
         case 11:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;  

         case 12:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;  

         case 13:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;  

         case 14:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

         case 15:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;  

         case 16:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;  

         case 17:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;   

         case 18:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;   

         case 19:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3. , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

         case 20:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

        //same as above but with TOFreq
         case 21:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;  

         case 22:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;  

         case 23:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;  

         case 24:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1.5, 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

         case 25:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3., 1.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;  

         case 26:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;  

         case 27:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;   

         case 28:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;   

         case 29:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.5 , 3. , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

         case 30:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 1. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3. , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

        case -1:
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.0 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 0.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3.0 , 3.0 , 0. ,100., kFALSE, AliDielectronPID::kRequire);
          break;  

        // loose and tight pid cuts
         case -98:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.5, 1.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -2., 3.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4.5 , 0. ,100., kTRUE);
          break;   

        case -99:
          pidFilterCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron, -3.0, 0.5 , 0. ,100., kFALSE);       
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -1., 2.5 , 0. ,100., kFALSE);
          pidFilterCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -2.0 , 2.5 , 0. ,100., kFALSE, AliDielectronPID::kIfAvailable);
          pidFilterCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0. ,100., kTRUE);
          break;        

    }
  }
  return(pidFilterCuts);       //Add nanoAODfilter PID cuts

}

AliDielectronCutGroup* LMEECutLib::GetPairCuts(int selPair) {  
  
  AliDielectronCutGroup* pairCutsGr = new AliDielectronCutGroup("PairCutsGr","PairCutsGr",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts* pairCuts =new AliDielectronVarCuts("pairCutsAOD","pairCutsAOD"); 
 
    switch (selPair) {     
      case 0:  
        pairCuts->AddCut(AliDielectronVarManager::kQnDeltaPhiTPCrpH2,     -TMath::Pi()/4., TMath::Pi()/4.,kTRUE);   //out of event plane
        pairCuts->AddCut(AliDielectronVarManager::kQnDeltaPhiTPCrpH2,     -TMath::Pi(),-3*TMath::Pi()/4.,kTRUE);    //out of event plane
        pairCuts->AddCut(AliDielectronVarManager::kQnDeltaPhiTPCrpH2,     3*TMath::Pi()/4., TMath::Pi(),kTRUE);     //out of event plane
        break;
      case 1:  
        pairCuts->AddCut(AliDielectronVarManager::kQnDeltaPhiTPCrpH2,     TMath::Pi()/4., 3*TMath::Pi()/4.,kTRUE);    //in event plane
        pairCuts->AddCut(AliDielectronVarManager::kQnDeltaPhiTPCrpH2,     -3*TMath::Pi()/4., -TMath::Pi()/4.,kTRUE);  //in event plane   
        break;
    }
    pairCutsGr->AddCut(pairCuts);
    return pairCutsGr;

}
AliDielectronCutGroup* LMEECutLib::GetTrackCuts(int selTr, int selPID,  Int_t MVACut, Bool_t useAODFilterCuts, TString TMVAweight, Bool_t isUPC) {
  
  
  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("CutsAna","CutsAna",AliDielectronCutGroup::kCompAND);
    
  ////Add nanoAOD filter cuts
  AliDielectronVarCuts *varCutsFilter   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCutsFilter = new AliDielectronTrackCuts("TrkCuts","TrkCuts");

  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD"); 
  
  cout<<"!!!!LMEE_CutLib_slehner: setting Track cuts "<<selTr<<endl;
  ::Info("LMEE_CutLib_slehner","setting Track cuts %d",selTr);
  
  //ITS tracklets
  if(selTr==-999){
    trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kITSonly+AliDielectronTrackCuts::kTPCqual);
    trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
    trackCuts->AddCut(GetPIDCutsAna(-999,kFALSE));     
  }
  else if(selTr==-998){
    trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kITSonly+AliDielectronTrackCuts::kTPCqual);
    trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
    trackCuts->AddCut(GetPIDCutsAna(-998,kFALSE));     
  }
  //normal global tracks
  else{
    trkCutsFilter->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
    trkCutsFilter->SetRequireITSRefit(kTRUE);
    trkCutsFilter->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

    varCutsFilter->AddCut(AliDielectronVarManager::kPt,           0.2, 8.0);
    varCutsFilter->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  //  varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
    varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
    varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  //varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster    // did not work on ESD when filtering nanoAODs
    varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  //  varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
    varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
    varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
    varCutsFilter->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  
    switch (selTr) {     
      case -1:     // Carsten's cuts from Physics Forum 12.04.18, except for conversion rejection cuts (commented out) with frac 0.6-1.05
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //          trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //          trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
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
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -0.3,   0.3);                     //<-
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
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
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.98, 1.02);          //<- 
              break;
        case 3:
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);                    
  //          trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      130.0, 150.0);     //<-
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
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);   
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.90, 1.10);   
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


        //cuts 11-20: Adopt Carsten's Cut Variation (AN 8.2.19) for the variables that are not in the MVA ()
        // these are combinations of kNclsTPC-min: 80 or 100 and kNFclsTPCr-min 80 or 100      
        case 11:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 12:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 13:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 14:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 15:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 16:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 17:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 18:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 19:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;
        case 20:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.95, 1.05);
              break;


        //cuts like 11-20 but with    kNFclsTPCfCross,     0.6, 1.05
        case 21:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 22:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 23:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 24:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 25:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 26:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 27:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 28:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 29:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;
        case 30:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              break;

   // track cut chosen for the results in 70-90% in 0.5-2.7 gev/c^2 (Track cut 30) and  1.1-2.7 gev/c^2 (Track cut 12) with additional event plane cuts  for tracks  
        case 300:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              //event plane cuts 
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     TMath::Pi()/4., 3*TMath::Pi()/4.,kTRUE);    //in event plane
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     -3*TMath::Pi()/4., -TMath::Pi()/4.,kTRUE);  //in event plane   
              break;              
         case 301:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      100.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              //event plane cuts 
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     -TMath::Pi()/4., TMath::Pi()/4.,kTRUE);   //out of event plane
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     -TMath::Pi(),-3*TMath::Pi()/4.,kTRUE);    //out of event plane
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     3*TMath::Pi()/4., TMath::Pi(),kTRUE);     //out of event plane
              break;
        case 220:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              //event plane cuts 
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     TMath::Pi()/4., 3*TMath::Pi()/4.,kTRUE);    //in event plane
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     -3*TMath::Pi()/4., -TMath::Pi()/4.,kTRUE);  //in event plane
              break;
        case 221:    
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      4.0, 100.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   4.0);
  //            trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
              trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.6, 1.05);
              //event plane cuts 
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     -TMath::Pi()/4., TMath::Pi()/4.,kTRUE);   //out of event plane
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     -TMath::Pi(),-3*TMath::Pi()/4.,kTRUE);    //out of event plane
              trackCutsAOD->AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2,     3*TMath::Pi()/4., TMath::Pi(),kTRUE);     //out of event plane       
              break;
    }
  
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);   //(1<<4) -> error
  trackCutsDiel->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
    
  trackCuts->AddCut(GetPIDCutsAna(selPID,useAODFilterCuts)); 
  
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
//  TMVACuts->AddTMVAInput("NCrossedRowsTPC",(Float_t) AliDielectronVarManager::kNclsCrTPC);
//  TMVACuts->AddTMVAInput("NClustersTPC",(Float_t) AliDielectronVarManager::kNclsTPC);
//  TMVACuts->AddTMVAInput("NTPCSignal",(Float_t) AliDielectronVarManager::kTPCsignalN);
  TMVACuts->AddTMVAInput("log(abs(DCAxy))",(Float_t) AliDielectronVarManager::kLogDCAXY);
  TMVACuts->AddTMVAInput("log(abs(DCAz))",(Float_t) AliDielectronVarManager::kLogDCAZ);  
  TMVACuts->AddTMVAInput("chi2GlobalPerNDF",(Float_t) AliDielectronVarManager::kChi2GlobalNDF);
  TMVACuts->AddTMVAInput("chi2ITS",(Float_t) AliDielectronVarManager::kITSchi2);
  TMVACuts->AddTMVAInput("eta",(Float_t) AliDielectronVarManager::kEta);
//  TMVACuts->AddTMVAInput("phi",(Float_t) AliDielectronVarManager::kPhi);
  TMVACuts->AddTMVAInput("pt",(Float_t) AliDielectronVarManager::kPt);  
  TMVACuts->AddTMVAInput("centrality",(Float_t) AliDielectronVarManager::kCentrality);

  TMVACuts->SetTMVAWeights("BDTG method", weightFile.Data());
  Printf("Use TMVA cut value = %f _________",-1. + 0.2*MVACut);
  TMVACuts->SetTMVACutValue(-1. + 0.2*MVACut);
  }

  if(MVACut!=0) trackCuts->AddCut(TMVACuts);

  if(useAODFilterCuts){
    trackCuts->AddCut(varCutsFilter);
    trackCuts->AddCut(trkCutsFilter);
  }
  trackCuts->AddCut(trackCutsDiel);
  if(!isUPC)trackCuts->AddCut(trackCutsAOD);    
  return trackCuts;
}
