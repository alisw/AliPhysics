//======================================================
// Disclaimer: applicable currently only for one cut set
//======================================================


class LMEECutLib {

  public:

  enum LMEECutSet{
    // Basic cuts used for QA stage
    kAllSpecies,
    kElectrons,
    // Cuts used when creating TTrees
    // TTrees used for MVA training (done locally)
    kTTreeCuts,
    // Current primary analysis cut
    kCutSet1,
    kCutSet2, // Remove fITSshared clsuter cut
    kCutSet3, // Lower pt cut to  10 GeV
    // Cut settings to obtain PID correction maps
    kV0_ITScorr,
    kV0_TPCcorr,
    kV0_TOFcorr,
    // Select V0 particles
    kV0_trackCuts,
    // Select electrons in MC via pdg code
    kPdgSel,
    // Loose cuts used in conjuction with kPdgSelf
    kMCsel,
    // Test setting (ignore)
    kV0_allAcc,
    // Loose cuts used to obtain resolution maps
    kResolutionTrackCuts,
    // Linearly increasing cuts over MVA output
    kPIDcut1,
    kPIDcut2,
    kPIDcut3,
    kPIDcut4,
    kPIDcut5,
    kPIDcut6,
    kPIDcut7,
    kPIDcut8,
    kPIDcut9,
    kPIDcut10,
    kPIDcut11,
    kPIDcut12,
    kPIDcut13,
    kPIDcut14,
    kPIDcut15,
    kPIDcut16,
    kPIDcut17,
    kPIDcut18,
    kPIDcut19,
    kPIDcut20,
    // Traditional ePID cut set taken from Run 1 pPb analysis
    kTheoPID,
    kTOFreq,
    // 20 "randomly" chosen cut variations (track+ePID)
    kCutVar1,
    kCutVar2,
    kCutVar3,
    kCutVar4,
    kCutVar5,
    kCutVar6,
    kCutVar7,
    kCutVar8,
    kCutVar9,
    kCutVar10,
    kCutVar11,
    kCutVar12,
    kCutVar13,
    kCutVar14,
    kCutVar15,
    kCutVar16,
    kCutVar17,
    kCutVar18,
    kCutVar19,
    kCutVar20,
    // Cut set used by Sebastian Scheid (also pPb at 5 TeV analysis)
    kScheidCuts,
    // Five different R factor binning schemes
    kMixScheme1,
    kMixScheme2,
    kMixScheme3,
    kMixScheme4,
    kMixScheme5,
    // Cut settings chosen to select regions in eta/phi that show similar
    // acceptance between run1 and run2
    // Select out regions in eta/phi where ITS coverage is "good"
    kGoodEtaPhi,
    // Invert the above selection
    kBadEtaPhi,
    // Vary fITSshared cut (essentially opy of kScheidCuts)
    kITSshared1,
    kITSshared2,
    kITSshared3,
    kITSshared4
  };



  LMEECutLib(Bool_t wSDD): wSDD(wSDD){

    ::Info("LMEECutLib_acapon", "Creating new LMEECutLib");
    fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);

  }

  // Getters
  AliDielectronEventCuts*     GetEventCuts(Bool_t reqAliEventCuts, Bool_t reqAliEventCutsCorrelated);
  AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);

  AliAnalysisCuts*       GetTrackCuts(Int_t cutSet, Int_t PIDcuts);
  AliAnalysisCuts*       GetPIDCuts(Int_t PIDcuts);
  AliDielectronCutGroup* GetPairCuts(Int_t cutSet);
  AliDielectronV0Cuts* GetV0finder();

  // Define signal at MCtruth level
  void SetSignalsMC(AliDielectron* die);

  // PID correction functions used within dielectron framework
  void SetEtaCorrectionTPC(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim);
  void SetEtaCorrectionITS(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC);
  void SetEtaCorrectionTOF(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC);

  // PID correction function used by SimpleTTreeMaker
  // i.e it doesn't need an AliDielectron object
  static TH3D SetEtaCorrectionTPCTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection );
  static TH3D SetEtaCorrectionITSTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC);
  static TH3D SetEtaCorrectionTOFTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC);

  static TBits* fUsedVars;           // Used Variables for correction
  static TH1* fPostPIDCntrdCorrTPC;  // Post PID correction object for electron sigma centroids in TPC
  static TH1* fPostPIDWdthCorrTPC;   // Post PID correction object for electron sigma widths in TPC
  static TH1* fPostPIDCntrdCorrITS;  // Post PID correction object for electron sigma centroids in ITS
  static TH1* fPostPIDWdthCorrITS;   // Post PID correction object for electron sigma widths in ITS
  static TH1* fPostPIDCntrdCorrTOF;  // Post PID correction object for electron sigma centroids in TOF
  static TH1* fPostPIDWdthCorrTOF;   //Post PID correction object for electron sigma widths in TOF

  private:
    Bool_t wSDD;
};

// Eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
void LMEECutLib::SetEtaCorrectionTPC(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim){

  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
  TString fileName;
  if(wSDD == kTRUE){
    fileName = "outputTPC.root";
  }else{
    fileName = "outputTPC_woSDD.root";
  }

  TFile* inFile = TFile::Open(localPath+fileName, "READ");
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
    std::cout << "Copy TPC correction from Alien" << std::endl;
    inFile = TFile::Open(fileName);
  }
  else {
    std::cout << "Correction loaded locally" << std::endl;
  }
  TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
  die->SetCentroidCorrFunction(mean, corrXdim, corrYdim, corrZdim);
  die->SetWidthCorrFunction(width, corrXdim, corrYdim, corrZdim);

}

// Eta correction for the centroid and width of electron sigmas in the ITS, can be one/two/three-dimensional
void LMEECutLib::SetEtaCorrectionITS(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC){

  if(!wSDD){
    return;
  }
  std::cout << "starting LMEECutLib::SetEtaCorrectionITS()\n";
  TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
  TString fileName = "outputITS";
  if(hasMC){
    fileName.Append("_MC.root");
  }else{
    fileName.Append(".root");
  }

  TFile* inFile = TFile::Open(localPath+fileName);
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
    std::cout << "Copy ITS correction from Alien" << std::endl;
    inFile = TFile::Open(fileName);
  }
  else{
    std::cout << "Correction loaded locally" << std::endl;
  }

  TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
  die->SetCentroidCorrFunctionITS(mean, corrXdim, corrYdim, corrZdim);
  die->SetWidthCorrFunctionITS(width, corrXdim, corrYdim, corrZdim);

}

// Eta correction for the centroid and width of electron sigmas in the TOF, can be one/two/three-dimensional
void LMEECutLib::SetEtaCorrectionTOF(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC){

  std::cout << "starting LMEECutLib::SetEtaCorrectionTOF()\n";
  TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
  TString fileName = "outputTOF";
  if(hasMC){
    if(wSDD == kTRUE){
      fileName.Append("_MC.root");
    }else{
      fileName.Append("_woSDD_MC.root");
    }
  }else{
    if(wSDD == kTRUE){
      fileName.Append(".root");
    }else{
      fileName.Append("_woSDD.root");
    }
  }

  TFile* inFile = TFile::Open(localPath+fileName);
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    inFile = TFile::Open(fileName);
  }
  else{
    std::cout << "Correction loaded localy" << std::endl;
  }

  TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
  TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
  die->SetCentroidCorrFunctionTOF(mean, corrXdim, corrYdim, corrZdim);
  die->SetWidthCorrFunctionTOF(width, corrXdim, corrYdim, corrZdim);

}


// Eta correction functions to be used with the SimpleTreeMakerClass
// Only differs from those above in that they don't require an AliDielectron object
// Eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
static TH3D LMEECutLib::SetEtaCorrectionTPCTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection){

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionTPC() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
  TString fileName;
  if(wSDD == kTRUE){
    fileName = "outputTPC.root";
  }else{
    fileName = "outputTPC_woSDD.root";
  }

  TFile* inFile = TFile::Open(localPath+fileName, "READ");
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
    std::cout << "Copy TPC correction from Alien" << std::endl;
    inFile = TFile::Open(fileName);
  }
  else {
    std::cout << "Correction loaded locally" << std::endl;
  }
  TH3D* mean  = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
  TH3D* width = dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
  if(!mean || !width){
    ::Error("LMEECutLib_acapon", "Recal histograms not found.");
    return 0x0;
  }

  // AliDielectron::SetCentroidCorrFunction
  UInt_t valType[20] = {0};
  valType[0] = corrXdim;     valType[1] = corrYdim;     valType[2] = corrZdim;

  AliDielectronHistos::StoreVariables(mean, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d", corrXdim, corrYdim, corrZdim);
  printf("%s", key);
  fPostPIDCntrdCorrTPC = (TH1*)mean->Clone(key.Data());
  if(!fPostPIDCntrdCorrTPC){
    ::Error("LMEECutLib_acapon", "CentroidCorr TH1 not cloned");
    return 0x0;
  }

  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorrTPC){
    printf("POST TPC PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTPC->GetDimension()){
      case 1: printf(" %s ",fPostPIDCntrdCorrTPC->GetXaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTPC->GetYaxis()->GetName());
      case 3: printf(" %s, ",fPostPIDCntrdCorrTPC->GetZaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }

  if(fPostPIDCntrdCorrTPC){
    AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorrTPC);
  }

  // AliDielectron::SetWidthCorrFunction
  {
  UInt_t valType[20] = {0};
  valType[0]=corrXdim;     valType[1]=corrYdim;     valType[2]=corrZdim;
  AliDielectronHistos::StoreVariables(width, valType);

  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",corrXdim,corrYdim,corrZdim);
  fPostPIDWdthCorrTPC = (TH1*)width->Clone(key.Data());

  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorrTPC)  {
    printf("POST TPC PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTPC->GetDimension()) {
      case 1: printf(" %s ",fPostPIDWdthCorrTPC->GetXaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTPC->GetYaxis()->GetName());
      case 3: printf(" %s, ",fPostPIDWdthCorrTPC->GetZaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }
  }

  if(fPostPIDWdthCorrTPC){
    AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorrTPC);
  }

  if(selection == 1){
    if(mean){
      ::Info("LMEECutLib::SetEtaCorrectionTPC","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
    }else{
      ::Info("LMEECutLib::SetEtaCorrectionTPC","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
      return 0;
    }
  return *mean;
  }
  else{
    if(width){
      ::Info("LMEECutLib::SetEtaCorrectionTPC","Width Correction Histo loaded, entries: %f",width->GetEntries());
    }else {
      ::Info("LMEECutLib::SetEtaCorrectionTPC","Width Correction Histo not loaded! entries: %f",width->GetEntries());
      return 0;
    }
    return *width;
  }

}

// Eta correction for the centroid and width of electron sigmas in the ITS, can be one/two/three-dimensional
TH3D LMEECutLib::SetEtaCorrectionITSTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC){

  if(!wSDD){
    return;
  }
  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionITSTTree() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionITSTTree()\n";
  TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
  TString fileName = "outputITS";
  if(hasMC){
    fileName.Append("_MC.root");
  }else{
    fileName.Append(".root");
  }

  TFile* inFile = TFile::Open(localPath+fileName);
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
    std::cout << "Copy ITS correction from Alien" << std::endl;
    inFile = TFile::Open(fileName);
  }
  else {
    std::cout << "Correction loaded locally" << std::endl;
  }
  TH3D* mean  = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
  TH3D* width = dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
  if(!mean || !width){
    ::Error("LMEECutLib_acapon", "Recal histograms not found.");
    return 0x0;
  }

  // AliDielectron::SetCentroidCorrFunction
  UInt_t valType[20] = {0};
  valType[0] = corrXdim;     valType[1] = corrYdim;     valType[2] = corrZdim;

  AliDielectronHistos::StoreVariables(mean, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d", corrXdim, corrYdim, corrZdim);
  printf("%s", key);
  fPostPIDCntrdCorrITS = (TH1*)mean->Clone(key.Data());
  if(!fPostPIDCntrdCorrITS){
    ::Error("LMEECutLib_acapon", "CentroidCorr TH1 not cloned");
    return 0x0;
  }

  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorrITS){
    printf("POST ITS PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrITS->GetDimension()){
      case 1: printf(" %s ",fPostPIDCntrdCorrITS->GetXaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrITS->GetYaxis()->GetName());
      case 3: printf(" %s, ",fPostPIDCntrdCorrITS->GetZaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }

  if(fPostPIDCntrdCorrITS){
    AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorrITS);
  }

  // AliDielectron::SetWidthCorrFunction
  {
  UInt_t valType[20] = {0};
  valType[0]=corrXdim;     valType[1]=corrYdim;     valType[2]=corrZdim;
  AliDielectronHistos::StoreVariables(width, valType);

  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",corrXdim,corrYdim,corrZdim);
  fPostPIDWdthCorrITS = (TH1*)width->Clone(key.Data());

  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorrITS)  {
    printf("POST ITS PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrITS->GetDimension()) {
      case 1: printf(" %s ",fPostPIDWdthCorrITS->GetXaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrITS->GetYaxis()->GetName());
      case 3: printf(" %s, ",fPostPIDWdthCorrITS->GetZaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }
  }

  if(fPostPIDWdthCorrITS){
    AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDWdthCorrITS);
  }

  if(selection == 1){
    if(mean){
      ::Info("LMEECutLib::SetEtaCorrectionITSTTree","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
    }else{
      ::Info("LMEECutLib::SetEtaCorrectionITSTTree","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
      return 0;
    }
  return *mean;
  }
  else{
    if(width){
      ::Info("LMEECutLib::SetEtaCorrectionITSTTree","Width Correction Histo loaded, entries: %f",width->GetEntries());
    }else {
      ::Info("LMEECutLib::SetEtaCorrectionITSTTree","Width Correction Histo not loaded! entries: %f",width->GetEntries());
      return 0;
    }
    return *width;
  }

}

// Eta correction for the centroid and width of electron sigmas in the TOF, can be one/two/three-dimensional
static TH3D LMEECutLib::SetEtaCorrectionTOFTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC){

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionTOFTTree() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionTOFTTree()\n";
  TString localPath = "/home/aaron/Data/diElec_framework_output/PIDcalibration/";
  TString fileName = "outputTOF";
  if(hasMC){
    if(wSDD == kTRUE){
      fileName.Append("_MC.root");
    }else{
      fileName.Append("_woSDD_MC.root");
    }
  }else{
    if(wSDD == kTRUE){
      fileName.Append(".root");
    }else{
      fileName.Append("_woSDD.root");
    }
  }

  TFile* inFile = TFile::Open(localPath+fileName);
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/"+fileName+" .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    inFile = TFile::Open(fileName);
  }
  else {
    std::cout << "Correction loaded locally" << std::endl;
  }
  TH3D* mean  = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
  TH3D* width = dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
  if(!mean || !width){
    ::Error("LMEECutLib_acapon", "Recal histograms not found.");
    return 0x0;
  }

  // AliDielectron::SetCentroidCorrFunction
  UInt_t valType[20] = {0};
  valType[0] = corrXdim;     valType[1] = corrYdim;     valType[2] = corrZdim;

  AliDielectronHistos::StoreVariables(mean, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d", corrXdim, corrYdim, corrZdim);
  printf("%s", key);
  fPostPIDCntrdCorrTOF = (TH1*)mean->Clone(key.Data());
  if(!fPostPIDCntrdCorrTOF){
    ::Error("LMEECutLib_acapon", "CentroidCorr TH1 not cloned");
    return 0x0;
  }

  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorrTOF){
    printf("POST TOF PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTOF->GetDimension()){
      case 1: printf(" %s ",fPostPIDCntrdCorrTOF->GetXaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTOF->GetYaxis()->GetName());
      case 3: printf(" %s, ",fPostPIDCntrdCorrTOF->GetZaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }

  if(fPostPIDCntrdCorrTOF){
    AliDielectronPID::SetCentroidCorrFunctionTOF(fPostPIDCntrdCorrTOF);
  }

  // AliDielectron::SetWidthCorrFunction
  {
  UInt_t valType[20] = {0};
  valType[0]=corrXdim;     valType[1]=corrYdim;     valType[2]=corrZdim;
  AliDielectronHistos::StoreVariables(width, valType);

  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",corrXdim,corrYdim,corrZdim);
  fPostPIDWdthCorrTOF = (TH1*)width->Clone(key.Data());

  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorrTOF)  {
    printf("POST TOF PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTOF->GetDimension()) {
      case 1: printf(" %s ",fPostPIDWdthCorrTOF->GetXaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTOF->GetYaxis()->GetName());
      case 3: printf(" %s, ",fPostPIDWdthCorrTOF->GetZaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(corrXdim, kTRUE);
    fUsedVars->SetBitNumber(corrYdim, kTRUE);
    fUsedVars->SetBitNumber(corrZdim, kTRUE);
  }
  }

  if(fPostPIDWdthCorrTOF){
    AliDielectronPID::SetWidthCorrFunctionTOF(fPostPIDWdthCorrTOF);
  }

  if(selection == 1){
    if(mean){
      ::Info("LMEECutLib::SetEtaCorrectionTOFTTree","Mean Correction Histo loaded, entries: %f",mean->GetEntries());
    }else{
      ::Info("LMEECutLib::SetEtaCorrectionTOFTTree","Mean Correction Histo not loaded! entries: %f",mean->GetEntries());
      return 0;
    }
  return *mean;
  }
  else{
    if(width){
      ::Info("LMEECutLib::SetEtaCorrectionTOFTTree","Width Correction Histo loaded, entries: %f",width->GetEntries());
    }else {
      ::Info("LMEECutLib::SetEtaCorrectionTOFTTree","Width Correction Histo not loaded! entries: %f",width->GetEntries());
      return 0;
    }
    return *width;
  }

}

// Note: event cuts are identical for all analysis 'cutDefinition's that run together!
// the selection is hardcoded in the AddTask
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Bool_t reqAliEventCuts, Bool_t reqAliEventCutsCorrelated){

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetEventCuts() >>>>>>>>>>>>>>>>>>>>>> ");
  AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCuts_acapon","Vertex Track && |vtxZ|<10 && ncontrib>0");

  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);

  // Possibility to use AliEventCuts class for additional pile up rejection
  if(reqAliEventCuts){
    std::cout << "------ Using AliEventCuts class -------" << std::endl;
    eventCuts->SetRequireAliEventCuts(reqAliEventCuts, reqAliEventCutsCorrelated);
  }

  eventCuts->Print();
  return eventCuts;
}

AliDielectronMixingHandler* LMEECutLib::GetMixingHandler(Int_t cutSet){
  AliDielectronMixingHandler* mixingHandler = 0x0;
  switch (cutSet) {
    case kAllSpecies:
    case kElectrons:
    case kCutSet1:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 10, 20, 30, 40, 60, 80,100");
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kMixScheme1:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 5, 10, 20, 30, 40, 60, 80,100");
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kMixScheme2:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -5., 0., 5., 10.");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 10, 20, 30, 40, 60, 80,100");
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kMixScheme3:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -5., 0., 5., 10.");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 5, 10, 20, 30, 40, 60, 80,100");
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kMixScheme4:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., 0., 10.");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 10, 20, 30, 40, 60, 80,100");
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    case kMixScheme5:
      mixingHandler = new AliDielectronMixingHandler;
      mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
      mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 20, 40, 60, 80,100");
      mixingHandler->SetDepth(60);
      mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
      break;
    default:
      std::cout << "No Mixer defined" << std::endl;
  }
  return mixingHandler;
}



// Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
AliDielectronCutGroup* LMEECutLib::GetPairCuts(Int_t cutSet)  {

  ::Info("LMEECutLibg_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetPairCuts() >>>>>>>>>>>>>>>>>>>>>> ");
  // Final OR cut group to incorporate the following cuts (below)
  AliDielectronCutGroup* allCuts    = new AliDielectronCutGroup("allCuts", "allCuts", AliDielectronCutGroup::kCompOR);

  // AND cut group to select low mass pairs with large opening angle
  AliDielectronCutGroup* convRejCut = new AliDielectronCutGroup("convRejCut", "convRejCut", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts* convMassCut = new AliDielectronVarCuts("convMassCut", "convMassCut");
  AliDielectronVarCuts* convPhiVCut = new AliDielectronVarCuts("convPhiVCut", "convPhiVCut");
  convMassCut->AddCut(AliDielectronVarManager::kM, 0.00, 0.14);
  convPhiVCut->AddCut(AliDielectronVarManager::kPhivPair, 0., 2.);
  convRejCut->AddCut(convMassCut);
  convRejCut->AddCut(convPhiVCut);

  // Mass cut to include any pairs with mass greater than 0.1 GeV
  AliDielectronVarCuts* pairMassCut = new AliDielectronVarCuts("pairMassCut", "pairMassCut");
  pairMassCut->AddCut(AliDielectronVarManager::kM, 0.14, 5.0);

  allCuts->AddCut(convRejCut);
  allCuts->AddCut(pairMassCut);
  allCuts->Print();

  return allCuts;

}


AliAnalysisCuts* LMEECutLib::GetPIDCuts(Int_t PIDcuts) {

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetPIDCuts() >>>>>>>>>>>>>>>>>>>>>> ");

  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts", "cuts", AliDielectronCutGroup::kCompAND);
  AliDielectronPID* cutsPID   = new AliDielectronPID("PID", "PID");
  //-----------------------------------------------
  // PID cuts depend on TPC_inner_p, if not specified
  // PID cut ranges correspond to global momentum P
  // check it again!!!
  //-----------------------------------------------

   switch(PIDcuts){
    case kScheidCuts:
      // "Hadron rejection band" PID scheme
      // PID with the TPC as per usual, then recover electron passing second
      // criteria
      AliDielectronCutGroup* hadBandRej = new AliDielectronCutGroup("hadBandRejs", "hadBandRej", AliDielectronCutGroup::kCompOR);
      AliDielectronPID* cutsTPC   = new AliDielectronPID("cutsTCP", "cutsTCP");
      cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3., 3., 0., 100., kFALSE);
      cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kPion, -100., 3.5, 0., 100., kTRUE);
      cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kKaon, -3., 3., 0.2, 100., kTRUE);
      cutsTPC->AddCut(AliDielectronPID::kTPC, AliPID::kProton, -3., 3., 0.2, 100., kTRUE);
      AliDielectronPID* recoverTOF = new AliDielectronPID("recoverTOF", "recoverTOF");
      recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3., 3., 0., 100., kFALSE);
      recoverTOF->AddCut(AliDielectronPID::kTPC, AliPID::kPion, -100., 4, 0., 100., kTRUE);
      recoverTOF->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3., 3., 0., 100., kFALSE, AliDielectronPID::kRequire);
      hadBandRej->AddCut(cutsTPC);
      hadBandRej->AddCut(recoverTOF);
      return hadBandRej;
    default:
      std::cout << "No Analysis PID Cut defined " << std::endl;
      return 0x0;
    }
    std::cout << "No Analysis PID Cut defined " << std::endl;
    return 0x0;
}

// Make/Tighten track Cuts that are *NOT* already //done in the AOD production
// **IMPORTANT**: For AODs, select FilterBit. The method is ignored for ESDs
AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet, Int_t PIDcuts){

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> ");
  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("trackCuts", "trackCuts", AliDielectronCutGroup::kCompAND);

  // Track cut groups
  AliDielectronVarCuts* varCutsFilter     = new AliDielectronVarCuts("varCutsFilter", "varCutsFilter");
  AliDielectronTrackCuts* trackCutsFilter = new AliDielectronTrackCuts("trackCutsFilter", "trackCutsFilter");

  // V0 track cut groups
  AliDielectronV0Cuts* gammaV0cuts = new AliDielectronV0Cuts("gammaV0cuts", "gammaV0cuts");
  AliDielectronVarCuts* trackCutsV0 = new AliDielectronVarCuts("trackCutsV0", "trackCutsV0");


  switch(cutSet){
    //----------
    // these MAIN settings just load the main track selection directly below:
    //----------
    case kScheidCuts:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);

      //trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      AliESDtrackCuts* esdTrackCuts = 0x0;
      esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
      esdTrackCuts->SetMaxDCAToVertexXY(2.4);
      esdTrackCuts->SetMaxDCAToVertexZ(3.2);
      esdTrackCuts->SetDCAToVertex2D(kTRUE);
      
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,   -1., 1.);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,    -3., 3.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100., 160.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80., 160.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.4, 1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,        3. , 10.);
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0, 5.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,       1.0, 6.0, kTRUE);

      trackCuts->AddCut(esdTrackCuts);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    default:
      std::cout << "No Analysis Track Cut defined" << std::endl;
    }
    std::cout << "Track cuts not applied...." << std::endl;
    return 0x0;
}

AliDielectronV0Cuts* LMEECutLib::GetV0finder(){

  // V0 finder to exclude conversions
  AliDielectronV0Cuts* rejConversions = new AliDielectronV0Cuts("IsGamma", "IsGamma");
  std::cout << "Adding V0 conversion cut!" << std::endl;
  // which V0 finder you want to use
  rejConversions->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
  // add some pdg codes (they are used then by the KF package and important for gamma conversions)
  rejConversions->SetPdgCodes(22,11,11); // mother, daughter1 and 2
  // add default PID cuts (defined in AliDielectronPID)
  // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
  rejConversions->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
  // add the pair cuts for V0 candidates
  rejConversions->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.00, kFALSE);
  rejConversions->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.00, kFALSE);
  rejConversions->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  rejConversions->AddCut(AliDielectronVarManager::kR,                             3.0,  90.00, kFALSE);
  rejConversions->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  rejConversions->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
  rejConversions->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  // selection or rejection of V0 tracks
  rejConversions->SetExcludeTracks(kTRUE);

  return rejConversions;
}

//______________________________________________________________________________________
//----------------------------- Define MC Signals --------------------------------------
void LMEECutLib::SetSignalsMC(AliDielectron* die){

  // Dielectrons originating from pion dalitz decays
  AliDielectronSignalMC* PiDalitz = new AliDielectronSignalMC("Pi0","di-electrons from Pi0 dalitz");
  PiDalitz->SetLegPDGs(11,-11);
  PiDalitz->SetMotherPDGs(111,111);
  PiDalitz->SetMothersRelation(AliDielectronSignalMC::kSame);
  PiDalitz->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PiDalitz->SetCheckBothChargesLegs(kTRUE,kTRUE);
  PiDalitz->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(PiDalitz);

  // Dielectron pairs from same mother (excluding conversions)
  AliDielectronSignalMC* pair_sameMother = new AliDielectronSignalMC("sameMother","sameMother");
  pair_sameMother->SetLegPDGs(11,-11);
  pair_sameMother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_sameMother->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  pair_sameMother->SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameMother->SetMotherPDGs(22,22,kTRUE,kTRUE); // Exclude conversion
  die->AddSignalMC(pair_sameMother);
  // Used pdg codes (defined in AliDielectronMC::ComparePDG)
  // 401: open charm meson
  // 404: charged open charmed mesons NO s quark
  // 405: neutral open charmed mesons
  // 406: charged open charmed mesons with s quark
  // 501: open beauty mesons
  // 503: all beauty hadrons
  // 504: charged open beauty mesons NO s quark
  // 505: neutral open beauty mesons
  // 506: charged open beauty mesons with s quark
  // all D mesons

  // decay channels
  // (1) D -> e X
  // (1) B -> e X
  // (2) B -> D X -> e X Y
  // (3) B -> e D X -> ee X Y always produces ULS pair

  // Electrons from open beauty mesons and baryons
  AliDielectronSignalMC* eleFinalStateFromB = new AliDielectronSignalMC("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB->SetLegPDGs(11,-11);
  eleFinalStateFromB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB->SetMotherPDGs(502, 502);
  eleFinalStateFromB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  eleFinalStateFromB->SetCheckCorrelatedHF(kTRUE);
  die->AddSignalMC(eleFinalStateFromB);

  // Electrons from open charm mesons and baryons
  AliDielectronSignalMC* eleFinalStateFromD = new AliDielectronSignalMC("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD->SetLegPDGs(11,-11);
  eleFinalStateFromD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD->SetMotherPDGs(402, 402);
  eleFinalStateFromD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  eleFinalStateFromD->SetCheckCorrelatedHF(kTRUE);
  die->AddSignalMC(eleFinalStateFromD);

  AliDielectronSignalMC* eleFromJPsi = new AliDielectronSignalMC("eleFromJPsi", "eleFromJPsi");
  eleFromJPsi->SetLegPDGs(11,-11);
  eleFromJPsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFromJPsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFromJPsi->SetMotherPDGs(443, 443);
  eleFromJPsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  eleFromJPsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(eleFromJPsi);

}
