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
  TString localPath = "~/Data/diElec_outputs/PIDcalibration/";
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
  TString localPath = "~/Data/diElec_outputs/PIDcalibration/";
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
  TString localPath = "~/Data/diElec_outputs/PIDcalibration/";
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
  TString localPath = "~/Data/diElec_outputs/PIDcalibration/";
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
  TString localPath = "~/Data/diElec_outputs/PIDcalibration/";
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
  TString localPath = "~/Data/diElec_outputs/PIDcalibration/";
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

  // Pdg code selection of electrons
  AliDielectronVarCuts* PdgElectron = new AliDielectronVarCuts("PdgElectron","PdgElectron");
  PdgElectron->AddCut(AliDielectronVarManager::kPdgCode, 11.);
  AliDielectronVarCuts* PdgPositron = new AliDielectronVarCuts("PdgElectron","PdgElectron");
  PdgPositron->AddCut(AliDielectronVarManager::kPdgCode, -11.);
  AliDielectronCutGroup* PdgLepton = new AliDielectronCutGroup("PdgLepton","PdgLepton",AliDielectronCutGroup::kCompOR);
  PdgLepton->AddCut(PdgElectron);
  PdgLepton->AddCut(PdgPositron);
  // Pdg selection cuts of V0 electrons

  AliDielectronVarCuts *MotherIsPhoton = new AliDielectronVarCuts("MotherIsPhoton","MotherIsPhoton");
  MotherIsPhoton->AddCut(AliDielectronVarManager::kPdgCodeMother, 22.);

  // TMVA weight file (BDT for ePID)
  // Standard BDT using only electron nSigma values and pt
  //TString weightFile = "alien:///alice/cern.ch/user/a/acapon/TMVAclassifiers/TMVAClassification_BDT.weights.xml";
  // BDT using electron nSigma values, pt, and the pion nSigma from the TPC
  TString weightFile;
  if(wSDD){
    weightFile = "alien:///alice/cern.ch/user/a/acapon/TMVAclassifiers/TMVAClassification_BDT_18f3_wPsigTPC.weights.xml";
  }else{
    weightFile = "alien:///alice/cern.ch/user/a/acapon/TMVAclassifiers/TMVAClassification_BDT_18f3_wPsigTPC_woSDD.weights.xml";
  }
  AliDielectronTMVACuts* pidCuts = new AliDielectronTMVACuts("PIDCutsTMVA","PIDCutsTMVA");
  pidCuts->AddTMVAInput("pt",          AliDielectronVarManager::kPt);
  pidCuts->AddTMVAInput("EsigTPC",     AliDielectronVarManager::kTPCnSigmaEle);
  if(wSDD){
    pidCuts->AddTMVAInput("EsigITScorr", AliDielectronVarManager::kITSnSigmaEle);
  }
  pidCuts->AddTMVAInput("EsigTOFcorr", AliDielectronVarManager::kTOFnSigmaEle);
  pidCuts->AddTMVAInput("PsigTPC",     AliDielectronVarManager::kTPCnSigmaPio);
  pidCuts->AddTMVASpectator("pdg",     AliDielectronVarManager::kPdgCode);
  pidCuts->SetTMVAWeights("BDT", weightFile.Data());

  switch(PIDcuts){
    case kElectrons:
      if(wSDD){
        cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -3.0,  1.0, 0.2, 100., kFALSE);
        cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -1.5,  4.0, 0.2, 100., kFALSE);
        cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kPion,     -100., 3.5, 0.2, 100., kTRUE);
        cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3.0,  3.0, 0.2, 100., kFALSE, AliDielectronPID::kIfAvailable);
        cuts->AddCut(cutsPID);
        return cuts;
      }else{
        cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3.0,  3.0, 0.2, 100., kFALSE);
        cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kPion,     -100., 4.0, 0.2, 100., kTRUE);
        cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3.0,  3.0, 0.4, 100., kFALSE, AliDielectronPID::kRequire);
        cuts->AddCut(cutsPID);
        return cuts;
      }
      break;
    case kAllSpecies:
      break;
    case kCutSet1:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      if(wSDD){
        Printf("Using BDT cut value: %f",0.15);
        pidCuts->SetTMVACutValue(0.15);
      }else{
        Printf("Using BDT cut value: %f",0.15);
        pidCuts->SetTMVACutValue(0.15);
      }
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kTTreeCuts:
      // PID cuts used during TTree creating
      // Momentum range relaxed as it cuts on P not Pt. Kinematic cuts applied separately.
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      //cuts->AddCut(cutsPID);
      return cutsPID;
    case kV0_TPCcorr:
      // PID cuts used to select out a very pure sample of V0 electrons using only ITS (if available) and TOF
      if(wSDD){
        cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -1., 1., 0.1, 100., kFALSE);
        cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1., 1., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
      }else{
        cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3., 3., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
      }
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3.0, 3.0, 0.1, 100., kFALSE, AliDielectronPID::kRequire);
      cuts->AddCut(cutsPID);
      cuts->Print();
      return cuts;
    case kV0_ITScorr:
      // PID cuts used to select out a very pure sample of V0 electrons using only TPC and TOF
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -0.6, 1., 0.1, 100., kFALSE);
      cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1.,  1., 0.1, 0.4,  kFALSE, AliDielectronPID::kIfAvailable);
      cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1.,  1., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
      cuts->AddCut(cutsPID);
      cuts->Print();
      return cuts;
    case kV0_TOFcorr:
      // PID cuts used to select out a very pure sample of V0 electrons using only ITS (if available) and TPC
      if(wSDD){
        cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -1., 1., 0.1, 100., kFALSE);
        cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -0.6, 1., 0.1, 100., kFALSE);
      }else{
        cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -3., 3., 0.1, 100., kFALSE);
        cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -3., 3., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
      }
      cuts->AddCut(cutsPID);
      cuts->Print();
      return cuts;
    case kPdgSel:
      cuts->AddCut(PdgLepton);
      return cuts;
    case kPIDcut1:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.03);
      pidCuts->SetTMVACutValue(0.03);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut2:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.05);
      pidCuts->SetTMVACutValue(0.05);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut3:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.07);
      pidCuts->SetTMVACutValue(0.07);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut4:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.09);
      pidCuts->SetTMVACutValue(0.09);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut5:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.11);
      pidCuts->SetTMVACutValue(0.11);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut6:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.13);
      pidCuts->SetTMVACutValue(0.13);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut7:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.15);
      pidCuts->SetTMVACutValue(0.15);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut8:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.17);
      pidCuts->SetTMVACutValue(0.17);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut9:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.19);
      pidCuts->SetTMVACutValue(0.19);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut10:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.21);
      pidCuts->SetTMVACutValue(0.21);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut11:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.23);
      pidCuts->SetTMVACutValue(0.23);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut12:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.25);
      pidCuts->SetTMVACutValue(0.25);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut13:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.27);
      pidCuts->SetTMVACutValue(0.27);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut14:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.29);
      pidCuts->SetTMVACutValue(0.29);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut15:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.31);
      pidCuts->SetTMVACutValue(0.31);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut16:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.33);
      pidCuts->SetTMVACutValue(0.33);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut17:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.35);
      pidCuts->SetTMVACutValue(0.35);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut18:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.37);
      pidCuts->SetTMVACutValue(0.37);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut19:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.40);
      pidCuts->SetTMVACutValue(0.40);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kPIDcut20:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.45);
      pidCuts->SetTMVACutValue(0.45);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kTheoPID:
      cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. , 4., 0.0, 100., kTRUE , AliDielectronPID::kRequire    , AliDielectronVarManager::kPt);
      cutsPID->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. , 3., 0.0, 100., kFALSE, AliDielectronPID::kIfAvailable, AliDielectronVarManager::kPt);
      if(wSDD){
        cutsPID->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. , 1., 0.0, 100., kFALSE, AliDielectronPID::kRequire    , AliDielectronVarManager::kPt);
      }
      cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5, 3., 0.0, 100., kFALSE, AliDielectronPID::kRequire    , AliDielectronVarManager::kPt);
      cuts->AddCut(cutsPID);
      cuts->Print();
      return cuts;
    case kTOFreq:
      cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. , 4. , 0.0, 100., kTRUE , AliDielectronPID::kRequire   , AliDielectronVarManager::kPt);
      cutsPID->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. , 3. , 0.0, 100., kFALSE, AliDielectronPID::kRequire   , AliDielectronVarManager::kPt);
      if(wSDD){
        cutsPID->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. , 1. , 0.0, 100., kFALSE, AliDielectronPID::kRequire   , AliDielectronVarManager::kPt);
      }
      cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5, 3. , 0.0, 100., kFALSE, AliDielectronPID::kRequire   , AliDielectronVarManager::kPt);
      cuts->AddCut(cutsPID);
      cuts->Print();
      return cuts;
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
    case kCutVar1:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.05);
      pidCuts->SetTMVACutValue(0.05);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar2:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.15);
      pidCuts->SetTMVACutValue(0.15);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar3:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.31);
      pidCuts->SetTMVACutValue(0.31);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar4:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.03);
      pidCuts->SetTMVACutValue(0.03);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar5:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.05);
      pidCuts->SetTMVACutValue(0.05);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar6:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.15);
      pidCuts->SetTMVACutValue(0.15);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar7:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.00);
      pidCuts->SetTMVACutValue(0.00);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar8:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.45);
      pidCuts->SetTMVACutValue(0.45);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar9:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.08);
      pidCuts->SetTMVACutValue(0.08);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar10:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.31);
      pidCuts->SetTMVACutValue(0.31);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar11:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.15);
      pidCuts->SetTMVACutValue(0.15);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar12:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.03);
      pidCuts->SetTMVACutValue(0.03);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar13:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.00);
      pidCuts->SetTMVACutValue(0.00);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar14:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.22);
      pidCuts->SetTMVACutValue(0.22);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar15:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.05);
      pidCuts->SetTMVACutValue(0.05);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar16:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.24);
      pidCuts->SetTMVACutValue(0.24);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar17:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.38);
      pidCuts->SetTMVACutValue(0.38);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar18:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.12);
      pidCuts->SetTMVACutValue(0.12);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar19:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.10);
      pidCuts->SetTMVACutValue(0.10);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
    case kCutVar20:
      cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
      Printf("Use TMVA cut value = %f",0.45);
      pidCuts->SetTMVACutValue(0.45);
      cuts->AddCut(cutsPID);
      cuts->AddCut(pidCuts);
      return cuts;
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
    case kAllSpecies:
    case kElectrons:
      varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 20.);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,  - 1.0, 1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,   - 3.0, 3.0);
      if(wSDD){
          varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0);
      }else{
          varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  2.0, 100.0);
      }
      //varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0, 0.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 4.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      60.0, 200.); //clusters
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    70.0, 200.); //findable
      //varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0, 6.0);
      //varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCrFrac,  0.3, 10.); //Number of found/findable
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1); //Crossed rows over findable

      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);

      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(varCutsFilter);

      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kTTreeCuts:
      varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 20.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      60.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,  - 3.0, 3.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,   - 4.0, 4.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0); // < 3
        varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 36.);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  2.0, 100.0); // < 2
        varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 40.);
      }
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutSet1:
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 3
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);

      // Select filterbit 4
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);

      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutSet2:
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 3
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);

      // Select filterbit 4
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);

      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutSet3:
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   10.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 3
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);

      // Select filterbit 4
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);

      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kV0_trackCuts: // Does not work for MC (checked 2019.05.08)
      // V0 specific track cuts
      gammaV0cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);
      // Cut on the angle between the total momentum vector of the daughter
      // tracks and a line connecting the primary and secondary vertices
      gammaV0cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02), 1.0);
      gammaV0cuts->AddCut(AliDielectronVarManager::kChi2NDF, 0.0, 10.0);
      // Restrict distance between legs
      gammaV0cuts->AddCut(AliDielectronVarManager::kLegDist, 0.0, 0.25);
      // Require minimum distance to secondary vertex
      gammaV0cuts->AddCut(AliDielectronVarManager::kR, 3.0, 90.0);
      // Angle between daughter momentum plane and plane perpendicular to magnetic field
      gammaV0cuts->AddCut(AliDielectronVarManager::kPsiPair, 0.0, 0.05);
      // Mass cut on V0 (mother) particle
      gammaV0cuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.05);
      // Armenteros-Podolanksi variables
      // Pt
      gammaV0cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0, 0.05);
      // Longitudinal momentum asymmentry between daughter particles
      gammaV0cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35, 0.35);
      // Default setting is to exclude V0 tracks
      gammaV0cuts->SetExcludeTracks(kFALSE);
      // Standard track cut variables
      trackCutsV0->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 4.0);
      trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCr, 100.0, 160.0);
      trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8, 1.1);

      trackCuts->AddCut(gammaV0cuts);
      trackCuts->AddCut(trackCutsV0);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kV0_allAcc: // Cut setting to check MC V0 features
      // No V0 track cuts applied as want to see distributions
      gammaV0cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);
      // Default setting is to exclude V0 tracks
      gammaV0cuts->SetExcludeTracks(kFALSE);
      /* // Standard track cut variables */
      trackCutsV0->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 4.0);
      trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCr, 100.0, 160.0);
      trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8, 1.1);
      trackCuts->AddCut(gammaV0cuts);
      trackCuts->AddCut(trackCutsV0);
      /* trackCuts->AddCut(GetPIDCuts(PIDcuts)); */
      trackCuts->Print();
      return trackCuts;
    case kMCsel:
      // Despite pdg selection, apply loose cuts to trim away edges of
      // distribution which has poor data/MC agreement
      varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 20.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      60.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,  - 3.0, 3.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,   - 4.0, 4.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0); // < 3
        varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 36.);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  2.0, 100.0); // < 2
        varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 40.);
      }
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);

      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    // ################## Cut setting used to obtain resolution files
    case kResolutionTrackCuts:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.05, 50.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta, -1.2, 1.2);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
        varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,  15.0);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      2.0, 100.0);
        varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,  20.0);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   4.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   5.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 1.1);

      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->Print();
      return trackCuts;
      // ################# Track cut variations start here
    case kCutVar1:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -2.0,  2.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar2:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        70.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      80.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.7,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.5,  1.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   3.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar3:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        90.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      110.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.6,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -0.5,  0.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      4.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   3.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar4:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      90.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -0.5,  0.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -4.0,  4.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   3.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar5:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        100.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar6:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar7:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        70.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      80.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.6,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -2.0,  2.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   7.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar8:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        90.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.9,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -2.0,  2.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   100.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar9:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        90.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      80.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.6,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -0.5,  0.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar10:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        90.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      110.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.5,  1.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar11:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        100.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      90.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.7,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -2.0,  2.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   3.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar12:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        100.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.6,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -0.5,  0.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      4.0,   100.0);
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar13:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      80.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.5,  1.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar14:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.5,  1.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   7.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar15:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        100.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      80.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.5,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.5,  1.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -4.0,  4.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   3.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar16:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        70.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      110.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.6,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   6.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar17:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        70.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      80.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.7,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.5,  1.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      4.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar18:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        100.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      90.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -2.0,  2.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   5.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar19:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        90.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -0.5,  0.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -2.0,  2.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   5.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kCutVar20:
      // Same cuts for all cut settings
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   20.);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      // Varied cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.9,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 5
        // Shared cluster cut the same for all cut settings
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   100.0);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kScheidCuts:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,   -1., 1.);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,    -3., 3.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100., 160.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        100., 160.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8, 1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0, 4.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,        3. , 10.);
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0, 5.5);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,       1.0, 6.0, kTRUE);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kGoodEtaPhi:

      // Standard track cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 3
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);

      // Select filterbit 4
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);

      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));

      // Select out eta/phi regions
      // Big OR cut group for all regions
      AliDielectronCutGroup* etaPhiRegions = new AliDielectronCutGroup("etaPhiRegions", "etaPhiRegions", AliDielectronCutGroup::kCompOR);

      // Define smaller regions with AND eta/phi selections
      AliDielectronCutGroup* etaPhiRegion1 = new AliDielectronCutGroup("etaPhiRegion1", "etaPhiRegion1", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut1        = new AliDielectronVarCuts("etaCut1", "etaCut1");
      AliDielectronVarCuts* phiCut1        = new AliDielectronVarCuts("phiCut1", "phiCut1");
      etaCut1->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
      phiCut1->AddCut(AliDielectronVarManager::kPhi,  0.0, 1.0);
      etaPhiRegion1->AddCut(etaCut1);
      etaPhiRegion1->AddCut(phiCut1);

      AliDielectronCutGroup* etaPhiRegion2 = new AliDielectronCutGroup("etaPhiRegion2", "etaPhiRegion2", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut2        = new AliDielectronVarCuts("etaCut2", "etaCut2");
      AliDielectronVarCuts* phiCut2        = new AliDielectronVarCuts("phiCut2", "phiCut2");
      etaCut2->AddCut(AliDielectronVarManager::kEta, -0.8, -0.5);
      phiCut2->AddCut(AliDielectronVarManager::kPhi,  3.2, 3.6);
      etaPhiRegion2->AddCut(etaCut2);
      etaPhiRegion2->AddCut(phiCut2);

      AliDielectronCutGroup* etaPhiRegion3 = new AliDielectronCutGroup("etaPhiRegion3", "etaPhiRegion3", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut3        = new AliDielectronVarCuts("etaCut3", "etaCut3");
      AliDielectronVarCuts* phiCut3        = new AliDielectronVarCuts("phiCut3", "phiCut3");
      etaCut3->AddCut(AliDielectronVarManager::kEta, 0.5, 0.8);
      phiCut3->AddCut(AliDielectronVarManager::kPhi, 3.0, 3.4);
      etaPhiRegion3->AddCut(etaCut3);
      etaPhiRegion3->AddCut(phiCut3);

      AliDielectronCutGroup* etaPhiRegion4 = new AliDielectronCutGroup("etaPhiRegion4", "etaPhiRegion4", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut4        = new AliDielectronVarCuts("etaCut4", "etaCut4");
      AliDielectronVarCuts* phiCut4        = new AliDielectronVarCuts("phiCut4", "phiCut4");
      etaCut4->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
      phiCut4->AddCut(AliDielectronVarManager::kPhi,  4.5, 4.9);
      etaPhiRegion4->AddCut(etaCut4);
      etaPhiRegion4->AddCut(phiCut4);

      etaPhiRegions->AddCut(etaPhiRegion1);
      etaPhiRegions->AddCut(etaPhiRegion2);
      etaPhiRegions->AddCut(etaPhiRegion3);
      etaPhiRegions->AddCut(etaPhiRegion4);

      // Final group contianing eta/phi regions AND standard track+PID cuts
      AliDielectronCutGroup* finalCutGroup = new AliDielectronCutGroup("finalCutGroup", "finalCutGroup", AliDielectronCutGroup::kCompAND);
      finalCutGroup->AddCut(trackCuts);
      finalCutGroup->AddCut(etaPhiRegions);

      finalCutGroup->Print();
      return finalCutGroup;
    case kBadEtaPhi:

      // Standard track cuts
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
      varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
      varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
      if(wSDD){
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      5.0,   100.0); // < 5
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }else{
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      3.0,   100.0); // < 3
        varCutsFilter->AddCut(AliDielectronVarManager::kNclsSFracITS, 0.0,   0.01);
      }
      varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);

      // Select filterbit 4
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
      trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
      // Refits
      trackCutsFilter->SetRequireITSRefit(kTRUE);
      trackCutsFilter->SetRequireTPCRefit(kTRUE);

      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));

      // Select out eta/phi regions
      // Big OR cut group for all regions
      AliDielectronCutGroup* etaPhiRegions = new AliDielectronCutGroup("etaPhiRegions", "etaPhiRegions", AliDielectronCutGroup::kCompOR);

      // Define smaller regions with AND eta/phi selections
      AliDielectronCutGroup* etaPhiRegion1 = new AliDielectronCutGroup("etaPhiRegion1", "etaPhiRegion1", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut1        = new AliDielectronVarCuts("etaCut1", "etaCut1");
      AliDielectronVarCuts* phiCut1        = new AliDielectronVarCuts("phiCut1", "phiCut1");
      etaCut1->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
      phiCut1->AddCut(AliDielectronVarManager::kPhi,  1.0, 3.0);
      etaPhiRegion1->AddCut(etaCut1);
      etaPhiRegion1->AddCut(phiCut1);

      AliDielectronCutGroup* etaPhiRegion2 = new AliDielectronCutGroup("etaPhiRegion2", "etaPhiRegion2", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut2        = new AliDielectronVarCuts("etaCut2", "etaCut2");
      AliDielectronVarCuts* phiCut2        = new AliDielectronVarCuts("phiCut2", "phiCut2");
      etaCut2->AddCut(AliDielectronVarManager::kEta, -0.5, 0.5);
      phiCut2->AddCut(AliDielectronVarManager::kPhi,  3.0, 4.5);
      etaPhiRegion2->AddCut(etaCut2);
      etaPhiRegion2->AddCut(phiCut2);

      AliDielectronCutGroup* etaPhiRegion3 = new AliDielectronCutGroup("etaPhiRegion3", "etaPhiRegion3", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut3        = new AliDielectronVarCuts("etaCut3", "etaCut3");
      AliDielectronVarCuts* phiCut3        = new AliDielectronVarCuts("phiCut3", "phiCut3");
      etaCut3->AddCut(AliDielectronVarManager::kEta, -0.8, -0.5);
      phiCut3->AddCut(AliDielectronVarManager::kPhi, 3.0, 3.2);
      etaPhiRegion3->AddCut(etaCut3);
      etaPhiRegion3->AddCut(phiCut3);

      AliDielectronCutGroup* etaPhiRegion4 = new AliDielectronCutGroup("etaPhiRegion4", "etaPhiRegion4", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut4        = new AliDielectronVarCuts("etaCut4", "etaCut4");
      AliDielectronVarCuts* phiCut4        = new AliDielectronVarCuts("phiCut4", "phiCut4");
      etaCut4->AddCut(AliDielectronVarManager::kEta, -0.8 -0.5);
      phiCut4->AddCut(AliDielectronVarManager::kPhi,  3.4, 4.5);
      etaPhiRegion4->AddCut(etaCut4);
      etaPhiRegion4->AddCut(phiCut4);

      AliDielectronCutGroup* etaPhiRegion5 = new AliDielectronCutGroup("etaPhiRegion5", "etaPhiRegion5", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut5        = new AliDielectronVarCuts("etaCut5", "etaCut5");
      AliDielectronVarCuts* phiCut5        = new AliDielectronVarCuts("phiCut5", "phiCut5");
      etaCut5->AddCut(AliDielectronVarManager::kEta, 0.5 0.8);
      phiCut5->AddCut(AliDielectronVarManager::kPhi,  3.4, 4.5);
      etaPhiRegion5->AddCut(etaCut5);
      etaPhiRegion5->AddCut(phiCut5);

      AliDielectronCutGroup* etaPhiRegion6 = new AliDielectronCutGroup("etaPhiRegion6", "etaPhiRegion6", AliDielectronCutGroup::kCompAND);
      AliDielectronVarCuts* etaCut6        = new AliDielectronVarCuts("etaCut6", "etaCut6");
      AliDielectronVarCuts* phiCut6        = new AliDielectronVarCuts("phiCut6", "phiCut6");
      etaCut6->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
      phiCut6->AddCut(AliDielectronVarManager::kPhi,  4.9, TMath::Pi()*2);
      etaPhiRegion6->AddCut(etaCut6);
      etaPhiRegion6->AddCut(phiCut6);

      etaPhiRegions->AddCut(etaPhiRegion1);
      etaPhiRegions->AddCut(etaPhiRegion2);
      etaPhiRegions->AddCut(etaPhiRegion3);
      etaPhiRegions->AddCut(etaPhiRegion4);
      etaPhiRegions->AddCut(etaPhiRegion5);
      etaPhiRegions->AddCut(etaPhiRegion6);

      // Final group contianing eta/phi regions AND standard track+PID cuts
      AliDielectronCutGroup* finalCutGroup = new AliDielectronCutGroup("finalCutGroup", "finalCutGroup", AliDielectronCutGroup::kCompAND);
      finalCutGroup->AddCut(trackCuts);
      finalCutGroup->AddCut(etaPhiRegions);

      finalCutGroup->Print();
      return finalCutGroup;
    case kITSshared1:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
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
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,       2.0, 6.0, kTRUE);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kITSshared2:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
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
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,       3.0, 6.0, kTRUE);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kITSshared3:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
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
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,       4.0, 6.0, kTRUE);
      trackCuts->AddCut(trackCutsFilter);
      trackCuts->AddCut(varCutsFilter);
      trackCuts->AddCut(GetPIDCuts(PIDcuts));
      trackCuts->Print();
      return trackCuts;
    case kITSshared4:
      varCutsFilter->AddCut(AliDielectronVarManager::kPt,  0.2, 100.);
      varCutsFilter->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);
      trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); //or 1<<4
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
      varCutsFilter->AddCut(AliDielectronVarManager::kNclsSITS,       5.0, 6.0, kTRUE);
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
