class LMEECutLib {

	public:

	enum LMEECutSet{
		kAllSpecies,
		kElectrons,
		kHighMult,
		kMidLowMult,
		kMidMult,
		kLowMult,
		kTTreeCuts,
		kCutSet1,
		kV0_ITScorr,
		kV0_TPCcorr,
		kV0_TOFcorr
	};


	LMEECutLib(Bool_t wSDD): wSDD(wSDD){
		
		::Info("LMEECutLib_acapon", "Creating new LMEECutLib");
		fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
	}

	//Getters
	AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
	AliAnalysisCuts*            GetCentralityCuts(Int_t centSel);
	AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);

	AliDielectronCutGroup* GetPairCuts(Int_t cutSet);
	AliAnalysisCuts* GetPIDCuts(Int_t PIDcuts);
	//AliDielectronPID* GetPIDCuts(Int_t PIDcuts);
	AliAnalysisCuts* GetTrackCuts(Int_t cutSet, Int_t PIDcuts);

	//PID correction functions used within dielectron framework
  void SetEtaCorrectionTPC(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void SetEtaCorrectionITS(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void SetEtaCorrectionTOF(AliDielectron *die, Int_t corrZdim, Int_t corrYdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/

	//PID correction function used by SimpleTTreeMaker
	//i.e it doesn't need an AliDielectron object
	static TH3D SetEtaCorrectionTPCTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise);
	static TH3D SetEtaCorrectionITSTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise);
	static TH3D SetEtaCorrectionTOFTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise);

	static TBits* fUsedVars; //Used Variables for correction
	TH1* fPostPIDCntrdCorrTPC; //Post PID correction object for electron sigma centroids in TPC
	TH1* fPostPIDWdthCorrTPC; //Post PID correction object for electron sigma widths in TPC
	TH1* fPostPIDCntrdCorrITS; //Post PID correction object for electron sigma centroids in ITS
	TH1* fPostPIDWdthCorrITS; //Post PID correction object for electron sigma widths in ITS
	TH1* fPostPIDCntrdCorrTOF; //Post PID correction object for electron sigma centroids in TOF
	TH1* fPostPIDWdthCorrTOF; //Post PID correction object for electron sigma widths in TOF

	private:
			Bool_t wSDD;
};

void LMEECutLib::SetEtaCorrectionTPC(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  std::string file_name = "/home/aaron/Data/PIDcalibration/outputTPC.root";

  TFile* inFile = TFile::Open(file_name.c_str(), "READ");
  std::cout << inFile << std::endl;
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTPC.root .");
    std::cout << "Copy TPC correction from Alien" << std::endl;
    inFile = TFile::Open("outputTPC.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }
  if(runwise){
    TObjArray* arr_mean = dynamic_cast<TObjArray*>(inFile->Get("mean_correction_arr"));
    TObjArray* arr_width =dynamic_cast<TObjArray*>( inFile->Get("width_correction_arr"));
    std::cout << arr_mean << " " << arr_width << std::endl;

    die->SetWidthCorrArr(arr_width, kTRUE, corrXdim, corrYdim, corrZdim);
    die->SetCentroidCorrArr(arr_mean, kTRUE, corrXdim, corrYdim, corrZdim);
  }
  else{
    TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
    die->SetCentroidCorrFunction(mean, corrXdim, corrYdim, corrZdim);
    die->SetWidthCorrFunction(width, corrXdim, corrYdim, corrZdim);
  }

}

void LMEECutLib::SetEtaCorrectionITS(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionITS()\n";
  std::string file_name = "/home/aaron/Data/PIDcalibration/outputITS.root";

  TFile* inFile = TFile::Open(file_name.c_str());
  std::cout << inFile << std::endl;
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputITS.root .");
    std::cout << "Copy ITS correction from Alien" << std::endl;
    inFile = TFile::Open("outputITS.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }
  if(runwise){
    TObjArray* arr_mean = dynamic_cast<TObjArray*>(inFile->Get("mean_correction_arr"));
    TObjArray* arr_width =dynamic_cast<TObjArray*>( inFile->Get("width_correction_arr"));
    std::cout << arr_mean << " " << arr_width << std::endl;

    die->SetWidthCorrArrITS(arr_width, kTRUE, corrXdim, corrYdim, corrZdim);
    die->SetCentroidCorrArrITS(arr_mean, kTRUE, corrXdim, corrYdim, corrZdim);
  }
  else{
    TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
    die->SetCentroidCorrFunctionITS(mean, corrXdim, corrYdim, corrZdim);
    die->SetWidthCorrFunctionITS(width, corrXdim, corrYdim, corrZdim);
  }

}
void LMEECutLib::SetEtaCorrectionTOF(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise) {
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionTOF()\n";
  std::string file_name = "/home/aaron/Data/PIDcalibration/outputTOF.root";

  TFile* inFile = TFile::Open(file_name.c_str());
  std::cout << inFile << std::endl;
  if(!inFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTOF.root .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    inFile = TFile::Open("outputTOF.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }
  if(runwise){
    TObjArray* arr_mean = dynamic_cast<TObjArray*>(inFile->Get("mean_correction_arr"));
    TObjArray* arr_width =dynamic_cast<TObjArray*>( inFile->Get("width_correction_arr"));
    std::cout << arr_mean << " " << arr_width << std::endl;

    die->SetWidthCorrArrTOF(arr_width, kTRUE, corrXdim, corrYdim, corrZdim);
    die->SetCentroidCorrArrTOF(arr_mean, kTRUE, corrXdim, corrYdim, corrZdim);
  }
  else{
    TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
    die->SetCentroidCorrFunctionTOF(mean, corrXdim, corrYdim, corrZdim);
    die->SetWidthCorrFunctionTOF(width, corrXdim, corrYdim, corrZdim);
  }

}
static TH3D LMEECutLib::SetEtaCorrectionTPCTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise, Int_t selection) {
	

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionTPC() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  std::string file_name = "/home/aaron/Data/PIDcalibration/outputTPC.root";

  TFile* recalFile = TFile::Open(file_name.c_str());
  std::cout << recalFile << std::endl;
  if(!recalFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTPC.root .");
    std::cout << "Copy TPC correction from Alien" << std::endl;
    recalFile = TFile::Open("outputTPC.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }  
  TH3D* mean  = dynamic_cast<TH3D*>(recalFile->Get("sum_mean_correction"));
  TH3D* width = dynamic_cast<TH3D*>(recalFile->Get("sum_width_correction"));
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

static TH3D LMEECutLib::SetEtaCorrectionITSTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise, Int_t selection) {
	

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionITSTTree() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionITSTTree()\n";
  std::string file_name = "/home/aaron/Data/PIDcalibration/outputITS.root";

  TFile* recalFile = TFile::Open(file_name.c_str());
  std::cout << recalFile << std::endl;
  if(!recalFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputITS.root .");
    std::cout << "Copy ITS correction from Alien" << std::endl;
    recalFile = TFile::Open("outputITS.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }  
  TH3D* mean  = dynamic_cast<TH3D*>(recalFile->Get("sum_mean_correction"));
  TH3D* width = dynamic_cast<TH3D*>(recalFile->Get("sum_width_correction"));
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

static TH3D LMEECutLib::SetEtaCorrectionTOFTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t runwise, Int_t selection) {
	

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionTOFTTree() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionTOFTTree()\n";
  std::string file_name = "/home/aaron/Data/PIDcalibration/outputTOF.root";

  TFile* recalFile = TFile::Open(file_name.c_str());
  std::cout << recalFile << std::endl;
  if(!recalFile){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTOF.root .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    recalFile = TFile::Open("outputTOF.root");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }  
  TH3D* mean  = dynamic_cast<TH3D*>(recalFile->Get("sum_mean_correction"));
  TH3D* width = dynamic_cast<TH3D*>(recalFile->Get("sum_width_correction"));
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
AliDielectronEventCuts* LMEECutLib::GetEventCuts(Int_t cutSet) {

    ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetEventCuts() >>>>>>>>>>>>>>>>>>>>>> ");
    AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCuts_acapon","Vertex Track && |vtxZ|<10 && ncontrib>0");

    switch(cutSet){
        case kAllSpecies:
        case kElectrons:
        case kHighMult:
        case kMidMult:
        case kMidLowMult:
        case kLowMult:
				case kTTreeCuts:
				case kV0_TPCcorr:
				case kV0_ITScorr:
				case kV0_TOFcorr:
				case kCutSet1:
            eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
            eventCuts->SetRequireVertex();
            eventCuts->SetMinVtxContributors(1);
            eventCuts->SetVertexZ(-10.,10.);
            break;
				default: cout << "No Event Cut defined" << endl;
    }
    return eventCuts;
}


//Centrality selection done in Event selection
AliAnalysisCuts* LMEECutLib::GetCentralityCuts(Int_t centSel) {
	AliDielectronVarCuts* centCuts = 0x0;
	switch(centSel){
		case kAllSpecies:
		case kElectrons:
		case kTTreeCuts:
		case kCutSet1:
		case kV0_TPCcorr:
		case kV0_ITScorr:
		case kV0_TOFcorr:
			centCuts = new AliDielectronVarCuts("centCutsHigh","cent0090");
			centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.,100.);
			break;
			break;
		case kHighMult:
			centCuts = new AliDielectronVarCuts("centCutsHigh","cent0020");
			centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.,20.);
			break;
	 case kMidMult:
			centCuts = new AliDielectronVarCuts("centCutsMid","cent2040");
			centCuts->AddCut(AliDielectronVarManager::kCentralityNew,20.,40.);
			break;
	 case kMidLowMult:
			centCuts = new AliDielectronVarCuts("centCutsMid","cent0460");
			centCuts->AddCut(AliDielectronVarManager::kCentralityNew,40.,60.);
			break;
		case kLowMult:
			centCuts  = new AliDielectronVarCuts("centCutsLow","cent60100");
			centCuts->AddCut(AliDielectronVarManager::kCentralityNew,60.,100.);
			break;
		default: cout << "No Centrality selected" << endl;
	}
	return centCuts;
}


AliDielectronMixingHandler* LMEECutLib::GetMixingHandler(Int_t cutSet) {
	AliDielectronMixingHandler* mixingHandler = 0x0;
	switch (cutSet) {
		case kAllSpecies:
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kMidLowMult:
		case kLowMult:
		case kCutSet1:
			mixingHandler = new AliDielectronMixingHandler;
			mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
			//mixingHandler->AddVariable(AliDielectronVarManager::kNacc,"0,500");
			mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100");
			mixingHandler->SetDepth(20);
			mixingHandler->SetMixType(AliDielectronMixingHandler::kAll);
			break;
		//[...]
		case kTTreeCuts:
			break;
		default: cout << "No Mixer defined" << endl;
	}
	return mixingHandler;
}



//Pair Cuts for Analysis step - take care of logic - inverted compared to other PairCuts!!
// cuts = SELECTION!!!
AliDielectronCutGroup* LMEECutLib::GetPairCuts(Int_t cutSet)  {

	::Info("LMEECutLibg_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetPairCuts() >>>>>>>>>>>>>>>>>>>>>> ");
	//Final OR cut group to incorporate the following cuts (below)
	AliDielectronCutGroup* allCuts    = new AliDielectronCutGroup("allCuts", "allCuts", AliDielectronCutGroup::kCompOR);

	case kCutSet1:
     //AND cut group to select low mass pairs with large opening angle
    AliDielectronCutGroup* convRejCut = new AliDielectronCutGroup("convRejCut", "convRejCut", AliDielectronCutGroup::kCompAND);
    AliDielectronVarCuts* convMassCut = new AliDielectronVarCuts("convMassCut", "convMassCut");
    AliDielectronVarCuts* convPhiVCut = new AliDielectronVarCuts("convPhiVCut", "convPhiVCut");
    convMassCut->AddCut(AliDielectronVarManager::kM, 0.00, 0.1);
    convPhiVCut->AddCut(AliDielectronVarManager::kPhivPair, 0., 2.);
    convRejCut->AddCut(convMassCut);
    convRejCut->AddCut(convPhiVCut);

    //Mass cut to include any pairs with mass greater than 0.1 GeV
    AliDielectronVarCuts* pairMassCut = new AliDielectronVarCuts("pairMassCut", "pairMassCut");
    pairMassCut->AddCut(AliDielectronVarManager::kM, 0.1, 5.0);

    allCuts->AddCut(convRejCut);
    allCuts->AddCut(pairMassCut);

    return allCuts;
		break;
	default:
		std::cout << "No pair cut applied" << std::endl;

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
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kMidLowMult:
		case kLowMult:
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

			// Copy xml file to local directory first
			// Classifier trained on reweighted data
			//TString weightFile = "alien:///alice/cern.ch/user/a/acapon/TMVAclassifiers/TMVAClassification_BDTunweighted.weights.xml";
			// Classifier trained on PID post calibratied data
			TString weightFile = "alien:///alice/cern.ch/user/a/acapon/TMVAclassifiers/TMVAClassification_BDT.weights.xml";
			
			AliDielectronTMVACuts *pidCuts = new AliDielectronTMVACuts("PIDCutsTMVA","PIDCutsTMVA");
			pidCuts->AddTMVAInput("pt", AliDielectronVarManager::kPt);
			pidCuts->AddTMVAInput("EsigTPC", AliDielectronVarManager::kTPCnSigmaEle);
			pidCuts->AddTMVAInput("EsigITS", AliDielectronVarManager::kITSnSigmaEle);
			pidCuts->AddTMVAInput("EsigTOF", AliDielectronVarManager::kTOFnSigmaEle);
			pidCuts->AddTMVAInput("PsigTPC", AliDielectronVarManager::kTPCnSigmaPio);
			pidCuts->AddTMVASpectator("pdg", AliDielectronVarManager::kPdgCode);
			pidCuts->SetTMVAWeights("BDT method", weightFile.Data());

			Printf("Use TMVA cut value = %f",0.1);
			pidCuts->SetTMVACutValue(0.1);

			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			cuts->Print();
			return cuts;
			break;
		case kTTreeCuts:
		  // PID cuts used during TTree creating
			// Momentum range relaxed as it cuts on P not Pt. Kinematic cuts applied
			// separately.
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			//cuts->AddCut(cutsPID);
			return cutsPID;
			break;
		case kV0_TPCcorr:
			// PID cuts used to select out a very pure sample of V0 electrons using only ITS and TOF
			cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -1., 1., 0.1, 100., kFALSE);
			cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1., 1., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
			cuts->AddCut(cutsPID);
			return cuts;
			break;
		case kV0_ITScorr:
			// PID cuts used to select out a very pure sample of V0 electrons using only TPC and TOF
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -0.6, 1., 0.1, 100., kFALSE);
			cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1.,  1., 0.1, 0.4,  kFALSE, AliDielectronPID::kIfAvailable);
			cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1.,  1., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
			cuts->AddCut(cutsPID);
			return cuts;
			break;
		case kV0_TOFcorr:
			// PID cuts used to select out a very pure sample of V0 electrons using only TPC and TOF
			cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -1., 1., 0.1, 100., kFALSE);
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -0.6, 1., 0.1, 100., kFALSE);
			cuts->AddCut(cutsPID);
			return cuts;
			break;
		default:
			std::cout << "No Analysis PID Cut defined " << std::endl;
			return 0x0;
    }
		std::cout << "No Analysis PID Cut defined " << std::endl;
    return 0x0;
}

//Make/Tighten track Cuts that are *NOT* already
//done in the AOD production
//**IMPORTANT**: For AODs, select FilterBit
//the method is ignored for ESDs

AliDielectronCutGroup* LMEECutLib::GetTrackCuts(Int_t cutSet, Int_t PIDcuts){

	::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> GetTrackCuts() >>>>>>>>>>>>>>>>>>>>>> ");
	AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("trackCuts", "trackCuts", AliDielectronCutGroup::kCompAND);

	AliDielectronVarCuts* varCutsFilter     = new AliDielectronVarCuts("varCutsFilter", "varCutsFilter");
	AliDielectronTrackCuts* trackCutsFilter = new AliDielectronTrackCuts("trackCutsFilter", "trackCutsFilter");

	switch(cutSet){
    //----------
    // these MAIN settings just load the main track selection directly below:
    //----------
		case kAllSpecies:
		case kElectrons:
		case kHighMult:
		case kMidMult:
		case kMidLowMult:
		case kLowMult:
			varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
			varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
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

			trackCutsFilter->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			trackCutsFilter->SetRequireITSRefit(kTRUE);
			trackCutsFilter->SetRequireTPCRefit(kTRUE);

			trackCuts->AddCut(trackCutsFilter);
			trackCuts->AddCut(varCutsFilter);
			
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			return trackCuts;
			break;
		case kTTreeCuts:
			varCutsFilter->AddCut(AliDielectronVarManager::kEta, -0.80, 0.80);
			varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.2, 10.);
			//TPC cuts
			//Clusters
			varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 200.); 
			//Crossed rows
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      60.0, 200.); 
			//Crossed rows over findable
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.3, 1.1); 
			//DCA
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,  - 3.0, 3.0);
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,   - 4.0, 4.0);
			//ITS cuts
			if(wSDD){
				varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  3.0, 100.0); // < 3
			}else{
				varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,  1.0, 100.0); // < 1
			}
			varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0, 36.);
			
			//Select filterbit 4
			trackCutsFilter->SetAODFilterBit(16);//or 1<<4
			//Refits	
			trackCutsFilter->SetRequireITSRefit(kTRUE);
			trackCutsFilter->SetRequireTPCRefit(kTRUE);

			trackCuts->AddCut(varCutsFilter);
			trackCuts->AddCut(trackCutsFilter);
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			return trackCuts;
			break;
		case kCutSet1:

			varCutsFilter->AddCut(AliDielectronVarManager::kEta,            -0.80, 0.80);
			varCutsFilter->AddCut(AliDielectronVarManager::kPt,             0.2,   10.);
			varCutsFilter->AddCut(AliDielectronVarManager::kNclsTPC,        80.0,  200.);
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCr,      100.0, 200.);
			varCutsFilter->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8,   1.1);
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,  1.0);
			varCutsFilter->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,  3.0);
			if(wSDD){
				varCutsFilter->AddCut(AliDielectronVarManager::kNclsITS,      4.0,   100.0); // < 4
			}else{
				varCutsFilter->AddCut(AliDielectronVarManager::kNclsSfracITS, 0.0,   0.01); 
			}
			varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);
			//Select filterbit 4
			trackCutsFilter->SetAODFilterBit(16);//or 1<<4
			//Refits	
			trackCutsFilter->SetRequireITSRefit(kTRUE);
			trackCutsFilter->SetRequireTPCRefit(kTRUE);
	
			trackCuts->AddCut(varCutsFilter);
			trackCuts->AddCut(trackCutsFilter);
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			return trackCuts;
			break;
		case kV0_TPCcorr:
		case kV0_ITScorr:
		case kV0_TOFcorr:
			// V0 specific track cuts
			AliDielectronV0Cuts* gammaV0cuts = new AliDielectronV0Cuts("gammaV0cuts", "gammaV0cuts");
			gammaV0cuts->SetV0finder(AliDielectronV0Cuts::kOnTheFly);
			// Cut on the angle between the total momentum vector of the daughter
			// tracks and a line connecting the primary and secondary vertices
			gammaV0cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02), 1.0,  kFALSE);
			gammaV0cuts->AddCut(AliDielectronVarManager::kChi2NDF, 0.0, 10.0, kFALSE);
			// Restrict distance between legs
			gammaV0cuts->AddCut(AliDielectronVarManager::kLegDist, 0.0, 0.25, kFALSE);
			// Require minimum distance to secondary vertex
			gammaV0cuts->AddCut(AliDielectronVarManager::kR, 3.0, 90.0, kFALSE);
			// Angle between daughter momentum plane and plane perpendicular to magnetic field
			gammaV0cuts->AddCut(AliDielectronVarManager::kPsiPair, 0.0, 0.05, kFALSE);
			// Mass cut on V0 (mother) particle
			gammaV0cuts->AddCut(AliDielectronVarManager::kM, 0.0, 0.05, kFALSE);
			// Armenteros-Podolanksi variables
			// Pt
			gammaV0cuts->AddCut(AliDielectronVarManager::kArmPt, 0.0, 0.05, kFALSE);
			// Longitudinal momentum asymmentry between daughter particles
			gammaV0cuts->AddCut(AliDielectronVarManager::kArmAlpha, -0.35, 0.35, kFALSE); 
			// Default setting is to exclude V0 tracks
			gammaV0cuts->SetExcludeTracks(kFALSE);
			// Standard track cut variables
			AliDielectronVarCuts* trackCutsV0 = new AliDielectronVarCuts("trackCutsV0", "trackCutsV0");
			trackCutsV0->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 4.0);
			trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCr, 100.0, 160.0);
			trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8, 1.1);

			trackCuts->AddCut(gammaV0cuts);
			trackCuts->AddCut(trackCutsV0);
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			return trackCuts;
			break;
		default:
			std::cout << "No Analysis Track Cut defined" << std::endl;
		}
		std::cout << "Track cuts not applied...." << std::endl;
		return 0x0;
}

