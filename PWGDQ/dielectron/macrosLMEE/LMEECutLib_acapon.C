class LMEECutLib {

	public:

	enum LMEECutSet{
		kAllSpecies,
		kElectrons,
		kTTreeCuts,
		kCutSet1,
		kV0_ITScorr,
		kV0_TPCcorr,
		kV0_TOFcorr,
		kV0_trackCuts,
		kPdgSel,
		kMCsel,
		kResolutionTrackCuts,
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
		kTheoPID
	};



	LMEECutLib(Bool_t wSDD): wSDD(wSDD){
		
		::Info("LMEECutLib_acapon", "Creating new LMEECutLib");
		fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);

		// Set the twenty different MVA cut settings here
		Float_t PIDcutRange[20] = {-0.03, 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36, 0.39, 0.42, 0.45, 0.48, 0.51};
	}

	//Getters
	AliDielectronEventCuts*     GetEventCuts(Int_t cutSet);
	AliAnalysisCuts*            GetCentralityCuts(Int_t centSel);
	AliDielectronMixingHandler* GetMixingHandler(Int_t cutSet);

	AliDielectronCutGroup* GetPairCuts(Int_t cutSet);
	AliAnalysisCuts* GetPIDCuts(Int_t PIDcuts);
	AliAnalysisCuts* GetTrackCuts(Int_t cutSet, Int_t PIDcuts);

	//PID correction functions used within dielectron framework
  void SetEtaCorrectionTPC(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void SetEtaCorrectionITS(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC); //giving default value fails: /* = AliDielectronVarManager::kEta*/
  void SetEtaCorrectionTOF(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC); //giving default value fails: /* = AliDielectronVarManager::kEta*/

	//PID correction function used by SimpleTTreeMaker
	//i.e it doesn't need an AliDielectron object
	static TH3D SetEtaCorrectionTPCTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection );
	static TH3D SetEtaCorrectionITSTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC);
	static TH3D SetEtaCorrectionTOFTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC);

	static TBits* fUsedVars; //Used Variables for correction
	static TH1* fPostPIDCntrdCorrTPC; //Post PID correction object for electron sigma centroids in TPC
	static TH1* fPostPIDWdthCorrTPC; //Post PID correction object for electron sigma widths in TPC
	static TH1* fPostPIDCntrdCorrITS; //Post PID correction object for electron sigma centroids in ITS
	static TH1* fPostPIDWdthCorrITS; //Post PID correction object for electron sigma widths in ITS
	static TH1* fPostPIDCntrdCorrTOF; //Post PID correction object for electron sigma centroids in TOF
	static TH1* fPostPIDWdthCorrTOF; //Post PID correction object for electron sigma widths in TOF

	private:
			Bool_t wSDD;
};

void LMEECutLib::SetEtaCorrectionTPC(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim){
  
	// eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
	
  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  std::string file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/outputTPC.root";

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
	TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
	TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
	die->SetCentroidCorrFunction(mean, corrXdim, corrYdim, corrZdim);
	die->SetWidthCorrFunction(width, corrXdim, corrYdim, corrZdim);

}

void LMEECutLib::SetEtaCorrectionITS(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC){
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionITS()\n";
  TString file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/output";
	if(hasMC){
		file_name.Append("ITS_MC.root");
	}else{
		file_name.Append("ITS.root");
	}

  TFile* inFile = TFile::Open(file_name.Data());
  std::cout << inFile << std::endl;
  if(!inFile){
		if(hasMC){
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputITS_MC.root .");
		}else{
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputITS.root .");
		}
    std::cout << "Copy ITS correction from Alien" << std::endl;
		if(hasMC){
			inFile = TFile::Open("outputITS_MC.root");
		}else{
			inFile = TFile::Open("outputITS.root");
		}
  }
  else{
    std::cout << "Correction loaded" << std::endl;
  }
 
	TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
	TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
	die->SetCentroidCorrFunctionITS(mean, corrXdim, corrYdim, corrZdim);
	die->SetWidthCorrFunctionITS(width, corrXdim, corrYdim, corrZdim);

}

void LMEECutLib::SetEtaCorrectionTOF(AliDielectron *die, Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Bool_t hasMC){
  //
  // eta correction for the centroid and width of electron sigmas in the TPC, can be one/two/three-dimensional
  //
  std::cout << "starting LMEECutLib::SetEtaCorrectionTOF()\n";
  TString file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/output";
	if(hasMC){
		file_name.Append("TOF_MC.root");
	}else{
		file_name.Append("TOF.root");
	}

  TFile* inFile = TFile::Open(file_name.Data());
  std::cout << inFile << std::endl;
  if(!inFile){
		if(hasMC){
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTOF_MC.root .");
		}else{
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTOF.root .");
		}
    std::cout << "Copy TOF correction from Alien" << std::endl;
		if(hasMC){
			inFile = TFile::Open("outputTOF_MC.root");
		}else{
			inFile = TFile::Open("outputTOF.root");
		}
  }
  else{
    std::cout << "Correction loaded" << std::endl;
  }
  
	TH3D* mean = dynamic_cast<TH3D*>(inFile->Get("sum_mean_correction"));
	TH3D* width= dynamic_cast<TH3D*>(inFile->Get("sum_width_correction"));
	die->SetCentroidCorrFunctionTOF(mean, corrXdim, corrYdim, corrZdim);
	die->SetWidthCorrFunctionTOF(width, corrXdim, corrYdim, corrZdim);

}


static TH3D LMEECutLib::SetEtaCorrectionTPCTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection){
	

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionTPC() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionTPC()\n";
  std::string file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/outputTPC.root";

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

TH3D LMEECutLib::SetEtaCorrectionITSTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC){
	

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionITSTTree() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionITSTTree()\n";
  std::string file_name;
	if(hasMC){
		file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/outputITS_MC.root";
	}else{
		file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/outputITS.root";
	}

  TFile* recalFile = TFile::Open(file_name.c_str());
  std::cout << recalFile << std::endl;
  if(!recalFile){
		if(hasMC){
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputITS_MC.root .");
		}else{
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputITS.root .");
		}
    std::cout << "Copy ITS correction from Alien" << std::endl;
    if(hasMC){
        recalFile = TFile::Open("outputITS_MC.root");
    }else{
        recalFile = TFile::Open("outputITS.root");
    }
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

static TH3D LMEECutLib::SetEtaCorrectionTOFTTree( Int_t corrXdim, Int_t corrYdim, Int_t corrZdim, Int_t selection, Bool_t hasMC) {
	

  ::Info("LMEECutLib_acapon", " >>>>>>>>>>>>>>>>>>>>>> SetEtaCorrectionTOFTTree() >>>>>>>>>>>>>>>>>>>>>> ");

  std::cout << "starting LMEECutLib::SetEtaCorrectionTOFTTree()\n";
  std::string file_name;
	if(hasMC){
		file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/outputTOF_MC.root";
	}else{
		file_name = "/home/aaron/Data/diElecOutput/PIDcalibration/outputTOF.root";
	}

  TFile* recalFile = TFile::Open(file_name.c_str());
  std::cout << recalFile << std::endl;
  if(!recalFile){
		if(hasMC){
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTOF_MC.root .");
		}else{
			gSystem->Exec("alien_cp alien:///alice/cern.ch/user/a/acapon/PIDcalibration/outputTOF.root .");
		}
    std::cout << "Copy TOF correction from Alien" << std::endl;
    if(hasMC){
        recalFile = TFile::Open("outputTOF_MC.root");
    }else{
        recalFile = TFile::Open("outputTOF.root");
    }
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
			centCuts = new AliDielectronVarCuts("centCutsHigh","cent00100");
			centCuts->AddCut(AliDielectronVarManager::kCentralityNew,0.,100.);
			break;
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
		case kCutSet1:
			mixingHandler = new AliDielectronMixingHandler;
			mixingHandler->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
			mixingHandler->AddVariable(AliDielectronVarManager::kCentralityNew,"0, 10, 20, 30, 40, 60, 80,100");
			mixingHandler->SetDepth(50);
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
	TString weightFile = "alien:///alice/cern.ch/user/a/acapon/TMVAclassifiers/TMVAClassification_BDT.weights.xml";
	AliDielectronTMVACuts* pidCuts = new AliDielectronTMVACuts("PIDCutsTMVA","PIDCutsTMVA");
	pidCuts->AddTMVAInput("pt",          AliDielectronVarManager::kPt);
	pidCuts->AddTMVAInput("EsigTPC",     AliDielectronVarManager::kTPCnSigmaEle);
	pidCuts->AddTMVAInput("EsigITScorr", AliDielectronVarManager::kITSnSigmaEle);
	pidCuts->AddTMVAInput("EsigTOFcorr", AliDielectronVarManager::kTOFnSigmaEle);
	//pidCuts->AddTMVAInput("PsigTPC",   AliDielectronVarManager::kTPCnSigmaPio);
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
		case kCutSet1:{
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);

			Printf("Use TMVA cut value = %f",0.1);
			pidCuts->SetTMVACutValue(0.1);

			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		}
		case kTTreeCuts:
		  // PID cuts used during TTree creating
			// Momentum range relaxed as it cuts on P not Pt. Kinematic cuts applied separately.
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			//cuts->AddCut(cutsPID);
			return cutsPID;
		case kV0_TPCcorr:
			// PID cuts used to select out a very pure sample of V0 electrons using only ITS and TOF
			cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -1., 1., 0.1, 100., kFALSE);
			cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1., 1., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
			cuts->AddCut(cutsPID);
			return cuts;
		case kV0_ITScorr:
			// PID cuts used to select out a very pure sample of V0 electrons using only TPC and TOF
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -0.6, 1., 0.1, 100., kFALSE);
			cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1.,  1., 0.1, 0.4,  kFALSE, AliDielectronPID::kIfAvailable);
			cutsPID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron, -1.,  1., 0.4, 100., kFALSE, AliDielectronPID::kRequire);
			cuts->AddCut(cutsPID);
			return cuts;
		case kV0_TOFcorr:
			// PID cuts used to select out a very pure sample of V0 electrons using only TPC and TOF
			cutsPID->AddCut(AliDielectronPID::kITS, AliPID::kElectron, -1., 1., 0.1, 100., kFALSE);
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -0.6, 1., 0.1, 100., kFALSE);
			cuts->AddCut(cutsPID);
			return cuts;
		case kPdgSel:
			cuts->AddCut(PdgLepton);
			return cuts;
		case kPIDcut1:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[0]);
			pidCuts->SetTMVACutValue(PIDcutRange[0]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut2:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[1]);
			pidCuts->SetTMVACutValue(PIDcutRange[1]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut3:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[2]);
			pidCuts->SetTMVACutValue(PIDcutRange[2]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut4:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[3]);
			pidCuts->SetTMVACutValue(PIDcutRange[3]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut5:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[4]);
			pidCuts->SetTMVACutValue(PIDcutRange[4]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut6:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[5]);
			pidCuts->SetTMVACutValue(PIDcutRange[5]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut7:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[6]);
			pidCuts->SetTMVACutValue(PIDcutRange[6]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut8:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[7]);
			pidCuts->SetTMVACutValue(PIDcutRange[7]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut9:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[8]);
			pidCuts->SetTMVACutValue(PIDcutRange[8]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut10:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[9]);
			pidCuts->SetTMVACutValue(PIDcutRange[9]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut11:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[10]);
			pidCuts->SetTMVACutValue(PIDcutRange[10]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut12:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[11]);
			pidCuts->SetTMVACutValue(PIDcutRange[11]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut13:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[12]);
			pidCuts->SetTMVACutValue(PIDcutRange[12]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut14:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[13]);
			pidCuts->SetTMVACutValue(PIDcutRange[13]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut15:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[14]);
			pidCuts->SetTMVACutValue(PIDcutRange[14]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut16:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[15]);
			pidCuts->SetTMVACutValue(PIDcutRange[15]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut17:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[16]);
			pidCuts->SetTMVACutValue(PIDcutRange[16]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut18:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[17]);
			pidCuts->SetTMVACutValue(PIDcutRange[17]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut19:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[18]);
			pidCuts->SetTMVACutValue(PIDcutRange[18]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kPIDcut20:
			cutsPID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4., 4., 0., 100., kFALSE);
			Printf("Use TMVA cut value = %f", PIDcutRange[19]);
			pidCuts->SetTMVACutValue(PIDcutRange[19]);
			cuts->AddCut(cutsPID);
			cuts->AddCut(pidCuts);
			return cuts;
		case kTheoPID:
			cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
			cutsPID->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
			cutsPID->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
			cutsPID->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
			cuts->AddCut(cutsPID);
			cuts->Print();
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
			break;
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
			break;
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
			}
			varCutsFilter->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0,   4.5);

			//Select filterbit 4
			trackCutsFilter->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);//or 1<<4
			trackCutsFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kFirst);
			//Refits	
			trackCutsFilter->SetRequireITSRefit(kTRUE);
			trackCutsFilter->SetRequireTPCRefit(kTRUE);
	
			trackCuts->AddCut(varCutsFilter);
			trackCuts->AddCut(trackCutsFilter);
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
			trackCuts->Print();
			return trackCuts;
		case kV0_trackCuts:
			// V0 specific track cuts
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
			trackCutsV0->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 4.0);
			trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCr, 100.0, 160.0);
			trackCutsV0->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.8, 1.1);

			trackCuts->AddCut(gammaV0cuts);
			trackCuts->AddCut(trackCutsV0);
			trackCuts->AddCut(GetPIDCuts(PIDcuts));
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
			break;
		case kResolutionTrackCuts:
			varCutsFilter->AddCut(AliDielectronVarManager::kPt, 0.1, 8.0);
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
			break;
		default:
			std::cout << "No Analysis Track Cut defined" << std::endl;
		}
		std::cout << "Track cuts not applied...." << std::endl;
		return 0x0;
}

