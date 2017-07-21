//=====================================================//
//Macro that reads the output of the flow analysis and 
//retrieves the QC and the MH containers.
//The macro produces an output called outputMH.root 
//that contains the four histograms related to the 
//3 particle correlator and the two integrated v2{2,QC} 
//and v2{4,QC}
//Author: Panos.Christakoglou@cern.ch

const Int_t nCentralityBins = 9;
TString strCentralityBins[nCentralityBins] = {"0To5","5To10","10To20",
					      "20To30","30To40","40To50",
					      "50To60","60To70","70To80"};

void plotMHCentrality(const char* filename = "finalAnalysisResults.root",
		      const char* systemType = "PbPb",
		      Bool_t isPID = kFALSE,
		      Int_t charge = -1) {
  gStyle->SetPalette(1,0);

  //----------------------------------------------------------
  // >>>>>>>>>>> Load libraries <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  gSystem->AddIncludePath("-I$ROOTSYS/include");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");

  AliFlowTrackCuts::PIDsource sourcePID = 0x0;
  AliPID::EParticleType particleType = 0x0;

  //----------------------------------------------------------
  // >>>>>>>>>>> Open file - Get objects <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  TFile *fInput = TFile::Open(filename);
  if(!fInput) {
    Printf("File not found!!!");
    break;
  }
  //fInput->ls();

  //Get the TList of the SP
  TList *listSPv2 = GetSPResults(fInput,systemType,2, isPID,
  sourcePID,particleType,charge);
  TList *listSPv4 = GetSPResults(fInput,systemType,4, isPID,
  sourcePID,particleType,charge);

  //Get the TList of the QC
  TList *listQCv2 = GetQCResults(fInput,systemType,2,isPID,
    sourcePID,particleType,charge);
  TList *listQCv4 = GetQCResults(fInput,systemType,4,isPID,
  sourcePID,particleType,charge);

  //Get the TList of the MH
  TList *listMHv1 = GetMHResults(fInput,systemType,1,isPID,
				 sourcePID,particleType,charge);
  TList *listMHv2 = GetMHResults(fInput,systemType,2,isPID,
				 sourcePID,particleType,charge);
  
  TFile *fOutput = TFile::Open("outputMH.root","recreate");
  listSPv2->Write("listSPv2",TObject::kSingleKey);
  listSPv4->Write("listSPv4",TObject::kSingleKey);
  listQCv2->Write("listQCv2",TObject::kSingleKey);
  listQCv4->Write("listQCv4",TObject::kSingleKey);
  listMHv1->Write("listMHv1",TObject::kSingleKey);
  listMHv2->Write("listMHv2",TObject::kSingleKey);
  fOutput->Close();
}

//____________________________________________________________//
TList *GetSPResults(TFile *fInput,
		    const char* systemType = 0x0,
		    Int_t iHarmonic = 2,
		    Bool_t isPID = kTRUE,
		    AliFlowTrackCuts::PIDsource sourcePID = AliFlowTrackCuts::kTOFbayesian,
		    AliPID::EParticleType particleType=AliPID::kPion,
		    Int_t charge = 0) {
  //Function that reads the TDirectoryFile of the MH
  //and returns a TList with the relevant plots.
  TList *listOutput = new TList();
  listOutput->SetName("listSP");

  //______________________________________________________________//
  //Global variables
  AliFlowCommonHistResults *commonHistRes;
  TString gAliFlowCommonHistName;

  //______________________________________________________________//
  //Get the TDirectoryFile
  TString directoryNameSP = 0;
  directoryNameSP = Form("outputSPv%danalysis",iHarmonic); 
  TDirectoryFile *outputSPanalysis = dynamic_cast<TDirectoryFile *>(fInput->Get(directoryNameSP.Data()));
  if(!outputSPanalysis) {
    Printf("SP directory not found!!!");
    break;
  }
  //outputSPanalysis->ls();

  //______________________________________________________________//
  //Get the TList
  TString listNameSP = 0;
  TList *cobjSP;
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++) {
    listNameSP = "Centrality"; 
    listNameSP += strCentralityBins[iCentralityBin];
    listNameSP += "_SP_v"; listNameSP += iHarmonic;
    if(isPID) {
      listNameSP += AliFlowTrackCuts::PIDsourceName(sourcePID);
      listNameSP += "_";
      listNameSP += AliPID::ParticleName(particleType);
    }
    cobjSP = dynamic_cast<TList *>(outputSPanalysis->Get(listNameSP.Data()));
    if(!cobjSP) {
      Printf("SP object list not found!!!");
      break;
    }
    // cobjSP->ls();

    //Common hist for SP
    Double_t nAnalyzedEvents = 0.;
    AliFlowCommonHist *commonHistSP = dynamic_cast<AliFlowCommonHist *>(cobjSP->FindObject("AliFlowCommonHist_SP"));
    if(commonHistSP) {
      TList *clistSPCommonHist = dynamic_cast<TList *>(commonHistSP->GetHistList());
      if(clistSPCommonHist) {
	TH1F *gHistMultRPSP = dynamic_cast<TH1F *>(clistSPCommonHist->FindObject("Control_Flow_MultRP_AliFlowCommonHist_SP"));
	if(gHistMultRPSP)
	  nAnalyzedEvents = gHistMultRPSP->GetEntries();
      }
    }

    //Integrated flow RP
    TH1F *histSPIntegratedVn = 0x0;
    TH1F *histSPPtDifferentialRPVn = 0x0;
    TH1F *histSPPtDifferentialPOIVn = 0x0;
    TH1F *histSPEtaDifferentialRPVn = 0x0;
    TH1F *histSPEtaDifferentialPOIVn = 0x0;
    gAliFlowCommonHistName = "AliFlowCommonHistResults_SP";
    commonHistRes = dynamic_cast<AliFlowCommonHistResults *>(cobjSP->FindObject(gAliFlowCommonHistName.Data()));
    if(commonHistRes) {
      TList *clistSPCommonHistRes = dynamic_cast<TList *>(commonHistRes->GetHistList());
      if(clistSPCommonHistRes) {
	histSPIntegratedVn = dynamic_cast<TH1F *>(clistSPCommonHistRes->FindObject("Flow_Integrated_AliFlowCommonHistResults_SP"));
	//Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[0].Data());
	if(histSPIntegratedVn) 
	  histSPIntegratedVn->SetName("fHistIntegratedFlowRPSP");

	//pT-differential plots
	histSPPtDifferentialRPVn = dynamic_cast<TH1F *>(clistSPCommonHistRes->FindObject("Flow_Differential_Pt_RP_AliFlowCommonHistResults_SP"));
	if(histSPPtDifferentialRPVn) 
	  histSPPtDifferentialRPVn->SetName("fHistPtDifferentialFlowRPSP");

	histSPPtDifferentialPOIVn = dynamic_cast<TH1F *>(clistSPCommonHistRes->FindObject("Flow_Differential_Pt_POI_AliFlowCommonHistResults_SP"));
	if(histSPPtDifferentialPOIVn) 
	  histSPPtDifferentialPOIVn->SetName("fHistPtDifferentialFlowPOISP");

	//eta-differential plots
	histSPEtaDifferentialRPVn = dynamic_cast<TH1F *>(clistSPCommonHistRes->FindObject("Flow_Differential_Eta_RP_AliFlowCommonHistResults_SP"));
	if(histSPEtaDifferentialRPVn) 
	  histSPEtaDifferentialRPVn->SetName("fHistEtaDifferentialFlowRPSP");

	histSPEtaDifferentialPOIVn = dynamic_cast<TH1F *>(clistSPCommonHistRes->FindObject("Flow_Differential_Eta_POI_AliFlowCommonHistResults_SP"));
	if(histSPEtaDifferentialPOIVn) 
	  histSPEtaDifferentialPOIVn->SetName("fHistEtaDifferentialFlowPOISP");
      }
    }

    //______________________________________________________________//
    //Print the results
    Printf("\n============================================================");
    Printf("Analyzed events: %d",(Int_t)nAnalyzedEvents);
    if(histSPIntegratedVn) 
      Printf(Form("v%d(SP): %lf +- %lf",iHarmonic,
		  histSPIntegratedVn->GetBinContent(1),
		  histSPIntegratedVn->GetBinError(1)));
    Printf("============================================================");

    TList *dirOutput = new TList();
    dirOutput->SetName(listNameSP.Data());
    if(histSPIntegratedVn) dirOutput->Add(histSPIntegratedVn);
    if(histSPPtDifferentialRPVn) dirOutput->Add(histSPPtDifferentialRPVn);
    if(histSPPtDifferentialPOIVn) dirOutput->Add(histSPPtDifferentialPOIVn);
    if(histSPEtaDifferentialRPVn) dirOutput->Add(histSPEtaDifferentialRPVn);
    if(histSPEtaDifferentialPOIVn) dirOutput->Add(histSPEtaDifferentialPOIVn);
    listOutput->Add(dirOutput);
  }//loop over centrality bins TLists

  return listOutput;
}

//____________________________________________________________//
TList *GetQCResults(TFile *fInput,
		    const char* systemType = 0x0,
		    Int_t iHarmonic = 2,
		    Bool_t isPID = kTRUE,
		    AliFlowTrackCuts::PIDsource sourcePID = AliFlowTrackCuts::kTOFbayesian,
		    AliPID::EParticleType particleType=AliPID::kPion,
		    Int_t charge = 0) {
  //Function that reads the TDirectoryFile of the QC
  //and returns a TList with the relevant plots.
  TList *listOutput = new TList();
  listOutput->SetName("listQC");

  //______________________________________________________________//
  //Global variables
  const Int_t nQCMethods = 4;
  AliFlowCommonHistResults *commonHistRes[nQCMethods];
  TString methodIntegratedFlowQC[nQCMethods] = {"2ndOrderQC",
						"4thOrderQC",
						"6thOrderQC",
						"8thOrderQC"};
  TString gAliFlowCommonHistName[nQCMethods];

  //______________________________________________________________//
  //Get the TDirectoryFile
  TString directoryNameQC = 0;
  directoryNameQC = Form("outputQCv%danalysis",iHarmonic); 
  TDirectoryFile *outputQCanalysis = dynamic_cast<TDirectoryFile *>(fInput->Get(directoryNameQC.Data()));
  if(!outputQCanalysis) {
    Printf("QC directory not found!!!");
    break;
  }
  //outputQCanalysis->ls();

  //______________________________________________________________//
  //Get the TList
  TString listNameQC = 0;
  TList *cobjQC;
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++) {
    listNameQC = "Centrality"; 
    listNameQC += strCentralityBins[iCentralityBin];
    listNameQC += "_QC_v"; listNameQC += iHarmonic;
    if(isPID) {
      listNameQC += AliFlowTrackCuts::PIDsourceName(sourcePID);
      listNameQC += "_";
      listNameQC += AliPID::ParticleName(particleType);
    }

    cobjQC = dynamic_cast<TList *>(outputQCanalysis->Get(listNameQC.Data()));
    if(!cobjQC) {
      Printf("QC object list not found!!!");
      break;
    }
    //cobjQC->ls();

    //Common hist for QC
    Double_t nAnalyzedEvents = 0.;
    AliFlowCommonHist *commonHistQC = dynamic_cast<AliFlowCommonHist *>(cobjQC->FindObject("AliFlowCommonHistQC"));
    if(commonHistQC) {
      TList *clistQCCommonHist = dynamic_cast<TList *>(commonHistQC->GetHistList());
      if(clistQCCommonHist) {
	TH1F *gHistMultRPQC = dynamic_cast<TH1F *>(clistQCCommonHist->FindObject("Control_Flow_MultRP_AliFlowCommonHistQC"));
	if(gHistMultRPQC)
	  nAnalyzedEvents = gHistMultRPQC->GetEntries();
      }
    }

    //2nd order QC
    gAliFlowCommonHistName[0] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[0] += methodIntegratedFlowQC[0].Data();
    commonHistRes[0] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[0].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[0].Data());
    TH1D *histQC2 = commonHistRes[0]->GetHistIntFlowRP();
    histQC2->SetName("fHistIntegratedFlowRPQC2");
    TH1D *histPtDifferentialQC2RP = commonHistRes[0]->GetHistDiffFlowPtRP();
    histPtDifferentialQC2RP->SetName("fHistPtDifferentialFlowRPQC2");
    TH1D *histPtDifferentialQC2POI = commonHistRes[0]->GetHistDiffFlowPtPOI();
    histPtDifferentialQC2POI->SetName("fHistPtDifferentialFlowPOIQC2");
    TH1D *histEtaDifferentialQC2RP = commonHistRes[0]->GetHistDiffFlowEtaRP();
    histEtaDifferentialQC2RP->SetName("fHistEtaDifferentialFlowRPQC2");
    TH1D *histEtaDifferentialQC2POI = commonHistRes[0]->GetHistDiffFlowEtaPOI();
    histEtaDifferentialQC2POI->SetName("fHistEtaDifferentialFlowPOIQC2");
    
    //4th order QC
    gAliFlowCommonHistName[1] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[1] += methodIntegratedFlowQC[1].Data();
    commonHistRes[1] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[1].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[1].Data());
    TH1D *histQC4 = commonHistRes[1]->GetHistIntFlowRP();
    histQC4->SetName("fHistIntegratedFlowRPQC4");
    TH1D *histPtDifferentialQC4RP = commonHistRes[1]->GetHistDiffFlowPtRP();
    histPtDifferentialQC4RP->SetName("fHistPtDifferentialFlowRPQC4");
    TH1D *histPtDifferentialQC4POI = commonHistRes[1]->GetHistDiffFlowPtPOI();
    histPtDifferentialQC4POI->SetName("fHistPtDifferentialFlowPOIQC4");
    TH1D *histEtaDifferentialQC4RP = commonHistRes[1]->GetHistDiffFlowEtaRP();
    histEtaDifferentialQC4RP->SetName("fHistEtaDifferentialFlowRPQC4");
    TH1D *histEtaDifferentialQC4POI = commonHistRes[1]->GetHistDiffFlowEtaPOI();
    histEtaDifferentialQC4POI->SetName("fHistEtaDifferentialFlowPOIQC4");
    
    //6th order QC
    gAliFlowCommonHistName[2] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[2] += methodIntegratedFlowQC[2].Data();
    commonHistRes[2] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[2].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[2].Data());
    TH1D *histQC6 = commonHistRes[2]->GetHistIntFlowRP();
    histQC6->SetName("fHistIntegratedFlowRPQC6");
    TH1D *histPtDifferentialQC6RP = commonHistRes[2]->GetHistDiffFlowPtRP();
    histPtDifferentialQC6RP->SetName("fHistPtDifferentialFlowRPQC6");
    TH1D *histPtDifferentialQC6POI = commonHistRes[2]->GetHistDiffFlowPtPOI();
    histPtDifferentialQC6POI->SetName("fHistPtDifferentialFlowPOIQC6");
    TH1D *histEtaDifferentialQC6RP = commonHistRes[2]->GetHistDiffFlowEtaRP();
    histEtaDifferentialQC6RP->SetName("fHistEtaDifferentialFlowRPQC6");
    TH1D *histEtaDifferentialQC6POI = commonHistRes[2]->GetHistDiffFlowEtaPOI();
    histEtaDifferentialQC6POI->SetName("fHistEtaDifferentialFlowPOIQC6");

    //8th order QC
    gAliFlowCommonHistName[3] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[3] += methodIntegratedFlowQC[3].Data();
    commonHistRes[3] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[3].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[3].Data());
    TH1D *histQC8 = commonHistRes[3]->GetHistIntFlowRP();
    histQC8->SetName("fHistIntegratedFlowRPQC8");
    TH1D *histPtDifferentialQC8RP = commonHistRes[3]->GetHistDiffFlowPtRP();
    histPtDifferentialQC8RP->SetName("fHistPtDifferentialFlowRPQC8");
    TH1D *histPtDifferentialQC8POI = commonHistRes[3]->GetHistDiffFlowPtPOI();
    histPtDifferentialQC8POI->SetName("fHistPtDifferentialFlowPOIQC8");
    TH1D *histEtaDifferentialQC8RP = commonHistRes[3]->GetHistDiffFlowEtaRP();
    histEtaDifferentialQC8RP->SetName("fHistEtaDifferentialFlowRPQC8");
    TH1D *histEtaDifferentialQC8POI = commonHistRes[3]->GetHistDiffFlowEtaPOI();
    histEtaDifferentialQC8POI->SetName("fHistEtaDifferentialFlowPOIQC8");

    TString gHistEventsName = Form("gHistEventsv%d",iHarmonic);
    gHistEventsName += strCentralityBins[iCentralityBin];
    TH1D *gHistEvents = new TH1D(gHistEventsName.Data(),
				 ";;N_{analyzed events}",
				 1,0.5,1.5);
    gHistEvents->SetBinContent(1,nAnalyzedEvents);

    //______________________________________________________________//
    //Print the results
    Printf("\n============================================================");
    Printf("Analyzed events: %d",(Int_t)nAnalyzedEvents);
    Printf(Form("v%d(2,QC): %lf +- %lf",iHarmonic,
		histQC2->GetBinContent(1),
		histQC2->GetBinError(1)));
    Printf(Form("v%d(4,QC): %lf +- %lf",iHarmonic,
		histQC4->GetBinContent(1),
		histQC4->GetBinError(1)));
    Printf(Form("v%d(6,QC): %lf +- %lf",iHarmonic,
		histQC6->GetBinContent(1),
		histQC6->GetBinError(1)));
    Printf(Form("v%d(8,QC): %lf +- %lf",iHarmonic,
		histQC8->GetBinContent(1),
		histQC8->GetBinError(1)));
    Printf("============================================================");

    TList *dirOutput = new TList();
    dirOutput->SetName(listNameQC.Data());
    dirOutput->Add(histQC2);
    dirOutput->Add(histQC4);
    dirOutput->Add(histQC6);
    dirOutput->Add(histQC8);
    dirOutput->Add(histPtDifferentialQC2RP);
    dirOutput->Add(histPtDifferentialQC2POI);
    dirOutput->Add(histEtaDifferentialQC2RP);
    dirOutput->Add(histEtaDifferentialQC2POI);
    dirOutput->Add(histPtDifferentialQC4RP);
    dirOutput->Add(histPtDifferentialQC4POI);
    dirOutput->Add(histEtaDifferentialQC4RP);
    dirOutput->Add(histEtaDifferentialQC4POI);
    dirOutput->Add(histPtDifferentialQC6RP);
    dirOutput->Add(histPtDifferentialQC6POI);
    dirOutput->Add(histEtaDifferentialQC6RP);
    dirOutput->Add(histEtaDifferentialQC6POI);
    dirOutput->Add(histPtDifferentialQC8RP);
    dirOutput->Add(histPtDifferentialQC8POI);
    dirOutput->Add(histEtaDifferentialQC8RP);
    dirOutput->Add(histEtaDifferentialQC8POI);
    dirOutput->Add(gHistEvents);
    listOutput->Add(dirOutput);
  }//loop over centrality bins TLists

  return listOutput;
}

//____________________________________________________________//
TList *GetMHResults(TFile *fInput,
		    const char* systemType = 0x0,
		    Int_t iHarmonic = 1,
		    Bool_t isPID = kTRUE,
		    AliFlowTrackCuts::PIDsource sourcePID = AliFlowTrackCuts::kTOFbayesian,
		    AliPID::EParticleType particleType=AliPID::kPion,
		    Int_t charge = 0) {
  //Function that reads the TDirectoryFile of the MH
  //and returns a TList with the relevant plots.
  TList *listOutput = new TList();
  listOutput->SetName("listMH");

  //______________________________________________________________//
  //Get the TDirectoryFile
  TString directoryNameMH = 0;
  if(charge == 0) {
    if(iHarmonic == 1)
      directoryNameMH = "outputMHUSanalysis"; 
    else if(iHarmonic == 2)
      directoryNameMH = Form("outputMHUS%danalysis",iHarmonic); 
  }
  else if((charge == -1)||(charge == 1)) {
    if(iHarmonic == 1)
      directoryNameMH = "outputMHLSanalysis"; 
    else if(iHarmonic == 2)
      directoryNameMH = Form("outputMHLS%danalysis",iHarmonic); 
  }

  TDirectoryFile *outputMHanalysis = dynamic_cast<TDirectoryFile *>(fInput->Get(directoryNameMH.Data()));
  if(!outputMHanalysis) {
    Printf("MH directory not found!!!");
    break;
  }
  //outputMHanalysis->ls();

  //______________________________________________________________//
  //Get the TList
  TString listNameMH = 0;
  TList *cobjMH;
  TCanvas *c[nCentralityBins];
  for(Int_t iCentralityBin = 0; iCentralityBin < nCentralityBins; iCentralityBin++) {
    listNameMH = "Centrality"; 
    listNameMH += strCentralityBins[iCentralityBin];
    if(charge == 0) 
      listNameMH += "_MHUS_v"; 
    else if((charge == -1)||(charge == 1)) 
      listNameMH += "_MHLS_v"; 
    listNameMH += iHarmonic;
    if(isPID) {
      listNameMH += AliFlowTrackCuts::PIDsourceName(sourcePID);
      listNameMH += "_";
      listNameMH += AliPID::ParticleName(particleType);
    }

    //Printf("%s",listNameMH.Data());
    cobjMH = dynamic_cast<TList *>(outputMHanalysis->Get(listNameMH.Data()));
    if(!cobjMH) {
      Printf("MH object list not found!!!");
      break;
    }
    //cobjMH->ls();
    
    //______________________________________________________________//
    //Get the daughter TLists
    TList *listWeights = dynamic_cast<TList *>(cobjMH->FindObject("Weights"));
    //listWeights->ls();
    TList *listProfiles = dynamic_cast<TList *>(cobjMH->FindObject("Profiles"));
    //listProfiles->ls();
    TList *listResults = dynamic_cast<TList *>(cobjMH->FindObject("Results"));
    //listResults->ls();

    if((!listWeights)||(!listProfiles)||(!listResults)) {
      Printf("MH output lists not found!!!");
      break;
    }
    
    //______________________________________________________________//
    TProfile *fAnalysisSettings = dynamic_cast<TProfile *>(cobjMH->FindObject("fAnalysisSettings"));

    //______________________________________________________________//
    //Get the objects from the Results list
    TH1D *f3pCorrelatorHist = dynamic_cast<TH1D *>(listResults->FindObject("f3pCorrelatorHist"));
    TH1D *fDetectorBiasHist = dynamic_cast<TH1D *>(listResults->FindObject("fDetectorBiasHist"));
    TH1D *f3pCorrelatorVsMHist = dynamic_cast<TH1D *>(listResults->FindObject("f3pCorrelatorVsMHist"));
    TH1D *fDetectorBiasVsMHist = dynamic_cast<TH1D *>(listResults->FindObject("fDetectorBiasVsMHist"));
    
    //______________________________________________________________//
    //Get the objects from the Profile list
    TProfile *f3pCorrelatorPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorPro"));

    //3-particle differential
    TProfile *f3pCorrelatorVsPtSumPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsPtSumPro"));
    TProfile *f3pCorrelatorVsPtDiffPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsPtDiffPro"));
    TProfile *f3pCorrelatorVsEtaSumPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsEtaSumPro"));
    TProfile *f3pCorrelatorVsEtaDiffPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsEtaDiffPro"));

    //2-particle differential: pT
    TProfile *f2pCorrelatorCosPsiDiffPtDiff = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiDiffPtDiff"));
    TProfile *f2pCorrelatorCosPsiSumPtDiff = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiSumPtDiff"));
    TProfile *f2pCorrelatorSinPsiSumPtDiff = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorSinPsiSumPtDiff"));
    TProfile *f2pCorrelatorCosPsiDiffPtSum = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiDiffPtSum"));
    TProfile *f2pCorrelatorCosPsiSumPtSum = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiSumPtSum"));
    TProfile *f2pCorrelatorSinPsiSumPtSum = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorSinPsiSumPtSum"));

    //2-particle differential: eta
    TProfile *f2pCorrelatorCosPsiDiffEtaDiff = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiDiffEtaDiff"));
    TProfile *f2pCorrelatorCosPsiSumEtaDiff = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiSumEtaDiff"));
    TProfile *f2pCorrelatorSinPsiSumEtaDiff = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorSinPsiSumEtaDiff"));
    TProfile *f2pCorrelatorCosPsiDiffEtaSum = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiDiffEtaSum"));
    TProfile *f2pCorrelatorCosPsiSumEtaSum = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorCosPsiSumEtaSum"));
    TProfile *f2pCorrelatorSinPsiSumEtaSum = dynamic_cast<TProfile *>(listProfiles->FindObject("f2pCorrelatorSinPsiSumEtaSum"));

    //integrated 3-particle correlator
    Double_t g3pCorrelatorValue = 0., g3pCorrelatorError = 0.;
    GetCorrelatorAndError(f3pCorrelatorVsPtSumPro,
			  g3pCorrelatorValue,
			  g3pCorrelatorError,1,17);
    TString g3pHistName = "g3pHist";
    g3pHistName += strCentralityBins[iCentralityBin];
    TH1D *g3pHist = new TH1D(g3pHistName,
			     ";;#LT#LTcos(#psi_{1}+#psi_{2}-2#phi_{3}#GT#GT",
			     1,0.5,1.5);
    g3pHist->SetBinContent(1,g3pCorrelatorValue);
    g3pHist->SetBinError(1,g3pCorrelatorError);

    //______________________________________________________________//
    //Draw the differential 3p correlator
    c[iCentralityBin] = new TCanvas(listNameMH.Data(),listNameMH.Data(),
				    iCentralityBin*50,
				    iCentralityBin*50,
				    500,500);
    c[iCentralityBin]->SetHighLightColor(10); 
    c[iCentralityBin]->SetFillColor(10);
    f3pCorrelatorVsPtDiffPro->DrawCopy("E");

    //if(iCentralityBin == 8) {
    //for(Int_t iBin = 1; iBin <= f3pCorrelatorVsPtDiffPro->GetNbinsX(); iBin++) {
	//if(f3pCorrelatorVsPtDiffPro->GetBinContent(iBin) != 0.) 
	//Printf("Entries: %lf - Error: %lf - Value: %lf",
	//f3pCorrelatorVsPtDiffPro->GetBinEntries(iBin),
	//f3pCorrelatorVsPtDiffPro->GetBinError(iBin),
	//f3pCorrelatorVsPtDiffPro->GetBinContent(iBin));
	//}
    //}

    //______________________________________________________________//
    //Print the results
    Printf("============================================================");
    Printf("=========================Bin: %s=========================",strCentralityBins[iCentralityBin].Data());
    cout<<"<cos(psi1 + psi2 - 2phi3)>: "<<
      g3pCorrelatorValue <<
      " +- " <<
      g3pCorrelatorError << endl;
    cout<<"<cos(phi1 + phi2 - 2phi3)>: "<<
      f3pCorrelatorPro->GetBinContent(1) <<
      " +- " <<
      f3pCorrelatorPro->GetBinError(1)<<endl;
    Printf("============================================================");

    TList *dirOutput = new TList();
    dirOutput->SetName(listNameMH.Data());
    dirOutput->Add(f3pCorrelatorPro);
    dirOutput->Add(f3pCorrelatorVsPtSumPro);
    dirOutput->Add(f3pCorrelatorVsPtDiffPro);
    dirOutput->Add(f3pCorrelatorVsEtaSumPro);
    dirOutput->Add(f3pCorrelatorVsEtaDiffPro);
    dirOutput->Add(g3pHist);

    dirOutput->Add(f2pCorrelatorCosPsiDiffPtDiff);
    dirOutput->Add(f2pCorrelatorCosPsiSumPtDiff);
    dirOutput->Add(f2pCorrelatorSinPsiSumPtDiff);
    dirOutput->Add(f2pCorrelatorCosPsiDiffPtSum);
    dirOutput->Add(f2pCorrelatorCosPsiSumPtSum);
    dirOutput->Add(f2pCorrelatorSinPsiSumPtSum);

    dirOutput->Add(f2pCorrelatorCosPsiDiffEtaDiff);
    dirOutput->Add(f2pCorrelatorCosPsiSumEtaDiff);
    dirOutput->Add(f2pCorrelatorSinPsiSumEtaDiff);
    dirOutput->Add(f2pCorrelatorCosPsiDiffEtaSum);
    dirOutput->Add(f2pCorrelatorCosPsiSumEtaSum);
    dirOutput->Add(f2pCorrelatorSinPsiSumEtaSum);
    listOutput->Add(dirOutput);
  }//loop over centrality bins TLists

  //listOutput->ls();
  return listOutput;
}

//____________________________________________________________//
void GetCorrelatorAndError(TProfile *g3pCorrelatorVsPt,
			   Double_t &g3pCorrelatorValue,
			   Double_t &g3pCorrelatorError,
			   Int_t iBinLow = 0,
			   Int_t iBinHigh = 0) {
  //Function to return the average value of the 3p correlator 
  //<cos(psi1 + psi2 - 2phi3)> and its error.
  //The first argument is one of the 3p TProfile objects vs pt.
  //The second and third argument give the integral and its error.
  //The fourth and fifth, if specified, indicate the lowest and 
  //highest bin the calculation should be performed for.
  Int_t gBinLow = 1, gBinHigh = g3pCorrelatorVsPt->GetNbinsX();
  if(iBinLow) gBinLow = iBinLow;
  if(iBinHigh) gBinHigh = iBinHigh;
  
  Double_t gSumXi = 0.;
  Double_t gSumYi = 0.;
  Double_t gSumXiYi = 0.;
  Double_t gSumXiYi2 = 0.;
  Double_t gSumXi2Yi2 = 0.;
  Double_t gSumDeltaXi2 = 0.;
  Double_t gSumYi2DeltaXi2 = 0.;
  Double_t dError = 0.; //Flow code driven error calculation

  Double_t kSumBi = 0., kSumBi2DeltaYi2 = 0.;
  Double_t kSumYiBi = 0., kSumYi2DeltaBi2 = 0.;
  Double_t kSumDeltaBi2 = 0.;

 for(Int_t iBin = gBinLow; iBin <= gBinHigh; iBin++) {
    gSumXi += g3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumYi += g3pCorrelatorVsPt->GetBinContent(iBin);
    gSumXiYi += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin);
    gSumXiYi2 += g3pCorrelatorVsPt->GetBinEntries(iBin)*TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2);
    gSumXi2Yi2 += TMath::Power(g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin),2);
    gSumDeltaXi2 += TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
    gSumYi2DeltaXi2 += TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2) + TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);

    dError += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinError(iBin)*g3pCorrelatorVsPt->GetBinError(iBin);  

    //new error calculation
    kSumBi += g3pCorrelatorVsPt->GetBinEntries(iBin);
    kSumYiBi += g3pCorrelatorVsPt->GetBinEntries(iBin)*g3pCorrelatorVsPt->GetBinContent(iBin);
    kSumBi2DeltaYi2 += TMath::Power(g3pCorrelatorVsPt->GetBinEntries(iBin),2)*TMath::Power(g3pCorrelatorVsPt->GetBinError(iBin),2);
    kSumYi2DeltaBi2 += TMath::Power(g3pCorrelatorVsPt->GetBinContent(iBin),2)*TMath::Power(TMath::Sqrt(g3pCorrelatorVsPt->GetBinEntries(iBin)),2);  
    kSumDeltaBi2 += TMath::Power(TMath::Sqrt(g3pCorrelatorVsPt->GetBinEntries(iBin)),2);
  }
  
  g3pCorrelatorValue = -1000.;
  g3pCorrelatorError = 1000.;
  
  if(gSumXi != 0.)
    g3pCorrelatorValue = gSumXiYi/gSumXi;
  if((gSumXi != 0.)&&(gSumXiYi != 0.))
    g3pCorrelatorError = TMath::Abs((gSumXiYi/gSumXi))*TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumYi2DeltaXi2)/gSumXiYi),2) + TMath::Power((gSumDeltaXi2/gSumXi),2));
    //g3pCorrelatorError = TMath::Sqrt((1./TMath::Power(kSumBi,2))*(kSumBi2DeltaYi2 + kSumYi2DeltaBi2 + TMath::Power(kSumYiBi,2)*kSumDeltaBi2/TMath::Power(kSumBi,2)));
  if(gSumXi != 0.)
    dError /= TMath::Power(gSumXi,2);
  dError = TMath::Sqrt(dError);
  g3pCorrelatorError = dError;

  return;
}

