//=====================================================//
//Macro that reads the output of the flow analysis and 
//retrieves the QC and the MH containers.
//The macro produces an output called outputMH.root 
//that contains the four histograms related to the 
//3 particle correlator and the two integrated v2{2,QC} 
//and v2{4,QC}
//Author: Panos.Christakoglou@cern.ch

const Int_t nCentralityBins = 9;
TString strCentralityBins[nCentralityBins] = {"0-5","5-10","10-20",
					      "20-30","30-40","40-50",
					      "50-60","60-70","70-80"};

void plotMHCentrality(const char* filename = "AnalysisResults.root",
		      const char* analysisType = "TPC only") {
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
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG2flowCommon");

  //----------------------------------------------------------
  // >>>>>>>>>>> Open file - Get objects <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  TFile *fInput = TFile::Open(filename);
  if(!fInput) {
    Printf("File not found!!!");
    break;
  }
  //fInput->ls();

  //Get the TList of the MH
  TList *listMH = GetMHResults(fInput,analysisType);

  //Get the TList of the QC
  TList *listQC = GetQCResults(fInput,analysisType);

  //Get the TList of the SP
  TList *listSP = GetSPResults(fInput,analysisType);

  TFile *fOutput = TFile::Open("outputMH.root","recreate");
  listMH->Write();
  listQC->Write();
  listSP->Write();
  fOutput->Close();
}

//____________________________________________________________//
TList *GetSPResults(TFile *fInput,
		    const char* analysisType = 0x0,
		    Int_t centrality = -1) {
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
  directoryNameSP = "outputSPanalysis"; 
  if(analysisType) directoryNameSP += analysisType;
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
    //for(Int_t iCentralityBin = 0; iCentralityBin < 1; iCentralityBin++) {
    listNameSP = "cobjSP_"; listNameSP += strCentralityBins[iCentralityBin];

    cobjSP = dynamic_cast<TList *>(outputSPanalysis->Get(listNameSP.Data()));
    if(!cobjSP) {
      Printf("SP object list not found!!!");
      break;
    }
    
    //Common hist for SP
    Double_t nAnalyzedEvents = 0.;
    AliFlowCommonHist *commonHistSP = dynamic_cast<AliFlowCommonHist *>(cobjSP->FindObject("AliFlowCommonHistSP"));
    if(commonHistSP) {
      TList *clistSPCommonHist = dynamic_cast<TList *>(commonHistSP->GetHistList());
      if(clistSPCommonHist) {
	TH1F *gHistMultRPSP = dynamic_cast<TH1F *>(clistSPCommonHist->FindObject("Control_Flow_MultRP_AliFlowCommonHistSP"));
	if(gHistMultRPSP)
	  nAnalyzedEvents = gHistMultRPSP->GetEntries();
      }
    }

    //Integrated flow RP
    gAliFlowCommonHistName = "AliFlowCommonHistResultsSP";
    commonHistRes = dynamic_cast<AliFlowCommonHistResults *>(cobjSP->FindObject(gAliFlowCommonHistName.Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[0].Data());
    TH1D *histSP = commonHistRes->GetHistIntFlowRP();
    histSP->SetName("fHistIntegratedFlowRPSP");
    
    //______________________________________________________________//
    //Print the results
    Printf("\n============================================================");
    Printf("Analyzed events: %d",(Int_t)nAnalyzedEvents);
    Printf("v2(SP): %lf +- %lf",histSP->GetBinContent(1),
	   histSP->GetBinError(1));
    Printf("============================================================");

    TList *dirOutput = new TList();
    dirOutput->SetName(listNameSP.Data());
    dirOutput->Add(histSP);
    listOutput->Add(dirOutput);
  }//loop over centrality bins TLists

  return listOutput;
}

//____________________________________________________________//
TList *GetQCResults(TFile *fInput,
		    const char* analysisType = 0x0,
		    Int_t centrality = -1) {
  //Function that reads the TDirectoryFile of the MH
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
  directoryNameQC = "outputQCanalysis"; 
  if(analysisType) directoryNameQC += analysisType;
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
    //for(Int_t iCentralityBin = 0; iCentralityBin < 1; iCentralityBin++) {
    listNameQC = "cobjQC_"; listNameQC += strCentralityBins[iCentralityBin];

    cobjQC = dynamic_cast<TList *>(outputQCanalysis->Get(listNameQC.Data()));
    if(!cobjQC) {
      Printf("QC object list not found!!!");
      break;
    }
    
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
    TH1D *histQC_2 = commonHistRes[0]->GetHistIntFlowRP();
    histQC_2->SetName("fHistIntegratedFlowRPQC_2");
    
    //4th order QC
    gAliFlowCommonHistName[1] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[1] += methodIntegratedFlowQC[1].Data();
    commonHistRes[1] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[1].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[1].Data());
    TH1D *histQC_4 = commonHistRes[1]->GetHistIntFlowRP();
    histQC_4->SetName("fHistIntegratedFlowRPQC_4");
    
    //6th order QC
    gAliFlowCommonHistName[2] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[2] += methodIntegratedFlowQC[2].Data();
    commonHistRes[2] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[2].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[2].Data());
    TH1D *histQC_6 = commonHistRes[2]->GetHistIntFlowRP();
    histQC_6->SetName("fHistIntegratedFlowRPQC_6");
    
    //8th order QC
    gAliFlowCommonHistName[3] = "AliFlowCommonHistResults";
    gAliFlowCommonHistName[3] += methodIntegratedFlowQC[3].Data();
    commonHistRes[3] = dynamic_cast<AliFlowCommonHistResults *>(cobjQC->FindObject(gAliFlowCommonHistName[3].Data()));
    //Printf("AliFlowCommonHist name %s",gAliFlowCommonHistName[3].Data());
    TH1D *histQC_8 = commonHistRes[3]->GetHistIntFlowRP();
    histQC_8->SetName("fHistIntegratedFlowRPQC_8");
    
    TString gHistEventsName = "gHistEvents";
    gHistEventsName += strCentralityBins[iCentralityBin];
    TH1D *gHistEvents = new TH1D(gHistEventsName.Data(),
				 ";;N_{analyzed events}",
				 1,0.5,1.5);
    gHistEvents->SetBinContent(1,nAnalyzedEvents);

    //______________________________________________________________//
    //Print the results
    Printf("\n============================================================");
    Printf("Analyzed events: %d",(Int_t)nAnalyzedEvents);
    Printf("v2(2,QC): %lf +- %lf",histQC_2->GetBinContent(1),
	   histQC_2->GetBinError(1));
    Printf("v2(4,QC): %lf +- %lf",histQC_4->GetBinContent(1),
	   histQC_4->GetBinError(1));
    Printf("v2(6,QC): %lf +- %lf",histQC_6->GetBinContent(1),
	   histQC_6->GetBinError(1));
    Printf("v2(8,QC): %lf +- %lf",histQC_8->GetBinContent(1),
	   histQC_8->GetBinError(1));
    Printf("============================================================");

    TList *dirOutput = new TList();
    dirOutput->SetName(listNameQC.Data());
    dirOutput->Add(histQC_2);
    dirOutput->Add(histQC_4);
    dirOutput->Add(histQC_6);
    dirOutput->Add(histQC_8);
    dirOutput->Add(gHistEvents);
    listOutput->Add(dirOutput);
  }//loop over centrality bins TLists

  return listOutput;
}

//____________________________________________________________//
TList *GetMHResults(TFile *fInput,
		    const char* analysisType = 0x0) {
  //Function that reads the TDirectoryFile of the MH
  //and returns a TList with the relevant plots.
  TList *listOutput = new TList();
  listOutput->SetName("listMH");

  //______________________________________________________________//
  //Get the TDirectoryFile
  TString directoryNameMH = 0;
  directoryNameMH = "outputMHanalysis"; 
  if(analysisType) directoryNameMH += analysisType;
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
    //for(Int_t iCentralityBin = 0; iCentralityBin < 1; iCentralityBin++) {
    listNameMH = "cobjMH_"; listNameMH += strCentralityBins[iCentralityBin];
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
    TProfile *fNonIsotropicTermsPro = dynamic_cast<TProfile *>(listProfiles->FindObject("fNonIsotropicTermsPro"));
    TProfile *f3pCorrelatorVsMPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsMPro"));
    TProfile2D *fNonIsotropicTermsVsMPro = dynamic_cast<TProfile2D *>(listProfiles->FindObject("fNonIsotropicTermsVsMPro"));
    TProfile *f3pCorrelatorVsPtSumPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsPtSumPro"));
    TProfile *f3pCorrelatorVsPtDiffPro = dynamic_cast<TProfile *>(listProfiles->FindObject("f3pCorrelatorVsPtDiffPro"));
    
    Double_t g3pCorrelatorValue = 0., g3pCorrelatorError = 0.;
    GetCorrelatorAndError(f3pCorrelatorVsPtSumPro,
			  g3pCorrelatorValue,
			  g3pCorrelatorError,1,17);
    TString g3pHistName = "g3pHistName";
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
      f3pCorrelatorHist->GetBinContent(1) <<
      " +- " <<
      f3pCorrelatorHist->GetBinError(1)<<endl;
    Printf("============================================================");

    TList *dirOutput = new TList();
    dirOutput->SetName(listNameMH.Data());
    dirOutput->Add(f3pCorrelatorHist);
    //dirOutput->Add(f3pCorrelatorVsMHist);
    dirOutput->Add(f3pCorrelatorVsPtSumPro);
    dirOutput->Add(f3pCorrelatorVsPtDiffPro);
    dirOutput->Add(g3pHist);
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

