//=====================================================//
//Macro that reads the output of the flow analysis and 
//retrieves the QC and the MH containers.
//The macro produces an output called outputMH.root 
//that contains the four histograms related to the 
//3 particle correlator and the two integrated v2{2,QC} 
//and v2{4,QC}
//Author: Panos.Christakoglou@cern.ch

void plotMH(const char* filename = "AnalysisResults.root",
	    const char* analysisType = "TPC only",
	    Int_t centrality = 0) {
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
  gSystem->Load("libPWGflowBase");

  //----------------------------------------------------------
  // >>>>>>>>>>> Open file - Get objects <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  TFile *fInput = TFile::Open(filename);
  if(!fInput) {
    Printf("File not found!!!");
    break;
  }
  //f->ls();

  //Get the TList of the MH
  TList *listMH = GetMHResults(fInput,analysisType,centrality);

  //Get the TList of the QC
  TList *listQC = GetQCResults(fInput,analysisType,centrality);

  TFile *fOutput = TFile::Open("outputMH.root","recreate");
  listMH->Write();
  listQC->Write();
  fOutput->Close();
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
  if(centrality > -1) {
    listNameQC = "cobjQC_outputCentrality"; 
    listNameQC += centrality; listNameQC += ".root";
  }
  else listNameQC = "cobjQC";
  TList *cobjQC = dynamic_cast<TList *>(outputQCanalysis->Get(listNameQC.Data()));
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

  TH1D *gHistEvents = new TH1D("gHistEvents",";;N_{analyzed events}",
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

  listOutput->Add(histQC_2);
  listOutput->Add(histQC_4);
  listOutput->Add(histQC_6);
  listOutput->Add(histQC_8);
  listOutput->Add(gHistEvents);

  return listOutput;
}

//____________________________________________________________//
TList *GetMHResults(TFile *fInput,
		    const char* analysisType = 0x0,
		    Int_t centrality = -1) {
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
  if(centrality > -1) {
    listNameMH = "cobjMH_outputCentrality"; 
    listNameMH += centrality; listNameMH += ".root";
  }
  else listNameMH = "cobjMH";
  TList *cobjMH = dynamic_cast<TList *>(outputMHanalysis->Get(listNameMH.Data()));
  if(!cobjMH) {
    Printf("MH object list not found!!!");
    break;
  }
  //cobjMH->ls();

  //______________________________________________________________//
  //Get the daughter TLists
  TList *listWeights = dynamic_cast<TList *>(cobjMH->At(0));
  //listWeights->ls();
  TList *listProfiles = dynamic_cast<TList *>(cobjMH->At(1));
  //listProfiles->ls();
  TList *listResults = dynamic_cast<TList *>(cobjMH->At(2));
  //listResults->ls();

  if((!listWeights)||(!listProfiles)||(!listResults)) {
    Printf("MH output lists not found!!!");
    break;
  }

  //______________________________________________________________//
  TProfile *fAnalysisSettings = dynamic_cast<TProfile *>(cobjMH->At(3));

  //______________________________________________________________//
  //Get the objects from the Results list
  TH1D *f3pCorrelatorHist = dynamic_cast<TH1D *>(listResults->At(0));
  TH1D *fDetectorBiasHist = dynamic_cast<TH1D *>(listResults->At(1));
  TH1D *f3pCorrelatorVsMHist = dynamic_cast<TH1D *>(listResults->At(2));
  TH1D *fDetectorBiasVsMHist = dynamic_cast<TH1D *>(listResults->At(3));

  //______________________________________________________________//
  //Get the objects from the Profile list
  TProfile *f3pCorrelatorPro = dynamic_cast<TProfile *>(listProfiles->At(0));
  TProfile *fNonIsotropicTermsPro = dynamic_cast<TProfile *>(listProfiles->At(1));
  TProfile *f3pCorrelatorVsMPro = dynamic_cast<TProfile *>(listProfiles->At(2));
  TProfile2D *fNonIsotropicTermsVsMPro = dynamic_cast<TProfile2D *>(listProfiles->At(3));
  TProfile *f3pCorrelatorVsPtSumPro = dynamic_cast<TProfile *>(listProfiles->At(4));
  TProfile *f3pCorrelatorVsPtDiffPro = dynamic_cast<TProfile *>(listProfiles->At(5));
  
  //______________________________________________________________//
  //Draw the histograms
  TCanvas *c0 = new TCanvas("c0","Analysis settings",0,0,500,500);
  c0->SetHighLightColor(10); c0->SetFillColor(10);
  fAnalysisSettings->Draw();

  TCanvas *c1 = new TCanvas("c1","3p correlator",50,50,500,500);
  c1->SetHighLightColor(10); c1->SetFillColor(10);
  f3pCorrelatorHist->Draw("E");

  TCanvas *c2 = new TCanvas("c2","Detector bias",100,100,500,500);
  c2->SetHighLightColor(10); c2->SetFillColor(10);
  fDetectorBiasHist->Draw("E");

  TCanvas *c3 = new TCanvas("c3","3p correlator vs M",150,150,500,500);
  c3->SetHighLightColor(10); c3->SetFillColor(10);
  f3pCorrelatorVsMHist->Draw("E");

  TCanvas *c4 = new TCanvas("c4","Detector bias vs M",200,200,500,500);
  c4->SetHighLightColor(10); c4->SetFillColor(10);
  fDetectorBiasVsMHist->Draw("E");

  TCanvas *c5 = new TCanvas("c5","3p correlator (pro)",500,0,500,500);
  c5->SetHighLightColor(10); c5->SetFillColor(10);
  f3pCorrelatorPro->Draw("E");

  TCanvas *c6 = new TCanvas("c6","Non isotropic terms",550,50,500,500);
  c6->SetHighLightColor(10); c6->SetFillColor(10);
  fNonIsotropicTermsPro->Draw("E");

  TCanvas *c7 = new TCanvas("c7","3p correlator vs M (pro)",600,100,500,500);
  c7->SetHighLightColor(10); c7->SetFillColor(10);
  f3pCorrelatorVsMPro->Draw("E");

  TCanvas *c8 = new TCanvas("c8","Non isotropic terms vs M (pro)",650,150,500,500);
  c8->SetHighLightColor(10); c8->SetFillColor(10);
  fNonIsotropicTermsVsMPro->Draw("colz");

  TCanvas *c9 = new TCanvas("c9","3p correlator vs sum pt",700,200,500,500);
  c9->SetHighLightColor(10); c9->SetFillColor(10);
  f3pCorrelatorVsPtSumPro->Draw("E");

  TCanvas *c10 = new TCanvas("c10","3p correlator vs diff pt",750,250,500,500);
  c10->SetHighLightColor(10); c10->SetFillColor(10);
  f3pCorrelatorVsPtDiffPro->Draw("E");

  Double_t g3pCorrelatorValue = 0., g3pCorrelatorError = 0.;
  GetCorrelatorAndError(f3pCorrelatorVsPtSumPro,
			//GetCorrelatorAndError(f3pCorrelatorVsPtDiffPro,
			g3pCorrelatorValue,
			g3pCorrelatorError,1,20);
			
  //______________________________________________________________//
  //Print the results
  Printf("============================================================");
  cout<<"<cos(psi1 + psi2 - 2phi3)>: "<<
    g3pCorrelatorValue <<
    " +- " <<
    g3pCorrelatorError << endl;
  cout<<"<cos(phi1 + phi2 - 2phi3)>: "<<
    f3pCorrelatorHist->GetBinContent(1) <<
    " +- " <<
    f3pCorrelatorHist->GetBinError(1)<<endl;
  Printf("============================================================");

  listOutput->Add(f3pCorrelatorHist);
  listOutput->Add(f3pCorrelatorVsMHist);
  listOutput->Add(f3pCorrelatorVsPtSumPro);
  listOutput->Add(f3pCorrelatorVsPtDiffPro);

  return listOutput;
}

//____________________________________________________________//
void GetCorrelatorAndError(TProfile *f3pCorrelatorVsPt,
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
  Int_t gBinLow = 1, gBinHigh = f3pCorrelatorVsPt->GetNbinsX();
  if(iBinLow) gBinLow = iBinLow;
  if(iBinHigh) gBinHigh = iBinHigh;
  
  Int_t iBinCounter = 0;
  Double_t gSumBinContentTimesWeight = 0., gSumWeight = 0.;
  Double_t gSumBinContentTimesWeightSquared = 0.;
  for(Int_t iBin = gBinLow; iBin <= gBinHigh; iBin++) {
    iBinCounter += 1;
    
    gSumBinContentTimesWeight += f3pCorrelatorVsPt->GetBinContent(iBin)*f3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumWeight += f3pCorrelatorVsPt->GetBinEntries(iBin);
    gSumBinContentTimesWeightSquared += TMath::Power(gSumBinContentTimesWeight,2);
  }

  Printf("%lf - %d",gSumWeight,iBinCounter);
  //Calculate the g3pCorrelatorValue and its error
  g3pCorrelatorValue = -1000.;
  g3pCorrelatorError = 1000.;
  if((gSumWeight)&&(iBinCounter)) {
    g3pCorrelatorValue = gSumBinContentTimesWeight/(gSumWeight*iBinCounter);
    g3pCorrelatorError = TMath::Sqrt(gSumBinContentTimesWeightSquared/TMath::Power(gSumWeight,2))/iBinCounter;
  }

  return;
}

