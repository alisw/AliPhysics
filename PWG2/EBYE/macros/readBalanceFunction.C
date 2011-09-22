const TString gBFAnalysisType[7] = {"y","eta","qlong","qout","qside","qinv","phi"};

void readBalanceFunction(Bool_t bHistos = kTRUE, TString inFile = "AnalysisResults.root") {
  // Macro to read the output of the BF analysis:  MW: CHANGE THIS!!!!
  //i) Prints and draws the final BF output
  //ii) Plots the QA part of the analysis
  //iii) store BF in output file
  //Author: Panos.Christakoglou@cern.ch, m.weber@cern.ch
  //Loading the needed libraries
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2ebye.so");

  //Draw BF              
  drawBF(bHistos,inFile);

  //Merge the output
  //mergeOutput("/alice/cern.ch/user/p/pchrist/Balance/pp/7TeV/LHC10b/output/");
}

//___________________________________________________________//
void drawBF(Bool_t bHistos = kTRUE, TString inFile = "AnalysisResults.root") {
  //Function to draw the BF objects and write them into the output file

  Int_t maximumCanvases = 10;
  Int_t iCanvas         = 0;
  TCanvas *cQA[10];
  TCanvas *cBF[10];
  TCanvas *cBFS[10];

  // get the file
  TFile *f = TFile::Open(inFile.Data());
  if(!f) {
    Printf("File not found!!!");
    break;
  }

  // get the BF output directory
  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWG2EbyE.outputBalanceFunctionAnalysis"));
  if(!dir) {
    Printf("Output directory not found!!!");
    break;
  }

  // loop over all lists and plot the BF and QA
  TList *list = NULL;
  TString listName;
  TIter nextkey( dir->GetListOfKeys() );
  TKey *key;

  AliBalance *bf[10][7];
  AliBalance *bfs[10][7];
  TGraphErrors *gbf[10][7];
  TGraphErrors *gbfs[10][7];

  TH1D *fHistP[7]; //N+
  TH1D *fHistN[7]; //N-
  TH1D *fHistPN[7]; //N+-
  TH1D *fHistNP[7]; //N-+
  TH1D *fHistPP[7]; //N++
  TH1D *fHistNN[7]; //N--

  while ( (key = (TKey*)nextkey())) {

    list = (TList*)key->ReadObj();
    listName = TString(list->GetName());

    cout<<"Processing list "<<listName<<endl;

    // ----------------------------------------------------
    // plot QA histograms
    if(listName.Contains("QA")){     

      iCanvas ++;
      if(iCanvas == 10) {cout<<"TOO MANY LISTS --> increase MAXIMUM"<<endl; return;}
      cQA[iCanvas] = new TCanvas(listName,listName);
      cQA[iCanvas]->Divide(3,3);

      cQA[iCanvas]->cd(1);
      TH1F* histVx = (TH1F*)list->FindObject("fHistVx");
      if(histVx){
	histVx->SetFillColor(9);
	histVx->Draw();
      }
      cQA[iCanvas]->cd(2);
      TH1F* histVy = (TH1F*)list->FindObject("fHistVy");
      if(histVy){
	histVy->SetFillColor(9);
	histVy->Draw();
      }
      
      cQA[iCanvas]->cd(3);
      TH1F* histVz = (TH1F*)list->FindObject("fHistVz");
      if(histVz){
	histVz->SetFillColor(9);
	histVz->Draw();
      }
      
      cQA[iCanvas]->cd(4);
      TH1F* histEventStats = (TH1F*)list->FindObject("fHistEventStats");
      if(histEventStats){
	histEventStats->SetFillColor(9);
	histEventStats->Draw();
      }
      
      cQA[iCanvas]->cd(5);
      cQA[iCanvas]->cd(5)->SetLogz();
      TH2F* histClus = (TH2F*)list->FindObject("fHistClus");
      if(histClus) histClus->Draw("colz");
      
      cQA[iCanvas]->cd(6);
      cQA[iCanvas]->cd(6)->SetLogz();
      TH2F* histDCA = (TH2F*)list->FindObject("fHistDCA");  
      if(histDCA) histDCA->Draw("colz");
      
      cQA[iCanvas]->cd(7);
      TH1F* histPt = (TH1F*)list->FindObject("fHistPt");
      if(histPt){
	histPt->SetFillColor(9);
	histPt->Draw();
      }
      
      cQA[iCanvas]->cd(8);
      TH1F* histEta = (TH1F*)list->FindObject("fHistEta");
      if(histEta){
	histEta->SetFillColor(9);
	histEta->Draw();
      }
      
      cQA[iCanvas]->cd(9);
      TH1F* histPhi = (TH1F*)list->FindObject("fHistPhi");
      if(histPhi){
	histPhi->SetFillColor(9);
	histPhi->Draw();
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF 
    if(listName.Contains("BF_")){
      cBF[iCanvas] = new TCanvas(listName,listName);
      cBF[iCanvas]->Divide(3,3);

      for(Int_t a = 0; a < 7; a++){

	cout<<"ANALYSE "<<gBFAnalysisType[a]<<endl;

	// create the BF object
	bf[iCanvas][a]  = new AliBalance();

	fHistP[a] = (TH1D*)list->FindObject(Form("fHistP%s",gBFAnalysisType[a].Data()));
	fHistN[a] = (TH1D*)list->FindObject(Form("fHistP%s",gBFAnalysisType[a].Data()));
	fHistPP[a] = (TH1D*)list->FindObject(Form("fHistPP%s",gBFAnalysisType[a].Data()));
	fHistPN[a] = (TH1D*)list->FindObject(Form("fHistPN%s",gBFAnalysisType[a].Data()));
	fHistNP[a] = (TH1D*)list->FindObject(Form("fHistNP%s",gBFAnalysisType[a].Data()));
	fHistNN[a] = (TH1D*)list->FindObject(Form("fHistNN%s",gBFAnalysisType[a].Data()));

	// set the binning (p1 doesn't play a role --> 0)
	bf[iCanvas][a]->SetNumberOfBins(a,fHistNN[a]->GetNbinsX());
	bf[iCanvas][a]->SetInterval(0,0,a,fHistNN[a]->GetBinCenter(1) - fHistNN[a]->GetBinWidth(1)/2,fHistNN[a]->GetBinCenter(fHistNN[a]->GetNbinsX()) + fHistNN[a]->GetBinWidth(1)/2);

	// set the members for each bin in histogram
	for(Int_t ibin = 1; ibin <= fHistNN[a]->GetNbinsX(); ibin++){
	  bf[iCanvas][a]->SetNpp(a,ibin-1,fHistPP[a]->GetBinContent(ibin));
	  bf[iCanvas][a]->SetNpn(a,ibin-1,fHistPN[a]->GetBinContent(ibin));
	  bf[iCanvas][a]->SetNnp(a,ibin-1,fHistNP[a]->GetBinContent(ibin));
	  bf[iCanvas][a]->SetNnn(a,ibin-1,fHistNN[a]->GetBinContent(ibin));
	}
	bf[iCanvas][a]->SetNp(a,fHistP[a]->GetEntries());
	bf[iCanvas][a]->SetNn(a,fHistN[a]->GetEntries());

	gbf[iCanvas][a] = bf[iCanvas][a]->drawBalance(a);
	gbf[iCanvas][a]->SetName(Form("%s_BF_%s",listName.Data(),gBFAnalysisType[a].Data()));

	cBF[iCanvas]->cd(a+1);
	gbf[iCanvas][a]->SetMarkerStyle(20);

	if(!bHistos){
	  gbf[iCanvas][a]->Draw("AP");
	}
	else{
	  fHistPP[a]->Draw();
	  fHistPN[a]->SetLineColor(2);
	  fHistPN[a]->Draw("same");
	  fHistNP[a]->SetLineColor(4);
	  fHistNP[a]->Draw("same");
	  fHistNN[a]->SetLineColor(8);
	  fHistNN[a]->Draw("same");
	}
      }
    }
    // ----------------------------------------------------

    // ----------------------------------------------------
    // calculate and plot BF (shuffled)
    if(listName.Contains("BFShuffled")){
  
      cBFS[iCanvas] = new TCanvas(listName,listName);
      cBFS[iCanvas]->Divide(3,3);

      for(Int_t a = 0; a < 7; a++){

	// create the BF object
	bfs[iCanvas][a]  = new AliBalance();

	fHistP[a] = (TH1D*)list->FindObject(Form("fHistP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistN[a] = (TH1D*)list->FindObject(Form("fHistP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistPP[a] = (TH1D*)list->FindObject(Form("fHistPP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistPN[a] = (TH1D*)list->FindObject(Form("fHistPN%s_shuffle",gBFAnalysisType[a].Data()));
	fHistNP[a] = (TH1D*)list->FindObject(Form("fHistNP%s_shuffle",gBFAnalysisType[a].Data()));
	fHistNN[a] = (TH1D*)list->FindObject(Form("fHistNN%s_shuffle",gBFAnalysisType[a].Data()));

	// set the binning (p1 doesn't play a role --> 0)
	bfs[iCanvas][a]->SetNumberOfBins(a,fHistNN[a]->GetNbinsX());
	bfs[iCanvas][a]->SetInterval(0,0,a,fHistNN[a]->GetBinCenter(1) - fHistNN[a]->GetBinWidth(1)/2,fHistNN[a]->GetBinCenter(fHistNN[a]->GetNbinsX()) + fHistNN[a]->GetBinWidth(1)/2);

	// set the members for each bin in histogram
	for(Int_t ibin = 1; ibin <= fHistNN[a]->GetNbinsX(); ibin++){
	  bfs[iCanvas][a]->SetNpp(a,ibin-1,fHistPP[a]->GetBinContent(ibin));
	  bfs[iCanvas][a]->SetNpn(a,ibin-1,fHistPN[a]->GetBinContent(ibin));
	  bfs[iCanvas][a]->SetNnp(a,ibin-1,fHistNP[a]->GetBinContent(ibin));
	  bfs[iCanvas][a]->SetNnn(a,ibin-1,fHistNN[a]->GetBinContent(ibin));
	}
	bfs[iCanvas][a]->SetNp(a,fHistP[a]->GetEntries());
	bfs[iCanvas][a]->SetNn(a,fHistN[a]->GetEntries());

	gbfs[iCanvas][a] = bf[iCanvas][a]->drawBalance(a);
	gbfs[iCanvas][a]->SetName(Form("%s_BF_%s",listName.Data(),gBFAnalysisType[a].Data()));

	cBFS[iCanvas]->cd(a+1);
	gbfs[iCanvas][a]->SetMarkerStyle(20);
	if(!bHistos){
	  gbf[iCanvas][a]->Draw("AP");
	}
	else{
	  fHistPP[a]->Draw();
	  fHistPN[a]->SetLineColor(2);
	  fHistPN[a]->Draw("same");
	  fHistNP[a]->SetLineColor(4);
	  fHistNP[a]->Draw("same");
	  fHistNN[a]->SetLineColor(8);
	  fHistNN[a]->Draw("same");
	}
      }
    }
    // ----------------------------------------------------
  }

  TFile *fOut = TFile::Open(Form("Histograms_%s",inFile.Data()),"RECREATE");
  fOut->cd();
  for(Int_t i = 0; i < iCanvas; i++){
    for(Int_t a = 0; a < 7; a++){
      gbf[iCanvas][a]->Write();
      gbfs[iCanvas][a]->Write();
    }
  }
}

//___________________________________________________________//
void mergeOutput(const char* outputDir) {
  //Function to merge the output of the sub-jobs
  //Create a BF object
  AliBalance *bf = new AliBalance();

  //connect to AliEn's API services
  TGrid::Connect("alien://"); 

  //Getting the output dir from the env. variable 
  //(JDL field: JDLVariables={"OutputDir"};)
  TGridResult* result = gGrid->Query(outputDir,"*/root_archive.zip","","-l 1000");
  
  Int_t nEntries = result->GetEntries();

  TString alienUrl;
  TDirectoryFile *dirSubJob;

  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  TH1F *fHistEventStats = new TH1F("fHistEventStats",
				   "Event statistics;;N_{events}",
				   4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());

  AliESDtrackCuts *bfTrackCuts = new AliESDtrackCuts("bfTrackCuts");
  for(Int_t i = 0; i < nEntries; i++) {
    alienUrl = result->GetKey(i,"turl");
    alienUrl += "#AnalysisResults.root";
    Printf("Opening file: %s",alienUrl.Data());
    TFile *file = TFile::Open(alienUrl.Data());
    dirSubJob = dynamic_cast<TDirectoryFile *>(file->Get("PWG2EbyE.outputBalanceFunctionAnalysis.root"));

    //merge BF
    AliBalance *bfSubJob = dynamic_cast<AliBalance *>(dirSubJob->Get("AliBalance"));
    //bfSubJob->PrintResults();
    bf->Merge(bfSubJob);
    //delete bfSubJob;

    //merge event stats
    TList *listSubJob = dynamic_cast<TList *>(dirSubJob->Get("listQA"));
    fHistEventStats->Add(dynamic_cast<TH1F *>(listSubJob->At(0)));

    bfTrackCuts = dynamic_cast<AliESDtrackCuts *>(listSubJob->At(1));
    delete listSubJob;
  }

  //Create the output file
  TString outputFile = "AnalysisResults.Merged.root";
  TFile *foutput = TFile::Open(outputFile.Data(),"recreate");
  TDirectoryFile *dirOutput = new TDirectoryFile();
  dirOutput->SetName("PWG2EbyE.outputBalanceFunctionAnalysis.root");
  //dirOutput->cd();
  dirOutput->Add(bf);
  TList *list = new TList();
  list->SetName("listQA");
  list->Add(fHistEventStats);
  list->Add(bfTrackCuts);
  dirOutput->Add(list);
  dirOutput->Write();
  bf->Write();
  list->Write();
  foutput->Close();

    //cout<<alienUrl.Data()<<endl;
}
