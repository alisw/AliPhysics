void readBalanceFunction() {
  //Macro to read the output of the BF analysis:
  //i) Merges the output of each sub-job
  //ii) Prints and draws the final output
  //iii) Printd the details of the BF analysis (i.e. type, bins,...)
  //iv) Reads the QA part of the analysis
  //Author: Panos.Christakoglou@cern.ch
  //Loading the needed libraries
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2ebye.so");


  //Draw BF              
  drawBF();

  //Print BF config
  //printBFConfig();     

  //Merge the output
  //mergeOutput("/alice/cern.ch/user/p/pchrist/Balance/pp/7TeV/LHC10b/output/");
}

//___________________________________________________________//
void drawBF() {
  //Function to draw the BF object
  TFile *f = TFile::Open("AnalysisResults.118558.root");
  if(!f) {
    Printf("File not found!!!");
    break;
  }

  TDirectoryFile *dir = dynamic_cast<TDirectoryFile *>(f->Get("PWG2EbyE.outputBalanceFunctionAnalysis.root"));
  if(!dir) {
    Printf("Output directory not found!!!");
    break;
  }

  AliBalance *bf = dynamic_cast<AliBalance *>(dir->Get("AliBalance"));
  if(!bf) {
    Printf("BF object not found!!!");
    break;
  }

  bf->PrintResults();
  TGraphErrors *gr = bf->DrawBalance();
  gr->SetMarkerStyle(24);

  TCanvas *c1 = new TCanvas("c1","Balance Function",0,0,500,500);
  c1->SetFillColor(10); c1->SetHighLightColor(10);

  TH2F *hEmpty = new TH2F("hEmpty","",1000,-1000,1000,1000,-0.2,1.0);
  hEmpty->SetStats(kFALSE);
  hEmpty->GetXaxis()->SetTitle(gr->GetXaxis()->GetTitle());
  hEmpty->GetYaxis()->SetTitle(gr->GetYaxis()->GetTitle());
  hEmpty->GetXaxis()->SetRangeUser(gr->GetXaxis()->GetXmin(),
                                   gr->GetXaxis()->GetXmax());
  hEmpty->GetYaxis()->SetRangeUser(gr->GetYaxis()->GetXmin(),
                                   gr->GetYaxis()->GetXmax());
  hEmpty->Draw();
  gr->Draw("P");
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
