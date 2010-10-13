void readBalanceFunctionF() {
  //Macro to read the output of the BF analysis:
  //i) Merges the output of each sub-job
  //ii) Prints and draws the final output
  //iii) Printd the details of the BF analysis (i.e. type, bins,...)
  //iv) Reads the QA part of the analysis
  //iv) 
  //Author: Panos.Christakoglou@cern.ch

  //Merge the output
  mergeOutput("/alice/cern.ch/user/p/pchrist/Balance/pp/7TeV/LHC10b/output/");
}

//___________________________________________________________//
void mergeOutput(const char* outputDir) {
  //Function to merge the output of the sub-jobs
 //Loading the needed libraries
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2ebye.so");

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
