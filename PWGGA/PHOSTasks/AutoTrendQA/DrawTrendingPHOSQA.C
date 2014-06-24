void DrawTrendingPHOSQA(TString mergedTrendFile = "trending.root")
{
  //Reads merged trending.root file and draws trending plots.
  
  TFile * fin = TFile::Open(mergedTrendFile.Data());
  if(!fin) { printf("File trending.root not exists.\n"); return; }
  
  TTree * ttree = (TTree*) fin->Get("trending"); 
  if (!ttree) { printf("Trending tree not found."); return; }

  Int_t nRuns = ttree->GetEntries();

  Int_t runNumber = 0; Float_t nEvents = 0;
  Float_t avCluEnergySM1=0., avCluMultSM1=0., avNcellPerCluSM1=0.; // Module 1
  Float_t avCluEnergySM2=0., avCluMultSM2=0., avNcellPerCluSM2=0.; // Module 2
  Float_t avCluEnergySM3=0., avCluMultSM3=0., avNcellPerCluSM3=0.; // Module 3
 
  ttree->SetBranchAddress("run",&runNumber);
  ttree->SetBranchAddress("nEvents",&nEvents);

  ttree->SetBranchAddress("avCluEnergySM1",&avCluEnergySM1);
  ttree->SetBranchAddress("avCluMultSM1",&avCluMultSM1);
  ttree->SetBranchAddress("avNcellPerCluSM1",&avNcellPerCluSM1);

  ttree->SetBranchAddress("avCluEnergySM2",&avCluEnergySM2);
  ttree->SetBranchAddress("avCluMultSM2",&avCluMultSM2);
  ttree->SetBranchAddress("avNcellPerCluSM2",&avNcellPerCluSM2);

  ttree->SetBranchAddress("avCluEnergySM3",&avCluEnergySM3);
  ttree->SetBranchAddress("avCluMultSM3",&avCluMultSM3);
  ttree->SetBranchAddress("avNcellPerCluSM3",&avNcellPerCluSM3);

  // booking histograms
  TH1F* havCluEnergySM1 = new TH1F("havCluEnergySM1","Average cluster energy in Module 4",nRuns,0.,nRuns);
  TH1F* havCluEnergySM2 = new TH1F("havCluEnergySM2","Average cluster energy in Module 3",nRuns,0.,nRuns);
  TH1F* havCluEnergySM3 = new TH1F("havCluEnergySM3","Average cluster energy in Module 2",nRuns,0.,nRuns);

  TH1F* havCluMultSM1 = new TH1F("havCluMultSM1","Average number of clusters per event in Module 4",nRuns,0.,nRuns);
  TH1F* havCluMultSM2 = new TH1F("havCluMultSM2","Average number of clusters per event in Module 3",nRuns,0.,nRuns);
  TH1F* havCluMultSM3 = new TH1F("havCluMultSM3","Average number of clusters per event in Module 2",nRuns,0.,nRuns);

  TH1F* havNcellPerCluSM1 = new TH1F("havNcellPerCluSM1","Average number of cells in cluster in Module 4",nRuns,0.,nRuns);
  TH1F* havNcellPerCluSM2 = new TH1F("havNcellPerCluSM2","Average number of cells in cluster in Module 3",nRuns,0.,nRuns);
  TH1F* havNcellPerCluSM3 = new TH1F("havNcellPerCluSM3","Average number of cells in cluster in Module 2",nRuns,0.,nRuns);

  // List of histograms for saving
  TList list;
  list.Add(havCluEnergySM1);
  list.Add(havCluEnergySM2);
  list.Add(havCluEnergySM3);
  list.Add(havCluMultSM1);
  list.Add(havCluMultSM2);
  list.Add(havCluMultSM3);
  list.Add(havNcellPerCluSM1);
  list.Add(havNcellPerCluSM2);
  list.Add(havNcellPerCluSM3);

  char runlabel[6];
  
  for (Int_t irun=0;irun<nRuns;irun++){
    ttree->GetEntry(irun);
    
    sprintf(runlabel,"%i",runNumber);
    
    havCluEnergySM1->SetBinContent(irun+1, avCluEnergySM1);
    havCluEnergySM1->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    havCluEnergySM2->SetBinContent(irun+1, avCluEnergySM2);
    havCluEnergySM2->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    havCluEnergySM3->SetBinContent(irun+1, avCluEnergySM3);
    havCluEnergySM3->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM1->SetBinContent(irun+1, avCluMultSM1);
    havCluMultSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM2->SetBinContent(irun+1, avCluMultSM2);
    havCluMultSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM3->SetBinContent(irun+1, avCluMultSM3);
    havCluMultSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM1->SetBinContent(irun+1, avNcellPerCluSM1);
    havNcellPerCluSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM2->SetBinContent(irun+1, avNcellPerCluSM2);
    havNcellPerCluSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM3->SetBinContent(irun+1, avNcellPerCluSM3);
    havNcellPerCluSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);	
  }
  
  
  TFile* fout = new TFile("ProductionQA.hist.root","recreate");
  list.Write(); // save selected trend histograms to file
  fout->Close();
}
