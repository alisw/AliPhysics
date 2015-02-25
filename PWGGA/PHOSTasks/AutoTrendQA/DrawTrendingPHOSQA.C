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
  Float_t avCluEnergySM4=0., avCluMultSM4=0., avNcellPerCluSM4=0.; // Module 4

  Float_t avPi0NumSM1=0., avPi0MassSM1=0., avPi0SigmaSM1=0.;
  Float_t avPi0NumErrSM1=0., avPi0MassErrSM1=0., avPi0SigmaErrSM1=0.;

  Float_t avPi0NumSM2=0., avPi0MassSM2=0., avPi0SigmaSM2=0.;
  Float_t avPi0NumErrSM2=0., avPi0MassErrSM2=0., avPi0SigmaErrSM2=0.;

  Float_t avPi0NumSM3=0., avPi0MassSM3=0., avPi0SigmaSM3=0.;
  Float_t avPi0NumErrSM3=0., avPi0MassErrSM3=0., avPi0SigmaErrSM3=0.;

  Float_t avPi0NumSM4=0., avPi0MassSM4=0., avPi0SigmaSM4=0.;
  Float_t avPi0NumErrSM4=0., avPi0MassErrSM4=0., avPi0SigmaErrSM4=0.;

  //clusters
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

  ttree->SetBranchAddress("avCluEnergySM4",&avCluEnergySM4);
  ttree->SetBranchAddress("avCluMultSM4",&avCluMultSM4);
  ttree->SetBranchAddress("avNcellPerCluSM4",&avNcellPerCluSM4);

  //pi0
  ttree->SetBranchAddress("avPi0NumSM1",&avPi0NumSM1);
  ttree->SetBranchAddress("avPi0NumSM2",&avPi0NumSM2);
  ttree->SetBranchAddress("avPi0NumSM3",&avPi0NumSM3);
  ttree->SetBranchAddress("avPi0NumSM4",&avPi0NumSM4);
  ttree->SetBranchAddress("avPi0NumErrSM1",&avPi0NumErrSM2);
  ttree->SetBranchAddress("avPi0NumErrSM2",&avPi0NumErrSM2);
  ttree->SetBranchAddress("avPi0NumErrSM3",&avPi0NumErrSM3);
  ttree->SetBranchAddress("avPi0NumErrSM4",&avPi0NumErrSM4);

  ttree->SetBranchAddress("avPi0MassSM1",&avPi0MassSM1);
  ttree->SetBranchAddress("avPi0MassSM2",&avPi0MassSM2);
  ttree->SetBranchAddress("avPi0MassSM3",&avPi0MassSM3);
  ttree->SetBranchAddress("avPi0MassSM4",&avPi0MassSM4);
  ttree->SetBranchAddress("avPi0MassErrSM1",&avPi0MassErrSM2);
  ttree->SetBranchAddress("avPi0MassErrSM2",&avPi0MassErrSM2);
  ttree->SetBranchAddress("avPi0MassErrSM3",&avPi0MassErrSM3);
  ttree->SetBranchAddress("avPi0MassErrSM4",&avPi0MassErrSM4);

  ttree->SetBranchAddress("avPi0SigmaSM1",&avPi0SigmaSM1);
  ttree->SetBranchAddress("avPi0SigmaSM2",&avPi0SigmaSM2);
  ttree->SetBranchAddress("avPi0SigmaSM3",&avPi0SigmaSM3);
  ttree->SetBranchAddress("avPi0SigmaSM4",&avPi0SigmaSM4);
  ttree->SetBranchAddress("avPi0SigmaErrSM1",&avPi0SigmaErrSM2);
  ttree->SetBranchAddress("avPi0SigmaErrSM2",&avPi0SigmaErrSM2);
  ttree->SetBranchAddress("avPi0SigmaErrSM3",&avPi0SigmaErrSM3);
  ttree->SetBranchAddress("avPi0SigmaErrSM4",&avPi0SigmaErrSM4);

  // booking histograms
  TH1F* havCluEnergySM1 = new TH1F("havCluEnergySM1","Average cluster energy in Module 4",nRuns,0.,nRuns);
  TH1F* havCluEnergySM2 = new TH1F("havCluEnergySM2","Average cluster energy in Module 3",nRuns,0.,nRuns);
  TH1F* havCluEnergySM3 = new TH1F("havCluEnergySM3","Average cluster energy in Module 2",nRuns,0.,nRuns);
  TH1F* havCluEnergySM4 = new TH1F("havCluEnergySM4","Average cluster energy in Module 1",nRuns,0.,nRuns);

  TH1F* havCluMultSM1 = new TH1F("havCluMultSM1","Average number of clusters per event in Module 4",nRuns,0.,nRuns);
  TH1F* havCluMultSM2 = new TH1F("havCluMultSM2","Average number of clusters per event in Module 3",nRuns,0.,nRuns);
  TH1F* havCluMultSM3 = new TH1F("havCluMultSM3","Average number of clusters per event in Module 2",nRuns,0.,nRuns);
  TH1F* havCluMultSM4 = new TH1F("havCluMultSM4","Average number of clusters per event in Module 1",nRuns,0.,nRuns);

  TH1F* havNcellPerCluSM1 = new TH1F("havNcellPerCluSM1","Average number of cells in cluster in Module 4",nRuns,0.,nRuns);
  TH1F* havNcellPerCluSM2 = new TH1F("havNcellPerCluSM2","Average number of cells in cluster in Module 3",nRuns,0.,nRuns);
  TH1F* havNcellPerCluSM3 = new TH1F("havNcellPerCluSM3","Average number of cells in cluster in Module 2",nRuns,0.,nRuns);
  TH1F* havNcellPerCluSM4 = new TH1F("havNcellPerCluSM4","Average number of cells in cluster in Module 1",nRuns,0.,nRuns);

  TH1F* havPi0NumSM1 = new TH1F("havPi0NumSM1","Average number of #pi^{0}s per event in Module 4",nRuns,0.,nRuns);
  TH1F* havPi0NumSM2 = new TH1F("havPi0NumSM2","Average number of #pi^{0}s per event in Module 3",nRuns,0.,nRuns);
  TH1F* havPi0NumSM3 = new TH1F("havPi0NumSM3","Average number of #pi^{0}s per event in Module 2",nRuns,0.,nRuns);
  TH1F* havPi0NumSM4 = new TH1F("havPi0NumSM4","Average number of #pi^{0}s per event in Module 1",nRuns,0.,nRuns);

  TH1F* havPi0MassSM1 = new TH1F("avPi0MassSM1","#pi^{0} mass position in Module 4",nRuns,0.,nRuns);
  TH1F* havPi0MassSM2 = new TH1F("avPi0MassSM2","#pi^{0} mass position in Module 3",nRuns,0.,nRuns);
  TH1F* havPi0MassSM3 = new TH1F("avPi0MassSM3","#pi^{0} mass position in Module 2",nRuns,0.,nRuns);
  TH1F* havPi0MassSM4 = new TH1F("avPi0MassSM4","#pi^{0} mass position in Module 1",nRuns,0.,nRuns);

  TH1F* havPi0SigmaSM1 = new TH1F("avPi0SigmaSM1","#pi^{0} peak width in Module 4",nRuns,0.,nRuns);
  TH1F* havPi0SigmaSM2 = new TH1F("avPi0SigmaSM2","#pi^{0} peak width in Module 3",nRuns,0.,nRuns);
  TH1F* havPi0SigmaSM3 = new TH1F("avPi0SigmaSM3","#pi^{0} peak width in Module 2",nRuns,0.,nRuns);
  TH1F* havPi0SigmaSM4 = new TH1F("avPi0SigmaSM4","#pi^{0} peak width in Module 1",nRuns,0.,nRuns);

  // List of histograms for saving
  TList list;
  list.Add(havCluEnergySM1);
  list.Add(havCluEnergySM2);
  list.Add(havCluEnergySM3);
  list.Add(havCluEnergySM4);
  list.Add(havCluMultSM1);
  list.Add(havCluMultSM2);
  list.Add(havCluMultSM3);
  list.Add(havCluMultSM4);
  list.Add(havNcellPerCluSM1);
  list.Add(havNcellPerCluSM2);
  list.Add(havNcellPerCluSM3);
  list.Add(havNcellPerCluSM4);
  list.Add(havPi0NumSM1);
  list.Add(havPi0NumSM2);
  list.Add(havPi0NumSM3);
  list.Add(havPi0NumSM4);
  list.Add(havPi0MassSM1);
  list.Add(havPi0MassSM2);
  list.Add(havPi0MassSM3);
  list.Add(havPi0MassSM4);
  list.Add(havPi0SigmaSM1);
  list.Add(havPi0SigmaSM2);
  list.Add(havPi0SigmaSM3);
  list.Add(havPi0SigmaSM4);

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

    havCluEnergySM4->SetBinContent(irun+1, avCluEnergySM4);
    havCluEnergySM4->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM1->SetBinContent(irun+1, avCluMultSM1);
    havCluMultSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM2->SetBinContent(irun+1, avCluMultSM2);
    havCluMultSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM3->SetBinContent(irun+1, avCluMultSM3);
    havCluMultSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havCluMultSM4->SetBinContent(irun+1, avCluMultSM4);
    havCluMultSM4->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM1->SetBinContent(irun+1, avNcellPerCluSM1);
    havNcellPerCluSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM2->SetBinContent(irun+1, avNcellPerCluSM2);
    havNcellPerCluSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM3->SetBinContent(irun+1, avNcellPerCluSM3);
    havNcellPerCluSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havNcellPerCluSM4->SetBinContent(irun+1, avNcellPerCluSM4);
    havNcellPerCluSM4->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0NumSM1->SetBinContent(irun+1,avPi0NumSM1);
    havPi0NumSM1->SetBinError(irun+1,avPi0NumErrSM1);
    havPi0NumSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0NumSM2->SetBinContent(irun+1,avPi0NumSM2);
    havPi0NumSM2->SetBinError(irun+1,avPi0NumErrSM2);
    havPi0NumSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0NumSM3->SetBinContent(irun+1,avPi0NumSM3);
    havPi0NumSM3->SetBinError(irun+1,avPi0NumErrSM3);
    havPi0NumSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0NumSM4->SetBinContent(irun+1,avPi0NumSM4);
    havPi0NumSM4->SetBinError(irun+1,avPi0NumErrSM4);
    havPi0NumSM4->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0MassSM1->SetBinContent(irun+1,avPi0MassSM1);
    havPi0MassSM1->SetBinError(irun+1,avPi0MassErrSM1);
    havPi0MassSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0MassSM2->SetBinContent(irun+1,avPi0MassSM2);
    havPi0MassSM2->SetBinError(irun+1,avPi0MassErrSM2);
    havPi0MassSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0MassSM3->SetBinContent(irun+1,avPi0MassSM3);
    havPi0MassSM3->SetBinError(irun+1,avPi0MassErrSM3);
    havPi0MassSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0MassSM4->SetBinContent(irun+1,avPi0MassSM4);
    havPi0MassSM4->SetBinError(irun+1,avPi0MassErrSM4);
    havPi0MassSM4->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0SigmaSM1->SetBinContent(irun+1,avPi0SigmaSM1);
    havPi0SigmaSM1->SetBinError(irun+1,avPi0SigmaErrSM1);
    havPi0SigmaSM1->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0SigmaSM2->SetBinContent(irun+1,avPi0SigmaSM2);
    havPi0SigmaSM2->SetBinError(irun+1,avPi0SigmaErrSM2);
    havPi0SigmaSM2->GetXaxis()->SetBinLabel(irun+1,runlabel);

    havPi0SigmaSM3->SetBinContent(irun+1,avPi0SigmaSM3);
    havPi0SigmaSM3->SetBinError(irun+1,avPi0SigmaErrSM3);
    havPi0SigmaSM3->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    havPi0SigmaSM4->SetBinContent(irun+1,avPi0SigmaSM4);
    havPi0SigmaSM4->SetBinError(irun+1,avPi0SigmaErrSM4);
    havPi0SigmaSM4->GetXaxis()->SetBinLabel(irun+1,runlabel);

  }
  
  
  TFile* fout = new TFile("ProductionQA.hist.root","recreate");
  list.Write(); // save selected trend histograms to file
  fout->Close();
}
