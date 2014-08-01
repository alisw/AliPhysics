void MakeTrendingPHOSQA(const char* file="QAresults.root", Int_t runNumber, Bool_t IsOnGrid = kFALSE)
{
  //const char* dirNames[] = {"PHOSCellsQA_AnyInt","PHOSCellsQA_PHI7","PHOSPbPbQAResults","PHOSTriggerQAResults"};
  if (IsOnGrid) TGrid::Connect("alien://");  

  TFile* fin = TFile::Open(file);
  if(!fin) {printf("Cannot open file %s. exit.\n",file); return; }

  const char* listNameAnyInt = "PHOSCellsQA_AnyInt";
  TObjArray* listAnyInt = (TObjArray*)fin->Get(listNameAnyInt);

  const char* listNamePHI7 = "PHOSCellsQA_PHI7";
  TObjArray* listPHI7 = (TObjArray*)fin->Get(listNamePHI7);

  const char* listNamePbPb = "PHOSPbPbQAResults";
  TList* listPbPb = (TList*)fin->Get(listNamePbPb);

  const char* listNameTrig = "PHOSTriggerQAResults";
  TList* listTrig = (TList*)fin->Get(listNameTrig);
  
  
  TFile * trendFile = new TFile("trending.root","recreate");
  TTree * ttree=new TTree("trending","tree of trending variables");
 
  Int_t   nEvents=0;
  Float_t avCluEnergySM1=-9999., avCluMultSM1=-9999., avNcellPerCluSM1=-9999.; // Module 1
  Float_t avCluEnergySM2=-9999., avCluMultSM2=-9999., avNcellPerCluSM2=-9999.; // Module 2
  Float_t avCluEnergySM3=-9999., avCluMultSM3=-9999., avNcellPerCluSM3=-9999.; // Module 3

  ttree->Branch("run",&runNumber,"run/I");
  ttree->Branch("nEvents",&nEvents,"nEvents/F");

  ttree->Branch("avCluEnergySM1",&avCluEnergySM1,"avCluEnergySM1/F");
  ttree->Branch("avCluMultSM1",&avCluMultSM1,"avCluMultSM1/F");
  ttree->Branch("avNcellPerCluSM1",&avNcellPerCluSM1,"avNcellPerCluSM1/F");

  ttree->Branch("avCluEnergySM2",&avCluEnergySM2,"avCluEnergySM2/F");
  ttree->Branch("avCluMultSM2",&avCluMultSM2,"avCluMultSM2/F");
  ttree->Branch("avNcellPerCluSM2",&avNcellPerCluSM2,"avNcellPerCluSM2/F");

  ttree->Branch("avCluEnergySM3",&avCluEnergySM3,"avCluEnergySM3/F");
  ttree->Branch("avCluMultSM3",&avCluMultSM3,"avCluMultSM3/F");
  ttree->Branch("avNcellPerCluSM3",&avNcellPerCluSM3,"avNcellPerCluSM3/F");
  
  char hnam[60]; 
  TH2* h;
  
  Float_t emin = 0.3, emax = 1000.; // minimum and maximum energy of the cluster

  //-------- Number of processed events --------------------------------------------------------------
  TH1* hNEventsProcessedPerRun = (TH1*)listAnyInt->FindObject("hNEventsProcessedPerRun");
  nEvents =  hNEventsProcessedPerRun->GetEntries();

  //-------- Mean cluster energy, number of cells in the cluster Mean number of clusters per event ---
  sprintf(hnam,"run%d_hNCellsInClusterSM1",runNumber); h = (TH2*)listAnyInt->FindObject(hnam); 
  h->GetXaxis()->SetRangeUser(emin,emax);
  avCluEnergySM1 = h->ProjectionX()->GetMean(); avNcellPerCluSM1 = h->ProjectionY()->GetMean();
  avCluMultSM1 = h->Integral()/nEvents;

  sprintf(hnam,"run%d_hNCellsInClusterSM2",runNumber); h = (TH2*)listAnyInt->FindObject(hnam); 
  h->GetXaxis()->SetRangeUser(emin,emax);
  avCluEnergySM2 = h->ProjectionX()->GetMean(); avNcellPerCluSM2 = h->ProjectionY()->GetMean();
  avCluMultSM2 = h->Integral()/nEvents;

  sprintf(hnam,"run%d_hNCellsInClusterSM3",runNumber); h = (TH2*)listAnyInt->FindObject(hnam); 
  h->GetXaxis()->SetRangeUser(emin,emax);
  avCluEnergySM3 = h->ProjectionX()->GetMean(); avNcellPerCluSM3 = h->ProjectionY()->GetMean();
  avCluMultSM3 = h->Integral()/nEvents;;

      
  //---------------------------------------------------------------------------------------------------

  ttree->Fill();
  trendFile->cd();

  ttree->Write();
  trendFile->Close();
  
}

