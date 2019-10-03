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
  Float_t avCluEnergySM4=-9999., avCluMultSM4=-9999., avNcellPerCluSM4=-9999.; // Module 4

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

  ttree->Branch("avCluEnergySM4",&avCluEnergySM4,"avCluEnergySM4/F");
  ttree->Branch("avCluMultSM4",&avCluMultSM4,"avCluMultSM4/F");
  ttree->Branch("avNcellPerCluSM4",&avNcellPerCluSM4,"avNcellPerCluSM4/F");
  
  char hnam[60]; 
  TH2* h;
  
  Float_t emin = 0.3, emax = 1000.; // minimum and maximum energy of the cluster

  //-------- Number of processed events --------------------------------------------------------------
  TH1* hNEventsProcessedPerRun = (TH1*)listAnyInt->FindObject("hNEventsProcessedPerRun");
  if(hNEventsProcessedPerRun) nEvents =  hNEventsProcessedPerRun->GetEntries();

  //-------- Mean cluster energy, number of cells in the cluster Mean number of clusters per event ---
  sprintf(hnam,"run%d_hNCellsInClusterSM1",runNumber); 
  h = (TH2*)listAnyInt->FindObject(hnam);
  
  if(h && nEvents ) {
    h->GetXaxis()->SetRangeUser(emin,emax);
    avCluEnergySM1 = h->ProjectionX()->GetMean(); avNcellPerCluSM1 = h->ProjectionY()->GetMean();
    avCluMultSM1 = h->Integral()/nEvents;
  }

  sprintf(hnam,"run%d_hNCellsInClusterSM2",runNumber); 
  h = (TH2*)listAnyInt->FindObject(hnam);
  
  if(h && nEvents) {
    h->GetXaxis()->SetRangeUser(emin,emax);
    avCluEnergySM2 = h->ProjectionX()->GetMean(); avNcellPerCluSM2 = h->ProjectionY()->GetMean();
    avCluMultSM2 = h->Integral()/nEvents;
  }

  sprintf(hnam,"run%d_hNCellsInClusterSM3",runNumber); 
  h = (TH2*)listAnyInt->FindObject(hnam);

  if(h && nEvents) {
    h->GetXaxis()->SetRangeUser(emin,emax);
    avCluEnergySM3 = h->ProjectionX()->GetMean(); avNcellPerCluSM3 = h->ProjectionY()->GetMean();
    avCluMultSM3 = h->Integral()/nEvents;
  }

  sprintf(hnam,"run%d_hNCellsInClusterSM4",runNumber); 
  h = (TH2*)listAnyInt->FindObject(hnam);

  if(h && nEvents) {
    h->GetXaxis()->SetRangeUser(emin,emax);
    avCluEnergySM4 = h->ProjectionX()->GetMean(); avNcellPerCluSM4 = h->ProjectionY()->GetMean();
    avCluMultSM4 = h->Integral()/nEvents;
  }
    
  //-------- Average pi0s number per event, pi0 mass, width -------------------------------------------
  Double_t nraw, enraw, mass, emass, sigma, esigma;
  Int_t sm; TH1* hm;
  char name[20],leaf[20];


  sm = 1;
  hm = (TH1*)listAnyInt->FindObject(Form("run%i_hPi0MassSM%iSM%i",runNumber,sm,sm));

  Float_t avPi0NumSM1   =0., avPi0MassSM1   =-9999., avPi0SigmaSM1   =-9999.;
  Float_t avPi0NumErrSM1=0., avPi0MassErrSM1=-9999., avPi0SigmaErrSM1=-9999.;
    
  sprintf(name,"avPi0NumSM%d",sm); sprintf(leaf,"avPi0NumSM%d/F",sm);
  ttree->Branch(name,&avPi0NumSM1,leaf);
    
  sprintf(name,"avPi0MassSM%d",sm); sprintf(leaf,"avPi0MassSM%d/F",sm);
  ttree->Branch(name,&avPi0MassSM1,leaf);
  
  sprintf(name,"avPi0SigmaSM%d",sm); sprintf(leaf,"avPi0SigmaSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaSM1,leaf);
  
  sprintf(name,"avPi0NumErrSM%d",sm); sprintf(leaf,"avPi0NumErrSM%d/F",sm);
  ttree->Branch(name,&avPi0NumErrSM1,leaf);
  
  sprintf(name,"avPi0MassErrSM%d",sm); sprintf(leaf,"avPi0MassErrSM%d/F",sm);
  ttree->Branch(name,&avPi0MassErrSM1,leaf);
  
  sprintf(name,"avPi0SigmaErrSM%d",sm); sprintf(leaf,"avPi0SigmaErrSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaErrSM1,leaf);
  
  FitPi0(hm, nraw, enraw, mass, emass, sigma, esigma);
  
  if(nEvents) avPi0NumSM1 = nraw/nEvents; avPi0MassSM1 = mass; avPi0SigmaSM1 = sigma;
  if(nEvents) avPi0NumErrSM1 = enraw/nEvents; avPi0MassErrSM1 = emass; avPi0SigmaErrSM1 = esigma; 
  

  sm = 2;
  hm = (TH1*)listAnyInt->FindObject(Form("run%i_hPi0MassSM%iSM%i",runNumber,sm,sm));
  
  Float_t avPi0NumSM2   =0., avPi0MassSM2   =-9999., avPi0SigmaSM2   =-9999.;
  Float_t avPi0NumErrSM2=0., avPi0MassErrSM2=-9999., avPi0SigmaErrSM2=-9999.;
  
  sprintf(name,"avPi0NumSM%d",sm); sprintf(leaf,"avPi0NumSM%d/F",sm);
  ttree->Branch(name,&avPi0NumSM2,leaf);
  
  sprintf(name,"avPi0MassSM%d",sm); sprintf(leaf,"avPi0MassSM%d/F",sm);
  ttree->Branch(name,&avPi0MassSM2,leaf);
  
  sprintf(name,"avPi0SigmaSM%d",sm); sprintf(leaf,"avPi0SigmaSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaSM2,leaf);
    
  sprintf(name,"avPi0NumErrSM%d",sm); sprintf(leaf,"avPi0NumErrSM%d/F",sm);
  ttree->Branch(name,&avPi0NumErrSM2,leaf);
    
  sprintf(name,"avPi0MassErrSM%d",sm); sprintf(leaf,"avPi0MassErrSM%d/F",sm);
  ttree->Branch(name,&avPi0MassErrSM2,leaf);
  
  sprintf(name,"avPi0SigmaErrSM%d",sm); sprintf(leaf,"avPi0SigmaErrSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaErrSM2,leaf);
  
  FitPi0(hm, nraw, enraw, mass, emass, sigma, esigma);
  
  if(nEvents) avPi0NumSM2 = nraw/nEvents; avPi0MassSM2 = mass; avPi0SigmaSM2 = sigma;
  if(nEvents) avPi0NumErrSM2 = enraw/nEvents; avPi0MassErrSM2 = emass; avPi0SigmaErrSM2 = esigma;
  

  sm = 3;
  hm = (TH1*)listAnyInt->FindObject(Form("run%i_hPi0MassSM%iSM%i",runNumber,sm,sm));
    
  Float_t avPi0NumSM3   =0., avPi0MassSM3   =-9999., avPi0SigmaSM3   =-9999.;
  Float_t avPi0NumErrSM3=0., avPi0MassErrSM3=-9999., avPi0SigmaErrSM3=-9999.;
    
  sprintf(name,"avPi0NumSM%d",sm); sprintf(leaf,"avPi0NumSM%d/F",sm);
  ttree->Branch(name,&avPi0NumSM3,leaf);
  
  sprintf(name,"avPi0MassSM%d",sm); sprintf(leaf,"avPi0MassSM%d/F",sm);
  ttree->Branch(name,&avPi0MassSM3,leaf);
  
  sprintf(name,"avPi0SigmaSM%d",sm); sprintf(leaf,"avPi0SigmaSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaSM3,leaf);
  
  sprintf(name,"avPi0NumErrSM%d",sm); sprintf(leaf,"avPi0NumErrSM%d/F",sm);
  ttree->Branch(name,&avPi0NumErrSM3,leaf);
    
  sprintf(name,"avPi0MassErrSM%d",sm); sprintf(leaf,"avPi0MassErrSM%d/F",sm);
  ttree->Branch(name,&avPi0MassErrSM3,leaf);
    
  sprintf(name,"avPi0SigmaErrSM%d",sm); sprintf(leaf,"avPi0SigmaErrSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaErrSM3,leaf);
  
  FitPi0(hm, nraw, enraw, mass, emass, sigma, esigma);
    
  if(nEvents) avPi0NumSM3 = nraw/nEvents; avPi0MassSM3 = mass; avPi0SigmaSM3 = sigma;
  if(nEvents) avPi0NumErrSM3 = enraw/nEvents; avPi0MassErrSM3 = emass; avPi0SigmaErrSM3 = esigma;
  
  sm = 4;
  hm = (TH1*)listAnyInt->FindObject(Form("run%i_hPi0MassSM%iSM%i",runNumber,sm,sm));
  
  Float_t avPi0NumSM4   =0., avPi0MassSM4   =-9999., avPi0SigmaSM4   =-9999.;
  Float_t avPi0NumErrSM4=0., avPi0MassErrSM4=-9999., avPi0SigmaErrSM4=-9999.;
  
  sprintf(name,"avPi0NumSM%d",sm); sprintf(leaf,"avPi0NumSM%d/F",sm);
  ttree->Branch(name,&avPi0NumSM4,leaf);
  
  sprintf(name,"avPi0MassSM%d",sm); sprintf(leaf,"avPi0MassSM%d/F",sm);
  ttree->Branch(name,&avPi0MassSM4,leaf);
  
  sprintf(name,"avPi0SigmaSM%d",sm); sprintf(leaf,"avPi0SigmaSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaSM4,leaf);
  
  sprintf(name,"avPi0NumErrSM%d",sm); sprintf(leaf,"avPi0NumErrSM%d/F",sm);
  ttree->Branch(name,&avPi0NumErrSM4,leaf);
  
  sprintf(name,"avPi0MassErrSM%d",sm); sprintf(leaf,"avPi0MassErrSM%d/F",sm);
  ttree->Branch(name,&avPi0MassErrSM4,leaf);
  
  sprintf(name,"avPi0SigmaErrSM%d",sm); sprintf(leaf,"avPi0SigmaErrSM%d/F",sm);
  ttree->Branch(name,&avPi0SigmaErrSM4,leaf);
  
  FitPi0(hm, nraw, enraw, mass, emass, sigma, esigma);
  
  if(nEvents) avPi0NumSM4 = nraw/nEvents; avPi0MassSM4 = mass; avPi0SigmaSM4 = sigma;
  if(nEvents) avPi0NumErrSM4 = enraw/nEvents; avPi0MassErrSM4 = emass; avPi0SigmaErrSM4 = esigma;
  
  
  //---------------------------------------------------------------------------------------------------
  
  ttree->Fill();
  trendFile->cd();
  
  ttree->Write();
  trendFile->Close();
  
}


//-----------------------------------------------------------------------------------------------------
void FitPi0(TH1* h, Double_t &nraw, Double_t &enraw,
            Double_t &mass, Double_t &emass,
            Double_t &sigma, Double_t &esigma,
            Double_t emin = 0.05, Double_t emax = 0.3, Int_t rebin = 1)
{
    // Fits the pi0 peak with crystal ball + pol2,
    // fills number of pi0s, mass, width and their errors.
    
    nraw = enraw = 0;
    mass = emass = 0;
    sigma = esigma = 0;

    if(!h) return;
    if (h->GetEntries() == 0) return;
    
    if (rebin > 1) h->Rebin(rebin);
    
    // crystal ball parameters (fixed by hand for EMCAL)
    Double_t alpha = 1.1;  // alpha >= 0
    Double_t n = 2.;       // n > 1
    
    // CB tail parameters
    Double_t a = pow((n/alpha), n) * TMath::Exp(-alpha*alpha/2.);
    Double_t b = n/alpha - alpha;
    
    // integral value under crystal ball with amplitude = 1, sigma = 1
    // (will be sqrt(2pi) at alpha = infinity)
    Double_t nraw11 = a * pow(b+alpha, 1.-n)/(n-1.) + TMath::Sqrt(TMath::Pi()/2.) * TMath::Erfc(-alpha/TMath::Sqrt(2.));
    
    // signal (crystal ball);
    new TF1("cball", Form("(x-[1])/[2] > -%f ?                        \
                          [0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))    \
                          : [0]*%f*(%f-(x-[1])/[2])^(-%f)", alpha, a, b, n));
    
    // background
    new TF1("mypol2", "[0] + [1]*(x-0.135) + [2]*(x-0.135)^2", emin, emax);
    
    // signal + background
    TF1 *fitfun = new TF1("fitfun", "cball + mypol2", emin, emax);
    fitfun->SetParNames("A", "M", "#sigma", "a_{0}", "a_{1}", "a_{2}");
    fitfun->SetLineColor(kRed);
    fitfun->SetLineWidth(2);
    
    // make a preliminary fit to estimate parameters
    TF1* ff = new TF1("fastfit", "gaus(0) + [3]");
    ff->SetParLimits(0, 0., h->GetMaximum()*1.5);
    ff->SetParLimits(1, 0.1, 0.2);
    ff->SetParLimits(2, 0.004,0.030);
    ff->SetParameters(h->GetMaximum()/3., 0.135, 0.010, 0.);
    h->Fit(ff, "0QL", "", 0.105, 0.165);
    
    fitfun->SetParLimits(0, 0., h->GetMaximum()*1.5);
    fitfun->SetParLimits(1, 0.12, 0.15);
    fitfun->SetParLimits(2, 0.004,0.030);
    fitfun->SetParameters(ff->GetParameter(0), ff->GetParameter(1), ff->GetParameter(2), ff->GetParameter(3));
    h->Fit(fitfun,"QLR", "");
    
    // calculate pi0 values
    mass = fitfun->GetParameter(1);
    emass = fitfun->GetParError(1);
    
    sigma = fitfun->GetParameter(2);
    esigma = fitfun->GetParError(2);
    
    Double_t A = fitfun->GetParameter(0);
    Double_t eA = fitfun->GetParError(0);
    
    nraw = nraw11 * A * sigma / h->GetBinWidth(1);
    enraw = nraw * (eA/A + esigma/sigma);
}

