
void readHit(){
  gROOT->SetStyle("Plain");

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");

  TH2F *xyGlob = new TH2F("xyGlob"," X - Y Global coordinates ",100,-50,50,100,-50,50);
  xyGlob->SetXTitle("cm");
  xyGlob->SetMarkerStyle(7);
  TH1F *zGlob  = new TH1F("zGlob", " Z Global coordinates ",200, -100,100 );
  zGlob->SetXTitle("cm");

  Int_t nbins=100;
  Double_t xmin=0;
  Double_t xmax=1e-04;

  const Int_t nLayers = 8;

  TH1D *hDeLoss[nLayers];
  for(Int_t i=0; i< nLayers; i++ ) {
    hDeLoss[i] = new TH1D(Form("hDeLossl%i",i),Form("E loss distribution [ Layer %i] ",i),nbins,xmin,xmax);
    hDeLoss[i]->SetXTitle("GeV");
  }

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadSDigits();
  runLoader->LoadHits();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  //HIT INIT

  TTree *hitTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITShit");

  AliITSsegmentationUpgrade *segmentation = new AliITSsegmentationUpgrade();
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
 
    runLoader->GetEvent(iEvent);
    
    hitTree=dl->TreeH();
    hitTree->SetBranchAddress("ITSupgrade",&hitList);
    for(Int_t iEnt=0;iEnt<hitTree->GetEntries();iEnt++){//entries loop degli hits
      hitTree->GetEntry(iEnt);
      for(Int_t iHit=0; iHit<hitList->GetEntries();iHit++){
	      
	AliITShit *pHit = (AliITShit*)hitList->At(iHit);
  
	if(pHit->GetParticle()->IsPrimary()){
	  Double_t xg,yg,zg=0.;
	  pHit->GetPositionG(xg,yg,zg);
	  xyGlob->Fill(xg,yg);
	  zGlob->Fill(zg);
	  Int_t module = pHit->GetModule();
	  cout<<module<<" "<<(module%100)<<endl;
	  hDeLoss[pHit->GetModule()%100]->Fill(pHit->GetIonization());
	} // is primary
      }//loop hit 
    }//entryloopHitList
	
  }//event loop

  TCanvas *xyCanv =  new TCanvas("xyCanv","Hit X-Y positions",500,500);
  xyCanv->cd();
  xyGlob->Draw();
  TCanvas *zCanv =  new TCanvas("zCanv","Hit Z positions",500,500);
  zCanv->cd();
  zGlob->Draw();

  TCanvas *c = new TCanvas("c","E loss  distribution per layer",1000,800);
  c->Divide(3,2);
  for(Int_t ip =1; ip<=nLayers; ip++){
    c->cd(ip);
    hDeLoss[ip-1]->Draw();
  }

}
