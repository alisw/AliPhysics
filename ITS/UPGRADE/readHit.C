
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

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadSDigits();
  runLoader->LoadHits();

  AliGeomManager::LoadGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");



  Int_t nLayers = gm->GetNLayers();

  TH1D **hDeLoss = new TH1D*[nLayers];
  for(Int_t i=0; i< nLayers; i++ ) {
    hDeLoss[i] = new TH1D(Form("hDeLossl%i",i),Form("E loss distribution [ Layer %i] ",i),nbins,xmin,xmax);
    hDeLoss[i]->SetXTitle("GeV");
  }


  //HIT INIT

  TTree *hitTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITSUHit");

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {

    printf("\nEvent\t%d\n",iEvent);
 
    runLoader->GetEvent(iEvent);
    
    hitTree=dl->TreeH();
    hitTree->SetBranchAddress("ITS",&hitList);
    for(Int_t iEnt=0;iEnt<hitTree->GetEntries();iEnt++){//entries loop degli hits
      hitTree->GetEntry(iEnt);
      for(Int_t iHit=0; iHit<hitList->GetEntries();iHit++){
	      
	AliITSUHit *pHit = (AliITSUHit*)hitList->At(iHit);
	int id = pHit->GetModule();
	int lr = gm->GetLayer(id);
	int ld = gm->GetLadder(id);
	//
	//	if(pHit->GetParticle()->IsPrimary()){
	Double_t xg,yg,zg=0.;
	Double_t xg0,yg0,zg0=0.,tg0;
	pHit->GetPositionG(xg,yg,zg);
	pHit->GetPositionG0(xg0,yg0,zg0,tg0);
	xyGlob->Fill(xg,yg);
	zGlob->Fill(zg);
	printf("Module %5d | Lr:%2d Ladder: %3d, X:[%+.5e:%+.5e] Y:[%+.5e:%+.5e] Z:[%+.5e %+.5e] DE: %.5e TrackID: %d\n",id,lr,ld,
	       xg0,xg,yg0,yg,zg0,zg,pHit->GetIonization(),pHit->GetTrack());
	hDeLoss[lr]->Fill(pHit->GetIonization());
	//	} // is primary
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
  c->Divide(nLayers&0x1 ?  (nLayers/2+1) : nLayers/2,2);
  for(Int_t ip =1; ip<=nLayers; ip++){
    c->cd(ip);
    hDeLoss[ip-1]->Draw();
  }

}
