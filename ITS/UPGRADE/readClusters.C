void readClusters(int nev=-1,int evStart=0)
{

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeRec");
  gROOT->SetStyle("Plain");

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  AliGeomManager::LoadGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  Int_t nLayers = gm->GetNLayers();
  AliITSUClusterPix::SetGeom(gm);
  AliITSURecoDet *its = new AliITSURecoDet(gm, "ITSinterface");
  its->CreateClusterArrays();
  //
  TH2F *xyGlob = new TH2F("xyGlob"," X - Y Global coordinates ",500,-50,50,500,-50,50);
  xyGlob->SetXTitle("cm"); 
  xyGlob->SetMarkerStyle(7); 
  TH1F *zGlob  = new TH1F("zGlob", " Z Global coordinates ",200, -50,50 );
  zGlob->SetXTitle("cm"); 
  //
  TH1F* rGlob = new TH1F("rGlob","R global", 5000, 0,50.);

  TTree * cluTree = 0x0;
  
  int nevTot = (Int_t)runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",nevTot);
  evStart = evStart<nevTot ? evStart : nevTot-1;
  if (evStart<0) evStart = 0;
  //
  int lastEv = nev<0 ? nevTot : evStart+nev;
  if (lastEv > nevTot) lastEv = nevTot;
  //
  for (Int_t iEvent = evStart; iEvent < lastEv; iEvent++) {
    printf("\n Event %i \n",iEvent);
    runLoader->GetEvent(iEvent);
    //   AliStack *stack = runLoader->Stack();
    cluTree=dl->TreeR();
    int nlr=0;
    while(1) {
      TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",nlr));
      if (!br) break;
      br->SetAddress(its->GetLayerActive(nlr)->GetClustersAddress());
      nlr++;
    }    
    printf(" tree entries: %d\n",cluTree->GetEntries());
    cluTree->GetEntry(0);      
    AliITSUClusterPix::SetSortMode( AliITSUClusterPix::SortModeIdTrkYZ());
    for (int ilr=0;ilr<nlr;ilr++) {its->GetLayerActive(ilr)->GetClusters()->Sort();}
    its->ProcessClusters();
    //
    for (int ilr=0;ilr<nlr;ilr++) {
      AliITSURecoLayer* lr = its->GetLayerActive(ilr);
      TClonesArray* clr = lr->GetClusters();
      
      int nClu = clr->GetEntries();
      printf("Layer %d : %d clusters\n",ilr,nClu);
      //
      for (int icl=0;icl<nClu;icl++) {
	AliITSUClusterPix *cl = (AliITSUClusterPix*)clr->At(icl);
	printf("#%4d | ",icl); cl->Print("glo p");
	Float_t loc[3];
	cl->GetLocalXYZ(loc);
	Float_t glob[3]; 
	cl->GetGlobalXYZ(glob);
	//printf("%d: mod %d: loc(%.4lf,%.4lf,%.4lf); glob(%.4lf,%.4lf,%.4lf); \n",icl,cl->GetVolumeId(),
	//       loc[0],loc[1],loc[2],glob[0],glob[1],glob[2]);
	xyGlob->Fill(glob[0],glob[1]);
	zGlob->Fill(glob[2]);
	rGlob->Fill(TMath::Sqrt(glob[0]*glob[0]+glob[1]*glob[1]));
      }
    }
  }//event loop

  Int_t size = 400;

  TCanvas *xyCanv =  new TCanvas("xvCanvClus","RecPoint X-Y Positions",10,10,size,size);
  xyCanv->cd();
  xyGlob->Draw("P");

  TCanvas *zCanv =  new TCanvas("zCanvClus","RecPoint Z Positions",size+20,10,size,size);
  zCanv->cd();
  zGlob->Draw();

  TCanvas *rCanv =  new TCanvas("rCanvClus","RecPoint R Positions",size+20,10,size,size);
  rCanv->cd();
  rGlob->Draw();

}
