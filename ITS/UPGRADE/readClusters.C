void readClusters(){

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

  TH2F *xyGlob = new TH2F("xyGlob"," X - Y Global coordinates ",500,-50,50,500,-50,50);
  xyGlob->SetXTitle("cm"); 
  xyGlob->SetMarkerStyle(7); 
  TH1F *zGlob  = new TH1F("zGlob", " Z Global coordinates ",200, -50,50 );
  zGlob->SetXTitle("cm"); 


  TTree * cluTree = 0x0;
  TObjArray layerClus;

  printf("N Events : %i \n",(Int_t)runLoader->GetNumberOfEvents());

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    printf("\n Event %i \n",iEvent);
    runLoader->GetEvent(iEvent);
    //   AliStack *stack = runLoader->Stack();
    cluTree=dl->TreeR();
    int nlr=0;
    while(1) {
      TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",nlr));
      if (!br) break;
      TClonesArray* clr = 0;
      br->SetAddress(&clr);
      layerClus.AddLast(clr);
      nlr++;
    }

    printf(" tree entries: %d\n",cluTree->GetEntries());
    cluTree->GetEntry(0);      
    //
    for (int ilr=0;ilr<nlr;ilr++) {
      TClonesArray* clr = (TClonesArray*)layerClus.At(ilr);
      int nClu = clr->GetEntries();
      printf("Layer %d : %d clusters\n",ilr,nClu);
      //
      for (int icl=0;icl<nClu;icl++) {
	AliITSUClusterPix *cl = (AliITSUClusterPix*)clr->At(icl);
	cl->Print("glo");
	Double_t loc[3]={cl->GetX(),cl->GetY(),cl->GetZ()}; 
	Double_t glob[3]; 
	gm->LocalToGlobal(cl->GetVolumeId(),loc,glob);
	//	printf("%d: mod %d: loc(%.4lf,%.4lf,%.4lf); glob(%.4lf,%.4lf,%.4lf); \n",icl,cl->GetVolumeId(),
	//	       loc[0],loc[1],loc[2],glob[0],glob[1],glob[2]);
	xyGlob->Fill(glob[0],glob[1]);
	zGlob->Fill(glob[2]);
 
      }
    }
    layerClus.Clear();
  }//event loop

  Int_t size = 400;

  TCanvas *xyCanv =  new TCanvas("xvCanvClus","RecPoint X-Y Positions",10,10,size,size);
  xyCanv->cd();
  xyGlob->Draw("P");

  TCanvas *zCanv =  new TCanvas("zCanvClus","RecPoint Z Positions",size+20,10,size,size);
  zCanv->cd();
  zGlob->Draw();

}
