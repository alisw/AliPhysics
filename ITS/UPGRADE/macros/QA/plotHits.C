void plotHits(Int_t ilay=3, Int_t ilrs=1, Int_t nev=-1, Int_t evStart=0)
{
  // Adapted from readHit.C macro - M.S. 10 jul 14

  enum {kLSLayer, kLSStave , kLSHalfStave, kLSModule,
	kLSChip , kLSSensor, kLSNumPar   };
  char lrsNames[kLSNumPar][10];
  Int_t i = -1;
  strcpy(lrsNames[++i],"Layer");
  strcpy(lrsNames[++i],"Stave");
  strcpy(lrsNames[++i],"HalfStave");
  strcpy(lrsNames[++i],"Module");
  strcpy(lrsNames[++i],"Chip");
  strcpy(lrsNames[++i],"Sensor");

  Int_t idWrapVol[7] = {0, 0, 0, 1, 1, 2, 2};

  Bool_t argerr = kFALSE;

  if (ilay < 0 || ilay > 6 ||
      ilrs < 0 || ilrs >= kLSNumPar ) argerr = kTRUE;

  if (argerr) {
    printf("Wrong parameters! ilay = %d - ilrs = %d\n\n",ilay,ilrs);
    printf("Usage:\n.x plotHit(ilay, ilrs, nev, evStart)\n");
    printf("   ilay : layer number (0 < ilay < 7, default 0)\n");
    printf("   ilrs : local reference system (0 layer, 1 stave, 2 half stave, 3 module, 4 chip, 5 sensor, default 1)\n");
    printf("   nev : number of events (default -1 = all)\n");
    printf("   evStart : first event (default 0)\n");
    return;
  }

  printf("--> Plotting hits at %s level for Layer %d <--\n",
	 lrsNames[ilrs],ilay);

  gROOT->SetStyle("Plain");

  // Load libraries, geometry, pointers, whatever
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");

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

  // Create histos
  char volName[15];
  sprintf(volName,"ITSU%s%d",lrsNames[ilrs],ilay);
  TGeoBBox *vShape = (TGeoBBox*)(gGeoManager->GetVolume(volName)->GetShape());
  Double_t xHist = 1.1*vShape->GetDX(); // 10% greater
  Double_t yHist = 1.1*vShape->GetDY();
  Double_t zHist = 1.1*vShape->GetDZ();

  // OB stave are strongly asymmetric in Y: dimension increased for better view
  if (ilay > 2 && ilrs == kLSStave) yHist*=1.20;

  char histTitle[50];

  sprintf(histTitle," X - Y %s Local coordinates",lrsNames[ilrs]);
  TH2F *xy = new TH2F("xy",histTitle,100,-xHist,xHist,100,-yHist,yHist);
  xy->SetXTitle("X (cm)");
  xy->SetYTitle("Y (cm)");
  xy->SetMarkerStyle(7);

  sprintf(histTitle," X - Z %s Local coordinates",lrsNames[ilrs]);
  TH2F *xz = new TH2F("xz",histTitle,100,-xHist,xHist,100,-zHist,zHist);
  xz->SetXTitle("X (cm)");
  xz->SetYTitle("Z (cm)");
  xz->SetMarkerStyle(7);

  sprintf(histTitle," Z - Y %s Local coordinates",lrsNames[ilrs]);
  TH2F *zy = new TH2F("zy",histTitle,100,-zHist,zHist,100,-yHist,yHist);
  zy->SetXTitle("Z (cm)");
  zy->SetYTitle("Y (cm)");
  zy->SetMarkerStyle(7);

  // Get number of events which to loop on
  Int_t nevTot = (Int_t)runLoader->GetNumberOfEvents();
  printf("Total events : %i\n\n",nevTot);
  evStart = evStart<nevTot ? evStart : nevTot-1;
  if (evStart<0) evStart = 0;

  Int_t lastEv = nev<0 ? nevTot : evStart+nev;
  if (lastEv > nevTot) lastEv = nevTot;

  // Loop on hit structure
  TTree *hitTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITSUHit");

  TString path;
  TGeoMatrix *mat = 0x0;
  char nodName[15], parName[15];

  for (Int_t iEvent = evStart; iEvent < lastEv; iEvent++) {

    printf("Event\t%d\n",iEvent);
 
    runLoader->GetEvent(iEvent);
    
    hitTree = dl->TreeH();
    hitTree->SetBranchAddress("ITS",&hitList);

    for(Int_t iEnt=0; iEnt<hitTree->GetEntries(); iEnt++){
      hitTree->GetEntry(iEnt);

      for(Int_t iHit=0; iHit<hitList->GetEntries();iHit++){
	AliITSUHit *pHit = (AliITSUHit*)hitList->At(iHit);

	Int_t id = pHit->GetChip();
	Int_t lr = gm->GetLayer(id);
	Int_t st = gm->GetStave(id);
	Int_t hs = gm->GetHalfStave(id);
	Int_t md = gm->GetModule(id);
	Int_t cp = gm->GetChipIdInModule(id);

	if (lr == ilay) {
	  Double_t xg, yg, zg=0.;
	  Double_t loc[3], mas[3];

	  pHit->GetPositionG(xg, yg, zg);
	  mas[0] = xg;
	  mas[1] = yg;
	  mas[2] = zg;
	  path.Clear();
	  path = Form("/ALIC_1/ITSV_2/");
	  path += Form("ITSUWrapVol%d_1/ITSULayer%d_1",idWrapVol[ilay],ilay);
	  if (ilrs > kLSLayer)
	    path += Form("/ITSUStave%d_%d",ilay,st);
	  if (ilrs > kLSStave)
	    path += Form("/ITSUHalfStave%d_%d",ilay,hs);
	  if (ilrs > kLSHalfStave)
	    path += Form("/ITSUModule%d_%d",ilay,md);
	  if (ilrs > kLSModule)
	    path += Form("/ITSUChip%d_%d",ilay,cp);
	  if (ilrs > kLSChip)
	    path += Form("/ITSUSensor%d_1",ilay);

	  // Get the matrix for the current volume
	  gGeoManager->PushPath();
	  if (!gGeoManager->cd(path.Data())) {
	    gGeoManager->PopPath();
	    printf("Error in cd-ing to %s!\n",path.Data());
	    printf("Chip id %5d | Lr:%2d Stave:%3d, HStave:%3d Module:%3d Chip:%3d\n",
		   id,lr,st,hs,md,cp);
	  
	    return 1;
	  }
	  mat = gGeoManager->GetCurrentMatrix();
	  gGeoManager->PopPath();

	  // Convert global coordinate to local reference system
	  mat->MasterToLocal(mas, loc);
	  xy->Fill(loc[0],loc[1]);
	  xz->Fill(loc[0],loc[2]);
	  zy->Fill(loc[2],loc[1]);

	} // lr == ilay
      }//loop hit 
    }//entryloopHitList
	
  }//event loop

  // Finally plot the histos
  TCanvas *xyCanv =  new TCanvas("xyCanv","Hit X-Y positions",500,500);
  xyCanv->cd();
  xy->Draw();

  TCanvas *xzCanv =  new TCanvas("xzCanv","Hit X-Z positions",500,500);
  xzCanv->cd();
  xz->Draw();

  TCanvas *zyCanv =  new TCanvas("zyCanv","Hit Z-Y positions",500,500);
  zyCanv->cd();
  zy->Draw();

}
