void AliITSspdLayer1Coverage(TString hfn="galice.root",Int_t mod=-1,
                     Int_t evnt=-1){
  // Macro to show the coverage of each spd layer 1 module

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    if(gAlice){
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
  }
  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");
  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number 0\n");
    man->SetDefaultStorage("local://$ALICE_ROOT");
    man->SetRun(0);
  }else {
    printf("Using deafult storage \n");
  }
  // retrives geometry 
  TString geof(gSystem->DirName(hfn));
  geof += "/geometry.root";
  TGeoManager::Import(geof.Data());
  if (!gGeoManager) {
    cout<<"geometry not found\n";
    return -1;
  }


  AliRunLoader *rl = AccessFile(hfn); // Set up to read in Data
  Int_t retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSPrintHits.C : LoadHeader returned error"<<endl;
    return;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");

  if(!ITSloader){
    cerr<<"AliITSPrintHits.C :  ITS loader not found"<<endl;
    return;
  }

  ITSloader->LoadHits("read");
  AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
  if(!ITS || ITS==0){
    cout << "Error: no ITS found. Aborting"<<endl;
    return;
  } // end if !ITS
  //cout << ITS << endl;
  AliITSInitGeometry *initgeom = new AliITSInitGeometry(ITS->GetMajorVersion(),ITS->GetMinorVersion());
  //cout << initgeom << endl;
  AliITSgeom *geom = initgeom->CreateAliITSgeom();
  if(!geom){
    cout << "Error: not AliITSgeom object found."<<endl;
    return;
  } // end if geom
  ITSloader->SetITSgeom(geom);

  Int_t evNumber1 = 0;
  Int_t evNumber2 = AliRunLoader::GetNumberOfEvents();
  if(evnt>=0){
    evNumber1 = evnt;
    evNumber2 = evnt+1;
  } // end if evnt>=0
  Int_t mod1 = 0;
  Int_t mod2 = geom->GetIndexMax();
  if(mod>=0){
    mod1 = mod;
    mod2 = mod+1;
  } // end if mod>=0
  AliITShit *hp = 0;

  Int_t i;
  Double_t zbin[161];
  for(i=0;i<161;i++) zbin[i] = 425.0E-4;
  zbin[32] = zbin[64] = zbin[96] = zbin[128] = 625.0E-4;
  zbin[0] = -3.536;
  for(i=1;i<161;i++) zbin[i] += zbin[i-1];
  Char_t name[20],title[40];
  TH2I *detCoverage[4*2*10];
  for(i=0;i<4*2*10;i++){
    sprintf(name,"detCoverage[%d]",i);
    sprintf(title,"Cell coverage for SPD module %d",i);
    detCoverage[i] = new TH2I(name,title,160,(Double_t *)zbin,256,-0.64,+0.64);
  } // end for i
  Double_t xg,yg,zg,xl,yl,zl,phi;
  Int_t nmodules,size=-1;
  Int_t event,m,i,i2,hit,trk,lay,lad,det;
  for(event = evNumber1; event < evNumber2; event++){
    //cout<<"Processing event "<<event<<endl;
    rl->GetEvent(event);
    ITS->InitModules(size,nmodules);
    ITS->FillModules(event,0,-1," "," ");
    for(m=mod1;m<mod2;m++)if(geom->GetModuleType(m)==(AliITSDetector)(0)){
      i2 = (ITS->GetModule(m))->GetNhits();
      //cout <<  "Event=" << event << " module=" << m <<
      //  " Number of Hits=" << i2 <<endl;
      for(i=0;i<i2;i++){
        trk = (ITS->GetModule(m))->GetHitTrackIndex(i);
        hit = (ITS->GetModule(m))->GetHitHitIndex(i);
        hp  = (ITS->GetModule(m))->GetHit(i);
	hp->GetPositionG(xg,yg,zg);
	hp->GetPositionL(xl,yl,zl);
	geom->GetModuleId(m,lay,lad,det);
	phi = TMath::ATan2(yg,xg)*TMath::RadToDeg();
	if(phi<0.0) phi+=360.0;
	switch(lay){
	case 1:
	  detCoverage[m]->Fill(zl,xl,1.0);
	  break;
	case 2:
	break;
	default:
	} // end switch
      } // end for i
    } // end for m
    ITS->ClearModules();
  } // end for event
  //
  cout << "Creating pad"<<endl;
  TVirtualPad *pad=0;
  TCanvas *c0 = new TCanvas("c0","SPD Hits Tests",262,50,700,534);
  c0->Divide(4,20,.0,0.0);
  //
  for(i=0;i<4*2*10;i++){
    c0->cd(i+1);
    detCoverage[i]->Draw();
  } // end for i
  c0->Update();
}
