void AliITSAnalizeSPDHits(TString hfn="galice.root",Int_t mod=-1,
                     Int_t evnt=-1){
  // Macro to analize hits in the SPD to chatch posible problems

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
  TH1I *phiDistribution1 = new TH1I("phiDistribution1",
		   "Phi Distribution for the inner most layer of the SPD",
				    5027,0.0,360.0);
  TH2I *detCoverage1 = new TH2I("detCoverage1",
	      "The SPD layer 1 cell coverage summed over the first ladder per sector",
				160,(Double_t *)zbin,256,-0.64,+0.64);
  TH2I *detCoverage1a = new TH2I("detCoverage1a",
	      "The SPD layer 1 cell coverage sum over the second ladder per sector",
				 160,(Double_t *)zbin,256,-0.64,+0.64);
  TH2I *detThicknessZ1 = new TH2I("detThicknessZ1",
                    "Hit local y distribution as a function of z",
				  160,(Double_t*)zbin,200,-100.E-4,+100.E-4);
  TH2I *detThicknessX1 = new TH2I("detThicknessX1",
                    "Hit local y distribution as a function of x",
				  256,-0.64,+0.64,200,-100.E-4,+100.E-4);
  TH1I *phiDistribution2 = new TH1I("phiDistribution2",
		   "Phi Distribution for the outer most layer of the SPD",
				    8800,0.0,360.0);
  TH2I *detCoverage2 = new TH2I("detCoverage2",
	      "The SPD layer 2 cell coverage summed over all but the forth ladder per sector",
				160,(Double_t *)zbin,256,-0.64,+0.64);
  TH2I *detCoverage2a = new TH2I("detCoverage2a",
	      "The SPD layer 2 cell coverage summed over all forth ladder per sector",
				 160,(Double_t *)zbin,256,-0.64,+0.64);
  TH2I *detThicknessZ2 = new TH2I("detThicknessZ2",
                    "Hit local y distribution as a function of z",
				  160,(Double_t*)zbin,200,-100.E-4,+100.E-4);
  TH2I *detThicknessX2 = new TH2I("detThicknessX2",
                    "Hit local y distribution as a function of x",
				  256,-0.64,+0.64,200,-100.E-4,+100.E-4);
  TH2I *xySpace = new TH2I("xySpace",
			   "ALICE global coorinate location of hits in X,y",
			   100,-8.0,+8.0,100,-8.0,+8.0);
  TH2I *zrSpace = new TH2I("zrSpace",
			   "ALICE global coorinate location of hits in z,r",
			   100,-15.0,+15.0,100,0.0,+8.0);
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
	  phiDistribution1->Fill(phi,1.0);
	  if(lad%2==1) detCoverage1->Fill(zl,xl,1.0);
	  if(lad%2==0) detCoverage1a->Fill(zl,xl,1.0);
	  detThicknessZ1->Fill(zl,yl,1.0);
	  detThicknessX1->Fill(xl,yl,1.0);
	  break;
	case 2:
	  phiDistribution2->Fill(phi,1.0);
	  if(lad%4!=0) detCoverage2->Fill(zl,xl,1.0);
	  if(lad%4==0) detCoverage2a->Fill(zl,xl,1.0);
	  detThicknessZ2->Fill(zl,yl,1.0);
	  detThicknessX2->Fill(xl,yl,1.0);
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
  c0->Range(0,0,1,1);
  c0->Divide(3,2);
  //
  pad = c0->cd(1);
  phiDistribution2->Draw();
  pad = c0->cd(4);
  phiDistribution1->Draw();
  TPad *pad2 = c0->cd(2);
  delete pad2;
  c0->cd(0);
  TPad *pad21 = new TPad("pad21","CoverageZ",0.34,0.875,0.65,1.0);
  pad21->Draw();
  pad21->cd();
  pad21->Range(0,0,1,1);
  pad21->Modified();
  detCoverage2->Draw();
  c0->Update();
  c0->cd(0);
  TPad *pad22 = new TPad("pad22","CoverageZa",0.34,0.75,0.65,0.875);
  pad22->Draw();
  pad22->cd();
  pad22->Range(0,0,1,1);
  pad22->Modified();
  detCoverage2a->Draw();
  c0->Update();
  c0->cd(0);
  TPad *pad23 = new TPad("pad23","ThicknessZ",0.34,0.625,0.65,0.75);
  pad23->Draw();
  pad23->cd();
  pad23->Range(0,0,1,1);
  pad23->Modified();
  detThicknessZ2->Draw();
  c0->Update();
  c0->cd(0);
  TPad *pad24 = new TPad("pad24","ThicknessX",0.34,0.50,0.65,0.625);
  pad24->Draw();
  pad24->cd();
  pad24->Range(0,0,1,1);
  pad24->Modified();
  detThicknessX2->Draw();
  c0->Update();
  //
  c0->cd(0);
  TPad *pad5 = c0->cd(5);
  delete pad5;
  c0->cd(0);
  TPad *pad51 = new TPad("pad51","CoverageZ",0.34,0.375,0.65,0.50);
  pad51->Draw();
  pad51->cd();
  pad51->Range(0,0,1,1);
  pad51->Modified();
  detCoverage1->Draw();
  c0->Update();
  c0->cd(0);
  TPad *pad52 = new TPad("pad52","CoverageZa",0.34,0.25,0.65,0.375);
  pad52->Draw();
  pad52->cd();
  pad52->Range(0,0,1,1);
  pad52->Modified();
  detCoverage1a->Draw();
  c0->Update();
  c0->cd(0);
  TPad *pad53 = new TPad("pad53","ThicknessZ",0.34,0.125,0.65,0.25);
  pad53->Draw();
  pad53->cd();
  pad53->Range(0,0,1,1);
  pad53->Modified();
  detThicknessZ1->Draw();
  c0->Update();
  c0->cd(0);
  TPad *pad54 = new TPad("pad54","ThicknessX",0.34,0.00,0.65,0.125);
  pad54->Draw();
  pad54->cd();
  pad54->Range(0,0,1,1);
  pad54->Modified();
  detThicknessX1->Draw();
  c0->Update();
  //
  TList *tvTreeList = new TList;
  TTree *tvTree = (TTree *) gROOT->FindObject("TreeH");
  tvTreeList->Add(tvTree);
  pad = c0->cd(3);
  tvTree->Draw("ITS.fY:ITS.fX","ITS.fModule<240","", 100000, 0);
  pad = c0->cd(6);
  tvTree->Draw("ITS.fX*ITS.fX+ITS.fY*ITS.fY:ITS.fZ","ITS.fModule<240","", 
		 100000, 0);

}
