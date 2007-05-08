void AliITSPrintRecPoints(Int_t outtype=1,TString rfn="galice.root",
                          Int_t mod=-1,Int_t evnt=-1){
  // Macro to print out the recpoints for all or a specific module

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } 
  else {
    if(gAlice){
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
  }

  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

  AliRunLoader *rl = AccessFile(rfn); // Set up to read in Data
  Int_t retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSPrintRecPoints.C : LoadHeader returned error"<<endl;
    return;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");

  if(!ITSloader){
    cerr<<"AliITSPrintRecPoints.C :  ITS loader not found"<<endl;
    return;
  }
  cout <<"ITSloader ok"<<endl;

  if(!gGeoManager){
      gGeoManger = new TGeoManager();
      gGeoManager->Import("geometry.root","");
  } // end if
  ITSloader->LoadHits("read");
  ITSloader->LoadDigits("read");
  ITSloader->LoadRecPoints("read");
  cout << "loaded hits, digits, and RecPoints"<< endl;
  AliITS *ITS = 0;
  ITS = (AliITS*)(gAlice->GetDetector("ITS"));
  if(!ITS){
    cout << "Error: no ITS found. Aborting"<<endl;
    return;
  } // end if !ITS
  cout <<"ITS="<<ITS<<endl;
  if(!(ITS->GetDetTypeSim())){
      cout <<"No AliITSDetTypeSim object found in ITS"<<endl;
      return;
  } // end if
  AliITSgeom *gm=0;
  gm = ITS->GetITSgeom();
  if(!gm){
      cout <<"No AliITSgeom object found in ITS"<<endl;
      if(!gGeoManager){
          cout <<"No gGeoManger. Aborting"<<endl;
          return;
      }else{
          ITS->UpdateInternalGeometry(gGeoManager);
      } // end if
  } // end if !AliITSgeom

  Int_t evNumber1 = 0;
  Int_t evNumber2 = gAlice->GetEventsPerRun();
  if(evnt>=0){
    evNumber1 = evnt;
    evNumber2 = evnt+1;
  } // end if evnt>=0
  Int_t mod1 = 0;
  Int_t mod2 = gm->GetIndexMax();
  if(mod>=0){
    mod1 = mod;
    mod2 = mod+1;
  } // end if mod>=0
  TClonesArray *rpa;
  AliITSRecPoint *rp = 0;
  AliITSDetTypeRec* rec = new AliITSDetTypeRec(ITSloader);
  //rec->SetITSgeom(gm);
  rec->SetDefaults();

  Int_t event,m,i,i2;
  Float_t xyz[3];
  for(event = evNumber1; event < evNumber2; event++){
    rl->GetEvent(event);
    rec->SetTreeAddress();
    for(m=mod1;m<mod2;m++){
      rec->ResetRecPoints();
      TTree *TR = ITSloader->TreeR();
      TR->GetEvent(m);
      rpa = rec->RecPoints();
      i2 = rpa->GetEntriesFast();
      cout <<  "Event=" << event << " module=" << m <<
	" Number of Recpoints=" << i2 <<endl;
      for(i=0;i<i2;i++){
          rp = (AliITSRecPoint*)(rpa->At(i));
          switch(outtype){
          case 1:
              rp->GetGlobalXYZ(xyz);
              cout << i << " lx=" << rp->GetDetLocalX() 
                        << " lz=" << rp->GetDetLocalZ() 
                   << " x=" << rp->GetX() 
                   << " y=" << rp->GetY()<< " z=" << rp->GetZ() 
                   <<" gx=" << xyz[0] << " gy="<< xyz[1] <<" gz="<<xyz[2]
                   << endl;
              break;
          default:
              cout << i << " ";
              rp->Print((ostream*)cout);
              cout << endl;
              break;
          } // end switch
      } // end for i
    } // end for m
  } // end for event

}
