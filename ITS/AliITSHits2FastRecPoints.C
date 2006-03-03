void AliITSHits2FastRecPoints (Int_t evNumber1=0,Int_t evNumber2=0, TString inFile = "galice.root", Int_t nsignal=25, Int_t size=-1) 
{
  /////////////////////////////////////////////////////////////////////////
  //   
  //   This macro creates fast recpoints, optionally on a separate file
  //   
  /////////////////////////////////////////////////////////////////////////


  // Dynamically link some shared libs

  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else if (gAlice){
     delete gAlice->GetRunLoader();
     delete gAlice; 
     gAlice=0;
  }

  // Connect the Root Galice file containing Geometry, Kine and Hits
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    if (rl == 0x0)
     {
       ::Error("AliITSHits2FastRecPoints.C","Can not open session RL=NULL");
       return;
     }
     
    Int_t retval = rl->LoadgAlice();
    if (retval)
     {
       ::Error("AliITSHits2FastRecPoints.C","LoadgAlice returned error");
       delete rl;
       return;
     }
    gAlice=rl->GetAliRun();
    rl->LoadHeader();
    retval = rl->LoadKinematics();
    if (retval)
     {
       ::Error("AliITSHits2FastRecPoints.C","LoadKinematics returned error");
       delete rl;
       return;
     }
    
    AliITSLoader* gime = (AliITSLoader*)rl->GetLoader("ITSLoader");
    if (gime == 0x0)
     {
       ::Error("AliITSHits2FastRecPoints.C","can not get ITS loader");
       delete rl;
       return;
     }
    retval = gime->LoadHits("read");
    if (retval)
     {
       ::Error("AliITSHits2FastRecPoints.C","LoadHits returned error");
       delete rl;
       return;
     }
    gime->SetRecPointsFileName("ITS.FastRecPoints.root"); 
    retval = gime->LoadRecPoints("update");
    if (retval)
     {
       ::Error("AliITSHits2FastRecPoints.C","LoadRecPoints returned error");
       delete rl;
       return;
     }
    
   
    TDirectory * olddir = gDirectory;
    rl->CdGAFile();
    AliITSgeom* geom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
    olddir->cd();
    if(!geom){
      Error("GetITSgeom","no ITS geometry available");
      return NULL;
    }

 
   AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
   if (!ITS) return;

  // Set the simulation model

   

  //
  // Event Loop
  //

  Int_t nbgr_ev=0;
  TStopwatch timer;

  cout << "Creating fast reconstructed points from hits for the ITS..." << endl;
  AliITSDetTypeSim* dettyp = new AliITSDetTypeSim();
  dettyp->SetITSgeom(geom);
  dettyp->SetLoader(gime);
  ITS->SetDetTypeSim(dettyp);
  for (Int_t i=0;i<3;i++) {
    ITS->SetSimulationModel(i,new AliITSsimulationFastPoints());
  }

  
  for (int ev=evNumber1; ev<= evNumber2; ev++) {
    cout << "...working on event "<< ev << " ..." << endl;
    Int_t nparticles = gAlice->GetEvent(ev);
    cout << "event         " <<ev<<endl;
    cout << "nparticles  " <<nparticles<<endl;
    rl->GetEvent(ev);
    //if(gime->TreeR() == 0x0) gime->MakeTree("R");
    
     if (ev < evNumber1) continue;
    if (nparticles <= 0) return;

    Int_t bgr_ev=Int_t(ev/nsignal);
    timer.Start();
    ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
    timer.Stop(); timer.Print();
  } // event loop 

  delete rl;
}

