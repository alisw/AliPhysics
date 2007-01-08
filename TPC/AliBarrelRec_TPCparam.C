void AliBarrelRec_TPCparam(Int_t firstEvent=0,Int_t lastEvent=0) {
  //
  // Macro to create AliESDs.root using parametrized TPC tracking
  // and AliITStrackerSA (MI + SA)
  //
  // Input files:
  // - galice.root
  // - Kinematics.root
  // - TrackReferences.root
  // - ITS.RecPoints.root (use AliRecontruction class)
  // - ITS.Vertex.root (use $ALICE_ROOT/ITS/AliITSVertexerZTest.C)
  //
  // A. Dainese - INFN Legnaro
  //

  Int_t  collcode = 1; // pp collisions
  Bool_t useMeanVtx = kTRUE;
  
  
  if (gAlice) {
    delete gAlice->GetRunLoader();
    delete gAlice; 
    gAlice=0;
  }
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval) {
    cerr<<"LoadgAlice returned error"<<endl;
    delete rl;
    return;
  }
  retval = rl->LoadHeader();
  if (retval) {
    cerr<<"LoadHeader returned error"<<endl;
    delete rl;
    return;
  }
  gAlice=rl->GetAliRun();
  
  
  
  // Get field from galice.root
  AliMagF *fiel = (AliMagF*)gAlice->Field();
  Double_t fieval=TMath::Abs((Double_t)fiel->SolenoidField()/10.);
  // Set the conversion constant between curvature and Pt
  AliTracker::SetFieldMap(fiel,kTRUE);
  
  /**** The TPC corner ********************/
  
  AliTPCtrackerParam tpcTrackerPar(collcode,fieval);
  tpcTrackerPar.Init();
  
  //**** Switch on the PID class (mandatory!)
  AliPID pid;
  
  /**** The ITS corner ********************/
  
  AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
  if (itsl == 0x0) {
    cerr<<"Cannot get the ITS loader"<<endl;
    return;
  }
  itsl->LoadRecPoints("read");
  
  AliITS *dITS = (AliITS*)gAlice->GetDetector("ITS");
  if (!dITS) {
    cerr<<"Cannot find the ITS detector !"<<endl;
    return;
  }
  AliITSgeom *geom = dITS->GetITSgeom();
  
  // Instance of the ITS tracker
  AliITStrackerSA itsTracker(geom);
  Int_t ITSclusters[6] = {1,1,1,1,1,1};
  itsTracker.SetLayersNotToSkip(ITSclusters);
  
  
  // Primary vertex reconstruction in pp
  AliESDVertex *initVertex = 0;
  TFile *invtx = new TFile("AliESDVertexMean.root");
  if(collcode==1 && useMeanVtx) {
    // open the mean vertex
    initVertex = (AliESDVertex*)invtx->Get("vtxmean");
  } else {
    Double_t pos[3]={0.5,0.5,0.}; 
    Double_t err[3]={3.,3.,5.};
    initVertex = new AliESDVertex(pos,err);
  }
  invtx->Close();
  delete invtx;
  AliVertexerTracks *vertexer = new AliVertexerTracks();
  vertexer->SetVtxStart(initVertex);
  vertexer->SetDebug(0);
  delete initVertex;
  initVertex=0;
  
  /***** The TREE for ESD is born *****/
  TTree *esdTree=new TTree("esdTree","Tree with ESD objects");
  AliESD *event=0;
  esdTree->Branch("ESD","AliESD",&event);
  
  if(firstEvent>rl->GetNumberOfEvents()) firstEvent=rl->GetNumberOfEvents()-1;
  if(lastEvent>rl->GetNumberOfEvents())  lastEvent=rl->GetNumberOfEvents()-1;
  cout<<" Number of events: "<<1+lastEvent-firstEvent<<endl;
  
  TFile *ppZ = TFile::Open("ITS.Vertex.root"); // z vertices from SPD
  AliESDVertex *vertexSPD = new AliESDVertex();
  Char_t zver[100];
  Double_t vtx[3]={0,0,0};
  Double_t sigmavtx[3]={0.07,0.07,0.1};


  //<---------------------------------- The Loop over events begins
  TStopwatch timer;
  Int_t trc;
  for(Int_t i=firstEvent; i<=lastEvent; i++) { 
    
    cerr<<" Processing event number : "<<i<<endl;
    AliESD *event = new AliESD(); 
    event->SetRunNumber(gAlice->GetRunNumber());
    event->SetEventNumber(i);
    event->SetMagneticField(gAlice->Field()->SolenoidField());
    rl->GetEvent(i);

    //***** Primary vertex from SPD from file 
    sprintf(zver,"Event%d/Vertex",i);
    vertexSPD = (AliESDVertex*)ppZ->Get(zver);
    if(!vertexSPD) {
      esdTree->Fill(); delete event;
      continue;
    }      
    event->SetVertex(vertexSPD);
    vertexSPD->GetXYZ(vtx);
    vertexSPD->GetSigmaXYZ(sigmavtx);

    //***** TPC tracking
    if ( (trc=tpcTrackerPar.BuildTPCtracks(event)) ) {
      printf("exiting TPC tracker with code %d in event %d\n",trc,i);
      esdTree->Fill(); delete event;
      continue;
    }

    //***** ITS tracking
    itsTracker.AliTracker::SetVertex(vtx,sigmavtx);
    TTree *itsTree=itsl->TreeR();
    if (!itsTree) {
      cerr<<"Can't get the ITS cluster tree !\n";
      esdTree->Fill(); delete event;
      return;
    }     
    itsTracker.UnloadClusters();
    itsTracker.LoadClusters(itsTree);
    if ( (trc=itsTracker.Clusters2Tracks(event)) ) {
      printf("exiting ITS tracker with code %d in event %d\n",trc,i);
      esdTree->Fill(); delete event;
      continue;
    }
    
    //***** Propagate kITSin tracks to local x = 0 (through beam pipe)
    ToLocalX0(event);

    //***** Vertex from ESD tracks
    if(collcode==1) { // pp
      AliESDVertex *vertexTrks = 
	(AliESDVertex*)vertexer->FindPrimaryVertex(event);
      event->SetPrimaryVertex(vertexTrks); 
    }
    
    esdTree->Fill();
    delete event;

  }//<-----------------------------------The Loop over events ends here
  timer.Stop(); timer.Print();
  
  // The AliESDs.root is born
  TFile *ef = TFile::Open("AliESDs.root","RECREATE"); 
  if (!ef || !ef->IsOpen()) {cerr<<"Can't open AliESDs.root !\n"; return;}

  //Write the TREE and close everything
  esdTree->Write();
  delete esdTree;
  ef->Close();

  delete vertexer;
  delete rl;

  return;
}
//--------------------------------------------------------------------------
void ToLocalX0(AliESD *esd) {

  Int_t ntracks = esd->GetNumberOfTracks();
  AliESDtrack *esdTrack = 0;

  // loop on tracks
  for(Int_t tr = 0; tr<ntracks; tr++) {
    esdTrack = (AliESDtrack*)esd->GetTrack(tr);
    // ask for kITSin
    if(!(esdTrack->GetStatus()&AliESDtrack::kITSin)) continue;
    AliITStrackV2 *itsTrack = 0;
    try {
      itsTrack = new AliITStrackV2(*esdTrack);
      esdTrack = 0;
    }
    catch (const Char_t *msg) {
        Warning("FindPrimaryVertexForCurrentEvent",msg);
        continue;
    }
    // propagate track to beam pipe (0.8 mm of Be) 
    if(itsTrack->GetX()>3.) itsTrack->PropagateTo(3.,0.0023,65.19);
    // propagate track to (0,0)
    itsTrack->PropagateTo(0.,0.,0.);
    itsTrack->UpdateESDtrack(AliESDtrack::kITSin);
    esdTrack = new AliESDtrack(*(itsTrack->GetESDtrack()));
    delete itsTrack;
  } // end loop on tracks

  return;
}
