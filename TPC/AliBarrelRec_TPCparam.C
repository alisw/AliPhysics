/// \file AliBarrelRec_TPCparam.C
/// \brief Macro to create AliESDs.root using parametrized TPC tracking and AliITStrackerSA (MI + SA)
///
/// Input files:
/// - galice.root
/// - Kinematics.root
/// - TrackReferences.root
/// - ITS.RecPoints.root (use AliRecontruction class)
/// - ITS.Vertex.root (use $ALICE_ROOT/ITS/AliITSVertexerZTest.C)
///
/// \author A. Dainese - INFN Legnaro

void AliBarrelRec_TPCparam(Int_t firstEvent=0,Int_t lastEvent=0) {

  Int_t  collcode = 1; // pp collisions
  Bool_t useMeanVtx = kFALSE;
  
  AliGeomManager::LoadGeometry("geometry.root");
  
  if (gAlice) {
    delete AliRunLoader::Instance();
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
  
  AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetLowFluxParam();
  AliITSReconstructor::SetRecoParam(itsRecoParam);
  
  // Instance of the ITS tracker
  AliITStrackerSA itsTracker(0);
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
  AliVertexerTracks *vertexer = new AliVertexerTracks(AliTracker::GetBz());
  vertexer->SetVtxStart(initVertex);
  vertexer->SetDebug(0);
  delete initVertex;
  initVertex=0;
  
  /***** The TREE for ESD is born *****/
  TTree *esdTree=new TTree("esdTree","Tree with ESD objects");
  AliESDEvent *event=0; AliESDEvent *eventTPCin=0;
  event = new AliESDEvent();
  event->CreateStdContent();
  event->WriteToTree(esdTree);
  
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
    
    cout<<" Processing event number : "<<i<<endl;
    //AliESDEvent *event = new AliESDEvent(); 
    event->SetRunNumber(gAlice->GetRunNumber());
    event->SetEventNumberInFile(i);
    event->SetMagneticField(gAlice->Field()->SolenoidField());
    rl->GetEvent(i);

    //***** Primary vertex from SPD from file 
    sprintf(zver,"Event%d/Vertex",i);
    vertexSPD = (AliESDVertex*)ppZ->Get(zver);
    if(!vertexSPD) {
      esdTree->Fill(); event->Reset();
      continue;
    }      
    event->SetVertex(vertexSPD);
    vertexSPD->GetXYZ(vtx);
    vertexSPD->GetSigmaXYZ(sigmavtx);

    //***** TPC tracking
    if ( (trc=tpcTrackerPar.BuildTPCtracks(event)) ) {
      printf("exiting TPC tracker with code %d in event %d\n",trc,i);
      esdTree->Fill(); event->Reset();
      continue;
    }

    // make a copy of the ESD at this stage
    eventTPCin = event;

    //***** ITS tracking
    itsTracker.AliTracker::SetVertex(vtx,sigmavtx);
    //    itsl->LoadRecPoints("read");
    TTree *itsTree=itsl->TreeR();
    if (!itsTree) {
      cerr<<"Can't get the ITS cluster tree !\n";
      esdTree->Fill(); event->Reset();
      return;
    }     
    itsTracker.UnloadClusters();
    itsTracker.LoadClusters(itsTree);
    if ( (trc=itsTracker.Clusters2Tracks(event)) ) {
      printf("exiting ITS tracker with code %d in event %d\n",trc,i);
      esdTree->Fill(); event->Reset();
      continue;
    }

    // Bring kTPCin-tracks back to the TPC inner wall
    BackToTPCInnerWall(event,eventTPCin);

    // refit inward in ITS:
    // - refit without vertex constraint
    // - propagate through beam pipe to local x = 0
    itsTracker.RefitInward(event);
    

    //***** Vertex from ESD tracks
    if(collcode==1) { // pp
      AliESDVertex *vertexTrks = 
	(AliESDVertex*)vertexer->FindPrimaryVertex(event);
      event->SetPrimaryVertex(vertexTrks); 
    }
    
    esdTree->Fill();
    event->Reset();

  }//<-----------------------------------The Loop over events ends here
  timer.Stop(); timer.Print();
  
  // The AliESDs.root is born
  TFile *ef = TFile::Open("AliESDs.root","RECREATE"); 
  if (!ef || !ef->IsOpen()) {cerr<<"Can't open AliESDs.root !\n"; return;}

  //Write the tree and close everything
  esdTree->Write();
  delete esdTree;
  ef->Close();

  delete vertexer;
  delete rl;

  return;
}
//--------------------------------------------------------------------------
void BackToTPCInnerWall(AliESDEvent *event,AliESDEvent *eventTPC) {

  Int_t ntracks = eventTPC->GetNumberOfTracks();
  AliESDtrack *esdTrackTPC = 0;

  // create relation between event and eventTPC
  Int_t labelsTPC[100000000];
  for(Int_t tr = 0; tr<ntracks; tr++) {
    esdTrackTPC = (AliESDtrack*)event->GetTrack(tr);
    labelsTPC[TMath::Abs(esdTrackTPC->GetLabel())] = tr;
  }

  ntracks = event->GetNumberOfTracks();
  AliESDtrack *esdTrack = 0;
  esdTrackTPC = 0;
  Int_t indexTPC;

  // loop on tracks
  for(tr = 0; tr<ntracks; tr++) {
    esdTrack = (AliESDtrack*)event->GetTrack(tr);
    // set to kITSout the tracks that don't have kTPCin
    // (they've been found by AliITStrackerSA)
    if(!(esdTrack->GetStatus()&AliESDtrack::kTPCin)) {
      esdTrack->SetStatus(AliESDtrack::kITSout);
      continue;
    }

    // skip tracks that don't have kITSin
    if(!(esdTrack->GetStatus()&AliESDtrack::kITSin)) continue;

    indexTPC = labelsTPC[TMath::Abs(esdTrack->GetLabel())];
    esdTrackTPC = (AliESDtrack*)eventTPC->GetTrack(indexTPC);

    AliITStrackMI *itsTrack = 0;
    try {
      itsTrack = new AliITStrackMI(*esdTrackTPC);
      esdTrack = 0;
    }
    catch (const Char_t *msg) {
        Warning("ToTPCInnerWall",msg);
        continue;
    }
    itsTrack->UpdateESDtrack(AliESDtrack::kITSout);
    esdTrack = new AliESDtrack(*(itsTrack->GetESDtrack()));

    delete itsTrack;
  } // end loop on tracks

  return;
}
//--------------------------------------------------------------------------
