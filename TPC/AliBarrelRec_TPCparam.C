void AliBarrelRec_TPCparam(Int_t firstEvent=0,Int_t lastEvent=0) {
  //
  // Macro to create AliESDs.root using parametrized TPC tracking
  // and AliITStrackerV2
  //
  // Input files:
  // - galice.root
  // - Kinematics.root
  // - TrackReferences.root
  // - ITS.RecPoints.root (use AliRecontruction class)
  // - ITS.Vertex.root (use $ALICE_ROOT/ITS/AliITSVertexerZTest.C)
  //
  // A. Dainese - LNL
  //


   /**** Initialization of the NewIO *******/

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
       

   TDatabasePDG *DataBase = TDatabasePDG::Instance();

   // Get field from galice.root
   AliMagF *fiel = (AliMagF*)gAlice->Field();
   Double_t fieval=TMath::Abs((Double_t)fiel->SolenoidField()/10.);
   // Set the conversion constant between curvature and Pt
   AliTracker::SetFieldMap(fiel,kTRUE);

   /**** The TPC corner ********************/

   Int_t collcode = 1; // pp collisions
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

   //An instance of the ITS tracker
   AliITStrackerV2 itsTracker(geom);
   Int_t ITSclusters[6] = {1,1,1,1,1,1};
   itsTracker.SetLayersNotToSkip(ITSclusters);

   /***** The TREE is born *****/
   
   TTree *esdTree=new TTree("esdTree","Tree with ESD objects");
   AliESD *event=0;
   esdTree->Branch("ESD","AliESD",&event);
   
   if(firstEvent>rl->GetNumberOfEvents()) firstEvent=rl->GetNumberOfEvents()-1;
   if(lastEvent>rl->GetNumberOfEvents())  lastEvent=rl->GetNumberOfEvents()-1;
   cout<<" Number of events: "<<1+lastEvent-firstEvent<<endl;
   
   TFile *ppZ = TFile::Open("ITS.Vertex.root"); // z vertices from SPD
   AliESDVertex *myvertex = new AliESDVertex();
   Char_t zver[100];
   Double_t vtx[3];
   Double_t cvtx[6];


   //<----------------------------------The Loop over events begins
   TStopwatch timer;
   Int_t trc;
   for(Int_t i=firstEvent; i<=lastEvent; i++) { 
     
     cerr<<" Processing event number : "<<i<<endl;
     AliESD *event = new AliESD(); 
     event->SetRunNumber(gAlice->GetRunNumber());
     event->SetEventNumber(i);
     event->SetMagneticField(gAlice->Field()->SolenoidField());
     rl->GetEvent(i);

    //***** Primary vertex 
     sprintf(zver,"Event%d/Vertex",i);
     myvertex = (AliESDVertex*)ppZ->Get(zver);
     if(!myvertex) {
       esdTree->Fill(); delete event;
       continue;
     }      
     event->SetVertex(myvertex);
     myvertex->GetXYZ(vtx);
     myvertex->GetCovMatrix(cvtx);

     //***** TPC tracking
     if ( (trc=tpcTrackerPar.BuildTPCtracks(event)) ) {
       printf("exiting TPC tracker with code %d in event %d\n",trc,i);
       esdTree->Fill(); delete event;
       continue;
     }

     //***** ITS tracking
     itsTracker.SetVertex(vtx,cvtx);
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


     esdTree->Fill();
     delete event;

   }//<-----------------------------------The Loop over events ends here
   timer.Stop(); timer.Print();

   //        The AliESDs.root is born
   TFile *ef = TFile::Open("AliESDs.root","RECREATE"); 
   if (!ef || !ef->IsOpen()) {cerr<<"Can't open AliESDs.root !\n"; return;}

   esdTree->Write(); //Write the TREE and close everything
   delete esdTree;
   ef->Close();

   delete rl;

   return;
}
