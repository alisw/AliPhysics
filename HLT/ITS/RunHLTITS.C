// The following macro runs the HLT ITS tracker over the HLT
// tracks stored in the ESD and stores the output tracks in a
// separate ESD file AliESDits.root

Int_t RunHLTITS(Int_t nev=1,Int_t run=0) {

  TStopwatch timer;
  timer.Start();

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }

   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"AliESDtest.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 1;
   }
   retval = rl->LoadHeader();
   if (retval) {
      cerr<<"AliESDtest.C : LoadHeader returned error"<<endl;
      delete rl;
      return 2;
   }
   gAlice=rl->GetAliRun();
       

   AliKalmanTrack::SetConvConst(
      1000/0.299792458/gAlice->Field()->SolenoidField()
   );

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliESDtest.C : Can not get the ITS loader"<<endl;
      return 3;
   }
   itsl->LoadRecPoints("read");

   AliITS *dITS = (AliITS*)gAlice->GetDetector("ITS");
   if (!dITS) {
      cerr<<"AliESDtest.C : Can not find the ITS detector !"<<endl;
      return 4;
   }
   AliITSgeom *geom = dITS->GetITSgeom();

   //An instance of the HLT ITS tracker
   AliL3ITStracker itsTracker(geom);

   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef || !ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   AliESD* event = new AliESD;
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   tree->SetBranchAddress("ESD", &event);

   TFile *itsf=TFile::Open("AliESDits.root","RECREATE");
   if ((!itsf)||(!itsf->IsOpen())) {
      cerr<<"Can't AliESDits.root !\n"; return 1;
   }

   Int_t rc=0;
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   //The loop over events
   for (Int_t i=0; i<nev; i++) {

     cerr<<"\n\nProcessing event number : "<<i<<endl;
     tree->GetEvent(i);
     
     rl->GetEvent(i);

     TArrayF v(3);     
     rl->GetHeader()->GenEventHeader()->PrimaryVertex(v);
     Double_t vtx[3]={v[0],v[1],v[2]};
     Double_t cvtx[3]={0.005,0.005,0.010};
     AliESDVertex vertex1(vtx,cvtx);
     event->SetVertex(&vertex1);

     Double_t vtx[3];
     Double_t cvtx[3];
     AliESDVertex *vertex = event->GetVertex();
     vertex->GetXYZ(vtx);
     vertex->GetSigmaXYZ(cvtx);
     itsTracker.SetVertex(vtx,cvtx);

     TTree *itsTree=itsl->TreeR();
     if (!itsTree) {
       cerr<<"Can't get the ITS cluster tree !\n";
       return 4;
     }     
     itsTracker.LoadClusters(itsTree);
     rc+=itsTracker.Clusters2Tracks(event);
     //     rc+=itsTracker.PropagateBack(event);
     itsTracker.UnloadClusters();

     if (rc==0) {
       TTree* tree = new TTree("esdTree", "Tree with ESD objects");
       tree->Branch("ESD", "AliESD", &event);
       tree->Fill();
       itsf->cd();
       tree->Write();
     } 
     if (rc) {
       cerr<<"Something bad happened...\n";
     }

   }
   delete event;

   itsf->Close();
   ef->Close();

   //   delete rl;

   timer.Stop();
   timer.Print();

   return rc;
}
