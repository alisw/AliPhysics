//
// Script to test the dN/dEta analysis using the dNdEtaAnalysis and
// dNdEtaCorrection classes. Note that there is a cut on the events,
// so the measurement will be biassed.
//

testAnalysis(Char_t* dataDir, Int_t nRuns=20) {

  Char_t str[256];

  gSystem->Load("../esdTrackCuts/libESDtrackQuality.so");
  gSystem->Load("libdNdEta.so");

  // ########################################################
  // selection of esd tracks
  ESDtrackQualityCuts* esdTrackCuts = new ESDtrackQualityCuts();    
  esdTrackCuts->DefineHistograms(1);
  
  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  
  esdTrackCuts->SetMinNsigmaToVertex(3);
  esdTrackCuts->SetAcceptKingDaughters(kFALSE);

  AliLog::SetClassDebugLevel("ESDtrackQualityCuts",1);

  // ########################################################
  // definition of dNdEta objects

  dNdEtaCorrection* dNdEtaMap = new dNdEtaCorrection();

  dNdEtaMap->LoadHistograms("correction_map.root","dndeta_correction");
  dNdEtaMap->RemoveEdges(2,0,2);

  dNdEtaAnalysis* analyse = new dNdEtaAnalysis("dndeta");

  // ########################################################
  // get the data dir  
  Char_t execDir[256];
  sprintf(execDir,gSystem->pwd());
  TSystemDirectory* baseDir = new TSystemDirectory(".",dataDir);
  TList* dirList            = baseDir->GetListOfFiles();
  Int_t nDirs               = dirList->GetEntries();
  // go back to the dir where this script is executed
  gSystem->cd(execDir);


  // ########################################################
  // definition of used pointers
  TFile* esdFile;
  TTree* esdTree;
  TBranch* esdBranch;

  AliESD* esd =0;

  // ########################################################
  // loop over runs
  Int_t nRunCounter = 0;
  for (Int_t r=1; r<=nDirs; r++) {

    TSystemFile* presentDir = (TSystemFile*)dirList->At(r);
    if (!presentDir->IsDirectory())
      continue;
    // check that the files are there
    sprintf(str,"%s/%s",dataDir, presentDir->GetName());
    if ((!gSystem->Which(str,"galice.root")) ||
        (!gSystem->Which(str,"AliESDs.root"))) 
      continue;
    
    if (nRunCounter++ >= nRuns)
      break;    
    
    cout << "run #" << nRunCounter << endl;
    
    // #########################################################
    // setup galice and runloader
    if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
    }

    sprintf(str,"%s/run%d/galice.root",dataDir,r);
    AliRunLoader* runLoader = AliRunLoader::Open(str);

    runLoader->LoadgAlice();
    gAlice = runLoader->GetAliRun();
    runLoader->LoadKinematics();
    runLoader->LoadHeader();

    // #########################################################
    // open esd file and get the tree

    sprintf(str,"%s/run%d/AliESDs.root",dataDir,r);
    // close it first to avoid memory leak
    if (esdFile)
      if (esdFile->IsOpen())
        esdFile->Close();

    esdFile = TFile::Open(str);
    esdTree = (TTree*)esdFile->Get("esdTree");
    if (!esdTree)
      continue;
    esdBranch = esdTree->GetBranch("ESD");
    esdBranch->SetAddress(&esd);
    if (!esdBranch)
      continue;

    // ########################################################
    // Magnetic field
    AliTracker::SetFieldMap(gAlice->Field(),kTRUE); // kTRUE means uniform magnetic field
    //AliKalmanTrack::SetConvConst(1000/0.299792458/5.);

    // ########################################################
    // getting number of events
    Int_t nEvents = (Int_t)runLoader->GetNumberOfEvents();
    Int_t nESDEvents = esdBranch->GetEntries();
    
    if (nEvents!=nESDEvents) {
      cout << " Different number of events from runloader and esdtree!!!" << nEvents << " / " << nESDEvents << endl;
      return;
    }
    // ########################################################
    // loop over number of events
    cout << " looping over events..." << endl;
    for(Int_t i=1; i<nEvents; i++) {

      esdBranch->GetEntry(i);
      runLoader->GetEvent(i);            
  

      // ########################################################
      // get the EDS vertex
      AliESDVertex* vtxESD = esd->GetVertex();

      Double_t vtx[3];
      Double_t vtx_res[3];
      vtxESD->GetXYZ(vtx);          
    
      vtx_res[0] = vtxESD->GetXRes();
      vtx_res[1] = vtxESD->GetYRes();
      vtx_res[2] = vtxESD->GetZRes();

      // we need a good vertex 
      //  => there will be a bias on the measurement, since this cuts away some events
      if (strcmp(vtxESD->GetName(),"default")==0) 
	continue;
      if (vtx_res[2]==0 || vtx_res[2]>0.1)
	continue;

      // ########################################################
      // loop over esd tracks      
      Int_t nTracks = esd->GetNumberOfTracks();            
      for (Int_t t=0; t<nTracks; t++) {
	AliESDtrack* esdTrack = esd->GetTrack(t);      
	
	if (!esdTrackCuts->AcceptTrack(esdTrack))
	  continue;
	
	AliTPCtrack* tpcTrack = new AliTPCtrack(*esdTrack);
	
	if (tpcTrack->GetAlpha()==0) {
	  cout << " Warning esd track alpha = 0" << endl;
	  continue;
	}

	Float_t theta = tpcTrack->Theta();
	Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
	Float_t correction = dNdEtaMap->GetCorrection(vtx[2], eta);
	
	dNdEtaMap->FillMeas(vtx[2], eta);

	analyse   ->FillTrack(vtx[2], eta, correction);
	
      } // end of track loop
      analyse->FillEvent(vtx[2]);

    } // end  of event loop
  } // end of run loop

  analyse->Finish();

  TFile* fout = new TFile("out.root","RECREATE");
  
  esdTrackCuts->SaveHistograms("esd_tracks_cuts");
  dNdEtaMap->SaveHistograms();
  analyse   ->SaveHistograms();
  
  fout->Write();
  fout->Close();

}
