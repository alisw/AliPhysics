/* $Id$ */

testESDtrackCuts(Char_t* dataDir=, Int_t nRuns=10) {

  Char_t str[256];

  gSystem->Load("libESDtrackCuts.so");

  // ########################################################
  // definition of ESD track cuts

  AliESDtrackCuts* trackCuts = new AliESDtrackCuts();
  trackCuts->DefineHistograms(4);

  trackCuts->SetMinNClustersTPC(50);
  trackCuts->SetMaxChi2PerClusterTPC(3.5);
  trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  trackCuts->SetRequireTPCRefit(kTRUE);

  trackCuts->SetMinNsigmaToVertex(3);
  trackCuts->SetAcceptKingDaughters(kFALSE);

  trackCuts->SetPRange(0.3);

  AliLog::SetClassDebugLevel("AliESDtrackCuts",1);

  // ########################################################
  // definition of used pointers
  TFile* esdFile;
  TTree* esdTree;
  TBranch* esdBranch;

  AliESD* esd = 0;

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
  // loop over runs
  Int_t nRunCounter = 0;
  for (Int_t r=1; r<=nDirs; r++) {

    TSystemFile* presentDir = (TSystemFile*)dirList->At(r);
    if (!presentDir->IsDirectory())
      continue;
    // first check that the files are there
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

    // ########################################################
    // getting number of events
    Int_t nEvents    = (Int_t)runLoader->GetNumberOfEvents();
    Int_t nESDEvents = esdBranch->GetEntries();
    
    if (nEvents!=nESDEvents) 
      cout << " Warning: Different number of events from runloader and esdtree!!!" << nEvents << " / " << nESDEvents << endl;
    
    // ########################################################
    // loop over number of events
    cout << " looping over events..." << endl;
    for(Int_t i=1; i<nEvents; i++) {
      
      esdBranch->GetEntry(i);
      runLoader->GetEvent(i);
      
      // ########################################################
      // get the EDS vertex
      AliESDVertex* vtxESD = esd->GetVertex();
      
      Double_t vtxSigma[3];
      vtxESD->GetSigmaXYZ(vtxSigma);
      
      // ########################################################
      // loop over esd tracks      
      Int_t nTracks = esd->GetNumberOfTracks();      

      for (Int_t t=0; t<nTracks; t++) {
	AliESDtrack* esdTrack = esd->GetTrack(t);      
	
	//trackCuts->AcceptTrack(esdTrack, vtxESD, esd->GetMagneticField());
	trackCuts->AcceptTrack(esdTrack);
	
      } // end of track loop
    } // end  of event loop
  } // end of run loop

  TFile* fout = new TFile("out.root","RECREATE");

  trackCuts->SaveHistograms("esd_track_cuts");
  
  fout->Write();
  fout->Close();

}
