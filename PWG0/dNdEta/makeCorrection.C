//
// Script to make correction maps for dndeta measurements using the
// dNdEtaCorrection class.
// 

makeCorrection(Char_t* dataDir, Int_t nRuns=20) {

  Char_t str[256];

  gSystem->Load("../libPWG0base.so");

  // ########################################################
  // selection of esd tracks
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();    
  esdTrackCuts->DefineHistograms(1);
  
  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  
  esdTrackCuts->SetMinNsigmaToVertex(3);
  esdTrackCuts->SetAcceptKingDaughters(kFALSE);

  AliLog::SetClassDebugLevel("AliESDtrackCuts",1);

  // ########################################################
  // definition of dNdEta correction object

  AlidNdEtaCorrection* dNdEtaMap = new AlidNdEtaCorrection();

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
  for (Int_t r=0; r<nDirs; r++) {

    TSystemFile* presentDir = (TSystemFile*)dirList->At(r);
    if (!presentDir || !presentDir->IsDirectory())
      continue;
    // check that the files are there
    TString currentDataDir;
    currentDataDir.Form("%s/%s",dataDir, presentDir->GetName());
    cout << "Processing directory " << currentDataDir.Data() << endl;
    if ((!gSystem->Which(currentDataDir,"galice.root")) ||
          (!gSystem->Which(currentDataDir,"AliESDs.root"))) 
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

    sprintf(str,"%s/galice.root",currentDataDir.Data());
    AliRunLoader* runLoader = AliRunLoader::Open(str);

    runLoader->LoadgAlice();
    gAlice = runLoader->GetAliRun();
    runLoader->LoadKinematics();
    runLoader->LoadHeader();

    // #########################################################
    // open esd file and get the tree

    sprintf(str,"%s/AliESDs.root",currentDataDir.Data());
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

    Int_t nEventsTriggered = 0;
    Int_t nEventsAll       = 0;
    
    if (nEvents!=nESDEvents) {
      cout << " Different number of events from runloader and esdtree!!!" << nEvents << " / " << nESDEvents << endl;
      return;
    }
   // ########################################################
    // loop over number of events
    cout << " looping over events..." << endl;
    for(Int_t i=0; i<nEvents; i++) {

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

      Bool_t vertexReconstructed = kTRUE;

      // the vertex should be reconstructed
      if (strcmp(vtxESD->GetName(),"default")==0) 
	vertexReconstructed = kFALSE;

      // the resolution should be reasonable???
      if (vtx_res[2]==0 || vtx_res[2]>0.01)
	vertexReconstructed = kFALSE;

      // ########################################################
      // get the trigger info
      
      Bool_t triggered = kFALSE;

      // MB should be 
      // ITS_SPD_GFO_L0  : 32
      // VZERO_OR_LEFT   : 1
      // VZERO_OR_RIGHT  : 2

      ULong64_t triggerMask = esd->GetTriggerMask();

      nEventsAll++;      

      if (triggerMask&32 && ((triggerMask&1) || (triggerMask&2))) {
	triggered = kTRUE;       
	nEventsTriggered++;
      }

      // ########################################################
      // get the MC vertex
      AliGenPythiaEventHeader* genHeader =
        (AliGenPythiaEventHeader*)runLoader->GetHeader()->GenEventHeader();
      
      TArrayF vtxMC(3);
      genHeader->PrimaryVertex(vtxMC);
      Double_t vtx[3];
      vtx[0] = vtxMC[0];
      vtx[1] = vtxMC[1];
      vtx[2] = vtxMC[2];

      // ########################################################
      // loop over mc particles
      AliStack* particleStack = runLoader->Stack();
      Int_t nPrim    = particleStack->GetNprimary();

      for (Int_t i_mc=0; i_mc<nPrim; i_mc++) {
      
	TParticle* part = particleStack->Particle(i_mc);      
	if (!part || strcmp(part->GetName(),"XXX")==0) 
	  continue;
      
	TParticlePDG* pdgPart = part->GetPDG();

	Bool_t prim = kFALSE;
	// check if it's a primary particle - is there a better way ???
	if ((part->GetFirstDaughter() >= nPrim) || (part->GetFirstDaughter()==-1)) {
	  if (TMath::Abs(pdgPart->PdgCode())>10 && pdgPart->PdgCode()!=21 && strcmp(pdgPart->ParticleClass(),"Unknown")!=0)
	    prim = kTRUE;
	}
	if (!prim)
	  continue;

	if (pdgPart->Charge()==0)
	  continue;	

	Float_t eta = part->Eta();
	Float_t pt  = part->Pt();

	if (prim) {

	  dNdEtaMap->FillParticleAllEvents(eta, pt);		  

	  if (triggered)
	    dNdEtaMap->FillParticleWhenEventTriggered(eta, pt);
	  
	  if (vertexReconstructed)
	    dNdEtaMap->FillParticle(vtx[2], eta, 1.);	
	}
	
      }// end of mc particle
      
      // ########################################################
      // loop over esd tracks      
      Int_t nGoodTracks = 0;

      Int_t nTracks = esd->GetNumberOfTracks();            
      for (Int_t t=0; t<nTracks; t++) {
	AliESDtrack* esdTrack = esd->GetTrack(t);      
	
	// cut the esd track?
	if (!esdTrackCuts->AcceptTrack(esdTrack))
	  continue;
	
	nGoodTracks++;	
	
	Double_t p[3];
	esdTrack->GetPxPyPz(p);
	Float_t momentum = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2) + TMath::Power(p[2],2));
	
	Float_t eta = -100.;
	if((momentum != TMath::Abs(p[2]))&&(momentum != 0))
	  eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
	
	// using the eta of the mc particle
	Int_t label = TMath::Abs(esdTrack->GetLabel());
	if (label<0) {
	  cout << " cannot find corresponding mc part !!! " << label << endl;
	  continue;
	}
	TParticle* mcPart = particleStack->Particle(label);	
	eta = mcPart->Eta();
	Float_t pt = mcPart->Pt();
	
	if (vertexReconstructed)
	  dNdEtaMap->FillParticleWhenMeasuredTrack(vtx[2], eta, pt);	
	
      } // end of track loop
      
      dNdEtaMap->FillEvent(vtxMC[2], nGoodTracks);
      
      if (vertexReconstructed)
	dNdEtaMap->FillEventWithReconstructedVertex(vtxMC[2], nGoodTracks);

    } // end  of event loop
  } // end of run loop
  
  dNdEtaMap->Finish(nEventsAll, nEventsTriggered);  

  TFile* fout = new TFile("correction_map.root","RECREATE");
  
  esdTrackCuts->SaveHistograms("esd_track_cuts");
  dNdEtaMap->SaveHistograms();
  
  fout->Write();
  fout->Close();
  
  dNdEtaMap->DrawHistograms();

}
