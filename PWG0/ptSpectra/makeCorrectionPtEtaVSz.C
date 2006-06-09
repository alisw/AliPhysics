/* $Id$ */

//
// Script to calculate PtvsEta correction map at different z of vtx.
// 

makeCorrectionPtEtaVSz(Char_t* dataDir, Int_t nRuns=20) {

  Char_t str[256];
static const Int_t  NZbin=5;
Float_t Zbin[]={-10., -3.,0., 3., 10.};

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
  AliLog::SetClassDebugLevel("ESDtrackQualityCuts",1);

  // ########################################################
  // definition of PtEta correction objects


CorrectionMatrix2D *PtEtaMap[NZbin];
for(Int_t i=0; i<NZbin-1; i++)
  PtEtaMap[i] = new CorrectionMatrix2D(Form("PtvsEta%d",i),Form("%f<z_{vtx}<%f",Zbin[i],Zbin[i+1]),80,0.,10.,120,-2.,2.);

  PtEtaMap[NZbin-1] = new CorrectionMatrix2D("xxx","xxx",80,0.,10.,120,-2.,2.);


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

      // the vertex should be reconstructed
      if (strcmp(vtxESD->GetName(),"default")==0) 
	continue;

      // the resolution should be reasonable???
      if (vtx_res[2]==0 || vtx_res[2]>0.1)
	continue;

Int_t IZbin=NZbin-1;
	  for(Int_t l=0; l<NZbin-1; l++)
            if((vtx[3]<Zbin[l+1])&&((vtx[3]>Zbin[l])))IZbin=l;



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

	if (prim)
            PtEtaMap[IZbin]->FillGene(part->Pt(), part->Eta());	
	
      }// end of mc particle

      // ########################################################
      // loop over esd tracks      
      Int_t nTracks = esd->GetNumberOfTracks();            
      for (Int_t t=0; t<nTracks; t++) {
	AliESDtrack* esdTrack = esd->GetTrack(t);      
	
	// cut the esd track?
	if (!esdTrackCuts->AcceptTrack(esdTrack))
	  continue;

        AliTPCtrack* tpcTrack = new AliTPCtrack(*esdTrack);
        if (tpcTrack->GetAlpha()==0) {
          cout << " Warning esd track alpha = 0" << endl;
          continue;
        }

        Float_t theta = tpcTrack->Theta();
        Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
	Float_t pt=  = tpcTrack->Pt();


	PtEtaMap[IZbin]->FillMeas(eta, pt);	

      } // end of track loop
    } // end  of event loop
  } // end of run loop

  
  for(Int_t l=0; l<NZbin-1; l++)
  {
  PtEtaMap[l]->Divide();  
  PtEtaMap[l]->RemoveEdges();  
  }
  TFile* fout = new TFile("PtEtaCorrectionMap.root","RECREATE");
  
  esdTrackCuts->SaveHistograms("EsdTrackCuts");
  for(Int_t l=0; l<NZbin-1; l++) PtEtaMap[l]->SaveHistograms();
  
  fout->Write();
  fout->Close();

}
