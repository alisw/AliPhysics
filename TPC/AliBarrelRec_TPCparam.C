/****************************************************************************
 * This macro performs track and vertex reconstruction in TPC and ITS.      *
 * The ITS Kalman tracker V2 is feeded "with" parameterized TPC tracks.     * 
 *                                                                          *
 * Reconstruction is performed in the following steps:                      *
 *             1) TPC tracking parameterization                             *
 *             2) ITS clusters: slow or fast                                *
 *             3) Primary vertex reconstruction                             *
 *                - read from event header for Pb-Pb events                 *
 *                - determined using points in pixels for pp/pA events      *
 *             4) ITS track finding V2                                      *
 *                - in pp/pA, redetermine the position of primary vertex    *
 *                  using the reconstructed tracks                          *
 *             5) Create a reference file with simulation info (p,PDG...)   *
 *                                                                          *
 * If mode='A' all 5 steps are executed                                     *
 * If mode='B' only steps 4-5 are executed                                  *
 *                                                                          *  
 *  Origin: A.Dainese, Padova,   andrea.dainese@pd.infn.it                  * 
 *  (from AliTPCtest.C & AliITStestV2.C by I.Belikov)                       *
 ****************************************************************************/

// structure for track references
typedef struct {
  Int_t lab;
  Int_t pdg;
  Int_t mumlab;
  Int_t mumpdg;
  Float_t Vx,Vy,Vz;
  Float_t Px,Py,Pz;
} RECTRACK;

//===== Functions definition ================================================= 

void CopyVtx(const Char_t *inName,const Char_t *outName);

void ITSFindClustersV2(Char_t SlowOrFast);

void ITSFindTracksV2(Int_t *skipEvt);

void ITSMakeRefFile(Int_t *skipEvt);

void MarkEvtsToSkip(const Char_t *evtsName,Int_t *skipEvt);

void PrimaryVertex(const Char_t *outName,Char_t vtxMode);

void TPCParamTracks(Int_t coll,Double_t Bfield);

Int_t UpdateEvtsToSkip(const Char_t *logName,const Char_t *evtsName);

void VtxFromHeader(const Char_t *outName,Bool_t smear);

void VtxFromTracks(const Char_t *outName);

void ZvtxFromSPD(const Char_t *outName);

//=============================================================================

// number of events to be processed
Int_t    gNevents;
// magnetic field
Double_t gBfieldValue;

void AliBarrelRec_TPCparam(Int_t n=-1,Char_t mode='A') {

  //---------------------------------------------------------------------
  //                    CONFIGURATION
  //
  // _Magnetic_field_
  gBfieldValue = 0.4;
  //
  // _Type_of_collision_ (needed for TPC tracking parameterization) 
  // Available choices:   !!! ONLY B = 0.4 TESLA !!!
  //    collcode = 0  ->   PbPb6000 (HIJING with b<2fm) 
  //    collcode = 1  ->   low multiplicity: pp or pA
  Int_t collcode = 1;  
  // 
  // _ITS_clusters_reconstruction_
  // Available choices:  (from AliITStestV2.C)
  //    SlowOrFast = 's'    slow points
  //    SlowOrFast = 'f'    fast points
  Char_t SlowOrFast = 'f';
  //
  // _Primary_vertex_for_ITS_tracking_
  // Available choices:
  //    Vtx4Tracking = 'H'   from event Header
  //    --- for Pb-Pb ---
  //    Vtx4Tracking = 'S'   from event header + Smearing 
  //                                           (x=15,y=15,z=10) micron
  //    --- for pp/pA ---
  //    Vtx4Tracking = 'P'   z from pixels, x,y in(0,0)
  Char_t Vtx4Tracking = 'P';
  // _Primary_vertex_for_analysis_ (AliESDVertex stored in tracks file)
  // Available choices:
  //    Vtx4Analysis = 'C'   Copy the same used for tracking
  //    --- for pp/pA ---
  //    Vtx4Analysis = 'T'   x,y,z from Tracks
  Char_t Vtx4Analysis = 'T';
  //
  //                  END CONFIGURATION
  //---------------------------------------------------------------------

  const Char_t *name=" AliBarrelRec_TPCparam";
  printf("\n %s\n",name);
  gBenchmark->Start(name);

  if(n==-1) { // read number of events to be processed from file
    TFile *f = new TFile("galice.root");
    gAlice = (AliRun*)f->Get("gAlice");
    n = gAlice->GetEventsPerRun();
    delete gAlice; 
    gAlice=0;
    f->Close(); 
    delete f;
    printf(" All %d events in file will be processed\n",n);
  }
  gNevents = n;


  // conversion constant for kalman tracks 
  AliKalmanTrack::SetConvConst(100/0.299792458/gBfieldValue);

  // Selection of execution mode
  switch(mode) {
  case 'A':
    // Build TPC tracks with parameterization
    TPCParamTracks(collcode,gBfieldValue);
  
    // ITS clusters
    ITSFindClustersV2(SlowOrFast);

    // Vertex for ITS tracking
    PrimaryVertex("Vtx4Tracking.root",Vtx4Tracking);

    break;
  
  case 'B':
    printf("       ---> only tracking in ITS <---\n");

    // Update list of events to be skipped
    if(!UpdateEvtsToSkip("itstracking.log","evtsToSkip.dat")) return;

    break;
  }

  // Mark events that have to be skipped (if any)
  Int_t *skipEvt = new Int_t[gNevents];
  for(Int_t i=0; i<gNevents; i++) skipEvt[i] = 0;
  if(!gSystem->AccessPathName("evtsToSkip.dat",kFileExists)) 
    MarkEvtsToSkip("evtsToSkip.dat",skipEvt);
    
  // Tracking in ITS
  ITSFindTracksV2(skipEvt); 

  // Vertex for analysis 
  PrimaryVertex("AliITStracksV2.root",Vtx4Analysis);

  // Create ITS tracks reference file
  ITSMakeRefFile(skipEvt);
  delete [] skipEvt;



  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//-----------------------------------------------------------------------------
void CopyVtx(const Char_t *inName,const Char_t *outName) {

  // Open input and output files
  TFile *inFile = new TFile(inName);
  TFile *outFile = new TFile(outName,"update");

  TDirectory *curdir;
  Char_t vname[20];


  for(Int_t ev=0; ev<gNevents; ev++) {
    sprintf(vname,"Vertex_%d",ev);
    AliESDVertex *vertex = (AliESDVertex*)inFile->Get(vname);
    if(!vertex) continue;
    curdir = gDirectory;
    outFile->cd();
    vertex->Write();
    curdir->cd();
    vertex = 0;
  }

  inFile->Close();
  outFile->Close();
  delete inFile;
  delete outFile;

  return;
}
//-----------------------------------------------------------------------------
void ITSFindClustersV2(Char_t SlowOrFast) {

  printf("\n------------------------------------\n");

  const Char_t *name="ITSFindClustersV2";
  printf("\n %s\n",name);
  gBenchmark->Start(name);

  //---  taken from AliITStestV2.C--------------------------------------
  //
  if (SlowOrFast=='f') {
    //cerr<<"Fast AliITSRecPoint(s) !\n";
    //gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2FastRecPoints.C");
    //AliITSHits2FastRecPoints();
  } else {
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2SDigits.C");
    AliITSHits2SDigits();
    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSSDigits2Digits.C");
    AliITSSDigits2Digits();
    //gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSDigits2RecPoints.C");
    //AliITSDigits2RecPoints();
  }
  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersV2.C");
  AliITSFindClustersV2(SlowOrFast,gNevents);
  //
  //--------------------------------------------------------------------


  gBenchmark->Stop(name);
  gBenchmark->Show(name);

   return;
}
//-----------------------------------------------------------------------------
Int_t ITSFindTracksV2(Int_t *skipEvt) {

  printf("\n------------------------------------\n");
  
  const Char_t *name="ITSFindTracksV2";
  printf("\n %s\n",name);
  gBenchmark->Start(name);

 
  TFile *outFile     = new TFile("AliITStracksV2.root","recreate");
  TFile *inTPCtrks   = new TFile("AliTPCtracksParam.root");
  TFile *inVertex    = new TFile("Vtx4Tracking.root");
  TFile *inClusters  = new TFile("AliITSclustersV2.root");
  
  AliITSgeom *geom=(AliITSgeom*)inClusters->Get("AliITSgeom");
  if(!geom) { printf("can't get ITS geometry !\n"); return;}
  
  Double_t vtx[3];
  Int_t flag1stPass,flag2ndPass;
  Char_t vname[20];

  // open logfile for done events
  FILE *logfile = fopen("itstracking.log","w");

  // Instantiate AliITStrackerV2
  AliITStrackerV2 tracker(geom);

  // loop on events
  for(Int_t ev=0; ev<gNevents; ev++){
    // write to logfile of begun events
    fprintf(logfile,"%d\n",ev);

    if(skipEvt[ev]) continue;
    printf(" --- Processing event %d ---\n",ev);

    // pass event number to the tracker
    tracker.SetEventNumber(ev);

    // set position of primary vertex
    sprintf(vname,"Vertex_%d",ev);
    AliESDVertex *vertex = (AliESDVertex*)inVertex->Get(vname);
    if(vertex) {
      vertex->GetXYZ(vtx);
      delete vertex;
    } else {
      printf(" AliESDVertex not found for event %d\n",ev);
      printf(" Using (0,0,0) for ITS tracking\n");
      vtx[0] = vtx[1] = vtx[2] = 0.;
    }

    flag1stPass=1; // vtx constraint
    flag2ndPass=0; // no vtx constraint

    // no vtx constraint if vertex not found
    if(vtx[2]<-999.) {
      flag1stPass=0;
      vtx[2]=0.;
    }

    tracker.SetVertex(vtx);  

    // setup vertex constraint in the two tracking passes
    Int_t flags[2];
    flags[0]=flag1stPass;
    tracker.SetupFirstPass(flags);
    flags[0]=flag2ndPass;
    tracker.SetupSecondPass(flags);
    
    // find the tracks
    tracker.Clusters2Tracks(inTPCtrks,outFile);

  } // loop on events

  fprintf(logfile,"%d\n",gNevents); //this means all evts are successfully completed
  fclose(logfile);

  delete geom;

  inTPCtrks->Close();
  inClusters->Close();
  inVertex->Close();
  outFile->Close();
 

  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  
  return;
}
//-----------------------------------------------------------------------------
void ITSMakeRefFile(Int_t *skipEvt) {

  printf("\n------------------------------------\n");

  const Char_t *name="ITSMakeRefFile";
  printf("\n %s\n",name);
  gBenchmark->Start(name);
  
  
  TFile *out = TFile::Open("ITStracksRefFile.root","recreate");
  TFile *trk = TFile::Open("AliITStracksV2.root");
  TFile *kin = TFile::Open("galice.root");

  
  // Get gAlice object from file
  gAlice=(AliRun*)kin->Get("gAlice");
  
  Int_t label;
  TParticle *Part;  
  TParticle *Mum;
  RECTRACK rectrk;
  

  for(Int_t ev=0; ev<gNevents; ev++){
    if(skipEvt[ev]) continue;
    printf(" --- Processing event %d ---\n",ev);

    gAlice->GetEvent(ev);  

    trk->cd();

    // Tree with ITS tracks
    char tname[100];
    sprintf(tname,"TreeT_ITS_%d",ev);

    TTree *tracktree=(TTree*)trk->Get(tname);
    if(!tracktree) continue;
    AliITStrackV2 *itstrack=new AliITStrackV2; 
    tracktree->SetBranchAddress("tracks",&itstrack);
    Int_t nentr=(Int_t)tracktree->GetEntries();

    // Tree for true track parameters
    char ttname[100];
    sprintf(ttname,"Tree_Ref_%d",ev);
    TTree *reftree = new TTree(ttname,"Tree with true track params");
    reftree->Branch("rectracks",&rectrk,"lab/I:pdg:mumlab:mumpdg:Vx/F:Vy:Vz:Px:Py:Pz");

    for(Int_t i=0; i<nentr; i++) {
      tracktree->GetEvent(i);
      label = TMath::Abs(itstrack->GetLabel());

      Part = (TParticle*)gAlice->Particle(label);
      rectrk.lab=label;
      rectrk.pdg=Part->GetPdgCode();
      rectrk.mumlab = Part->GetFirstMother();
      if(Part->GetFirstMother()>=0) {
	Mum = (TParticle*)gAlice->Particle(Part->GetFirstMother());
	rectrk.mumpdg=Mum->GetPdgCode();
      } else {
	rectrk.mumpdg=-1;
      }
      rectrk.Vx=Part->Vx();
      rectrk.Vy=Part->Vy();
      rectrk.Vz=Part->Vz();
      rectrk.Px=Part->Px();
      rectrk.Py=Part->Py();
      rectrk.Pz=Part->Pz();
      
      reftree->Fill();
    } // loop on tracks   

    out->cd();
    reftree->Write();

    delete itstrack;
    delete reftree;
  } // loop on events

  trk->Close();
  kin->Close();
  out->Close();
  
  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  

  return;
}
//-----------------------------------------------------------------------------
void MarkEvtsToSkip(const Char_t *evtsName,Int_t *skipEvt) {

  printf("\n------------------------------------\n");
  printf("\nChecking for events to skip...\n");

  Int_t evt,ncol;

  FILE *f = fopen(evtsName,"r");
  while(1) {
    ncol = fscanf(f,"%d",&evt);
    if(ncol<1) break;
    skipEvt[evt] = 1;
    printf(" event %d will be skipped\n",evt);
  }
  fclose(f);

  return;
}
//-----------------------------------------------------------------------------
void PrimaryVertex(const Char_t *outName,Char_t vtxMode) {

  printf("\n------------------------------------\n");

  const Char_t *name="PrimaryVertex";
  printf("\n %s\n",name);
  gBenchmark->Start(name);

  switch(vtxMode) {
  case 'H':
    printf(" ... from event header\n");
    VtxFromHeader(outName,kFALSE);
    break;
  case 'S':
    printf(" ... from event header + smearing\n");
    VtxFromHeader(outName,kTRUE);
    break;
  case 'P':
    printf(" ... z from pixels for pp/pA\n");
    ZvtxFromSPD(outName);
    break;
  case 'T':
    printf(" ... from tracks for pp/pA\n");
    VtxFromTracks(outName);
    break;
  case 'C':
    printf(" ... copied from Vtx4Tracking.root to AliITStracksV2.root\n");
    CopyVtx("Vtx4Tracking.root",outName);
    break;
  }

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//-----------------------------------------------------------------------------
void TPCParamTracks(Int_t coll,Double_t Bfield) {

  printf("\n------------------------------------\n");

  const Char_t *name="TPCParamTracks";
  printf("\n %s\n",name);
  gBenchmark->Start(name);

  TFile *outFile=TFile::Open("AliTPCtracksParam.root","recreate");
  TFile *inFile =TFile::Open("galice.root");
 
  AliTPCtrackerParam tracker(coll,Bfield,gNevents);
  tracker.BuildTPCtracks(inFile,outFile);

  delete gAlice; gAlice=0;

  inFile->Close();
  outFile->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//-----------------------------------------------------------------------------
Int_t UpdateEvtsToSkip(const Char_t *logName,const Char_t *evtsName) {

    if(!gSystem->AccessPathName(logName,kFileExists)) { 
      FILE *ifile = fopen(logName,"r");
      Int_t lEvt=0,nCol=1;
      while(nCol>0) {
	nCol = fscanf(ifile,"%d",&lEvt);
      }
      fclose(ifile);
      if(lEvt==gNevents) { 
	printf(" All events already reconstructed\n"); 
	return 0;  
      } else {
	FILE *ofile = fopen("evtsToSkip.dat","a");
	fprintf(ofile,"%d\n",lEvt);
	fclose(ofile);
      }
    } else { 
      printf("File itstracking.log not found\n");
    }

    return 1;
}
//-----------------------------------------------------------------------------
void VtxFromHeader(const Char_t *outName,Bool_t smear) {

  TDatime t;
  UInt_t seed = t.Get();
  gRandom->SetSeed(seed);

  TFile *galice  = new TFile("galice.root");  
  TFile *outFile = new TFile(outName,"update");

  TDirectory *curdir;
  Double_t pos[3],sigma[3];
  if(smear) {
    sigma[0]=15.e-4;
    sigma[1]=15.e-4;
    sigma[2]=10.e-4;
  } else {
    sigma[0]=0.;
    sigma[1]=0.;
    sigma[2]=0.;
  }
  Char_t vname[20];

  galice->cd();

  for(Int_t ev=0; ev<gNevents; ev++){
    printf(" event %d\n",ev);
    sprintf(vname,"Vertex_%d",ev);
    TArrayF o = 0;
    o.Set(3);
    AliHeader* header = 0;
    TTree* treeE = (TTree*)gDirectory->Get("TE");
    treeE->SetBranchAddress("Header",&header);
    treeE->GetEntry(ev);
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if(genHeader) {
      // get primary vertex position
      genHeader->PrimaryVertex(o);
      pos[0] = (Double_t)o[0];
      pos[1] = (Double_t)o[1];
      pos[2] = (Double_t)o[2];
      if(smear) {
	pos[0] = gRandom->Gaus(pos[0],sigma[0]);
	pos[1] = gRandom->Gaus(pos[1],sigma[1]);
	pos[2] = gRandom->Gaus(pos[2],sigma[2]);
      }
      // create AliESDVertex
      AliESDVertex *vertex = new AliESDVertex(pos,sigma,vname);
    } else {
      printf(" ! event header not found : setting vertex to (0,0,0) !");
      pos[0] = 0.;
      pos[1] = 0.;
      pos[2] = 0.;
      // create AliESDVertex
      AliESDVertex *vertex = new AliESDVertex(pos,sigma,vname);
    }    
    delete header;
    // write AliESDVertex to file
    curdir = gDirectory;
    outFile->cd();
    if(smear) {
      vertex->SetTitle("vertex from header, smeared");
    } else {
      vertex->SetTitle("vertex from header");
    }
    vertex->Write();
    curdir->cd();
    vertex = 0;
  }

  outFile->Close();
  galice->Close();

  delete outFile;
  delete galice;

  return;
}
//-----------------------------------------------------------------------------
void VtxFromTracks(const Char_t *outName) {

  // Open input and output files
  TFile *inFile  = new TFile("AliITStracksV2.root");
  TFile *outFile = new TFile(outName,"update");

  // set AliRun object to 0
  if(gAlice) gAlice = 0;

  // Create vertexer
  AliITSVertexerTracks *vertexer = 
    new AliITSVertexerTracks(inFile,outFile,gBfieldValue);
  vertexer->SetFirstEvent(0);
  vertexer->SetLastEvent(gNevents-1);
  vertexer->SetDebug(0);
  vertexer->PrintStatus();
  // Find vertices
  vertexer->FindVertices();

  delete vertexer;

  inFile->Close();
  outFile->Close();
  delete inFile;
  delete outFile;

  return;
}
//-----------------------------------------------------------------------------
void ZvtxFromSPD(const Char_t *outName) {

  // create fast RecPoints, which are used for vertex finding
  cerr<<"Fast AliITSRecPoint(s) !\n";
  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2FastRecPoints.C");
  AliITSHits2FastRecPoints(0,gNevents-1);

  // delphi ---> azimuthal range to accept tracklets
  // window ---> window in Z around the peak of tracklets proj. in mm
  Float_t delphi=0.05;
  Float_t window=3.;
  Float_t initx=0.;
  Float_t inity=0.;

  TFile *infile = new TFile("galice.root");
  TFile *outfile = new TFile(outName,"update");

  AliITSVertexerPPZ *vertexer = new AliITSVertexerPPZ(infile,outfile,initx,inity);
  vertexer->SetFirstEvent(0);
  vertexer->SetLastEvent(gNevents-1);
  vertexer->SetDebug(0);
  vertexer->SetDiffPhiMax(delphi);
  vertexer->SetWindow(window);
  vertexer->PrintStatus();
  vertexer->FindVertices();
  delete vertexer;
  vertexer=0;

  outfile->Close();
  infile->Close();
  delete infile;
  delete outfile;


  return;
}
//-----------------------------------------------------------------------------









