//========================================================//
//Macro that reads the output of the hydro calculations 
//(Therminator) using the format that can be found at
//arXiv:1102.0273
//Author: Panos.Christakoglou@nikhef.nl
//========================================================//

//========================================================//
//Event structure
struct StructEvent {
  UInt_t eventID;//unique identifier of the event (random number)
  UInt_t entries;//number of particles of each event
  UInt_t entriesprev;//set to 0 by default
};
//========================================================//

//========================================================//
//Particle structure
//Primordial particles: fatherid==-1
struct StructParticle {
  Float_t mass;//the mass of the particle (in GeV)
  Float_t t;//the time coordinate of the particle in fm/c
  Float_t x;//the spacial coordinate x in fm/c
  Float_t y;//the spacial coordinate y in fm/c
  Float_t z;//the spacial coordinate z in fm/c
  Float_t e;//the energy in GeV
  Float_t px;//the x-coordinate of the particle's momentum in GeV
  Float_t py;//the y-coordinate of the particle's momentum in GeV
  Float_t pz;//the z-coordinate of the particle's momentum in GeV
  Int_t decayed;//decay flag (1: decayed, 0: no decay)
  Int_t pid;//PDG of particle
  Int_t fatherpid;//PDG of parent
  Int_t rootpid;//root (primordial) particle PDG number
  Int_t eid;//sequence number in the event
  Int_t fathereid;//parent sequence number in the event
  UInt_t eventid;//unique identifier of the event (random number)
};
//========================================================//

//========================================================//
//Balance function analysis variables
Bool_t gRunShuffling=kFALSE;
Bool_t gRunMixing=kFALSE;
Bool_t gRunMixingWithEventPlane=kFALSE;

Double_t gEtaMin = -0.8;
Double_t gEtaMax = 0.8;
Double_t gPtMin = 0.2;
Double_t gPtMax = 20.0;
//========================================================//


void runBalanceFunctionOnHydro(TString aEventDir = "/glusterfs/alice1/alice2/pchrist/Therminator/pA/ChargeBalancing/Centrality1/lhyquid3v-LHCpPb5020s3.1Ti319t0.20Tf150e001/", Int_t aEventFiles = 2) {
  //Macro that reads the themrinator events
  //Author: Panos.Christakoglou@nikhef.nl
  // Time:
  TStopwatch timer;
  timer.Start();

  //========================================================//
  //Load the aliroot libraries
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFebye.so");
  //========================================================//

  //========================================================//
  //Configure the bf objects
  AliBalancePsi *bf = new AliBalancePsi();
  bf->SetAnalysisLevel("MC");
  bf->SetShuffle(gRunShuffling);
  bf->SetEventClass("EventPlane");
  bf->SetDeltaEtaMax(TMath::Abs(gEtaMax-gEtaMin));
  bf->InitHistograms();

  AliBalancePsi *bfs = 0x0;
  if(gRunShuffling) {
    bfs = new AliBalancePsi();
    bfs->SetAnalysisLevel("MC");
    bfs->SetEventClass("EventPlane");
    bfs->SetDeltaEtaMax(TMath::Abs(gEtaMax-gEtaMin));
    bfs->InitHistograms();
  }

  AliBalancePsi *bfm = 0x0;
  if(gRunMixing) {
    bfm = new AliBalancePsi();
    bfm->SetAnalysisLevel("MC");
    bfm->SetShuffle(gRunShuffling);
    bfm->SetEventClass("EventPlane");
    bfm->SetDeltaEtaMax(TMath::Abs(gEtaMax-gEtaMin));
    bfm->InitHistograms();
  }
  //========================================================//

  //========================================================//
  //Output objects
  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  TList *fListBF = new TList();
  fListBF->SetName("listBF");
  fListBF->SetOwner();
  
  //Balance function list: shuffling
  TList *fListBFS = 0x0;
  if(gRunShuffling) {
    fListBFS = new TList();
    fListBFS->SetName("listBFShuffled");
    fListBFS->SetOwner();
    fListBFS->Add(bfs->GetHistNp());
    fListBFS->Add(bfs->GetHistNn());
    fListBFS->Add(bfs->GetHistNpn());
    fListBFS->Add(bfs->GetHistNnn());
    fListBFS->Add(bfs->GetHistNpp());
    fListBFS->Add(bfs->GetHistNnp());
  }

  //Balance function list: event mixing
  TList *fListBFM = 0x0;
  if(gRunMixing) {
    fListBFM = new TList();
    fListBFM->SetName("listBFMixed");
    fListBFM->SetOwner();
    fListBFM->Add(bfm->GetHistNp());
    fListBFM->Add(bfm->GetHistNn());
    fListBFM->Add(bfm->GetHistNpn());
    fListBFM->Add(bfm->GetHistNnn());
    fListBFM->Add(bfm->GetHistNpp());
    fListBFM->Add(bfm->GetHistNnp());
  }
  //========================================================//

  //========================================================//
  //Event Mixing
  if(gRunMixing){
    Int_t trackDepth = 50000;
    Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
    
    // centrality bins
    Double_t centralityBins[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* centbins        = centralityBins;
    Int_t nCentralityBins     = sizeof(centralityBins) / sizeof(Double_t) - 1;
    
    // Zvtx bins
    Double_t vertexBins[] = {-10., -7., -5., -3., -1., 1., 3., 5., 7., 10.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* vtxbins     = vertexBins;
    Int_t nVertexBins     = sizeof(vertexBins) / sizeof(Double_t) - 1;
    
    // Event plane angle (Psi) bins
    Double_t psiBins[] = {0.,45.,135.,215.,305.,360.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* psibins     = psiBins;
    Int_t nPsiBins     = sizeof(psiBins) / sizeof(Double_t) - 1;

    AliEventPoolManager *fPoolMgr = 0x0;
    // run the event mixing also in bins of event plane (statistics!)
    if(gRunMixingEventPlane){
      fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins, nPsiBins, psibins);
    }
    else{
      fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);
    }
  }
  //========================================================//

  //========================================================//
  //Create the TChain object: events
  TChain *eventChain = new TChain("events");
  StructEvent tStructEvents;
  eventChain->SetBranchAddress("event",&tStructEvents);
  //========================================================//

  //========================================================//
  //Create the TChain object: particles 
  TChain *particleChain = new TChain("particles");
  StructParticle tStructParticles;
  particleChain->SetBranchAddress("particle",&tStructParticles);
  //========================================================//

  //========================================================//
  //Fill the TChain with the files in the directory
  for(Int_t iFile = 1; iFile <= 9; iFile++) {
    TString filename = aEventDir.Data();
    filename += "/event00"; filename += iFile;
    filename += ".root";
    cout<<"Adding file "<<filename.Data()<<" to the chain..."<<endl;
    eventChain->Add(filename.Data());
    particleChain->Add(filename.Data());
  }
  //========================================================//
 
  //========================================================//
  //Histograms
  //Event stats.
  TString gCutName[5] = {"Total","Offline trigger",
                         "Vertex","Analyzed","sel. Centrality"};
  TH2F *fHistEventStats = new TH2F("fHistEventStats",
				   "Event statistics;;Centrality percentile;N_{events}",
				   5,0.5,5.5,220,-5,105);
  for(Int_t i = 1; i <= 5; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  //Number of accepted particles
  TH2F *fHistNumberOfAcceptedTracks = new TH2F("fHistNumberOfAcceptedTracks",";N_{acc.};Centrality percentile;Entries",4001,-0.5,4000.5,220,-5,105);
  fList->Add(fHistNumberOfAcceptedTracks);
  //========================================================//

  //========================================================//
  //loop over the events
  Int_t nEvents = eventChain->GetEntries();
  cout<<"========================================="<<endl;
  cout<<"Number of events in the chain: "<<nEvents<<endl;
  cout<<"========================================="<<endl;
  Int_t iParticleCounter = 0;
  Int_t nTotalParticles = 0;
  
  //for(Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
  for(Int_t  // for local changed BF configuration
  //gROOT->LoadMacro("./configBalanceFunctionPsiAnalysis.C");
 iEvent = 0; iEvent < 1; iEvent++) {
    eventChain->GetEntry(iEvent);

    Int_t nParticles = tStructEvents.entries;
    cout<<"Event: "<<iEvent+1<<" - ID: "<<tStructEvents.eventID<<" - Entries: "<<tStructEvents.entries<<" - Entries(prev.): "<<tStructEvents.entriesprev<<endl;

    //========================================================//
    //loop over particles
    for(Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      particleChain->GetEntry(nTotalParticles+iParticle);
      
      Double_t gPt = TMath::Sqrt(TMath::Power(tStructParticles.px,2) + TMath::Power(tStructParticles.py,2));
      hPt->Fill(gPt);
      	
      iParticleCounter += 1;
      cout<<"\t Particle counter: "<<iParticleCounter<<" - Particle in the event: "<<iParticle+1<<" - eventID: "<<tStructParticles.eventid<<" - pid: "<<tStructParticles.pid<<" - fatherpid: "<<tStructParticles.fatherpid<<" - rootpid: "<<tStructParticles.rootpid<<" - fathereid: "<<tStructParticles.fathereid<<" - eid: "<<tStructParticles.eid<<endl;
    }//particle loop
    nTotalParticles += nParticles;
  }

  //========================================================//
  //Output file
  TFile *f = TFile::Open("therminator.root","recreate");
  hPdgCode->Write();
  hEta->Write();
  hPt->Write();
  f->Close();  
  //========================================================//

  // Print real and CPU time used for analysis:  
  timer.Stop();
  timer.Print();
}
