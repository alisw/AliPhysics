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
//Primordial particles: fathereid==-1 ==> 
//pid==fatherpid==rootpid
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
    
    AliEventPoolManager *fPoolMgr = 0x0;
    fPoolMgr = new AliEventPoolManager(poolsize, trackDepth);
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
  TH1F *fHistEventStats = new TH1F("fHistEventStats",
				   "Event statistics;;N_{events}",
				   5,0.5,5.5);
  for(Int_t i = 1; i <= 5; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  //Number of accepted particles
  TH1F *fHistNumberOfAcceptedTracks = new TH1F("fHistNumberOfAcceptedTracks",";N_{acc.};Entries",4001,-0.5,4000.5);
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
  for(Int_t iEvent = 0; iEvent < 1; iEvent++) {
    eventChain->GetEntry(iEvent);

    //========================================================//
    //Create the TObjArray object to store the particles
    TObjArray* tracksAccepted = new TObjArray;
    tracksAccepted->SetOwner(kTRUE);
    //========================================================//

    Int_t nParticles = tStructEvents.entries;
    if(((iEvent+1)%100)==0)
      cout<<"Event: "<<iEvent+1<<" - ID: "<<tStructEvents.eventID<<" - Entries: "<<tStructEvents.entries<<endl;

    //========================================================//
    //Filling event level histos
    fHistEventStats->Fill(1);
    fHistEventStats->Fill(2);
    fHistEventStats->Fill(3);
    fHistEventStats->Fill(4);
    fHistEventStats->Fill(5);

    Int_t gNumberOfAcceptedParticles = 0;
    Double_t gReactionPlane = 0.;
    //========================================================//
    //loop over particles
    for(Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      particleChain->GetEntry(nTotalParticles+iParticle);
      
      //========================================================//
      //consider only primordial particles
      if(tStructParticles.fathereid != -1) continue;

      //========================================================//
      //Calculate kinematic variables
      Double_t gPt = TMath::Sqrt(TMath::Power(tStructParticles.px,2) + 
				 TMath::Power(tStructParticles.py,2));
      Double_t gP = TMath::Sqrt(TMath::Power(tStructParticles.px,2) + 
				TMath::Power(tStructParticles.py,2) + 
				TMath::Power(tStructParticles.pz,2) );
      Double_t gEta = -100.;
      if(gP != tStructParticles.pz)
	gEta = 0.5*TMath::Log((gP + tStructParticles.pz)/(gP - tStructParticles.pz));
      Double_t gPhi = 0.;
      Double_t gCharge = 0.;
      
      //========================================================//
      //Apply cuts
      if((gEta > gEtaMax)||(gEta < gEtaMin)) continue;
      if((gPt > gPtMax)||(gPt < gPtMin)) continue;
      
      tracksAccepted->Add(new AliBFBasicParticle(gEta,gPhi,gPt,gCharge, 1.));
      gNumberOfAcceptedParticles += 1;
      
      iParticleCounter += 1;
      //cout<<"\t Particle counter: "<<iParticleCounter<<" - Particle in the event: "<<iParticle+1<<" - eventID: "<<tStructParticles.eventid<<" - pid: "<<tStructParticles.pid<<" - fatherpid: "<<tStructParticles.fatherpid<<" - rootpid: "<<tStructParticles.rootpid<<" - fathereid: "<<tStructParticles.fathereid<<" - eid: "<<tStructParticles.eid<<endl;
    }//particle loop
    
    //========================================================//
    // Event mixing (borrowed code from the task) 
    /*if (gRunMixing) {
      AliEventPool* pool = fPoolMgr->GetEventPool(gCentrality, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane);
      
      if (!pool) {
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f, psi = %f", gCentrality, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane));
      }
      else {
      //pool->SetDebug(1);
      if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5){ 
      
      Int_t nMix = pool->GetCurrentNEvents();
      //cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
      
      //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
      //((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
      //if (pool->IsReady())
      //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3);
      
      // Fill mixed-event histos here  
      for (Int_t jMix=0; jMix<nMix; jMix++) {
      TObjArray* tracksMixed = pool->GetEvent(jMix);
      fMixedBalance->CalculateBalance(gReactionPlane,tracksMain,tracksMixed,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
      }
      }
      
      // Update the Event pool
      pool->UpdatePool(tracksMain);
      //pool->PrintInfo();
      
      }//pool NULL check  
      }*///run mixing
      //========================================================//

    //========================================================//
    // calculate balance function
    //fBalance->CalculateBalance(gReactionPlane,tracksMain,NULL,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());

    fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedParticles);
    nTotalParticles += nParticles;
  }

  //========================================================//
  //Output file
  TFile *f = TFile::Open("AnalysisResults.root","recreate");
  fList->Write();
  fListBF->Write();
  if(gRunMixing) fListBFM->Write();
  f->Close();  
  //========================================================//

  // Print real and CPU time used for analysis:  
  timer.Stop();
  timer.Print();
}
