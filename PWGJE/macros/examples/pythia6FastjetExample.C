void run(Int_t nEvent = 50, Float_t e_cms = 2760) {

  // example macro for on-the-fly generation of PYTHIA6 events
  // and analysis with new reader interface and FastJet
  // M. van Leeuwen

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libSISConePlugin");
  gSystem->Load("libJETANdev");
  gSystem->Load("libFASTJETANdev");

  gSystem->Load("libpythia6.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");


  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliPythia6 *pythia=AliPythia6::Instance();

  pythia->SetCKIN(3,10);   // minimum hard pt
  pythia->SetCKIN(4,1000);  // maximum hard pt

  pythia->SetMDCY(pythia->Pycomp(111),1,0);  // switch off pi0 decay

  pythia->Initialize("CMS","p","p",e_cms);


  AliFastJetHeaderV1 *header = new AliFastJetHeaderV1;
  header->SetBGMode(0);
  header->SetRadius(0.4);
  //header->SetGhostEtaMax(2);
  //header->SetGhostArea(0.05);
  header->SetAlgorithm(2); // antikt_algorithm = 2, kt = 0 (see fastjet/fastjet/JetDefinition.hh

  AliFastJetFinder *FastJet = new AliFastJetFinder;
  FastJet->SetJetHeader(header);

  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();
  FastJet->ConnectAOD(aod);

  AliJetCalTrkEvent JetFinderEvent(0,1);
  TClonesArray *plist = new TClonesArray("TParticle");
  TClonesArray aliplist("AliMCParticle",1000);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    pythia->GenerateEvent();

    pythia->GetParticles(plist);

    aliplist.Clear();
    JetFinderEvent.Clear();

    Int_t n_part = plist->GetEntries();
    for (Int_t i_part = 0; i_part < n_part; i_part++) {
      part=(TParticle*) plist->At(i_part);

      if (part->GetStatusCode() >= 10)  // Not a final state particle
	continue;

      new (aliplist[i_part]) AliMCParticle(part);

      JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)aliplist[i_part],1,1);
    }

    aod->ClearStd();
    FastJet->Reset();
    FastJet->SetCalTrkEvent(JetFinderEvent);
    FastJet->ProcessEvent();
    if (aod->GetNJets() > 0) {
      cout << "event " << iEvent << " " << aod->GetNJets() << " jets found" << endl;
      for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
	jet = aod->GetJet(iJet);
	cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl; 
      }
    }
  }
}
