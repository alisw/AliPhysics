void pythia6FastjetExample(Int_t nEvent = 50, Float_t e_cms = 2760, Float_t pthardmin = 10, Float_t pthardmax = 50) {

  // example macro for on-the-fly generation of PYTHIA6 events
  // and analysis with new reader interface and FastJet
  // M. van Leeuwen

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");

  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");


  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliGenPythia *pythia=new AliGenPythia(1);
  pythia->SetProcess(kPyJets);
  pythia->SetPtHard(pthardmin, pthardmax);
  pythia->SetEnergyCMS(e_cms);

  pythia->Init();

  //
  // Jet finder settings go via the FastJetHeader
  //
  AliFastJetHeaderV1 *header = new AliFastJetHeaderV1;
  header->SetBGMode(0);
  //  header->SetRadius(0.4);
  header->SetRparam(0.4); 
  //header->SetGhostEtaMax(2);
  //header->SetGhostArea(0.05);
  header->SetAlgorithm(2); // antikt_algorithm = 2, kt = 0 (see fastjet/fastjet/JetDefinition.hh

  AliFastJetFinder *FastJet = new AliFastJetFinder;
  FastJet->SetJetHeader(header);

  // Infrastructure needed for AliGenPythia
  AliStack *stack = new AliStack();
  // There is a big mess with the connection between runloader, AliRun and the gAlice pointer. 
  // This order of things seems to work...
  AliRunLoader *dummyrl = new AliRunLoader("dummyEvtFolder");
  dummyrl->MakeHeader();
  dummyrl->SetEventNumber(0);
  gAlice->SetRunLoader(dummyrl);
  pythia->SetStack(stack);

  // Set up dummy AOD event for JETAN output
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();
  FastJet->ConnectAOD(aod);

  // Set up input structures for FastJet
  AliJetCalTrkEvent JetFinderEvent(0,1);
  TClonesArray aliplist("AliMCParticle",1000);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    stack->Reset();
    pythia->Generate();

    aliplist.Clear();
    JetFinderEvent.Clear();

    Int_t n_part = stack->GetNtrack();
    for (Int_t i_part = 0; i_part < n_part; i_part++) {
      part=(TParticle*) stack->Particle(i_part);

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

// kept for backward compatibility; this was the old function name
void run(Int_t nEvent = 50, Float_t e_cms = 2760) {
  pythia6FastjetExample(nEvent, e_cms);
}
