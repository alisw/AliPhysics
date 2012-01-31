typedef enum {kPhojet = -1,kPyTuneCDFA=100,kPyTuneAtlasCSC=306, kPyTuneCMS6D6T=109, kPyTunePerugia0=320 } Tune_t;

AliGenerator*  CreateGenerator(Tune_t tune = kPyTuneCDFA , Float_t energy);

void fastGen(Tune_t tune = kPyTuneCDFA , Float_t energy, Int_t nev = 1, TString process)
{
  // Add all particles to the PDG database
  AliPDG::AddParticlesToPdgDataBase();

  // set the random seed
  TDatime date;
  UInt_t seed    = date.Get()+gSystem->GetPid();
  gRandom->SetSeed(seed);
  cout<<"Seed for random number generation= "<<seed<<endl; 


  //  Runloader  
  AliRunLoader* rl = AliRunLoader::Open("galice.root", "FASTRUN","recreate");
    
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);

  //  Create stack
  rl->MakeStack();
  AliStack* stack      = rl->Stack();
 
  //  Header
  AliHeader* header = rl->GetHeader();
  //
  //  Create and Initialize Generator
  AliGenerator *gener = CreateGenerator(tune,energy);
  gener->Init();
  // if nsd switch off single diffraction
  if ( process == "NSD"){
    if(tune != kPhojet) {
      AliPythia::Instance()->	SetMSUB(92,0);             // single diffraction AB-->XB
      AliPythia::Instance()-> SetMSUB(93,0);             // single diffraction AB-->AX
    }
    else {
      cout << "NSD not yet implemented in the phojet case" << endl;
      exit(1);
    }
  }
  gener->SetStack(stack);
    
  //
  //                        Event Loop
  //
  Int_t iev;
     
  for (iev = 0; iev < nev; iev++) {

    if(!(iev%500)) printf("\n \n Event number %d \n \n", iev);
	
    //  Initialize event
    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
    //	stack->ConnectTree();
    
    //  Generate event
    gener->Generate();
    //  Analysis
    // 	Int_t npart = stack->GetNprimary();
    // 	printf("Analyse %d Particles\n", npart);
    // 	for (Int_t part=0; part<npart; part++) {
    // 	    TParticle *MPart = stack->Particle(part);
    // 	    Int_t mpart  = MPart->GetPdgCode();
    // 	    printf("Particle %d\n", mpart);
    // 	}
	
    //  Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  
    //      I/O
    //	
    stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");

  } // event loop
    //
    //                         Termination
    //  Generator
  gener->FinishRun();
  //  Write file
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();
    
}


AliGenerator*  CreateGenerator(Tune_t tune, Float_t energy)
{

  

  if (tune == -1) {

    // phojet
    AliGenDPMjet* gener = new AliGenDPMjet(1);
    gener->SetProcess(kDpmMb);
    gener->SetProjectile("P", 1, 1);
    gener->SetTarget("P", 1, 1);

    gener->SetEnergyCMS(energy);
    return gener;
  }

  if (tune != kPyTuneAtlasCSC && tune != kPyTuneCDFA && tune != kPyTuneCMS6D6T && tune != kPyTunePerugia0) {
    
    Printf("Unknown pythia tune, quitting");
    exit(1);

  } else {

    AliGenPythia * gener =  new AliGenPythia(1);
    //
    //  set pythia tune
    gener->SetTune(tune);

    //   structure function  
    if(tune == kPyTuneAtlasCSC) {
      cout << "Setting structure function" << endl;      
      gener->SetStrucFunc(kCTEQ61);
    }
    if(tune == kPyTuneCMS6D6T) {
      cout << "Setting structure function" << endl;      
      gener->SetStrucFunc(kCTEQ6l);
    }
    if(tune == kPyTunePerugia0) {
      cout << "Setting new parton shower" << endl;      
      gener->UseNewMultipleInteractionsScenario();
    }
    //   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
    gener->SetProcess(kPyMb);
    //   Centre of mass energy 
    gener->SetEnergyCMS(energy);

    // Set target/projectile // Is this needded?
    gener->SetProjectile("P", 1, 1);
    gener->SetTarget("P", 1, 1);

    return gener;
  }
}












