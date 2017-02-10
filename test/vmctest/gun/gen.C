// $Id$
//
// Generation of kinematics tree which can be then used as an extrenal
// event generator.
// According to: $ALICE_ROOT/test/genkine/gen/fastGen.C
//
// By I. Hrivnacova, IPN Orsay

void gen(Int_t nev = 1)
{

  AliPDG::AddParticlesToPdgDataBase();
  TDatabasePDG::Instance();

  // Run loader
  AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");
  
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);
  
  //  Create stack
  rl->MakeStack();
  AliStack* stack = rl->Stack();
  
  //  Header
  AliHeader* header = rl->GetHeader();
  
  //  Create and Initialize Generator
  AliGenerator* gener = genGunConfig();

  // Go to galice.root
  rl->CdGAFile();

  // Forbid some decays. Do it after gener->Init(0, because
  // the initialization of the generator includes reading of the decay table.
  // ...

  //
  // Event Loop
  //
  
  TStopwatch timer;
  timer.Start();
  for (Int_t iev = 0; iev < nev; iev++) {
    
    cout <<"Event number "<< iev << endl;
    
    // Initialize event
    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
    
    // Generate event
    stack->Reset();
    stack->ConnectTree(rl->TreeK());
    gener->Generate();
    cout << "Number of particles " << stack->GetNprimary() << endl;
    
    // Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  
    
    // I/O
    stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");
    
  } // event loop
  timer.Stop();
  timer.Print();
  
  //                         Termination
  //  Generator
  gener->FinishRun();
  //  Write file
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();
}
