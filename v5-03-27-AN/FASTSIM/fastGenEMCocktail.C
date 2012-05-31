void fastGenEMCocktail(Int_t nev = 1, char* filename = "galice.root")
{
// load libraries

  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6.so");  // ALICE specific implementations

// Runloader
    
  AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");
    
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(nev);
  rl->LoadKinematics("RECREATE");
  rl->MakeTree("E");
  gAlice->SetRunLoader(rl);

// Create stack
  rl->MakeStack();
  AliStack* stack      = rl->Stack();
 
// Header
  AliHeader* header = rl->GetHeader();

// Create and Initialize Generator

  AliGenEMCocktail *gener = new AliGenEMCocktail();

//=======================================================================
// Set External decayer
  TVirtualMCDecayer *decayer = new AliDecayerPythia();

  gener->SetNPart(1000);               // source multiplicity per event
  gener->SetPtRange(0.,20.);
  gener->SetYRange(-1.,1.);
  gener->SetPhiRange(0., 360.);
  gener->SetOrigin(0.,0.,0.); 
  gener->SetSigma(0.,0.,0.);
  gener->SetVertexSmear(kPerEvent);
  gener->SetTrackingFlag(0);
  gener->SetDecayMode(kElectronEM);    // select cocktail:
                                       // kElectronEM   => single electron
                                       // kDiElectronEM => electron-positron
                                       // kGammaEM      => single photon
  gener->SetDecayer(decayer);
  gener->SetWeightingMode(kNonAnalog); // select weighting:
                                       // kNonAnalog => weight ~ dN/dp_T
                                       // kAnalog    => weight ~ 1
  gener->CreateCocktail();
  gener->Init();
  gener->SetStack(stack);

//
//                        Event Loop
//
  Int_t iev;
     
  for (iev = 0; iev < nev; iev++) {

    printf("\n \n Event number %d \n \n", iev);
	
// Initialize event
    header->Reset(0,iev);
    rl->SetEventNumber(iev);
    stack->Reset();
    rl->MakeTree("K");
//  stack->ConnectTree();
    
// Generate event
    gener->Generate();
// Analysis
    Int_t npart = stack->GetNprimary();
    printf("Analyse %d Particles\n", npart);
    for (Int_t part=0; part<npart; part++) {
      TParticle *MPart = stack->Particle(part);
      Int_t mpart  = MPart->GetPdgCode();
//      printf("Particle %d\n", mpart);
    }
	
// Finish event
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());  

// I/O
//	
    stack->FinishEvent();
    header->SetStack(stack);
    rl->TreeE()->Fill();
    rl->WriteKinematics("OVERWRITE");

  } // event loop

//
//                         Termination
// Generator
  gener->FinishRun();
// Write file
  rl->WriteHeader("OVERWRITE");
  gener->Write();
  rl->Write();
}












