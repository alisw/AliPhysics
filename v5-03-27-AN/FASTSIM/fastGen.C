AliGenerator*  CreateGenerator();

void fastGen(Int_t nev = 1, char* filename = "galice.root")
{
//  Load libraries
  gSystem->Load("liblhapdf.so");
  gSystem->Load("libpythia6.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");

//  Runloader
    
    AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");
    
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
    AliGenerator *gener = CreateGenerator();
    gener->Init();
    gener->SetStack(stack);
    
//
//                        Event Loop
//
    Int_t iev;
     
    for (iev = 0; iev < nev; iev++) {

	printf("\n \n Event number %d \n \n", iev);
	
//  Initialize event
	header->Reset(0,iev);
	rl->SetEventNumber(iev);
	stack->Reset();
	rl->MakeTree("K");
//	stack->ConnectTree();
    
//  Generate event
	gener->Generate();
//  Analysis
	Int_t npart = stack->GetNprimary();
	printf("Analyse %d Particles\n", npart);
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = stack->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
	    printf("Particle %d\n", mpart);
	}
	
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


AliGenerator*  CreateGenerator()
{
    gener = new AliGenPythia(1);
//
//
//   vertex position and smearing 
    gener->SetVertexSmear(kPerEvent);
//   structure function
    gener->SetStrucFunc(kCTEQ6);
//   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
    gener->SetProcess(kPyJets);
//   Centre of mass energy 
    gener->SetEnergyCMS(5500.);
//   Pt transfer of the hard scattering
    gener->SetPtHard(50.,50.2);
//   Initialize generator    
    return gener;
}












