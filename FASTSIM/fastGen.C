AliGenerator*  CreateGenerator();

void fastGen(Int_t nev = 1, char* filename = "galice.root")
{
//
//                        Construction
//
//  Output file
    TFile*  file         = new TFile(filename, "recreate");
//  Create stack
    AliStack* stack      = new AliStack(10000);
    stack->MakeTree(0, filename);

//  Create Header
    AliHeader* header    = new AliHeader();
//  Create Header Tree
    TTree* treeE         = new TTree("TE","Headers");
    treeE->Branch("Header", "AliHeader", &header, 4000, 0);
    treeE->Write();
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
	stack->BeginEvent(iev);

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
//	stack->FinishEvent();
//	header->SetStack(stack);
//	treeE->Fill();
//	Reset stack

	stack->Reset();
    } // event loop
//
//                         Termination
//  Generator
    gener->FinishRun();
//  Header
    treeE->Write(0,TObject::kOverwrite);
    delete treeE;   treeE = 0;
//  Stack
    stack->FinishRun();
//  Write file
    gener->Write();
    file->Write();
}


AliGenerator*  CreateGenerator()
{
    gener = new AliGenPythia(1);
//
//
//   vertex position and smearing 
    gener->SetVertexSmear(kPerEvent);
//   structure function
    gener->SetStrucFunc(kGRVHO);
//   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
    gener->SetProcess(kPyJets);
//   Centre of mass energy 
    gener->SetEnergyCMS(14000.);
//   Pt transfer of the hard scattering
    gener->SetPtHard(5.,5.1);
//   Initialize generator    
    return gener;
}












