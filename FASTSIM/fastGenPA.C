AliGenerator*  CreateGenerator();

void fastGenPA(Int_t nev = 1, char* filename = "galice.root")
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
    AliPythia* pyth = AliPythia::Instance();
    
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
	stack->Reset();
	stack->BeginEvent(iev);

//  Generate event
	gener->Generate();
		
//  Analysis
	Int_t npart = stack->GetNprimary();
//	printf("Analyse %d Particles\n", npart);
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
//	(stack->TreeK())->Write(0,TObject::kOverwrite);
    } // event loop
//
//                         Termination
//  Generator
    printf("Calling Finish Run \n");
    
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
    AliGenCocktail *gener  = new AliGenCocktail();
    gener->SetTrackingFlag(0);
    AliGenHijing   *hijing = new AliGenHijing(-1);
// centre of mass energy 
    hijing->SetEnergyCMS(TMath::Sqrt(82./208.) * 14000.);
// impact parameter range
    hijing->SetImpactParameterRange(0., 6.);
// reference frame
    hijing->SetReferenceFrame("CMS");
    hijing->SetBoostLHC(1);
// projectile
    hijing->SetProjectile("P", 1, 1);
    hijing->SetTarget    ("A", 208, 82);
// tell hijing to keep the full parent child chain
    hijing->KeepFullEvent();
// enable jet quenching
    hijing->SetJetQuenching(4);
// enable shadowing
    hijing->SetShadowing(1);
// Don't track spectators
    hijing->SetSpectators(0);
// kinematic selection
    hijing->SetSelectAll(0);
//
    AliGenSlowNucleons*  gray    = new AliGenSlowNucleons(1);
    AliSlowNucleonModel* model   = new AliSlowNucleonModelExp();
    gray->SetSlowNucleonModel(model);
    gray->SetDebug(1);
    
    
    gener->AddGenerator(hijing,"Hijing pPb", 1);
    gener->AddGenerator(gray,  "Gray Particles",1);
    
    return gener;
}












