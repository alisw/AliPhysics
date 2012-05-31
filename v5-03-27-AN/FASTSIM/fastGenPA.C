AliGenerator*  CreateGenerator();

void fastGenPA(Int_t nev = 1, char* filename = "galice.root")
{
//  Runloader
    AliRunLoader* rl = AliRunLoader::Open("galice.root", "FASTRUN", "recreate");

    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(10000);
    rl->LoadKinematics("RECREATE");
    rl->MakeTree("E");
    gAlice->SetRunLoader(rl);


//  Create stack
    rl->MakeStack();
    AliStack* stack      = rl->Stack();
 
//  Header
    AliHeader* header = rl->GetHeader();


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

//  Generate event

	gener->Generate();

//  Analysis
	Int_t npart = stack->GetNprimary();
	printf("Analyse %d Particles\n", npart);
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = stack->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
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












