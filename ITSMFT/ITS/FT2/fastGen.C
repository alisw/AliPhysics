
AliGenerator *Hijing();
AliGenerator *Hijing2500();
AliGenerator *Hijing2500HF();
static Float_t energy = 5500.; // energy in CMS
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static TString comment;
AliDecayer *decayer = 0;

void fastGen(Int_t nev = 1,UInt_t seed = 0)
{
	gSystem->Load("libITSUpgradeBase.so");
	gSystem->Load("libITSUpgradeSim.so");
	gSystem->Load("libEVGEN");
	//  Load libraries
	gSystem->Load("liblhapdf");
	gSystem->Load("libpythia6");
	gSystem->Load("libEGPythia6");
	//	gSystem->Load("libpythia6_4_21");   // Pythia 6.4
	gSystem->Load("libAliPythia6");
	//	gSystem->Load("libhijing");  already loaded
	gSystem->Load("libTHijing");
	//  Runloader
	gRandom->SetSeed(seed);
	printf("### Seed of random number generation ###\n###            %0.f             ###\n",seed);
	
	AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");
	
	rl->SetCompressionLevel(2);
	rl->SetNumberOfEventsPerFile(nev);
	rl->LoadKinematics("RECREATE");
	rl->MakeTree("E");
	gAlice->SetRunLoader(rl);
	
	decayer = new AliDecayerPythia();
	//	decayer->SetForceDecay(kAll);
	decayer->SetForceDecay(kHadronicDWithout4Bodies);
	decayer->Init();
	
	//  Create stack
	rl->MakeStack();
	AliStack* stack      = rl->Stack();
	//
	//  Header
	AliHeader* header = rl->GetHeader();
	//
	//  Create and Initialize Generator
	AliGenerator *gener = Hijing2500HF();
	//
	//
	// Size of the interaction diamond
	// Longitudinal
	Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
	//
	// Transverse
	Float_t betast  = 3.5;												// beta* [m]
	Float_t eps     = 3.75e-6;										// emittance [m]
	Float_t gamma   = energy / 2.0 / 0.938272;		// relativistic gamma [1]
	Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
	printf("###         Diamond size             ###\n###   x-y: %10.3e z: %10.3e  ###\n \n", sigmaxy, sigmaz);
	gener->SetSigma(sigmaxy, sigmaxy, sigmaz);    // Sigma in (X,Y,Z) (cm) on IP position
	gener->SetVertexSmear(kPerEvent);
	gener->Init();
	gener->SetStack(stack);
	//
	//	Event Loop
	//
	Int_t iev;
	
	for (iev = 0; iev < nev; iev++) {
		//
		printf("\n \n Event number %d \n \n", iev);
		//
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
			//	printf("Particle %d\n", mpart);
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

AliGenerator* Hijing()
{
	AliGenHijing *gener = new AliGenHijing(-1);
	// centre of mass energy
	gener->SetEnergyCMS(energy);
	gener->SetImpactParameterRange(0., 0.);//5.); // as in LHC13d19
	// reference frame
	gener->SetReferenceFrame("CMS");
	// projectile
	gener->SetProjectile("A", 208, 82);
	gener->SetTarget    ("A", 208, 82);
	// tell hijing to keep the full parent child chain
	gener->KeepFullEvent();
	// enable jet quenching
	gener->SetJetQuenching(1);
	// enable shadowing
	gener->SetShadowing(1);
	// Don't track spectators
	gener->SetSpectators(0);
	// kinematic selection
	gener->SetSelectAll(0);
	return gener;
}
AliGenerator* Hijing2500()
{
	AliGenHijing *gener = (AliGenHijing*) Hijing();
	gener->SetJetQuenching(0);
	gener->SetPtHardMin (3.7);
	return gener;
}

AliGenerator* Hijing2500HF()
{
	comment = comment.Append(" PbPb: Hjing2500 at 5.5 + ITS Upgrade signals");
	AliPDG::AddParticlesToPdgDataBase();
	
	AliGenCocktail *cocktail = new AliGenCocktail();
	
	cocktail->SetProjectile("A", 208, 82);
	cocktail->SetTarget    ("A", 208, 82);
	cocktail->SetEnergyCMS(energy);
	//
	// 1 Hijing event
	TFormula* one    = new TFormula("one",    "1.");
	// provides underlying event and collision geometry
	AliGenHijing *hijing = Hijing2500();
	cocktail->AddGenerator(hijing,"hijing",1);
	Float_t thminH = (180./TMath::Pi())*2.*atan(exp(-2.5));
	Float_t thmaxH = (180./TMath::Pi())*2.*atan(exp( 2.5));
	hijing->SetChildThetaRange(thminH,thmaxH);
	hijing->SetDecaysOff(0);  // HIJIJNG decays all particles, including weak strange decayss
	//
	// Set pseudorapidity range from -1. to 1.
	//
	
	Float_t thmin          = (180./TMath::Pi())*2.*atan(exp(-1.));
	Float_t thmax          = (180./TMath::Pi())*2.*atan(exp( 1.));
	//
// as in LHC13d19
	AliGenParam *gen[14];
	UInt_t partId[7] = {AliGenITSULib::kLc,AliGenITSULib::kLb,AliGenITSULib::kXi_c,AliGenITSULib::kBplus, AliGenITSULib::kBzero, AliGenITSULib::kDs, AliGenITSULib::kDplus};
	for(Int_t iPart=0; iPart<14 ; iPart++){
		if(iPart%2==0) gen[iPart] = new AliGenParam(15,new AliGenITSULib(),partId[iPart/2],"DIST");
		if(iPart%2==1) gen[iPart]= new AliGenParam(15,new AliGenITSULib(),-partId[iPart/2],"DIST");
		gen[iPart]->SetDecayer(decayer);
		gen[iPart]->SetPtRange(0.,999.);
		gen[iPart]->SetPhiRange(0., 360.);
		gen[iPart]->SetYRange(-1.,1.);
		gen[iPart]->SetCutOnChild(1);
		gen[iPart]->SetChildThetaRange(thmin,thmax);
		gen[iPart]->SetSelectAll(kTRUE);
		gen[iPart]->SetForceDecay(kBeautyUpgrade);
		gen[iPart]->SetPreserveFullDecayChain(kTRUE);
		cocktail->AddGenerator(gen[iPart], Form("Generator_%i_%i",partId[iPart/2],iPart%2), 1);
	}
 
	// Hypernuclei: 10 per type for 3LH, 4LH, 4LHe
	AliGenBox *pG1=new AliGenBox(5);
	pG1->SetPart(1010010030);
	pG1->SetPtRange(0,10);
	pG1->SetPhiRange(0,360);
	pG1->SetYRange(-1,1);
	cocktail->AddGenerator(pG1,"g1",1);
	
	AliGenBox *pG2=new AliGenBox(5);
	pG2->SetPart(-1010010030);
	pG2->SetPtRange(0,10);
	pG2->SetPhiRange(0,360);
	pG2->SetYRange(-1,1);
	cocktail->AddGenerator(pG2,"g2",1);
	
	AliGenBox *pG3=new AliGenBox(5);
	pG3->SetPart(1010010040);
	pG3->SetPtRange(0,10);
	pG3->SetPhiRange(0,360);
	pG3->SetYRange(-1,1);
	cocktail->AddGenerator(pG3,"g3",1);
	
	AliGenBox *pG4=new AliGenBox(5);
	pG4->SetPart(-1010010040);
	pG4->SetPtRange(0,10);
	pG4->SetPhiRange(0,360);
	pG4->SetYRange(-1,1);
	cocktail->AddGenerator(pG4,"g4",1);
	
	AliGenBox *pG5=new AliGenBox(5);
	pG5->SetPart(1010020040);
	pG5->SetPtRange(0,10);
	pG5->SetPhiRange(0,360);
	pG5->SetYRange(-1,1);
	cocktail->AddGenerator(pG5,"g5",1);
	
	AliGenBox *pG6=new AliGenBox(5);
	pG6->SetPart(-1010020040);
	pG6->SetPtRange(0,10);
	pG6->SetPhiRange(0,360);
	pG6->SetYRange(-1,1);
	cocktail->AddGenerator(pG6,"g6",1);
	
	/*
	 AliGenParam *gen[2];
	 UInt_t partId[1] = {AliGenITSULib::kBplus};
	 for(Int_t iPart=0; iPart<2 ; iPart++){
		if(iPart%2==0) gen[iPart] = new AliGenParam(100,new AliGenITSULib(),partId[iPart/2],"DIST");
		if(iPart%2==1) gen[iPart]= new AliGenParam(100,new AliGenITSULib(),-partId[iPart/2],"DIST");
		gen[iPart]->SetDecayer(decayer);
		gen[iPart]->SetPtRange(0.,999.);
		gen[iPart]->SetPhiRange(0., 360.);
		gen[iPart]->SetYRange(-1.,1.);
		gen[iPart]->SetCutOnChild(1);
		gen[iPart]->SetChildThetaRange(thmin,thmax);
		gen[iPart]->SetSelectAll(kTRUE);
		gen[iPart]->SetForceDecay(kBeautyUpgrade);
		gen[iPart]->SetPreserveFullDecayChain(kTRUE);
		cocktail->AddGenerator(gen[iPart], Form("Generator_%i_%i",partId[iPart/2],iPart%2), 1);
	 }
	 */

/*
	AliGenBox *pG1=new AliGenBox(50);
	pG1->SetPart(521);
	pG1->SetPtRange(0.5,24);
	pG1->SetPhiRange(0,360);
	pG1->SetYRange(-1,1);
	pG1->Init();
	
	AliGenBox *pG2=new AliGenBox(50);
	pG2->SetPart(-521);
	pG2->SetPtRange(0.5,24);
	pG2->SetPhiRange(0,360);
	pG2->SetYRange(-1,1);
	pG2->Init();
	
	printf("###     Adding particle decayer      ###\n");
	
	AliGenEvtGen *gene = new AliGenEvtGen();
	gene->SetUserDecayTable("./BTODPI.DEC");
	gene->SetParticleSwitchedOff(AliGenEvtGen::kBeautyPart);
	cocktail->AddGenerator(pG1,"OnTheFly_gbx_Bp",1);
	cocktail->AddGenerator(pG2,"OnTheFly_gbx_Bm",1);
	cocktail->AddGenerator(gene,"OnTheFly_evtGen",1.);
	*/
	return cocktail;
}
