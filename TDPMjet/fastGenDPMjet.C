AliGenerator*  CreateGenerator();

void fastGenDPMjet(Int_t nev = 1, char* filename = "dpmjet.root")
{
//  Runloader
#if defined(__CINT__)
    gSystem->Load("liblhapdf"); 
    gSystem->Load("libpythia6");     
    gSystem->Load("libdpmjet.so");
    gSystem->Load("libTDPMjet.so");
#endif

    AliPDG::AddParticlesToPdgDataBase();
    
    AliRunLoader* rl = AliRunLoader::Open(filename,"FASTRUN","recreate");
    AliPythiaRndm::SetPythiaRandom(new TRandom());

    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(nev);
    rl->LoadKinematics("RECREATE");
    rl->MakeTree("E");
    gAlice->SetRunLoader(rl);

    //  Create stack
    rl->MakeStack();
    AliStack* stack   = rl->Stack();
 
    //  Header
    AliHeader* header = rl->GetHeader();
    //AliGenDPMjetEventHEader* dpmHeader;

    //  Create and Initialize Generator
    AliGenDPMjet *gener = CreateGenerator();
    //AliCollisionGeometry *coll = gener->CollisionGeometry();
    gener->Init();
    gener->SetStack(stack);
    
    // Event Loop
    Int_t iev;
     
    for (iev = 0; iev < nev; iev++) {

	if(iev%500==0) printf("\n  Event number %d  \n", iev);
	
	//  Initialize event
	header->Reset(0,iev);
	rl->SetEventNumber(iev);
	stack->Reset();
	rl->MakeTree("K");
    
	//  Generate event
	gener->Generate();
	
	Int_t npart = stack->GetNprimary();
        
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

    gener->FinishRun();

    rl->WriteHeader("OVERWRITE");
    gener->Write();
    rl->Write();
    
}


AliGenerator*  CreateGenerator()
{ 
  AliGenDPMjet *gener = new AliGenDPMjet(-1);
  
  //  p-A
  /*Float_t pEnergy = 4000.;
  gener->SetProjectile("P", 1, 1);
  gener->SetTarget("A", 208, 82);
  //
  gener->SetEnergyCMS(TMath::Sqrt(82./208.) * 2* pEnergy);
  gener->SetProjectileBeamEnergy(pEnergy); */ 
  
  //  A-p
  Float_t pEnergy = 1577.;
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget("P", 1, 1);
  //
  gener->SetEnergyCMS(TMath::Sqrt(208./82.) * 2* pEnergy);
  gener->SetProjectileBeamEnergy(pEnergy);

  // Pb-Pb
  /*  gener->SetProjectile("A", 208, 82);
  gener->SetTarget("A", 208, 82);
  gener->SetEnergyCMS(2760.);
  //
  gener->SetProcess(kDpmMb);
  gener->SetImpactParameterRange(0., 100.);
  //gener->SetFragmentProd(kTRUE);*/

  gener->SetTrackingFlag(0);
  
  return gener;
}












