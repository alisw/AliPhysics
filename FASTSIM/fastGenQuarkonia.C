AliGenerator*  CreateGenerator();

void fastGenQuarkonia(Int_t nev = 1, char* filename = "galice.root")
{
//  Update data base
    AliPDG::AddParticlesToPdgDataBase();
    
// Create the fast tracker
    AliFastMuonTriggerEff *trigeff = new AliFastMuonTriggerEff();
    AliFastMuonTrackingAcc *acc = new AliFastMuonTrackingAcc();
    AliFastMuonTrackingEff *eff = new AliFastMuonTrackingEff();
    AliFastMuonTrackingRes *res = new AliFastMuonTrackingRes();
    acc->SetBackground(0);
    eff->SetBackground(0);
    res->SetBackground(0);  
    acc->Init(); 
    eff->Init(); 
    res->Init(); 
    AliFastDetector* tracker = new AliFastDetector("Tracker", "Muon Tracker");
    tracker->AddResponse(trigeff);
    tracker->AddResponse(acc);
    tracker->AddResponse(eff);
    tracker->AddResponse(res);
    tracker->Init();
    tracker->Dump();
//  Histos
    TH1F* massHU = new TH1F("massHU", "Mass Spectrum:                ", 500, 5, 15.);
    TH1F* massHS = new TH1F("massHS", "Mass Spectrum Smeared:        ", 500, 5, 15.);
    TH2F* etaH   = new TH2F("etaH",   "eta vs etas        ", 100, 0., 5., 100, 0., 5.);
    TH1F* etadH  = new TH1F("etaH",   "eta vs etas        ", 100, -1., 1.);
//
//                        Construction
//
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
//  Generate event
	gener->Generate();
		
//  Analysis
	Int_t npart = stack->GetNprimary();
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = stack->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
	    if (mpart != 553 && mpart != 100553 && mpart != 200553) continue;
	    Int_t ch1 = MPart->GetFirstDaughter();
	    Int_t ch2 = MPart->GetLastDaughter();
	    
	    TParticle *muon1 = stack->Particle(ch1);
	    TParticle *muon2 = stack->Particle(ch2);
	    
	    Float_t theta1 = muon1->Theta();
	    Float_t thetad1 = theta1 * 180./TMath::Pi();
	    Float_t eta1   = muon1->Eta();
	    Float_t pt1    = muon1->Pt();
	    Float_t phi1   = muon1->Phi();
	    Float_t phid1  = phi1 * 180./TMath::Pi() - 180.;
	    Float_t p1     = muon1->P();

	    if (thetad1 > 9. || thetad1 < 2.) continue;

	    Float_t theta2 = muon2->Theta();
	    Float_t thetad2 = theta2 * 180./TMath::Pi();
	    Float_t eta2   = muon2->Eta();
	    Float_t pt2    = muon2->Pt();
	    Float_t phi2   = muon2->Phi();
	    Float_t phid2  = phi2 * 180./TMath::Pi() - 180.;
	    Float_t p2     = muon2->P();


	    if (thetad2 > 9. || thetad2 < 2.) continue;	    
	    
	    Float_t dphi   = phi1 - phi2;
	    Float_t deta   = eta1 - eta2;
	    Float_t m      =  TMath::Sqrt(2. * pt1 * pt2 * (TMath::CosH(deta) - TMath::Cos(dphi)));

	    massHU->Fill(m);
// Smear
	    // the mu+
	    Float_t thetas1, phis1, ps1, thetas2, phis2, ps2, pts1, pts2, etas1, etas2;
	    Float_t trigeffpL, trigeffpH, trigeffnL, trigeffnH, trigeffnA;
	    
	    res->SetCharge(1);
	    eff->SetCharge(1);
	    acc->SetCharge(1);
	    res->Evaluate(p1, thetad1, phid1, ps1, thetas1, phis1);

	    
	    Float_t effp     = eff->Evaluate(pt1, thetad1, phid1);
	    Float_t accp     = acc->Evaluate(pt1, thetad1, phid1);
	    trigeff->Evaluate(1, pt1, thetad1, phid1, trigeffpL, trigeffpH, trigeffnA);
	    thetas1 *= TMath::Pi()/180.;
	    phis1 *= TMath::Pi()/180.;
	    Float_t etas;
	    
	    if (TMath::Tan(thetas1/2.) > 0.) {
		etas  = -TMath::Log(TMath::Tan(thetas1/2.));
		etaH->Fill(etas, eta1);
		etadH->Fill(etas-eta1);
	    }
	    
	    // the mu- 
	    res->SetCharge(-1);
	    acc->SetCharge(-1);
	    eff->SetCharge(-1);
	    res->Evaluate(p2, thetad2, phid2, ps2, thetas2, phis2);
	    Float_t effn     = eff->Evaluate(pt2, thetad2, phid2);
	    Float_t accn     = acc->Evaluate(pt2, thetad2, phid2);
	    trigeff->Evaluate(-1, pt2, thetad2, phid2, trigeffnL, trigeffnH, trigeffnA);
	    thetas2 *= TMath::Pi()/180.;
	    phis2 *= TMath::Pi()/180.;
	    if (TMath::Tan(thetas2/2.) > 0) {
		etas  = -TMath::Log(TMath::Tan(thetas2/2.));
		etaH->Fill(etas, eta2);
		etadH->Fill(etas-eta2);
	    }
	    
	    dphi   = phis1 - phis2;
	    etas1  = - TMath::Log(TMath::Abs(TMath::Tan(thetas1/2.)+1e-4));
	    etas2  = - TMath::Log(TMath::Abs(TMath::Tan(thetas2/2.)+1e-4));	    
	    deta   = etas1 - etas2;
	    pts1   = ps1 * TMath::Sin(thetas1);
	    pts2   = ps2 * TMath::Sin(thetas2);	    
	    
	    m      =  TMath::Sqrt(2. * pts1 * pts2 * (TMath::CosH(deta) - TMath::Cos(dphi)));
	    Float_t wgt = effn * effp * accn * accp * trigeffpH * trigeffnH;
	    printf("Mass %f\n", m);
	    if (pts1 > 0. && pts2 > 0.)
	    massHS->Fill(m, wgt);
	}
	
//  Finish event
//      I/O
//	
	stack->FinishEvent();
	header->SetStack(stack);
	rl->TreeE()->Fill();
	header->SetNprimary(stack->GetNprimary());
	header->SetNtrack(stack->GetNtrack());  
	rl->WriteKinematics("OVERWRITE");
    } // event loop
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    massHU->Draw();
    massHS->Draw("same");

    TCanvas *c2 = new TCanvas("c2","Canvas 2",400,10,600,700);
    etaH->Draw();
    
//
//                         Termination
//  Generator
    gener->FinishRun();
//  Stack
    stack->FinishRun();
//  Write file
    rl->WriteHeader("OVERWRITE");
    gener->Write();
    rl->Write();
}


AliGenerator*  CreateGenerator()
{
    AliGenParam *gener 
	= new AliGenParam(2000, AliGenMUONlib::kUpsilon, "");
        
    gener->SetMomentumRange(0,999);
    gener->SetPtRange(0,100.);
    gener->SetPhiRange(0., 360.);
    gener->SetYRange(2.5,4);
    gener->SetCutOnChild(1);
    gener->SetChildThetaRange(2.0,9.);
    gener->SetForceDecay(kDiMuon);

    AliDecayer *decayer = new AliDecayerPythia();
    decayer->SetForceDecay(kAll);
    decayer->Init();
    gener->SetDecayer(decayer);
    gener->Init();
    
    gener->Draw();
    
    return gener;
}












