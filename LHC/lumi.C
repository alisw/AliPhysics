void lumi () 
{
// Dynamically link some shared libs
    gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libLHC");

//
//  LHC 
//    
    AliLHC* lhc = new AliLHC();
    lhc->SetRadius(2665887./2./TMath::Pi());
    lhc->SetAverageBeta(5050.);
    lhc->SetAverageDisp(129.);
    lhc->SetSetUpTime(1./3. * 3600.);
    lhc->SetFillingTime(3.  * 3600.);
//
//  The Interaction regions
//
//

    Float_t bstar = 50.;
//
    AliLhcIRegion* alice =  new AliLhcIRegion(lhc, "ALICE", "ALICE");
    alice->SetBetaStar(bstar);
    AliLhcIRegion* cms   =  new AliLhcIRegion(lhc, "CMS"  , "CMS"  );
    cms->SetBetaStar(bstar);
    AliLhcIRegion* atlas =  new AliLhcIRegion(lhc, "ATLAS", "ATLAS"); 
    atlas->SetBetaStar(bstar);

    lhc->AddIRegion(alice);
//    lhc->AddIRegion(cms);
//    lhc->AddIRegion(atlas); 
//
//  The beams
//
    Float_t n         = 6.8e7;
    Float_t epsH      = 1.5e-4;
    Float_t epsL      = 2.5e-9;
    Float_t de        = 1.14e-4;
    Float_t energy    = 7000;
    Int_t a1          = 208;
    Int_t z1          =  82;
    Int_t a2          = 208;
    Int_t z2          =  82;
    
    
//
    AliLhcBeam* beam1 =  new AliLhcBeam(lhc);
    beam1->SetN(n);
    beam1->SetNEmittance(epsH);            // (cm)
    beam1->SetLongEmittance(epsL);         // (GeV s)
    beam1->SetEnergy(energy);                 // (GeV)
    beam1->SetParticle(a1,z1); 
    beam1->SetEnergySpread(de);

    AliLhcBeam* beam2 =  new AliLhcBeam(lhc);
    beam2->SetN(n);
    beam2->SetNEmittance(epsH);            // (cm)
    beam2->SetLongEmittance(epsL);         // (GeV s)
    beam2->SetEnergy(energy);                 // (GeV)
    beam2->SetParticle(a2,z2);
    beam2->SetEnergySpread(de);

    
    lhc->SetBeams(beam1, beam2);
//
//  The Processes
//    
    AliLhcProcessBB*  bb  = new AliLhcProcessBB(lhc, "BB", "Beam-Beam Losses");
    bb->SetCrossSection(505.);
    AliLhcProcessIBS* ibs = new AliLhcProcessIBS(lhc, "IBS", "Intra Beam Scattering");
    AliLhcProcessBT*  bt  = new AliLhcProcessBT(lhc, "BT", "Beta* Tuning");
    bt->SetBetaMin(50.);

    lhc->AddProcess(ibs);
    lhc->AddProcess(bb);
//    lhc->AddProcess(bt);
//
//  Run the collider
//
    lhc->SetTime(100., 20.*3600.);
    lhc->Init();
    lhc->EvolveTime();
    lhc->Evaluate();
}










