/**
 * @file   EGConfig.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 16 11:01:38 2014
 * 
 * @brief  Specific configuration for event generator. 
 * 
 * 
 */
struct EGCfg : public VirtualEGCfg
{
  Int_t   hftype;     // Heavy flavour type (random)

  EGCfg() 
    : hftype(-1)
  {
    hftype = HFType();
  }
protected:
  /** 
   * Make a random Heavy Flavour type 
   * 
   * @return Heavy flavour type 
   */
  Int_t HFType() const 
  {
    Int_t type = -1;
    Int_t r = gRandom->Rndm();
    if      (r < 0.16) type = 0; 
    else if (r < 0.32) type = 1; 
    else if (r < 0.48) type = 2; 
    else if (r < 0.64) type = 3; 
    else if (r < 0.72) type = 4; 
    else if (r < 0.80) type = 5; 
    else               type = 6; 
    return r;
  }
  /** 
   * Make the generator 
   * 
   * @param rt   Generator ID
   * @param b1   Least impact parameter 
   * @param b2   Largest impact parameter 
   *
   * @return Point to newly allocated generator or null 
   */
  AliGenerator* CreateGenerator(const TString& rt,
				Float_t b1, Float_t b2)
  {
    Bool_t        asym = grp->IsPA()||grp->IsAP();
    AliGenerator* g    = 0;
    if      (rt.EndsWith("perugia0chadr"))     g=PythiaHF(0);
    else if (rt.EndsWith("perugia0bchadr"))    g=PythiaHF(1);
    else if (rt.EndsWith("perugia0cele"))      g=PythiaHF(2);
    else if (rt.EndsWith("perugia0bele"))      g=PythiaHF(3);
    else if (rt.EndsWith("perugia0jspi2e"))    g=PythiaHF(4);
    else if (rt.EndsWith("perugia0btojspi2e")) g=PythiaHF(5);
    else if (rt.BeginsWith("pythia"))          g=Pythia(rt);
    else if (rt.BeginsWith("hijing2000hf"))    g=HFCocktail(rt, b1, b2);
    else if (rt.BeginsWith("hijing2000"))      g=Hijing(b1, b2, asym, 
							false, 2.3);
    else if (rt.BeginsWith("hijing"))          g=Hijing(b1, b2, asym, 
							grp->IsAA(), 0);
    else if (rt.BeginsWith("ampthf"))          g=HFCocktail(rt,b1,b2);
    else if (rt.BeginsWith("ampt"))            g=Ampt(b1, b2);
    else if (rt.BeginsWith("dpmjet"))          g=Dpmjet(b1, b2);
    else if (rt.BeginsWith("phojet"))          g=Dpmjet(b1, b2);
    else if (rt.BeginsWith("hydjet"))          g=Hydjet(b1, b2);

    if (!g)
      Fatal("", "Invalid run type \"%s\" specified", rt.Data());
    return g;
  }
  /** 
   * Make our decayer 
   * 
   * @return Newly allocated decayer or null
   */
  TVirtualMCDecayer* CreateDecayer(const TString& runType)
  {
    if (runType.BeginsWith("hydjet")) return 0;

    LoadPythia();
    TVirtualMCDecayer* decayer = new AliDecayerPythia();
    if (runType.EqualTo("hijing2000hf") && hftype < 2) 
      decayer->SetForceDecay(kHadronicD);
    else 
      decayer->SetForceDecay(kAll);
    decayer->Init();
    return decayer;
  }

  // === PYTHIA ========================================================
  // Normal 
  /** 
   * Greate a pythia6 event generator 
   * 
   * @param tune Possible tune 
   * 
   * @return newly allocated generator or null
   */
  AliGenerator* Pythia(const TString & tune)
  {
    // Int_t kCTEQ6l = 8;
    if (!grp->IsPP()) Fatal("Setup", "Pythia6 only works for pp");

    TString t(tune);
    t.ToUpper();
    t.ReplaceAll("PYTHIA6", "");
    t.ReplaceAll("PYTHIA", "");
    Info("Setup", "Making Pythia6 event generator (tune: %s)", t.Data());

    LoadPythia();
    AliGenPythia* pythia = new AliGenPythia(-1); 
    pythia->SetMomentumRange(0, 999999.);
    pythia->SetThetaRange(0., 180.);
    pythia->SetYRange(-12.,12.);
    pythia->SetPtRange(0,1000.);
    pythia->SetProcess(kPyMb);
    pythia->SetEnergyCMS(grp->energy);

    if (t == "D6T") {
      //    Tune
      //    109     D6T : Rick Field's CDF Tune D6T 
      //                  (NB: needs CTEQ6L pdfs externally)
      pythia->SetTune(109); // F I X 
      pythia->SetStrucFunc(kCTEQ6l);
    }
    else if (t == "PERUGIA0") { 
      //    Tune
      //    320     Perugia 0
      pythia->SetTune(320); 
      pythia->UseNewMultipleInteractionsScenario();
    }
    else if (t == "ATLAS") {
      //    Tune
      //    C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune 
      //                        (needs CTEQ6L externally)
      pythia->SetTune(306);
      pythia->SetStrucFunc(kCTEQ6l);
    }
    else if (t == "JETS") { 
      pythia->SetProcess(kPyJets);
      pythia->SetStrucFunc(kCTEQ6l);
      pythia->SetJetEtaRange(-1.5, 1.5); 
      pythia->SetJetEtRange(50., 800.);
      pythia->SetPtHard(45., 1000.);
      pythia->SetPycellParameters(2.2, 300, 432, 0., 4., 5., 0.7);
    }
    else if (t == "ATLAS_FLAT") {
      // set high multiplicity trigger
      // this weight achieves a flat multiplicity distribution
      TH1 *weight = new TH1D("weight","weight",201,-0.5,200.5);
      weight->SetBinContent(1,5.49443);
      weight->SetBinContent(2,8.770816);
      weight->SetBinContent(6,0.4568624);
      weight->SetBinContent(7,0.2919915);
      weight->SetBinContent(8,0.6674189);
      weight->SetBinContent(9,0.364737);
      weight->SetBinContent(10,0.8818444);
      weight->SetBinContent(11,0.531885);
      weight->SetBinContent(12,1.035197);
      weight->SetBinContent(13,0.9394057);
      weight->SetBinContent(14,0.9643193);
      weight->SetBinContent(15,0.94543);
      weight->SetBinContent(16,0.9426507);
      weight->SetBinContent(17,0.9423649);
      weight->SetBinContent(18,0.789456);
      weight->SetBinContent(19,1.149026);
      weight->SetBinContent(20,1.100491);
      weight->SetBinContent(21,0.6350525);
      weight->SetBinContent(22,1.351941);
      weight->SetBinContent(23,0.03233504);
      weight->SetBinContent(24,0.9574557);
      weight->SetBinContent(25,0.868133);
      weight->SetBinContent(26,1.030998);
      weight->SetBinContent(27,1.08897);
      weight->SetBinContent(28,1.251382);
      weight->SetBinContent(29,0.1391099);
      weight->SetBinContent(30,1.192876);
      weight->SetBinContent(31,0.448944);
      for (Int_t i = 32; i <= 201; i++) weight->SetBinContent(i,1);
      weight->SetEntries(526);
      
      Int_t limit = weight->GetRandom();
      pythia->SetTriggerChargedMultiplicity(limit, 1.4);
    }
    return pythia;
  }
  /** 
   * Create a Pythia6 generator for high-flavor physics 
   * 
   * @param type    Which kind 
   * @param harder  If true, make harder processes 
   * 
   * @return Newly allocated generator or null
   */
  AliGenerator* PythiaHF(Int_t type, Bool_t harder=0) 
  { 
    LoadPythia();
    if (type == 6) return Pythia("jets");
    if (type == 4) { 
      AliGenParam *jpsi =  AliGenParam(1, AliGenMUONlib::kJpsi,
				       (harder?"CDF pp 8.8":"CDF pp 7"),"Jpsi");
      jpsi->SetPtRange(0.,999.);
      jpsi->SetYRange(-1.0, 1.0);
      jpsi->SetPhiRange(0.,360.);
      jpsi->SetForceDecay(kDiElectron);
      return jpsi;
    }
    AliGenPythia* pythia = static_cast<AliGenPythia*>(Pythia("PERUGIA0"));
    switch (type) { 
    case 0: // chadr
      pythia->SetProcess(kPyCharmppMNRwmi);
      pythia->SetForceDecay(kHadronicD);
      break;
    case 1: // bchadr
      pythia->SetProcess(kPyBeautyppMNRwmi);
      pythia->SetForceDecay(kHadronicD);
      break;
    case 2: // cele
      pythia->SetProcess(kPyCharmppMNRwmi);
      pythia->SetCutOnChild(1);
      pythia->SetPdgCodeParticleforAcceptanceCut(11);
      pythia->SetChildYRange(-1.2,1.2);
      pythia->SetChildPtRange(0,10000.);
      break;
    case 3: // bele
      pythia->SetProcess(kPyBeautyppMNRwmi);
      pythia->SetCutOnChild(1);
      pythia->SetPdgCodeParticleforAcceptanceCut(11);
      pythia->SetChildYRange(-1.2,1.2);
      pythia->SetChildPtRange(0,10000.);
      break;
    case 5:
      pythia->SetProcess(kPyBeautyppMNRwmi);
      pythia->SetCutOnChild(1);
      pythia->SetPdgCodeParticleforAcceptanceCut(443);
      pythia->SetChildYRange(-2,2);
      pythia->SetChildPtRange(0,10000.);
    }
    return pythia;
  }
  /** 
   * Make a Min-Bias AA, pA, or Ap Hijing generator 
   * 
   * @param minB    Least impact parameter 
   * @param maxB    Largest impact parameter 
   * @param quench  If true, enable quenching 
   * @param slowN   If true, make a cocktail with slow neutrons 
   * @param ptHard  Hard pT cut-off 
   * 
   * @return Generator 
   */
  AliGenerator* Hijing(Float_t minB, 
		       Float_t maxB, 
		       Bool_t  slowN=false, 
		       Bool_t  quench=1, 
		       Float_t ptHard=0) 
  {
    LoadHijing();
    AliGenHijing *gener = new AliGenHijing(-1);
    // centre of mass energy 
    gener->SetEnergyCMS(grp->energy);
    gener->SetImpactParameterRange(minB, maxB);	
    // reference frame
    gener->SetReferenceFrame("CMS");
    // projectil
    gener->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    gener->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    // tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
    // enable jet quenching
    gener->SetJetQuenching(quench);
    // enable shadowing
    gener->SetShadowing(slowN);
    // Don't track spectators
    gener->SetSpectators(!slowN);
    // 
    if (ptHard > 0) hi->SetPtHardMin(ptHard);

    // kinematic selection
    gener->SetSelectAll(0);
    // Boosted CMS 
    gener->SetBoostLHC(grp->IsPA() || grp->IsAP());
    // No need for cocktail 
    if (!slowN || !grp->IsPA() || !grp->IsAP()) return gener;

    AliGenCocktail* cocktail = new AliGenCocktail();
    cocktail->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    cocktail->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    cocktail->SetEnergyCMS(grp->energy);

    AliGenSlowNucleons*     gray  = new AliGenSlowNucleons(1);
    AliCollisionGeometry*   coll  = gener->CollisionGeometry();
    AliSlowNucleonModelExp* model = new AliSlowNucleonModelExp();
    //  Not yet in the release...
    //      model->SetSaturation(kTRUE);
    gray->SetSlowNucleonModel(model);
    gray->SetTarget(grp->beam2.a, grp->beam2.z);
    gray->SetThetaDist(1);
    gray->SetProtonDirection(grp->beam1.IsP() ? 1 : 2);
    //      gray->SetDebug(1);
    gray->SetNominalCmsEnergy(2*grp->beamEnergy);
    gray->NeedsCollisionGeometry();
    gray->SetCollisionGeometry(coll);

    cocktail->AddGenerator(gener, "Hijing pPb", 1);
    cocktail->AddGenerator(gray, "Gray Particles", 1);
    
    return cocktail;
  }
  /** 
   * Make a DPMJet generator for pp, AA, pA, or Ap. 
   * 
   * @param fragments If true, make fragments 
   * 
   * @return Generator 
   */
  AliGenerator* Dpmjet(Float_t minB, Float_t maxB, 
		       Bool_t fragments=0)
  {
    LoadDpmjet();
    AliGenDPMjet* dpmjet = new AliGenDPMjet(-1); 
    dpmjet->SetEnergyCMS(grp->energy);
    dpmjet->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    dpmjet->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    dpmjet->SetImpactParameterRange(minB, maxB);
    dpmjet->SetProjectileBeamEnergy(grp->beam1.z*grp->beamEnergy/grp->beam1.a);
    if (grp->IsAA()) { 
      dpmjet->SetPi0Decay(0);
    }
    else if (grp->IsPA() || grp->IsAP()) { 
      dpmjet->SetTriggerParticle(3312, 1.2, 2.0);
      dpmjet->SetFragmentProd(fragments); // Alwas disabled 
    }
    else if (grp->IsPP()) { // PhoJet
      dpmjet->SetMomentumRange(0, 999999.);
      dpmjet->SetThetaRange(0., 180.);
      dpmjet->SetYRange(-12.,12.);
      dpmjet->SetPtRange(0,1000.);
    }
    return dpmjet;
  }
  /** 
   * Make an AMPT generator for AA collisions 
   * 
   * @return Generator 
   */
  AliGenerator* Ampt(Float_t minB, Float_t maxB)
  {
    LoadAmpt();
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    genHi->SetEnergyCMS(grp->energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    genHi->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    genHi->SetPtHardMin (2);
    genHi->SetImpactParameterRange(minB,maxB);
    // disable jet quenching
    genHi->SetJetQuenching(0); 
    // enable shadowing
    genHi->SetShadowing(1);    
    // neutral pion and heavy particle decays switched off
    genHi->SetDecaysOff(1);
    genHi->SetSpectators(0);   // track spectators 
    genHi->KeepFullEvent();
    genHi->SetSelectAll(0);
    return genHi;
  }
  /** 
   * Make an HydJet generator for A-A
   * 
   * @return Generator 
   */
  AliGenerator* Hydjet(Float_t minB, Float_t maxB)
  {
    LoadHydjet();
    AliGenUHKM *genHi = new AliGenUHKM(-1);
    genHi->SetAllParametersLHC();
    genHi->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    genHi->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    genHi->SetEcms(grp->energy);
    genHi->SetEnergyCMS(grp->energy);
    genHi->SetBmin(minB);
    genHi->SetBmax(maxB);
    genHi->SetPyquenPtmin(9);
    return genHi;
  }
  /** 
   * Make a heavy flavour cocktail 
   * 
   * @param base Underlying event. 
   * 
   * @return Generator 
   */
  AliGeneator* HFCocktail(const TString& base, Float_t minB, Float_t maxB) 
  {
    
    AliGenCocktail *cocktail = new AliGenCocktail();
    cocktail->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    cocktail->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    cocktail->SetEnergyCMS(grp->energy);
    
    // Add underlying event
    if (base.BeginsWith("ampt", TString::kIgnoreCase)) { 
      hi = Ampt(minB, maxB);
      cocktail->AddGenerator(hi,"ampt",1);
    }
    else { 
      hi = Hijing(minB, maxB, grp->IsPA() || grp->IsAP(), false, 2.3);
      cocktail->AddGenerator(hi,"hijing",1);

    }

    // --- Default formula -------------------------------------------
    TForumla* one = new TFormula("one", "1.");

    // --- Pythia ----------------------------------------------------
    AliGenerator* pythia = PythiaHF(hftype);
    switch (hftype) { 
    case 6: 
      cocktail->AddGenerator(pythia, "pythiaJets", 1, one);
      break;
    defualt:
      cocktail
	->AddGenerator(pythia, "pythiaHF", 1, 
		       new TFormula("Signals", 
				    "20.*(x<5.)+80./3.*(1.-x/20.)*(x>5.)"));
      break;
    }
    // --- Dummy -----------------------------------------------------
    AliGenParam* param = 0;

    // --- Phos stuff ------------------------------------------------
    AliGenPHOSlib *plib = new AliGenPHOSlib();
    Double_t lower[] = { 0, 3, 6, 9, 12, -1 };
    Double_t *pLow   = lower;
    for (Int_t i = 0; i < 5; i++) {
      param = new AliGenParam(5, plib, AliGenPHOSlib::kPi0);
      param->SetPhiRange(0., 360.) ;
      param->SetYRange(-1.2, 1.2) ;
      param->SetPtRange(lower[i], 30.) ;
      cocktail->AddGenerator(param,Form("Pi0HagPt%d", i), 1., one);
    }

    // --- Jpsi->mu+ mu-----------------------------------------------
    param = new AliGenParam(1, AliGenMUONlib::kJpsi, "CDF pp 3.94", "Jpsi");
    param->SetPtRange(0.,999.);
    param->SetYRange(-4.2, -2.3);
    param->SetPhiRange(0.,360.);
    param->SetForceDecay(kDiMuon);
    param->SetCutOnChild(1);
    param->SetChildPhiRange(0.,360.);
    param->SetChildThetaRange(168.5,178.5);
    cocktail->AddGenerator(param, "Jpsi2M", 1, one); 

    // --- Chi_c -> J/Psi + gamma, J/Psi -> e+e- ---------------------
    Float_t thmin  = (180./TMath::Pi())*2.*atan(exp(-1.2));  
    Float_t thmax  = (180./TMath::Pi())*2.*atan(exp( 1.2));  
    param  = new AliGenParam(1, AliGenMUONlib::kChic,"default","Chic");
    param->SetMomentumRange(0, 999.);        // Wide cut on the momentum
    param->SetPtRange(0, 100.);              // Wide cut on Pt
    param->SetYRange(-2.5, 2.5);
    param->SetCutOnChild(1);                 // Enable cuts on decay products
    param->SetChildPhiRange(0., 360.);
    // In the acceptance of the Central Barrel
    param->SetChildThetaRange(thmin, thmax); 
    // Chi_c -> J/Psi + gamma, J/Psi -> e+e-
    param->SetForceDecay(kChiToJpsiGammaToElectronElectron); 
    cocktail->AddGenerator(param, "Chi_c", 1, one); 

    // --- Dummy -----------------------------------------------------
    AliGenBox* box = 0;

    // --- Some particles --------------------------------------------
    Double_t    boxR   = gRandom->Integer(3)+1;
    Int_t       boxP[] = { 310, 3122,   3312, 3322, 
			   (boxR==1 ?3334: boxR==2 ?-3334: -3312) };
    Int_t       boxN[] = {   1,    1,      3,    3, 1 }
    const char* boxT[] = { "K0s", "Lambda", "Xi-", "Xi0", 
			   (boxR==1 ? "Omega-": boxR==2 ? "Omega+": "Xi+") };
    for (Int_t i = 0; i < 5; i++) {
      box = new AliGenBox(boxN[i]);
      box->SetPart(boxP[i]);
      box->SetPtRange(0,13);
      box->SetThetaRange(45, 135);
      cocktail->AddGenerator(box, boxT[i], 1, one);
    }
    
    // --- High pT charged particle ----------------------------------
    TFormula* hptF = new TFormula("High Pt", 
				  "5.*(x<5.)+20./3.*(1.-x/20.)*(x > 5.)");
    Int_t       hptP[] = { 211, 321, 2212 };
    const char* hptT[] = { "pi", "K", "p" };
    for (Int_t i = 0; i < 3; i++) { 
      for (Int_t j = -1; j <= 1; j++) {
	box->SetPart(j*hptP[i]);
	box->SetPtRange(2., 50.);
	box->SetThetaRange(thmin, thmax);
	cocktail->AddGenerator(box, Form("%s%c",hptT[i],(j<0,'-','+')),1,hptF);
      }
    }
    return cocktail;
  }

};

void 
EGConfig()
{
  ::Info("EGConfig", "Creating EG factory");
  egCfg = new EGCfg;
}

// 
// EOF
// 


