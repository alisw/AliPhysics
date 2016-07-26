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
  Bool_t  fIsLego;    //
  
  EGCfg() 
    : hftype(-1), fIsLego(false);
  {
    hftype = HFType();
  }
  virtual Bool_t IsLego() const { return fIsLego; }
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
    TString t(rt);
    if      (t.EqualTo("default"))            t = DeduceRunType();
    if      (t.EndsWith("perugia0chadr"))     g=PythiaHF(0);
    else if (t.EndsWith("perugia0bchadr"))    g=PythiaHF(1);
    else if (t.EndsWith("perugia0cele"))      g=PythiaHF(2);
    else if (t.EndsWith("perugia0bele"))      g=PythiaHF(3);
    else if (t.EndsWith("perugia0jspi2e"))    g=PythiaHF(4);
    else if (t.EndsWith("perugia0btojspi2e")) g=PythiaHF(5);
    else if (t.BeginsWith("pythia"))          g=Pythia(rt);
    else if (t.BeginsWith("hijing2000hf"))    g=HFCocktail(rt, b1, b2);
    else if (t.BeginsWith("hijing2000"))      g=Hijing(b1, b2, rt);
    else if (t.BeginsWith("hijing"))          g=Hijing(b1, b2, rt);
    else if (t.BeginsWith("ampthf"))          g=HFCocktail(rt,b1,b2);
    else if (t.BeginsWith("ampt"))            g=Ampt(b1, b2, rt);
    else if (t.BeginsWith("dpmjet"))          g=Dpmjet(b1, b2);
    else if (t.BeginsWith("phojet"))          g=Dpmjet(b1, b2);
    else if (t.BeginsWith("hydjet"))          g=Hydjet(b1, b2);
    else if (t.BeginsWith("epos-lhc"))        g=EposLHC(b1, b2);
    else if (t.BeginsWith("epos"))            g=Epos(b1, b2);
    else if (t.BeginsWith("therminator"))     g=Therminator(b1, b2);
    else if (t.BeginsWith("lego"))            g=Lego(rt);
    if (!g && !fIsLego)
      Fatal("", "Invalid run type \"%s\" specified", t.Data());
    AliPDG::AddParticlesToPdgDataBase();
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

    TString tit(Form("Pythia6 %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		     grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		     grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		     Int_t(grp->energy)));

    TString t(tune);
    t.ToUpper();
    t.ReplaceAll("PYTHIA6", "");
    t.ReplaceAll("PYTHIA", "");
    Info("Setup", "Making Pythia6 event generator (tune: %s)", t.Data());

    tit.Append(Form(" tune=%s", t.Data()));
    
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
    pythia->SetTitle(tit);
    return pythia;
  }
  /** 
   * Create a Pythia6 generator for heavy-flavor physics 
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
		       const TString& rt) 
  {
    LoadHijing();
    TString opt(rt);
    opt.ToLower();
    opt.Remove(0,6); // Remove hijing prefix 

    // When no options are passed, we
    //  Enable shadowing
    //  Disable spectators
    //  Disable quenching
    // Which corresponds to the defaults for PbPb MB.  Note, if we
    // simulate pA/Ap then, we always
    //  Disable shadowing
    //  Enable spectators
    // and add a slow-nucleon model afterburner 
    Bool_t   quench = opt.Contains("quench");
    Bool_t   spec   = ((grp->IsPA()||grp->IsAP() || opt.Contains("spectators"))
		       && !opt.Contains("nospectators"));
    Bool_t   slow   = (grp->IsPA() || grp->IsAP()) && !opt.Contains("noslow");
    Bool_t   shadow = !spec && !slow && !opt.Contains("noshadow");
    UInt_t   sNN    = (grp->energy/10)*10;
    Double_t ptCut  = sNN < 5000 ? 2.3 : 2.8;
    TString  tit(Form("Hijing%s %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		      (slow ? "+SNM" : ""), 
		      grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		      grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		      Int_t(grp->energy), minB, maxB));
    tit.Append(Form(" %squench",    quench ? "" : "no"));
    tit.Append(Form(" %sspectator", spec   ? "" : "no"));
    tit.Append(Form(" %sshadow",    shadow ? "" : "no"));
    tit.Append(Form(" ptCut=%5.2f", ptCut));
    
    AliGenHijing *gener = new AliGenHijing(-1);
    // --- Centre of mass energy -------------------------------------
    gener->SetEnergyCMS(sNN);
    // --- Impact parameter range ------------------------------------
    gener->SetImpactParameterRange(minB, maxB);	
    // --- Reference frame -------------------------------------------
    gener->SetReferenceFrame("CMS");
    // --- projectile ------------------------------------------------
    gener->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    // -- Target -----------------------------------------------------
    gener->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    // --- tell hijing to keep the full parent child chain - default 0
    gener->KeepFullEvent();
    // --- enable jet quenching - off for PbPb - 1 by default --------
    gener->SetJetQuenching(quench);
    // --- enable shadowing - on for PbPb - 1 by default -------------
    gener->SetShadowing(shadow); 
    // --- Don't track spectators - off for PbPb - 1 by default ------
    gener->SetSpectators(spec);
    // --- Possibly Pt cut-off - 2.3 for PbPb ------------------------
    gener->SetPtHardMin(ptCut);
    // --- Do not disable decays -- 3 for PbPb -----------------------
    // gener->SetDecaysOff(3); 
    // --- kinematic selection - 0 by default ------------------------
    // gener->SetSelectAll(0);
    // Boosted CMS 
    gener->SetBoostLHC(grp->IsPA() || grp->IsAP());

    // Debug
    // gener->GetTHijing()->SetIHPR2(10,1);
    
    // No need for cocktail
    if (!grp->IsPA() && !grp->IsAP()) {
      gener->SetTitle(tit);
      return gener;
    }


    AliGenCocktail* cocktail = new AliGenCocktail();
    cocktail->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    cocktail->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    cocktail->SetEnergyCMS(grp->energy);
    cocktail->SetName("HIJINGSNM");
    
    AliGenSlowNucleons*     gray  = new AliGenSlowNucleons(1);
    AliCollisionGeometry*   coll  = gener->CollisionGeometry();
    AliSlowNucleonModelExp* model = new AliSlowNucleonModelExp();
    //  Not yet in the release...
    //      model->SetSaturation(kTRUE);
    gray->SetSlowNucleonModel(model);
    gray->SetTarget(grp->beam1.a, grp->beam1.z);
    gray->SetThetaDist(1);
    gray->SetProtonDirection(grp->beam1.IsP() ? 1 : 2);
    //      gray->SetDebug(1);
    gray->SetNominalCmsEnergy(2*grp->beamEnergy);
    gray->NeedsCollisionGeometry();
    gray->SetCollisionGeometry(coll);

    cocktail->AddGenerator(gener, "Hijing pPb", 1);
    cocktail->AddGenerator(gray, "Gray Particles", 1);

    tit.Append(" + slow nucleon");
    cocktail->SetTitle(tit);
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
    TString tit(Form("DpmJet %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		     grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		     grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		     Int_t(grp->energy), minB, maxB));

    dpmjet->SetEnergyCMS(grp->energy);
    dpmjet->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    dpmjet->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    dpmjet->SetImpactParameterRange(minB, maxB);
    dpmjet->SetProjectileBeamEnergy(grp->beam2.z*grp->beamEnergy/grp->beam2.a);
    if (grp->IsPA() || grp->IsAP()) { 
      // dpmjet->SetTriggerParticle(3312, 1.2, 2.0);
      dpmjet->SetFragmentProd(false/*fragments*/); // Alwas disabled 
      dpmjet->SetSpectators(true); // Spectators
    }
    else if (grp->IsPP()) { // PhoJet
      dpmjet->SetMomentumRange(0, 999999.);
      dpmjet->SetThetaRange(0., 180.);
      dpmjet->SetYRange(-12.,12.);
      dpmjet->SetPtRange(0,1000.);
      tit.Append(" (Phojet)");
    }
    dpmjet->SetTitle(tit);
    return dpmjet;
  }
  /** 
   * Make an AMPT generator for AA collisions 
   * 
   * @return Generator 
   */
  AliGenerator* Ampt(Float_t minB, Float_t maxB, const TString& rt)
  {
    LoadAmpt();
    TString opt(rt);
    opt.ToLower();
    opt.Remove(0,4); // Remove AMPT prefix 

    TString tit(Form("AMPT %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		     grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		     grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		     Int_t(grp->energy), minB, maxB));

    // When no options are passed, we
    //  Enable decayer
    //  Turn on screening mass
    //  Dislable string melting
    //  Enable shadowing
    //  Disable spectators
    //  Disable quenching
    // Which corresponds to the defaults for PbPb MB 
    Bool_t decayer = !opt.Contains("nodecay"); // Default for PbPb 
    Bool_t screen  = !opt.Contains("noscreen");
    Bool_t melt    = opt.Contains("melt");
    Bool_t shadow  = !opt.Contains("noshadow");
    Bool_t spec    = opt.Contains("spectators");
    Bool_t quench  = opt.Contains("quench");
    tit.Append(Form(" %sdecay",     decayer ? "" : "no"));
    tit.Append(Form(" %sscreen",    screen  ? "" : "no"));
    tit.Append(Form(" %smelt",      melt    ? "" : "no"));
    tit.Append(Form(" %sshadow",    shadow  ? "" : "no"));
    tit.Append(Form(" %sspectator", spec    ? "" : "no"));
    tit.Append(Form(" %squench",    quench  ? "" : "no"));
    
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    genHi->SetEnergyCMS(grp->energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    genHi->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    genHi->SetImpactParameterRange(minB,maxB);
    // --- Least hard Pt ---------------------------------------------
    genHi->SetPtHardMin (3);
    // --- disable jet quenching -------------------------------------
    genHi->SetJetQuenching(quench); 
    // --- enable shadowing ------------------------------------------
    genHi->SetShadowing(shadow);    
    // neutral pion and heavy particle decays switched off
    genHi->SetDecaysOff(decayer);
    // --- track spectators ------------------------------------------
    genHi->SetSpectators(spec);
    // --- Keep everything -------------------------------------------
    genHi->KeepFullEvent();
    // --- Do not use all --------------------------------------------
    genHi->SetSelectAll(0);

    // -- String melting: 1 default ----------------------------------
    genHi->SetIsoft(melt ? 4 : 1);
    // --- Lund string fragmentation parameters ----------------------
    genHi->SetStringFrag(0.5, 0.9);
    // --- MAx number of time steps ----------------------------------
    genHi->SetNtMax(150);
    // --- Boost according to LHC parameters -------------------------
    genHi->SetBoostLHC(1);
    // --- Create random reaction plane ------------------------------
    genHi->SetRandomReactionPlane(true);
    // --- parton screening mass in fm^(-1) (D=3.2264d0) -------------
    if (screen) genHi->SetXmu(3.2264);
    
    // -- Process tunes ----------------------------------------------
    if (decayer) {
      // Add a decayer here
      hftype = 1;
      genHi->SetDecayer(CreateDecayer("hijing2000hf"));
    }
    genHi->SetTitle(tit);
    
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
    genHi->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    genHi->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    genHi->SetEcms(grp->energy);
    genHi->SetEnergyCMS(grp->energy);
    genHi->SetBmin(minB);
    genHi->SetBmax(maxB);
    genHi->SetPyquenPtmin(9);

    TString tit(Form("Hydjet %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		     grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		     grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		     Int_t(grp->energy), minB, maxB));
    genHi->SetTitle(tit);
    return genHi;
  }
  /** 
   * Make an Epos=LHC generator for p-p, A-A, p-A, or A-p
   * 
   * @return Generator 
   */
  AliGenerator* EposLHC(Float_t minB, Float_t maxB)
  {
    LoadEposLHC();
    AliGenEposLHC* gen = new AliGenEposLHC(-1);
    gen->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    gen->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    // We should set the beam momenta instead of the CM energy 
    //  gen->SetEnergyCMS(grp->energy);
    // Note, one of the beams should have negative momentum 
    gen->SetPTarget    (-grp->BeamMomentum(1));
    gen->SetPProjectile(+grp->BeamMomentum(2));
    gen->SetImpactParameterRange(minB, maxB);
    TString tit(Form("EPOS-LHC %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		     grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		     grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		     Int_t(grp->energy), minB, maxB));
    gen->SetTitle(tit);
    return gen;
  }    
  /** 
   * Make an Epos generator for p-p (A-A, or p-A possible?)
   * 
   * @return Generator 
   */
  AliGenerator* Epos(Float_t,Float_t)
  {
    LoadEpos();
    AliGenEpos* gen = new AliGenEpos();
    gen->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    gen->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    gen->SetEnergyCMS(grp->energy);
    TString tit(Form("EPOS %s(%d,%d)+%s(%d,%d) @ %5d b in[%4.1f,%4.1f]",
		     grp->beam1.Name(), grp->beam1.a, grp->beam1.z, 
		     grp->beam2.Name(), grp->beam2.a, grp->beam2.z,
		     Int_t(grp->energy), minB, maxB));
    gen->SetTitle(tit);
    return gen;
  }    
  AliGenerator* Therminator(Float_t,Float_t)
  {
    LoadTherminator();
    AliGenTherminator* gen = new AliGenTherminator();
    gen->SetFileName("event.out");
    gen->SetEventNumberInFile(1);
    gen->SetTemperature(.145);
    gen->SetMiuI(-0.0009);
    gen->SetMiuS(0.000);
    gen->SetMiuB(0.0008);
    gen->SetAlfaRange(8.0);
    gen->SetRapRange(4.0);
    gen->SetRhoMax(7.74);
    gen->SetTau(9.74);
    gen->SetModel("Lhyquid3D");
    gen->SetLhyquidSet("LHC500C2030");
    gen->SetTitle("Therminator");
    return gen;

  }
  // === Lego ========================================================
  /** 
   * Greate a lego event generator 
   * 
   * @param tune Possible tune 
   * 
   * @return newly allocated generator or null
   */
  AliGenerator* Lego(const TString & variant)
  {
    fIsLego = true;
    return 0;
#if 0
    TString v(variant);
    v.ToUpper();
    v.ReplaceAll("LEGO", "");
    Info("Setup", "Making Lego event generator (variant: %s)", v.Data());

    AliLegoGenerator* ret = 0;
    // XYZ varies origin of the particles in two dimensions:
    //  X:  o=(0,t1,t2), p=(1,0,0)
    //  Y:  o=(t1,0,t2), p=(0,1,0)
    //  Z:  o=(t1,t2,0), p=(0,0,1)
    // PhiZ varies the momentum in two dimensions
    //  o=(0,0,t1) p=(cos(t2),sin(t2),0)
    // Eta varies momentum in two dimensions
    //  phi=t1
    //  theta=2*atan(exp(-t2))
    //  o=(0,0,0) p=(cos(phi)*sin(theta),sin(phi)*cos(theta),cos(theta))
    // Base varies in two dimensions
    //  phi=t1
    //  theta=t2
    //  o=(0,0,0) p=(cos(phi)*sin(theta),sin(phi)*cos(theta),cos(theta))    
    if (v.BeginsWith("X") || v.BeginsWith) {
      const char* c[] = { v(0), '\0' };
      ret = new AliLegoGeneratorXYZ(c);
      ret->SetCoor1Range(10,-2,2);   // Y, X
      ret->SetCoor2Range(10,-10,10); // Z
    }
    else if (v.BeginsWith("Z")) {
      ret = new AliLegoGeneratorXYZ("Z");
      ret->SetCoor1Range(10,-2,2);   // X
      ret->SetCoor2Range(10,-2,2);   // Y
    }
    else if (v.BeginsWith("PHIZ")) {
      ret = new AliLegoGeneratorPhiZ();
      ret->SetCoor1Range(10,-10,10); // Z
      ret->SetCoor2Range(360,0,360); // phi
    }
    else if (v.BeginsWith("ETA")) {
      ret = new AliLegoGeneratorEta();
      ret->SetCoor1Range(360,0,360); // phi
      Double_t aEta = 7;
      Double_t dEta = (6--4)/200;
      ret->SetCoor2Range(2*aEta/dEta,-aEta,+aEta); // Eta
    }
    else {
      ret = new AliLegoGenerator();
      ret->SetCoor1Range(180,0,180); // theta
      ret->SetCoor1Range(360,0,360); // phi
    }
    return ret;
#endif
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
    cocktail->SetTarget    (grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    cocktail->SetProjectile(grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    cocktail->SetEnergyCMS(grp->energy);
    
    // Add underlying event
    if (base.BeginsWith("ampt", TString::kIgnoreCase)) { 
      hi = Ampt(minB, maxB, base);
      cocktail->AddGenerator(hi,"ampt",1);
    }
    else { 
      hi = Hijing(minB, maxB, base);
      cocktail->AddGenerator(hi,"hijing",1);      
    }
    TString tit(hi->GetTitle());
    
    // --- Default formula -------------------------------------------
    TForumla* one = new TFormula("one", "1.");

    // --- Pythia ----------------------------------------------------
    AliGenerator* pythia = PythiaHF(hftype);
    tit.Append(Form(" + %s", pythia->GetTitle()));
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
    cocktail->SetTitle(tit);
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


