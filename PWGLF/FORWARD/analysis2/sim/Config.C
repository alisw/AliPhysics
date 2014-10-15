// -------------------------------------------------------------------
struct Setup
{
  TString runType;    // Event generator chosen
  UInt_t  seed;       // Random number seed - (env)
  Float_t minB;       // Least imp.param. (env)
  Float_t maxB;       // Largest imp.param. (env)  
  Int_t   hftype;     // Heavy flavour type (random)
  /** 
   * Get the value of an environment variable as a float 
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As float, or default
   */
  static Float_t Env2Float(const char* envName, Float_t def) 
  { 
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return val.Atof();
  }
  /** 
   * Get the value of an environment variable as a unsigned int 
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As unsigned int, or default
   */
  static UInt_t Env2UInt(const char* envName, UInt_t def)
  {
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return UInt_t(val.Atoll());
  }
  /** 
   * Get the value of an environment variable as a int
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As int, or default
   */
  static UInt_t Env2Int(const char* envName, Int_t def)
  {
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return val.Atoi();
  }

  /** 
   * Constructor - retrieves needed information from CDB manager and
   * environment.
   */
  Setup() 
    : runType(""),
      seed(0),
      minB(0), 
      maxB(100),
      hftype(-1)
  {
    TDatime now;
    runType = gSystem->Getenv("CONFIG_RUN_TYPE");
    seed    = Env2UInt("CONFIG_SEED", now.Get());
    minB    = Env2Float("CONFIG_BMIN", 0);
    maxB    = Env2Float("CONFIG_BMAX", 100);
    if (runType[0] == 'k') runType.Remove(0,1);
    runType.ToLower();

    // gROOT->Macro("GetGRP.C");

    Int_t pytR = gRandom->Rndm();
    if      (pytR < 0.16) hftype = 0; 
    else if (pytR < 0.32) hftype = 1; 
    else if (pytR < 0.48) hftype = 2; 
    else if (pytR < 0.64) hftype = 3; 
    else if (pytR < 0.72) hftype = 4; 
    else if (pytR < 0.80) hftype = 5; 
    else                  hftype = 6; 

    if (runType.IsNull() || runType == "default") DeduceRunType();

    Print();
  }
  void Print()
  {
    Printf("=======================================================\n"
	   "           Set-up of the simulation\n");
    grp->Print();
    Printf("Run type: '%s'", runType.Data());
    Printf("Seed:     %d", seed);
    Printf("\n"
	   "=======================================================");
  }
  /** 
   * Set the default generator based on the beam type 
   *
   * - p-p PYTHIA
   * - p-A or A-p DPMJet
   * - A-A Hijing 
   */
  void DeduceRunType()
  {
    if      (grp->IsPP())                runType = "pythia";
    else if (grp->IsPA() || grp->IsAP()) runType = "dpmjet";
    else if (grp->IsAA())                runType = "hijing";
  }
  void LoadGen() {
    if (!gROOT->GetClass("AliStructFuncType")) 
      gSystem->Load("liblhapdf");      // Parton density functions
    if (!gROOT->GetClass("TPythia6"))
      gSystem->Load("libEGPythia6")   // TGenerator interface
    if (!runType.EqualTo("hydjet")) 
      LoadPythia(false);
  }
  /** 
   * Load the pythia libraries 
   * 
   * @param vers Optional version post-fix
   */
  void LoadPythia(Bool_t gen=true, const char* vers="6.4.21")
  {
    if (gen) LoadGen();
    char m = vers[0];
    if (gROOT->GetClass(Form("AliPythia6%c", m))) return;
    gSystem->Load(Form("libpythia%s", vers));
    gSystem->Load(Form("libAliPythia%c", m));
  }
  /** 
   * Load HIJING libraries 
   */
  void LoadHijing() 
  {
    LoadPythia();
    if (gROOT->GetClass("THijing")) return;
    gSystem->Load("libhijing");
    gSystem->Load("libTHijing");
    AliPDG::AddParticlesToPdgDataBase();
  }
  /** 
   * Load HydJet libraries 
   */
  void LoadHydjet()
  {
    if (gROOT->GetClass("TUHKMgen")) return;
    gSystem->Load("libTUHKMgen");
  }
  /** 
   * Load DPMJet libraries
   */
  void LoadDpmjet()
  {
    LoadPythia();
    if (gROOT->GetClass("TDPMjet")) return;
    gSystem->Load("libdpmjet");
    gSystem->Load("libTDPMjet");
  }
  /** 
   * Load AMPT libraries
   */
  void LoadAmpt() 
  {
    LoadPythia();
    if (gROOT->GetClass("TAmpt")) return;
    gSystem->Load("libampt");
    gSystem->Load("libTAmpt");
  }
  /** 
   * Make the generator 
   * 
   * @return Point to newly allocated generator or null 
   */
  AliGenerator* MakeGenerator()
  {
    Bool_t asym = grp->IsPA()||grp->IsAP();
    TString& rt = runType;
    if (rt.EndsWith("perugia0chadr"))     return PythiaHF(0);
    if (rt.EndsWith("perugia0bchadr"))    return PythiaHF(1);
    if (rt.EndsWith("perugia0cele"))      return PythiaHF(2);
    if (rt.EndsWith("perugia0bele"))      return PythiaHF(3);
    if (rt.EndsWith("perugia0jspi2e"))    return PythiaHF(4);
    if (rt.EndsWith("perugia0btojspi2e")) return PythiaHF(5);
    if (rt.BeginsWith("pythia"))          return Pythia(rt);
    if (rt.BeginsWith("hijing2000hf"))    return HFCocktail(rt);
    if (rt.BeginsWith("hijing2000"))      return Hijing(asym, 
							false, 2.3);
    if (rt.BeginsWith("hijing"))          return Hijing(asym, 
							grp->IsAA(), 0);
    if (rt.BeginsWith("ampthf"))          return HFCocktail(rt);
    if (rt.BeginsWith("ampt"))            return Ampt();
    if (rt.BeginsWith("dpmjet"))          return Dpmjet();
    if (rt.BeginsWith("phojet"))          return Dpmjet();
    if (rt.BeginsWith("hydjet"))          return Hydjet();

    Fatal("", "Invalid run type \"%s\" specified", runType.Data());
    return 0;
  }
  TVirtualMCDecayer* MakeDecayer()
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
   * @param slowN If true, make a cocktail with slow neutrons 
   * 
   * @return Generator 
   */
  AliGenerator* Hijing(Bool_t slowN=false, Bool_t quench=1, Float_t ptHard=0) 
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
   * Make a DPMJet generator for AA, pA, or Ap. 
   * 
   * @param fragments If true, make fragments 
   * 
   * @return Generator 
   */
  AliGenerator* Dpmjet(Bool_t fragments=0)
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
  AliGenerator* Ampt()
  {
    LoadAmpt();
    AliGenAmpt *genHi = new AliGenAmpt(-1);
    genHi->SetEnergyCMS(grp->energy);
    genHi->SetReferenceFrame("CMS");
    genHi->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    genHi->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    genHi->SetPtHardMin (2);
    genHi->SetImpactParameterRange(bMin,bMax);
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
  AliGenerator* Hydjet()
  {
    LoadHydjet();
    AliGenUHKM *genHi = new AliGenUHKM(-1);
    genHi->SetAllParametersLHC();
    genHi->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    genHi->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    genHi->SetEcms(grp->energy);
    genHi->SetEnergyCMS(grp->energy);
    genHi->SetBmin(bMin);
    genHi->SetBmax(bMax);
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
  AliGeneator* HFCocktail(const TString& base) 
  {
    
    AliGenCocktail *cocktail = new AliGenCocktail();
    cocktail->SetProjectile(grp->beam1.Name(), grp->beam1.a, grp->beam1.z);
    cocktail->SetTarget    (grp->beam2.Name(), grp->beam2.a, grp->beam2.z);
    cocktail->SetEnergyCMS(grp->energy);
    
    // Add underlying event
    if (base.BeginsWith("ampt", TString::kIgnoreCase)) { 
      hi = Ampt();
      cocktail->AddGenerator(hi,"ampt",1);
    }
    else { 
      hi = Hijing(grp->IsPA() || grp->IsAP(), false, 2.3);
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





void Config()
{
  // --- Get settings from environment variables --------------------
  Setup    s;

  // ---- Seed random number generator -------------------------------
  gRandom->SetSeed(s.seed);
  std::cerr << "Seed for random number generation= " << s.seed << std::endl; 

  //------------------------------------------------------------------
  // 
  // Geometry and tracking 
  //
  // --- Libraries required by geant321 ------------------------------
  s.LoadGen();
  gSystem->Load("libgeant321");
  new TGeant3TGeo("C++ Interface to Geant3");

  // -----------------------------------------------------------------
  //  Create the output file
  std::cout<< "Config.C: Creating Run Loader ..." << std::endl;
  AliRunLoader* rl = AliRunLoader::Open("galice.root",
					AliConfig::GetDefaultEventFolderName(),
					"recreate");
  if (!rl) Fatal("Config","Can not instatiate the Run Loader");

  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);

  // --- Trigger configuration ---------------------------------------
  // AliSimulation::Instance()->SetTriggerConfig(grp->IsAA() ? "Pb-Pb" : "p-p");

  //
  //=======================================================================
  // Steering parameters for ALICE simulation
  // 
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  // --- Process switches --------------------------------------------
  gMC->SetProcess("DCAY",1);
  gMC->SetProcess("PAIR",1);
  gMC->SetProcess("COMP",1);
  gMC->SetProcess("PHOT",1);
  gMC->SetProcess("PFIS",0);
  gMC->SetProcess("DRAY",0);
  gMC->SetProcess("ANNI",1);
  gMC->SetProcess("BREM",1);
  gMC->SetProcess("MUNU",1);
  gMC->SetProcess("CKOV",1);
  gMC->SetProcess("HADR",1);
  gMC->SetProcess("LOSS",2);
  gMC->SetProcess("MULS",1);
  gMC->SetProcess("RAYL",1);
  
  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;

  // --- Tracking cuts -----------------------------------------------
  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE",  cut); 
  gMC->SetCut("BCUTM",  cut); 
  gMC->SetCut("DCUTE",  cut); 
  gMC->SetCut("DCUTM",  cut); 
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax); 
  
  // --- Set External decayer ----------------------------------------
  TVirtualMCDecayer* decayer = s.MakeDecayer();
  if (decayer) gMC->SetExternalDecayer(decayer);  

  //------------------------------------------------------------------
  // 
  // Generator Configuration 
  //
  // --- Make the generator - this loads libraries 
  AliGenerator* gener = s.MakeGenerator();
  gener->Init();

  // --- Go back to galice.root --------------------------------------
  rl->CdGAFile();
  
  // --- Switch on and off detectors ---------------------------------
  Int_t iABSO  = 1;
  Int_t iACORDE= 0;
  Int_t iDIPO  = 1;
  Int_t iEMCAL = 1;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
  Int_t iMAG   = 1;
  Int_t iMUON  = 1;
  Int_t iPHOS  = 1;
  Int_t iPIPE  = 1;
  Int_t iPMD   = 1;
  Int_t iHMPID = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 1;
  Int_t iTOF   = 1;
  Int_t iTPC   = 1;
  Int_t iTRD   = 1;
  Int_t iVZERO = 1;
  Int_t iZDC   = 1;
  

  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  
  
  if (iMAG)    new AliMAG("MAG", "Magnet");
  if (iABSO)   new AliABSOv3("ABSO", "Muon Absorber");
  if (iDIPO)   new AliDIPOv3("DIPO", "Dipole version 3");
  if (iHALL)   new AliHALLv3("HALL", "Alice Hall");
  if (iFRAME)  (new AliFRAMEv2("FRAME", "Space Frame"))->SetHoles(1);
  if (iSHIL)   new AliSHILv3("SHIL", "Shielding Version 3");
  if (iPIPE)   new AliPIPEv3("PIPE", "Beam Pipe");
  if (iITS)    new AliITSv11("ITS","ITS v11");
  // if (iITS)   new AliITSv11Hybrid("ITS","ITS v11Hybrid");
  if (iTPC)    new AliTPCv2("TPC", "Default");
  if (iTOF)    new AliTOFv6T0("TOF", "normal TOF");
  if (iHMPID)  new AliHMPIDv3("HMPID", "normal HMPID");
  if (iZDC) {
    AliZDC *ZDC = 0;
    if (grp->period.EqualTo("LHC10h")) {
      // Need to use older ZDC for PbPb 
      ZDC = new AliZDCv3("ZDC", "normal ZDC");
      ZDC->SetSpectatorsTrack();
    }
    else 
      ZDC = new AliZDCv4("ZDC", "normal ZDC");
    if (grp->Year() < 2011) { //?
      // What are these? Do they need to be set properly? 
      //Collimators aperture
      ZDC->SetVCollSideCAperture(0.85);
      ZDC->SetVCollSideCCentre(0.);
      ZDC->SetVCollSideAAperture(0.75);
      ZDC->SetVCollSideACentre(0.);
      //Detector position
      ZDC->SetYZNC(1.6);
      ZDC->SetYZNA(1.6);
      ZDC->SetYZPC(1.6);
      ZDC->SetYZPA(1.6);
    }
    ZDC->SetLumiLength(0.);
    if (grp->IsPA() || grp->IsAP()) { 
      ZDC->SetpAsystem();
      ZDC->SetBeamEnergy(82.*grp->beamEnergy/208.);
    }
  }
  if (iTRD) {
    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
    AliTRDgeometry *geoTRD = TRD->GetGeometry();
    // Total of 18 super modules. We turn them all off by default 
    for (Int_t i = 0; i < 18; i++) geoTRD->SetSMstatus(i, 0);

    // '09-'10 had 7 super modules 
    geoTRD->SetSMstatus( 0,1);
    geoTRD->SetSMstatus( 1,1);
    geoTRD->SetSMstatus( 7,1);
    geoTRD->SetSMstatus( 8,1);
    geoTRD->SetSMstatus( 9,1);
    geoTRD->SetSMstatus(14,1);//?
    geoTRD->SetSMstatus(17,1);

    // In jan '11 3 more were added 
    if (grp->Year() > 2010) { 
      geoTRD->SetSMstatus(11, 1);
      geoTRD->SetSMstatus(15, 1);
      geoTRD->SetSMstatus(16, 1);//?
    }

    // In the 2012 shutdow 3 more were added 
    if (grp->Year() > 2012) { 
      geoTRD->SetSMstatus( 2,1);
      geoTRD->SetSMstatus( 3,1);
      geoTRD->SetSMstatus( 6,1);
    }
    if (grp->Year() > 2014) {
      geoTRD->SetSMstatus( 4,1);
      geoTRD->SetSMstatus( 5,1);
      geoTRD->SetSMstatus(10,1);
      geoTRD->SetSMstatus(12,1);
      geoTRD->SetSMstatus(13,1);
    }      
  }
  if (iFMD)    new AliFMDv1("FMD", "normal FMD");
  if (iMUON) {
    AliMUON *MUON = new AliMUONv1("MUON", "default");
    MUON->SetTriggerEffCells(1); // not needed if raw masks
    MUON->SetTriggerResponseV1(2);
  }
  if (iPHOS)   new AliPHOSv1("PHOS", "noCPV_Modules123");
  if (iPMD)    new AliPMDv1("PMD", "normal PMD");
  if (iT0)     new AliT0v1("T0", "T0 Detector");
  if (iEMCAL)  new AliEMCALv2("EMCAL", "EMCAL_COMPLETE12SMV1");
  if (iACORDE) new AliACORDEv1("ACORDE", "normal ACORDE");
  if (iVZERO)  new AliVZEROv7("VZERO", "normal VZERO");
}




// 
// EOF
// 
