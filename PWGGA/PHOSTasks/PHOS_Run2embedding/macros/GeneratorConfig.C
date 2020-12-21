/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Generator configuration script
 *
 */

enum EGenerator_t {
  kGeneratorDefault,
  kGeneratorPythia6,
  kGeneratorPythia6_Perugia2011,
  kGeneratorPythia8,
  kGeneratorPythia8_Monash2013,
  kGeneratorPythia8_Monash2013_Rsn001, // [ALIROOT-6685]
  kGeneratorPhojet,
  kGeneratorEPOSLHC_pp,
  kGeneratorHijing,
  kGeneratorHijing_Rsn002a, kGeneratorHijing_Rsn002b, kGeneratorHijing_Rsn002c, // [ALIROOT-6721] [ALIROOT-6722]
  kGeneratorCustom,
  kNGenerators
};

const Char_t *GeneratorName[kNGenerators] = {
  "Default",
  "Pythia6",
  "Pythia6_Perugia2011",
  "Pythia8",
  "Pythia8_Monash2013",
  "Pythia8_Monash2013_Rsn001",
  "Phojet",
  "EPOSLHC_pp",
  "Hijing",
  "Hijing_Rsn002a", "Hijing_Rsn002b", "Hijing_Rsn002c",
  "Custom"
};

enum ETrigger_t {
  kTriggerDefault,
  kTriggerPP,
  kTriggerPbPb,
  kNTriggers
};

const Char_t *TriggerName[kNTriggers] = {
  "ocdb",
  "p-p",
  "Pb-Pb"
};

enum EPythiaTune_t {
  kPerugia2011 = 350,
  kMonash2013  = 14
};

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

// functions
AliGenerator *GeneratorCocktail(TString projN, Int_t projA, Int_t projZ, TString targN, Int_t targA, Int_t targZ);
AliGenerator *GeneratorInjector(Int_t ninj, Int_t pdg, Float_t ptmin, Float_t ptmax, Float_t ymin, Float_t ymax); 
AliGenerator *GeneratorPythia6(Int_t tune = 0, Int_t ntrig = 0, Int_t *trig = NULL);
AliGenerator *GeneratorPythia8(Int_t tune = 0, Int_t ntrig = 0, Int_t *trig = NULL);
AliGenerator *GeneratorPhojet();
AliGenerator *GeneratorEPOSLHC(TString system);
AliGenerator *GeneratorHijing();

/*****************************************************************/

// global variables
static TString comment;

/*****************************************************************/

GeneratorConfig(Int_t tag, Int_t run)
{

  AliGenerator *gen = NULL;
  
  switch (tag) {

    // Default
  case kGeneratorDefault:
    abort();
    return;
    
    // Pythia6 (Perugia2011)
  case kGeneratorPythia6:
  case kGeneratorPythia6_Perugia2011:
    gen = GeneratorPythia6(kPerugia2011);
    break;

    // Pythia8 (Monash2013)
  case kGeneratorPythia8:
  case kGeneratorPythia8_Monash2013:
    gen = GeneratorPythia8(kMonash2013);
    break;
    
    // Pythia8 (Monash2013) - Rsn001
  case kGeneratorPythia8_Monash2013_Rsn001:
    AliGenCocktail *ctl = GeneratorCocktail("p", 1, 1, "p", 1, 1);
    // pythia8
    AliGenerator   *py8 = GeneratorPythia8(kMonash2013);
    ctl->AddGenerator(py8, "Pythia8 (Monash2013)", 1.);
    // randomly injected particles
    Int_t pdglist[] = {225, 3124, -3124, 9010221}; // f2(1270), Lambda(1520), Lambda_bar(1520), f0(980)
    Int_t pdg = pdglist[uidConfig % (sizeof(pdglist) / 4)]; // select according to unique ID
    inj = GeneratorInjector(1, pdg, 0., 15., -0.6, 0.6);
    ctl->AddGenerator(inj, "Injector (Rsn001)", 1.);
    //
    gen = ctl;
    break;
    
    // Phojet
  case kGeneratorPhojet:
    gen = GeneratorPhojet();
    break;

    // EPOSLHC (pp)
  case kGeneratorEPOSLHC_pp:
    gen = GeneratorEPOSLHC("pp");
    break;

    // Hijing
  case kGeneratorHijing:
    gen = GeneratorHijing();
    break;

    // Hijing - Rsn002
  case kGeneratorHijing_Rsn002a:
  case kGeneratorHijing_Rsn002b:
  case kGeneratorHijing_Rsn002c:
    Int_t ninjlist[3] = {25, 7, 3};
    Int_t ninj = ninjlist[tag - kGeneratorHijing_Rsn002a];
    AliGenCocktail *ctl  = GeneratorCocktail("A", 208, 82, "A", 208, 82);
    AliGenerator   *hij  = GeneratorHijing();
    AliGenerator   *inj1 = GeneratorInjector(ninj,  3124, 0., 10., -0.6, 0.6);
    AliGenerator   *inj2 = GeneratorInjector(ninj, -3124, 0., 10., -0.6, 0.6);
    ctl->AddGenerator(hij, "Hijing (central)", 1.);
    ctl->AddGenerator(inj1, "Injector (Rsn002)", 1.);
    ctl->AddGenerator(inj2, "Injector (Rsn002)", 1.);
    gen = ctl;
    break;
    
    // Custom
  case kGeneratorCustom:
    if ((gROOT->LoadMacro("Sim/GeneratorCustom.C")) != 0) {
      printf("ERROR: cannot find GeneratorCustom.C\n");
      abort();
      return;
    }
    gen = GeneratorCustom();

  }
  
  //  gener->SetOrigin(0.075, 0.522, -0.884); // R+HACK                        
  //  gener->SetSigma(65e-4, 65e-4, 5.); // R+HACK                                           
  gen->SetVertexSmear(kPerEvent);
  gen->Init();
  printf(">>>>> Generator Configuration: %s \n", comment.Data());
  // Set the trigger configuration: proton-proton
  AliSimulation::Instance()->SetTriggerConfig(TriggerName[triggerConfig]);
  printf(">>>>> Trigger configuration:   %s \n", TriggerName[triggerConfig]);
 
}

/*** PYTHIA 6 ****************************************************/

AliGenerator *
GeneratorPythia6(Int_t tune, Int_t ntrig, Int_t *trig)
{
  comment = comment.Append(" | Pythia6 low-pt");
  //
  // Pythia
  AliGenPythia* pythia = new AliGenPythia(-1); 
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-12.,12.);
  pythia->SetPtRange(0,1000.);
  pythia->SetProcess(kPyMb);
  pythia->SetEnergyCMS(energyConfig);
  pythia->SetCrossingAngle(0,crossingConfig);
  //
  // Tune
  if (tune > 0) {
    comment = comment.Append(Form(" | tune %d", tune));
    pythia->SetTune(tune); 
    //    pythia->UseNewMultipleInteractionsScenario();
  }
  //
  // Trigger particles
  if (ntrig > 0) {
    Int_t pdg = trig[gRandom->Integer(ntrig)];
    comment = comment.Append(Form(" | %s enhanced", TDatabasePDG::Instance()->GetParticle(pdg)->GetName()));
    pythia->SetTriggerParticle(pdg, 1.2);
  }
  //
  return pythia;
}

/*** PYTHIA 8 ****************************************************/

AliGenerator *
GeneratorPythia8(Int_t tune, Int_t ntrig, Int_t *trig)
{
  //
  // Libraries
  gSystem->Load("libpythia8.so");
  gSystem->Load("libAliPythia8.so");
  //
  //
  comment = comment.Append(" | Pythia8 low-pt");
  //
  // Pythia
  AliGenPythiaPlus *pythia = new AliGenPythiaPlus(AliPythia8::Instance()); 
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-12.,12.);
  pythia->SetPtRange(0,1000.);
  pythia->SetProcess(kPyMbDefault); // pythia->SetProcess(kPyMb);
  pythia->SetEnergyCMS(energyConfig);
  pythia->SetCrossingAngle(0, crossingConfig);
  //
  // Initialize
  pythia->SetEventListRange(-1, 2); 
  (AliPythia8::Instance())->ReadString("Random:setSeed = on");
  (AliPythia8::Instance())->ReadString(Form("Random:seed = %ld", seedConfig % 900000000)); 
  //
  // Tune
  if (tune > 0) {
    comment = comment.Append(Form(" | tune %d", tune));
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d", tune));
  }
  //
  // Trigger particles
  if (ntrig > 0) {
    Int_t pdg = trig[gRandom->Integer(ntrig)];
    comment = comment.Append(Form(" | %s enhanced", DatabasePDG::Instance()->GetParticle(pdg)->GetName()));
    pythia->SetTriggerParticle(pdg, 1.2);
  }
  //
  return pythia;
}

/*** PHOJET ****************************************************/

AliGenerator *
GeneratorPhojet()
{
  //
  // Libraries
  gSystem->Load("libDPMJET");
  gSystem->Load("libTDPMjet");
  //
  comment = comment.Append(" | Phojet low-pt");
  //                                                                                      
  //    DPMJET                                                                            
  AliGenDPMjet* dpmjet = new AliGenDPMjet(-1);
  dpmjet->SetMomentumRange(0, 999999.);
  dpmjet->SetThetaRange(0., 180.);
  dpmjet->SetYRange(-12.,12.);
  dpmjet->SetPtRange(0,1000.);
  dpmjet->SetProcess(kDpmMb);
  dpmjet->SetEnergyCMS(energyConfig);
  dpmjet->SetCrossingAngle(0,crossingConfig);
  return dpmjet;
}

/*** EPOSLHC ****************************************************/

AliGenerator *
GeneratorEPOSLHC(TString system)
{
  
  if (system.EqualTo("pp")) {
    comment = comment.Append(" | EPOS-LHC");
    Float_t beamEnergy = energyConfig / 2.;
    TString fifoname = Form("/tmp/crmceventfifo%d", gRandom->Integer(kMaxInt));
    gROOT->ProcessLine(Form(".! rm -rf %s", fifoname.Data()));
    gROOT->ProcessLine(Form(".! mkfifo %s", fifoname.Data()));
    gROOT->ProcessLine(Form(".! sh $ALIDPG_ROOT/MC/EXTRA/gen_eposlhc_pp.sh %s %d %f %f &",
			    fifoname.Data(), neventsConfig, beamEnergy, beamEnergy));
  }
  else {
    printf("EPOSLHC not implemented for %s system\n", system.Data());
    abort();
  }
  
  AliGenReaderHepMC *reader = new AliGenReaderHepMC();
  reader->SetFileName("crmceventfifo");
  AliGenExtFile *gener = new AliGenExtFile(-1);
  gener->SetReader(reader);
  
  return gener;
}

/*** HIJING ****************************************************/

AliGenerator * 
GeneratorHijing()
{
  //
  // Libraries
  gSystem->Load("libHIJING");
  gSystem->Load("libTHijing");

  comment = comment.Append(Form(" | HIJING (b = %f-%f fm)", bminConfig, bmaxConfig));
  AliGenHijing *gener = new AliGenHijing(-1);
  // centre of mass energy
  gener->SetEnergyCMS(energyConfig);
  gener->SetImpactParameterRange(bminConfig, bmaxConfig);
  // reference frame
  gener->SetReferenceFrame("CMS");
  // projectile
  gener->SetProjectile("A", 208, 82);
  gener->SetTarget    ("A", 208, 82);
  // tell hijing to keep the full parent child chain
  gener->KeepFullEvent();
  // enable jet quenching
  gener->SetJetQuenching(0);
  // enable shadowing
  gener->SetShadowing(1);
  // neutral pion and heavy particle decays switched off
  gener->SetDecaysOff(1);
  // Don't track spectators
  gener->SetSpectators(0);
  // kinematic selection
  gener->SetSelectAll(0);
  gener->SetPtHardMin (2.3);
  return gener;
}

/*** COCKTAIL ****************************************************/

AliGenerator * 
GeneratorCocktail(TString projN, Int_t projA, Int_t projZ,
		  TString targN, Int_t targA, Int_t targZ)
{
  comment = comment.Append(Form(" | cocktail (%s-%s)", projN.Data(), targN.Data()));
  //
  AliGenCocktail *ctl = new AliGenCocktail();
  ctl->SetProjectile(projN, projA, projZ);
  ctl->SetTarget(targN, targA, targZ);
  ctl->SetEnergyCMS(energyConfig);
  return ctl;
}

/*** INJECTOR ****************************************************/

AliGenerator * 
GeneratorInjector(Int_t ninj, Int_t pdg, Float_t ptmin, Float_t ptmax, Float_t ymin, Float_t ymax)
{
  comment = comment.Append(Form(" | injected (pdg=%d, %d particles)", pdg, ninj));
  //
  // Injected particles
  AliGenBox *box = new AliGenBox(ninj);
  box->SetPart(pdg);
  box->SetPtRange(ptmin, ptmax);
  box->SetYRange(ymin, ymax);
  box->SetPhiRange(0., 360.);
  return box;
}

