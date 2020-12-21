/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Central configuration script
 *
 */

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

// global variables

static Int_t   runNumber       = 0;         // run number
static Int_t   neventsConfig   = 4;       // number of events
static Int_t   magnetConfig    = 0;         // magnetic field
static Int_t   detectorConfig  = 0;         // detector
static Int_t   generatorConfig = 0 ;         // MC generator
static Float_t energyConfig    = 0.;        // CMS energy
static Float_t triggerConfig   = 0.;        // trigger
static Float_t bminConfig      = 20.;       // impact parameter min
static Float_t bmaxConfig      = 0.;        // impact parameter max
static Float_t crossingConfig  = 0.;        // 2.8e-4 // crossing angle
static Int_t   seedConfig      = 123456789; // random seed
static Int_t   uidConfig       = 1;         // unique ID

/*****************************************************************/

void
Config()
{

  /* initialise */
  gROOT->LoadMacro("Sim/DetectorConfig.C");
  gROOT->LoadMacro("Sim/GeneratorConfig.C");
  ProcessEnvironment();

  /* verbose */
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  printf(">>>>>       run number: %d \n", runNumber);
  printf(">>>>> number of events: %d \n", neventsConfig);
  printf(">>>>>   magnetic field: %s \n", MagnetName[magnetConfig]);
  printf(">>>>>         detector: %s \n", DetectorName[detectorConfig]);
  printf(">>>>>     MC generator: %s \n", GeneratorName[generatorConfig]);
  printf(">>>>>       CMS energy: %f \n", energyConfig);
  printf(">>>>>          trigger: %s \n", TriggerName[triggerConfig]);
  printf(">>>>>            b-min: %f \n", bminConfig);
  printf(">>>>>            b-max: %f \n", bmaxConfig);
  printf(">>>>>   crossing angle: %f \n", crossingConfig);
  printf(">>>>>      random seed: %d \n", seedConfig);
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");

  /* load libraries */
  LoadLibraries();

  /* setup geant */
  new TGeant3TGeo("C++ Interface to Geant3");

  /* create galice.root */
  CreateGAlice();

  /* configure detector */
  DetectorConfig(detectorConfig, runNumber);

  /* configure MC generator */
  GeneratorConfig(generatorConfig, runNumber);
  GeneratorOptions();
}

/*****************************************************************/

Float_t
SetEnergyFromGRP()
{
  AliCDBEntry *cdbe = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject *grpd = dynamic_cast<AliGRPObject*>(cdbe->GetObject()); 
  return (grpd->GetBeamEnergy() * 2.);
}

/*****************************************************************/

void
ProcessEnvironment()
{

  // run number
  if (gSystem->Getenv("CONFIG_RUN"))
    runNumber = atoi(gSystem->Getenv("CONFIG_RUN"));

  // number of events configuration
  neventsConfig = 200;
  if (gSystem->Getenv("CONFIG_NEVENTS"))
    neventsConfig = atoi(gSystem->Getenv("CONFIG_NEVENTS"));

  // magnetic field configuration
  magnetConfig = kMagnetDefault;
  if (gSystem->Getenv("CONFIG_MAGNET")) {
    Bool_t valid = kFALSE;
    for (Int_t imag = 0; imag < kNMagnets; imag++)
      if (strcmp(gSystem->Getenv("CONFIG_MAGNET"), MagnetName[imag]) == 0) {
	magnetConfig = imag;
	valid = kTRUE;
	break;
      }
    if (!valid) {
      printf(">>>>> Unknown magnetic field configuration: %s \n", gSystem->Getenv("CONFIG_MAGNET"));
      abort();
    }
  }
	
  // detector configuration
  detectorConfig = kDetectorDefault;
  if (gSystem->Getenv("CONFIG_DETECTOR")) {
    Bool_t valid = kFALSE;
    for (Int_t idet = 0; idet < kNDetectors; idet++)
      if (strcmp(gSystem->Getenv("CONFIG_DETECTOR"), DetectorName[idet]) == 0) {
	detectorConfig = idet;
	valid = kTRUE;
	break;
      }
    if (!valid) {
      printf(">>>>> Unknown detector configuration: %s \n", gSystem->Getenv("CONFIG_DETECTOR"));
      abort();
    }
  }
	
//   // generator configuration
   generatorConfig = kGeneratorCustom;
//   if (gSystem->Getenv("CONFIG_GENERATOR")) {
//     Bool_t valid = kFALSE;
//     for (Int_t igen = 0; igen < kNGenerators; igen++)
//       if (strcmp(gSystem->Getenv("CONFIG_GENERATOR"), GeneratorName[igen]) == 0) {
// 	generatorConfig = igen;
// 	valid = kTRUE;
// 	break;
//       }
//     if (!valid) {
//       printf(">>>>> Unknown MC generator configuration: %s \n", gSystem->Getenv("CONFIG_GENERATOR"));
//       abort();
//     }
//   }
  
  // energy configuration
  energyConfig = SetEnergyFromGRP();
  if (gSystem->Getenv("CONFIG_ENERGY"))
    energyConfig = atoi(gSystem->Getenv("CONFIG_ENERGY"));
  if (energyConfig <= 0) {
    printf(">>>>> Invalid CMS energy: %f \n", energyConfig);
    abort();
  }

  // trigger configuration
  triggerConfig = kGeneratorDefault;
  if (gSystem->Getenv("CONFIG_TRIGGER")) {
    Bool_t valid = kFALSE;
    for (Int_t itrg = 0; itrg < kNTriggers; itrg++)
      if (strcmp(gSystem->Getenv("CONFIG_TRIGGER"), TriggerName[itrg]) == 0) {
	triggerConfig = itrg;
	valid = kTRUE;
	break;
      }
    if (!valid) {
      printf(">>>>> Unknown trigger configuration: %s \n", gSystem->Getenv("CONFIG_TRIGGER"));
      abort();
    }
  }
  
  // impact parameter configuration
  bminConfig = 0.;
  if (gSystem->Getenv("CONFIG_BMIN"))
    bminConfig = atoi(gSystem->Getenv("CONFIG_BMIN"));
  if (bminConfig < 0) {
    printf(">>>>> Invalid min impact parameter: %f \n", bminConfig);
    abort();
  }
  bmaxConfig = 20.;
  if (gSystem->Getenv("CONFIG_BMAX"))
    bmaxConfig = atoi(gSystem->Getenv("CONFIG_BMAX"));
  if (bmaxConfig <= bminConfig) {
    printf(">>>>> Invalid max impact parameter: %f \n", bmaxConfig);
    abort();
  }

  // seed configuration
  seedConfig = TDatime().Get();
  if (gSystem->Getenv("CONFIG_SEED"))
    seedConfig = atoi(gSystem->Getenv("CONFIG_SEED"));
  
  // unique ID configuration
  uidConfig = 1;
  if (gSystem->Getenv("CONFIG_UID"))
    uidConfig = atoi(gSystem->Getenv("CONFIG_UID"));
  
}

/*****************************************************************/

void
LoadLibraries()
{
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libgeant321");
  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libAliPythia6");

}

/*****************************************************************/

void
CreateGAlice() 
{
  //=======================================================================
  //  Create the output file
   
  AliRunLoader* rl=0x0;

  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
  if (!rl) {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();
  rl->CdGAFile();
}

/*****************************************************************/

void
GeneratorOptions()
{
  //======================//
  //    Set MC options    //
  //======================//
  
  //
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
  //
  //======================//
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  //
}
