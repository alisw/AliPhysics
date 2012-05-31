///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TTherminator: global interface class to the Therminator model.            //
// Initialized global variables and runs the event generation                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
/******************************************************************************
 *                      T H E R M I N A T O R                                 *
 *                   THERMal heavy-IoN generATOR                              *
 *                           version 1.0                                      *
 *                                                                            *
 * Authors of the model: Wojciech Broniowski, Wojciech.Broniowski@ifj.edu.pl, *
 *                       Wojciech Florkowski, Wojciech.Florkowski@ifj.edu.pl  *
 * Authors of the code:  Adam Kisiel, kisiel@if.pw.edu.pl                     *
 *                       Tomasz Taluc, ttaluc@if.pw.edu.pl                    *
 * Code designers: Adam Kisiel, Tomasz Taluc, Wojciech Broniowski,            *
 *                 Wojciech Florkowski                                        *
 *                                                                            *
 * For the detailed description of the program and furhter references         * 
 * to the description of the model plesase refer to: nucl-th/0504047,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0504047                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-th/0504047 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#include <TTherminator.h>
#include <fstream>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TClass.h>

ReadPar *sRPInstance;
STR      sRPFileName;
int      sNumEvents;
int      sRunType; 
int      sTables;
int      sModel;
int      sIntegrateSample;

ClassImp(TTherminator)

TTherminator::TTherminator():
  fCalka(0),
  fEvent(0),
  fPartDB(0)
{
  // Default constructor
  fPartDB = new ParticleDB();
}
TTherminator::TTherminator(const TTherminator & therm) :
  TGenerator(therm), 
  fCalka(0),
  fEvent(0),
  fPartDB(0)
{
  // Copy constructor
  //  fPartDB = new ParticleDB();
  if (fCalka) delete fCalka;
  fCalka = new Integrator(*therm.fCalka);
  if (fEvent) delete fEvent;
  fEvent = new Event(*therm.fEvent);
  if (fPartDB) delete fPartDB;
  fPartDB = new ParticleDB();
}
TTherminator& TTherminator::operator=(const TTherminator & therm)
{
  if (this != &therm) {
    fCalka = therm.fCalka;
    fEvent = therm.fEvent;
    delete fPartDB;
    fPartDB = new ParticleDB();
  }

  return *this;
}

TTherminator::~TTherminator()
{
  // Destructor
  if (sRPInstance) delete sRPInstance;
  if (fEvent) delete fEvent;
  if (fCalka) delete fCalka;
  if (fPartDB) delete fPartDB;
}

void TTherminator::ReadParameters()
{
  // Read in global model parameters
  STR tModel;
  STR tTable;
  
  try {
    sNumEvents = atoi(sRPInstance->getPar("NumberOfEvents").Data());
    tModel = sRPInstance->getPar("FreezeOutModel");
    if (tModel.Contains("SingleFreezeOut"))
      sModel = 0;
/*MCH begin*/
    else if (tModel.Contains("Lhyquid3D"))
      sModel = 10;
    else if (tModel.Contains("Lhyquid2D"))
      sModel = 11;
/*MCH end*/      
    else if (tModel.Contains("BlastWaveVTDelay"))
      sModel = 6;
    else if (tModel.Contains("BlastWaveVT"))
      sModel = 2;
    else if (tModel.Contains("BlastWaveVLinearFormation"))
      sModel = 7;
    else if (tModel.Contains("BlastWaveVLinearDelay"))
      sModel = 8;
    else if (tModel.Contains("BlastWaveVLinear"))
      sModel = 4;
    else if (tModel.Contains("FiniteHydro"))
      sModel = 5;
    else {
      PRINT_MESSAGE("Unknown model type: " << tModel.Data());
      PRINT_MESSAGE("Please provide the proper name");
      exit(0);
    }
    if (atoi(sRPInstance->getPar("Randomize").Data()) == 1)
      sRunType = 4;
    else
      sRunType = 3;
    tTable = sRPInstance->getPar("TableType");
    if (tTable.Contains("Mathematica"))
      sTables = 0;
    else if (tTable.Contains("SHARE"))
      sTables = 1;
    else {
      PRINT_MESSAGE("Unknown table type: " << tTable.Data());
      exit(0);
    }
    sIntegrateSample = atoi(sRPInstance->getPar("NumberOfIntegrateSamples").Data());
    fInputDir = sRPInstance->getPar("InputDirSHARE");
  }
  catch (STR tError) {
    PRINT_DEBUG_1("RunBFPW::ReadParameters - Caught exception " << tError);
    PRINT_MESSAGE("Did not find one of the neccessary parameters in the parameters file.");
    exit(0);
  }
}

void        TTherminator::Initialize(){
  // Initialize global variables
  // and particle and decay tables
  Parser       *tParser;
  
  sRPFileName = "therminator.in";
  sRPInstance = new ReadPar(sRPFileName.Data());

  ReadParameters();
  
  tParser = new Parser(fPartDB);

  tParser->ReadShare();
    
  fCalka = new Integrator(sIntegrateSample);
  fCalka->ReadMultInteg(fPartDB);


  // Read in particle table from share
  // and add missing particle type to TDatabasePDG

  char str[50];
  char str1[200];
    
  //  ParticleType *tPartBuf;

  TDatabasePDG *tInstance = TDatabasePDG::Instance();
  TParticlePDG *tParticleType;

  //  AliWarning(Form("Reading particle types from particles.data"));

  ifstream in((fInputDir+"/"+"particles.data").Data());
  //  ifstream in("particles.data");
  
  int charge;
    
  int number=0;
  if ((in) && (in.is_open()))
    {
      //START OF HEAD-LINE
      in.ignore(200,'\n');
      in.ignore(200,'\n');
      in.ignore(200,'\n');
      //END OF HEAD-LINE
      
      while (in>>str)
	{
	  if (/*(*str == '#')||*/(*str<65)||(*str>122))
	    {
	      in.getline(str1,200);
	      continue;
	    }
	  double mass, gamma, spin, tI3, tI, q, s, aq, as, c, ac, mc;
	  
	  in>>mass>>gamma>>spin>>tI>>tI3>>q>>s>>aq>>as>>c>>ac>>mc;
	  number++;
	  tParticleType = tInstance->GetParticle((int) mc);
	  if (!tParticleType) {
	    charge = 0;
	    if (strstr(str, "plu")) charge = 1;
	    if (strstr(str, "min")) charge = -1;
	    if (strstr(str, "plb")) charge = -1;
	    if (strstr(str, "mnb")) charge = 1;
	    if (strstr(str, "plp")) charge = 2;
	    if (strstr(str, "ppb")) charge = -2;
	    tInstance->AddParticle(str, str, mass, gamma == 0.0 ? 1:0, gamma, charge , "meson", (int) mc);
	    //	    printf("Added particle %s with PDG PID %d charge %d", str, (int) mc, charge);
	    //AliWarning(Form("Added particle %s with PDG PID %d charge %d", str, (int) mc, charge));
	    //	    AliWarning(Form("Quantum numbers q s c aq as ac tI3 %lf %lf %lf %lf %lf %lf %lf", q, s, c, aq, as, ac, tI3));

	  }
	}
      in.close();
    }

  delete tParser;
}

void        TTherminator::GenerateEvent()
{
  // Run single event generation
  if ((sRunType == 3) || (sRunType == 4))
    {
      if (sRunType==4)
	fCalka->Randomize();

      if (!fEvent) fEvent = new Event(fPartDB, fCalka);
      if (sRunType == 4) fEvent->Randomize();
      
      fEvent->GenerateEvent();
      fEvent->DecayParticles();
      //      fEvent->WriteEvent(0);
    }  
}

Int_t       TTherminator::ImportParticles(TClonesArray *particles, Option_t */*option*/)
{
  // Import particles from a generated event into an external array
  const double kFmToGev = 0.197327;

  if (particles == 0) return 0;
  TClonesArray &particlesR = *particles;
  particlesR.Clear();
  Int_t nump = 0;
  if (!fEvent) return 0;
  Int_t numpart = fEvent->GetParticleCount();
  printf("\n TTherminator: Therminator stack contains %d particles.\n", numpart);
  for (Int_t iPart=0; iPart<numpart; iPart++) {
    Particle *tPart = fEvent->GetParticleOfCount(iPart);
    Int_t tFather;
    tFather = tPart->GetFather();

    if (tFather> -1) {
      TParticle *mother = (TParticle*) (particlesR.UncheckedAt(tFather+1));	   
      mother->SetLastDaughter(iPart);
      if (mother->GetFirstDaughter()==-1)
 	mother->SetFirstDaughter(iPart);
    }
    nump++;
    //    printf("Putting %d %d %lf %d\n", tPart->GetParticleType()->GetPDGCode(), iPart, tPart->px, tPart);
    new (particlesR[iPart]) TParticle(tPart->GetParticleType()->GetPDGCode(), tPart->HadDecayed(),
 				      tFather, -1, -1, -1,
 				      tPart->px, tPart->py, tPart->pz, tPart->GetEnergy() ,
 				      tPart->rx*1.e-13*kFmToGev, tPart->ry*1.e-13*kFmToGev, tPart->rz*1.e-13*kFmToGev, tPart->rt*1.e-13*kFmToGev/3e10);
    particlesR[iPart]->SetUniqueID(iPart);
  }
  
  return nump;
}

TObjArray*  TTherminator::ImportParticles(Option_t */*option*/)
{
  // Import generated particles into an internal array
  const double kFmToGev = 0.197327;

  Int_t nump = 0;
  fParticles->Clear();
  if (!fEvent) return 0;
  Int_t numpart = fEvent->GetParticleCount();
  printf("\n TTherminator: Therminator stack contains %d particles.", numpart);
  fEvent->InitParticleScan();
  for (Int_t iPart=0; iPart<numpart; iPart++) {
    Particle *tPart = fEvent->GetNextParticle();
    Int_t tFather;
    tFather = tPart->GetFather();
    
    if (tFather> -1) {
      TParticle *mother = (TParticle*) (fParticles->UncheckedAt(tFather));	   
      mother->SetLastDaughter(iPart);
      if (mother->GetFirstDaughter()==-1)
 	mother->SetFirstDaughter(iPart);
    }
    TParticle* p = new TParticle(tPart->GetParticleType()->GetPDGCode(), tPart->HadDecayed(),
				 tFather, -1, -1, -1,
				 tPart->px, tPart->py, tPart->pz, tPart->GetEnergy() ,
				 tPart->rx*1.e-13*kFmToGev, tPart->ry*1.e-13*kFmToGev, tPart->rz*1.e-13*kFmToGev, tPart->rt*1.e-13*kFmToGev/3e10);
    p->SetUniqueID(iPart);
    fParticles->Add(p);
  }
  nump = fEvent->GetParticleCount();
  
  return fParticles;
  
}
