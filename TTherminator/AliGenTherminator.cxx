// ALICE event generator based on the THERMINATOR model
// It reads the test output of the model and puts it onto
// the stack
// Author: Adam.Kisiel@cern.ch

#include <iostream>
#include <fstream>
#include <sstream>
#include <TClonesArray.h>
#include <TMCProcess.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliConst.h"
#include "AliDecayer.h"
#include "AliGenEventHeader.h"
#include "AliGenTherminator.h"
#include "AliLog.h"
#include "AliRun.h"

ClassImp(AliGenTherminator)

using namespace std;

AliGenTherminator::AliGenTherminator():
  AliGenMC(),
  fNt(0),
  fEventNumber(0),
  fFileName(""),
  fFreezeOutModel(""),
  fTemperature(0.1656),
  fMiuI(-0.0009),
  fMiuS(0.0),
  fMiuB(0.0008),
  fAlfaRange(6.0),
  fRapRange(4.0),
  fRhoMax(8.0),
  fTau(8.0),
  fBWA(0.0),
  fBWVt(1.41),
  fBWDelay(0.0)
{
  // Default constructor
  fParticles = new TClonesArray("TParticle",20000);    
}
AliGenTherminator::AliGenTherminator(Int_t npart):
  AliGenMC(npart),
  fNt(0),
  fEventNumber(0),
  fFileName(""),
  fFreezeOutModel(""),
  fTemperature(0.1656),
  fMiuI(-0.0009),
  fMiuS(0.0),
  fMiuB(0.0008),
  fAlfaRange(6.0),
  fRapRange(4.0),
  fRhoMax(8.0),
  fTau(8.0),
  fBWA(0.0),
  fBWVt(1.41),
  fBWDelay(0.0)
{
  // Constructor specifying the size of the particle table
  fParticles = new TClonesArray("TParticle",20000);    
  fNprimaries = 0;
}

AliGenTherminator::~AliGenTherminator()
{
  //  AliGenMC::~AliGenMC();
  //  if (fTherminator) delete fTherminator;
}

void AliGenTherminator::Generate()
{
  // Run single event generation with the Therminator model
  AliWarning("Generating event from AliGenTherminator");

  Float_t polar[3]    =   {0,0,0};
  Float_t origin[3]   =   {0,0,0};
  Float_t origin0[3]  =   {0,0,0};
  Float_t p[3];
  Float_t mass, energy;

  Int_t nt  = 0;
  Int_t j, kf, ks, imo;
  kf = 0;

  Vertex();
  for (j=0; j < 3; j++) origin0[j] = fVertex[j];

  // Generate one event

  ((TTherminator *) fMCEvGen)->GenerateEvent();
  AliWarning("Generated");
  ((TTherminator *) fMCEvGen)->ImportParticles(fParticles);

  Int_t np = fParticles->GetEntriesFast();
  AliWarning(Form("Imported %d particles", np));

  TParticle *iparticle;
  
  for (int i = 0; i < np; i++) {
    iparticle = (TParticle *) fParticles->At(i);
    Bool_t  hasMother   = (iparticle->GetFirstMother()     >=0);
    Bool_t  hasDaughter = (iparticle->GetFirstDaughter()   >=0);
    
    kf   = iparticle->GetPdgCode();
    ks   = iparticle->GetStatusCode();
    p[0] = iparticle->Px();
    p[1] = iparticle->Py();
    p[2] = iparticle->Pz();
    mass = TDatabasePDG::Instance()->GetParticle(kf)->Mass();
    energy = sqrt(mass*mass + p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    origin[0] = origin0[0]+iparticle->Vx();
    origin[1] = origin0[1]+iparticle->Vy();
    origin[2] = origin0[2]+iparticle->Vz();
	      
    imo = -1;
    TParticle* mother = 0;
    if (hasMother) {
      imo = iparticle->GetFirstMother();
      mother = (TParticle *) fParticles->At(imo);
    } // if has mother   
    Bool_t tFlag = (!hasDaughter);

    printf("Pushing Track %d with status %d mother %d\n", kf, tFlag, imo);
    PushTrack(tFlag,imo,kf,
	      p[0],p[1],p[2],energy,
	      origin[0],origin[1],origin[2],iparticle->T()*0.197327*1e-13/300000000,
	      polar[0],polar[1],polar[2],
	      hasMother ? kPDecay:kPNoProcess,nt);
    
    fNprimaries++;
    KeepTrack(nt);
  }

  SetHighWaterMark(fNprimaries);

  TArrayF eventVertex;
  eventVertex.Set(3);
  eventVertex[0] = origin0[0];
  eventVertex[1] = origin0[1];
  eventVertex[2] = origin0[2];

  // Header
  AliGenEventHeader* header = new AliGenEventHeader("Therminator");
  // Event Vertex
  header->SetPrimaryVertex(eventVertex);
  header->SetNProduced(fNprimaries);
  gAlice->SetGenEventHeader(header); 
}

void AliGenTherminator::Init()
{
  // Initialize global variables and
  // particle and decay tables
  if (fFileName.Length() == 0)
    fFileName = "event.out";
  ReadShareParticleTable();

  SetMC(new TTherminator());
  
  AliGenMC::Init();
  ((TTherminator *) fMCEvGen)->Initialize();
}
void AliGenTherminator::SetFileName(const char *infilename)
{
  // Set parameter filename
  fFileName = infilename;
}
void AliGenTherminator::SetEventNumberInFile(int evnum)
{
  // Set number of events to generate - default: 1
  fEventNumber = evnum;
}

void AliGenTherminator::ReadShareParticleTable()
{
  // Read in particle table from share
  // and add missing particle type to TDatabasePDG

  char str[50];
  char str1[200];
    
  TDatabasePDG *tInstance = TDatabasePDG::Instance();
  TParticlePDG *tParticleType;

  AliWarning(Form("Reading particle types from particles.data"));
  ifstream in("particles.data");
  
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
	    AliWarning(Form("Added particle %s with PDG PID %d charge %d", str, (int) mc, charge));
	    //	    AliWarning(Form("Quantum numbers q s c aq as ac tI3 %lf %lf %lf %lf %lf %lf %lf", q, s, c, aq, as, ac, tI3));

	  }
	}
      in.close();
    }
  CreateTherminatorInputFile();
}

void AliGenTherminator::CreateTherminatorInputFile()
{
  // Create Therminator input file
  ofstream *ostr = new ofstream("therminator.in");
  (*ostr) << "NumberOfEvents = 1" << endl;
  (*ostr) << "Randomize = 1" << endl;
  (*ostr) << "TableType = SHARE" << endl;
  (*ostr) << "InputDirSHARE = ." << endl;
  (*ostr) << "EventOutputFile = " << fFileName.Data() << endl;
  (*ostr) << "FreezeOutModel = " << fFreezeOutModel.Data() << endl;
  (*ostr) << "BWVt = " << fBWVt << endl;
  (*ostr) << "Tau = " << fTau << endl;
  (*ostr) << "RhoMax = " << fRhoMax << endl;
  (*ostr) << "Temperature = " << fTemperature << endl;
  (*ostr) << "MiuI = " << fMiuI << endl;
  (*ostr) << "MiuS = " << fMiuS << endl;
  (*ostr) << "MiuB = " << fMiuB << endl;
  (*ostr) << "AlphaRange = " << fAlfaRange << endl;
  (*ostr) << "RapidityRange = " << fRapRange << endl;
  (*ostr) << "NumberOfIntegrateSamples = 1000000" << endl;
}

void AliGenTherminator::SetModel(const char *model)
{
  // Set the freeze-out model to use
  fFreezeOutModel = model;
}
