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
#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include "Particle.h"
#include "Event.h"
#include "Integrator.h"
#include "ParticleDecayer.h"
#include "DecayTable.h"
#include "ReadPar.h"

extern ReadPar* sRPInstance;

Event::Event(ParticleDB *aDB, Integrator* aInteg) :
  mCurIter(),
  mScanStarted(0)
{
  mDB = aDB;
  mInteg = aInteg;
  mRandom = new TRandom2();
#ifdef _ROOT_4_
  mRandom->SetSeed2(31851, 14327);
#else
  mRandom->SetSeed(31851);
#endif
  mAverageMultiplicities.resize(mDB->GetParticleTypeCount());
  mMultiplicities.resize(mDB->GetParticleTypeCount());
  ReadParameters();
  ReadMultiplicities();


}

Event::~Event()
{
  
  mParticles.clear();
  delete mRandom;
}

void 
Event::GenerateEvent(int aSeed)
{
  Particle **tPartBuf;

  mParticles.clear();
#ifdef _ROOT_4_
  if (aSeed) mRandom->SetSeed2(aSeed, (aSeed*2) % (7*11*23*31));
#else
  if (aSeed) mRandom->SetSeed(aSeed);
#endif
  GenerateMultiplicities();
  for (unsigned int tIter=0; tIter<mAverageMultiplicities.size(); tIter++) {
    mInteg->Generate(mDB->GetParticleType(tIter), mMultiplicities[tIter], &tPartBuf);
    for (int tBufIter=0; tBufIter<mMultiplicities[tIter]; tBufIter++) {
      AddParticle(tPartBuf[tBufIter]);
      delete tPartBuf[tBufIter];
    }
    free(tPartBuf);
  }
}

int  
Event::GetParticleCount()
{
  return mParticles.size();
}

Particle* 
Event::GetNextParticle()
{
  if (mScanStarted)
    InitParticleScan();
  if (mCurIter == mParticles.end())
    return 0;

  Particle *tBuf;
  tBuf = &(*mCurIter);
  mCurIter++;
  
  return tBuf;
}

void 
Event::InitParticleScan()
{
  mCurIter = mParticles.begin();
  mScanStarted = 1;
}

Particle* 
Event::GetParticleOfCount(int aCount)
{
  if (aCount == 0) {
    InitParticleScan();
  }
  else if (aCount == 1) {
    InitParticleScan();
    mCurIter++;
  }
  else {
    mCurIter--;
    mCurIter++;
    mCurIter++;
  }

  if (mCurIter == mParticles.end())
    return 0;
  
  Particle *tBuf;
  tBuf = &(*mCurIter);
  return tBuf;
}

void 
Event::WriteEvent(int nEvent)
{
  ParticleListIter tPLIter;
  tPLIter = mParticles.begin();
  int tCount = 0;
  ofstream os;
  
  if (nEvent == 0) {
    try {
      os.open(sRPInstance->getPar("EventOutputFile").Data());
    }
    catch (STR e) {
      PRINT_MESSAGE("Very strange: No event output filename ???");
      exit(0);
    }
  }
  else {
    try {
      os.open(sRPInstance->getPar("EventOutputFile").Data(), ios::app);
    }
    catch (STR e) {
      PRINT_MESSAGE("Very strange: No event output filename ???");
      exit(0);
    }
  }
    
  

  if (nEvent == 0) 
    {
      os << "Therminator text output\nall_particles_p_x\nTherminator " << endl;
    }
  os << nEvent+1 << "	" << GetParticleCount() << "	" << "0.000" << "	" << "0.000" << endl;

  while (tPLIter != mParticles.end())
    {
      os.flags(ios::right);
      os.precision(6);
      os.width(6);
      os << tCount << "   ";
      (*tPLIter).WriteParticle(&os);
      os << endl;

      tCount++;
      tPLIter++;
    }
}

  
void 
Event::AddParticle(Particle* aParticle)
{
  if (0) {
    ParticleListIter tPLIter;
    tPLIter = mParticles.begin();
    while (tPLIter != mParticles.end())
      {
	if (((*tPLIter).GetMass() == aParticle->GetMass()) &&
	    ((*tPLIter).GetI3() == aParticle->GetI3()) &&
	    ((*tPLIter).GetBarionN() == aParticle->GetBarionN()) &&
	    ((*tPLIter).GetStrangeness() == aParticle->GetStrangeness())) {
	  break; }
	if ((*tPLIter).GetMass() < aParticle->GetMass()) {
	  break; }
	else if ((*tPLIter).GetMass() == aParticle->GetMass()) {
	  if ((*tPLIter).GetI3() < aParticle->GetI3()) {
	    break; }
	  else if ((*tPLIter).GetI3() == aParticle->GetI3()) {
	    if ((*tPLIter).GetBarionN() < aParticle->GetBarionN()) {
	      break; }
	    else if ((*tPLIter).GetBarionN() == aParticle->GetBarionN()) {
	      if ((*tPLIter).GetStrangeness() < aParticle->GetStrangeness()) {
		break; }
	    }
	  }
	}
	tPLIter++;
      }
    mParticles.insert(tPLIter, *aParticle);
  }
  else
    mParticles.push_back(*aParticle);
}

void           
Event::ReadMultiplicities()
{
  char   tName[200];
  double tMult;
  
  char *tHash;
  int   multsize;
  char *tMultName = 0;
  ifstream *fin = NULL;
  
  tHash = mInteg->ParameterHash();
  multsize = strlen(mFOHSlocation.Data()) + 25 + strlen(tHash);
  tMultName = (char *) malloc(sizeof(char) * multsize);

  if (!tMultName) {
    printf("Cannot allocate memory!\n");
    exit(0);
  }

  if (mFOHSlocation != "") {
    strcpy(tMultName, mFOHSlocation.Data());
    strcat(tMultName, "/");
    strcat(tMultName, "fmultiplicity_");
    strcat(tMultName, tHash);
    strcat(tMultName, ".txt");
    fin = new ifstream(tMultName);
  }
  else if (!((fin) && (fin->is_open()))) {
    strcpy(tMultName, "fmultiplicity_");
    strcat(tMultName, tHash);
    strcat(tMultName, ".txt");
    fin = new ifstream(tMultName);
  }

//   strcpy(tMultName, "fmultiplicity_");
//   strcat(tMultName, tHash);
//   strcat(tMultName, ".txt");

//  fin = new ifstream(tMultName);
  if ((fin) && (fin->is_open())) {
    PRINT_MESSAGE("Reading Multiplicities from " << tMultName);

    while (!fin->eof())
    {
      (*fin) >> tName >> tMult;
      PRINT_DEBUG_2(tName << " " <<  mDB->GetParticleTypeIndex(tName) << " " << tMult);
      mAverageMultiplicities[mDB->GetParticleTypeIndex(tName)] = tMult;
    }
    fin->close();
  }

  delete tHash;
  delete tMultName;
}

void           
Event::GenerateMultiplicities()
{
  for (unsigned int tIter=0; tIter<mAverageMultiplicities.size(); tIter++)
    mMultiplicities[tIter] = mRandom->Poisson(mAverageMultiplicities[tIter]);
}

void 
Event::DecayParticles()
{
  Particle *tPart1, *tPart2, *tPart3, *tFather;
  ParticleDecayer* tDecayer;
  int tCount = 0;
  
  tDecayer = new ParticleDecayer(mDB);
  tDecayer->SeedSet(mRandom->Integer(100000000));
  
  //  InitParticleScan();
  while ((tFather = GetParticleOfCount(tCount))) 
    {
      if (tFather->GetParticleType()->GetGamma() >= 0.0) {
	if ((tFather->GetParticleType()->GetTable()) && (((DecayTable *) tFather->GetParticleType()->GetTable())->GetChannelCount()+1 > 0))
	  {
	    tDecayer->DecayParticle(tFather, &tPart1, &tPart2, &tPart3);
#ifndef _RESCALE_CHANNELS_
	    if (!tPart1) {
	      tCount++;
	      continue;
	    }
#endif
#ifdef _NO_THREE_BODY_DECAYS_
	    if (tPart3){
	      tCount++;
	      if (tPart1) delete tPart1;
	      if (tPart2) delete tPart2;
	      continue;
	    }
#endif
#ifdef _OMIT_TWO_BODY_
	    if (tPart1 && tPart2 && !tPart3){
	      tCount++;
	      delete tPart1;
	      delete tPart2;
	      continue;
	    }
#endif

	    tPart1->SetFather(tCount);
	    tPart2->SetFather(tCount);
	    if (tPart3)
	      tPart3->SetFather(tCount);
	    
	    AddParticle(tPart1);
	    AddParticle(tPart2);
	    
	    delete tPart1;
	    delete tPart2;

	    if (tPart3) {
	      AddParticle(tPart3);
	      delete tPart3;
	    }
	  }
	else
	  {
	  }
      }
      tCount++;
    }

  delete tDecayer;
}

void 
Event::Randomize()
{
  TDatime dat;
  
#ifdef _ROOT_4_
  mRandom->SetSeed2(dat.Get() / 2 * 3, dat.Get() / 11 * 9);
#else
  mRandom->SetSeed(dat.Get() / 2 * 3);
#endif
}

void           
Event::ReadParameters()
{
  STR tEvFile;
  STR tUseNegativeBinomial;

  try {
    tEvFile = sRPInstance->getPar("EventOutputFile");
  }
  catch (STR e) {
    PRINT_DEBUG_1("Event::ReadParameters - Caught exception " << e);
    PRINT_MESSAGE("Did not find one of the neccessary parameters in the parameters file.");
    exit(0);
  }
  try {
    tUseNegativeBinomial = sRPInstance->getPar("MultiplicityDistribution");
    if (tUseNegativeBinomial.Contains("NegativeBinomial"))
      mNegBin = 1;
    else
      mNegBin = 0;
  }
  catch (STR e) {
    PRINT_MESSAGE("Using default multiplicty distribution: Poissonian");
    mNegBin = 0;
  }
  try {
    mFOHSlocation = sRPInstance->getPar("FOHSLocation");
  }
  catch (STR e) {
  }
	
}

