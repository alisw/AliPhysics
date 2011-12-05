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
#include "ParticleType.h"
#include "DecayTable.h"

ParticleType::ParticleType()
{
  mName = "";
  mNumber=0;
  mMass=-1.;
  mStrangeness=-1;
  mBarionN=-1;
  mCharmN=-1;
  mSpin=-1.;
  mI=-1.;
  mI3=-1.;
  mGamma=-1.;
  mDecayChannelCount2=0;
  mDecayChannelCount3=0;
  mTable = new DecayTable();
  mPDGCode = 0;
  mFMax = 0.0;
}

ParticleType::ParticleType(const ParticleType& aParticleType)
{
  mName = aParticleType.GetName();
  mNumber = aParticleType.GetNumber();
  mMass = aParticleType.GetMass();
  mStrangeness = aParticleType.GetStrangeness();
  mBarionN=aParticleType.GetBarionN();
  mCharmN=aParticleType.GetCharmN();
  mSpin=aParticleType.GetSpin();
  mI=aParticleType.GetI();
  mI3=aParticleType.GetI3();
  mGamma=aParticleType.GetGamma();
  mDecayChannelCount2=aParticleType.GetDecayChannelCount2();
  mDecayChannelCount3=aParticleType.GetDecayChannelCount3();
  mPDGCode = aParticleType.GetPDGCode();
  mFMax = aParticleType.GetFMax();
  mTable = new DecayTable(*(aParticleType.GetTable()));
}

ParticleType& ParticleType::operator=(const ParticleType& aParticleType)
{
  if (this != &aParticleType) {
    mName = aParticleType.GetName();
    mNumber = aParticleType.GetNumber();
    mMass = aParticleType.GetMass();
    mStrangeness = aParticleType.GetStrangeness();
    mBarionN=aParticleType.GetBarionN();
    mCharmN=aParticleType.GetCharmN();
    mSpin=aParticleType.GetSpin();
    mI=aParticleType.GetI();
    mI3=aParticleType.GetI3();
    mGamma=aParticleType.GetGamma();
    mDecayChannelCount2=aParticleType.GetDecayChannelCount2();
    mDecayChannelCount3=aParticleType.GetDecayChannelCount3();
    mPDGCode = aParticleType.GetPDGCode();
    mFMax = aParticleType.GetFMax();
    delete mTable;
    mTable = new DecayTable(*(aParticleType.GetTable()));
  }

  return *this;
}


ParticleType::~ParticleType()
{
  if (mTable)
    delete mTable;
}

void ParticleType::WriteParticle(int /*i*/)
{
    mMass=-1;
    mStrangeness=-1;
    mBarionN=-1;
    mCharmN=-1;
    mSpin=-1;
    mName="1";
}

int    ParticleType::GetNumber() const { return mNumber; }
float  ParticleType::GetMass() const { return mMass; }
int    ParticleType::GetStrangeness() const { return mStrangeness; }
int    ParticleType::GetBarionN() const { return mBarionN; }
int    ParticleType::GetCharmN() const { return mCharmN; }
float  ParticleType::GetSpin() const { return mSpin; }
float  ParticleType::GetI() const { return mI; }
float  ParticleType::GetI3() const { return mI3; }
float  ParticleType::GetGamma() const { return mGamma; }
std::string ParticleType::GetName() const { return mName; }
int    ParticleType::GetDecayChannelCount2() const { return mDecayChannelCount2; }
int    ParticleType::GetDecayChannelCount3() const { return mDecayChannelCount3; }
int    ParticleType::GetCharge() { return int(mI3+(mBarionN+mStrangeness)/2.); }	/*MCH int() added */
int    ParticleType::GetPDGCode() const { return mPDGCode; }

  

void ParticleType::SetNumber(int aNumber) { mNumber = aNumber; }
void ParticleType::SetMass(float aMass) { mMass = aMass; }
void ParticleType::SetStrangeness(int aStrangeness) { mStrangeness = aStrangeness; }
void ParticleType::SetBarionN(int aBarionN) { mBarionN = aBarionN; }
void ParticleType::SetCharmN(int aCharmN) { mCharmN = aCharmN; }
void ParticleType::SetSpin(float aSpin) { mSpin = aSpin; }
void ParticleType::SetI(float aI) { mI = aI; }
void ParticleType::SetI3(float aI3) { mI3 = aI3; }
void ParticleType::SetGamma(float aGamma) { mGamma = aGamma; }
void ParticleType::SetName(char *aName) { mName = aName; }
void ParticleType::SetDecayChannelCount2(int aDCCount2) { mDecayChannelCount2 = aDCCount2; }
void ParticleType::SetDecayChannelCount3(int aDCCount3) { mDecayChannelCount3 = aDCCount3; }
void ParticleType::SetPDGCode(int aCode) { mPDGCode = aCode; }

DecayTable* 
ParticleType::GetTable() const
{
  if (mTable)
    return mTable;
  else
    return NULL;
}

void        
ParticleType::AddDecayChannel(DecayChannel aChannel)
{
  if (!mTable)
    mTable = new DecayTable();
  mTable->AddDecayChannel(aChannel);
}

float  
ParticleType::GetFMax() const
{
  return mFMax;
}

void 
ParticleType::SetFMax(float aFMax)
{
  mFMax = aFMax;
}
