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
#ifndef _BFPW_EVENT_
#define _BFPW_EVENT_

#include <vector>
#include <list>
#include <TRandom2.h>
#include "THGlobal.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "Integrator.h"

using namespace std;

typedef list<Particle> ParticleList;
typedef list<Particle>::iterator ParticleListIter;

class Particle;

class Event
{
 public:
  Event(ParticleDB *aDB, Integrator *aInteg);
  ~Event();

  void GenerateEvent(int aSeed=0);
  void DecayParticles();
  
  int  GetParticleCount();
  void InitParticleScan();
  Particle* GetNextParticle();
  Particle* GetParticleOfCount(int aCount);
  void WriteEvent(int nEvent=0);
  void AddParticle(Particle* aParticle);
  void Randomize();
  
 private:
  ParticleList   mParticles;
  vector<int>    mMultiplicities;
  vector<double> mAverageMultiplicities;
  void           GenerateMultiplicities();
  void           ReadMultiplicities();
  void           ReadParameters();
  ParticleDB    *mDB;
  Integrator    *mInteg;
  TRandom2      *mRandom;
  ParticleListIter mCurIter;
  Int_t          mNegBin;
  Int_t          mScanStarted;
  TString        mFOHSlocation;
};


#endif
