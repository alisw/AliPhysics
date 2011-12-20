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
#include "DecayChannel.h"

DecayChannel::DecayChannel() :
  mParticleType1(0), mParticleType2(0), mParticleType3(-1), mBranchRatio(0.0)
{
}

DecayChannel::DecayChannel(const DecayChannel& aChannel) 
{
  mBranchRatio = aChannel.GetBranchingRatio();
  mParticleType1 = aChannel.GetParticle1();
  mParticleType2 = aChannel.GetParticle2();
  mParticleType3 = aChannel.GetParticle3();
}

DecayChannel::DecayChannel(double aBranchRatio, int aPartType1, int aPartType2, int aPartType3) :
  mParticleType1(aPartType1), mParticleType2(aPartType2), mParticleType3(aPartType3), mBranchRatio(aBranchRatio)
{
}

DecayChannel::~DecayChannel()
{
}

DecayChannel& DecayChannel::operator=(const DecayChannel& aChannel)
{
  if (this != &aChannel) {
    mBranchRatio = aChannel.GetBranchingRatio();
    mParticleType1 = aChannel.GetParticle1();
    mParticleType2 = aChannel.GetParticle2();
    mParticleType3 = aChannel.GetParticle3();
  }

  return *this;
}

int    
DecayChannel::GetParticle1() const
{
  return mParticleType1;
}

int    
DecayChannel::GetParticle2() const
{
  return mParticleType2;
}

int    
DecayChannel::GetParticle3() const
{
  return mParticleType3;
}

double 
DecayChannel::GetBranchingRatio() const
{
  return mBranchRatio;
}

int    
DecayChannel::Is3Particle() const
{
  return (mParticleType3 != -1);
}


void   DecayChannel::SetParticle1(int aPartType1)
{
  mParticleType1 = aPartType1;
}

void   DecayChannel::SetParticle2(int aPartType2)
{
  mParticleType2 = aPartType2;
}

void   
DecayChannel::SetParticle3(int aPartType3)
{
  mParticleType3 = aPartType3;
}

void   DecayChannel::SetBranchingRatio(double aRatio)
{
  mBranchRatio = aRatio;
}
