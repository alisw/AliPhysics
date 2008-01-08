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

#ifndef _BFPW_DECAY_CHANNEL_
#define _BFPW_DECAY_CHANNEL_
#include "ParticleType.h"

class DecayChannel
{
 public:
  DecayChannel();
  DecayChannel(double aBranchRatio, int aPartType1, int aPartType2, int aPartType3 = -1);
  DecayChannel(const DecayChannel& aChannel);
  ~DecayChannel();
  
  int    GetParticle1() const;
  int    GetParticle2() const;
  int    GetParticle3() const;
  double GetBranchingRatio() const;
  int    Is3Particle() const;
  
  void   SetParticle1(int aPartType1);
  void   SetParticle2(int aPartType2);
  void   SetParticle3(int aPartType3);
  void   SetBranchingRatio(double aRatio);
  
 protected:
  int    mParticleType1;
  int    mParticleType2;
  int    mParticleType3;
  double mBranchRatio;
  
};

#endif
