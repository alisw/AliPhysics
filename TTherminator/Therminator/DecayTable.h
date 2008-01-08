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
#ifndef _BFPW_DECAY_TABLE_
#define _BFPW_DECAY_TABLE_

#include <vector>
#include "THGlobal.h"
#include "DecayChannel.h"

class DecayTable
{
 public:
  DecayTable();
  DecayTable(const DecayTable& aTable);
  ~DecayTable();
  
  void AddDecayChannel(DecayChannel aChannel);

  int  GetChannelCount() const;
  const DecayChannel* GetDecayChannel(int aIndex) const;
  float GetDecayStep(int aIndex);
  int   ChooseDecayChannel(double aProb);
  int   ChooseDecayChannelOrNot(double aProb);
  
 private:
  std::vector<DecayChannel> mDecayChannels;
  std::vector<float>        mBranchingRatios;
  
  void RecalculateBranchingRatios();
};



#endif
