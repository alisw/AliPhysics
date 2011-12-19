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
#include "DecayTable.h"

using namespace std;

DecayTable::DecayTable()
{
  mDecayChannels.clear();
  mBranchingRatios.clear();
}

DecayTable::DecayTable(const DecayTable& aTable)
{
  
  mDecayChannels.clear();
  mBranchingRatios.clear();
  for (int iter=0; iter<aTable.GetChannelCount(); iter++)
    AddDecayChannel(*(aTable.GetDecayChannel(iter)));
}

DecayTable::~DecayTable()
{
}

DecayTable& DecayTable::operator=(const DecayTable& aTable)
{
  if (this != &aTable) {
    mDecayChannels.clear();
    mBranchingRatios.clear();
    for (int iter=0; iter<aTable.GetChannelCount(); iter++)
      AddDecayChannel(*(aTable.GetDecayChannel(iter)));
  }
  
  return *this;
}

void 
DecayTable::AddDecayChannel(DecayChannel aChannel)
{
  mDecayChannels.push_back(aChannel);
  RecalculateBranchingRatios();
}

void 
DecayTable::RecalculateBranchingRatios()
{
  float tSumRatio = 0.0;
  float tCurRatio = 0.0;
  
  for (unsigned int tIter=0; tIter<mDecayChannels.size(); tIter++)
    tSumRatio += mDecayChannels[tIter].GetBranchingRatio();
  
  for (unsigned int tIter=0; tIter<mDecayChannels.size(); tIter++)
    {
      tCurRatio += mDecayChannels[tIter].GetBranchingRatio()/tSumRatio;
      if (mBranchingRatios.size() <= tIter)
	mBranchingRatios.push_back(tCurRatio);
      else
	mBranchingRatios[tIter] = tCurRatio;
    }
}

int  
DecayTable::GetChannelCount() const
{
  return mDecayChannels.size()-1;
}


const DecayChannel* 
DecayTable::GetDecayChannel(int aIndex) const
{
  return &(mDecayChannels[aIndex]);
}

float 
DecayTable::GetDecayStep(int aIndex)
{
  return mBranchingRatios[aIndex];
}

int   
DecayTable::ChooseDecayChannel(double aProb)
{
  unsigned int tChanIndex = 0;
  while ((mBranchingRatios[tChanIndex] < aProb) && (tChanIndex < mDecayChannels.size()))
    tChanIndex++;
  
  return tChanIndex;
}

int   
DecayTable::ChooseDecayChannelOrNot(double aProb)
{
  float tSumRatio = 0.0;
  
  for (unsigned int tIter=0; tIter<mDecayChannels.size(); tIter++) {
    if ((aProb > tSumRatio) && (aProb <= tSumRatio+mDecayChannels[tIter].GetBranchingRatio()))
      return tIter;
    tSumRatio += mDecayChannels[tIter].GetBranchingRatio();
  }

  return -1;
}

