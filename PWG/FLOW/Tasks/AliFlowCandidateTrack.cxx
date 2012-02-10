/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

////////////////////////////////////////////////////
// AliFlowCandidateTrack:
// Class for reconstructed particles to be used in flow analysis
// Author: Carlos Perez (cperez@cern.ch)
////////////////////////////////////////////////////

#include "AliFlowCandidateTrack.h"

ClassImp(AliFlowCandidateTrack)

AliFlowCandidateTrack::AliFlowCandidateTrack():
    AliFlowTrack(),
    fNDaughters(0)
{
  // ctor
  for(int i=0; i!=5; ++i) {
    fDaughter[i] = -1;
    fTrack[i] = NULL;
  }
}

AliFlowCandidateTrack::AliFlowCandidateTrack(const AliFlowCandidateTrack& aTrack):
  AliFlowTrack(aTrack),
  fNDaughters(aTrack.fNDaughters)
{
  // ctor
  for(int i=0; i!=5; ++i) {
    fDaughter[i] = aTrack.fDaughter[i];
    fTrack[i] = aTrack.fTrack[i];
  }
}

AliFlowCandidateTrack&  AliFlowCandidateTrack::operator=(const AliFlowCandidateTrack& aTrack)
{
  // assignment
  if (this == &aTrack) return *this; //handles self assignment

  AliFlowTrack::operator=(aTrack);
  fNDaughters = aTrack.fNDaughters;
  for(int i=0; i!=5; ++i) {
    fDaughter[i] = aTrack.fDaughter[i];
    fTrack[i] = aTrack.fTrack[i];
  }
  return *this;
}

AliFlowCandidateTrack::~AliFlowCandidateTrack()
{
  // dtor
}

