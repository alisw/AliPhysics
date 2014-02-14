/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisNonMuonTrackCuts.h"
#include "AliAODTrack.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisNonMuonTrackCuts)
/// \endcond

/////////////////////////////////////////////////////////////////////////

AliAnalysisNonMuonTrackCuts::AliAnalysisNonMuonTrackCuts()
{
  /// default ctor
}

Bool_t AliAnalysisNonMuonTrackCuts::IsSelected(TObject* obj)
{
  /// Returns true if the object is a muon track
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);
  //  if (track && track->IsMuonTrack()) return kTRUE;
  if (track && (track->IsMuonTrack() || track->IsMuonGlobalTrack())) return kTRUE;
  return kFALSE;
}
