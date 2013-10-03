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

#include "AliAnalysisNonPrimaryVertices.h"
#include "AliAODVertex.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisNonPrimaryVertices)
/// \endcond

AliAnalysisNonPrimaryVertices::AliAnalysisNonPrimaryVertices()
{
  /// default ctor
}

Bool_t AliAnalysisNonPrimaryVertices::IsSelected(TObject* obj)
{
  /// Returns true if the object is a primary or pileup vertex
  
  AliAODVertex* vertex = dynamic_cast<AliAODVertex*>(obj);
  if (vertex)
  {
    if ( vertex->GetType() == AliAODVertex::kPrimary ||
        vertex->GetType() == AliAODVertex::kMainSPD ||
        vertex->GetType() == AliAODVertex::kPileupSPD ||
        vertex->GetType() == AliAODVertex::kPileupTracks ||
        vertex->GetType() == AliAODVertex::kMainTPC )
    {
      return kTRUE;
    }
  }
  
  //  enum AODVtx_t {kUndef=-1, kPrimary, kKink, kV0, kCascade, kMulti, kMainSPD, kPileupSPD, kPileupTracks,kMainTPC};
  
  return kFALSE;
}

