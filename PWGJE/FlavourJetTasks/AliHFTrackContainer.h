/// \class AliHFTrackContainer
/// \brief Select tracks based on specific prescriptions of HF analysis
///
/// This class derives from AliParticleContainer. It allows
/// to select tracks based on specific prescriptions of HF analysis.
/// In particular it will reject tracks that are daughters of a
/// specified D meson candidate
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb 10, 2016

#ifndef ALIHFTRACKCONTAINER_H
#define ALIHFTRACKCONTAINER_H

/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

class AliAODRecoDecay;
class AliAODTrack;

#include "AliParticleContainer.h"
#include "AliAnalysisTaskDmesonJets.h"

class AliHFTrackContainer : public AliParticleContainer {
 public:
  AliHFTrackContainer();
  AliHFTrackContainer(const char *name);

  void                 SetDMesonCandidate(AliAODRecoDecay* c);

  Bool_t               AcceptParticle(AliVParticle* vp);
  Bool_t               AcceptParticle(Int_t i);
  void                 GenerateDaughterList();
  const TObjArray&     GetDaughterList() const                         { return fDaughterList            ; }
  
 protected:
  void                 AddDaughters(AliAODRecoDecay* cand);
  Bool_t               IsDMesonDaughter(AliAODTrack* track);
 
  AliAODRecoDecay*     fDMesonCandidate;             ///<  Exclude daughters of this D meson candidate
  TObjArray            fDaughterList;                ///<  Daughters of the D meson candidate
  
 private:
  AliHFTrackContainer(const AliHFTrackContainer&);            // not implemented
  AliHFTrackContainer &operator=(const AliHFTrackContainer&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliHFTrackContainer, 1);
  /// \endcond
};
#endif
