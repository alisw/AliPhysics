/// \class AliTrackContainerToyModel
/// \brief Allows to modify the tracks to implement toy models
///
/// This class derives from AliTrackContainer.
/// It allows to select tracks based and modify their momenta according to some
/// toy modeling.
/// At the moment only pt scaling is implemented
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \author Leticia Conqueiro Mendez <leticia.cunqueiro.mendez@cern.ch>
/// \date Jul 22, 2016

#ifndef ALITRACKCONTAINERTOYMODEL_H
#define ALITRACKCONTAINERTOYMODEL_H

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

class AliAODEvent;

#include "AliAODEvent.h"
#include "AliTrackContainer.h"


class AliTrackContainerToyModel : public AliTrackContainer {
public:
  AliTrackContainerToyModel();
  AliTrackContainerToyModel(const char *name);

  void SetTrackScalePt(Double_t t) { fTrackScalePt  = t; }

  virtual Bool_t GetMomentumFromTrack(TLorentzVector &mom, const AliVTrack* track, Double_t mass) const;
  virtual Bool_t GetMomentum(TLorentzVector &mom, Int_t i) const;
  virtual Bool_t GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;
  virtual Bool_t GetNextMomentum(TLorentzVector &mom);
  virtual Bool_t GetNextAcceptMomentum(TLorentzVector &mom);

  void ScalePtOfLorentzVector(TLorentzVector &mom) const;

protected:

  Double_t               fTrackScalePt;           //scaling of the track pT by given fraction (0....1)

private:
  /// \cond CLASSIMP
  ClassDef(AliTrackContainerToyModel, 1);
  /// \endcond
};
#endif
