/**
 * \file AliReducedPatchContainer.h
 * \brief Declaration of class AliReducedPatchContainer, a container for reduced trigger patches
 *
 * In this file, the class AliReducedPatchContainer, a container for reduced trigger patches, is defined
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDPATCHCONTAINER_H
#define ALIREDUCEDPATCHCONTAINER_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

class TObjArray;

/**
 * \namespace HighPtTracks
 * \brief Namespace for classes creating trees of events with jets
 *
 * This namespace contains classes descibing reduced events with high
 * jets. A jet event consists of the following classes.
 *    - AliReducedJetEvent
 *    - AliReducedJetParticle
 *    - AliReducedJetConstituent
 *    - AliReducedMatchedTrack
 *  The task AliHighPtReconstructionEfficiency produces the reduced jet events
 *  and is included in the namespace as well. Also the helper class AliParticleMap
 *  is part of this namespace.
 */
namespace HighPtTracks {

/**
 * \class AliReducedPatchContainer
 * \brief Container structure for reduced trigger patches
 *
 * This class stores reduced trigger patch information, depending on whether they are online or offline patches,
 * and depending on the trigger class
 */
class AliReducedPatchContainer: public TObject {
public:
  /**
   * \enum Definition of jet trigger patches
   */
  enum PatchType_t{
    kEMCGammaLow = 0,         ///< Gamma trigger, low threshold
    kEMCGammaHigh = 1,        ///< Gamma trigger, high threshold
    kEMCJetLow = 2,           ///< Jet trigger, low threshold
    kEMCJetHigh = 3           ///< Jet trigger, high threshold
  };
  AliReducedPatchContainer(Bool_t doAlloc = kFALSE);
  AliReducedPatchContainer(const AliReducedPatchContainer &cont);
  AliReducedPatchContainer &operator=(const AliReducedPatchContainer &cont);
  void Copy(TObject &target) const;
  virtual ~AliReducedPatchContainer();

  void AddTriggerPatch(Bool_t isOffline, PatchType_t patchtype, Float_t energy, Int_t amplitude, Float_t eta, Float_t phi);
  TObjArray *GetTriggerPatches(Bool_t isOffline, PatchType_t patchtype);

protected:
  TObjArray                   *fOnlinePatches[4];             ///< Trigger type dependent container for online patches
  TObjArray                   *fOfflinePatches[4];            ///< Trigger type dependent container for offline patches

  /// \cond CLASSIMP
  ClassDef(AliReducedPatchContainer, 1);;
  /// \endcond CLASSIMP
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDPATCHCONTAINER_H */
