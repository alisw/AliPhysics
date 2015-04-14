/**
 * \file AliReducedReconstuctedTrack.h
 * \brief Declaration of class AliReducedReconstuctedTrack
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDRECONSTRUCTEDTRACK_H
#define ALIREDUCEDRECONSTRUCTEDTRACK_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TArrayI;
class TList;
class TVector3;

/**
 * \namespace HighPtTracks
 * \brief Namespace for classes creating trees of events with jets
 *
 * This namespace contains classes describing reduced events with high
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
 * \class AliReducedReconstructedTrack
 * \brief Structure for reconstructed track information
 */
class AliReducedReconstructedTrack: public TObject {
public:
  AliReducedReconstructedTrack();
  virtual ~AliReducedReconstructedTrack();

  void FillMomentumVector(TVector3 &pvec);
  /**
   * Set the Momentum vector of the track at the vertex
   * \param px x-component of the primary vertex
   * \param py y-component of the primary vertex
   * \param pz z-component of the primary vertex
   */
  void SetMomentumVector(Double_t px, Double_t py, Double_t pz){
    fPVec[0] = px;
    fPVec[1] = py;
    fPVec[2] = pz;
  }
  /**
   * Check whether track was selected by a certain track cut combination during filtering
   * \param index Index of the track cuts in the map
   * \return True if the track was selected
   */
  Bool_t TestTrackCuts(Int_t index) { return TESTBIT(fTrackCutsMap, index); }
  /**
   * Get the index of the matched cluster
   * \return Index of the cluster
   */
  Int_t GetClusterIndex() const { return fClusterIndex; }
  /**
   * Get the index of the matched particle
   * \return Index of the particle
   */
  Int_t GetMatchedParticleIndex() const { return fParticleIndex; }
  /**
   * Check if track is a good track according to the label
   * \return True if track is good, false otherwise
   */
  Bool_t IsGoodTrack() const {return fGoodMCTrack; }

  /**
   * Mark track as selected by a given track cuts using the index
   * \param index Index of the track cuts
   */
  void SetTrackCuts(Int_t index) { SETBIT(fTrackCutsMap, index); }
  /**
   * Set index of matched cluster
   * \param index Index of the cluster
   */
  void SetMatchedClusterIndex(Int_t index) { fClusterIndex = index; }
  /**
   * Set index of the matched particle
   * \param index Index of the particle
   */
  void SetMatchedParticleIndex(Int_t index) { fParticleIndex = index; }

  /**
   * Mark track as MC-good (according to the label)
   * \param isGood Boolean for good tracks
   */
  void SetGoodTrackLabel(Bool_t isGood) { fGoodMCTrack = isGood; }


protected:
  Double_t                      fPVec[3];             ///< Momentum vector of the particle
  Long_t                        fTrackCutsMap;        ///< Map of Track cuts
  Int_t                         fClusterIndex;        ///< Index of cluster matched to the track
  Int_t                         fParticleIndex;       ///< Index of true particle matched to this track
  Bool_t                        fGoodMCTrack;         ///< Mark track as good (according to the track label)

  /// \cond CLASSIMP
  ClassDef(AliReducedReconstructedTrack, 1);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDRECONSTRUCTEDTRACK_H */
