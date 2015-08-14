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

  void FillMomentumVector(TVector3 &pvec) const;
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

  Double_t GetPx() const { return fPVec[0]; }
  Double_t GetPy() const { return fPVec[1]; }
  Double_t GetPz() const { return fPVec[2]; }

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
   * Get the number of TPC clusters
   * \return Number of TPC clustes
   */
  Int_t GetNumberOfTPCClusters() const { return fClustersTPC; }
  /**
   * Get the number of TPC crossed rows
   * \return Number of TPC crossed rows
   */
  Int_t GetNumberOfTPCCrossedRows() const { return fCrossedRowsTPC; }
  /**
   * Get the number of shared clusters in the TPC
   * \return Number of shared clusters
   */
  Int_t GetNumberOfTPCSharedClusters() const { return fSharedClustersTPC; }
  /**
   * Get the number of findable clusters in the TPC
   * \return Number of findable clusters
   */
  Int_t GetNumberOfTPCFindableClusters() const { return fFindableClustersTPC; }

  Double_t Pt() const;
  Double_t Eta() const;
  Double_t Phi() const;

  /**
   * Access to components of the  3-momentum vector.
   * \param px x-component
   * \param py y-component
   * \param pz z-component
   */
  void GetMomentumVector(Double_t &px, Double_t &py, Double_t &pz) const {
    px = fPVec[0];
    py = fPVec[1];
    pz = fPVec[2];
  }

  /**
   * Get the charge of the reconstructed track
   * \return Charge of the reconstructed track
   */
  Char_t Charge() const { return fCharge; }

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
  /**
   * Set the number of clusters in the TPC
   * \param nclusters Number of clusters
   */
  void SetTPCClusters(Int_t nclusters) { fClustersTPC = nclusters; }
  /**
   * Set the number of crossed rows in the TPC
   * \param crossedRows Number of crossed rows
   */
  void SetTPCCrossedRows(Int_t crossedRows) { fCrossedRowsTPC = crossedRows; }
  /**
   * Set Number of shared clusters in the TPC
   * @param nshared Number of shared clusters
   */
  void SetTPCSharedClusters(Int_t nshared) { fSharedClustersTPC = nshared; }
  /**
   * Set number of findable clusters in the TPC
   * \param nfindable Number of findable clusters
   */
  void SetTPCFindableClusters(Int_t nfindable) { fFindableClustersTPC = nfindable; }
  /**
   * Set the charge of the reconstructed track
   * \param charge Charge of the reconstructed track
   */
  void SetCharge(Char_t charge) { fCharge = charge; }


protected:
  Double_t                      fPVec[3];             ///< Momentum vector of the particle
  Char_t                        fCharge;              ///< Charge
  Long_t                        fTrackCutsMap;        ///< Map of Track cuts
  Int_t                         fClusterIndex;        ///< Index of cluster matched to the track
  Int_t                         fParticleIndex;       ///< Index of true particle matched to this track
  Bool_t                        fGoodMCTrack;         ///< Mark track as good (according to the track label)
  Char_t                        fClustersTPC;         ///< Number of TPC clusters
  Char_t                        fCrossedRowsTPC;      ///< Number of TPC crossed rows
  Char_t                        fSharedClustersTPC;   ///< Number of shared clusters in the TPC
  Char_t                        fFindableClustersTPC; ///< Number of findable clusters in the TPC

  /// \cond CLASSIMP
  ClassDef(AliReducedReconstructedTrack, 2);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDRECONSTRUCTEDTRACK_H */
