/**
 * \file AliReducedHighPtEvent.h
 * \brief Reduced event structure for high-\f$ p_{t} \f$ analysis
 *
 * In this file a simple reduced event structure containing
 * - event information
 * - reconstructed EMCAL trigger information
 * - reconstructed EMCAL cluster
 * - generated particles (if created from Monte-Carlo events)
 * - Monte-Carlo event information (if created from Monte-Carlo events)
 * - reconstructed tracks
 * is created
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDHIGHPTEVENT_H
#define ALIREDUCEDHIGHPTEVENT_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#include <vector>

#include <TObject.h>
#include "AliReducedPatchContainer.h"

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

class AliReducedEmcalCluster;
class AliReducedGeneratedParticle;
class AliReducedReconstructedTrack;
class AliReducedMCHeader;

/**
 * \class AliReducedHighPtEvent
 * \brief Event structure for high-pt analysis
 */
class AliReducedHighPtEvent : public TObject {
public:
  AliReducedHighPtEvent(Bool_t doAlloc = kFALSE);
  AliReducedHighPtEvent(const AliReducedHighPtEvent &ref);
  AliReducedHighPtEvent &operator=(const AliReducedHighPtEvent &ref);
  virtual ~AliReducedHighPtEvent();
  void Copy(TObject &target) const;

  /**
   * Get the cluster container
   * \return The cluster container
   */
  TObjArray *GetClusterContainer() { return fReducedClusterInfo; }
  std::vector<HighPtTracks::AliReducedEmcalCluster *> GetClusterVector() const;
  /**
   * Get the particle container (at generator level)
   * \return The particle container
   */
  TObjArray *GetParticleContainer() { return fReducedParticleInfo; }
  std::vector<HighPtTracks::AliReducedGeneratedParticle *> GetParticleVector() const;
  /**
   * Get the container with reconstructed tracks
   * \return Container with reconstructed tracks
   */
  TObjArray *GetTrackContainer() { return fReducedTrackInfo; }
  std::vector<HighPtTracks::AliReducedReconstructedTrack *> GetTrackVector() const;
  AliReducedEmcalCluster * GetClusterForIndex(Int_t index);
  AliReducedGeneratedParticle *GetParticleForIndex(Int_t index);
  /**
   * Get the container with trigger patches
   * \return container with trigger patches
   */
  AliReducedPatchContainer *GetPatchContainer() { return fReducedPatchInfo; }
  /**
   * Get the centrality percentile of the event
   * \return centrality percentile of the event
   */
  Float_t GetCentralityPercentile() const { return fCentralityPercentile; }
  /**
   * Get the z-position of the primary vertex
   * \return z-position of the primary vertex
   */
  Float_t GetVertexZ() const { return fVertexZ; }
  /**
   * Get the Monte-Carlo header
   * \return The Monte-Carlo header (NULL if not set)
   */
  AliReducedMCHeader *GetMonteCarloHeader() { return fMCHeader; }
  /**
   * Check if event is a min. bias event
   * \return true if event is a min. bias event, false otherwise
   */
  Bool_t IsMinBias() const { return fIsMinBias; }
  /**
   * Check if event is triggered as Gamma low
   * \return true if event is a gamma low event, false otherwise
   */
  Bool_t IsGammaLowFromString() const { return fGammaTriggerString[0]; }
  /**
   * Check if event is triggered as Gamma high
   * \return true if event is a gamma high event, false otherwise
   */
  Bool_t IsGammaHighFromString() const { return fGammaTriggerString[1]; }
  /**
   * Check if event is triggered as Jet low
   * \return true if event is a jet low event, false otherwise
   */
  Bool_t IsJetLowFromString() const { return fJetTriggerString[0]; }
  /**
   * Check if event is triggered as Jet high
   * \return true if event is a jet high event, false otherwise
   */
  Bool_t IsJetHighFromString() const { return fJetTriggerString[1]; }
  /**
   * Get the run number
   * \return Run number
   */
  Int_t GetRunNumber() const { return fRunNumber; }

  void AddReducedCluster(AliReducedEmcalCluster *cluster);
  void AddReducedGeneratedParticle(AliReducedGeneratedParticle *part);
  void AddReducedReconstructedParticle(AliReducedReconstructedTrack *trk);
  /**
   * Set the z-Position of the primary vertex
   * \param vz z-Position of the primary vertex
   */
  void SetVertexZ(Float_t vz)                   { fVertexZ = vz; }
  /**
   * Set the event centrality percentile
   * \param cent Event centrality percentile
   */
  void SetCentralityPercentile(Float_t cent)    { fCentralityPercentile = cent; }
  /**
   * Set the trigger decision from the trigger string
   * \param isGammaLow If true event is triggered as gamma low event
   * \param isGammaHigh If true event is triggered as gamma high event
   * \param isJetLow If true event is triggered as jet low event
   * \param isJetHigh If true event is triggered as jet high event
   */
  void SetDecisionFromTriggerString(Bool_t isGammaLow, Bool_t isGammaHigh, Bool_t isJetLow, Bool_t isJetHigh) {
    fGammaTriggerString[0] = isGammaLow;
    fGammaTriggerString[1] = isGammaHigh;
    fJetTriggerString[0] = isJetLow;
    fJetTriggerString[1] = isJetHigh;
  }
  /**
   * Flag event as min. bias event
   * \param isMinBias If true evnet is a min. bias event
   */
  void SetMinBiasEvent(Bool_t isMinBias) { fIsMinBias = isMinBias; }
  /**
   * Set the reduced MC event header
   * \param header The Monte-Carlo event header
   */
  void SetMonteCarloHeader(AliReducedMCHeader *header) { fMCHeader = header; }
  /**
   * Set the run number
   * \param runnumber The run number
   */
  void SetRunNumber(Int_t runnumber) { fRunNumber = runnumber; }

protected:
  Int_t                                   fRunNumber;                         ///< Run number
  Float_t                                 fCentralityPercentile;              ///< Centrality percentile
  Float_t                                 fVertexZ;                           ///< z-position of the primary vertex
  AliReducedMCHeader                      *fMCHeader;                         ///< Reduced Monte-Carlo header
  Bool_t                                  fJetTriggerString[2];               ///< jet trigger selection from trigger string
  Bool_t                                  fGammaTriggerString[2];             ///< gamma trigger selection from trigger string
  Bool_t                                  fIsMinBias;                         ///< Flag event as min. bias event
  AliReducedPatchContainer                *fReducedPatchInfo;                 ///< Container for reduced trigger patches
  TObjArray                               *fReducedClusterInfo;               ///< Container for reduced EMCAL clusters
  TObjArray                               *fReducedParticleInfo;              ///< Container for reduced true particles
  TObjArray                               *fReducedTrackInfo;                 ///< Container for reduced reconstructed tracks

  /// \cond CLASSIMP
  ClassDef(AliReducedHighPtEvent, 2);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDHIGHPTEVENT_H */
