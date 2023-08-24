/**
 * \file AliEMCalTriggerEventData.h
 * \brief Event Data used in exchange to the different analysis components
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGEREVENTDATA_H
#define ALIEMCALTRIGGEREVENTDATA_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TObject.h>

class TClonesArray;
class AliJetContainer;
class AliMCEvent;
class AliVEvent;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * \class AliEMCalTriggerEventData
 * \brief Simple event container within the high-\f$ p_{t} \f$ track analysis
 *
 * This class combines the reconstructed event (and if available the Monte-Carlo evnet) with the
 * containers from the EMCAL jet framework relevant for this taks and the
 */
class AliEMCalTriggerEventData : public TObject {
public:
  AliEMCalTriggerEventData();
  AliEMCalTriggerEventData(const AliEMCalTriggerEventData &ref);
  AliEMCalTriggerEventData &operator=(const AliEMCalTriggerEventData &ref);
  /**
   * Destructor
   */
  virtual ~AliEMCalTriggerEventData() {}

  /**
   * Get the reconstructed event
   * \return The reconstructed event
   */
  AliVEvent    * GetRecEvent() const { return fRecEvent; }
  /**
   * Get the trigger bit selection
   * \return The trigger bit selection of the event
   */
  UInt_t        GetTriggerBitSelection() const { return fTriggerBitSelection; }
  /**
   * Get the MC event (if available)
   * \return The corresponding Monte-Carlo event
   */
  AliMCEvent   * GetMCEvent() const { return fMCEvent; }
  /**
   * Get the cluster container
   * \return the container with corrected clusters
   */
  const TClonesArray * GetClusterContainer() const { return fClusterContainer; }
  /**
   * Get the container with tracks used for cluster-track matching and jet finding
   * \return The container with selected tracks
   */
  const TClonesArray * GetMatchedTrackContainer() const { return fTrackContainer; }
  /**
   * Get the container with MC-true particles selected
   * \return The container with true particles
   */
  const TClonesArray * GetParticleContainer() const { return fParticleContainer; }
  /**
   * Get the container with EMCAL trigger patches
   * \return The container with EMCAL trigger patches
   */
  const TClonesArray * GetTriggerPatchContainer() const { return fTriggerPatchContainer; }
  /**
   * Get the container with jets in the reconstructed event
   * \return The container with jets (if available - null otherwise)
   */
  AliJetContainer * GetJetContainerData() const { return fJetContainerData; }
  /**
   * Get the container with jets in the Monte-Carlo event
   * \return The container with jets (if available - null otherwise)
   */
  AliJetContainer * GetJetContainerMC() const { return fJetContainerMC; }

  /**
   * Set the reconstructed event
   * \param ev The reconstructed event (ESD or AOD event)
   */
  void SetRecEvent(AliVEvent * const ev) { fRecEvent = ev; }
  /**
   * Set the trigger bit selection for the given event
   * \param triggerbits
   */
  void SetTriggerBitSelection(Int_t triggerbits)  { fTriggerBitSelection = triggerbits; }
  /**
   * Set the corresponding Monte-Carlo event
   * \param ev The corresponding Monte-Carlo Event
   */
  void SetMCEvent(AliMCEvent * const ev) { fMCEvent = ev; }
  /**
   * Set the container with corrected EMCAL clusters
   * \param cont The container with corrected EMCAL clusters
   */
  void SetClusterContainer(const TClonesArray *const cont) { fClusterContainer = cont; }
  /**
   * Set the container with tracks selected for cluster-track matching and jet finding
   * \param cont The container with selected tracks
   */
  void SetTrackContainer(const TClonesArray * const cont) { fTrackContainer = cont; }
  /**
   * Set the container with selected MC-true particles
   * \param cont The container with selected particles
   */
  void SetParticleContainer(const TClonesArray * const cont) { fParticleContainer = cont ;}
  /**
   * Set the container with reconstructed trigger patches
   * \param cont The container with reconstructed trigger patches
   */
  void SetTriggerPatchContainer(const TClonesArray *const cont) { fTriggerPatchContainer = cont; }
  /**
   * Set the data with jets found in the Monte-Carlo event
   * \param cont The container with found jets
   */
  void SetMCJetContainer(AliJetContainer * const cont) { fJetContainerMC = cont; }
  /**
   * Set the data with jets found in the reconstructed event
   * \param cont The container with found jets
   */
  void SetDataJetContainer(AliJetContainer * const cont) { fJetContainerData = cont; }

protected:
  AliVEvent             *fRecEvent;                     ///< Reconstructed event
  AliMCEvent            *fMCEvent;                      ///< Monte-Carlo event
  Int_t                 fTriggerBitSelection;           ///< Event trigger bit selection
  const TClonesArray    *fClusterContainer;             ///< Container with calibrated clusters
  const TClonesArray    *fTrackContainer;               ///< Container with matched tracks
  const TClonesArray    *fParticleContainer;            ///< Container with MC-true filtered particles
  const TClonesArray    *fTriggerPatchContainer;        ///< Container with trigger patches
  AliJetContainer       *fJetContainerMC;               ///< Container with reconstructed jets
  AliJetContainer       *fJetContainerData;             ///< Container with reconstructed jets

  ClassDef(AliEMCalTriggerEventData, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGEREVENTDATA_H */
