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

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData : public TObject {
public:
  AliEMCalTriggerEventData();
  AliEMCalTriggerEventData(const AliEMCalTriggerEventData &ref);
  AliEMCalTriggerEventData &operator=(const AliEMCalTriggerEventData &ref);
  virtual ~AliEMCalTriggerEventData() {}

  AliVEvent    * GetRecEvent() const { return fRecEvent; }
  const AliMCEvent   * GetMCEvent() const { return fMCEvent; }
  const TClonesArray * GetClusterContainer() const { return fClusterContainer; }
  const TClonesArray * GetMatchedTrackContainer() const { return fTrackContainer; }
  const TClonesArray * GetParticleContainer() const { return fParticleContainer; }
  const TClonesArray * GetTriggerPatchContainer() const { return fTriggerPatchContainer; }
  const AliJetContainer * GetJetContainerData() const { return fJetContainerData; }
  const AliJetContainer * GetJetContainerMC() const { return fJetContainerMC; }

  void SetRecEvent(AliVEvent * const ev) { fRecEvent = ev; }
  void SetMCEvent(const AliMCEvent * const ev) { fMCEvent = ev; }
  void SetClusterContainer(const TClonesArray *const cont) { fClusterContainer = cont; }
  void SetTrackContainer(const TClonesArray * const cont) { fTrackContainer = cont; }
  void SetParticleContainer(const TClonesArray * const cont) { fParticleContainer = cont ;}
  void SetTriggerPatchContainer(const TClonesArray *const cont) { fTriggerPatchContainer = cont; }
  void SetMCJetContainer(const AliJetContainer * const cont) { fJetContainerMC = cont; }
  void SetDataJetContainer(const AliJetContainer * const cont) { fJetContainerData = cont; }

protected:
  AliVEvent             *fRecEvent;                     // Reconstructed event
  const AliMCEvent      *fMCEvent;                      // Monte-Carlo event
  const TClonesArray    *fClusterContainer;             // Container with calibrated clusters
  const TClonesArray    *fTrackContainer;               // Container with matched tracks
  const TClonesArray    *fParticleContainer;            // Container with MC-true filtered particles
  const TClonesArray    *fTriggerPatchContainer;        // Container with trigger patches
  const AliJetContainer *fJetContainerMC;               // Container with reconstructed jets
  const AliJetContainer *fJetContainerData;             // Container with reconstructed jets

  ClassDef(AliEMCalTriggerEventData, 1);          // Data structure exchanged to trigger analysis components
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGEREVENTDATA_H */
