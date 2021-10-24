/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliEMCalTriggerEventData.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerEventData)

using namespace PWGJE::EMCALJetTasks;

/**
 * Default constructor
 */
AliEMCalTriggerEventData::AliEMCalTriggerEventData() :
    TObject(),
    fRecEvent(NULL),
    fMCEvent(NULL),
    fTriggerBitSelection(0),
    fClusterContainer(NULL),
    fTrackContainer(NULL),
    fParticleContainer(NULL),
    fTriggerPatchContainer(NULL),
    fJetContainerMC(NULL),
    fJetContainerData(NULL)
{
}

/**
 * Copy constructor
 */
AliEMCalTriggerEventData::AliEMCalTriggerEventData(const AliEMCalTriggerEventData &ref) :
    TObject(ref),
    fRecEvent(ref.fRecEvent),
    fMCEvent(ref.fMCEvent),
    fTriggerBitSelection(ref.fTriggerBitSelection),
    fClusterContainer(ref.fClusterContainer),
    fTrackContainer(ref.fTrackContainer),
    fParticleContainer(ref.fParticleContainer),
    fTriggerPatchContainer(ref.fTriggerPatchContainer),
    fJetContainerMC(ref.fJetContainerMC),
    fJetContainerData(ref.fJetContainerData)
{
}

/**
 * assignment operator
 */
AliEMCalTriggerEventData &AliEMCalTriggerEventData::operator=(const AliEMCalTriggerEventData &ref) {
  TObject::operator=(ref);
  if(this != &ref){
    fRecEvent = ref.fRecEvent;
    fTriggerBitSelection = ref.fTriggerBitSelection;
    fMCEvent = ref.fMCEvent;
    fClusterContainer = ref.fClusterContainer;
    fTrackContainer = ref.fTrackContainer;
    fParticleContainer = ref.fParticleContainer;
    fTriggerPatchContainer = ref.fTriggerPatchContainer;
    fJetContainerMC = ref.fJetContainerMC;
    fJetContainerData = ref.fJetContainerData;
  }
  return *this;
}
