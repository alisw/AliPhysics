/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <TObjArray.h>

#include "AliLog.h"

#include "AliReducedEmcalCluster.h"
#include "AliReducedGeneratedParticle.h"
#include "AliReducedMCHeader.h"
#include "AliReducedPatchContainer.h"
#include "AliReducedReconstructedTrack.h"
#include "AliReducedHighPtEvent.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedHighPtEvent)
/// \endcond

namespace HighPtTracks {

/**
 * Constructor. If doAlloc is true, then the containers are allocated. Otherwise the constructor becomes
 * the ROOT dummy constructor needed for I/O
 * \param doAlloc if true then the containers are allocated
 */
AliReducedHighPtEvent::AliReducedHighPtEvent(Bool_t doAlloc):
  TObject(),
  fRunNumber(0),
  fCentralityPercentile(100.),
  fVertexZ(-999.9),
  fMCHeader(NULL),
  fIsMinBias(kFALSE),
  fReducedPatchInfo(NULL),
  fReducedClusterInfo(NULL),
  fReducedParticleInfo(NULL),
  fReducedTrackInfo(NULL)
{
  for(int i = 0; i < 2; i++){
    fJetTriggerString[i] = kFALSE;
    fGammaTriggerString[i] = kFALSE;
  }
  if(doAlloc){
    fReducedPatchInfo = new AliReducedPatchContainer(doAlloc);
    fReducedClusterInfo = new TObjArray;
    fReducedClusterInfo->SetOwner(kTRUE);
    fReducedTrackInfo = new TObjArray;
    fReducedTrackInfo->SetOwner(kTRUE);
  }
}

/**
 * Copy constructor, taking ownership over pointer members. For this a deep copy is performed. Copy
 * functionality is implemented in the function Copy
 * \param ref The reference for the copy
 */
AliReducedHighPtEvent::AliReducedHighPtEvent(const AliReducedHighPtEvent& ref):
  TObject(ref),
  fRunNumber(0),
  fCentralityPercentile(ref.fCentralityPercentile),
  fVertexZ(ref.fVertexZ),
  fMCHeader(NULL),
  fIsMinBias(ref.fIsMinBias),
  fReducedPatchInfo(NULL),
  fReducedClusterInfo(NULL),
  fReducedParticleInfo(NULL),
  fReducedTrackInfo(NULL)
{
  ref.Copy(*this);
}

/**
 * Assignment operator, taking ownership over pointer members. For this a deep copy is performed. Copy
 * functionality is implemented in the function Copy
 * \param ref The reference for the copy
 */
AliReducedHighPtEvent& AliReducedHighPtEvent::operator=(const AliReducedHighPtEvent& ref) {
  if(this != &ref){
    this->~AliReducedHighPtEvent();
    TObject::operator=(ref);
    ref.Copy(*this);

  }
  return *this;
}

/**
 * Destructor, cleaning up allocated memory
 */
AliReducedHighPtEvent::~AliReducedHighPtEvent() {
  if(fReducedPatchInfo) delete fReducedPatchInfo;
  if(fMCHeader) delete fMCHeader;
  if(fReducedClusterInfo) delete fReducedClusterInfo;
  if(fReducedParticleInfo) delete fReducedParticleInfo;
  if(fReducedTrackInfo) delete fReducedTrackInfo;
}

/**
 * Implementation of the copy functionality. Copies information from this object
 * into the target object. Performs a deep copy.
 * @param target
 */
void AliReducedHighPtEvent::Copy(TObject& target) const {
  AliReducedHighPtEvent *targetevent = dynamic_cast<AliReducedHighPtEvent *>(&target);
  if(!targetevent) return;
  targetevent->fCentralityPercentile = fCentralityPercentile;
  targetevent->fVertexZ = fVertexZ;
  targetevent->fIsMinBias = fIsMinBias;
  memcpy(targetevent->fGammaTriggerString, fGammaTriggerString, sizeof(Bool_t) *2);
  memcpy(targetevent->fJetTriggerString, fJetTriggerString, sizeof(Bool_t) *2);
  if(fMCHeader) targetevent->fMCHeader = new AliReducedMCHeader(*fMCHeader);
  if(fReducedPatchInfo) targetevent->fReducedPatchInfo = new AliReducedPatchContainer(*fReducedPatchInfo);
  if(fReducedClusterInfo){
    targetevent->fReducedClusterInfo = new TObjArray;
    targetevent->fReducedClusterInfo->SetOwner(kTRUE);
    for(TIter clustiter = TIter(fReducedClusterInfo).Begin(); clustiter != TIter::End(); ++clustiter)
      targetevent->fReducedClusterInfo->Add(new AliReducedEmcalCluster(*(static_cast<AliReducedEmcalCluster *>(*clustiter))));
  }
  if(fReducedParticleInfo){
    targetevent->fReducedParticleInfo = new TObjArray;
    targetevent->fReducedParticleInfo->SetOwner(kTRUE);
    for(TIter partiter = TIter(fReducedParticleInfo).Begin(); partiter != TIter::End(); ++partiter)
      targetevent->fReducedParticleInfo->Add(new AliReducedGeneratedParticle(*(static_cast<const AliReducedGeneratedParticle *>(*partiter))));
  }
  if(fReducedTrackInfo){
    targetevent->fReducedTrackInfo = new TObjArray;
    targetevent->fReducedTrackInfo->SetOwner(kTRUE);
    for(TIter partiter = TIter(fReducedTrackInfo).Begin(); partiter != TIter::End(); ++partiter)
      targetevent->fReducedTrackInfo->Add(new AliReducedReconstructedTrack(*(static_cast<const AliReducedReconstructedTrack *>(*partiter))));
  }
}

/**
 * Find EMCAL cluster by internal cluster ID
 * \param index ID of the EMCAL cluster
 * \return The found EMCAL cluster
 */
AliReducedEmcalCluster* AliReducedHighPtEvent::GetClusterForIndex(Int_t index) {
  if(!fReducedClusterInfo) return NULL;
  AliReducedEmcalCluster *foundCluster = NULL, *tmpcluster = NULL;
  for(TIter clustIter = TIter(fReducedClusterInfo).Begin(); clustIter != TIter::End(); ++clustIter){
    tmpcluster = static_cast<AliReducedEmcalCluster *>(*clustIter);
    if(tmpcluster->GetClusterID() == index){
      foundCluster = tmpcluster;
      break;
    }
  }
  return foundCluster;
}

/**
 * Find generated particle in the event with a given unique ID
 * \param index unique ID of the particle
 * \return the matching particle at generator level (NULL if no particle is available or particle is not found)
 */
AliReducedGeneratedParticle* AliReducedHighPtEvent::GetParticleForIndex(Int_t index) {
  if(!fReducedParticleInfo) return NULL;
  AliReducedGeneratedParticle *foundparticle(NULL), *tmpparticle(NULL);
  for(TIter partiter = TIter(fReducedParticleInfo).Begin(); partiter != TIter::End(); ++partiter){
    tmpparticle = static_cast<AliReducedGeneratedParticle *>(*partiter);
    if(tmpparticle->GetID() == index){
      foundparticle =  tmpparticle;
      break;
    }
  }
  return foundparticle;
}

/**
 * Add cluster information to the reduced event
 * \param cluster The reduced cluster to be added
 */
void AliReducedHighPtEvent::AddReducedCluster(AliReducedEmcalCluster* cluster) {
  if(!fReducedClusterInfo){
    AliError("Cluster container not allocated");
    return;
  }
  fReducedClusterInfo->Add(cluster);
}

/**
 * Add reduced generated particle to the list of particles
 * \param part The reduced particle to be added
 */
void AliReducedHighPtEvent::AddReducedGeneratedParticle(AliReducedGeneratedParticle* part) {
  if(!fReducedParticleInfo){
    fReducedParticleInfo = new TObjArray;
    fReducedParticleInfo->SetOwner(kTRUE);
  }
  fReducedParticleInfo->Add(part);
}

/**
 * Add reduced reconstructed particle to the list of tracks
 * \param trk The track to be added
 */
void AliReducedHighPtEvent::AddReducedReconstructedParticle(AliReducedReconstructedTrack* trk) {
  if(!fReducedTrackInfo){
    AliError("Track container not allocated");
    return;
  }
  fReducedTrackInfo->Add(trk);
}

/**
 * Create a stl vector with reduced Emcal cluster information
 * \return vector with reduced Emcal clusters
 */
std::vector<HighPtTracks::AliReducedEmcalCluster*> AliReducedHighPtEvent::GetClusterVector() const {
  std::vector<HighPtTracks::AliReducedEmcalCluster*> resultvec;
  for(TIter clustIter = TIter(fReducedClusterInfo).Begin(); clustIter != TIter::End(); ++clustIter){
    resultvec.push_back(static_cast<AliReducedEmcalCluster *>(*clustIter));
  }
  return resultvec;
}

/**
 * Create a stl vector with reduced generated particle information
 * \return vector with reduced generated particles
 */
std::vector<HighPtTracks::AliReducedGeneratedParticle*> AliReducedHighPtEvent::GetParticleVector() const {
  std::vector<HighPtTracks::AliReducedGeneratedParticle*> resultvec;
  for(TIter partIter = TIter(fReducedParticleInfo).Begin(); partIter != TIter::End(); ++partIter){
    resultvec.push_back(static_cast<AliReducedGeneratedParticle *>(*partIter));
  }
  return resultvec;
}

/**
 * Create a stl vector with reduced reconstructed
 * \return vector with reduced reconstructed tracks
 */
std::vector<HighPtTracks::AliReducedReconstructedTrack*> HighPtTracks::AliReducedHighPtEvent::GetTrackVector() const {
  std::vector<HighPtTracks::AliReducedReconstructedTrack*> resultvec;
  for(TIter trackIter = TIter(fReducedTrackInfo).Begin(); trackIter != TIter::End(); ++trackIter){
    resultvec.push_back(static_cast<AliReducedReconstructedTrack *>(*trackIter));
  }
  return resultvec;
}


} /* namespace HighPtTracks */

