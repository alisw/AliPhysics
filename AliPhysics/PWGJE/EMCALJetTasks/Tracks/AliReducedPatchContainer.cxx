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

#include "AliReducedPatchInfo.h"
#include "AliReducedPatchContainer.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedPatchContainer)
/// \endcond

namespace HighPtTracks {

/**
 * Constructor. In case doAlloc is true, the containers are allocated. Otherwise the
 * constructor becomes the dummy constructor, needed for ROOT I/O
 * \param doAlloc If true containers are allocated
 */
AliReducedPatchContainer::AliReducedPatchContainer(Bool_t doAlloc) :
  TObject()
{
  if(doAlloc){
    for(int icont = 0; icont < 4; icont++){
      fOnlinePatches[icont] = new TObjArray;
      fOnlinePatches[icont]->SetOwner(kTRUE);
      fOfflinePatches[icont] = new TObjArray;
      fOfflinePatches[icont]->SetOwner(kTRUE);
    }
  } else {
    memset(fOnlinePatches, 0, sizeof(TObjArray *) *4);
    memset(fOfflinePatches, 0, sizeof(TObjArray *) *4);
  }
}

/**
 * Copy constructor, taking ownership over the objects
 * \param ref Reference for the copy
 */
AliReducedPatchContainer::AliReducedPatchContainer(const AliReducedPatchContainer &ref):
  TObject(ref)
{
  ref.Copy(*this);
}

/**
 * Assignment operator. First cleans memory assigns to this object, then performs a deep copy.
 * \param ref Reference for the assignment
 * \return this object
 */
AliReducedPatchContainer &AliReducedPatchContainer::operator=(const AliReducedPatchContainer &ref){
  if(&ref != this){
    this->~AliReducedPatchContainer();
    TObject::operator=(ref);
    ref.Copy(*this);
  }
  return *this;
}

/**
 * Implementation of the copy functionality for the patch container. Called in assignment operator and in the
 * copy constructor. Copies entries from this object into the target object.
 * \param target Target where to copy to
 */
void AliReducedPatchContainer::Copy(TObject &target) const{
  AliReducedPatchContainer *targetcont = dynamic_cast<AliReducedPatchContainer *>(&target);
  if(!targetcont) return;
  for(int icont = 0; icont < 4; icont++){
    targetcont->fOnlinePatches[icont] = new TObjArray;
    targetcont->fOnlinePatches[icont]->SetOwner(kTRUE);
    for(TIter refIter = TIter(fOnlinePatches[icont]).Begin(); refIter != TIter::End(); ++ refIter){
      targetcont->fOnlinePatches[icont]->Add(new AliReducedPatchInfo(*(static_cast<AliReducedPatchInfo *>(*refIter))));
    }
    targetcont->fOfflinePatches[icont] = new TObjArray;
    targetcont->fOfflinePatches[icont]->SetOwner(kTRUE);
    for(TIter refIter = TIter(fOfflinePatches[icont]).Begin(); refIter != TIter::End(); ++ refIter){
      targetcont->fOfflinePatches[icont]->Add(new AliReducedPatchInfo(*(static_cast<AliReducedPatchInfo *>(*refIter))));
    }
  }
}

/**
 * Detructor, cleaning up containers
 */
AliReducedPatchContainer::~AliReducedPatchContainer() {
  for(int icont = 0; icont < 4; icont++){
    if(fOnlinePatches[icont]) delete fOnlinePatches[icont];
    if(fOfflinePatches[icont]) delete fOfflinePatches[icont];
  }
}

/**
 * Add new patch to the trigger patch container
 * \param isOffline Online or offline patch
 * \param patchtype Type of the patch
 * \param energy Patch energy
 * \param amplitude Patch amplitude
 * \param eta Patch \f$ \eta \f$
 * \param phi Patch \f$ \phi \f$
 */
void AliReducedPatchContainer::AddTriggerPatch(Bool_t isOffline, PatchType_t patchtype, Float_t energy, Int_t amplitude, Float_t eta, Float_t phi){
  TObjArray *contpointer = isOffline ? fOfflinePatches[patchtype] : fOnlinePatches[patchtype];
  contpointer->Add(new AliReducedPatchInfo(energy, amplitude, eta, phi));
}

/**
 * Get list of reconstructed patches of a given patch type
 * @param isOffline Online or offline patches
 * @param patchtype Type of the trigger patch
 * @return list of reconstructed patches
 */
TObjArray *AliReducedPatchContainer::GetTriggerPatches(Bool_t isOffline, PatchType_t patchtype){
  return isOffline ? fOfflinePatches[patchtype] : fOnlinePatches[patchtype];
}


} /* namespace HighPtTracks */
