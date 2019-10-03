// $Id$
//
// Jet model task to merge to existing branches
// only implemented for track branches
//
// Author: M. Verweij

#include "AliJetModelMergeBranches.h"

#include <TClonesArray.h>
#include <TFolder.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include <TProfile.h>
#include "AliAnalysisManager.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

ClassImp(AliJetModelMergeBranches)

//________________________________________________________________________
AliJetModelMergeBranches::AliJetModelMergeBranches() : 
  AliJetModelBaseTask("AliJetModelMergeBranches"),
  fTracksMergeName(),
  fTracksMerge(0)
{
  // Default constructor.
  SetSuffix("Merged");
  SetCopyArray(kTRUE);
}

//________________________________________________________________________
AliJetModelMergeBranches::AliJetModelMergeBranches(const char *name) : 
  AliJetModelBaseTask(name),
  fTracksMergeName(),
  fTracksMerge(0)
{
  // Standard constructor.
  SetSuffix("Merged");
  SetCopyArray(kTRUE);
}

//________________________________________________________________________
AliJetModelMergeBranches::~AliJetModelMergeBranches()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliJetModelMergeBranches::ExecOnce() 
{
  // Exec only once.

  AliJetModelBaseTask::ExecOnce();
  
  if (!fTracksMergeName.IsNull()) {
    fTracksMerge = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksMergeName));
    if (!fTracksMerge) {
      AliWarning(Form("%s: Couldn't retrieve tracks with name %s!", GetName(), fTracksMergeName.Data()));
    }
    else if (!fTracksMerge->GetClass()->GetBaseClass("AliPicoTrack")) {
      AliError(Form("%s: Collection %s does not contain AliPicoTrack objects!", GetName(), fTracksMergeName.Data())); 
      fTracksMerge = 0;
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
void AliJetModelMergeBranches::Run() 
{
  // Merge two branches into a new one

  CopyTracks();
  MergeTracks();
}

//________________________________________________________________________
void AliJetModelMergeBranches::MergeTracks()
{
  // Copy all the tracks in the new collection

  if (!fTracksMerge)
    return;

  const Int_t nTracks = fTracksMerge->GetEntriesFast();
  Int_t nCopiedTracks = fOutTracks->GetEntriesFast();
  for (Int_t i = 0; i < nTracks; ++i) {
    AliPicoTrack *picotrack = static_cast<AliPicoTrack*>(fTracksMerge->At(i));
    if (!picotrack)
      continue;
    AliPicoTrack *track = new ((*fOutTracks)[nCopiedTracks]) AliPicoTrack(*picotrack);
    track->SetBit(TObject::kBitMask,1);
    nCopiedTracks++;
  }
}
