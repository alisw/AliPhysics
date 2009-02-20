// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveKineTools.h"
#include "AliEveTrack.h"

#include <AliStack.h>
#include <AliTrackReference.h>

#include <TTree.h>
#include <TBranchElement.h>
#include <TClonesArray.h>

#include <map>

//______________________________________________________________________________
// AliEveKineTools
//
// Tools for import of kinematics.
//

ClassImp(AliEveKineTools)

namespace {

  // Map to store label-to-track association.
  //
  // multimap is used as there are cases when initial particles (in
  // particular resonances) are not assigned proper status-codes
  // and can thus be found several times in the eve-track-list.

  typedef std::multimap<Int_t, AliEveTrack*>                 TrackMap_t;
  typedef std::multimap<Int_t, AliEveTrack*>::const_iterator TrackMap_ci;

  void MapTracks(TrackMap_t& map, TEveElement* cont, Bool_t recurse)
  {
    TEveElement::List_i i = cont->BeginChildren();
    while (i != cont->EndChildren()) {
      AliEveTrack* track = dynamic_cast<AliEveTrack*>(*i);
      map.insert(std::make_pair(track->GetLabel(), track));
      if (recurse)
        MapTracks(map, track, recurse);
      ++i;
    }
  }
}

/**************************************************************************/

void AliEveKineTools::SetDaughterPathMarks(TEveElement* cont, AliStack* stack, Bool_t recurse)
{
  // Import daughters birth points.

  TEveElement::List_i  iter = cont->BeginChildren();
  while(iter != cont->EndChildren())
  {
    AliEveTrack* track = dynamic_cast<AliEveTrack*>(*iter);
    TParticle* p = stack->Particle(track->GetLabel());
    if (p->GetNDaughters())
    {
      Int_t d0 = p->GetDaughter(0), d1 = p->GetDaughter(1);
      for(int d = d0; d > 0 && d <= d1; ++d)
      {
	TParticle* dp = stack->Particle(d);
        track->AddPathMark(TEvePathMark(TEvePathMark::kDaughter,
                                        TEveVector(dp->Vx(), dp->Vy(), dp->Vz()),
                                        TEveVector(dp->Px(), dp->Py(), dp->Pz()),
                                        dp->T()));
        // printf("Daughter path-mark for %d, %d, t=%e, r=%f,%f,%f\n",
        //        track->GetLabel(), d, dp->T(), dp->Vx(), dp->Vy(), dp->Vz());
      }

      // Check last process, set decay if needed.
      Int_t lp = stack->Particle(d1)->GetUniqueID();
      if (lp != kPBrem && lp != kPDeltaRay && lp < kPCerenkov)
      {
        TParticle* dp = stack->Particle(d1);
        track->AddPathMark(TEvePathMark(TEvePathMark::kDecay,
                                        TEveVector(dp->Vx(), dp->Vy(), dp->Vz()),
                                        TEveVector(0, 0,0),  dp->T()));
      }

      if (recurse)
	SetDaughterPathMarks(track, stack, recurse);
    }
    ++iter;
  }
}

/**************************************************************************/

void AliEveKineTools::SetTrackReferences(TEveElement* cont, TTree* treeTR, Bool_t recurse)
{
  // Set decay and track reference path-marks.

  static const TEveException kEH("AliEveKineTools::ImportPathMarks");

  TrackMap_t map;
  MapTracks(map, cont, recurse);

  TClonesArray* arr = 0;
  treeTR->SetBranchAddress("TrackReferences", &arr);

  Int_t nTreeEntries = (Int_t) treeTR->GetEntries();

  for (Int_t te = 0; te < nTreeEntries; ++te)
  {
    treeTR->GetEntry(te);

    Int_t last_label = -1;
    std::pair<TrackMap_ci, TrackMap_ci> range;
    Int_t nArrEntries = arr->GetEntriesFast();

    // printf("tree-entry %d, n-arr-entries %d\n", te, nArrEntries);

    for (Int_t ae = 0; ae < nArrEntries; ++ae)
    {
      AliTrackReference* atr = (AliTrackReference*)arr->UncheckedAt(ae);
      Bool_t isRef = (atr->DetectorId() != -1);
      Int_t  label = atr->GetTrack();

      // printf("    arr-entry %d, label %d, detid %d, len=%f, t=%e r=%f,%f,%f\n",
      //        ae, label, atr->DetectorId(), 
      //        atr->GetLength(), atr->GetTime(), atr->X(), atr->Y(), atr->Z());

      if (label < 0)
        throw(kEH + Form("negative label for array-entry %d in tree-entry %d.",
                         ae, te));

      if (label != last_label)
      {
        range      = map.equal_range(label);
        last_label = label;
      }

      for (TrackMap_ci i = range.first; i != range.second; ++i)
      {
	i->second->AddPathMark
          (TEvePathMark(isRef ? TEvePathMark::kReference : TEvePathMark::kDecay,
                        TEveVector(atr->X(),  atr->Y(),  atr->Z()),
                        TEveVector(atr->Px(), atr->Py(), atr->Pz()),
                        atr->GetTime()));
      }
    }
  }
  delete arr;
}

/**************************************************************************/

void AliEveKineTools::SortPathMarks(TEveElement* el, Bool_t recurse)
{
  // Sort path-marks for track by time.
  // If recurse is true, descends down all track children.

  AliEveTrack* track = dynamic_cast<AliEveTrack*>(el);
  if (track)
    track->SortPathMarksByTime();

  if (recurse)
  {
    for (TEveElement::List_i i = el->BeginChildren(); i != el->EndChildren(); ++i)
      SortPathMarks(*i, kTRUE);
  }
}
