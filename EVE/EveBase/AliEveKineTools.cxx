// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveKineTools.h"

#include <TObject.h>
#include <TTree.h>
#include <TBranchElement.h>
#include <TClonesArray.h>

#include <AliStack.h>
#include <AliTrackReference.h>

#include <TEveTrack.h>
#include <TEveElement.h>

#include <map>

//______________________________________________________________________________
// AliEveKineTools
//
// Tools for import of kinematics. Preliminary version.
//

using namespace std;

ClassImp(AliEveKineTools)


namespace {

  typedef std::map<Int_t, TEveTrack*> TrackMap_t;

  void MapTracks(TrackMap_t& map, TEveElement* cont, Bool_t recurse)
  {
    TEveElement::List_i i = cont->BeginChildren();
    while (i != cont->EndChildren()) {
      TEveTrack* track = dynamic_cast<TEveTrack*>(*i);
      map[track->GetLabel()] = track;
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
    TEveTrack* track = dynamic_cast<TEveTrack*>(*iter);
    TParticle* p = stack->Particle(track->GetLabel());
    if(p->GetNDaughters()) {
      Int_t d0 = p->GetDaughter(0), d1 = p->GetDaughter(1);
      for(int d=d0; d>0 && d<=d1; ++d)
      {
	TParticle* dp = stack->Particle(d);
	TEvePathMark* pm = new TEvePathMark(TEvePathMark::kDaughter);
        pm->fV.Set(dp->Vx(), dp->Vy(), dp->Vz());
	pm->fP.Set(dp->Px(), dp->Py(), dp->Pz());
        pm->fTime = dp->T();
        track->AddPathMark(pm);
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

  static const TEveException eH("AliEveKineTools::ImportPathMarks");

  TrackMap_t map;
  MapTracks(map, cont, recurse);

  Int_t nPrimaries = (Int_t) treeTR->GetEntries();
  TIter next(treeTR->GetListOfBranches());
  TBranchElement* el;
  Bool_t isRef = kTRUE;

  while ((el = (TBranchElement*) next()))
  {
    if (strcmp("AliRun",el->GetName()) == 0)
      isRef = kFALSE;

    TClonesArray* arr = 0;
    el->SetAddress(&arr);
    for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++)
    {
      el->GetEntry(iPrimPart);

      Int_t last_label = -1;
      TrackMap_t::iterator iter = map.end();
      Int_t Nent =  arr->GetEntriesFast();
      for (Int_t iTrackRef = 0; iTrackRef < Nent; iTrackRef++)
      {
	AliTrackReference* atr = (AliTrackReference*)arr->UncheckedAt(iTrackRef);

	Int_t label = atr->GetTrack();
	if (label < 0)
	  throw(eH + Form("negative label for entry %d in branch %s.",
			  iTrackRef, el->GetName()));

        if(label != last_label) {
	  iter = map.find(label);
	  last_label = label;
	}

	if (iter != map.end()) {
	  TEvePathMark* pm = new TEvePathMark(isRef ? TEvePathMark::kReference : TEvePathMark::kDecay);
	  pm->fV.Set(atr->X(),atr->Y(), atr->Z());
	  pm->fP.Set(atr->Px(),atr->Py(), atr->Pz());
	  pm->fTime = atr->GetTime();
          TEveTrack* track  = iter->second;
          track->AddPathMark(pm);
	}
      } // loop track refs
    } // loop primaries in clones arrays
    delete arr;
  } // end loop through top branches
}

/**************************************************************************/

void AliEveKineTools::SortPathMarks(TEveElement* el, Bool_t recurse)
{
  // Sort path-marks for all tracks by time.

  TEveTrack* track = dynamic_cast<TEveTrack*>(el);
  if(track) track->SortPathMarksByTime();

  TEveElement::List_i i = el->BeginChildren();
  while (i != el->EndChildren() && recurse) {
    track = dynamic_cast<TEveTrack*>(el);
    if (track) track->SortPathMarksByTime();
    i++;
  }
}
