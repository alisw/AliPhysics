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

#include <algorithm>
#include <map>

//______________________________________________________________________
// AliEveKineTools
//

using namespace std;

ClassImp(AliEveKineTools)

AliEveKineTools::AliEveKineTools()
{}

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

namespace {
struct cmp_pathmark
{
  bool operator()(TEvePathMark* const & a, TEvePathMark* const & b)
  { return a->fTime < b->fTime; }
};

void slurp_tracks(map<Int_t, TEveTrack*>& tracks, TEveElement* cont, Bool_t recurse)
{
  TEveElement::List_i citer = cont->BeginChildren();
  while(citer != cont->EndChildren())
  { 
    TEveTrack* track = dynamic_cast<TEveTrack*>(*citer); 
    tracks[track->GetLabel()] = track;
    if (recurse)
      slurp_tracks(tracks, track, recurse);
    ++citer;
  }
}

}

void AliEveKineTools::SetTrackReferences(TEveElement* cont, TTree* treeTR, Bool_t recurse)
{
  // set decay and reference points

  static const TEveException eH("AliEveKineTools::ImportPathMarks");

  // Fill map
  map<Int_t, TEveTrack*> tracks;
  slurp_tracks(tracks, cont, recurse);
 
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
      map<Int_t, TEveTrack*>::iterator iter = tracks.end(); 
      Int_t Nent =  arr->GetEntriesFast();
      for (Int_t iTrackRef = 0; iTrackRef < Nent; iTrackRef++) 
      {
	AliTrackReference* atr = (AliTrackReference*)arr->UncheckedAt(iTrackRef);

	Int_t label = atr->GetTrack();
	if (label < 0)
	  throw(eH + Form("negative label for entry %d in branch %s.",
			  iTrackRef, el->GetName()));
	
        if(label != last_label) {
	  iter = tracks.find(label);
	  last_label = label;
	}

	if (iter != tracks.end()) {
	  TEvePathMark* pm = new TEvePathMark(isRef ? TEvePathMark::kReference : TEvePathMark::kDecay);
	  pm->fV.Set(atr->X(),atr->Y(), atr->Z());
	  pm->fP.Set(atr->Px(),atr->Py(), atr->Pz());  
	  pm->fTime = atr->GetTime();
          TEveTrack* track  = iter->second;
          track->AddPathMark(pm);
	}
      } // loop track refs 
    } // loop primaries, clones arrays
    delete arr;
  } // end loop through top branches
}

void AliEveKineTools::SortPathMarks(TEveElement* cont, Bool_t recurse)
{
  // Sort path-marks for all tracks by time.

  // Fill map
  map<Int_t, TEveTrack*> tracks;
  slurp_tracks(tracks, cont, recurse);

  // sort 
  for(map<Int_t, TEveTrack*>::iterator j=tracks.begin(); j!=tracks.end(); ++j)
  {
    j->second->SortPathMarksByTime();
  }
}
