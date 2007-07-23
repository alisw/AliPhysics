// $Header$

#include "KineTools.h"

#include <TObject.h>
#include <TTree.h>
#include <TBranchElement.h>
#include <TClonesArray.h>

#include <AliStack.h>
#include <AliTrackReference.h>

#include "Reve/Track.h"
#include "Reve/RenderElement.h"

#include <algorithm>
#include <map>

//______________________________________________________________________
// KineTools
//

using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(KineTools)

KineTools::KineTools()
{}

/**************************************************************************/

void KineTools::SetDaughterPathMarks(RenderElement* cont, AliStack* stack, Bool_t recurse)
{
  // Import daughters birth points.

  RenderElement::List_i  iter = cont->BeginChildren();

  while(iter != cont->EndChildren())
  {
    Track* track = dynamic_cast<Track*>(*iter); 
    TParticle* p = stack->Particle(track->GetLabel());
    if(p->GetNDaughters()) {
      Int_t d0 = p->GetDaughter(0), d1 = p->GetDaughter(1);
      for(int d=d0; d>0 && d<=d1; ++d) 
      {	
	TParticle* dp = stack->Particle(d);
	Reve::PathMark* pm = new PathMark(PathMark::Daughter);
        pm->V.Set(dp->Vx(), dp->Vy(), dp->Vz());
	pm->P.Set(dp->Px(), dp->Py(), dp->Pz()); 
        pm->time = dp->T();
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
  bool operator()(PathMark* const & a, PathMark* const & b)
  { return a->time < b->time; }
};

void slurp_tracks(map<Int_t, Track*>& tracks, RenderElement* cont, Bool_t recurse)
{
  RenderElement::List_i citer = cont->BeginChildren();
  while(citer != cont->EndChildren())
  { 
    Track* track = dynamic_cast<Track*>(*citer); 
    tracks[track->GetLabel()] = track;
    if (recurse)
      slurp_tracks(tracks, track, recurse);
    ++citer;
  }
}

}

void KineTools::SetTrackReferences(RenderElement* cont, TTree* treeTR, Bool_t recurse)
{
  // set decay and reference points

  static const Exc_t eH("KineTools::ImportPathMarks");

  // Fill map
  map<Int_t, Track*> tracks;
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
      map<Int_t, Track*>::iterator iter = tracks.end(); 
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
	  PathMark* pm = new PathMark(isRef ? PathMark::Reference : PathMark::Decay);
	  pm->V.Set(atr->X(),atr->Y(), atr->Z());
	  pm->P.Set(atr->Px(),atr->Py(), atr->Pz());  
	  pm->time = atr->GetTime();
          Track* track  = iter->second;
          track->AddPathMark(pm);
	}
      } // loop track refs 
    } // loop primaries, clones arrays
    delete arr;
  } // end loop through top branches
}

void KineTools::SortPathMarks(RenderElement* cont, Bool_t recurse)
{
  // Sort path-marks for all tracks by time.

  // Fill map
  map<Int_t, Track*> tracks;
  slurp_tracks(tracks, cont, recurse);

  // sort 
  for(map<Int_t, Track*>::iterator j=tracks.begin(); j!=tracks.end(); ++j)
  {
    j->second->SortPathMarksByTime();
  }
}
