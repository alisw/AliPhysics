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
{

}

/**************************************************************************/
void KineTools::SetDaughterPathMarks(TrackList* cont,  AliStack* stack )
{
  // import daughters birth points 
  RenderElement::List_i  iter = cont->BeginChildren();

  while(iter != cont->EndChildren())
  {
    Track* track = dynamic_cast<Track*>(*iter); 
    TParticle* p = stack->Particle(track->GetLabel());
    if(p->GetNDaughters()) {
      Int_t d0 = p->GetDaughter(0), d1 = p->GetDaughter(1);
      for(int d=d0; d>0 && d<=d1;++d) 
      {	
	TParticle* dp = stack->Particle(d);
	Reve::PathMark* pm = new PathMark( PathMark::Daughter);
        pm->V.Set(dp->Vx(),dp->Vy(), dp->Vz());
	pm->P.Set(dp->Px(),dp->Py(), dp->Pz()); 
        pm->time = dp->T();
        track->AddPathMark(pm);
      }
    }
    ++iter;
  }
}

/**************************************************************************/

namespace {
struct cmp_pathmark {
  bool operator()(PathMark* const & a, PathMark* const & b)
  { return a->time < b->time; }
};
}

void KineTools::SetPathMarks(TrackList* cont, AliStack* stack , TTree* treeTR)
{
  // set decay and reference points

  static const Exc_t eH("KineTools::ImportPathMarks");

  if(treeTR == 0) {
    SetDaughterPathMarks(cont, stack);
    return;
  }

  map<Int_t, Track::vpPathMark_t > refs;

  Int_t nPrimaries = (Int_t) treeTR->GetEntries();
  TIter next(treeTR->GetListOfBranches());
  TBranchElement* el;
  Bool_t isRef = kTRUE;
  treeTR->SetBranchStatus("*",0);

  while ((el = (TBranchElement*) next()))
  {
    if (strcmp("AliRun",el->GetName()) == 0)
      isRef = kFALSE;

    treeTR->SetBranchStatus(Form("%s*", el->GetName()), 1);
    for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) 
    {
     
      TClonesArray* arr = 0;
      treeTR->SetBranchAddress(el->GetName(), &arr);
      treeTR->GetEntry(iPrimPart);
    
      for (Int_t iTrackRef = 0; iTrackRef < arr->GetEntriesFast(); iTrackRef++) 
      {
	AliTrackReference* atr = (AliTrackReference*)arr->At(iTrackRef);
	Int_t track = atr->GetTrack();
        if(atr->TestBit(TObject::kNotDeleted)) {
	  if(track > 0)
	  { 
            PathMark* pm;
	    if(isRef) 
	      pm = new PathMark(PathMark::Reference);
	    else
 	      pm = new PathMark(PathMark::Decay);
	      
	    pm->V.Set(atr->X(),atr->Y(), atr->Z());
	    pm->P.Set(atr->Px(),atr->Py(), atr->Pz());  
	    pm->time = atr->GetTime();

	    Track::vpPathMark_t& v = refs[track];
            v.push_back(pm);
	  }
	  else
	    throw(eH + "negative label for entry " + Form("%d",iTrackRef) + " in branch " + el->GetName()+ ".");
	}
      } // loop track refs 
      treeTR->SetBranchAddress(el->GetName(), 0);
    } // loop primaries, clones arrays
    treeTR->SetBranchStatus(Form("%s*", el->GetName()), 0);
  } // end loop through top branches


  // sort references and add it to tracks
  RenderElement::List_i  cit = cont->BeginChildren();
  while(cit != cont->EndChildren())
  {
    Track* track = dynamic_cast<Track*>(*cit);

    // add daughters path marks in the map 
    TParticle* p = stack->Particle(track->GetLabel());
    if(p->GetNDaughters()) {
      for(int d=p->GetFirstDaughter(); d>0 && d<=p->GetLastDaughter();++d) 
      {	
	TParticle* dp = stack->Particle(d);
	Reve::PathMark* pm = new PathMark( PathMark::Daughter);
	pm->V.Set(dp->Vx(),dp->Vy(), dp->Vz());
	pm->P.Set(dp->Px(),dp->Py(), dp->Pz()); 
	pm->time = dp->T();
	Track::vpPathMark_t& v = refs[track->GetLabel()];
	v.push_back(pm);
      }
    }
    
    map<Int_t, Track::vpPathMark_t > ::iterator mi = refs.find(track->GetLabel());
    if(mi != refs.end()) {
      Track::vpPathMark_t& v = mi->second;
      sort(v.begin(), v.end(), cmp_pathmark());
      for(Track::vpPathMark_i i=v.begin(); i!=v.end(); ++i){
	track->AddPathMark(*i);
      }
    }
    ++cit;
  }
}
