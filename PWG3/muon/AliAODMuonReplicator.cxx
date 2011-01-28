/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// $Id$

#include "AliAODMuonReplicator.h"
#include "AliAnalysisCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include <cassert>

ClassImp(AliAODMuonReplicator)

//_____________________________________________________________________________
AliAODMuonReplicator::AliAODMuonReplicator(const char* name, const char* title,
                                                 AliAnalysisCuts* trackCut,
                                                 AliAnalysisCuts* vertexCut)
: AliAODBranchReplicator(name,title), 
fTrackCut(trackCut), fTracks(0x0), 
fVertexCut(vertexCut), fVertices(0x0),
fList(0x0)
{
  // default ctor
}

//_____________________________________________________________________________
AliAODMuonReplicator::~AliAODMuonReplicator()
{
  // dtor
  delete fTrackCut;
  delete fVertexCut;
  delete fList;
}

//_____________________________________________________________________________
TList* AliAODMuonReplicator::GetList() const
{
  // return (and build if not already done) our internal list of managed objects
  if (!fList)
  {
    fList = new TList;
    fList->SetOwner(kTRUE);
    fTracks = new TClonesArray("AliAODTrack",100);
		fTracks->SetName("tracks");    
    fVertices = new TClonesArray("AliAODVertex",10);
		fVertices->SetName("vertices");    
    fList->Add(fTracks);
    fList->Add(fVertices);
  }
  return fList;
}

//_____________________________________________________________________________
void AliAODMuonReplicator::ReplicateAndFilter(const AliAODEvent& source)
{
  // Replicate (and filter if filters are there) the relevant parts we're interested in AODEvent
  
  static int n(0);
  
  ++n;
  
  assert(fTracks!=0x0);
  fTracks->Clear("C");
  TIter next(source.GetTracks());
  AliAODTrack* t;
  Int_t ntracks(0);
  
  
  while ( ( t = static_cast<AliAODTrack*>(next()) ) )
  {
    if ( !fTrackCut || fTrackCut->IsSelected(t) ) 
    {
      new((*fTracks)[ntracks++]) AliAODTrack(*t);
    }
  }
  
  assert(fVertices!=0x0);
  fVertices->Clear("C");
  TIter nextV(source.GetVertices());
  AliAODVertex* v;
  Int_t nvertices(0);
  
  while ( ( v = static_cast<AliAODVertex*>(nextV()) ) )
  {
    if ( !fVertexCut || fVertexCut->IsSelected(v) ) 
    {
      AliAODVertex* tmp = v->CloneWithoutRefs();
      new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
      delete tmp;
    }
  }
  
  AliInfo(Form("n=%d tracks=%d vertices=%d",n,fTracks->GetEntries(),fVertices->GetEntries()));
}

