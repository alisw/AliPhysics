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

//
// Implementation of a branch replicator 
// to produce slim muon and dimuon aods.
//
// This replicator is in charge of replicating the tracks,vertices and dimuons
// branches of the standard AOD into muon AODs (AliAOD.Muons.root and
// AliAOD.Dimuons.root)
// 
// The tracks are filtered so that only muon tracks (and only muon tracks
// that pass the trackCut if present) make it to the output aods
//
// The vertices are filtered so that only the primary vertices make it
// to the output aods.
//
// The dimuons are recreated here, according to the set of tracks
// that pass the trackCut (that set may be the same as the input set,
// but to be 100% safe, we always recreate the dimuons).
// 
// Author: L. Aphecetche (Subatech)

#include "AliAODMuonReplicator.h"
#include "AliAnalysisCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODDimuon.h"
#include <cassert>

ClassImp(AliAODMuonReplicator)

//_____________________________________________________________________________
AliAODMuonReplicator::AliAODMuonReplicator(const char* name, const char* title,
                                                 AliAnalysisCuts* trackCut,
                                                 AliAnalysisCuts* vertexCut)
: AliAODBranchReplicator(name,title), 
fTrackCut(trackCut), fTracks(0x0), 
fVertexCut(vertexCut), fVertices(0x0), 
fDimuons(0x0),
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
    fTracks = new TClonesArray("AliAODTrack",30);
		fTracks->SetName("tracks");    
    
    fVertices = new TClonesArray("AliAODVertex",2);
		fVertices->SetName("vertices");    

    fDimuons = new TClonesArray("AliAODDimuon",2);
    fDimuons->SetName("dimuons");
    
    fList = new TList;
    fList->SetOwner(kTRUE);
    
    fList->Add(fTracks);
    fList->Add(fVertices);
    fList->Add(fDimuons);
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
  
  fDimuons->Clear("C");
  
  // as there might be a track cut (different from the one of the main filtering),
  // we must recreate the dimuon completely from scratch to be 100% safe...

  Int_t ndimuons(0);

  for ( Int_t i = 0; i < ntracks; ++i ) 
  {
    for ( Int_t j = i+1; j < ntracks; ++j ) 
    {
      new((*fDimuons)[ndimuons++]) AliAODDimuon(fTracks->At(i),fTracks->At(j));
    }
  }
  
  AliDebug(1,Form("n=%d tracks=%d vertices=%d ndimuons=%d",n,fTracks->GetEntries(),fVertices->GetEntries(),fDimuons->GetEntries()));
}

