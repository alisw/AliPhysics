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
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include <cassert>

ClassImp(AliAODMuonReplicator)

//_____________________________________________________________________________
AliAODMuonReplicator::AliAODMuonReplicator(const char* name, const char* title,
                                           AliAnalysisCuts* trackCut,
                                           AliAnalysisCuts* vertexCut,
                                           Int_t mcMode)
:AliAODBranchReplicator(name,title), 
fTrackCut(trackCut), fTracks(0x0), 
fVertexCut(vertexCut), fVertices(0x0), 
fDimuons(0x0),
fList(0x0),
fMCParticles(0x0),
fMCHeader(0x0),
fMCMode(mcMode),
fLabelMap(),
fParticleSelected()
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
void AliAODMuonReplicator::SelectParticle(Int_t i)
{
  // taking the absolute values here, need to take care 
  // of negative daughter and mother
  // IDs when setting!
  
  if (!IsParticleSelected(TMath::Abs(i)))
  {
    fParticleSelected.Add(TMath::Abs(i),1);    
  }
}

//_____________________________________________________________________________
Bool_t AliAODMuonReplicator::IsParticleSelected(Int_t i)  
{
  // taking the absolute values here, need to take 
  // care with negative daughter and mother
  // IDs when setting!
  return (fParticleSelected.GetValue(TMath::Abs(i))==1);
}


//_____________________________________________________________________________
void AliAODMuonReplicator::CreateLabelMap(const AliAODEvent& source)
{  
  //
  // this should be called once all selections are done 
  //
  
  fLabelMap.Delete();
  
  TClonesArray* mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  Int_t i(0);
  Int_t j(0);
  
  TIter next(mcParticles);
  
  while ( next() ) 
  {
    if (IsParticleSelected(i))
    {
      fLabelMap.Add(i,j++);
    }
    ++i;
  }
}

//_____________________________________________________________________________
Int_t AliAODMuonReplicator::GetNewLabel(Int_t i) 
{
  // Gets the label from the new created Map
  // Call CreatLabelMap before
  // otherwise only 0 returned
  return fLabelMap.GetValue(TMath::Abs(i));
}

//_____________________________________________________________________________
void AliAODMuonReplicator::FilterMC(const AliAODEvent& source)
{
  // Filter MC information

  fMCHeader->Reset();
  fMCParticles->Clear("C");

  AliAODMCHeader* mcHeader(0x0);
  TClonesArray* mcParticles(0x0);
  
  fParticleSelected.Delete();
  
  if ( fMCMode>=2 && !fTracks->GetEntries() ) return;
  // for fMCMode==1 we only copy MC information for events where there's at least one muon track
    
  mcHeader = static_cast<AliAODMCHeader*>(source.FindListObject(AliAODMCHeader::StdBranchName()));
  
  if ( mcHeader ) 
  {
    *fMCHeader = *mcHeader;
  }
  
  mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  if ( mcParticles && fMCMode>=2 )
  {
    // loop on (kept) muon tracks to find their ancestors
    TIter nextMT(fTracks);
    AliAODTrack* mt;
    
    while ( ( mt = static_cast<AliAODTrack*>(nextMT()) ) )
    {
      Int_t label = mt->GetLabel();
      
      while ( label >= 0 ) 
      {
        SelectParticle(label);
        AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(mcParticles->UncheckedAt(label));
        if (!mother)
        {
          AliError("Got a null mother ! Check that !");
          label = -1;
        }
        else
        {
          label = mother->GetMother();
        }
      }
    }
    
    CreateLabelMap(source);
    
    // Actual filtering and label remapping (shamelessly taken for the implementation of AliAODHandler::StoreMCParticles)
    TIter nextMC(mcParticles);
    AliAODMCParticle* p;
    Int_t nmc(0);
    Int_t nmcout(0);
    
    while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
    {
      AliAODMCParticle c(*p);
      
      if ( IsParticleSelected(nmc) )
      {
        // 
        Int_t d0 =  p->GetDaughter(0);
        Int_t d1 =  p->GetDaughter(1);
        Int_t m =   p->GetMother();
        
        // other than for the track labels, negative values mean
        // no daughter/mother so preserve it
        
        if(d0<0 && d1<0)
        {
          // no first daughter -> no second daughter
          // nothing to be done
          // second condition not needed just for sanity check at the end
          c.SetDaughter(0,d0);
          c.SetDaughter(1,d1);
        } else if(d1 < 0 && d0 >= 0) 
        {
          // Only one daughter
          // second condition not needed just for sanity check at the end
          if(IsParticleSelected(d0))
          {
            c.SetDaughter(0,GetNewLabel(d0));
          } else 
          {
            c.SetDaughter(0,-1);
          }
          c.SetDaughter(1,d1);
        }
        else if (d0 > 0 && d1 > 0 )
        {
          // we have two or more daughters loop on the stack to see if they are
          // selected
          Int_t d0tmp = -1;
          Int_t d1tmp = -1;
          for (int id = d0; id<=d1;++id)
          {
            if (IsParticleSelected(id))
            {
              if(d0tmp==-1)
              {
                // first time
                d0tmp = GetNewLabel(id);
                d1tmp = d0tmp; // this is to have the same schema as on the stack i.e. with one daugther d0 and d1 are the same 
              }
              else d1tmp = GetNewLabel(id);
            }
          }
          c.SetDaughter(0,d0tmp);
          c.SetDaughter(1,d1tmp);
        } else 
        {
          AliError(Form("Unxpected indices %d %d",d0,d1));
        }
        
        if ( m < 0 )
        {
          c.SetMother(m);
        } else 
        {
          if (IsParticleSelected(m)) 
          {
            c.SetMother(GetNewLabel(m));              
          }
          else 
          {
            AliError(Form("PROBLEM Mother not selected %d", m));              
          }
        }
        
        new ((*fMCParticles)[nmcout++]) AliAODMCParticle(c);
      }
      
      ++nmc;        
    }      
    
    // now remap the tracks...
    
    TIter nextTrack(fTracks);
    AliAODTrack* t;
    
    while ( ( t = static_cast<AliAODTrack*>(nextTrack()) ) )
    {
      t->SetLabel(GetNewLabel(t->GetLabel()));
    }
    
  }
  else if ( mcParticles ) 
  {
    // simple copy of input MC particles to ouput MC particles
    TIter nextMC(mcParticles);
    AliAODMCParticle* p;
    Int_t nmcout(0);
    
    while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
    {
      new ((*fMCParticles)[nmcout++]) AliAODMCParticle(*p);
    }
  }
  
  AliDebug(1,Form("input mc %d output mc %d",
                  mcParticles ? mcParticles->GetEntries() : 0,
                  fMCParticles ? fMCParticles->GetEntries() : 0));
  
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
    
    if ( fMCMode > 0 )
    {
      fMCHeader = new AliAODMCHeader;    
      fMCParticles = new TClonesArray("AliAODMCParticle",1000);
      fMCParticles->SetName(AliAODMCParticle::StdBranchName());
      fList->Add(fMCHeader);
      fList->Add(fMCParticles);
    }
  }
  return fList;
}

//_____________________________________________________________________________
void AliAODMuonReplicator::ReplicateAndFilter(const AliAODEvent& source)
{
  // Replicate (and filter if filters are there) the relevant parts we're interested in AODEvent
  
  assert(fTracks!=0x0);
  
  fTracks->Clear("C");
  TIter next(source.GetTracks());
  AliAODTrack* t;
  Int_t ntracks(0);
  Int_t input(0);
  
  while ( ( t = static_cast<AliAODTrack*>(next()) ) )
  {
    if ( t->IsMuonTrack() ) 
    {
      ++input;
    }
    
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
      AliAODVertex* copiedVertex = new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
      // to insure the main vertex retains the ncontributors information
      // (which is otherwise computed dynamically from
      // references to tracks, which we do not keep in muon aods...)
      // we set it here
      copiedVertex->SetNContributors(v->GetNContributors()); 
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

  AliDebug(1,Form("input mu tracks=%d tracks=%d vertices=%d ndimuons=%d",
                  input,fTracks->GetEntries(),fVertices->GetEntries(),fDimuons->GetEntries()));

  // Finally, deal with MC information, if needed

  if ( fMCMode > 0 )
  {
    FilterMC(source);      
  }
}

