/**************************************************************************
* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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

///
/// Implementation of an AliAODBranchReplicator to produce slim muon aods.
///
/// This replicator is in charge of replicating the {tracks,vertices,
/// vzero, tzero, zdc, ad, and, optionally, the SPD tracklets} branches
/// of the standard AOD into muon AODs (AliAOD.Muons.root)
///
/// The tracks are filtered so that only muon tracks (and only muon tracks
/// that pass the trackCut if present) make it to the output aods
///
/// The vertices are filtered so that only the primary (and pileup) vertices make it
/// to the output aods.
///
/// \author L. Aphecetche (Subatech)

#include "AliAODMuonReplicator.h"

#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTZERO.h"
#include "AliAODTrack.h"
#include "AliAODVZERO.h"
#include "AliAnalysisCuts.h"
#include <cassert>

/// \cond CLASSIMP
ClassImp(AliAODMuonReplicator)
/// \endcond

//_____________________________________________________________________________
AliAODMuonReplicator::AliAODMuonReplicator(const char* name, const char* title,
                                           AliAnalysisCuts* trackCut,
                                           AliAnalysisCuts* vertexCut,
                                           Int_t mcMode,
                                           Bool_t replicateHeader,
                                           Bool_t replicateTracklets)
:AliAODBranchReplicator(name,title), 
fTrackCut(trackCut), fTracks(0x0), 
fVertexCut(vertexCut), fVertices(0x0), 
fVZERO(0x0),
fTZERO(0x0),
fAD(0x0),
fHeader(0x0),
fTracklets(0x0),
fZDC(0x0),
fList(0x0),
fMCParticles(0x0),
fMCHeader(0x0),
fMCMode(mcMode),
fLabelMap(),
fParticleSelected(),
fReplicateHeader(replicateHeader),
fReplicateTracklets(replicateTracklets)
{
  /// default ctor
  ///
  /// \param trackCut if present will filter out tracks
  /// \param vertexCut if present will filter out vertices
  /// \param mcMode what to do with MC information (0: skip it, 1: copy all,
  /// 2: copy only for events with at least one muon )
  /// \param replicateHeader whether or not to handle the replication of the
  /// AOD header branch
  /// \param replicateTracklets whether or not to include the SPD tracklets branch
  ///
}

//_____________________________________________________________________________
AliAODMuonReplicator::~AliAODMuonReplicator()
{
  /// dtor
  delete fTrackCut;
  delete fVertexCut;
  delete fList;
}

//_____________________________________________________________________________
void AliAODMuonReplicator::SelectParticle(Int_t i)
{
  /// taking the absolute values here, need to take care
  /// of negative daughter and mother
  /// IDs when setting!
  
  if (!IsParticleSelected(TMath::Abs(i)))
  {
    fParticleSelected.Add(TMath::Abs(i),1);    
  }
}

//_____________________________________________________________________________
Bool_t AliAODMuonReplicator::IsParticleSelected(Int_t i)  
{
  /// taking the absolute values here, need to take
  /// care with negative daughter and mother
  /// IDs when setting!
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
  /// Gets the label from the new created Map
  /// Call CreatLabelMap before
  /// otherwise only 0 returned
  if ( i < 0 ) {
    AliError(Form("Searching for new label of particle with invalid label %i",i));
    return i;
  }
  return fLabelMap.GetValue(i);
}

//_____________________________________________________________________________
void AliAODMuonReplicator::FilterMC(const AliAODEvent& source)
{
  /// Filter MC information

  fMCHeader->Reset();
  fMCParticles->Clear("C");

  AliAODMCHeader* mcHeader(0x0);
  TClonesArray* mcParticles(0x0);
  
  fParticleSelected.Delete();
  
  if ( fMCMode==2 && !fTracks->GetEntries() ) return;
  // for fMCMode==2 we only copy MC information for events where there's at least one muon track
    
  mcHeader = static_cast<AliAODMCHeader*>(source.FindListObject(AliAODMCHeader::StdBranchName()));
  
  if ( mcHeader ) 
  {
    *fMCHeader = *mcHeader;
  }
  
  mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  if ( mcParticles && mcParticles->GetLast() >= 0 && fMCMode>=2 )
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
    
    if ( mcParticles && fMCMode==3 )
    {
      // loop on MC muon tracks to find their ancestors
      for ( Int_t ipart=0; ipart<mcParticles->GetEntries(); ipart++ )
      {
        AliAODMCParticle* mcp = static_cast<AliAODMCParticle*>(mcParticles->UncheckedAt(ipart));
        if ( TMath::Abs(mcp->PdgCode()) != 13 ) continue;
        Int_t label = ipart;

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
        Int_t d0 =  p->GetDaughterLabel(0);
        Int_t d1 =  p->GetDaughterLabel(1);
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
      if ( t->GetLabel() >= 0 ) t->SetLabel(GetNewLabel(t->GetLabel()));
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
  /// return (and build if not already done) our internal list of managed objects
  if (!fList)
  {
    if ( fReplicateHeader )
    {
      fHeader = new AliAODHeader;
    }
    
    if ( fReplicateTracklets )
    {
      fTracklets = new AliAODTracklets;
      fTracklets->SetName("tracklets");
    }

    fTracks = new TClonesArray("AliAODTrack",30);
    fTracks->SetName("tracks");    
    
    fVertices = new TClonesArray("AliAODVertex",2);
    fVertices->SetName("vertices");    
    
    fVZERO = new AliAODVZERO;
    
    fTZERO = new AliAODTZERO;
    
    fZDC = new AliAODZDC;
    
    fAD = new AliAODAD;
    
    fList = new TList;
    fList->SetOwner(kTRUE);
    
    fList->Add(fHeader);
    fList->Add(fTracks);
    fList->Add(fVertices);
    fList->Add(fTracklets);
    fList->Add(fVZERO);
    fList->Add(fTZERO);
    fList->Add(fZDC);
    fList->Add(fAD);
    
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
  /// Replicate (and filter if filters are there) the relevant
  /// parts we're interested in AODEvent
  
  assert(fTracks!=0x0);
  
  if (fReplicateHeader)
  {
    AliAODHeader * header = dynamic_cast<AliAODHeader*>(source.GetHeader());
    if(!header) AliFatal("Not a standard AOD");
    *fHeader = *(header);
  }

  if (fReplicateTracklets)
  {
    *fTracklets = *(source.GetTracklets());
  }
  
  *fVZERO = *(source.GetVZEROData());

  *fTZERO = *(source.GetTZEROData());

  *fZDC = *(source.GetZDCData());
  
  if ( source.GetADData() )
  {
    *fAD = *(source.GetADData());
  }
  
  fTracks->Clear("C");
  TIter next(source.GetTracks());
  AliAODTrack* t;
  Int_t nMuons=0;
  Int_t inputMuons=0;
  
  while (( t = static_cast<AliAODTrack*>(next()) )) {

    if (t->IsMuonTrack() || t->IsMuonGlobalTrack()) ++inputMuons;    // just a counter: MUON and MUON+MFT tracks before track cuts are applied
     
    // MUON and MUON+MFT tracks selected                    // AU
    if (!fTrackCut || fTrackCut->IsSelected(t)) {
      new ((*fTracks)[nMuons++]) AliAODTrack(*t);
    }

  }
  
  assert(fVertices!=0x0);
  fVertices->Clear("C");
  TIter nextV(source.GetVertices());
  AliAODVertex* v;
  Int_t nvertices(0);
  
  while ( ( v = static_cast<AliAODVertex*>(nextV()) ) ) {
    if ( !fVertexCut || fVertexCut->IsSelected(v) ) {
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
  
  AliDebug(1,Form("input mu tracks=%d tracks=%d vertices=%d",
                  inputMuons,fTracks->GetEntries(),fVertices->GetEntries()));

  // Finally, deal with MC information, if needed

  if (fMCMode>0) FilterMC(source);

}

