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
/// Implementation of an AliAODBranchReplicator to produce slim UPC aods.
///
/// This replicator is in charge of replicating the {tracks,vertices,
/// vzero, zdc, ad, and, optionally, the SPD tracklets} branches
/// of the standard AOD into UPC AODs.
///
/// The tracks are filtered so that only filter bits 0 and 1 make it to the output aods
///
///
/// \author Michal Broz

#include "AliAODUPCReplicator.h"

#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTZERO.h"
#include "AliAODTrack.h"
#include "AliAODVZERO.h"
#include "AliAnalysisCuts.h"
#include <cassert>

/// \cond CLASSIMP
ClassImp(AliAODUPCReplicator)
/// \endcond

//_____________________________________________________________________________
AliAODUPCReplicator::AliAODUPCReplicator(const char* name, const char* title,
                                           Bool_t replicateHeader,
                                           Bool_t replicateTracklets,
					   Bool_t replicateMuonTracks)
:AliAODBranchReplicator(name,title), 
fTracks(0x0), 
fVertices(0x0), 
fVZERO(0x0),
fAD(0x0),
fHeader(0x0),
fTracklets(0x0),
fZDC(0x0),
fList(0x0),
fTOFHeader(0x0),
fReplicateHeader(replicateHeader),
fReplicateTracklets(replicateTracklets),
fReplicateMuonTracks(replicateMuonTracks)
{

}

//_____________________________________________________________________________
AliAODUPCReplicator::~AliAODUPCReplicator()
{
  /// dtor
  delete fList;
}

//_____________________________________________________________________________
TList* AliAODUPCReplicator::GetList() const
{
  /// return (and build if not already done) our internal list of managed objects
  if (!fList)
  {
    if ( fReplicateHeader )
    {
      fHeader = new AliAODHeader;
      fTOFHeader = new AliTOFHeader;
    }
    
    if ( fReplicateTracklets )
    {
      fTracklets = new AliAODTracklets;
      fTracklets->SetName("tracklets");
    }

    fTracks = new TClonesArray("AliAODTrack",30);
    fTracks->SetName("tracks");
    
    fVertices = new TClonesArray("AliAODVertex",3);
    fVertices->SetName("vertices");     
      
    fVZERO = new AliAODVZERO;
    
    fZDC = new AliAODZDC;
    
    fAD = new AliAODAD;
    
    fList = new TList;
    fList->SetOwner(kTRUE);
    
    fList->Add(fHeader);
    fList->Add(fTracks);
    fList->Add(fVertices);
    if ( fReplicateTracklets ) fList->Add(fTracklets);
    fList->Add(fVZERO);
    fList->Add(fZDC);
    fList->Add(fAD);
    fList->Add(fTOFHeader);
  }
  return fList;
}

//_____________________________________________________________________________
void AliAODUPCReplicator::ReplicateAndFilter(const AliAODEvent& source)
{
  /// Replicate (and filter if filters are there) the relevant
  /// parts we're interested in AODEvent
  
  assert(fTracks!=0x0);
  
  if (fReplicateHeader)
  {
    AliAODHeader * header = dynamic_cast<AliAODHeader*>(source.GetHeader());
    if(!header) AliFatal("Not a standard AOD");
    *fHeader = *(header);
    
    AliTOFHeader * tofheader = (AliTOFHeader*)(source.GetTOFHeader());
    if(!tofheader) AliFatal("Not a standard AOD");
    *fTOFHeader = *(tofheader);
  }

  if (fReplicateTracklets)
  {
    *fTracklets = *(source.GetTracklets());
  }
  
  *fVZERO = *(source.GetVZEROData());

  *fZDC = *(source.GetZDCData());
  
  if ( source.GetADData() )
  {
    *fAD = *(source.GetADData());
  }
  
  fTracks->Clear("C");
  TIter next(source.GetTracks());
  AliAODTrack* trk;
  Int_t nGlobalTracks=0;
  Int_t nITSsaTracks=0;
  Int_t nMuonTracks=0;
  
  while (( trk = static_cast<AliAODTrack*>(next()) )) {

    if(IsGoodGlobalTrack(trk)) new ((*fTracks)[nGlobalTracks++]) AliAODTrack(*trk);
    if(IsGoodITSsaTrack(trk)) new ((*fTracks)[nITSsaTracks++]) AliAODTrack(*trk);
    if(IsGoodMuonTrack(trk)) new ((*fTracks)[nMuonTracks++]) AliAODTrack(*trk);

  }
  
  assert(fVertices!=0x0);
  fVertices->Clear("C");
  TIter nextV(source.GetVertices());
  AliAODVertex* v;
  Int_t nvertices(0);
  
  while ( ( v = static_cast<AliAODVertex*>(nextV()) ) ) {
    if (v->GetType() == AliAODVertex::kPrimary || v->GetType() == AliAODVertex::kMainSPD) {
      AliAODVertex* tmp = v->CloneWithoutRefs();
      AliAODVertex* copiedVertex = new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
      // to insure the main vertex retains the ncontributors information
      // which is otherwise computed dynamically from
      // references to tracks we set it here
      copiedVertex->SetNContributors(v->GetNContributors()); 
      delete tmp;
    }
  }
  //AliInfo(Form("%d good global tracks",nGlobalTracks));
  //AliInfo(Form("%d good ITSsa tracks",nITSsaTracks));

  
}

//_____________________________________________________________________________
Bool_t AliAODUPCReplicator::IsGoodGlobalTrack(const AliAODTrack *trk)
{

  if(!(trk->TestFilterBit(1<<0))) return kFALSE;  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAODUPCReplicator::IsGoodITSsaTrack(const AliAODTrack *trk)
{

  if(fReplicateMuonTracks){
    Bool_t goodMuonCaloTrack = kTRUE;
    if(!(trk->GetStatus() & AliAODTrack::kITSrefit)) goodMuonCaloTrack=kFALSE;
    if(trk->GetITSNcls() < 4)goodMuonCaloTrack=kFALSE;
    if(trk->Chi2perNDF() > 2.5)goodMuonCaloTrack=kFALSE;
    if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1)))goodMuonCaloTrack=kFALSE;
    
    if(goodMuonCaloTrack)return kTRUE;
    }

  if(!(trk->TestFilterBit(1<<1))) return kFALSE;
  return kTRUE;
} 

//_____________________________________________________________________________
Bool_t AliAODUPCReplicator::IsGoodMuonTrack(const AliAODTrack *trk)
{

  if(!fReplicateMuonTracks)return kFALSE;
  if(!trk->IsMuonTrack())return kFALSE;
  if( trk->GetRAtAbsorberEnd() < 17.5 || trk->GetRAtAbsorberEnd() > 89.5 ) return kFALSE;
  if( trk->Eta() < -4.0 || trk->Eta() > -2.5 ) return kFALSE;
  
  return kTRUE;
}   
      
