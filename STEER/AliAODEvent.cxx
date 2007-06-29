/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TTree.h>

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"

ClassImp(AliAODEvent)

//______________________________________________________________________________
AliAODEvent::AliAODEvent() :
  fAODObjects(new TList()),
  fHeader(0),
  fTracks(0),
  fVertices(0),
  fClusters(0),
  fJets(0)
{
  // default constructor
}

//______________________________________________________________________________
AliAODEvent::~AliAODEvent() 
{
// destructor
    delete fAODObjects;
}

//______________________________________________________________________________
void AliAODEvent::AddObject(TObject* obj) 
{
  // Add an object to the list of object.
  // Please be aware that in order to increase performance you should
  // refrain from using TObjArrays (if possible). Use TClonesArrays, instead.
  
  fAODObjects->AddLast(obj);
}

//______________________________________________________________________________
TObject *AliAODEvent::GetObject(const char *objName) const 
{
  // Return the pointer to the object with the given name.

  return fAODObjects->FindObject(objName);
}

//______________________________________________________________________________
void AliAODEvent::CreateStdContent() 
{
  // create the standard AOD content and set pointers

  // create standard objects and add them to the TList of objects
  AddObject(new AliAODHeader());
  AddObject(new TClonesArray("AliAODTrack", 0));
  AddObject(new TClonesArray("AliAODVertex", 0));
  AddObject(new TClonesArray("AliAODCluster", 0));
  AddObject(new TClonesArray("AliAODJet", 0));

  // read back pointers
  GetStdContent();

  // set names
  fTracks->SetName("tracks");
  fVertices->SetName("vertices");
  fClusters->SetName("clusters");
  fJets->SetName("jets");

}

//______________________________________________________________________________
void AliAODEvent::GetStdContent()
{
  // set pointers for standard content

  fHeader   = (AliAODHeader*)fAODObjects->At(0);
  fTracks   = (TClonesArray*)fAODObjects->At(1);
  fVertices = (TClonesArray*)fAODObjects->At(2);
  fClusters = (TClonesArray*)fAODObjects->At(3);
  fJets     = (TClonesArray*)fAODObjects->At(4);
}

//______________________________________________________________________________
void AliAODEvent::ResetStd(Int_t trkArrSize, Int_t vtxArrSize)
{
  // deletes content of standard arrays and resets size
  fTracks->Delete();
  if (trkArrSize > fTracks->GetSize()) 
    fTracks->Expand(trkArrSize);

  fVertices->Delete();
  if (vtxArrSize > fVertices->GetSize()) 
    fVertices->Expand(vtxArrSize);
}

void AliAODEvent::ClearStd()
{
  // clears the standard arrays
    fTracks   ->Clear();
    fVertices ->Clear();
    fClusters ->Clear();
    fJets     ->Clear();
}

//______________________________________________________________________________
Int_t AliAODEvent::GetMuonTracks(TRefArray *muonTracks) const
{
  // fills the provided TRefArray with all found muon tracks

  muonTracks->Clear();

  AliAODTrack *track = 0;
  for (Int_t iTrack = 0; iTrack < GetNTracks(); iTrack++) {
    if ((track = GetTrack(iTrack))->GetMostProbablePID() == AliAODTrack::kMuon) {
      muonTracks->Add(track);
    }
  }
  
  return muonTracks->GetSize();
}



void AliAODEvent::ReadFromTree(TTree *tree)
{
    // connects aod event to tree

    fAODObjects = (TList*)((AliAODEvent*)tree->GetTree()->GetUserInfo()->FindObject("AliAODEvent"))->GetList(); 
    TIter next(fAODObjects);
    TNamed *el;
    while((el=(TNamed*)next())){
	TString bname(el->GetName());
	tree->SetBranchAddress(bname.Data(),fAODObjects->GetObjectRef(el));
    }
    GetStdContent();
}

