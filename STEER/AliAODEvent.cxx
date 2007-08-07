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

// definition of std AOD member names
  const char* AliAODEvent::fAODListName[kAODListN] = {"header",
						      "tracks",
						      "vertices",
						      "clusters",
						      "jets",
						      "tracklets"};
//______________________________________________________________________________
AliAODEvent::AliAODEvent() :
  AliVEvent(),
  fAODObjects(new TList()),
  fHeader(0),
  fTracks(0),
  fVertices(0),
  fClusters(0),
  fJets(0),
  fTracklets(0)
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
TObject *AliAODEvent::FindListObject(const char *objName)
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
  AddObject(new AliAODTracklets());

  // set names
  SetStdNames();

  // read back pointers
  GetStdContent();

  return;
}

//______________________________________________________________________________
void AliAODEvent::SetStdNames()
{
  // introduce the standard naming

  if(fAODObjects->GetEntries()==kAODListN){
    for(int i = 0;i < fAODObjects->GetEntries();i++){
      TObject *fObj = fAODObjects->At(i);
      if(fObj->InheritsFrom("TNamed")){
	((TNamed*)fObj)->SetName(fAODListName[i]);
      }
      else if(fObj->InheritsFrom("TClonesArray")){
	((TClonesArray*)fObj)->SetName(fAODListName[i]);
      }
    }
  }
  else{
    printf("%s:%d SetStdNames() Wrong number of Std Entries \n",(char*)__FILE__,__LINE__);
  }
} 

//______________________________________________________________________________
void AliAODEvent::GetStdContent()
{
  // set pointers for standard content

  fHeader    = (AliAODHeader*)fAODObjects->FindObject("header");
  fTracks    = (TClonesArray*)fAODObjects->FindObject("tracks");
  fVertices  = (TClonesArray*)fAODObjects->FindObject("vertices");
  fClusters  = (TClonesArray*)fAODObjects->FindObject("clusters");
  fJets      = (TClonesArray*)fAODObjects->FindObject("jets");
  fTracklets = (AliAODTracklets*)fAODObjects->FindObject("tracklets");
  
  /* old implementation
  fHeader    = (AliAODHeader*)fAODObjects->At(0);
  fTracks    = (TClonesArray*)fAODObjects->At(1);
  fVertices  = (TClonesArray*)fAODObjects->At(2);
  fClusters  = (TClonesArray*)fAODObjects->At(3);
  fJets      = (TClonesArray*)fAODObjects->At(4);
  fTracklets = (AliAODTracklets*)fAODObjects->At(5);
  */
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
    fTracklets->DeleteContainer();
}

//______________________________________________________________________________
Int_t AliAODEvent::GetMuonTracks(TRefArray *muonTracks) const
{
  // fills the provided TRefArray with all found muon tracks

  muonTracks->Clear();

  AliAODTrack *track = 0;
  for (Int_t iTrack = 0; iTrack < GetNTracks(); iTrack++) {
    if ((track = GetTrack(iTrack))->IsMuonTrack()) {
      muonTracks->Add(track);
    }
  }
  
  return muonTracks->GetSize();
}


void AliAODEvent::ReadFromTree(TTree *tree)
{
  // connects aod event to tree
  
  // load the TTree
  tree->LoadTree(0);

  fAODObjects = (TList*)((AliAODEvent*)tree->GetTree()->GetUserInfo()->FindObject("AliAODEvent"))->GetList(); 
  TIter next(fAODObjects);
  TNamed *el;
  while((el=(TNamed*)next())){
    TString bname(el->GetName());
    tree->SetBranchAddress(bname.Data(),fAODObjects->GetObjectRef(el));
  }
  GetStdContent();
}

//______________________________________________________________________________
void AliAODEvent::Print(Option_t *) const
{
  // Something meaningful should be implemented here.
  
  return;
}
