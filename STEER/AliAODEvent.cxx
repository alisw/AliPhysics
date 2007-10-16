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
						      "v0s",
						      "tracklets",
						      "jets",
						      "caloCells",
						      "caloClusters",
						      "fmdClusters",
						      "pmdClusters"
};
//______________________________________________________________________________
AliAODEvent::AliAODEvent() :
  AliVEvent(),
  fAODObjects(new TList()),
  fHeader(0),
  fTracks(0),
  fVertices(0),
  fV0s(0),
  fTracklets(0),
  fJets(0),
  fCaloCells(0),
  fCaloClusters(0),
  fFmdClusters(0),
  fPmdClusters(0)
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
  AddObject(new TClonesArray("AliAODv0", 0));
  AddObject(new AliAODTracklets());
  AddObject(new TClonesArray("AliAODJet", 0));
  AddObject(new AliAODCaloCells());
  AddObject(new TClonesArray("AliAODCaloCluster", 0));
  AddObject(new TClonesArray("AliAODFmdCluster", 0));
  AddObject(new TClonesArray("AliAODPmdCluster", 0));
  
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

  fHeader        = (AliAODHeader*)fAODObjects->FindObject("header");
  fTracks        = (TClonesArray*)fAODObjects->FindObject("tracks");
  fVertices      = (TClonesArray*)fAODObjects->FindObject("vertices");
  fV0s           = (TClonesArray*)fAODObjects->FindObject("v0s");
  fTracklets     = (AliAODTracklets*)fAODObjects->FindObject("tracklets");
  fJets          = (TClonesArray*)fAODObjects->FindObject("jets");
  fCaloCells     = (AliAODCaloCells*)fAODObjects->FindObject("caloCells");
  fCaloClusters  = (TClonesArray*)fAODObjects->FindObject("caloClusters");
  fFmdClusters   = (TClonesArray*)fAODObjects->FindObject("fmdClusters");
  fPmdClusters   = (TClonesArray*)fAODObjects->FindObject("pmdClusters");
}

//______________________________________________________________________________
void AliAODEvent::ResetStd(Int_t trkArrSize, 
			   Int_t vtxArrSize, 
			   Int_t v0ArrSize, 
			   Int_t jetSize, 
			   Int_t caloClusSize, 
			   Int_t fmdClusSize, 
			   Int_t pmdClusSize)
{
  // deletes content of standard arrays and resets size 

  fTracks->Delete();
  if (trkArrSize > fTracks->GetSize()) 
    fTracks->Expand(trkArrSize);

  fVertices->Delete();
  if (vtxArrSize > fVertices->GetSize()) 
    fVertices->Expand(vtxArrSize);
 
  fV0s->Delete();
  if (v0ArrSize > fV0s->GetSize()) 
    fV0s->Expand(v0ArrSize);

  fJets->Delete();
  if (jetSize > fJets->GetSize()) 
    fJets->Expand(jetSize);

  fCaloClusters->Delete();
  if (caloClusSize > fCaloClusters->GetSize()) 
    fCaloClusters->Expand(caloClusSize);

  fFmdClusters->Delete();
  if (fmdClusSize > fFmdClusters->GetSize()) 
    fFmdClusters->Expand(fmdClusSize);

  fPmdClusters->Delete();
  if (pmdClusSize > fPmdClusters->GetSize()) 
    fPmdClusters->Expand(pmdClusSize);
}

void AliAODEvent::ClearStd()
{
  // clears the standard arrays
  fTracks        ->Clear();
  fVertices      ->Clear();
  fV0s           ->Clear();
  fTracklets     ->DeleteContainer();
  fJets          ->Clear();
  fCaloCells     ->DeleteContainer();
  fCaloClusters ->Clear();
  fFmdClusters  ->Clear();
  fPmdClusters  ->Clear();
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
  
  return muonTracks->GetEntriesFast();
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
