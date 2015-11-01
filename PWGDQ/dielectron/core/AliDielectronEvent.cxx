/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                Dielectron Event                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TExMap.h>
#include <TProcessID.h>

#include <AliVTrack.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>

#include "AliDielectronEvent.h"

ClassImp(AliDielectronEvent)

AliDielectronEvent::AliDielectronEvent() :
  TNamed(),
  fArrTrackP(),
  fArrTrackN(),
  fArrVertex("AliAODVertex",10),
  fArrPairs("AliKFParticle",0),
  fNTracksP(0),
  fNTracksN(0),
  fIsAOD(kFALSE),
  fEventData(),
  fPID(0x0),
  fPIDIndex(0)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronEvent::AliDielectronEvent(const char* name, const char* title) :
  TNamed(name, title),
  fArrTrackP(),
  fArrTrackN(),
  fArrVertex("AliAODVertex",10),
  fArrPairs("AliKFParticle",0),
  fNTracksP(0),
  fNTracksN(0),
  fIsAOD(kFALSE),
  fEventData(),
  fPID(0x0),
  fPIDIndex(0)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronEvent::~AliDielectronEvent()
{
  //
  // Default Destructor
  //
  fArrTrackP.Delete();
  fArrTrackN.Delete();
  fArrVertex.Delete();
  fArrPairs.Delete();
}

//______________________________________________
void AliDielectronEvent::SetTracks(const TObjArray &arrP, const TObjArray &arrN, const TObjArray &/*arrPairs*/)
{
  //
  // Setup AliKFParticles
  // assumes that the objects in arrP and arrN are AliVTracks
  //

  //Clear out old entries before filling new ones
  Clear();
  // we keep the tracks buffered to minimise new / delete operations
  fNTracksN=0;
  fNTracksP=0;

  //check size of the arrays
  if (fArrTrackP.GetSize()<arrP.GetSize()) {
    fArrTrackP.Expand(arrP.GetSize());
  }
  if (fArrTrackN.GetSize()<arrN.GetSize()) {
    fArrTrackN.Expand(arrN.GetSize());
  }

  TExMap mapStoredVertices;
  fPIDIndex=TProcessID::GetPIDs()->IndexOf(fPID);
  // fill particles
  Int_t tracks=0;
  for (Int_t itrack=0; itrack<arrP.GetEntriesFast(); ++itrack){
    if (!fIsAOD){
      AliESDtrack *track=dynamic_cast<AliESDtrack*>(arrP.At(itrack));
      if (!track) continue;
      new (fArrTrackP[tracks]) AliESDtrack(*track);
      ++tracks;
    } else {
      AliAODTrack *track=dynamic_cast<AliAODTrack*>(arrP.At(itrack));
      if (!track) continue;

      // buffer track
      AliAODTrack   *ctrack = new (fArrTrackP[tracks]) AliAODTrack(*track);

      // buffer vertex, don't duplicate
      // most particles will be assiciated to the primary vertex ...
      AliAODVertex *vtx=track->GetProdVertex();
      AliAODVertex *cvertex = 0x0;
      if (vtx){
        cvertex = reinterpret_cast<AliAODVertex*>(mapStoredVertices.GetValue(reinterpret_cast<ULong64_t>(vtx)));
        if (!cvertex) {
          if (mapStoredVertices.Capacity()<=mapStoredVertices.GetSize()) mapStoredVertices.Expand(2*mapStoredVertices.GetSize());
          if (fArrVertex.GetSize()<=fArrVertex.GetEntriesFast()) fArrVertex.Expand(2*fArrVertex.GetSize());
          cvertex = new (fArrVertex[fArrVertex.GetEntriesFast()]) AliAODVertex(*vtx);
          AssignID(cvertex);
          mapStoredVertices.Add(reinterpret_cast<ULong64_t>(vtx),reinterpret_cast<ULong64_t>(cvertex));
        }
      }
      ctrack->SetProdVertex(cvertex);
      ++tracks;
    }
  }
  fNTracksP=tracks;

  tracks=0;
  for (Int_t itrack=0; itrack<arrN.GetEntriesFast(); ++itrack){
    if (!fIsAOD){
      AliESDtrack *track=dynamic_cast<AliESDtrack*>(arrN.At(itrack));
      if (!track) continue;
      new (fArrTrackN[tracks]) AliESDtrack(*track);
      ++tracks;
    } else {
      AliAODTrack *track=dynamic_cast<AliAODTrack*>(arrN.At(itrack));
      if (!track) continue;

      //buffer track
      AliAODTrack   *ctrack = new (fArrTrackN[tracks]) AliAODTrack(*track);
      
      // buffer vertex, don't duplicate
      // most particles will be assiciated to the primary vertex ...
      AliAODVertex *vtx=track->GetProdVertex();
      AliAODVertex *cvertex = 0x0;
      if (vtx){
        cvertex = reinterpret_cast<AliAODVertex*>(mapStoredVertices.GetValue(reinterpret_cast<ULong64_t>(vtx)));
        if (!cvertex) {
          if (mapStoredVertices.Capacity()<=mapStoredVertices.GetSize()) mapStoredVertices.Expand(2*mapStoredVertices.GetSize());
          if (fArrVertex.GetSize()<=fArrVertex.GetEntriesFast()) fArrVertex.Expand(2*fArrVertex.GetSize());
          cvertex = new (fArrVertex[fArrVertex.GetEntriesFast()]) AliAODVertex(*vtx);
          AssignID(cvertex);
          mapStoredVertices.Add(reinterpret_cast<ULong64_t>(vtx),reinterpret_cast<ULong64_t>(cvertex));
        }
      }
      ctrack->SetProdVertex(cvertex);
      ++tracks;
    }
  }
  fNTracksN=tracks;

  //TODO: pair arrays
}

//______________________________________________
void AliDielectronEvent::Clear(Option_t *opt)
{
  //
  // clear arrays
  //
//   fArrTrackP.Clear(opt);
//   fArrTrackN.Clear(opt);

  for (Int_t i=fArrTrackP.GetEntriesFast()-1; i>=0; --i){
    delete fArrTrackP.RemoveAt(i);
  }
  
  for (Int_t i=fArrTrackN.GetEntriesFast()-1; i>=0; --i){
    delete fArrTrackN.RemoveAt(i);
  }

  for (Int_t i=0; i<fArrVertex.GetEntriesFast(); ++i){
    delete fArrVertex.RemoveAt(i);
  }

  fArrPairs.Clear(opt);
}

//______________________________________________
void AliDielectronEvent::SetAOD(Int_t sizeP, Int_t sizeN)
{
  //
  // use AOD as input
  //

  //overwrite with fixed sizes
  sizeP=sizeN=1000;
  fArrTrackP.SetClass("AliAODTrack",sizeP);
  fArrTrackN.SetClass("AliAODTrack",sizeN);
  fIsAOD=kTRUE;
}

//______________________________________________
void AliDielectronEvent::SetESD(Int_t sizeP, Int_t sizeN)
{
  //
  // use ESD as input
  //

  //overwrite with fixed sizes
  sizeP=sizeN=1000;
  fArrTrackP.SetClass("AliESDtrack",sizeP);
  fArrTrackN.SetClass("AliESDtrack",sizeN);
  fIsAOD=kFALSE;
}

//______________________________________________
void AliDielectronEvent::SetEventData(const Double_t data[AliDielectronVarManager::kNMaxValues])
{
  //
  // copy only evnet variables
  //
  for (Int_t i=AliDielectronVarManager::kPairMax; i<AliDielectronVarManager::kNMaxValues;++i) fEventData[i]=data[i];
}

//______________________________________________
void AliDielectronEvent::AssignID(TObject *obj)
{
  //
  // Custom function to assign a uid to an object with an own process id
  // to avoid problems buffering the vertices
  //
  UInt_t uid=1;
  if (fPID->GetObjects()) uid=fPID->GetObjects()->GetEntriesFast();
  uid+=(fPIDIndex<<24);
  obj->SetBit(kIsReferenced);
  obj->SetUniqueID(uid);
  fPID->PutObjectWithID(obj);
}
