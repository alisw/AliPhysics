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

#include <AliVTrack.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>

#include "AliDielectronEvent.h"

ClassImp(AliDielectronEvent)

AliDielectronEvent::AliDielectronEvent() :
  TNamed(),
  fArrTrackP(),
  fArrTrackN(),
  fArrVertexP("AliAODVertex",1000),
  fArrVertexN("AliAODVertex",1000),
  fArrPairs("AliKFParticle",0),
  fNTracksP(0),
  fNTracksN(0),
  fIsAOD(kFALSE),
  fEventData()
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
  fArrVertexP("AliAODVertex",1000),
  fArrVertexN("AliAODVertex",1000),
  fArrPairs("AliKFParticle",0),
  fNTracksP(0),
  fNTracksN(0),
  fIsAOD(kFALSE),
  fEventData()
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
  fArrVertexP.Delete();
  fArrVertexN.Delete();
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
    fArrVertexP.Expand(arrP.GetSize());
  }
  if (fArrTrackN.GetSize()<arrN.GetSize()) {
    fArrTrackN.Expand(arrN.GetSize());
    fArrVertexN.Expand(arrN.GetSize());
  }

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

      new (fArrTrackP[tracks]) AliAODTrack(*track);
      new (fArrVertexP[tracks]) AliAODVertex(*(track->GetProdVertex()));
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
      
      new (fArrTrackN[tracks]) AliAODTrack(*track);
      new (fArrVertexN[tracks]) AliAODVertex(*(track->GetProdVertex()));
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
    delete fArrVertexP.RemoveAt(i);
  }
  
  for (Int_t i=fArrTrackN.GetEntriesFast()-1; i>=0; --i){
    delete fArrTrackN.RemoveAt(i);
    delete fArrVertexN.RemoveAt(i);
  }
  
  fArrPairs.Clear(opt);
  
}

//______________________________________________
void AliDielectronEvent::SetAOD()
{
  //
  // use AOD as input
  //
  fArrTrackP.SetClass("AliAODTrack",1000);
  fArrTrackN.SetClass("AliAODTrack",1000);
  fIsAOD=kTRUE;
}

//______________________________________________
void AliDielectronEvent::SetESD()
{
  //
  // use ESD as input
  //
  fArrTrackP.SetClass("AliESDtrack",1000);
  fArrTrackN.SetClass("AliESDtrack",1000);
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
