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
//   Cut class providing cuts to all infomation                          //
//     available for the AliVParticle interface                          //
//                                                                       //
// Authors:                                                              //
//   Julian Book <Julian.Book@cern.ch>                                  //
/*



*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "AliDielectronV0Cuts.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTrackCuts.h"
#include "AliESDv0.h"

ClassImp(AliDielectronV0Cuts)


AliDielectronV0Cuts::AliDielectronV0Cuts() :
  AliDielectronVarCuts(),
  fV0TrackArr(0),
  fExcludeTracks(kTRUE),
  fMotherPdg(0),
  fNegPdg(0),
  fPosPdg(0),
  fOrbit(0),
  fPeriod(0),
  fBunchCross(0)
{
  //
  // Default costructor
  //
}

//________________________________________________________________________
AliDielectronV0Cuts::AliDielectronV0Cuts(const char* name, const char* title) :
  AliDielectronVarCuts(name,title),
  fV0TrackArr(0),
  fExcludeTracks(kTRUE),
  fMotherPdg(0),
  fNegPdg(0),
  fPosPdg(0),
  fOrbit(0),
  fPeriod(0),
  fBunchCross(0)
{
  //
  // Named contructor
  //
}

//________________________________________________________________________
AliDielectronV0Cuts::~AliDielectronV0Cuts()
{
  //
  // Destructor
  //

}

//________________________________________________________________________
void AliDielectronV0Cuts::InitEvent(AliVTrack *trk)
{
  //
  // Init the V0 candidates
  //

  // take current event from the track
  // TODO: this should be simplyfied by AliVTrack::GetEvent() as soon as implemented
  const AliVEvent *ev=0;
  if(trk->IsA() == AliAODTrack::Class())
    ev=static_cast<const AliVEvent*>((static_cast<const AliAODTrack*>(trk))->GetAODEvent());
  else if(trk->IsA() == AliESDtrack::Class())
    ev=static_cast<const AliVEvent*>((static_cast<const AliESDtrack*>(trk))->GetESDEvent());
  else
    return;


  // IsNewEvent
  if(!ev) return;
  if(!IsNewEvent(ev)) return;
  //  printf(" Build V0 candidates according to the applied cuts \n");

  // TODO think about MCevent
  //  Print();

  // rest booleans
  fV0TrackArr.ResetAllBits();

  // basic quality cut, at least one of the V0 daughters has to fullfill
  AliDielectronVarCuts dauQAcuts1;
  dauQAcuts1.AddCut(AliDielectronVarManager::kPt,           0.3,  1e30);
  dauQAcuts1.AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  dauQAcuts1.AddCut(AliDielectronVarManager::kNclsTPC,     50.0, 160.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  AliDielectronTrackCuts dauQAcuts2;
  //  dauQAcuts2.SetRequireITSRefit(kTRUE);
  dauQAcuts2.SetRequireTPCRefit(kTRUE);

  Int_t nV0s = 0;
  AliDielectronPair candidate;
  candidate.SetPdgCode(fMotherPdg);

  // ESD or AOD event
  if(ev->IsA() == AliESDEvent::Class()) {
    const AliESDEvent *esdev = static_cast<const AliESDEvent*>(ev);

    // loop over V0s
    for (Int_t iv=0; iv<esdev->GetNumberOfV0s(); ++iv){
      AliESDv0 *v = esdev->GetV0(iv);
      if(!v) continue;

      // should we make use of AliESDv0Cuts::GetPdgCode() to preselect candiadtes, e.g.:
      // if(fMotherPdg!=v->GetPdgCode()) continue;

      AliESDtrack *trNeg=esdev->GetTrack(v->GetIndex(0));
      AliESDtrack *trPos=esdev->GetTrack(v->GetIndex(1));
      if(!trNeg || !trPos){
	printf("Error: Couldn't get V0 daughter: %p - %p\n",trNeg,trPos);
	continue;
      }

      // reject tracks with neative ID
      if(trNeg->GetID()<0 || trPos->GetID()) continue;

      // at least one of the daughter has to pass basic QA cuts
      if(!(dauQAcuts1.IsSelected(trNeg) && dauQAcuts2.IsSelected(trNeg)) ||
	 !(dauQAcuts1.IsSelected(trPos) && dauQAcuts2.IsSelected(trPos))  ) continue;

      if(fMotherPdg==22) candidate.SetGammaTracks(trNeg, 11, trPos, 11);
      else candidate.SetTracks(trNeg, fNegPdg, trPos, fPosPdg);
      // eventually take the external trackparam and build the KFparticles by hand (see AliESDv0::GetKFInfo)
      // the folowing is not needed, because the daughters where used in the v0 vertex fit (I guess)
      //      AliKFVertex v0vtx = *v;
      //      candidate.SetProductionVertex(v0vtx);

      if(AliDielectronVarCuts::IsSelected(&candidate)) {
	nV0s++;
	//printf(" gamma found for vtx %p dau1id %d dau2id %d \n",v,trNeg->GetID(),trPos->GetID());
	fV0TrackArr.SetBitNumber(trNeg->GetID());
	fV0TrackArr.SetBitNumber(trPos->GetID());
      }
    }

  }
  else if(ev->IsA() == AliAODEvent::Class()) {
    const AliAODEvent *aodEv = static_cast<const AliAODEvent*>(ev);

    // loop over vertices
    for (Int_t ivertex=0; ivertex<aodEv->GetNumberOfVertices(); ++ivertex){
      AliAODVertex *v=aodEv->GetVertex(ivertex);
      if(v->GetType()!=AliAODVertex::kV0) continue;
      if(v->GetNDaughters()!=2) continue;

      AliAODTrack *trNeg=dynamic_cast<AliAODTrack*>(v->GetDaughter(0));
      AliAODTrack *trPos=dynamic_cast<AliAODTrack*>(v->GetDaughter(1));
      if(!trNeg || !trPos){
	printf("Error: Couldn't get V0 daughter: %p - %p\n",trNeg,trPos);
	continue;
      }

      // reject tracks with neative ID
      if(trNeg->GetID()<0 || trPos->GetID()) continue;

      // at least one of the daughter has to pass basic QA cuts
      if(!(dauQAcuts1.IsSelected(trNeg) && dauQAcuts2.IsSelected(trNeg)) ||
	 !(dauQAcuts1.IsSelected(trPos) && dauQAcuts2.IsSelected(trPos))  ) continue;

      AliKFVertex v0vtx = *v;
      if(fMotherPdg==22) candidate.SetGammaTracks(trNeg, 11, trPos, 11);
      else candidate.SetTracks(trNeg, (trNeg->Charge()<0?fNegPdg:fPosPdg), trPos, (trPos->Charge()<0?fNegPdg:fPosPdg));
      candidate.SetProductionVertex(v0vtx);

      if(AliDielectronVarCuts::IsSelected(&candidate)) {
	nV0s++;
	//printf(" gamma found for vtx %p dau1id %d dau2id %d \n",v,trNeg->GetID(),trPos->GetID());
	fV0TrackArr.SetBitNumber(trNeg->GetID());
	fV0TrackArr.SetBitNumber(trPos->GetID());
      }
    }

  }
  else
    return;

  //  printf(" Number of V0s candiates found %d \n",nV0s);

}
//________________________________________________________________________
Bool_t AliDielectronV0Cuts::IsSelected(TObject* track)
{
  //
  // Make cut decision
  //

  if(!track) return kFALSE;

  AliVTrack *vtrack = static_cast<AliVTrack*>(track);
  InitEvent(vtrack);
  //printf(" track ID %d selected result %d %d \n",vtrack->GetID(),(fV0TrackArr.TestBitNumber(vtrack->GetID())),fExcludeTracks);
  return ( (fV0TrackArr.TestBitNumber(vtrack->GetID()))^fExcludeTracks );
}

//________________________________________________________________________
Bool_t AliDielectronV0Cuts::IsNewEvent(const AliVEvent *ev)
{
  //
  // check weather we process a new event
  //

  //  printf(" current ev %d %d %d \n",fBunchCross, fOrbit, fPeriod);
  //  printf(" new event %p %d %d %d \n",ev, ev->GetBunchCrossNumber(), ev->GetOrbitNumber(), ev->GetPeriodNumber());

  if( fBunchCross == ev->GetBunchCrossNumber() ) {
    if( fOrbit == ev->GetOrbitNumber() )         {
      if( fPeriod == ev->GetPeriodNumber() )     {
	return kFALSE;
      }
    }
  }

  fBunchCross = ev->GetBunchCrossNumber();
  fOrbit      = ev->GetOrbitNumber();
  fPeriod     = ev->GetPeriodNumber();
  return kTRUE;
}
