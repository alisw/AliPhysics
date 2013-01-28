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
  fV0TrackArr(0x0)
{
  //
  // Default costructor
  //
}

//________________________________________________________________________
AliDielectronV0Cuts::AliDielectronV0Cuts(const char* name, const char* title) :
  AliDielectronVarCuts(name,title),
  fV0TrackArr(0x0)
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
  if(fV0TrackArr) {
    delete fV0TrackArr;
    fV0TrackArr=0;
  }
}


//________________________________________________________________________
void AliDielectronV0Cuts::Init()
{
  //
  // Init the V0 candidates
  //

  // TODO think about MCevent
  //  Print();

  // take current event from the varmanager
  AliVEvent *ev  =   AliDielectronVarManager::GetCurrentEvent();
  if(!ev) return;
  fV0TrackArr = new TArrayC(ev->GetNumberOfTracks());

  // basic quality cut, at least one of the V0 daughters has to fullfill
  AliDielectronVarCuts *dauQAcuts1 = new AliDielectronVarCuts();
  dauQAcuts1->AddCut(AliDielectronVarManager::kPt,           0.5,  1e30);
  dauQAcuts1->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  dauQAcuts1->AddCut(AliDielectronVarManager::kNclsTPC,     50.0, 160.0);
  AliDielectronTrackCuts *dauQAcuts2 = new AliDielectronTrackCuts();
  dauQAcuts2->SetRequireITSRefit(kTRUE);
  dauQAcuts2->SetRequireTPCRefit(kTRUE);

  Int_t nV0s = 0;

  // ESD or AOD event
  if(ev->IsA() == AliESDEvent::Class()) {
    const AliESDEvent *esdev = static_cast<const AliESDEvent*>(ev);
    // loop over V0s
    for (Int_t iv=0; iv<esdev->GetNumberOfV0s(); ++iv){
      AliESDv0 *v = esdev->GetV0(iv);
      if(!v) continue;

      AliESDtrack *tr1=esdev->GetTrack(v->GetIndex(0));
      AliESDtrack *tr2=esdev->GetTrack(v->GetIndex(1));
      if(!tr1 || !tr2){
	printf("Error: Couldn't get V0 daughter: %p - %p\n",tr1,tr2);
	continue;
      }

      // at least one of the daughter has to pass basic QA cuts
      if(!(dauQAcuts1->IsSelected(tr1) && dauQAcuts2->IsSelected(tr1)) ||
	 !(dauQAcuts1->IsSelected(tr1) && dauQAcuts2->IsSelected(tr1))  ) continue;
      //      printf(" One or both V0 daughters pass the qa cuts \n");

      // eventually take the external trackparam and build the KFparticles by hand (see AliESDv0::GetKFInfo)
      AliDielectronPair candidate;
      candidate.SetPdgCode(22);
      candidate.SetGammaTracks(tr1, 11, tr2, 11);
      // this is not needed, because the daughters where used in the v0 vertex fit (I guess)
      //      AliKFVertex v0vtx = *v;
      //      candidate.SetProductionVertex(v0vtx);

      //Fill values
      Double_t values[AliDielectronVarManager::kNMaxValues];
      AliDielectronVarManager::Fill(&candidate,values);

      fSelectedCutsMask=0;
      for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
	Int_t cut=fActiveCuts[iCut];
	SETBIT(fSelectedCutsMask,iCut);
	if ( ((values[cut]<fCutMin[iCut]) || (values[cut]>fCutMax[iCut]))^fCutExclude[iCut] ) CLRBIT(fSelectedCutsMask,iCut);
      }

      Bool_t isSelected=(fSelectedCutsMask==fActiveCutsMask);
      if ( fCutType==kAny ) isSelected=(fSelectedCutsMask>0);

      // store boolean at index=trackID
      if(isSelected) {
	nV0s++;
	if(tr1->GetID()>fV0TrackArr->GetSize()) {/* printf(" size of array %d too small expand to %d \n",fV0TrackArr->GetSize(), tr1->GetID());*/ fV0TrackArr->Set(tr1->GetID()+1); }
	if(tr2->GetID()>fV0TrackArr->GetSize()) {/* printf(" size of array %d too small expand to %d \n",fV0TrackArr->GetSize(), tr2->GetID());*/ fV0TrackArr->Set(tr2->GetID()+1); }
	//	printf(" gamma found for vtx %p dau1id %d dau2id %d \n",v,tr1->GetID(),tr2->GetID());
	fV0TrackArr->AddAt(1,tr1->GetID());
	fV0TrackArr->AddAt(1,tr2->GetID());
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

      AliAODTrack *tr1=dynamic_cast<AliAODTrack*>(v->GetDaughter(0));
      AliAODTrack *tr2=dynamic_cast<AliAODTrack*>(v->GetDaughter(1));
      if(!tr1 || !tr2){
	printf("Error: Couldn't get V0 daughter: %p - %p\n",tr1,tr2);
	continue;
      }

      // at least one of the daughter has to pass basic QA cuts
      if(!(dauQAcuts1->IsSelected(tr1) && dauQAcuts2->IsSelected(tr1)) ||
	 !(dauQAcuts1->IsSelected(tr1) && dauQAcuts2->IsSelected(tr1))  ) continue;

      //      printf(" One or both V0 daughters pass the qa cuts \n");

      AliKFVertex v0vtx = *v;
      AliDielectronPair candidate;
      candidate.SetPdgCode(22); // TODO setter?
      candidate.SetGammaTracks(tr1, 11, tr2, 11);
      candidate.SetProductionVertex(v0vtx);

      //Fill values
      Double_t values[AliDielectronVarManager::kNMaxValues];
      AliDielectronVarManager::Fill(&candidate,values);

      fSelectedCutsMask=0;
      for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
	Int_t cut=fActiveCuts[iCut];
	SETBIT(fSelectedCutsMask,iCut);
	if ( ((values[cut]<fCutMin[iCut]) || (values[cut]>fCutMax[iCut]))^fCutExclude[iCut] ) CLRBIT(fSelectedCutsMask,iCut);
      }

      Bool_t isSelected=(fSelectedCutsMask==fActiveCutsMask);
      if ( fCutType==kAny ) isSelected=(fSelectedCutsMask>0);

      // store boolean at index=trackID
      if(isSelected) {
	nV0s++;
	if(tr1->GetID()>fV0TrackArr->GetSize()) {/* printf(" size of array %d too small expand to %d \n",fV0TrackArr->GetSize(), tr1->GetID());*/ fV0TrackArr->Set(tr1->GetID()+1); }
	if(tr2->GetID()>fV0TrackArr->GetSize()) {/* printf(" size of array %d too small expand to %d \n",fV0TrackArr->GetSize(), tr2->GetID());*/ fV0TrackArr->Set(tr2->GetID()+1); }
	//	printf(" gamma found for vtx %p dau1id %d dau2id %d \n",v,tr1->GetID(),tr2->GetID());
	fV0TrackArr->AddAt(1,tr1->GetID());
	fV0TrackArr->AddAt(1,tr2->GetID());
      }
    }
  }
  else
    return;

  delete dauQAcuts1;
  delete dauQAcuts2;
  //  printf(" Number of V0s candiates found %d \n",nV0s);

}
//________________________________________________________________________
Bool_t AliDielectronV0Cuts::IsSelected(TObject* track)
{
  //
  // Make cut decision
  //
  if(!fV0TrackArr->GetSum()) return kTRUE;
  if(!track) return kFALSE;
  //what about VParticles MC, better to store pointers??
  AliVTrack *vtrack = static_cast<AliVTrack*>(track);
  if(!vtrack) return kFALSE;

  //  printf(" tttttrackID %d \n",vtrack->GetID());

  if(fV0TrackArr->At(vtrack->GetID()) != 1) return kTRUE;
  else {
    //printf(" track belongs to a V0 candidate \n");
    return kFALSE;
  }

}
