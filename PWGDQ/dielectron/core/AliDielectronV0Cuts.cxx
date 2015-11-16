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
//   Cut class providing cuts to V0 candidates                           //
//   Selection or deselection of V0 candiates can be done.               //
//                                                                       //
// Authors:                                                              //
//   Julian Book <Julian.Book@cern.ch>                                   //
/*

Class to provide some V0 selection using exactley the same cuts in ESDs as in AODs.
It implements the PID cut class AliDielectronPID and the standard AliDielectronVarCuts for
the configuration of leg respective pair cuts. These pair cuts are applied on the KFparticle
build by the legs.

Some QA cuts for the tracks are applied before the V0 pair is build. These cuts are:
AliDielectronVarCuts dauQAcuts1;
dauQAcuts1.AddCut(AliDielectronVarManager::kPt,            0.05, 100.0);
dauQAcuts1.AddCut(AliDielectronVarManager::kEta,          -0.9,    0.9);
dauQAcuts1.AddCut(AliDielectronVarManager::kNclsTPC,      50.0,  160.0);
dauQAcuts1.AddCut(AliDielectronVarManager::kTPCchi2Cl,     0.0,    4.0);
AliDielectronTrackCuts dauQAcuts2;
dauQAcuts2.SetRequireTPCRefit(kTRUE);



Example configuration:

  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");

  // which V0 finder you want to use
  gammaV0Cuts->SetV0finder(kOnTheFly);  // kAll(default), kOffline or kOnTheFly

  // add some pdg codes (they are used then by the KF package and important for gamma conversions)
  gammaV0Cuts->SetPdgCodes(22,11,11); // mother, daughter1 and 2

  // add default PID cuts (defined in AliDielectronPID)
  // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
  gammaV0Cuts->SetDefaultPID(16, AliDielectronV0Cuts::kAny);

  // add the pair cuts for V0 candidates
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);

  // selection or rejection of V0 tracks
  gammaV0Cuts->SetExcludeTracks(kTRUE);

  // add the V0cuts directly to the track filter or to some cut group of it

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "AliDielectronV0Cuts.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronPID.h"
#include "AliESDv0.h"
#include "AliAODv0.h"

ClassImp(AliDielectronV0Cuts)


AliDielectronV0Cuts::AliDielectronV0Cuts() :
  AliDielectronVarCuts(),
  fV0TrackArr(0),
  fExcludeTracks(kTRUE),
  fV0finder(kAll),
  fMotherPdg(0),
  fNegPdg(0),
  fPosPdg(0),
  fPID(-1),
  fPIDCutType(kBoth),
  fOrbit(0),
  fPeriod(0),
  fBunchCross(0),
  fEventId(-1)
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
  fV0finder(kAll),
  fMotherPdg(0),
  fNegPdg(0),
  fPosPdg(0),
  fPID(-1),
  fPIDCutType(kBoth),
  fOrbit(0),
  fPeriod(0),
  fBunchCross(0),
  fEventId(-1)
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

  // rest booleans
  fV0TrackArr.ResetAllBits();

  // basic quality cut, /*at least one*/ both of the V0 daughters has to fullfill
  // always update ::Print accordingly
  AliDielectronVarCuts dauQAcuts1;
  dauQAcuts1.AddCut(AliDielectronVarManager::kNclsTPC,      70.0,  160.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kTPCchi2Cl,     0.0,    4.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kKinkIndex0,            0.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kEta,          -0.9,    0.9);
  dauQAcuts1.AddCut(AliDielectronVarManager::kPt,            0.05, 100.0);
  AliDielectronTrackCuts dauQAcuts2;
  //  dauQAcuts2.SetRequireITSRefit(kTRUE);
  dauQAcuts2.SetRequireTPCRefit(kTRUE);
  AliDielectronPID dauPIDcuts;
  if(fPID>=0) dauPIDcuts.SetDefaults(fPID);

  Int_t nV0s      = 0;
  Int_t nV0stored = 0;
  AliDielectronPair candidate;
  candidate.SetPdgCode(fMotherPdg);

  // ESD or AOD event
  if(ev->IsA() == AliESDEvent::Class()) {
    const AliESDEvent *esdev = static_cast<const AliESDEvent*>(ev);

    //printf("there are %d V0s in the event \n",esdev->GetNumberOfV0s());
    // loop over V0s
    for (Int_t iv=0; iv<esdev->GetNumberOfV0s(); ++iv){
      AliESDv0 *v = esdev->GetV0(iv);
      if(!v) continue;

      // check the v0 finder
      if( v->GetOnFlyStatus() && fV0finder==AliDielectronV0Cuts::kOffline  ) continue;
      if(!v->GetOnFlyStatus() && fV0finder==AliDielectronV0Cuts::kOnTheFly ) continue;

      // should we make use of AliESDv0Cuts::GetPdgCode() to preselect candiadtes, e.g.:
      // if(fMotherPdg!=v->GetPdgCode()) continue;

      AliESDtrack *trNeg=esdev->GetTrack(v->GetIndex(0));
      AliESDtrack *trPos=esdev->GetTrack(v->GetIndex(1));
      if(!trNeg || !trPos){
	printf("Error: Couldn't get V0 daughter: %p - %p\n",trNeg,trPos);
	continue;
      }

      // protection against LS v0s
      if(trNeg->Charge() == trPos->Charge()) continue;

      // PID default cuts
      if(fPID>=0) {
	Bool_t selected=kFALSE;
	selected=dauPIDcuts.IsSelected(trNeg);
	if(fPIDCutType==kBoth) selected &= dauPIDcuts.IsSelected(trPos);
	if(fPIDCutType==kAny)  selected |= dauPIDcuts.IsSelected(trPos);
	if(!selected) continue;
      }

      // basic track cuts
      if( !dauQAcuts2.IsSelected(trNeg) ) continue;
      if( !dauQAcuts2.IsSelected(trPos) ) continue;
      if( !dauQAcuts1.IsSelected(trNeg) ) continue;
      if( !dauQAcuts1.IsSelected(trPos) ) continue;

      if(fMotherPdg==22) candidate.SetGammaTracks(trNeg, 11, trPos, 11);
      else candidate.SetTracks(trNeg, (trNeg->Charge()<0?fNegPdg:fPosPdg), trPos, (trPos->Charge()<0?fNegPdg:fPosPdg));
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
    if(!aodEv->GetV0s()) return; // protection for nano AODs

    // loop over vertices
    for (Int_t ivertex=0; ivertex<aodEv->GetNumberOfV0s(); ++ivertex){
      AliAODv0 *v=aodEv->GetV0(ivertex);
      if(!v) continue;

      // check the v0 finder
      if( v->GetOnFlyStatus() && fV0finder==AliDielectronV0Cuts::kOffline  ) continue;
      if(!v->GetOnFlyStatus() && fV0finder==AliDielectronV0Cuts::kOnTheFly ) continue;

      AliAODTrack *trNeg=dynamic_cast<AliAODTrack*>(v->GetDaughter(0));
      AliAODTrack *trPos=dynamic_cast<AliAODTrack*>(v->GetDaughter(1));
      if(!trNeg || !trPos){
	printf("Error: Couldn't get V0 daughter: %p - %p\n",trNeg,trPos);
	continue;
      }
      nV0stored++;

      // protection against LS v0s
      if(trNeg->Charge() == trPos->Charge()) continue;

      // PID default cuts
      if(fPID>=0) {
	Bool_t selected=kFALSE;
	selected=dauPIDcuts.IsSelected(trNeg);
	if(fPIDCutType==kBoth) selected &= dauPIDcuts.IsSelected(trPos);
	if(fPIDCutType==kAny)  selected |= dauPIDcuts.IsSelected(trPos);
	if(!selected) continue;
      }

      // basic track cuts
      if( !dauQAcuts2.IsSelected(trNeg) ) continue;
      if( !dauQAcuts2.IsSelected(trPos) ) continue;
      if( !dauQAcuts1.IsSelected(trNeg) ) continue;
      if( !dauQAcuts1.IsSelected(trPos) ) continue;

      AliKFVertex v0vtx = *(v->GetSecondaryVtx());
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
    //    printf("there are %d V0s in the event \n",nV0stored);
  }
  else
    return;

  //  printf(" Number of V0s candiates found %d/%d \n",nV0s,nV0stored);

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
 
  // NOTE: if event number in file is not enough then use in addition the inputfilename
  //  (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->GetInputFileName()

  if( fEventId == ev->GetEventNumberInFile() )     {
    if( fBunchCross == ev->GetBunchCrossNumber() ) {
      if( fOrbit == ev->GetOrbitNumber() )         {
	if( fPeriod == ev->GetPeriodNumber() )     {
	return kFALSE;
	}
      }
    }
  }

  fBunchCross = ev->GetBunchCrossNumber();
  fOrbit      = ev->GetOrbitNumber();
  fPeriod     = ev->GetPeriodNumber();
  fEventId    = ev->GetEventNumberInFile();

  return kTRUE;
}

//________________________________________________________________________
void AliDielectronV0Cuts::Print(const Option_t* /*option*/) const
{
  //
  // Print cuts and the range
  //
  printf(" V0 cuts:\n");
  printf(" V0 finder mode: %s \n",(fV0finder ==kOnTheFly ? "One-The-Fly":
				   (fV0finder==kOffline  ? "Offline":
				    "One-The-Fly+Offline") ) );
  AliDielectronVarCuts::Print();

  printf(" V0 daughter cuts (applied to both):\n");
  AliDielectronVarCuts dauQAcuts1;
  dauQAcuts1.AddCut(AliDielectronVarManager::kNclsTPC,      70.0,  160.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kTPCchi2Cl,     0.0,    4.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kKinkIndex0,            0.0);
  dauQAcuts1.AddCut(AliDielectronVarManager::kEta,          -0.9,    0.9);
  dauQAcuts1.AddCut(AliDielectronVarManager::kPt,            0.05, 100.0);
  dauQAcuts1.Print();
  AliDielectronTrackCuts dauQAcuts2;
  //  dauQAcuts2.SetRequireITSRefit(kTRUE);
  dauQAcuts2.SetRequireTPCRefit(kTRUE);
  //  dauQAcuts2.Print(); //TODO activate as soon as implemented

  if(fPID>=0) {
    printf(" V0 daughter PID cuts (applied to %s):\n",(fPIDCutType==kBoth?"both":"any"));
    AliDielectronPID dauPIDcuts;
    dauPIDcuts.SetDefaults(fPID);
    dauPIDcuts.Print(); //TODO activate as soon as implemented
  }

}
