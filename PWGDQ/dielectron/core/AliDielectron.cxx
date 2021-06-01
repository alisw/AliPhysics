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
//                Dielectron Analysis Main class                         //
//                                                                       //
/*
Framework to perform event selectoin, single track selection and track pair
selection.

Convention for the signs of the pair in fPairCandidates:
The names are available via the function PairClassName(Int_t i)

0: ev1+ ev1+  (same event like sign +)
1: ev1+ ev1-  (same event unlike sign)
2: ev1- ev1-  (same event like sign -)

3: ev1+ ev2+  (mixed event like sign +)
4: ev1- ev2+  (mixed event unlike sign -+)
6: ev1+ ev2-  (mixed event unlike sign +-)
7: ev1- ev2-  (mixed event like sign -)

5: ev2+ ev2+  (same event like sign +)
8: ev2+ ev2-  (same event unlike sign)
9: ev2- ev2-  (same event like sign -)

10: ev1+ ev1- (same event track rotation)
11: ev1+ ev1+ (same event track rotation)
12: ev1- ev1- (same event track rotation)

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TList.h>
#include <TMath.h>
#include <TObject.h>
#include <TGrid.h>
#include <TSystem.h>

#include <AliKFParticle.h>

#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliEPSelectionTask.h>
#include <AliEventplane.h>
#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliVTrack.h>
#include <AliLog.h>
#include "AliDielectronPair.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTrackRotator.h"
#include "AliDielectronDebugTree.h"
#include "AliDielectronSignalMC.h"
#include "AliDielectronMixingHandler.h"
#include "AliDielectronPairLegCuts.h"
#include "AliDielectronV0Cuts.h"
#include "AliDielectronPID.h"
#include "AliDielectronHistos.h"

#include "AliDielectron.h"

ClassImp(AliDielectron)

const char* AliDielectron::fgkTrackClassNames[6] = {
  "ev1+",
  "ev1-",
  "ev2+",
  "ev2-",
  "ev1_TR+",
  "ev1_TR-"
};

const char* AliDielectron::fgkPairClassNames[13] = {
  "ev1+_ev1+",
  "ev1+_ev1-",
  "ev1-_ev1-",
  "ev1+_ev2+",
  "ev1-_ev2+",
  "ev2+_ev2+",
  "ev1+_ev2-",
  "ev1-_ev2-",
  "ev2+_ev2-",
  "ev2-_ev2-",
  "ev1+_ev1-_TR",
  "ev1+_ev1+_TR",
  "ev1-_ev1-_TR"
};

//________________________________________________________________
AliDielectron::AliDielectron() :
  TNamed("AliDielectron","AliDielectron"),
  fCutQA(kFALSE),
  fQAmonitor(0x0),
  fPostPIDCntrdCorrArr(0x0),
  fPostPIDWdthCorrArr(0x0),
  fPostPIDCntrdCorr(0x0),
  fPostPIDWdthCorr(0x0),
  fPostPIDCntrdCorrITS(0x0),
  fPostPIDWdthCorrITS(0x0),
  fPostPIDCntrdCorrTOF(0x0),
  fPostPIDWdthCorrTOF(0x0),
  fRotateTrackCorrectionMap(),
	fPIDCalibinPU(kFALSE),
  fLegEffMap(0x0),
  fPairEffMap(0x0),
  fEventFilter("EventFilter"),
  fTrackFilter("TrackFilter"),
  fPairPreFilter1("PairPreFilter1"),
  fPairPreFilter2("PairPreFilter2"),
  fPairPreFilterLegs1("PairPreFilterLegs1"),
  fPairPreFilterLegs2("PairPreFilterLegs2"),
  fPairFilter("PairFilter"),
  fEventPlanePreFilter("EventPlanePreFilter"),
  fEventPlanePOIPreFilter("EventPlanePOIPreFilter"),
  fQnTPCACcuts(0x0),
  fQnVectorNorm(""),
  fPdgMother(443),
  fPdgLeg1(11),
  fPdgLeg2(11),
  fSignalsMC(0x0),
  fNoPairing(kFALSE),
  fProcessLS(kTRUE),
  fUseKF(kTRUE),
  fHistoArray(0x0),
  fHistos(0x0),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fPairCandidates(new TObjArray(13)),
  fCfManagerPair(0x0),
  fTrackRotator(0x0),
  fRotatePP(kFALSE),
  fRotateMM(kFALSE),
  fDebugTree(0x0),
  fMixing(0x0),
  fEvtVsTrkHist(0x0),
  fEvtVsTrkHistExists(kFALSE),
  fPreFilterEventPlane(kFALSE),
  fACremovalIsSetted(kFALSE),
  fLikeSignSubEvents(kFALSE),
  fPreFilterUnlikeOnly1(kFALSE),
  fPreFilterLikeOnly1(kFALSE),
  fPreFilterAllSigns1(kFALSE),
  fPreFilterPhotons1(kFALSE),
  fPreFilterOnlyOnePair1(kFALSE),
  fPreFilterUnlikeOnly2(kFALSE),
  fPreFilterLikeOnly2(kFALSE),
  fPreFilterAllSigns2(kFALSE),
  fPreFilterPhotons2(kFALSE),
  fPreFilterOnlyOnePair2(kFALSE),
  fHasMC(kFALSE),
  fStoreRotatedPairs(kFALSE),
  fDontClearArrays(kFALSE),
  fEventProcess(kTRUE),
  fUseGammaTracks(kTRUE),
  fEstimatorFilename(""),
  fEstimatorObjArray(0x0),
  fTRDpidCorrectionFilename(""),
  fQnCalibrationFilepath(""),
  fDoQnV0GainEqualization(kFALSE),
  fDoQnV0Recentering(kFALSE),
  fDoQnTPCRecentering(kFALSE),
  fVZEROCalibrationFilename(""),
  fVZERORecenteringFilename(""),
  fZDCRecenteringFilename(""),
  fUseAccMap(kTRUE),
  fIterations(1)


{
  //
  // Default constructor
  //

	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<15;j++){
			fPostPIDCntrdCorrPU[i][j] = 0x0;
			fPostPIDWdthCorrPU[i][j]  = 0x0;

		}
	}


}

//________________________________________________________________
AliDielectron::AliDielectron(const char* name, const char* title) :
  TNamed(name,title),
  fCutQA(kFALSE),
  fQAmonitor(0x0),
  fPostPIDCntrdCorrArr(0x0),
  fPostPIDWdthCorrArr(0x0),
  fPostPIDCntrdCorr(0x0),
  fPostPIDWdthCorr(0x0),
  fPostPIDCntrdCorrITS(0x0),
  fPostPIDWdthCorrITS(0x0),
  fPostPIDCntrdCorrTOF(0x0),
  fPostPIDWdthCorrTOF(0x0),
	fPIDCalibinPU(kFALSE),
  fLegEffMap(0x0),
  fPairEffMap(0x0),
  fRotateTrackCorrectionMap(),
  fEventFilter("EventFilter"),
  fTrackFilter("TrackFilter"),
  fPairPreFilter1("PairPreFilter1"),
  fPairPreFilter2("PairPreFilter2"),
  fPairPreFilterLegs1("PairPreFilterLegs1"),
  fPairPreFilterLegs2("PairPreFilterLegs2"),
  fPairFilter("PairFilter"),
  fEventPlanePreFilter("EventPlanePreFilter"),
  fEventPlanePOIPreFilter("EventPlanePOIPreFilter"),
  fQnTPCACcuts(0x0),
  fQnVectorNorm(""),
  fPdgMother(443),
  fPdgLeg1(11),
  fPdgLeg2(11),
  fSignalsMC(0x0),
  fNoPairing(kFALSE),
  fProcessLS(kTRUE),
  fUseKF(kTRUE),
  fHistoArray(0x0),
  fHistos(0x0),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fPairCandidates(new TObjArray(13)),
  fCfManagerPair(0x0),
  fTrackRotator(0x0),
  fRotatePP(kFALSE),
  fRotateMM(kFALSE),
  fDebugTree(0x0),
  fMixing(0x0),
  fEvtVsTrkHist(0x0),
  fEvtVsTrkHistExists(kFALSE),
  fPreFilterEventPlane(kFALSE),
  fACremovalIsSetted(kFALSE),
  fLikeSignSubEvents(kFALSE),
  fPreFilterUnlikeOnly1(kFALSE),
  fPreFilterLikeOnly1(kFALSE),
  fPreFilterAllSigns1(kFALSE),
  fPreFilterPhotons1(kFALSE),
  fPreFilterOnlyOnePair1(kFALSE),
  fPreFilterUnlikeOnly2(kFALSE),
  fPreFilterLikeOnly2(kFALSE),
  fPreFilterAllSigns2(kFALSE),
  fPreFilterPhotons2(kFALSE),
  fPreFilterOnlyOnePair2(kFALSE),
  fHasMC(kFALSE),
  fStoreRotatedPairs(kFALSE),
  fDontClearArrays(kFALSE),
  fEventProcess(kTRUE),
  fUseGammaTracks(kTRUE),
  fEstimatorFilename(""),
  fEstimatorObjArray(0x0),
  fTRDpidCorrectionFilename(""),
  fQnCalibrationFilepath(""),
  fDoQnV0GainEqualization(kFALSE),
  fDoQnV0Recentering(kFALSE),
  fDoQnTPCRecentering(kFALSE),
  fVZEROCalibrationFilename(""),
  fVZERORecenteringFilename(""),
  fZDCRecenteringFilename(""),
  fUseAccMap(kTRUE),
  fIterations(1)
{
  //
  // Named constructor
  //

	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<15;j++){
			fPostPIDCntrdCorrPU[i][j] = 0x0;
			fPostPIDWdthCorrPU[i][j]  = 0x0;
		}
	}

}

//________________________________________________________________
AliDielectron::~AliDielectron()
{
  //
  // Default destructor
  //
  if (fQAmonitor) delete fQAmonitor;
  if (fPostPIDCntrdCorrArr) delete fPostPIDCntrdCorrArr;
  if (fPostPIDWdthCorrArr) delete fPostPIDWdthCorrArr;
  if (fPostPIDCntrdCorr) delete fPostPIDCntrdCorr;
  if (fPostPIDWdthCorr) delete fPostPIDWdthCorr;
  if (fPostPIDCntrdCorrITS) delete fPostPIDCntrdCorrITS;
  if (fPostPIDWdthCorrITS) delete fPostPIDWdthCorrITS;
  if (fPostPIDCntrdCorrTOF) delete fPostPIDCntrdCorrTOF;
  if (fPostPIDWdthCorrTOF) delete fPostPIDWdthCorrTOF;
  if (fLegEffMap) delete fLegEffMap;
  if (fPairEffMap) delete fPairEffMap;
  if (fHistos) delete fHistos;
  if (fUsedVars) delete fUsedVars;
  if (fPairCandidates && fEventProcess) delete fPairCandidates;
  if (fDebugTree) delete fDebugTree;
  if (fMixing) delete fMixing;
  if (fEvtVsTrkHist) delete fEvtVsTrkHist;
  if (fSignalsMC) delete fSignalsMC;
  if (fCfManagerPair) delete fCfManagerPair;
  if (fHistoArray) delete fHistoArray;



	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<15;j++){
			if(fPostPIDCntrdCorrPU[i][j]) delete fPostPIDCntrdCorrPU[i][j];
			if(fPostPIDWdthCorrPU[i][j] ) delete fPostPIDWdthCorrPU[i][j] ;
		}
	}
}
//________________________________________________________________
void AliDielectron::Init()
{
  //
  // Initialise objects
  //
  if(GetHasMC()) AliDielectronMC::Instance()->SetHasMC(GetHasMC());

  if(fEventProcess) InitPairCandidateArrays();

  if (fCfManagerPair) {
    fCfManagerPair->SetSignalsMC(fSignalsMC);
    fCfManagerPair->InitialiseContainer(fPairFilter);
  }
  if (fTrackRotator)  {
    if(fRotatePP){
      fTrackRotator->SetTrackArrays(&fTracks[0],&fTracks[0]);
    }
    else if(fRotateMM){
      fTrackRotator->SetTrackArrays(&fTracks[1],&fTracks[1]);
    }
    else{
      fTrackRotator->SetTrackArrays(&fTracks[0],&fTracks[1]);
    }
    fTrackRotator->SetPdgLegs(fPdgLeg1,fPdgLeg2);
  }
  if (fDebugTree) fDebugTree->SetDielectron(this);

  if(fEstimatorFilename.Contains(".root"))        AliDielectronVarManager::InitEstimatorAvg(fEstimatorFilename.Data());
  if(fEstimatorObjArray)			  AliDielectronVarManager::InitEstimatorObjArrayAvg(fEstimatorObjArray);
  if(fTRDpidCorrectionFilename.Contains(".root")) AliDielectronVarManager::InitTRDpidEffHistograms(fTRDpidCorrectionFilename.Data());
  if(fVZEROCalibrationFilename.Contains(".root")) AliDielectronVarManager::SetVZEROCalibrationFile(fVZEROCalibrationFilename.Data());
  if(fVZERORecenteringFilename.Contains(".root")) AliDielectronVarManager::SetVZERORecenteringFile(fVZERORecenteringFilename.Data());
  if(fZDCRecenteringFilename.Contains(".root"))   AliDielectronVarManager::SetZDCRecenteringFile(fZDCRecenteringFilename.Data());

  if(fQnCalibrationFilepath != "") AliDielectronVarManager::SetQnCalibrationFilePath(fQnCalibrationFilepath.Data(), fDoQnV0GainEqualization, fDoQnV0Recentering, fDoQnTPCRecentering);

  if (fMixing) fMixing->Init(this);
  if (fHistoArray) {
    fHistoArray->SetSignalsMC(fSignalsMC);
    fHistoArray->Init();
  }

  if(!fEventProcess) {
    AliDielectronPairLegCuts *trk2leg = new AliDielectronPairLegCuts("trk2leg","trk2leg");
    // move all track cuts (if any) into pair leg cuts
    TIter listIterator(fTrackFilter.GetCuts());
    while (AliAnalysisCuts *thisCut = (AliAnalysisCuts*) listIterator()) {
      trk2leg->GetLeg1Filter().AddCuts((AliAnalysisCuts*)thisCut->Clone());
      trk2leg->GetLeg2Filter().AddCuts((AliAnalysisCuts*)thisCut->Clone());
    }
    // add pair leg cuts to pair filter
    fPairFilter.AddCuts(trk2leg);
  }

  if (fCutQA) {
    fQAmonitor = new AliDielectronCutQA(Form("QAcuts_%s",GetName()),"QAcuts");
    fQAmonitor->AddTrackFilter(&fTrackFilter);
    if(!fNoPairing) fQAmonitor->AddPairFilter(&fPairFilter);
    fQAmonitor->AddEventFilter(&fEventFilter);
    fQAmonitor->Init();
  }

  if(fHistos) {
    (*fUsedVars)|= (*fHistos->GetUsedVars());

    // Initialisation of AliDielectronEvtVsTrkHist
    if(fHistos->GetHistogramList()->FindObject("EvtVsTrk")){
      fEvtVsTrkHist = new AliDielectronEvtVsTrkHist("EvtVsTrkHistos", "EvtVsTrkHistos");
      fEvtVsTrkHist->SetHistogramList(fHistos);
    }
  }
}

//________________________________________________________________

void AliDielectron::Process(TObjArray *arr)
{
  //
  // Process the pair array
  //

  // set pair arrays
  fPairCandidates = arr;

  //fill debug tree if a manager is attached
  //  if (fDebugTree) FillDebugTree();
  //in case there is a histogram manager, fill the QA histograms
  //  if (fHistos && fSignalsMC) FillMCHistograms(ev1);

  // apply cuts and fill output
  if (fHistos) FillHistogramsFromPairArray();

  // never clear arrays !!!!


}

//________________________________________________________________
Bool_t AliDielectron::Process(AliVEvent *ev1, AliVEvent *ev2)
{
  //
  // Process the events
  //
  //at least first event is needed!
  if (!ev1){
    AliError("At least first event must be set!");
    return 0;
  }

  // modify event numbers in MC so that we can identify new events
  // in AliDielectronV0Cuts (not neeeded for collision data)
  if(GetHasMC()) {
    ev1->SetBunchCrossNumber(1);
    ev1->SetOrbitNumber(1);
    ev1->SetPeriodNumber(1);
  }

  // set qn vector normalisation to var manager
  AliDielectronVarManager::SetQnVectorNormalisation(fQnVectorNorm);

  // set pid correction function to var manager
  // from an array of objects for run-wise differences
  if(fPostPIDCntrdCorrArr){
    TString key;
    key.Form("Centroid_%d",ev1->GetRunNumber());
    TH1* hCent = (TH1*) fPostPIDCntrdCorrArr->FindObject(key.Data());
    AliDielectronPID::SetCentroidCorrFunction((TH1*) hCent->Clone());
  }
  if(fPostPIDWdthCorrArr){
    TString key;
    key.Form("Width_%d",ev1->GetRunNumber());
    TH1* hWidth = (TH1*) fPostPIDWdthCorrArr->FindObject(key.Data());
    AliDielectronPID::SetWidthCorrFunction((TH1*) hWidth->Clone());
  }
  // from one object for the full period
  if(fPostPIDCntrdCorr)     AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorr);
  if(fPostPIDWdthCorr)      AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorr);
  if(fPostPIDCntrdCorrITS)  AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorrITS);
  if(fPostPIDWdthCorrITS)   AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDWdthCorrITS);
  if(fPostPIDCntrdCorrTOF)  AliDielectronPID::SetCentroidCorrFunctionTOF(fPostPIDCntrdCorrTOF);
  if(fPostPIDWdthCorrTOF)   AliDielectronPID::SetWidthCorrFunctionTOF(fPostPIDWdthCorrTOF);

	for(Int_t id=0;id<15;id++){//detector loop TPC/ITS/TOF
		for(Int_t ip=0;ip<15;ip++){//particle loop e/mu/pi/k/p
			if(fPostPIDCntrdCorrPU[id][ip])   AliDielectronPID::SetCentroidCorrFunctionPU(id,ip,fPostPIDCntrdCorrPU[id][ip]);
			if(fPostPIDWdthCorrPU[id][ip])    AliDielectronPID::SetWidthCorrFunctionPU(   id,ip,fPostPIDWdthCorrPU[id][ip] );
		}
	}

	AliDielectronPID::SetPIDCalibinPU(fPIDCalibinPU);

  // set event
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::SetEvent(ev1);

  if (fMixing){
    //set mixing bin to event data
    Int_t bin=fMixing->FindBin(AliDielectronVarManager::GetData());
    AliDielectronVarManager::SetValue(AliDielectronVarManager::kMixingBin,bin);
  }

  // set efficiency maps
  AliDielectronVarManager::SetLegEffMap(fLegEffMap);
  AliDielectronVarManager::SetPairEffMap(fPairEffMap);

  //in case we have MC load the MC event and process the MC particles
  // why do not apply the event cuts first ????
  if (AliDielectronMC::Instance()->ConnectMCEvent()){
    ProcessMC(ev1);
  }

  //if candidate array doesn't exist, create it
  if (!fPairCandidates->UncheckedAt(0)) {
    InitPairCandidateArrays();
  } else {
    ClearArrays();
  }

  //mask used to require that all cuts are fulfilled
  UInt_t selectedMask=(1<<fEventFilter.GetCuts()->GetEntries())-1;

  //apply event cuts
  UInt_t cutmask = fEventFilter.IsSelected(ev1);
  if(fCutQA) fQAmonitor->FillAll(ev1);
  if(fCutQA) fQAmonitor->Fill(cutmask,ev1);
  if ((ev1&&cutmask!=selectedMask) ||
      (ev2&&fEventFilter.IsSelected(ev2)!=selectedMask)) return 0;

  if(fEvtVsTrkHist){
    fEvtVsTrkHist->SetPIDResponse(AliDielectronVarManager::GetPIDResponse());
    fEvtVsTrkHist->FillHistograms(ev1);
  }

  //fill track arrays for the first event
  if (ev1){
    FillTrackArrays(ev1);
    if (((fPreFilterAllSigns1)||(fPreFilterUnlikeOnly1)) && !fPreFilterLikeOnly1 && ( fPairPreFilter1.GetCuts()->GetEntries()>0 )) PairPreFilter(0, 1, fTracks[0], fTracks[1], ev1, 1);
    if (((fPreFilterAllSigns2)||(fPreFilterUnlikeOnly2)) && !fPreFilterLikeOnly2 && ( fPairPreFilter2.GetCuts()->GetEntries()>0 )) PairPreFilter(0, 1, fTracks[0], fTracks[1], ev1, 2);

    if(fPreFilterLikeOnly1 && fPairPreFilter1.GetCuts()->GetEntries()>0 ){
       PairPreFilter(0, 0, fTracks[0], fTracks[0], ev1, 1);
       PairPreFilter(1, 1, fTracks[1], fTracks[1], ev1, 1);
    }

    if(fPreFilterLikeOnly2 && fPairPreFilter2.GetCuts()->GetEntries()>0 ){
       PairPreFilter(0, 0, fTracks[0], fTracks[0], ev1, 2);
       PairPreFilter(1, 1, fTracks[1], fTracks[1], ev1, 2);
    }
  }


  //fill track arrays for the second event
  if (ev2) {
    FillTrackArrays(ev2,1);
    if (((fPreFilterAllSigns1)||(fPreFilterUnlikeOnly1)) && !fPreFilterLikeOnly1 && ( fPairPreFilter1.GetCuts()->GetEntries()>0 )) PairPreFilter(2, 3, fTracks[2], fTracks[3], ev2, 1);
    if (((fPreFilterAllSigns2)||(fPreFilterUnlikeOnly2)) && !fPreFilterLikeOnly2 && ( fPairPreFilter2.GetCuts()->GetEntries()>0 )) PairPreFilter(2, 3, fTracks[2], fTracks[3], ev2, 2);

    if(fPreFilterLikeOnly1 && fPairPreFilter1.GetCuts()->GetEntries()>0 ){
       PairPreFilter(2, 2, fTracks[2], fTracks[2], ev2, 1);
       PairPreFilter(3, 3, fTracks[3], fTracks[3], ev2, 1);
    }

    if(fPreFilterLikeOnly2 && fPairPreFilter2.GetCuts()->GetEntries()>0 ){
       PairPreFilter(2, 2, fTracks[2], fTracks[2], ev2, 2);
       PairPreFilter(3, 3, fTracks[3], fTracks[3], ev2, 2);
    }
  }

  // TPC event plane correction
  if (ev1 && fPreFilterEventPlane && ( fEventPlanePreFilter.GetCuts()->GetEntries()>0 || fEventPlanePOIPreFilter.GetCuts()->GetEntries()>0))
    EventPlanePreFilter(0, 1, fTracks[0], fTracks[1], ev1);
  // QnFramework est. 2016 auto-correlation removal
  if(fACremovalIsSetted){
    AliDielectronVarManager::SetTPCEventPlaneACremoval(fQnTPCACcuts);
  }
  if (!fNoPairing){
    // create pairs and fill pair candidate arrays
    for (Int_t itrackArr1=0; itrackArr1<4; ++itrackArr1){
      for (Int_t itrackArr2=itrackArr1; itrackArr2<4; ++itrackArr2){
	if(!fProcessLS && GetPairIndex(itrackArr1,itrackArr2)!=kEv1PM) continue;
	FillPairArrays(itrackArr1, itrackArr2, ev1);
      }
    }

    //track rotation
    if (fTrackRotator) {
      fTrackRotator->SetEvent(ev1);
      FillPairArrayTR();
    }
  }

  //fill debug tree if a manager is attached
  if (fDebugTree) FillDebugTree();

  //process event mixing
  if (fMixing) {
    fMixing->Fill(ev1,this);
    //     FillHistograms(0x0,kTRUE);
  }

  // fill candidate variables
  Double_t ntracks = fTracks[0].GetEntriesFast() + fTracks[1].GetEntriesFast();
  Double_t npairs  = PairArray(AliDielectron::kEv1PM)->GetEntriesFast();
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kTracks, ntracks);
  AliDielectronVarManager::SetValue(AliDielectronVarManager::kPairs,  npairs);

  //in case there is a histogram manager, fill the QA histograms
  if (fHistos && fSignalsMC) FillMCHistograms(ev1);
  if (fHistos) FillHistograms(ev1);
  // fill histo array with event information only
  if (fHistoArray && fHistoArray->IsEventArray())
    fHistoArray->Fill(0,const_cast<Double_t *>(AliDielectronVarManager::GetData()),0x0,0x0);

  // clear arrays
  if (!fDontClearArrays) ClearArrays();

  // reset TPC EP and unique identifiers for v0 cut class
  AliDielectronVarManager::SetTPCEventPlane(0x0);
  if(GetHasMC()) { // only for MC needed
    for (Int_t iCut=0; iCut<fTrackFilter.GetCuts()->GetEntries();++iCut) {
      if ( fTrackFilter.GetCuts()->At(iCut)->IsA() == AliDielectronV0Cuts::Class() )
	((AliDielectronV0Cuts*)fTrackFilter.GetCuts()->At(iCut))->ResetUniqueEventNumbers();
    }
  }

  if (fTrackRotator) {
    fTrackRotator->ClearRotatedTrackPool();
    fTrackRotator->ClearRotatedPairPool();
  }

  return 1;

}

//________________________________________________________________
void AliDielectron::ProcessMC(AliVEvent *ev1)
{
  //
  // Process the MC data
  //

  AliDielectronMC *dieMC=AliDielectronMC::Instance();

  if (fHistos) FillHistogramsMC(dieMC->GetMCEvent(), ev1);

  // mc tracks
  if(!dieMC->GetNMCTracks()) return;

  // signals to be studied
  if(!fSignalsMC) return;
  Int_t nSignals = fSignalsMC->GetEntries();
  if(!nSignals) return;

  //loop over all MC data and Fill the HF, CF containers and histograms if they exist
  if(fCfManagerPair) fCfManagerPair->SetPdgMother(fPdgMother);

  Bool_t bFillCF   = (fCfManagerPair ? fCfManagerPair->GetStepForMCtruth()  : kFALSE);
  Bool_t bFillHF   = (fHistoArray    ? fHistoArray->GetStepForMCGenerated() : kFALSE);
  Bool_t bFillHist = kFALSE;
  if(fHistos) {
    const THashList *histlist =  fHistos->GetHistogramList();
    for(Int_t isig=0;isig<nSignals;isig++) {
      TString sigName = fSignalsMC->At(isig)->GetName();
      bFillHist |= histlist->FindObject(Form("Pair_%s_MCtruth",sigName.Data()))!=0x0;
      bFillHist |= histlist->FindObject(Form("Track_Leg_%s_MCtruth",sigName.Data()))!=0x0;
      bFillHist |= histlist->FindObject(Form("Track_%s_%s_MCtruth",fgkPairClassNames[1],sigName.Data()))!=0x0;
      if(bFillHist) break;
    }
  }
  // check if there is anything to fill
  if(!bFillCF && !bFillHF && !bFillHist) return;


  // initialize 2D arrays of labels for particles from each MC signal
  Int_t** labels1;      // labels for particles satisfying branch 1
  Int_t** labels2;      // labels for particles satisfying branch 2
  Int_t** labels12;     // labels for particles satisfying both branches
  labels1 = new Int_t*[nSignals];
  labels2 = new Int_t*[nSignals];
  labels12 = new Int_t*[nSignals];
  Int_t* indexes1=new Int_t[nSignals];
  Int_t* indexes2=new Int_t[nSignals];
  Int_t* indexes12=new Int_t[nSignals];
  for(Int_t isig=0;isig<nSignals;++isig) {
    *(labels1+isig) = new Int_t[dieMC->GetNMCTracks()];
    *(labels2+isig) = new Int_t[dieMC->GetNMCTracks()];
    *(labels12+isig) = new Int_t[dieMC->GetNMCTracks()];
    for(Int_t ip=0; ip<dieMC->GetNMCTracks();++ip) {
      labels1[isig][ip] = -1;
      labels2[isig][ip] = -1;
      labels12[isig][ip] = -1;
    }
    indexes1[isig]=0;
    indexes2[isig]=0;
    indexes12[isig]=0;
  }

  Bool_t truth1=kFALSE;
  Bool_t truth2=kFALSE;
  // loop over the MC tracks
  for(Int_t ipart=0; ipart<dieMC->GetNMCTracks(); ++ipart) {
    for(Int_t isig=0; isig<nSignals; ++isig) {       // loop over signals
      // Proceed only if this signal is required in the pure MC step
      // NOTE: Some signals can be satisfied by many particles and this leads to high
      //       computation times (e.g. secondary electrons from the GEANT transport). Be aware of this!!
      if(!((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetFillPureMCStep()) continue;

      truth1 = dieMC->IsMCTruth(ipart, (AliDielectronSignalMC*)fSignalsMC->At(isig), 1);
      truth2 = dieMC->IsMCTruth(ipart, (AliDielectronSignalMC*)fSignalsMC->At(isig), 2);

      // particles satisfying both branches are treated separately to avoid double counting during pairing
      if(truth1 && truth2) {
	labels12[isig][indexes12[isig]] = ipart;
	++indexes12[isig];
      }
      else {
	if(truth1) {
	  labels1[isig][indexes1[isig]] = ipart;
	  ++indexes1[isig];
	}
	if(truth2) {
	  labels2[isig][indexes2[isig]] = ipart;
	  ++indexes2[isig];
	}
      }
    }
  }  // end loop over MC particles

  // Do the pairing and fill the CF container with pure MC info
  for(Int_t isig=0; isig<nSignals; ++isig) {
    //    printf("INDEXES: %d-%d both%d\n",indexes1[isig],indexes2[isig],indexes12[isig]);
    // mix the particles which satisfy only one of the signal branches
    for(Int_t i1=0;i1<indexes1[isig];++i1) {
      if(!indexes2[isig]) FillMCHistograms(labels1[isig][i1], -1, isig); // (e.g. single electrons only, no pairs)
      for(Int_t i2=0;i2<indexes2[isig];++i2) {
	// add pair cuts on mc truth level
	if(bFillCF) fCfManagerPair->FillMC(labels1[isig][i1], labels2[isig][i2], isig);
	if(bFillHF) fHistoArray->Fill(labels1[isig][i1], labels2[isig][i2], isig);
	FillMCHistograms(labels1[isig][i1], labels2[isig][i2], isig);
      }
    }
    // mix the particles which satisfy both branches
    for(Int_t i1=0;i1<indexes12[isig];++i1) {
      for(Int_t i2=0; i2<i1; ++i2) {
	// add pair cuts on mc truth level
	if(bFillCF) fCfManagerPair->FillMC(labels12[isig][i1], labels12[isig][i2], isig);
	if(bFillHF) fHistoArray->Fill(labels12[isig][i1], labels12[isig][i2], isig);
	FillMCHistograms(labels12[isig][i1], labels12[isig][i2], isig);
      }
    }
  }    // end loop over signals

  // release the memory
  for(Int_t isig=0;isig<nSignals;++isig) {
    delete [] *(labels1+isig);
    delete [] *(labels2+isig);
    delete [] *(labels12+isig);
  }
  delete [] labels1;
  delete [] labels2;
  delete [] labels12;
  delete [] indexes1;
  delete [] indexes2;
  delete [] indexes12;
}

//________________________________________________________________
void AliDielectron::FillHistogramsTracks(TObjArray **tracks)
{
  //
  // Fill Histogram information for tracks after prefilter
  // ignore mixed events - for prefilter, only single tracks +/- are relevant
  //

  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);

  //Fill track information, separately for the track array candidates
  for (Int_t i=0; i<2; ++i){
    className.Form("Pre_%s",fgkTrackClassNames[i]);
    if (!fHistos->GetHistogramList()->FindObject(className.Data())) continue;
    Int_t ntracks=tracks[i]->GetEntriesFast();
    for (Int_t itrack=0; itrack<ntracks; ++itrack){
      AliDielectronVarManager::Fill(tracks[i]->UncheckedAt(itrack), values);
      fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
    }
  }
}


//________________________________________________________________
void AliDielectron::FillHistogramsMC(const AliMCEvent *ev, AliVEvent *ev1)
{
  //
  // Fill Histogram information for MCEvents
  //

  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars);

  // Fill event information
  AliDielectronVarManager::Fill(ev1, values);    // ESD/AOD information
  AliDielectronVarManager::Fill(ev, values);     // MC truth info
  if (fHistos->GetHistogramList()->FindObject("MCEvent"))
    fHistos->FillClass("MCEvent", AliDielectronVarManager::kNMaxValues, values);
}


//________________________________________________________________
void AliDielectron::FillHistograms(const AliVEvent *ev, Bool_t pairInfoOnly)
{
  //
  // Fill Histogram information for tracks and pairs
  //

  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars);

  //Fill event information
  if (ev){
    if (fHistos->GetHistogramList()->FindObject("Event")) {
      fHistos->FillClass("Event", AliDielectronVarManager::kNMaxValues, AliDielectronVarManager::GetData());
    }
  }

  //Fill track information, separately for the track array candidates
  if (!pairInfoOnly){
    className2.Form("Track_%s",fgkPairClassNames[1]);  // unlike sign, SE only
    for (Int_t i=0; i<6; ++i){
      className.Form("Track_%s",fgkTrackClassNames[i]);
      Bool_t mergedtrkClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;
      Bool_t trkClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
      if (!trkClass && !mergedtrkClass) continue;

      Double_t ntracks; 
      Double_t nPos = fTracks[0].GetEntriesFast();
      Double_t nNeg = fTracks[1].GetEntriesFast();
      //Int_t numberLS_PP = nPos*(nPos-1);
      //Int_t numberLS_MM = nNeg*(nNeg-1);
      if (i < 4) ntracks = fTracks[i].GetEntriesFast();
      else if ( i == 4) ntracks = fTrackRotator->GetRotatedTrackPSize();
      else if ( i == 5) ntracks = fTrackRotator->GetRotatedTrackNSize();
      

      for (Int_t itrack=0; itrack<ntracks; ++itrack){
        if (i < 4) {
          AliDielectronVarManager::Fill(fTracks[i].UncheckedAt(itrack), values);
        }
        else if (i == 4){
          AliKFParticle* part = fTrackRotator->GetRotatedTrackP(itrack);
          //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks, fTrackRotator->GetRotatedTrackWeightP(itrack) * GetWeightFromRotation(part) * nPos / ntracks);
          AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks,  GetWeightFromRotation(part) * nPos / ntracks);

          AliDielectronVarManager::Fill(part, values);
        }
        else if (i == 5) {
          AliKFParticle* part = fTrackRotator->GetRotatedTrackN(itrack);
          //AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks, fTrackRotator->GetRotatedTrackWeightN(itrack) * GetWeightFromRotation(part) * nNeg / ntracks);
          AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks,  GetWeightFromRotation(part) * nNeg / ntracks);

          AliDielectronVarManager::Fill(part, values);
        }
        if(trkClass)
          fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
        if(mergedtrkClass && i<2)
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values); //only ev1
      }
    }
  }

  //Fill Pair information, separately for all pair candidate arrays and the legs
  TObjArray arrLegs(100);
  for (Int_t i=0; i<13; ++i){
    className.Form("Pair_%s",fgkPairClassNames[i]);
    className2.Form("Track_Legs_%s",fgkPairClassNames[i]);
    Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
    Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;
    if (!pairClass&&!legClass) continue;
    Int_t ntracks=PairArray(i)->GetEntriesFast();
    for (Int_t ipair=0; ipair<ntracks; ++ipair){
      AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(i)->UncheckedAt(ipair));

      //fill pair information
      if (pairClass){
        AliDielectronVarManager::Fill(pair, values);
        fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
      }

      //fill leg information, don't fill the information twice
      if (legClass){
        AliVParticle *d1=pair->GetFirstDaughterP();
        AliVParticle *d2=pair->GetSecondDaughterP();
        if (!arrLegs.FindObject(d1)){
          AliDielectronVarManager::Fill(d1, values);
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
          arrLegs.Add(d1);
        }
        if (!arrLegs.FindObject(d2)){
          AliDielectronVarManager::Fill(d2, values);
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
          arrLegs.Add(d2);
        }
      }
    }
    if (legClass) arrLegs.Clear();
  }

}

//________________________________________________________________
void AliDielectron::FillHistogramsPair(AliDielectronPair *pair,Bool_t fromPreFilter/*=kFALSE*/)
{
  //
  // Fill Histogram information for pairs and the track in the pair
  // NOTE: in this funtion the leg information may be filled multiple
  //       times. This funtion is used in the track rotation pairing
  //       and those legs are not saved!
  //
  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);

  //Fill Pair information, separately for all pair candidate arrays and the legs
  TObjArray arrLegs(100);
  const Int_t type=pair->GetType();
  if (fromPreFilter) {
    className.Form("RejPair_%s",fgkPairClassNames[type]);
    className2.Form("RejTrack_%s",fgkPairClassNames[type]);
  } else {
    className.Form("Pair_%s",fgkPairClassNames[type]);
    className2.Form("Track_Legs_%s",fgkPairClassNames[type]);
  }

  Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
  Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;

  //fill pair information
  if (pairClass){
    AliDielectronVarManager::Fill(pair, values);
    fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
  }

  if (legClass){
    AliVParticle *d1=pair->GetFirstDaughterP();
    AliDielectronVarManager::Fill(d1, values);
    fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);

    AliVParticle *d2=pair->GetSecondDaughterP();
    AliDielectronVarManager::Fill(d2, values);
    fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
  }
}

//________________________________________________________________
void AliDielectron::FillTrackArrays(AliVEvent * const ev, Int_t eventNr)
{
  //
  // select tracks and fill track candidate arrays
  // eventNr = 0: First  event, use track arrays 0 and 1
  // eventNr = 1: Second event, use track arrays 2 and 3
  //

  Int_t ntracks=ev->GetNumberOfTracks();

  UInt_t selectedMask=(1<<fTrackFilter.GetCuts()->GetEntries())-1;
  for (Int_t itrack=0; itrack<ntracks; ++itrack){
    //get particle
    AliVParticle *particle=ev->GetTrack(itrack);

    //apply track cuts
    UInt_t cutmask=fTrackFilter.IsSelected(particle);
    //fill cut QA
    if(fCutQA) fQAmonitor->FillAll(particle);
    if(fCutQA) fQAmonitor->Fill(cutmask,particle);

    if (cutmask!=selectedMask) continue;

    //fill selected particle into the corresponding track arrays
    Short_t charge=particle->Charge();
    if (charge>0)      fTracks[eventNr*2].Add(particle);
    else if (charge<0) fTracks[eventNr*2+1].Add(particle);


  }
}

//________________________________________________________________
void AliDielectron::EventPlanePreFilter(Int_t arr1, Int_t arr2, TObjArray arrTracks1, TObjArray arrTracks2, const AliVEvent *ev)
{
  //
  // Prefilter tracks and tracks from pairs
  // Needed for rejection in the Q-Vector of the event plane
  // remove contribution of all tracks to the Q-vector that are in invariant mass window
  //

  // Run1 eventplane correction obsolete for run2 (due to qn framework) variables commented out in VarManager for memory savior, therefore also this function is commented out 20180312 PD
  // AliEventplane *evplane = const_cast<AliVEvent *>(ev)->GetEventplane();
  // if(!evplane) { // nanoAODs , here we do NOT have sub event reaction planes
  //   //  if(1) {
  //   // get the EPselectionTask for recalculation of weighting factors
  //   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //   AliEPSelectionTask *eptask = dynamic_cast<AliEPSelectionTask *>(man->GetTask("EventplaneSelection"));
  //   if(!eptask) return;
  //
  //   // get recentering values for Qx and Qy (only rms is needed for correction)
  //   // Double_t mean[2]={0.,0.}
  //   // eptask->Recenter(0, mean);
  //   Double_t rms[2] ={1.,1.};
  //   eptask->Recenter(1, rms);
  //
  //
  //   // track mapping
  //   TMap mapRemovedTracks;
  //
  //
  //   Double_t cQX=0., cQY=0.;
  //   // apply cuts to the tracks, e.g. etagap
  //   if(fEventPlanePreFilter.GetCuts()->GetEntries()) {
  //     UInt_t selectedMask=(1<<fEventPlanePreFilter.GetCuts()->GetEntries())-1;
  //     Int_t ntracks=ev->GetNumberOfTracks();
  //     for (Int_t itrack=0; itrack<ntracks; ++itrack){
  //       AliVParticle *particle=ev->GetTrack(itrack);
  //       AliVTrack *track= static_cast<AliVTrack*>(particle);
  //       if (!track) continue;
  //       //event plane cuts
  //       UInt_t cutMask=fEventPlanePreFilter.IsSelected(track);
  //       //apply cut
  //       if (cutMask==selectedMask) continue;
  //
  //       mapRemovedTracks.Add(track,track);
  //       cQX += (eptask->GetWeight(track) * TMath::Cos(2*track->Phi()) / rms[0]);
  //       cQY += (eptask->GetWeight(track) * TMath::Sin(2*track->Phi()) / rms[1]);
  //     }
  //   }
  //
  //   // POI (particle of interest) rejection
  //   Int_t pairIndex=GetPairIndex(arr1,arr2);
  //
  //   Int_t ntrack1=arrTracks1.GetEntriesFast();
  //   Int_t ntrack2=arrTracks2.GetEntriesFast();
  //   AliDielectronPair candidate;
  //   candidate.SetKFUsage(fUseKF);
  //
  //   UInt_t selectedMask=(1<<fEventPlanePOIPreFilter.GetCuts()->GetEntries())-1;
  //   for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
  //     Int_t end=ntrack2;
  //     if (arr1==arr2) end=itrack1;
  //     Bool_t accepted=kFALSE;
  //     for (Int_t itrack2=0; itrack2<end; ++itrack2){
  //       TObject *track1=arrTracks1.UncheckedAt(itrack1);
  //       TObject *track2=arrTracks2.UncheckedAt(itrack2);
  //       if (!track1 || !track2) continue;
  //       //create the pair
  //       candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
  //                           static_cast<AliVTrack*>(track2), fPdgLeg2);
  //       candidate.SetType(pairIndex);
  //       candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
  //
  //       //event plane pair cuts
  //       UInt_t cutMask=fEventPlanePOIPreFilter.IsSelected(&candidate);
  //       //apply cut
  //       if (cutMask==selectedMask) continue;
  //
  //       accepted=kTRUE;
  //       //remove the tracks from the Track arrays
  //       arrTracks2.AddAt(0x0,itrack2);
  //     }
  //     if ( accepted ) arrTracks1.AddAt(0x0,itrack1);
  //   }
  //   //compress the track arrays
  //   arrTracks1.Compress();
  //   arrTracks2.Compress();
  //
  //   //Modify the components: subtract the tracks
  //   ntrack1=arrTracks1.GetEntriesFast();
  //   ntrack2=arrTracks2.GetEntriesFast();
  //   // remove leg1 contribution
  //   for (Int_t itrack=0; itrack<ntrack1; ++itrack){
  //     AliVTrack *track= static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack));
  //     if (!track) continue;
  //     // track contribution was already removed
  //     if (mapRemovedTracks.FindObject(track)) continue;
  //     else mapRemovedTracks.Add(track,track);
  //
  //     cQX += (eptask->GetWeight(track) * TMath::Cos(2*track->Phi()) / rms[0]);
  //     cQY += (eptask->GetWeight(track) * TMath::Sin(2*track->Phi()) / rms[1]);
  //   }
  //   // remove leg2 contribution
  //   for (Int_t itrack=0; itrack<ntrack2; ++itrack){
  //     AliVTrack *track= static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack));
  //     if (!track) continue;
  //     // track contribution was already removed
  //     if (mapRemovedTracks.FindObject(track)) continue;
  //     else mapRemovedTracks.Add(track,track);
  //
  //     cQX += (eptask->GetWeight(track) * TMath::Cos(2*track->Phi()) / rms[0]);
  //     cQY += (eptask->GetWeight(track) * TMath::Sin(2*track->Phi()) / rms[1]);
  //   }
  //
  //   // build a corrected alieventplane using the values from the var manager
  //   // these uncorrected values are filled using the stored magnitude and angle  in the header
  //   TVector2 qcorr;
  //   qcorr.Set(AliDielectronVarManager::GetValue(AliDielectronVarManager::kTPCxH2uc)-cQX,
	//       AliDielectronVarManager::GetValue(AliDielectronVarManager::kTPCyH2uc)-cQY);
  //   // fill alieventplane
  //   AliEventplane cevplane;
  //   cevplane.SetQVector(&qcorr);
  //   AliDielectronVarManager::SetTPCEventPlane(&cevplane);
  //   cevplane.SetQVector(0);
  //   return;
  // } //end: nanoAODs
  // else
  //   {
  //   // this is done in case of ESDs or AODs
  //   Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());
  //   // copy event plane object
  //   AliEventplane cevplane(*evplane);
  //   //    Int_t nMaxID=cevplane->GetQContributionXArray()->GetSize();
  //
  //   TVector2 *qcorr  = cevplane.GetQVector();
  //   if(!qcorr) return;
  //   TVector2 *qcsub1 = 0x0;
  //   TVector2 *qcsub2 = 0x0;
  //
  //   // eta gap ?
  //   Bool_t etagap = kFALSE;
  //   for (Int_t iCut=0; iCut<fEventPlanePreFilter.GetCuts()->GetEntries();++iCut) {
  //     TString cutName=fEventPlanePreFilter.GetCuts()->At(iCut)->GetName();
  //     if(cutName.Contains("eta") || cutName.Contains("Eta"))  etagap=kTRUE;
  //   }
  //
  //   // subevent configuration for eta gap or LS (default is rndm)
  //   if(fLikeSignSubEvents && etagap) {
  //     // start with the full Qvector/event in both sub events
  //     qcsub1 = new TVector2(*qcorr);
  //     qcsub2 = new TVector2(*qcorr);
  //     cevplane.SetQsub(qcsub1,qcsub2);
  //
  //     Int_t ntracks=ev->GetNumberOfTracks();
  //     // track removals
  //     for (Int_t itrack=0; itrack<ntracks; ++itrack){
  //       AliVParticle *particle=ev->GetTrack(itrack);
  //       AliVTrack *track= static_cast<AliVTrack*>(particle);
  //       if (!track) continue;
  //       if (track->GetID()>=0 && !isESD) continue;
  //       Int_t tmpID = isESD ? track->GetID() : track->GetID()*-1 - 1;
  //
  //       // set contributions to zero
  //       // charge sub1+ sub2-
  //       if(fLikeSignSubEvents) {
  //         Short_t charge=track->Charge();
  //         if (charge<0) {
  //           cevplane.GetQContributionXArraysub1()->SetAt(0.0, tmpID);
  //           cevplane.GetQContributionYArraysub1()->SetAt(0.0, tmpID);
  //         }
  //         if (charge>0) {
  //           cevplane.GetQContributionXArraysub2()->SetAt(0.0, tmpID);
  //           cevplane.GetQContributionYArraysub2()->SetAt(0.0, tmpID);
  //         }
  //       }
  //       // eta sub1+ sub2-
  //       if(etagap) {
  //         Double_t eta=track->Eta();
  //         if (eta<0.0) {
  //           cevplane.GetQContributionXArraysub1()->SetAt(0.0, tmpID);
  //           cevplane.GetQContributionYArraysub1()->SetAt(0.0, tmpID);
  //         }
  //         if (eta>0.0) {
  //           cevplane.GetQContributionXArraysub2()->SetAt(0.0, tmpID);
  //           cevplane.GetQContributionYArraysub2()->SetAt(0.0, tmpID);
  //         }
  //       }
  //     } // end: loop over tracks
  //   } // end: sub event configuration
  //
  //   // apply cuts, e.g. etagap
  //   if(fEventPlanePreFilter.GetCuts()->GetEntries()) {
  //     UInt_t selectedMask=(1<<fEventPlanePreFilter.GetCuts()->GetEntries())-1;
  //     Int_t ntracks=ev->GetNumberOfTracks();
  //     for (Int_t itrack=0; itrack<ntracks; ++itrack){
  //       AliVParticle *particle=ev->GetTrack(itrack);
  //       AliVTrack *track= static_cast<AliVTrack*>(particle);
  //       if (!track) continue;
  //       if (track->GetID()>=0 && !isESD) continue;
  //       Int_t tmpID = isESD ? track->GetID() : track->GetID()*-1 - 1;
  //
  //       //event plane cuts
  //       UInt_t cutMask=fEventPlanePreFilter.IsSelected(track);
  //       //apply cut
  //       if (cutMask==selectedMask) continue;
  //
  //       // set contributions to zero
  //       cevplane.GetQContributionXArray()->SetAt(0.0, tmpID);
  //       cevplane.GetQContributionYArray()->SetAt(0.0, tmpID);
  //       cevplane.GetQContributionXArraysub1()->SetAt(0.0, tmpID);
  //       cevplane.GetQContributionYArraysub1()->SetAt(0.0, tmpID);
  //       cevplane.GetQContributionXArraysub2()->SetAt(0.0, tmpID);
  //       cevplane.GetQContributionYArraysub2()->SetAt(0.0, tmpID);
  //     }
  //   } // end: track cuts
  //
  //   // POI (particle of interest) rejection
  //   Int_t pairIndex=GetPairIndex(arr1,arr2);
  //   Int_t ntrack1=arrTracks1.GetEntriesFast();
  //   Int_t ntrack2=arrTracks2.GetEntriesFast();
  //   AliDielectronPair candidate;
  //   candidate.SetKFUsage(fUseKF);
  //
  //   UInt_t selectedMask=(1<<fEventPlanePOIPreFilter.GetCuts()->GetEntries())-1;
  //   for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
  //     Int_t end=ntrack2;
  //     if (arr1==arr2) end=itrack1;
  //     Bool_t accepted=kFALSE;
  //     for (Int_t itrack2=0; itrack2<end; ++itrack2){
  //       TObject *track1=arrTracks1.UncheckedAt(itrack1);
  //       TObject *track2=arrTracks2.UncheckedAt(itrack2);
  //       if (!track1 || !track2) continue;
  //       //create the pair
  //       candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
  //                           static_cast<AliVTrack*>(track2), fPdgLeg2);
  //
  //       candidate.SetType(pairIndex);
  //       candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
  //
  //       //event plane cuts
  //       UInt_t cutMask=fEventPlanePOIPreFilter.IsSelected(&candidate);
  //       //apply cut
  //       if (cutMask==selectedMask) continue;
  //
  //       accepted=kTRUE;
  //       //remove the tracks from the Track arrays
  //       arrTracks2.AddAt(0x0,itrack2);
  //     }
  //     if ( accepted ) arrTracks1.AddAt(0x0,itrack1);
  //   }
  //   //compress the track arrays
  //   arrTracks1.Compress();
  //   arrTracks2.Compress();
  //
  //   //Modify the components: subtract the tracks
  //   ntrack1=arrTracks1.GetEntriesFast();
  //   ntrack2=arrTracks2.GetEntriesFast();
  //   // remove leg1 contribution
  //   for (Int_t itrack=0; itrack<ntrack1; ++itrack){
  //     AliVTrack *track= static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack));
  //     if (!track) continue;
  //     if (track->GetID()>=0 && !isESD) continue;
  //     Int_t tmpID = isESD ? track->GetID() : track->GetID()*-1 - 1;
  //     // set contributions to zero
  //     cevplane.GetQContributionXArray()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionYArray()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionXArraysub1()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionYArraysub1()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionXArraysub2()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionYArraysub2()->SetAt(0.0, tmpID);
  //   }
  //   // remove leg2 contribution
  //   for (Int_t itrack=0; itrack<ntrack2; ++itrack){
  //     AliVTrack *track= static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack));
  //     if (!track) continue;
  //     if (track->GetID()>=0 && !isESD) continue;
  //     Int_t tmpID = isESD ? track->GetID() : track->GetID()*-1 - 1;
  //     // set contributions to zero
  //     cevplane.GetQContributionXArray()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionYArray()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionXArraysub1()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionYArraysub1()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionXArraysub2()->SetAt(0.0, tmpID);
  //     cevplane.GetQContributionYArraysub2()->SetAt(0.0, tmpID);
  //   }
  //
  //   // set corrected AliEventplane and fill variables with corrected values
  //   AliDielectronVarManager::SetTPCEventPlane(&cevplane);
  //   delete qcsub1;
  //   delete qcsub2;
  // } // end: ESD or AOD case

}

//________________________________________________________________
void AliDielectron::PairPreFilter(Int_t arr1, Int_t arr2, TObjArray &arrTracks1, TObjArray &arrTracks2, const AliVEvent *ev, Int_t prefilterN )
{
  //
  // Prefilter tracks from pairs
  // Needed for Dalitz rejections
  // remove all tracks from the Single track arrays that pass the cuts in this filter
  //

  AliAnalysisFilter * pairPreFilter = prefilterN == 1 ? &fPairPreFilter1 : &fPairPreFilter2;
  Bool_t prefilterAllSigns = prefilterN == 1 ? fPreFilterAllSigns1 : fPreFilterAllSigns2;
  Bool_t prefilterPhotons = prefilterN == 1 ? fPreFilterPhotons1 : fPreFilterPhotons2;
  Bool_t prefilterOnlyOnePair = prefilterN == 1 ? fPreFilterOnlyOnePair1 : fPreFilterOnlyOnePair2;

  Int_t ntrack1=arrTracks1.GetEntriesFast();
  Int_t ntrack2=arrTracks2.GetEntriesFast();
  AliDielectronPair candidate;
  candidate.SetKFUsage(fUseKF);
  // flag arrays for track removal
  Bool_t *bTracks1 = new Bool_t[ntrack1];
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1) bTracks1[itrack1]=kFALSE;
  Bool_t *bTracks2 = new Bool_t[ntrack2];
  for (Int_t itrack2=0; itrack2<ntrack2; ++itrack2) bTracks2[itrack2]=kFALSE;

  UInt_t selectedMask= (1<<pairPreFilter->GetCuts()->GetEntries())-1 ;
  UInt_t selectedMaskPair=(1<<fPairFilter.GetCuts()->GetEntries())-1;

  Int_t nRejPasses = 1; //for fPreFilterUnlikeOnly and no set flag
  if (prefilterAllSigns) nRejPasses = 3;

  for (Int_t iRP=0; iRP < nRejPasses; ++iRP) {
    Int_t arr1RP=arr1, arr2RP=arr2;
    TObjArray *arrTracks1RP=&arrTracks1;
    TObjArray *arrTracks2RP=&arrTracks2;
    Bool_t *bTracks1RP = bTracks1;
    Bool_t *bTracks2RP = bTracks2;
    switch (iRP) {
      case 1: arr1RP=arr1;arr2RP=arr1;
				arrTracks1RP=&arrTracks1;
				arrTracks2RP=&arrTracks1;
				bTracks1RP = bTracks1;
				bTracks2RP = bTracks1;
				break;
      case 2: arr1RP=arr2;arr2RP=arr2;
				arrTracks1RP=&arrTracks2;
				arrTracks2RP=&arrTracks2;
				bTracks1RP = bTracks2;
				bTracks2RP = bTracks2;
				break;
      default: ;//nothing to do
    }
    Int_t ntrack1RP=(*arrTracks1RP).GetEntriesFast();
    Int_t ntrack2RP=(*arrTracks2RP).GetEntriesFast();

    Int_t pairIndex=GetPairIndex(arr1RP,arr2RP);

    if( prefilterOnlyOnePair ){
      Double_t maxLikelihood1[ntrack1RP];
      Double_t maxLikelihood2[ntrack2RP];
      Int_t partner1[ntrack1RP];
      Int_t partner2[ntrack2RP];
      for(Int_t itrack2=0; itrack2<ntrack2RP; ++itrack2){
        maxLikelihood2[itrack2] = -999.;
        partner2[itrack2] = -1;
      }

      for (Int_t itrack1=0; itrack1<ntrack1RP; ++itrack1){
        maxLikelihood1[itrack1] = -999.;
        partner1[itrack1] = -1;
        Int_t end=ntrack2RP;
        if (arr1RP==arr2RP) end=itrack1;
        for (Int_t itrack2=0; itrack2<end; ++itrack2){
          TObject *track1=(*arrTracks1RP).UncheckedAt(itrack1);
          TObject *track2=(*arrTracks2RP).UncheckedAt(itrack2);
          if (!track1 || !track2) continue;
          maxLikelihood2[itrack2] = -999.;
          //create the pair
          if(prefilterPhotons){
            candidate.SetGammaTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
                              static_cast<AliVTrack*>(track2), fPdgLeg2);
          }
          else{
            candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
                              static_cast<AliVTrack*>(track2), fPdgLeg2);
          }

          candidate.SetType(pairIndex);
          candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
          //relate to the production vertex
          //       if (AliDielectronVarManager::GetKFVertex()) candidate.SetProductionVertex(*AliDielectronVarManager::GetKFVertex());

          //pair cuts
          UInt_t cutMask=pairPreFilter->IsSelected(&candidate);

          //apply cut
          if (cutMask!=selectedMask) continue;
      //    Double_t likelihood = prefilterN == 1 ?  candidate.PhivPair(ev->GetMagneticField()) - 21. * candidate.M() : -1. * candidate.M();
            Double_t likelihood =  candidate.GetKFNdf() / candidate.GetKFChi2();
          if( likelihood > maxLikelihood1[itrack1]   && likelihood > maxLikelihood2[itrack2]  ){
            maxLikelihood1[itrack1] = likelihood;
            maxLikelihood2[itrack2] = likelihood;
            partner1[itrack1] = itrack2;
            partner2[itrack2] = itrack1;
          }

          if (fCfManagerPair) fCfManagerPair->Fill(selectedMaskPair+1 ,&candidate);
          if (fHistos) FillHistogramsPair(&candidate,kTRUE);
          //set flags for track removal
        }
      }
      for (Int_t itrack1=0; itrack1<ntrack1RP; ++itrack1){
        if(   partner1[itrack1] != -1     ){
          Int_t itrack2 = partner1[itrack1];
          if(  partner2[ itrack2  ] == itrack1    ){
            bTracks1RP[itrack1]=kTRUE;
            bTracks2RP[itrack2]=kTRUE;
          }
        }
      }
    }
    else{
      for (Int_t itrack1=0; itrack1<ntrack1RP; ++itrack1){
        Int_t end=ntrack2RP;
        if (arr1RP==arr2RP) end=itrack1;
        for (Int_t itrack2=0; itrack2<end; ++itrack2){
          TObject *track1=(*arrTracks1RP).UncheckedAt(itrack1);
          TObject *track2=(*arrTracks2RP).UncheckedAt(itrack2);
          if (!track1 || !track2) continue;
          //create the pair
          if(prefilterPhotons){
            candidate.SetGammaTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
                              static_cast<AliVTrack*>(track2), fPdgLeg2);
          }
          else{
            candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
                              static_cast<AliVTrack*>(track2), fPdgLeg2);
          }

          candidate.SetType(pairIndex);
          candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
          //relate to the production vertex
          //       if (AliDielectronVarManager::GetKFVertex()) candidate.SetProductionVertex(*AliDielectronVarManager::GetKFVertex());

          //pair cuts
          UInt_t cutMask=pairPreFilter->IsSelected(&candidate);

          //apply cut
          if (cutMask!=selectedMask) continue;
          if (fCfManagerPair) fCfManagerPair->Fill(selectedMaskPair+1 ,&candidate);
          if (fHistos) FillHistogramsPair(&candidate,kTRUE);
          //set flags for track removal
          bTracks1RP[itrack1]=kTRUE;
          bTracks2RP[itrack2]=kTRUE;
        }
      }
    }
  }

  //remove the tracks from the Track arrays
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
    if(bTracks1[itrack1]) {
      arrTracks1.AddAt(0x0, itrack1);
      if(arr1 == arr2) arrTracks2.AddAt(0x0, itrack1);
    }
  }
  for (Int_t itrack2=0; itrack2<ntrack2; ++itrack2){
    if(bTracks2[itrack2]) {
      arrTracks2.AddAt(0x0, itrack2);
      if(arr1 == arr2) arrTracks1.AddAt(0x0, itrack2);
    }
  }

  // clean up
  delete [] bTracks1;
  delete [] bTracks2;

  //compress the track arrays
  arrTracks1.Compress();
  arrTracks2.Compress();

  AliAnalysisFilter * pairPreFilterLegs = prefilterN == 1 ? &fPairPreFilterLegs1 : &fPairPreFilterLegs2;
  //apply leg cuts after the pre filter
  if ( pairPreFilterLegs->GetCuts()->GetEntries()>0 ) {
    selectedMask=(1<<pairPreFilterLegs->GetCuts()->GetEntries())-1;
    //loop over tracks from array 1
    for (Int_t itrack=0; itrack<arrTracks1.GetEntriesFast();++itrack){
      //test cuts
      UInt_t cutMask=pairPreFilterLegs->IsSelected(arrTracks1.UncheckedAt(itrack));

      //apply cut
      if (cutMask!=selectedMask) arrTracks1.AddAt(0x0,itrack);
    }
    arrTracks1.Compress();

    //in case of like sign don't loop over second array
    if (arr1==arr2) {
      arrTracks2=arrTracks1;
    } else {

      //loop over tracks from array 2
      for (Int_t itrack=0; itrack<arrTracks2.GetEntriesFast();++itrack){
      //test cuts
        UInt_t cutMask=pairPreFilterLegs->IsSelected(arrTracks2.UncheckedAt(itrack));
      //apply cut
        if (cutMask!=selectedMask) arrTracks2.AddAt(0x0,itrack);
      }
      arrTracks2.Compress();

    }
  }
  //For unlike-sign monitor track-cuts:
  if (arr1!=arr2&&fHistos) {
    TObjArray *unlikesignArray[2] = {&arrTracks1,&arrTracks2};
    FillHistogramsTracks(unlikesignArray);
  }
}

//________________________________________________________________
void AliDielectron::FillPairArrays(Int_t arr1, Int_t arr2, const AliVEvent *ev)
{
  //
  // select pairs and fill pair candidate arrays
  //

  TObjArray arrTracks1=fTracks[arr1];
  TObjArray arrTracks2=fTracks[arr2];

  //process pre filter if set
  if ((!fPreFilterAllSigns1) && (!fPreFilterUnlikeOnly1) && (!fPreFilterLikeOnly1) && ( fPairPreFilter1.GetCuts()->GetEntries()>0 ))  PairPreFilter(arr1, arr2, arrTracks1, arrTracks2, ev, 1);

  if ((!fPreFilterAllSigns2) && (!fPreFilterUnlikeOnly2) && (!fPreFilterLikeOnly2) && ( fPairPreFilter2.GetCuts()->GetEntries()>0 ))  PairPreFilter(arr1, arr2, arrTracks1, arrTracks2, ev, 2);

  Int_t pairIndex=GetPairIndex(arr1,arr2);

  Int_t ntrack1=arrTracks1.GetEntriesFast();
  Int_t ntrack2=arrTracks2.GetEntriesFast();

  AliDielectronPair *candidate=new AliDielectronPair;
  candidate->SetKFUsage(fUseKF);

  UInt_t selectedMask=(1<<fPairFilter.GetCuts()->GetEntries())-1;

  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1){
    Int_t end=ntrack2;
    if (arr1==arr2) end=itrack1;
    for (Int_t itrack2=0; itrack2<end; ++itrack2){
      //create the pair (direct pointer to the memory by this daughter reference are kept also for ME)
      candidate->SetTracks(&(*static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack1))), fPdgLeg1,
                           &(*static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack2))), fPdgLeg2);
      candidate->SetType(pairIndex);

      Int_t label=AliDielectronMC::Instance()->GetLabelMotherWithPdg(candidate,fPdgMother);
      candidate->SetLabel(label);
      if (label>-1) candidate->SetPdgCode(fPdgMother);
      else candidate->SetPdgCode(0);

      // check for gamma kf particle
      label=AliDielectronMC::Instance()->GetLabelMotherWithPdg(candidate,22);
      if (label>-1 && fUseGammaTracks) {
        candidate->SetGammaTracks(static_cast<AliVTrack*>(arrTracks1.UncheckedAt(itrack1)), fPdgLeg1,
                                  static_cast<AliVTrack*>(arrTracks2.UncheckedAt(itrack2)), fPdgLeg2);
      // should we set the pdgmothercode and the label
      }

      //pair cuts
      UInt_t cutMask=fPairFilter.IsSelected(candidate);

      //CF manager for the pair
      if (fCfManagerPair) fCfManagerPair->Fill(cutMask,candidate);

      // cut qa
      if(pairIndex==kEv1PM && fCutQA) {
        fQAmonitor->FillAll(candidate);
        fQAmonitor->Fill(cutMask,candidate);
      }

      //apply cut
      if (cutMask!=selectedMask) continue;

      //histogram array for the pair
      if (fHistoArray) fHistoArray->Fill(pairIndex,candidate);

      //add the candidate to the candidate array
      PairArray(pairIndex)->Add(candidate);
      //get a new candidate
      candidate=new AliDielectronPair;
      candidate->SetKFUsage(fUseKF);
    }
  }
  //delete the surplus candidate
  delete candidate;
}

//________________________________________________________________
void AliDielectron::FillPairArrayTR()
{
  //
  // select pairs and fill pair candidate arrays
  //
  UInt_t selectedMask=(1<<fPairFilter.GetCuts()->GetEntries())-1;
  while ( fTrackRotator->NextCombination() ){
    if(fTrackRotator->SameTracks() ) continue;
    AliDielectronPair candidate;
    candidate.SetKFUsage(fUseKF);
    candidate.SetTracks(&fTrackRotator->GetKFTrack1(), &fTrackRotator->GetKFTrack2(),
                        fTrackRotator->GetVTrack1(),fTrackRotator->GetVTrack2());

    Int_t label=AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother);
    candidate.SetLabel(label);
    if (label>-1) candidate.SetPdgCode(fPdgMother);
      else candidate.SetPdgCode(0);

    if ((fTrackRotator->GetChargeTrack1() == +1 && fTrackRotator->GetChargeTrack2() == -1) ||
        (fTrackRotator->GetChargeTrack1() == -1 && fTrackRotator->GetChargeTrack2() == +1)){
      candidate.SetType(kEv1PMRot);
      if(fStoreRotatedPairs) PairArray(kEv1PMRot)->Add(new AliDielectronPair(candidate));
    }
    else if (fTrackRotator->GetChargeTrack1() == +1 && fTrackRotator->GetChargeTrack2() == +1){
      candidate.SetType(kEv1PPRot);
      if(fStoreRotatedPairs) PairArray(kEv1PPRot)->Add(new AliDielectronPair(candidate));
    }
    else if (fTrackRotator->GetChargeTrack1() == -1 && fTrackRotator->GetChargeTrack2() == -1){
      candidate.SetType(kEv1MMRot);
      if(fStoreRotatedPairs) PairArray(kEv1MMRot)->Add(new AliDielectronPair(candidate));
    }
    else 
      continue;

    //pair cuts
    UInt_t cutMask=fPairFilter.IsSelected(&candidate);

    //CF manager for the pair
    if (fCfManagerPair) fCfManagerPair->Fill(cutMask,&candidate);

    //apply cut
    if (cutMask==selectedMask) {

      //histogram array for the pair
      if (fHistoArray) fHistoArray->Fill(candidate.GetType() ,&candidate);

      //if(fHistos) FillHistogramsPair(&candidate);
      //if(fStoreRotatedPairs) PairArray(kEv1PMRot)->Add(new AliDielectronPair(candidate));
    }
  }
}

//________________________________________________________________
void AliDielectron::FillDebugTree()
{
  //
  // Fill Histogram information for tracks and pairs
  //

  //Fill Debug tree
  AliDielectronVarManager::SetFillMap(fDebugTree->GetUsedVars());
  for (Int_t i=0; i<10; ++i){
    Int_t ntracks=PairArray(i)->GetEntriesFast();
    for (Int_t ipair=0; ipair<ntracks; ++ipair){
      fDebugTree->Fill(static_cast<AliDielectronPair*>(PairArray(i)->UncheckedAt(ipair)));
    }
  }
}

//________________________________________________________________
void AliDielectron::SaveDebugTree()
{
  //
  // delete the debug tree, this will also write the tree
  //
  if (fDebugTree) fDebugTree->DeleteStreamer();
}


//__________________________________________________________________
void AliDielectron::AddSignalMC(AliDielectronSignalMC* signal) {
  //
  //  Add an MC signal to the signals list
  //
  if(!fSignalsMC) {
    fSignalsMC = new TObjArray();
    fSignalsMC->SetOwner();
  }
  fSignalsMC->Add(signal);
}

//________________________________________________________________
void AliDielectron::FillMCHistograms(Int_t label1, Int_t label2, Int_t nSignal) {
  //
  // fill QA MC TRUTH histograms for pairs and legs of all added mc signals
  //

  TString className,className2,className3;
  className.Form("Pair_%s_MCtruth",fSignalsMC->At(nSignal)->GetName());
  className2.Form("Track_Legs_%s_MCtruth",fSignalsMC->At(nSignal)->GetName());
  className3.Form("Track_%s_%s_MCtruth",fgkPairClassNames[1],fSignalsMC->At(nSignal)->GetName());
  Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
  Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;
  Bool_t trkClass=fHistos->GetHistogramList()->FindObject(className3.Data())!=0x0;
  //  printf("fill signal %d: pair %d legs %d trk %d \n",nSignal,pairClass,legClass,trkClass);
  if(!pairClass && !legClass && !trkClass) return;

  //  printf("leg labels: %d-%d \n",label1,label2);
  AliVParticle* part1 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label1);
  AliVParticle* part2 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label2);
  if(!part1 && !part2) return;
  if(part1&&part2) {
    // fill only unlike sign (and only SE)
    if(part1->Charge()*part2->Charge()>=0) return;
  }


  AliDielectronMC* dieMC = AliDielectronMC::Instance();

  Int_t mLabel1 = dieMC->GetMothersLabel(label1);    // should work for both ESD and AOD
  Int_t mLabel2 = dieMC->GetMothersLabel(label2);

  // check the same mother option
  AliDielectronSignalMC* sigMC = (AliDielectronSignalMC*)fSignalsMC->At(nSignal);
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kSame && mLabel1!=mLabel2) return;
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent && mLabel1==mLabel2) return;

  // fill event values
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(dieMC->GetMCEvent(), values); // get event informations
  // @TODO: check if this Fill() is even needed. It might modify the fill map (fUsedVars).

  // fill the leg variables
  //  printf("leg:%d trk:%d part1:%p part2:%p \n",legClass,trkClass,part1,part2);
  if (legClass || trkClass) {
    if(part1) AliDielectronVarManager::Fill(part1,values);
    if(part1 && trkClass)          fHistos->FillClass(className3, AliDielectronVarManager::kNMaxValues, values);
    if(part1 && part2 && legClass) fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
    if(part2) AliDielectronVarManager::Fill(part2,values);
    if(part2 && trkClass)          fHistos->FillClass(className3, AliDielectronVarManager::kNMaxValues, values);
    if(part1 && part2 && legClass) fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
  }

  //fill pair information
  if (pairClass && part1 && part2) {
    AliDielectronVarManager::FillVarMCParticle2(part1,part2,values);
    fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
  }

}

//________________________________________________________________
void AliDielectron::FillMCHistograms(const AliVEvent *ev) {
  //
  // fill QA MC histograms for pairs and legs of all added mc signals
  //
  if (!fSignalsMC) return;
  TString className,className2,className3;
  TString className4,className5,className6;
  TString className4_2,className5_2,className6_2;
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars);
  // AliDielectronVarManager::Fill(ev, values);
  // not needed to get event information here, because done in FillVarVParticle() [and FillVarDielectronPair()].

  //loop over all added mc signals
  for(Int_t isig=0; isig<fSignalsMC->GetEntries(); isig++) {

    //check if and what to fill
    className.Form("Pair_%s",fSignalsMC->At(isig)->GetName());
    className2.Form("Track_Legs_%s",fSignalsMC->At(isig)->GetName());
    className3.Form("Track_%s_%s",fgkPairClassNames[1],fSignalsMC->At(isig)->GetName());  // unlike sign, SE only

    Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
    Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;
    Bool_t mergedtrkClass=fHistos->GetHistogramList()->FindObject(className3.Data())!=0x0;
    if(!pairClass && !legClass && !mergedtrkClass) continue;

    className4.Form("Pair_%s_ev1+_ev1-_TR",fSignalsMC->At(isig)->GetName());
    className5.Form("Pair_%s_ev1+_ev1+_TR",fSignalsMC->At(isig)->GetName());
    className6.Form("Pair_%s_ev1-_ev1-_TR",fSignalsMC->At(isig)->GetName());

    className4_2.Form("Track_Legs_%s_ev1+_ev1-_TR",fSignalsMC->At(isig)->GetName());
    className5_2.Form("Track_Legs_%s_ev1+_ev1+_TR",fSignalsMC->At(isig)->GetName());
    className6_2.Form("Track_Legs_%s_ev1-_ev1-_TR",fSignalsMC->At(isig)->GetName());

    Bool_t pairClass_ULS_TR=fHistos->GetHistogramList()->FindObject(className4.Data())!=0x0;
    Bool_t pairClass_LSpp_TR=fHistos->GetHistogramList()->FindObject(className5.Data())!=0x0;
    Bool_t pairClass_LSmm_TR=fHistos->GetHistogramList()->FindObject(className6.Data())!=0x0;
    Bool_t legClass_ULS_TR=fHistos->GetHistogramList()->FindObject(className4_2.Data())!=0x0;
    Bool_t legClass_LSpp_TR=fHistos->GetHistogramList()->FindObject(className5_2.Data())!=0x0;
    Bool_t legClass_LSmm_TR=fHistos->GetHistogramList()->FindObject(className6_2.Data())!=0x0;

    // fill pair and/or their leg variables
    if(pairClass || legClass) {
      if(((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetCheckUnlikeSign()){
        Int_t npairs=PairArray(AliDielectron::kEv1PM)->GetEntriesFast(); // SE +-
        for (Int_t ipair=0; ipair<npairs; ++ipair){
          AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(AliDielectron::kEv1PM)->UncheckedAt(ipair));
          Bool_t isMCtruth = AliDielectronMC::Instance()->IsMCTruth(pair, (AliDielectronSignalMC*)fSignalsMC->At(isig));
          if(isMCtruth) {
            //fill pair information
            if (pairClass){
              AliDielectronVarManager::Fill(pair, values);
              fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
            }
            //fill leg information, both + and - in the same histo
            if (legClass){
              AliDielectronVarManager::Fill(pair->GetFirstDaughterP(),values);
              fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
              AliDielectronVarManager::Fill(pair->GetSecondDaughterP(),values);
              fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
            }
          } //is signal
        } //loop: pairs
      } // fill +- pairs
      if(((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetCheckLikeSignPP()){
        Int_t npairsLS1=PairArray(AliDielectron::kEv1PP)->GetEntriesFast(); // SE ++
        for (Int_t ipair=0; ipair<npairsLS1; ++ipair){
          AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(AliDielectron::kEv1PP)->UncheckedAt(ipair));
          Bool_t isMCtruth = AliDielectronMC::Instance()->IsMCTruth(pair, (AliDielectronSignalMC*)fSignalsMC->At(isig));
          if(isMCtruth){
            //fill pair information
            if (pairClass){
              AliDielectronVarManager::Fill(pair, values);
              fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
            }
            //fill leg information, both + and - in the same histo
            if (legClass){
              AliDielectronVarManager::Fill(pair->GetFirstDaughterP(),values);
              fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
              AliDielectronVarManager::Fill(pair->GetSecondDaughterP(),values);
              fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
            }
          } //is signal
        } //loop: pairs
      } // fill ++ pairs
      if(((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetCheckLikeSignMM()){
        Int_t npairsLS2=PairArray(AliDielectron::kEv1MM)->GetEntriesFast(); // SE --
        for (Int_t ipair=0; ipair<npairsLS2; ++ipair){
          AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(AliDielectron::kEv1MM)->UncheckedAt(ipair));
          Bool_t isMCtruth = AliDielectronMC::Instance()->IsMCTruth(pair, (AliDielectronSignalMC*)fSignalsMC->At(isig));
          if(isMCtruth){
            //fill pair information
            if (pairClass){
              AliDielectronVarManager::Fill(pair, values);
              fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
            }
            //fill leg information, both + and - in the same histo
            if (legClass){
              AliDielectronVarManager::Fill(pair->GetFirstDaughterP(),values);
              fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
              AliDielectronVarManager::Fill(pair->GetSecondDaughterP(),values);
              fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
            }
          } //is signal
        } //loop: pairs
      } // fill -- pairs


      if (fTrackRotator) {
        if(((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetCheckUnlikeSign()){
          Int_t npairs=PairArray(AliDielectron::kEv1PMRot)->GetEntriesFast(); // SE +-
          for (Int_t ipair=0; ipair<npairs; ++ipair){
            AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(AliDielectron::kEv1PMRot)->UncheckedAt(ipair));
            Bool_t isMCtruth = AliDielectronMC::Instance()->IsMCTruth(pair, (AliDielectronSignalMC*)fSignalsMC->At(isig));
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kRotationAngle, fTrackRotator->GetRotatedPairPM(ipair)->rotAng); 
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks, fTrackRotator->GetRotatedPairPM(ipair)->weight); 
            if(isMCtruth) {
              //fill pair information
              if (pairClass_ULS_TR){
                AliDielectronVarManager::Fill(pair, values);
                fHistos->FillClass(className4, AliDielectronVarManager::kNMaxValues, values);
              }
              //fill leg information, both + and - in the same histo
              if (legClass_ULS_TR){
                AliDielectronVarManager::Fill(&(pair->GetKFFirstDaughter()),values);
                fHistos->FillClass(className4_2, AliDielectronVarManager::kNMaxValues, values);
                AliDielectronVarManager::Fill(&(pair->GetKFSecondDaughter()),values);
                fHistos->FillClass(className4_2, AliDielectronVarManager::kNMaxValues, values);
              }
            } //is signal
          } //loop: pairs
        } // fill +- pairs
        if(((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetCheckLikeSignPP()){
          Int_t npairsLS1=PairArray(AliDielectron::kEv1PPRot)->GetEntriesFast(); // SE ++
          for (Int_t ipair=0; ipair<npairsLS1; ++ipair){
            AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(AliDielectron::kEv1PPRot)->UncheckedAt(ipair));
            Bool_t isMCtruth = AliDielectronMC::Instance()->IsMCTruth(pair, (AliDielectronSignalMC*)fSignalsMC->At(isig));
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kRotationAngle, fTrackRotator->GetRotatedPairPP(ipair)->rotAng);
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks, fTrackRotator->GetRotatedPairPP(ipair)->weight); 
            if(isMCtruth){
              //fill pair information
              if (pairClass_LSpp_TR){
                AliDielectronVarManager::Fill(pair, values);
                fHistos->FillClass(className5, AliDielectronVarManager::kNMaxValues, values);
              }
              //fill leg information, both + and - in the same histo
              if (legClass_LSpp_TR){
                AliDielectronVarManager::Fill(&(pair->GetKFFirstDaughter()),values);
                fHistos->FillClass(className5_2, AliDielectronVarManager::kNMaxValues, values);
                AliDielectronVarManager::Fill(&(pair->GetKFSecondDaughter()),values);
                fHistos->FillClass(className5_2, AliDielectronVarManager::kNMaxValues, values);
              }
            } //is signal
          } //loop: pairs
        } // fill ++ pairs
        if(((AliDielectronSignalMC*)fSignalsMC->At(isig))->GetCheckLikeSignMM()){
          Int_t npairsLS2=PairArray(AliDielectron::kEv1MMRot)->GetEntriesFast(); // SE --
          for (Int_t ipair=0; ipair<npairsLS2; ++ipair){
            AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(AliDielectron::kEv1MMRot)->UncheckedAt(ipair));
            Bool_t isMCtruth = AliDielectronMC::Instance()->IsMCTruth(pair, (AliDielectronSignalMC*)fSignalsMC->At(isig));
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kRotationAngle, fTrackRotator->GetRotatedPairMM(ipair)->rotAng);
            AliDielectronVarManager::SetValue(AliDielectronVarManager::kWeightFromRotationSingleTracks, fTrackRotator->GetRotatedPairMM(ipair)->weight); 
            if(isMCtruth){
              //fill pair information
              if (pairClass_LSmm_TR){
                AliDielectronVarManager::Fill(pair, values);
                fHistos->FillClass(className6, AliDielectronVarManager::kNMaxValues, values);
              }
              //fill leg information, both + and - in the same histo
              if (legClass_LSmm_TR){
                AliDielectronVarManager::Fill(&(pair->GetKFFirstDaughter()),values);
                fHistos->FillClass(className6_2, AliDielectronVarManager::kNMaxValues, values);
                AliDielectronVarManager::Fill(&(pair->GetKFSecondDaughter()),values);
                fHistos->FillClass(className6_2, AliDielectronVarManager::kNMaxValues, values);
              }
            } //is signal
          } //loop: pairs
        } // fill -- pairs
      }

    }

    // fill single tracks of signals
    if(!mergedtrkClass) continue;
    // loop over SE track arrays
    for (Int_t i=0; i<2; ++i){
      Int_t ntracks=fTracks[i].GetEntriesFast();
      for (Int_t itrack=0; itrack<ntracks; ++itrack){
        Int_t label=((AliVParticle*)fTracks[i].UncheckedAt(itrack))->GetLabel();
        Bool_t isMCtruth1 = AliDielectronMC::Instance()->IsMCTruth(label, (AliDielectronSignalMC*)fSignalsMC->At(isig), 1);
        Bool_t isMCtruth2 = AliDielectronMC::Instance()->IsMCTruth(label, (AliDielectronSignalMC*)fSignalsMC->At(isig), 2);
        // skip if track does not correspond to the signal
        if(!isMCtruth1 && !isMCtruth2) continue;
        AliDielectronVarManager::Fill(fTracks[i].UncheckedAt(itrack), values);
        fHistos->FillClass(className3, AliDielectronVarManager::kNMaxValues, values);
      } //loop: tracks
    } //loop: arrays

  } //loop: MCsignals

}



//______________________________________________
void AliDielectron::SetCentroidCorrArr(TObjArray *arrFun, Bool_t bHisto, UInt_t varx, UInt_t vary, UInt_t varz)
{
  // Store the maps with the runnumber as name in the arrays
  // bHisto has to be setted to true if histograms and not functions are used for correction
  fPostPIDCntrdCorrArr = new TObjArray();
  Int_t nEntries = arrFun->GetEntriesFast();
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key;
  TH1 *histo;

  for (Int_t i = 0; i < nEntries; i++) {
    key = arrFun->At(i)->GetName();
    if(!key.IsDigit()){
      AliFatal("The stored functions in the PID correction array have a wrong naming scheme - please store them with the according run number as the name!");
      return;
    }
    key.Form("Centroid_%d",key.Atoi());
    // if TF1s are stored in the array
    if(!bHisto){
      TF1 *fun = (TF1*) arrFun->At(i);
      AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
      // clone temporare histogram, otherwise it will not be streamed to file!
      fPostPIDCntrdCorrArr->Add((TH1*)fun->GetHistogram()->Clone(key.Data()));
      histo = (TH1*) fPostPIDCntrdCorrArr->At(i);
      histo->GetListOfFunctions()->AddAt(fun, 0);
    }
    // if TH1s are stored in the array
    if(bHisto){
      TH1 *fun = (TH1*) arrFun->At(i);
      AliDielectronHistos::StoreVariables(fun, valType);
      // clone temporare histogram, otherwise it will not be streamed to file!
      fPostPIDCntrdCorrArr->Add((TH1*)fun->Clone(key.Data()));
      histo = (TH1*) fPostPIDCntrdCorrArr->At(i);
    }
    if(histo){
      // check for corrections and add their variables to the fill map
      printf("POST TPC PID CORRECTION run added for %s:  ",key.Data());
      switch(histo->GetDimension()) {
      case 3: printf(" %s, ",histo->GetZaxis()->GetName());
      case 2: printf(" %s, ",histo->GetYaxis()->GetName());
      case 1: printf(" %s ",histo->GetXaxis()->GetName());
      }
      printf("\n");
    }
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}

//______________________________________________
void AliDielectron::SetWidthCorrArr(TObjArray *arrFun, Bool_t bHisto, UInt_t varx, UInt_t vary, UInt_t varz)
{
  // Store the maps with the runnumber as name in the arrays
  // bHisto has to be setted to true if histograms and not functions are used for correction
  fPostPIDWdthCorrArr = new TObjArray();
  Int_t nEntries = arrFun->GetEntriesFast();
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key;
  TH1 *histo;
  for (Int_t i = 0; i < nEntries; i++) {
    key = arrFun->At(i)->GetName();
    if(!key.IsDigit()){
      AliFatal("The stored functions in the PID correction array have a wrong naming scheme - please store them with the according run number as the name!");
      return;
    }
    key.Form("Width_%d",key.Atoi());
    // if TF1s are stored in the array
    if(!bHisto){
      TF1 *fun = (TF1*) arrFun->At(i);
      AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
      // clone temporare histogram, otherwise it will not be streamed to file!
      fPostPIDWdthCorrArr->Add((TH1*)fun->GetHistogram()->Clone(key.Data()));
      histo = (TH1*) fPostPIDWdthCorrArr->At(i);
      histo->GetListOfFunctions()->AddAt(fun, 0);
    }
    // if TH1s are stored in the array
    if(bHisto){
      TH1 *fun = (TH1*) arrFun->At(i);
      AliDielectronHistos::StoreVariables(fun, valType);
      // clone temporare histogram, otherwise it will not be streamed to file!
      fPostPIDWdthCorrArr->Add((TH1*)fun->Clone(key.Data()));
      histo = (TH1*) fPostPIDWdthCorrArr->At(i);
    }
    if(histo){
      // check for corrections and add their variables to the fill map
      printf("POST TPC PID CORRECTION added for %s:  ",key.Data());
      switch(histo->GetDimension()) {
      case 3: printf(" %s, ",histo->GetZaxis()->GetName());
      case 2: printf(" %s, ",histo->GetYaxis()->GetName());
      case 1: printf(" %s ",histo->GetXaxis()->GetName());
      }
      printf("\n");
    }
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}

//______________________________________________
void AliDielectron::SetCentroidCorrFunction(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorr = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDCntrdCorr)  {
    fPostPIDCntrdCorr->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST TPC PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDCntrdCorr->GetZaxis()->GetName());
    case 2: printf(" %s, ",fPostPIDCntrdCorr->GetYaxis()->GetName());
    case 1: printf(" %s ",fPostPIDCntrdCorr->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetCentroidCorrFunction(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorr = (TH1*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorr)  {
    printf("POST TPC PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDCntrdCorr->GetZaxis()->GetName());
    case 2: printf(" %s, ",fPostPIDCntrdCorr->GetYaxis()->GetName());
    case 1: printf(" %s ",fPostPIDCntrdCorr->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunction(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorr = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDWdthCorr)  {
    fPostPIDWdthCorr->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST TPC PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDWdthCorr->GetZaxis()->GetName());
    case 2: printf(" %s, ",fPostPIDWdthCorr->GetYaxis()->GetName());
    case 1: printf(" %s ",fPostPIDWdthCorr->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunction(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorr = (TH1*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorr)  {
    printf("POST TPC PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorr->GetDimension()) {
    case 3: printf(" %s, ",fPostPIDWdthCorr->GetZaxis()->GetName());
    case 2: printf(" %s, ",fPostPIDWdthCorr->GetYaxis()->GetName());
    case 1: printf(" %s ",fPostPIDWdthCorr->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}

/// @TODO: implement smart way to use generic functions for all detectors.
//______________________________________________
void AliDielectron::SetCentroidCorrFunctionITS(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorrITS = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDCntrdCorrITS)  {
    fPostPIDCntrdCorrITS->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST ITS PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetCentroidCorrFunctionITS(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorrITS = (TH1*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorrITS)  {
    printf("POST ITS PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunctionITS(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorrITS = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDWdthCorrITS)  {
    fPostPIDWdthCorrITS->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST ITS PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunctionITS(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorrITS = (TH1*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorrITS)  {
    printf("POST ITS PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetCentroidCorrFunctionTOF(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorrTOF = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDCntrdCorrTOF)  {
    fPostPIDCntrdCorrTOF->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST TOF PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTOF->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrTOF->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTOF->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrTOF->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetCentroidCorrFunctionTOF(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d",varx,vary,varz);
  fPostPIDCntrdCorrTOF = (TH1*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorrTOF)  {
    printf("POST TOF PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTOF->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrTOF->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTOF->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrTOF->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunctionTOF(TF1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun->GetHistogram(), valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorrTOF = (TH1*)fun->GetHistogram()->Clone(key.Data());
  if(fPostPIDWdthCorrTOF)  {
    fPostPIDWdthCorrTOF->GetListOfFunctions()->AddAt(fun,0);
    // check for corrections and add their variables to the fill map
    printf("POST TOF PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTOF->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrTOF->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTOF->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrTOF->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunctionTOF(TH1 *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d",varx,vary,varz);
  fPostPIDWdthCorrTOF = (TH1*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorrTOF)  {
    printf("POST TOF PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTOF->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrTOF->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTOF->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrTOF->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE);
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}


//______________________________________________
void AliDielectron::SetCentroidCorrFunctionPU(UInt_t detID, UInt_t parID, THnBase *fun, UInt_t var0, UInt_t var1, UInt_t var2, UInt_t var3, UInt_t var4)
{
  UInt_t valType[20] = {0};
  valType[0]=var0;     valType[1]=var1;     valType[2]=var2;     valType[3]=var3;     valType[4]=var4;

  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("cntrd%d%d%d%d%d_%d%d",var0,var1,var2,var3,var4,detID,parID);

  fPostPIDCntrdCorrPU[detID][parID] = (THnBase*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDCntrdCorrPU[detID][parID])  {
    printf("detID = %u , parID = %u, POST PID CORRECTION in PU added for centroids:  ",detID,parID);
    switch(fPostPIDCntrdCorrPU[detID][parID]->GetNdimensions()) {
      case 5: printf(" %s, ",fPostPIDCntrdCorrPU[detID][parID]->GetAxis(4)->GetName());
      case 4: printf(" %s, ",fPostPIDCntrdCorrPU[detID][parID]->GetAxis(3)->GetName());
      case 3: printf(" %s, ",fPostPIDCntrdCorrPU[detID][parID]->GetAxis(2)->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrPU[detID][parID]->GetAxis(1)->GetName());
      case 1: printf(" %s " ,fPostPIDCntrdCorrPU[detID][parID]->GetAxis(0)->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(var0, kTRUE);
    fUsedVars->SetBitNumber(var1, kTRUE);
    fUsedVars->SetBitNumber(var2, kTRUE);
    fUsedVars->SetBitNumber(var3, kTRUE);
    fUsedVars->SetBitNumber(var4, kTRUE);
  }
}
//______________________________________________
void AliDielectron::SetWidthCorrFunctionPU(UInt_t detID, UInt_t parID, THnBase *fun, UInt_t var0, UInt_t var1, UInt_t var2, UInt_t var3, UInt_t var4)
{
  UInt_t valType[20] = {0};
  valType[0]=var0;     valType[1]=var1;     valType[2]=var2;     valType[3]=var3;     valType[4]=var4;

  AliDielectronHistos::StoreVariables(fun, valType);
  // clone temporare histogram, otherwise it will not be streamed to file!
  TString key = Form("wdth%d%d%d%d%d_%d%d",var0,var1,var2,var3,var4,detID,parID);

  fPostPIDWdthCorrPU[detID][parID] = (THnBase*)fun->Clone(key.Data());
  // check for corrections and add their variables to the fill map
  if(fPostPIDWdthCorrPU[detID][parID])  {
    printf("detID = %u , parID = %u, POST PID CORRECTION IN PU added for widths:  ",detID,parID);
    switch(fPostPIDWdthCorrPU[detID][parID]->GetNdimensions()) {
      case 5: printf(" %s, ",fPostPIDWdthCorrPU[detID][parID]->GetAxis(4)->GetName());
      case 4: printf(" %s, ",fPostPIDWdthCorrPU[detID][parID]->GetAxis(3)->GetName());
      case 3: printf(" %s, ",fPostPIDWdthCorrPU[detID][parID]->GetAxis(2)->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrPU[detID][parID]->GetAxis(1)->GetName());
      case 1: printf(" %s " ,fPostPIDWdthCorrPU[detID][parID]->GetAxis(0)->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(var0, kTRUE);
    fUsedVars->SetBitNumber(var1, kTRUE);
    fUsedVars->SetBitNumber(var2, kTRUE);
    fUsedVars->SetBitNumber(var3, kTRUE);
    fUsedVars->SetBitNumber(var4, kTRUE);
  }
}
//______________________________________________
TObject* AliDielectron::InitEffMap(TString filename, TString generatedname, TString foundname)
{
  // init an efficiency object for on-the-fly correction calculations
  if(filename.Contains("alien://") && !gGrid) TGrid::Connect("alien://",0,0,"t");

  TFile* file=TFile::Open(filename.Data());
  if(!file) return 0x0;
	else printf("[I]  AliDielectron::InitEffMap efficiency maps file %s loaded! \n",filename.Data());

  // NOTE: the spline must have the 'variable name' stored in its fHistogram
  TSpline3 *hEff = (TSpline3*) file->Get("hEfficiency");
  //if(hEff) printf("we use a TSpline!!!!!!!!!!! \n");
  if(hEff){
    printf("[II] AliDielectron::InitEffMap TSpline3 loaded! \n");
    return (hEff->Clone("effMap"));
  }
  THnBase *hGen = (THnBase*) file->Get(generatedname.Data());
  THnBase *hFnd = (THnBase*) file->Get(foundname.Data());
  if(!hFnd || !hGen) return 0x0;
  printf("[II] AliDielectron::InitEffMap THnBase generated: %s and found: %s loaded! \n",generatedname.Data(),foundname.Data());
  hFnd->Divide(hGen);
  return (hFnd->Clone("effMap"));
}

//________________________________________________________________
void AliDielectron::FillHistogramsFromPairArray(Bool_t pairInfoOnly/*=kFALSE*/)
{
  //
  // Fill Histogram information for tracks and pairs
  //

  TString  className,className2;
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::SetLegEffMap(fLegEffMap);
  AliDielectronVarManager::SetPairEffMap(fPairEffMap);

  //Fill event information
  if(!pairInfoOnly) {
    if(fHistos->GetHistogramList()->FindObject("Event")) {
      fHistos->FillClass("Event", AliDielectronVarManager::kNMaxValues, AliDielectronVarManager::GetData());
    }
  }

  UInt_t selectedMask=(1<<fPairFilter.GetCuts()->GetEntries())-1;

  //Fill Pair information, separately for all pair candidate arrays and the legs
  TObjArray arrLegs(100);
  for (Int_t i=0; i<10; ++i){ // ROT pairs??
    Int_t npairs=PairArray(i)->GetEntriesFast();
    if(npairs<1) continue;

    className.Form("Pair_%s",fgkPairClassNames[i]);
    className2.Form("Track_Legs_%s",fgkPairClassNames[i]);
    Bool_t pairClass=fHistos->GetHistogramList()->FindObject(className.Data())!=0x0;
    Bool_t legClass=fHistos->GetHistogramList()->FindObject(className2.Data())!=0x0;

    //    if (!pairClass&&!legClass) continue;
    for (Int_t ipair=0; ipair<npairs; ++ipair){
      AliDielectronPair *pair=static_cast<AliDielectronPair*>(PairArray(i)->UncheckedAt(ipair));

      // apply cuts
      UInt_t cutMask=fPairFilter.IsSelected(pair);

      // cut qa
      if(i==kEv1PM && fCutQA) {
        fQAmonitor->FillAll(pair);
        fQAmonitor->Fill(cutMask,pair);
      }

      //CF manager for the pair (TODO: check steps and if they are properly filled)
      //      if (fCfManagerPair) fCfManagerPair->Fill(cutMask,pair);

      //apply cut
      if (cutMask!=selectedMask) continue;

      //histogram array for the pair
      if (fHistoArray) fHistoArray->Fill(i,pair);

      // fill map
      AliDielectronVarManager::SetFillMap(fUsedVars);

      //fill pair information
      if (pairClass){
        AliDielectronVarManager::Fill(pair, values);
        fHistos->FillClass(className, AliDielectronVarManager::kNMaxValues, values);
      }

      //fill leg information, don't fill the information twice
      if (legClass){
        AliVParticle *d1=pair->GetFirstDaughterP();
        AliVParticle *d2=pair->GetSecondDaughterP();
        if (!arrLegs.FindObject(d1)){
          AliDielectronVarManager::Fill(d1, values);
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
          arrLegs.Add(d1);
        }
        if (!arrLegs.FindObject(d2)){
          AliDielectronVarManager::Fill(d2, values);
          fHistos->FillClass(className2, AliDielectronVarManager::kNMaxValues, values);
          arrLegs.Add(d2);
        }
      }
    }
    if (legClass) arrLegs.Clear();
  }

}

//______________________________________________
void AliDielectron::FinishEvtVsTrkHistoClass()
{
  if(fEvtVsTrkHist) fEvtVsTrkHist->CalculateMatchingEfficiency();
}


void AliDielectron::SetRotatedTrackWeightMap(TString filename, TString histoname){
  TFile* file = TFile::Open(filename.Data(), "READ");
  printf("%p\n", file);
  if (file == 0x0){
    gSystem->Exec(Form("alien_cp alien://%s .",filename.Data()));
    printf("Copy rotated track map from Alien\n");
    TObjArray *arrNames=filename.Tokenize("/");
    TString name = arrNames->Last()->GetName();
    file = TFile::Open(Form("%s", name.Data()));
  }
  else {
    printf("Track Correction Map loaded\n");
  }
  if (file == nullptr){
    AliFatal(Form("Rotated-Track-Weighting file %s not found!", filename.Data()));
  }
  fRotateTrackCorrectionMap = *(dynamic_cast<TH3F*>(file->Get(histoname.Data())));
  if (&fRotateTrackCorrectionMap == nullptr){
    AliFatal(Form("Weighting histogram %s not found!", histoname.Data()));
  }
  fRotateTrackCorrectionMap.SetDirectory(0);
  file->Close();
}


//______________________________________________
Double_t AliDielectron::GetWeightFromRotation(AliKFParticle* part){
  if(!fUseAccMap){
    return 1;
  }
  else{
    int bin_pt  = fRotateTrackCorrectionMap.GetXaxis()->FindBin(part->GetPt() );
    const int bin_eta = fRotateTrackCorrectionMap.GetYaxis()->FindBin(part->GetEta());

    if(bin_pt < fRotWeight_minPtBin){
      if(fRotWeight_minPtBin <= 1) bin_pt = 1;
      else bin_pt = fRotWeight_minPtBin;
    }
    if(bin_pt > fRotWeight_maxPtBin){
      if(fRotWeight_maxPtBin >= fRotateTrackCorrectionMap.GetXaxis()->GetNbins()) bin_pt = fRotateTrackCorrectionMap.GetXaxis()->GetNbins();
      else bin_pt = fRotWeight_maxPtBin;
    }


    double phi = part->GetPhi();
    if (phi < 0) phi += TMath::TwoPi();
    const int bin_phi = fRotateTrackCorrectionMap.GetZaxis()->FindBin(phi);

    Double_t weight = fRotateTrackCorrectionMap.GetBinContent(bin_pt, bin_eta, bin_phi);

    Int_t i = 1;
    while(weight <= 0){
      weight = fRotateTrackCorrectionMap.GetBinContent(bin_pt-i, bin_eta, bin_phi);
      i++;
      if (bin_pt-i <= 1)
        break;
    }

    return weight;  

  }
}



void AliDielectron::SetRotatedPairWeightMap(TString filename, TString histoname){
  
  TFile* file = TFile::Open(filename.Data(), "READ");
  printf("%p\n", file);
  if (file == 0x0){
    gSystem->Exec(Form("alien_cp alien://%s .",filename.Data()));
    printf("Copy rotated pair map from Alien\n");
    TObjArray *arrNames=filename.Tokenize("/");
    TString name = arrNames->Last()->GetName();
    file = TFile::Open(Form("%s", name.Data()));
  }
  else {
    printf("Pair Correction Map loaded\n");
  }
  if (file == nullptr){
    AliFatal(Form("Rotated-Pair-Weighting file %s not found!", filename.Data()));
  }
  fRotatePairCorrectionMap = *(dynamic_cast<TH2F*>(file->Get(histoname.Data())));
  if (&fRotatePairCorrectionMap == nullptr){
    AliFatal(Form("Weighting histogram %s not found!", histoname.Data()));
  }
  fRotatePairCorrectionMap.SetDirectory(0);
  file->Close();
}


Double_t AliDielectron::GetWeightFromOpeningAngle(AliKFParticle* KFpos, AliKFParticle* KFneg){
  // if(!fUseAccMap){
  //  return 1;
  //}
  //else{
    static const double electron_mass = AliPID::ParticleMass(AliPID::kElectron);
    TLorentzVector LvecPos;
    LvecPos.SetPtEtaPhiM(KFpos->GetPt(), KFpos->GetEta(), KFpos->GetPhi(), electron_mass);
    TLorentzVector LvecNeg;
    LvecNeg.SetPtEtaPhiM(KFneg->GetPt(), KFneg->GetEta(), KFneg->GetPhi(), electron_mass);
    TLorentzVector LvecMother = LvecPos + LvecNeg;
    Double_t ptee         = LvecMother.Pt();
    Double_t openingAngle = LvecPos.Angle(LvecNeg.Vect());

    if(openingAngle > 0.5) return 1.;

    int bin_ptee  = fRotatePairCorrectionMap.GetYaxis()->FindBin(ptee);
    const int bin_opAng = fRotatePairCorrectionMap.GetXaxis()->FindBin(openingAngle);


    double weight = fRotatePairCorrectionMap.GetBinContent(bin_opAng, bin_ptee);
    Int_t i = 1;
    while(weight <= 0){
      weight = fRotatePairCorrectionMap.GetBinContent(bin_opAng, bin_ptee-i);
      i++;
      if (bin_ptee-i <= 1)
        return 1.;
    }
    return weight;
  //}
}
 
