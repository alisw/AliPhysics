#include "AliLog.h"

// analysis framework
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// ESD stuff
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

// GTU simulation
#include "AliTRDgtuParam.h"
#include "AliTRDgtuSim.h"

#include "AliAnalysisTaskTRDgtuSim.h"
#include "AliTRDonlineTrackMatching.h"

AliAnalysisTaskTRDgtuSim::AliAnalysisTaskTRDgtuSim(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputList(0x0),
  fHistStat(0x0),
  fHistDeltaA(0x0),
  fHistDeltaB(0x0),
  fHistDeltaC(0x0),
  fGtuSim(new AliTRDgtuSim()),
  fTrackletLabel(-2),
  fLabel(-10),
  fDeltaY(12),
  fDeltaAlpha(16),
  fLimitNoTracklets(kTRUE),
  fMaxNoTracklets(62)
{
  // ctor

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskTRDgtuSim::~AliAnalysisTaskTRDgtuSim()
{
  // dtor

  delete fGtuSim;
}

Bool_t AliAnalysisTaskTRDgtuSim::Notify()
{
  // Implemented Notify()

  return AliAnalysisTaskSE::Notify();
}

void AliAnalysisTaskTRDgtuSim::UserCreateOutputObjects()
{
  // create output objects

  OpenFile(1);
  fOutputList = new TList;
  fOutputList->SetOwner(kTRUE);

  fHistStat = new TH1F("stat", "stat", 4, .5, 4.5);
  fHistStat->GetXaxis()->SetBinLabel(1, "identical");
  fHistStat->GetXaxis()->SetBinLabel(2, "mismatch");
  fHistStat->GetXaxis()->SetBinLabel(3, "sim only");
  fHistStat->GetXaxis()->SetBinLabel(4, "raw only");
  fOutputList->Add(fHistStat);

  fHistDeltaA = new TH1F("deltaA", "delta A", 100, -100., 100.);
  fOutputList->Add(fHistDeltaA);

  fHistDeltaB = new TH1F("deltaB", "delta B", 100, -4500., 4500.);
  fOutputList->Add(fHistDeltaB);

  fHistDeltaC = new TH1F("deltaC", "delta C", 100, -100., 100.);
  fOutputList->Add(fHistDeltaC);

  PostData(1, fOutputList);
}

void AliAnalysisTaskTRDgtuSim::UserExec(Option_t * /* option */)
{
  // perform actual analysis

  AliESDEvent *esdEvent =
    dynamic_cast<AliESDEvent*>(this->InputEvent());

  if (!esdEvent)
    return;

  AliTRDgtuParam::SetDeltaY(fDeltaY);
  AliTRDgtuParam::SetDeltaAlpha(fDeltaAlpha);
  AliTRDgtuParam::SetLimitNoTracklets(fLimitNoTracklets);
  AliTRDgtuParam::SetMaxNoTracklets(fMaxNoTracklets);

  fGtuSim->RunGTU(0x0, esdEvent, fTrackletLabel, fLabel);

  AliTRDonlineTrackMatching trdMatch;
  trdMatch.ProcessEvent(esdEvent, kTRUE, fLabel);

  // compare to hardware tracks
  Check(fLabel, -3);

  PostData(1, fOutputList);
}


void AliAnalysisTaskTRDgtuSim::Terminate(const Option_t*)
{

}

void AliAnalysisTaskTRDgtuSim::Check(Int_t label, Int_t labelRef)
{
  // compare tracks with label label with tracks with label labelRef

  AliESDEvent *esdEvent =
    dynamic_cast<AliESDEvent*>(this->InputEvent());

  if (!esdEvent)
    return;

  TList tracksCheck;
  TList tracksRef;

  Int_t nTracks = esdEvent->GetNumberOfTrdTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    AliESDTrdTrack *trk = esdEvent->GetTrdTrack(iTrack);
    if (trk->GetLabel() == label)
      tracksCheck.Add(trk);
    if (trk->GetLabel() == labelRef)
      tracksRef.Add(trk);
  }

  TIter iterTrack(&tracksCheck);
  while (AliESDTrdTrack *trk = (AliESDTrdTrack*) iterTrack()) {
    TIter iterTrackRef(&tracksRef);
    Bool_t foundMatch = kFALSE;
    while (AliESDTrdTrack *trkRef = (AliESDTrdTrack*) iterTrackRef()) {
      // check for common tracklets
      Bool_t commonTracklet = kFALSE;
      Bool_t allEqual = kTRUE;
      for (Int_t iLayer = 0; iLayer < 6; ++iLayer) {
	if (trk->GetTracklet(iLayer)) {
	  if (trk->GetTracklet(iLayer) == trkRef->GetTracklet(iLayer))
	    commonTracklet = kTRUE;
	  else
	    allEqual = kFALSE;
	}
	else if (trkRef->GetTracklet(iLayer))
	  allEqual = kFALSE;
      }
      if (commonTracklet) {
	// tracks with a common tracklet should be identical
	if (allEqual) {
	  // identical track composition
	  fHistStat->Fill(1);
	  Int_t deltaA = trk->GetA() - trkRef->GetA();
	  fHistDeltaA->Fill(deltaA);
	  Int_t deltaB = trk->GetB() - trkRef->GetB();
	  fHistDeltaB->Fill(deltaB);
	  Int_t deltaC = trk->GetC() - trkRef->GetC();
	  fHistDeltaC->Fill(deltaC);
	}
	else {
	  // track with different tracklets
	  fHistStat->Fill(2);
	}
	tracksRef.Remove(trkRef);
	foundMatch = kTRUE;
	break;
      }
    }
    if (!foundMatch) {
      // unmatched sim track
      fHistStat->Fill(3);
    }
  }
  TIter iterTrackRef(&tracksRef);
  while (iterTrackRef()) {
    // unmatched raw track
    fHistStat->Fill(4);
  }
}
