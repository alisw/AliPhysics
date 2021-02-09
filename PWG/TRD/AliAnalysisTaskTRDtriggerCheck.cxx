// ROOT
#include "TFile.h"
#include "TList.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom.h"

// analysis framework
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrdTrack.h"

// MC stuff
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliTrackReference.h"

// ESD stuff
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDTrdTrack.h"
#include "AliESDTrdTracklet.h"
#include "AliESDTrdTrigger.h"

// AOD stuff
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODTrack.h"

#include "AliTRDTriggerAnalysis.h"

#include "AliAnalysisTaskTRDtriggerCheck.h"

AliAnalysisTaskTRDtriggerCheck::AliAnalysisTaskTRDtriggerCheck(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputList(),
  fHist(),
  fShortTaskId("trd_trg_check")
{
  // default ctor

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskTRDtriggerCheck::~AliAnalysisTaskTRDtriggerCheck()
{
  // dtor

}

void AliAnalysisTaskTRDtriggerCheck::UserCreateOutputObjects()
{
  // create user output objects

  // setup list
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  TH2 *h = AddHistogram(ID(kHistTrgStat),
			"trigger contributions",
			kTrgLast, -.5, kTrgLast - .5,
			kLastComb - 1, .5, kLastComb - .5);
  h->GetXaxis()->SetBinLabel(kTrgHCO + 1, "HCO");
  h->GetXaxis()->SetBinLabel(kTrgHJT + 1, "HJT");
  h->GetXaxis()->SetBinLabel(kTrgHSE + 1, "HSE");
  h->GetXaxis()->SetBinLabel(kTrgHQU + 1, "HQU");
  h->GetXaxis()->SetBinLabel(kTrgHEE + 1, "HEE");
  h->GetYaxis()->SetBinLabel(ID(kNothing));
  h->GetYaxis()->SetBinLabel(ID(kTriggered));
  h->GetYaxis()->SetBinLabel(ID(kFired));
  h->GetYaxis()->SetBinLabel(ID(kCondition));
  h->GetYaxis()->SetBinLabel(ID(kTriggeredFired));
  h->GetYaxis()->SetBinLabel(ID(kTriggeredCondition));
  h->GetYaxis()->SetBinLabel(ID(kFiredCondition));
  h->GetYaxis()->SetBinLabel(ID(kTriggeredFiredCondition));

  h = AddHistogram(ID(kHistTrgStatSec),
		   "SM header flag per sector;;global stack",
		   kTrgLast + 1, -.5, kTrgLast + .5,
		   18, -.5, 17.5);
  h->GetXaxis()->SetBinLabel(kTrgHCO + 1, "HCO");
  h->GetXaxis()->SetBinLabel(kTrgHJT + 1, "HJT");
  h->GetXaxis()->SetBinLabel(kTrgHSE + 1, "HSE");
  h->GetXaxis()->SetBinLabel(kTrgHQU + 1, "HQU");
  h->GetXaxis()->SetBinLabel(kTrgHEE + 1, "HEE");
  h->GetXaxis()->SetBinLabel(kTrgLast + 1, "TO");

  h = AddHistogram(ID(kHistTrgStatStackCond),
		   "condition per stack;;global stack",
		   kTrgLast, -.5, kTrgLast - .5,
		   90, -.5, 89.5);
  h->GetXaxis()->SetBinLabel(kTrgHCO + 1, "HCO");
  h->GetXaxis()->SetBinLabel(kTrgHJT + 1, "HJT");
  h->GetXaxis()->SetBinLabel(kTrgHSE + 1, "HSE");
  h->GetXaxis()->SetBinLabel(kTrgHQU + 1, "HQU");
  h->GetXaxis()->SetBinLabel(kTrgHEE + 1, "HEE");

  h = AddHistogram(ID(kHistTrgStatStackCondNotFired),
		   "(condition + !fired) per stack;;global stack",
		   kTrgLast, -.5, kTrgLast - .5,
		   90, -.5, 89.5);
  h->GetXaxis()->SetBinLabel(kTrgHCO + 1, "HCO");
  h->GetXaxis()->SetBinLabel(kTrgHJT + 1, "HJT");
  h->GetXaxis()->SetBinLabel(kTrgHSE + 1, "HSE");
  h->GetXaxis()->SetBinLabel(kTrgHQU + 1, "HQU");
  h->GetXaxis()->SetBinLabel(kTrgHEE + 1, "HEE");

  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskTRDtriggerCheck::Notify()
{
  // actions to be taken upon notification about input file change

  return AliAnalysisTaskSE::Notify();
}

void AliAnalysisTaskTRDtriggerCheck::UserExec(Option_t * /* option */)
{
  // actual work

  // setup pointers to input data (null if unavailable)
  // mcEvent:  MC input
  // esdEvent: ESD input
  // outEvent: AOD output
  // aodEvent: AOD input if available, otherwise AOD output
  // AliMCEvent  *mcEvent   = this->MCEvent();
  AliESDEvent *esdEvent  = dynamic_cast<AliESDEvent*>(this->InputEvent()); // could also be AOD input
  // AliAODEvent* outEvent  = this->AODEvent();
  // AliAODEvent *aodEvent  = outEvent;
  // if (dynamic_cast<AliAODEvent*>(this->InputEvent()))
  //   aodEvent = (AliAODEvent*) (this->InputEvent());

  if ((fDebug > 0) && esdEvent)
    printf("event: %s-%06i\n", CurrentFileName(), esdEvent->GetEventNumberInFile());

//  if (!InputEvent()->GetFiredTriggerClasses().Contains("WU")) {
//    if (fDebug > 0)
//      printf("no TRD, returning\n");
//    return;
//  }

  // reproduce hardware decision
  AliTRDTriggerAnalysis trgAnalysis;
  // trgAnalysis.SetVerbosity(1);
  trgAnalysis.SetRequireInTime(kTRUE);
  trgAnalysis.SetRequireMatchElectron(kFALSE);
  trgAnalysis.CalcTriggers(InputEvent());

  // now what we rather want
  AliTRDTriggerAnalysis trgAnalysisGoal;
  trgAnalysisGoal.SetRequireInTime(kFALSE);
  trgAnalysisGoal.SetRequireMatchElectron(kTRUE);
  trgAnalysisGoal.CalcTriggers(InputEvent());

  Int_t nTracks = InputEvent()->GetNumberOfTrdTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    AliVTrdTrack *trk = InputEvent()->GetTrdTrack(iTrack);
    if (fDebug > 0)
      printf("trk %2i: %02i_%i - pt = %6.2f, PID = %3i, %10s, %11s\n",
	     iTrack, trk->GetSector(), trk->GetStack(),
	     trk->Pt(), trk->GetPID(),
	     trk->GetTrackInTime() ? "in-time" : "not in-time",
	     trk->GetTrackMatch() ? "matched" : "not matched");
  }

  const Int_t trgBits[] = { 0, 4, 5, 6, 1, 7 };
  for (Int_t iSector = 0; iSector < 18; ++iSector)
    for (Int_t iTrg = 0; iTrg <= kTrgLast; ++iTrg)
      if (trgAnalysis.CheckTrgFlags(trgBits[iTrg], iSector))
	FillH2(kHistTrgStatSec, iTrg, iSector);

  for (Int_t iTrg = 0; iTrg < kTrgLast; ++iTrg) {
    Bool_t trgd  = trgAnalysis.HasTriggered((AliTRDTriggerAnalysis::TRDTrigger_t) iTrg);
    Bool_t fired = trgAnalysis.HasFired((AliTRDTriggerAnalysis::TRDTrigger_t) iTrg);
    Bool_t cond  = trgAnalysis.CheckCondition((AliTRDTriggerAnalysis::TRDTrigger_t) iTrg);
    if (fDebug > 0)
      printf("trigger %i: %10s %10s %10s\n", iTrg,
	     trgd ? "triggered" : "",
	     fired ? "fired" : "",
	     cond ? "condition" : "");

    if (trgd && fired && cond)
      FillH2(kHistTrgStat, iTrg, kTriggeredFiredCondition);
    else if (trgd && fired && !cond)
      FillH2(kHistTrgStat, iTrg, kTriggeredFired);
    else if (trgd && !fired && cond)
      FillH2(kHistTrgStat, iTrg, kTriggeredCondition);
    else if (trgd && !fired && !cond)
      FillH2(kHistTrgStat, iTrg, kTriggered);
    else if (!trgd && fired && cond)
      FillH2(kHistTrgStat, iTrg, kFiredCondition);
    else if (!trgd && fired && !cond)
      FillH2(kHistTrgStat, iTrg, kFired);
    else if (!trgd && !fired && cond)
      FillH2(kHistTrgStat, iTrg, kCondition);
    else if (!trgd && !fired && !cond)
      FillH2(kHistTrgStat, iTrg, kNothing);

    for (Int_t iStack = 0; iStack < 90; ++iStack) {
      if (trgAnalysis.CheckCondition((AliTRDTriggerAnalysis::TRDTrigger_t) iTrg, iStack)) {
	FillH2(kHistTrgStatStackCond, iTrg, iStack);
	if (!fired)
	  FillH2(kHistTrgStatStackCondNotFired, iTrg, iStack);
      }
    }
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskTRDtriggerCheck::Terminate(const Option_t * /* option */)
{
  // actions at task termination
}



// ----- histogram management -----
TH1* AliAnalysisTaskTRDtriggerCheck::AddHistogram(Hist_t hist, const char *hid, TString title,
						 Int_t xbins, Float_t xmin, Float_t xmax,
						 Int_t binType)
{
  TString hName;
  hName.Form("%s_%s", fShortTaskId, hid);
  hName.ToLower();
  TH1 *h = 0x0;
  if (binType == 0)
    h = new TH1I(hName.Data(), title,
                 xbins, xmin, xmax);
  else
    h = new TH1F(hName.Data(), title,
                 xbins, xmin, xmax);
  GetHistogram(hist) = h;
  fOutputList->Add(h);
  return h;
}

TH2* AliAnalysisTaskTRDtriggerCheck::AddHistogram(Hist_t hist, const char *hid, TString title,
						 Int_t xbins, Float_t xmin, Float_t xmax,
						 Int_t ybins, Float_t ymin, Float_t ymax,
						 Int_t binType)
{
  TString hName;
  hName.Form("%s_%s", fShortTaskId, hid);
  hName.ToLower();
  TH1 *h = 0x0;
  if (binType == 0)
    h = GetHistogram(hist) = new TH2I(hName.Data(), title,
                                     xbins, xmin, xmax,
                                     ybins, ymin, ymax);
  else
    h = GetHistogram(hist) = new TH2F(hName.Data(), title,
                                     xbins, xmin, xmax,
                                     ybins, ymin, ymax);
  fOutputList->Add(h);
  return (TH2*) h;
}

TH3* AliAnalysisTaskTRDtriggerCheck::AddHistogram(Hist_t hist, const char *hid, TString title,
						 Int_t xbins, Float_t xmin, Float_t xmax,
						 Int_t ybins, Float_t ymin, Float_t ymax,
						 Int_t zbins, Float_t zmin, Float_t zmax,
						 Int_t binType)
{
  TString hName;
  hName.Form("%s_%s", fShortTaskId, hid);
  hName.ToLower();
  TH1 *h = 0x0;
  if (binType == 0)
    h = GetHistogram(hist) = new TH3I(hName.Data(), title,
                                     xbins, xmin, xmax,
                                     ybins, ymin, ymax,
                                     zbins, zmin, zmax);
  else
    h = GetHistogram(hist) = new TH3F(hName.Data(), title,
                                     xbins, xmin, xmax,
                                     ybins, ymin, ymax,
                                     zbins, zmin, zmax);
  fOutputList->Add(h);
  return (TH3*) h;
}
