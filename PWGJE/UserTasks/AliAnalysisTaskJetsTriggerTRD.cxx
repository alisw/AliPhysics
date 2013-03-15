// ROOT
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// analysis framework
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"

// MC stuff
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"

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

// jet tasks
#include "AliAnalysisTaskJetServices.h"
#include "AliAnalysisHelperJetTasks.h"

#include "AliAnalysisTaskJetsTriggerTRD.h"

const char * const histos[] = { "stat", "jetpt", "jetpt_emc", "jetpt_trd" };

AliAnalysisTaskJetsTriggerTRD::AliAnalysisTaskJetsTriggerTRD(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputList(),
  fHist(),
  fShortTaskId("jets_trg_trd"),
  fNoJetPtBins(40),
  fJetPtBinMax(400)
{
  // default ctor

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetsTriggerTRD::~AliAnalysisTaskJetsTriggerTRD()
{
  // dtor

}

void AliAnalysisTaskJetsTriggerTRD::UserCreateOutputObjects()
{
  // create user output objects

  // setup list
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  // setup histograms
  AddHistogram(HISTID(kHistStat), "event statistics",
  	       10, -.5, 9.5);
  AddHistogram(HISTID(kHistJetPt), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(HISTID(kHistJetPtEMC), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(HISTID(kHistJetPtTRD), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);

  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskJetsTriggerTRD::Notify()
{
  // actions to be taken upon notification about input file change

  return AliAnalysisTaskSE::Notify();
}

void AliAnalysisTaskJetsTriggerTRD::UserExec(Option_t * /* option */)
{
  // actual work

  // setup pointers to input data (null if unavailable)
  // mcEvent:  MC input
  // esdEvent: ESD input
  // outEvent: AOD output
  // aodEvent: AOD input if available, otherwise AOD output
  // AliMCEvent  *mcEvent   = this->MCEvent();
  AliESDEvent *esdEvent  = dynamic_cast<AliESDEvent*>(this->InputEvent()); // could also be AOD input
  AliAODEvent* outEvent  = this->AODEvent();
  AliAODEvent *aodEvent  = outEvent;
  if (dynamic_cast<AliAODEvent*>(this->InputEvent()))
    aodEvent = (AliAODEvent*) (this->InputEvent());

  if ((fDebug > 0) && esdEvent)
    printf("event: %s-%06i\n", CurrentFileName(), esdEvent->GetEventNumberInFile());

  // Int_t nTracksMC  = 0; // no. of MC tracks
  // Int_t nTracksESD = 0; // no. of global ESD tracks
  // Int_t nTracksGTU = 0; // no. of GTU tracks

  // Int_t nTracks[6][90]; // tracks above lower pt threshold, counted stack-wise
  // memset(nTracks, 0, sizeof(Int_t)*6*90);
  // Int_t nMax[6] = { 0 };

  TList partMC;
  TList partEsd;
  TList partGtu;

  // Float_t leadingJetPtMC = 0.; // leading jet energy from MC information
  Float_t leadingJetPtRec = 0.; // leading jet energy from AOD information

  // number of sampled events
  FillH1(kHistStat, 1);

  // event selection
  AliVEvent::EOfflineTriggerTypes physSel = (AliVEvent::EOfflineTriggerTypes) ((AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if (fDebug > 2)
    printf("trg: %8s %8s %8s %8s %8s %8s (%s)\n",
           (physSel & AliVEvent::kAnyINT)  ? "kAnyINT"  : " ",
           (physSel & AliVEvent::kCINT5)   ? "kCINT5" : " ",
           (physSel & AliVEvent::kINT7)    ? "kINT7" : " ",
           (physSel & AliVEvent::kINT8)    ? "kINT8" : " ",
           (physSel & AliVEvent::kMB)      ? "kMB"   : " ",
           (physSel & AliVEvent::kEMC7)    ? "kEMC7" : " ",
           esdEvent ? esdEvent->GetFiredTriggerClasses().Data() : "<no ESD event>"
           );
  // Int_t triggerMask = 0;
  if (!(physSel & (AliVEvent::kAnyINT)))
    return;

  // store no. of sampled events
  FillH1(kHistStat, 2);

  AliAODJet *leadJet = 0x0;
  AliAODJet *subleadJet = 0x0;

  if (aodEvent) {
    TClonesArray *jetArray = dynamic_cast<TClonesArray*> (aodEvent->FindListObject(fJetBranchName));
    if (jetArray) {
      Int_t nJets = jetArray->GetEntriesFast();
      for (Int_t iJet = 0; iJet < nJets; ++iJet) {
        AliAODJet *jet = (AliAODJet*) (*jetArray)[iJet];

        // check contributing tracks
        Int_t nJetTracks = jet->GetRefTracks()->GetEntriesFast();
        Int_t iLeadingTrack = -1;
        Float_t ptLeadingTrack = 0.;
        for (Int_t iTrack=0; iTrack < nJetTracks; ++iTrack) {
          AliAODTrack *track = (AliAODTrack*) jet->GetRefTracks()->At(iTrack);
          if (track->Pt() > ptLeadingTrack) {
            ptLeadingTrack = track->Pt();
            iLeadingTrack = iTrack;
          }
        }

        // check tracking flags for leading track
        // discard the jet if the leading track has no ITS contribution with a hit in any SPD layer
        AliAODTrack *leadingTrack = (AliAODTrack*) jet->GetRefTracks()->At(iLeadingTrack);
        if (!(leadingTrack->GetFlags() & AliVTrack::kITSrefit) ||
            !(leadingTrack->HasPointOnITSLayer(0) || leadingTrack->HasPointOnITSLayer(1)))
          continue;

        if (TMath::Abs(jet->Eta()) < 0.5) {
          if (TMath::Abs(jet->Pt()) > leadingJetPtRec) {
            leadingJetPtRec = TMath::Abs(jet->Pt());
            subleadJet = leadJet;
            leadJet = jet;
          }

  	  FillH1(kHistJetPt, jet->Pt());
          // if (esdEvent && esdEvent->IsTriggerClassFired("CEMC7WU-B-NOPF-ALL"))
  	  //   FillH1("jetpt_emc", jet->Pt());
          // // if (esdEvent && esdEvent->IsTriggerClassFired("CINT7WU-I-NOPF-ALL"))
          // if (triggerMask & (1 << kMBHJT))
  	  //   FillH1("jetpt_trd", jet->Pt());
        }
      }
    }
    else {
      printf("no jet array found\n");
    }
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskJetsTriggerTRD::Terminate(const Option_t * /* option */)
{
  // actions at task termination

}

// ----- histogram management -----
TH1* AliAnalysisTaskJetsTriggerTRD::AddHistogram(Hist_t hist, const char *hid, TString title,
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

TH2* AliAnalysisTaskJetsTriggerTRD::AddHistogram(Hist_t hist, const char *hid, TString title,
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

TH3* AliAnalysisTaskJetsTriggerTRD::AddHistogram(Hist_t hist, const char *hid, TString title,
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
