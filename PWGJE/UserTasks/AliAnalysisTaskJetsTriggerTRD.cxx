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
#include "AliVTrdTrack.h"

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

AliAnalysisTaskJetsTriggerTRD::AliAnalysisTaskJetsTriggerTRD(const char *name) :
  AliAnalysisTaskSE(name),
  fTriggerMask(0),
  fOutputList(),
  fHist(),
  fShortTaskId("jets_trg_trd"),
  fNoJetPtBins(80),
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
  TH1 *histStat = AddHistogram(ID(kHistStat), "event statistics;;counts",
                               kStatLast-1, .5, kStatLast-.5);
  histStat->GetXaxis()->SetBinLabel(ID(kStatSeen));
  histStat->GetXaxis()->SetBinLabel(ID(kStatUsed));
  histStat->GetXaxis()->SetBinLabel(ID(kStatMB));

  AddHistogram(ID(kHistNoJets), "number of jets;N^{jet};counts",
               100, -.5, 99.5);

  AddHistogram(ID(kHistTrackGTU), "GTU track p_{T};p_{T} (GeV/c);counts",
               100, 0., 25.);

  AddHistogram(ID(kHistJetPt), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistJetPtITS), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistJetPt3x3), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);

  AddHistogram(ID(kHistJetPtEMC), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistJetPtHJT), "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax);

  AddHistogram(ID(kHistJetPtNoTracks3), "number of tracks above 3 GeV;p_{T}^{jet};no. of tracks",
               fNoJetPtBins, 0., fJetPtBinMax,
               20, -.5, 19.5);

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
  Int_t nTracksGTU = 0; // no. of GTU tracks

  // Int_t nTracks[6][90]; // tracks above lower pt threshold, counted stack-wise
  // memset(nTracks, 0, sizeof(Int_t)*6*90);
  // Int_t nMax[6] = { 0 };

  TList partMC;
  TList partEsd;
  TList partGtu;

  // Float_t leadingJetPtMC = 0.; // leading jet energy from MC information
  Float_t leadingJetPtRec = 0.; // leading jet energy from AOD information

  // record number of sampled events and detect trigger contributions
  FillH1(kHistStat, kStatSeen);
  if (!DetectTriggers()) {
    AliError("Failed to detect the triggers");
    return;
  }

  if (!IsTrigger(kTrgInt))
    return;
  FillH1(kHistStat, kStatUsed);

  for (Int_t iTrack = 0; iTrack < nTracksGTU; ++iTrack) {
    AliVTrdTrack *trk = 0x0;
    FillH1(kHistTrackGTU, trk->Pt());
  }

  AliAODJet *leadJet = 0x0;
  AliAODJet *subleadJet = 0x0;

  if (aodEvent) {
    TClonesArray *jetArray = dynamic_cast<TClonesArray*> (aodEvent->FindListObject(fJetBranchName));
    if (jetArray) {
      Int_t nJets = jetArray->GetEntriesFast();
      FillH1(kHistNoJets, nJets);
      for (Int_t iJet = 0; iJet < nJets; ++iJet) {
        AliAODJet *jet = (AliAODJet*) (*jetArray)[iJet];

        // check contributing tracks
        Int_t nJetTracks = jet->GetRefTracks()->GetEntriesFast();
        Int_t nJetTracks3 = 0;
        Int_t iLeadingTrack = -1;
        Float_t ptLeadingTrack = 0.;
        for (Int_t iTrack=0; iTrack < nJetTracks; ++iTrack) {
          AliAODTrack *track = (AliAODTrack*) jet->GetRefTracks()->At(iTrack);
          // count tracks above 3 GeV/c
          if (track->Pt() > 3.)
            ++nJetTracks3;
          // find the leading track
          if (track->Pt() > ptLeadingTrack) {
            ptLeadingTrack = track->Pt();
            iLeadingTrack = iTrack;
          }
        }

        // retrieve the leading track
        AliAODTrack *leadingTrack = (AliAODTrack*) jet->GetRefTracks()->At(iLeadingTrack);

        if (TMath::Abs(jet->Eta()) < 0.5) {

          // find leading jet
          if (TMath::Abs(jet->Pt()) > leadingJetPtRec) {
            leadingJetPtRec = TMath::Abs(jet->Pt());
            subleadJet = leadJet;
            leadJet = jet;
          }

  	  FillH1(kHistJetPt, jet->Pt());

          // discard the jet if the leading track has no ITS contribution with a hit in any SPD layer
          if (leadingTrack &&
              (leadingTrack->GetFlags() & AliVTrack::kITSrefit) &&
              (leadingTrack->HasPointOnITSLayer(0) || leadingTrack->HasPointOnITSLayer(1)))
            FillH1(kHistJetPtITS, jet->Pt());

          // discard the jet if it does not have 3 tracks above 3 GeV/c
          FillH2(kHistJetPtNoTracks3, jet->Pt(), nJetTracks3);
          if (nJetTracks3 > 2)
            FillH1(kHistJetPt3x3, jet->Pt());

          if (IsTrigger(kTrgEMC))
            FillH1(kHistJetPtEMC, jet->Pt());
          if (IsTrigger(kTrgHJT))
            FillH1(kHistJetPtHJT, jet->Pt());
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

Bool_t AliAnalysisTaskJetsTriggerTRD::DetectTriggers()
{
  fTriggerMask = 0;

  AliVEvent::EOfflineTriggerTypes physSel = (AliVEvent::EOfflineTriggerTypes) ((AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  TString trgClasses = InputEvent()->GetFiredTriggerClasses();

  if (fDebug > 2)
    printf("trg: %8s %8s %8s %8s %8s %8s %8s (%s)\n",
           (physSel & AliVEvent::kAnyINT)  ? "kAnyINT"  : " ",
           (physSel & AliVEvent::kCINT5)   ? "kCINT5" : " ",
           (physSel & AliVEvent::kINT7)    ? "kINT7" : " ",
           (physSel & AliVEvent::kINT8)    ? "kINT8" : " ",
           (physSel & AliVEvent::kMB)      ? "kMB"   : " ",
           (physSel & AliVEvent::kEMC7)    ? "kEMC7" : " ",
           (physSel & AliVEvent::kTRD)     ? "kTRD" : " ",
           trgClasses.Data()
           );

  // physics selection
  if (physSel & (AliVEvent::kAnyINT))
    MarkTrigger(kTrgInt);

  // trigger classes
  if (trgClasses.Contains("CINT7WUHJT") || trgClasses.Contains("CINT8WUHJT"))
    MarkTrigger(kTrgHJT);
  if (trgClasses.Contains("CEMC7WU-") || trgClasses.Contains("CEMC8WU-"))
    MarkTrigger(kTrgEMC);

  return kTRUE;
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
