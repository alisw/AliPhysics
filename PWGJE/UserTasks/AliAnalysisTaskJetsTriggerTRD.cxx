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
  fJetPtBinMax(400),
  fXsection(0.),
  fAvgTrials(0.),
  fPtHard(0.)
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
  histStat->GetXaxis()->SetBinLabel(ID(kStatTrg));
  histStat->GetXaxis()->SetBinLabel(ID(kStatUsed));
  histStat->GetXaxis()->SetBinLabel(ID(kStatEvCuts));

  AddHistogram(ID(kHistJetPtMC), "leading jet spectrum (MC, |#eta| < 0.5);p_{T} (GeV/c);counts",
	       fNoJetPtBins, 0., fJetPtBinMax);

  AddHistogram(ID(kHistNoJets), "number of jets;N^{jet};counts",
               400, -.5, 399.5);

  AddHistogram(ID(kHistTrackGTU), "GTU track p_{T};p_{T} (GeV/c);counts",
               100, 0., 25.);

  AddHistogram(ID(kHistNPtMin),
	       "rejection;p_{T}^{min};N_{trk};trigger",
  	       100, 0., 10.,
	       20, 0., 20.,
	       kTrgLast - 1, .5, kTrgLast - .5);

  AddHistogram(ID(kHistLeadJetPt),
	       "leading jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       kTrgLast - 1, .5, kTrgLast - .5);
  AddHistogram(ID(kHistJetPt),
	       "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);trigger",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       kTrgLast - 1, .5, kTrgLast - .5);
  AddHistogram(ID(kHistJetPtITS),
	       "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);trigger",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       kTrgLast - 1, .5, kTrgLast - .5);
  AddHistogram(ID(kHistJetPt3x3),
	       "jet spectrum (|#eta| < 0.5);p_{T} (GeV/c);trigger",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       kTrgLast - 1, .5, kTrgLast - .5);

  for (Int_t iHist = kHistLeadJetPt; iHist <= kHistJetPt3x3; ++iHist) {
    TH1 *h = GetHistogram(Hist_t (iHist));
    h->GetYaxis()->SetBinLabel(ID(kTrgMinBias));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt7));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt8));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC7));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC8));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt7WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt8WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC7WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC8WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMCEJE));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMCEGA));
  }

  AddHistogram(ID(kHistJetPtNoTracks3),
	       "number of tracks above 3 GeV;p_{T}^{jet};no. of tracks",
               fNoJetPtBins, 0., fJetPtBinMax,
               40, -.5, 39.5);

  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskJetsTriggerTRD::Notify()
{
  // actions to be taken upon notification about input file change

  // ??? check ???

  fXsection = 0.;
  fAvgTrials = 1.;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();

  if (tree) {
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }

    Float_t nEntries = (Float_t) tree->GetTree()->GetEntries();
    if (!AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(), fXsection, fAvgTrials)) {
      AliError("retrieval of cross section failed");
      // ??? what to set as cross section?
    }

    if (nEntries > 0.)
      fAvgTrials /= nEntries;
  }

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
  AliMCEvent  *mcEvent   = this->MCEvent();
  AliESDEvent *esdEvent  = dynamic_cast<AliESDEvent*>(this->InputEvent()); // could also be AOD input
  AliAODEvent* outEvent  = this->AODEvent();
  AliAODEvent *aodEvent  = outEvent;
  if (dynamic_cast<AliAODEvent*>(this->InputEvent()))
    aodEvent = (AliAODEvent*) (this->InputEvent());

  if ((fDebug > 0) && esdEvent)
    printf("event: %s-%06i\n", CurrentFileName(), esdEvent->GetEventNumberInFile());

  // Int_t nTracksMC  = mcEvent ? mcEvent->GetNumberOfTracks() : 0; // no. of MC tracks
  // Int_t nTracks    = InputEvent()->GetNumberOfTracks(); // no. of global tracks
  Int_t nTracksTRD = InputEvent()->GetNumberOfTrdTracks(); // no. of GTU tracks

  TList partMC;
  TList partEsd;
  TList partGtu;

  Float_t leadingJetPtMC = 0.; // leading jet energy from MC information
  Float_t leadingJetPtRec = 0.; // leading jet energy from AOD information

  // record number of sampled events and detect trigger contributions
  FillH1(kHistStat, kStatSeen);
  if (!DetectTriggers()) {
    AliError("Failed to detect the triggers");
    return;
  }

  // only continue for events from interesting triggers
  if (fTriggerMask == 0)
    return;
  FillH1(kHistStat, kStatTrg);

  // no further technical requirements for the event at the moment
  FillH1(kHistStat, kStatUsed);

  // apply event cuts
  const AliVVertex *vtx = InputEvent()->GetPrimaryVertex();
  if (!vtx ||
      (vtx->GetNContributors() < 3.) ||
      (vtx->GetZ() > 10.))
    return;
  FillH1(kHistStat, kStatEvCuts);

  // extract MC information
  if (mcEvent) {
    // check for PYTHIA event header
    AliGenPythiaEventHeader *pythiaHeader =
      dynamic_cast<AliGenPythiaEventHeader*> (mcEvent->GenEventHeader());
    if (!pythiaHeader) {
      AliWarning("MC event without PYTHIA event header!\n");
    }
    else {
      fPtHard  = pythiaHeader->GetPtHard();
      // Int_t nTrials = pythiaHeader->Trials();

      // loop over jets from PYTHIA
      for (Int_t iJet = 0; iJet < pythiaHeader->NTriggerJets(); iJet++) {
        Float_t p[4];
        pythiaHeader->TriggerJet(iJet, p);
	TLorentzVector pJet(p);
        Float_t pt  = pJet.Pt();
	// only consider jets with |eta| < 0.5
	Float_t eta = pJet.Eta();
	if (TMath::Abs(eta) > 0.5)
	  continue;
        if (pt > leadingJetPtMC)
          leadingJetPtMC = pt;
      }
      // fill histogram for leading jet pt spectrum
      FillH1(kHistJetPtMC, leadingJetPtMC, fXsection);
    }
  }

  // loop over GTU tracks
  for (Int_t iTrack = 0; iTrack < nTracksTRD; ++iTrack) {
    AliVTrdTrack *trk = InputEvent()->GetTrdTrack(iTrack);
    FillH1(kHistTrackGTU, TMath::Abs(trk->Pt()));
    partGtu.Add(trk);
  }
  partGtu.Sort(kSortAscending);

  Int_t nTracksPerStack[90] = { 0 };
  Int_t nTracksPerStackMax = 0;

  TIter nextPartGtu(&partGtu);
  while (AliVTrdTrack *trdTrack = (AliVTrdTrack*) nextPartGtu()) {
    // count no. of tracks in stack,
    // check whether this number was reached before,
    // if not store pt^min(n),
    // i.e. pt of current track because of sorting

    Int_t sec    = trdTrack->GetSector();
    Int_t stack  = trdTrack->GetStack();

    if ((sec > -1) && (sec < 18) &&
        (stack > -1) && (stack < 5)) {
      ++nTracksPerStack[5*sec + stack];
      if (nTracksPerStack[5*sec + stack] > nTracksPerStackMax) {
        ++nTracksPerStackMax;

	for (Int_t iTrigger = 1; iTrigger < kTrgLast; ++iTrigger)
	  if (IsTrigger(Trigger_t (iTrigger)))
	      FillH3(kHistNPtMin,
		     TMath::Abs(trdTrack->Pt()), nTracksPerStackMax, iTrigger);
      }
    }
    else {
      AliError(Form("Invalid sector or stack: %i %i",
                    sec, stack));
    }
    
  }

  // loop over jets from AOD event
  if (aodEvent) {
    TClonesArray *jetArray =
      dynamic_cast<TClonesArray*> (aodEvent->FindListObject(fJetBranchName));
    if (jetArray) {
      Int_t nJets = jetArray->GetEntriesFast();
      FillH1(kHistNoJets, nJets);

      AliAODJet *leadJet = 0x0;
      // AliAODJet *subleadJet = 0x0;

      for (Int_t iJet = 0; iJet < nJets; ++iJet) {
        AliAODJet *jet = (AliAODJet*) (*jetArray)[iJet];
        if (TMath::Abs(jet->Eta()) < 0.5) {
	  // check contributing tracks
	  Int_t nJetTracks = jet->GetRefTracks()->GetEntriesFast();
	  Int_t nJetTracks3 = 0;
	  AliAODTrack *leadingTrack = 0x0;
	  for (Int_t iTrack=0; iTrack < nJetTracks; ++iTrack) {
	    AliAODTrack *track = (AliAODTrack*) jet->GetRefTracks()->At(iTrack);

	    // count constituents above 3 GeV/c
	    if (track->Pt() > 3.)
	      ++nJetTracks3;

	    // find the leading track
	    if (!leadingTrack ||
		(track->Pt() > leadingTrack->Pt()))
	      leadingTrack = track;
	  }

          // find leading jet
          if (TMath::Abs(jet->Pt()) > leadingJetPtRec) {
            leadingJetPtRec = TMath::Abs(jet->Pt());
            // subleadJet = leadJet;
            leadJet = jet;
          }

	  // jet pt spectrum
	  for (Int_t iTrigger = 1; iTrigger < kTrgLast; ++iTrigger)
	    if (IsTrigger(Trigger_t (iTrigger)))
	      FillH1(kHistJetPt, jet->Pt(), iTrigger);

	  // integrated over all triggers
          FillH2(kHistJetPtNoTracks3, jet->Pt(), nJetTracks3);

          // limit to jets with leading track having an ITS contribution
          // with a hit in any SPD layer
          if (leadingTrack &&
              (leadingTrack->GetFlags() & AliVTrack::kITSrefit) &&
              (leadingTrack->HasPointOnITSLayer(0) || leadingTrack->HasPointOnITSLayer(1)))
	    for (Int_t iTrigger = 1; iTrigger < kTrgLast; ++iTrigger)
	      if (IsTrigger(Trigger_t (iTrigger)))
		FillH1(kHistJetPtITS, jet->Pt(), iTrigger);

          // limit to jets having 3 tracks above 3 GeV/c
          if (nJetTracks3 > 2)
	    for (Int_t iTrigger = 1; iTrigger < kTrgLast; ++iTrigger)
	      if (IsTrigger(Trigger_t (iTrigger)))
		FillH1(kHistJetPt3x3, jet->Pt(), iTrigger);
        }
      }
      // fill leading jet information
      for (Int_t iTrigger = 1; iTrigger < kTrgLast; ++iTrigger)
	if (IsTrigger(Trigger_t (iTrigger)))
	  FillH1(kHistLeadJetPt, leadJet ? leadJet->Pt() : 0., iTrigger);
    }
    else {
      printf("no jet array found as branch %s\n", fJetBranchName);
      aodEvent->Print();
    }
  } else {
    printf("no AOD event found\n");
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

  AliInputEventHandler *inputHandler =
    (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

  AliVEvent::EOfflineTriggerTypes physSel = (AliVEvent::EOfflineTriggerTypes) inputHandler->IsEventSelected();
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
  if ((physSel & (AliVEvent::kMB)))
    MarkTrigger(kTrgMinBias);

  if ((physSel & (AliVEvent::kINT7)))
    MarkTrigger(kTrgInt7);

  if ((physSel & (AliVEvent::kINT8)))
    MarkTrigger(kTrgInt8);

  if ((physSel & (AliVEvent::kEMC7)) &&
      trgClasses.Contains("CEMC7"))
    MarkTrigger(kTrgEMC7);

  if ((physSel & (AliVEvent::kEMC8)) &&
      trgClasses.Contains("CEMC8"))
    MarkTrigger(kTrgEMC8);

  if ((physSel & (AliVEvent::kEMCEJE)))
    MarkTrigger(kTrgEMCEJE);

  if ((physSel & (AliVEvent::kEMCEGA)))
    MarkTrigger(kTrgEMCEGA);

  // for the TRD-triggered events we use the classes
  if (trgClasses.Contains("CINT7WUHJT-"))
    MarkTrigger(kTrgInt7WUHJT);

  if (trgClasses.Contains("CINT8WUHJT-"))
    MarkTrigger(kTrgInt8WUHJT);

  if (trgClasses.Contains("CEMC7WUHJT-"))
    MarkTrigger(kTrgEMC7WUHJT);

  if (trgClasses.Contains("CEMC8WUHJT-"))
    MarkTrigger(kTrgEMC8WUHJT);

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
