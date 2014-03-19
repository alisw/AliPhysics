// ROOT
#include "TFile.h"
#include "TList.h"
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
  fTrdTrg(),
  fOutputList(),
  fHist(),
  fShortTaskId("jets_trg_trd"),
  fNoTriggers(kTrgLast),
  fNoJetPtBins(80),
  fJetPtBinMax(400),
  fAvgXsection(0.),
  fAvgTrials(0.),
  fPtHard(0.),
  fNTrials(0),
  fGlobalEfficiencyGTU(.8),
  fGtuLabel(-1) // -3 for hw, -1821 for re-simulation
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

  if (HasMC()) {
    fNoTriggers = kTrgMCLast;
    AddHistogram(ID(kHistXsection), "xsection stats",
		 1, 0., 1.);
    AddHistogram(ID(kHistPtHard), "pt hard;#hat{p}_{T} (GeV/c);counts",
		 fNoJetPtBins, 0., fJetPtBinMax);
    AddHistogram(ID(kHistJetPtMC), "leading jet spectrum (MC, |#eta| < 0.5);p_{T}^{jet} (GeV/c);counts",
		 fNoJetPtBins, 0., fJetPtBinMax);
  }

  AddHistogram(ID(kHistJetEtaAvg), "average eta",
	       100, -1., 1.);

  AddHistogram(ID(kHistNoJets), "number of jets;N^{jet};counts",
               400, -.5, 399.5);

  AddHistogram(ID(kHistTrackGTU), "GTU track p_{T};p_{T} (GeV/c);counts",
               100, 0., 25.);

  AddHistogram(ID(kHistTrackEffGTU), "found GTU tracks;p_{T} (GeV/c);#eta;#varphi",
               100, 0., 25.,
               100, -1., 1.,
               100, 0., 2.*TMath::Pi());
  AddHistogram(ID(kHistTrackEffMC), "wanted GTU tracks;p_{T} (GeV/c);#eta;#varphi",
               100, 0., 25.,
               100, -1., 1.,
               100, 0., 2.*TMath::Pi());

  AddHistogram(ID(kHistNPtMin),
	       "rejection;p_{T}^{min};N_{trk};trigger",
  	       100, 0., 10.,
	       20, 0., 20.,
	       fNoTriggers - 1, .5, fNoTriggers - .5);

  // leading jet
  AddHistogram(ID(kHistLeadJetPt),
	       "leading jet spectrum (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);counts",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5);
  AddHistogram(ID(kHistLeadJetPtEta),
	       "leading jet #eta (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger;#eta",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
	       40, -.5, .5);
  AddHistogram(ID(kHistLeadJetPtPhi),
	       "leading jet #phi (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger;#varphi",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
	       40, 0., 2. * TMath::Pi());
  AddHistogram(ID(kHistLeadJetEtaPhi),
	       "leading jet #eta - #varphi (|#eta| < 0.5);#eta;trigger;#varphi",
  	       40, -.5, .5,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
	       40, 0., 2. * TMath::Pi());

  AddHistogram(ID(kHistLeadJetPtTrackPt),
	       "leading jet pt vs track pt;p_{T} (GeV/c);trigger;p_{T}^{jet,ch} (GeV/c)",
	       100., 0., 20.,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistLeadJetPtZ),
	       "leading jet pt vs z;z;trigger;p_{T}^{jet,ch} (GeV/c)",
	       100., 0., 1.,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistLeadJetPtXi),
	       "leading jet pt vs #xi;#xi;trigger;p_{T}^{jet,ch} (GeV/c)",
	       100., 0., 10.,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
  	       fNoJetPtBins, 0., fJetPtBinMax);

  // inclusive jets
  AddHistogram(ID(kHistJetPt),
	       "jet spectrum (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5);
  AddHistogram(ID(kHistJetPtEta),
	       "jet #eta (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger;#eta",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
	       40, -.5, .5);
  AddHistogram(ID(kHistJetPtPhi),
	       "jet #phi (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger;#varphi",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
	       40, 0., 2. * TMath::Pi());
  AddHistogram(ID(kHistJetEtaPhi),
	       "jet #eta - #varphi (|#eta| < 0.5);#eta;trigger;#varphi",
  	       40, -.5, .5,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
	       40, 0., 2. * TMath::Pi());

  AddHistogram(ID(kHistJetPtITS),
	       "jet spectrum (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5);
  AddHistogram(ID(kHistJetPt3x3),
	       "jet spectrum (|#eta| < 0.5);p_{T}^{jet,ch} (GeV/c);trigger",
  	       fNoJetPtBins, 0., fJetPtBinMax,
	       fNoTriggers - 1, .5, fNoTriggers - .5);

  AddHistogram(ID(kHistJetPtTrackPt),
	       "jet pt vs track pt;p_{T} (GeV/c);trigger;p_{T}^{jet,ch} (GeV/c)",
	       100., 0., 20.,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistJetPtZ),
	       "jet pt vs z;z;trigger;p_{T}^{jet,ch} (GeV/c)",
	       100., 0., 1.,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
  	       fNoJetPtBins, 0., fJetPtBinMax);
  AddHistogram(ID(kHistJetPtXi),
	       "jet pt vs #xi;#xi;trigger;p_{T}^{jet,ch} (GeV/c)",
	       100., 0., 10.,
	       fNoTriggers - 1, .5, fNoTriggers - .5,
  	       fNoJetPtBins, 0., fJetPtBinMax);

  for (Int_t iHist = kHistLeadJetPt; iHist <= kHistJetPtXi; ++iHist) {
    TH1 *h = GetHistogram(Hist_t (iHist));
    h->GetYaxis()->SetBinLabel(ID(kTrgMinBias));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt7));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt8));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC7));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC8));
    h->GetYaxis()->SetBinLabel(ID(kTrgPbPb));
    h->GetYaxis()->SetBinLabel(ID(kTrgCentral));
    h->GetYaxis()->SetBinLabel(ID(kTrgSemiCentral));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt7WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt8WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC7WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMC8WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMCEJE));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMCEGA));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt7_WU));
    h->GetYaxis()->SetBinLabel(ID(kTrgInt7_WUHJT));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMCEJE_WU));
    h->GetYaxis()->SetBinLabel(ID(kTrgEMCEJE_WUHJT));
    if (HasMC()) {
      h->GetYaxis()->SetBinLabel(ID(kTrgMC3x3Vtx));
      h->GetYaxis()->SetBinLabel(ID(kTrgMC3x3TRD));
      h->GetYaxis()->SetBinLabel(ID(kTrgMC3x3TRDeff));
      h->GetYaxis()->SetBinLabel(ID(kTrgMC3x3TRDeffmap));
    }
  }

  AddHistogram(ID(kHistJetPtNoTracks3),
	       "number of tracks above 3 GeV;p_{T}^{jet,ch};no. of tracks",
               fNoJetPtBins, 0., fJetPtBinMax,
               40, -.5, 39.5);

  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskJetsTriggerTRD::Notify()
{
  // actions to be taken upon notification about input file change

  // ??? check ???
  // we should only count the cross section we see in the analysis,
  // i.e. fXsection / nTrials

  if (HasMC()) {
    fAvgXsection = 0.;
    Float_t xSection = 0.;
    Float_t nTrials = 0.;
    Float_t nEntries = 0.;

    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();

    if (tree) {
      TFile *curfile = tree->GetCurrentFile();
      if (!curfile) {
	AliError("No current file");
	return kFALSE;
      }

      if (!AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(), xSection, nTrials)) {
	AliError("retrieval of cross section failed");
      }

      nEntries = (Float_t) tree->GetTree()->GetEntries();
      if (nEntries > 0.) {
	fAvgTrials = nTrials / nEntries;
	if (nTrials > 0.)
	  fAvgXsection = xSection / (nTrials / nEntries);
      }
    }
    printf("n_trials = %f\n", nTrials);
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

  Int_t nTracksMC  = mcEvent ? mcEvent->GetNumberOfTracks() : 0; // no. of MC tracks
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
  if (HasMC()) {
    if (fAvgXsection == 0.) {
      AliError("zero cross section");
      return;
    }
    FillH1(kHistXsection, .5, fAvgXsection);
  }

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
      Int_t nTrials = pythiaHeader->Trials();
      fNTrials += nTrials;

      FillH1(kHistPtHard, fPtHard);

      // loop over jets from PYTHIA
      Float_t eta1 = -10.;
      Float_t eta2 = -10.;
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

	// for eta averge determination consider only jets above 60 GeV
	if (pt < 60.)
	  continue;
	if (eta1 < -1.)
	  eta1 = eta;
	else if (eta2 < -1.)
	  eta2 = eta;
      }
      // fill histogram for leading jet pt spectrum
      FillH1(kHistJetPtMC, leadingJetPtMC);
      // fill histogram for eta average
      if ((eta1 > -1.) && (eta2 > -1.))
	FillH1(kHistJetEtaAvg, (eta1 + eta2)/2.);
    }

    for (Int_t iTrack = 0; iTrack < nTracksMC; ++iTrack) {
      AliVParticle *part = mcEvent->GetTrack(iTrack);
      if (AcceptTrackMC(iTrack)) {
	FillH3(kHistTrackEffMC, part->Pt(), part->Eta(), part->Phi());
      }
    }
  }

  // loop over GTU tracks
  for (Int_t iTrack = 0; iTrack < nTracksTRD; ++iTrack) {
    AliVTrdTrack *trk = InputEvent()->GetTrdTrack(iTrack);
    // printf("trk %p has pt %5.2f and label %i\n",
    // 	   trk, trk->Pt(), trk->GetLabel());
    if ((fGtuLabel != -1) && (trk->GetLabel() != fGtuLabel))
      continue;
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

	for (Int_t iTrigger = 1; iTrigger < fNoTriggers; ++iTrigger)
	  if (IsTrigger(Trigger_t (iTrigger)))
	      FillH3(kHistNPtMin,
		     TMath::Abs(trdTrack->Pt()), nTracksPerStackMax, iTrigger);
      }
      if (HasMC()) {
	Int_t label = trdTrack->GetLabel();
	if (label > -1) {
	  if (AcceptTrackMC(label)) {
	    AliVParticle *part = MCEvent()->GetTrack(label);
	    FillH3(kHistTrackEffGTU, part->Pt(), part->Eta(), part->Phi());
	  }
	}
	else {
	  AliWarning(Form("GTU track at %p with no label", trdTrack));
	  const Int_t nLayers = 6;
	  for (Int_t iLayer = 0; iLayer < nLayers; ++iLayer) {
	    AliVTrdTracklet *trkl = trdTrack->GetTracklet(iLayer);
	    if (trkl)
	      AliWarning(Form("tracklet in layer %i has label %i\n",
			      iLayer, trkl->GetLabel()));
	  }
	}
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

	  // jet pt spectrum and
	  // fragmentation function
	  for (Int_t iTrigger = 1; iTrigger < fNoTriggers; ++iTrigger)
	    if (IsTrigger(Trigger_t (iTrigger))) {
	      Float_t jetPt = jet->Pt();

	      FillH1(kHistJetPt, jetPt, iTrigger);

	      FillH3(kHistJetPtEta, jet->Pt(), iTrigger, jet->Eta());
	      FillH3(kHistJetPtPhi, jet->Pt(), iTrigger, jet->Phi());
	      if (jet->Pt() > 50.)
		FillH3(kHistJetEtaPhi, jet->Eta(), iTrigger, jet->Phi());

	      Int_t nTracks = jet->GetRefTracks()->GetEntriesFast();
	      for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
		AliVParticle *track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(iTrack));
		if (track) {
		  Float_t trackPt = track->Pt();
		  Float_t z = trackPt / jetPt;
		  Float_t xi = - TMath::Log(z);
		  FillH3(kHistJetPtTrackPt, trackPt, iTrigger, jet->Pt());
		  FillH3(kHistJetPtZ, z, iTrigger, jet->Pt());
		  FillH3(kHistJetPtXi, xi, iTrigger, jet->Pt());
		}
	      }
	    }

	  // integrated over all triggers
          FillH2(kHistJetPtNoTracks3, jet->Pt(), nJetTracks3);

          // limit to jets with leading track having an ITS contribution
          // with a hit in any SPD layer
          if (leadingTrack &&
              (leadingTrack->GetFlags() & AliVTrack::kITSrefit) &&
              (leadingTrack->HasPointOnITSLayer(0) || leadingTrack->HasPointOnITSLayer(1)))
	    for (Int_t iTrigger = 1; iTrigger < fNoTriggers; ++iTrigger)
	      if (IsTrigger(Trigger_t (iTrigger)))
		FillH1(kHistJetPtITS, jet->Pt(), iTrigger);

          // limit to jets having 3 tracks above 3 GeV/c
          if (nJetTracks3 > 2)
	    for (Int_t iTrigger = 1; iTrigger < fNoTriggers; ++iTrigger)
	      if (IsTrigger(Trigger_t (iTrigger)))
		FillH1(kHistJetPt3x3, jet->Pt(), iTrigger);
        }
      }

      // fill leading jet information
      for (Int_t iTrigger = 1; iTrigger < fNoTriggers; ++iTrigger)
	if (IsTrigger(Trigger_t (iTrigger))) {
	  FillH2(kHistLeadJetPt, leadJet ? leadJet->Pt() : 0., iTrigger);
	  if (leadJet) {
	    Float_t leadJetPt = leadJet->Pt();

	    FillH3(kHistLeadJetPtEta, leadJet->Pt(), iTrigger, leadJet->Eta());
	    FillH3(kHistLeadJetPtPhi, leadJet->Pt(), iTrigger, leadJet->Phi());
	    if (leadJet->Pt() > 50.)
	      FillH3(kHistLeadJetEtaPhi, leadJet->Eta(), iTrigger, leadJet->Phi());

	    Int_t nTracks = leadJet->GetRefTracks()->GetEntriesFast();
	    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
	      AliVParticle *track = dynamic_cast<AliVParticle*>(leadJet->GetRefTracks()->At(iTrack));
	      if (track) {
		Float_t trackPt = track->Pt();
		Float_t z = trackPt / leadJetPt;
		Float_t xi = - TMath::Log(z);
		FillH3(kHistLeadJetPtTrackPt, trackPt, iTrigger, leadJet->Pt());
		FillH3(kHistLeadJetPtZ, z, iTrigger, leadJet->Pt());
		FillH3(kHistLeadJetPtXi, xi, iTrigger, leadJet->Pt());
	      }
	    }
	  }
	}
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

  printf("total trials: %d\n", fNTrials);
}

Bool_t AliAnalysisTaskJetsTriggerTRD::DetectTriggers()
{
  fTriggerMask = 0;

  AliInputEventHandler *inputHandler =
    (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

  AliVEvent::EOfflineTriggerTypes physSel = (AliVEvent::EOfflineTriggerTypes) inputHandler->IsEventSelected();
  TString trgClasses = InputEvent()->GetFiredTriggerClasses();

  fTrdTrg.CalcTriggers(InputEvent());

  UInt_t inputMaskL1 = 0;
  if (AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*> (InputEvent()))
    inputMaskL1 = esdEvent->GetHeader()->GetL1TriggerInputs();
  else if (AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*> (InputEvent()))
    inputMaskL1 = aodEvent->GetHeader()->GetL1TriggerInputs();

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

  if (fDebug > 2)
    printf("trg: %8s %8s %8s %8s (%s)\n",
	   fTrdTrg.IsFired(AliTRDTriggerAnalysis::kHJT) ? "HJT"  : " ",
	   fTrdTrg.IsFired(AliTRDTriggerAnalysis::kHSE) ? "HSE"  : " ",
	   fTrdTrg.IsFired(AliTRDTriggerAnalysis::kHQU) ? "HQU"  : " ",
	   fTrdTrg.IsFired(AliTRDTriggerAnalysis::kHEE) ? "HEE"  : " ",
	   trgClasses.Data()
	   );

  // physics selection
  if ((physSel & (AliVEvent::kMB)))
    MarkTrigger(kTrgMinBias);

  if ((physSel & (AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral)))
    MarkTrigger(kTrgPbPb);

  if ((physSel & (AliVEvent::kINT7))) {
    MarkTrigger(kTrgInt7);
    if (trgClasses.Contains("WU")) {
      MarkTrigger(kTrgInt7_WU);
      if (inputMaskL1 & (1 << 9))
	MarkTrigger(kTrgInt7_WUHJT);
    }
  }

  if ((physSel & (AliVEvent::kINT8)))
    MarkTrigger(kTrgInt8);

  if ((physSel & (AliVEvent::kEMC7)) &&
      trgClasses.Contains("CEMC7"))
    MarkTrigger(kTrgEMC7);

  if ((physSel & (AliVEvent::kEMC8)) &&
      trgClasses.Contains("CEMC8"))
    MarkTrigger(kTrgEMC8);

  if ((physSel & (AliVEvent::kEMCEJE))) {
    MarkTrigger(kTrgEMCEJE);
    if (trgClasses.Contains("WU")) {
      MarkTrigger(kTrgEMCEJE_WU);
      if (inputMaskL1 & (1 << 9))
	MarkTrigger(kTrgEMCEJE_WUHJT);
    }
  }

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

  if (HasMC())
    DetectMCTriggers();

  return kTRUE;
}

Bool_t AliAnalysisTaskJetsTriggerTRD::DetectMCTriggers()
{
  AliMCEvent  *mcEvent   = MCEvent();

  Int_t nTracksMC  = mcEvent ? mcEvent->GetNumberOfTracks() : 0; // no. of MC tracks
  TList partMC; // list of particles

  // fill tracks passing cuts to list
  for (Int_t iTrackMC = 0; iTrackMC < nTracksMC; ++iTrackMC) {
    AliVParticle *part = mcEvent->GetTrack(iTrackMC);

    // only consider primaries
    if (!mcEvent->IsPhysicalPrimary(iTrackMC))
      continue;

    // only consider charged particles
    if (part->Charge() == 0)
      continue;

    // only look at particles in eta-acceptance of the central barrel
    if (TMath::Abs(part->Eta()) >= 0.9)
      continue;

    partMC.Add(part);
  }

  // sort, starting from highest pt
  partMC.Sort(kSortAscending);

  Int_t nTracksInStack[4][90] = { { 0 } };

  // iterate over accepted tracks
  TIter nextMcPart(&partMC);
  while (AliVParticle *part = (AliMCParticle*) (nextMcPart())) {
    Float_t pt  = part->Pt();
    Float_t eta = part->Eta();
    Float_t phi = part->Phi();

    Int_t phiBin = (Int_t) (phi / (2 / 18. *TMath::Pi()));
    // ??? use actual geometry
    Int_t etaBin = (Int_t) ((eta + .9) / (1.8 / 5.));

    Int_t pdgCode = ((AliMCParticle*) part)->Particle()->GetPdgCode();
    if (fDebug > 2)
      printf("pt = %f, eta = %f, phi = %f, phiBin = %d, etaBin = %d, pdgCode = %d\n",
	     pt, eta, phi, phiBin, etaBin, pdgCode);

    // increment number of tracks in stack
    Int_t bin = phiBin * 5 + etaBin;
    ++nTracksInStack[0][bin];

    // if minimum number reached and
    // the track has pt above threshold
    // trigger fired
    if ((nTracksInStack[0][bin] >= 3) &&
	(pt >= 3.)) 
      MarkTrigger(kTrgMC3x3Vtx);

    // only propagate particles that do not decay before ???

    // propagate to TRD
    // why not use track references?
    Float_t rTrack = pt / .15;
    Float_t rTrd   = 3.;
    Float_t phiTrd = phi + TMath::ASin(rTrd / 2 / rTrack);
    if (phiTrd < 0)
      phiTrd += 2*TMath::Pi();
    else if (phiTrd > 2*TMath::Pi())
      phiTrd -= 2*TMath::Pi();
    Float_t etaTrd = eta;

    Int_t secTrd   = (Int_t) (phiTrd / (2.*TMath::Pi()/18.));
    Int_t stackTrd = (Int_t) ((etaTrd+0.9) / (1.8/5.));
    Int_t binTrd   = secTrd * 5 + stackTrd;

    // increment number of tracks in stack
    ++nTracksInStack[1][binTrd];
    if (gRandom->Uniform() < fGlobalEfficiencyGTU)
      ++nTracksInStack[2][binTrd];
    if (gRandom->Uniform() < GetEfficiencyTRD(pt, eta, phi))
      ++nTracksInStack[3][binTrd];

    // if minimum number reached and
    // the track has pt above threshold
    // trigger fired
    if (pt >= 3.) {
      if (nTracksInStack[1][binTrd] >= 3)
	MarkTrigger(kTrgMC3x3TRD);
      if (nTracksInStack[2][binTrd] >= 3)
	MarkTrigger(kTrgMC3x3TRDeff);
      if (nTracksInStack[3][binTrd] >= 3)
	MarkTrigger(kTrgMC3x3TRDeffmap);
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskJetsTriggerTRD::AcceptTrackMC(Int_t track) const
{
  // only consider primaries
  if (!MCEvent()->IsPhysicalPrimary(track))
    return kFALSE;

  const AliVParticle *part = MCEvent()->GetTrack(track);

  // only consider charged particles
  if (part->Charge() == 0)
    return kFALSE;

  // only look at particles in eta-acceptance of the central barrel
  if (TMath::Abs(part->Eta()) >= 0.9)
    return kFALSE;

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
