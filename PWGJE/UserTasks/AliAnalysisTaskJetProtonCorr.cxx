// ROOT
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFormula.h"

// analysis framework
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrdTrack.h"
#include "AliVVertex.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"

// MC stuff
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"

// ESD stuff
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
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

#include "AliAnalysisTaskJetProtonCorr.h"

AliAnalysisTaskJetProtonCorr::AliAnalysisTaskJetProtonCorr(const char *name) :
  AliAnalysisTaskSE(name),
  fMCEvent(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTriggerMask(0),
  fClassMask(0),
  fCentrality(100.),
  fCentralityCheck(100.),
  fZvtx(0.),
  fPIDResponse(0x0),
  fEventPlane(5.),
  fPrimTrackArray(0x0),
  fJetArray(0x0),
  fPoolMgr(),
  fPool(),
  fHistCorr(0x0),
  fOutputList(),
  fHist(),
  fShortTaskId("jet_prot_corr"),
  fCutsPrim(0x0),
  fTrgPartPtMin(6.),
  fTrgPartPtMax(8.),
  fTrgJetPtMin(50.),
  fTrgJetPtMax(80.),
  fAssPartPtMin(2.),
  fAssPartPtMax(4.)
{
  // default ctor

  fkCorrTypeName[kCorrHadHad]  = "hh";
  fkCorrTypeName[kCorrHadProt] = "hp";
  fkCorrTypeName[kCorrJetHad]  = "jh";
  fkCorrTypeName[kCorrJetProt] = "jp";

  fkClassName[kClCentral]      = "cent";
  fkClassName[kClSemiCentral]  = "semi";
  fkClassName[kClDijet]        = "dijet";

  fkEvName[kEvSame] = "same";
  fkEvName[kEvMix]  = "mixed";

  // track cuts
  fCutsPrim = new AliESDtrackCuts();

  // this is taken from PWGJE track cuts
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  fCutsPrim->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  fCutsPrim->SetMinNClustersTPC(70);
  fCutsPrim->SetMaxChi2PerClusterTPC(4);
  fCutsPrim->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  fCutsPrim->SetAcceptKinkDaughters(kFALSE);
  fCutsPrim->SetRequireTPCRefit(kTRUE);
  fCutsPrim->SetMaxFractionSharedTPCClusters(0.4);
  // ITS
  fCutsPrim->SetRequireITSRefit(kTRUE);
  //accept secondaries
  fCutsPrim->SetMaxDCAToVertexXY(2.4);
  fCutsPrim->SetMaxDCAToVertexZ(3.2);
  fCutsPrim->SetDCAToVertex2D(kTRUE);
  //reject fakes
  fCutsPrim->SetMaxChi2PerClusterITS(36);
  fCutsPrim->SetMaxChi2TPCConstrainedGlobal(36);

  fCutsPrim->SetRequireSigmaToVertex(kFALSE);

  fCutsPrim->SetEtaRange(-0.9,0.9);
  fCutsPrim->SetPtRange(0.15, 1E+15);

  fCutsPrim->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

  // event mixing pool
  Double_t centralityBins[] = { 0., 10., 40., 80.};
  Int_t nCentralityBins = sizeof(centralityBins)/sizeof(centralityBins[0]);

  Double_t vertexBins[] = { -10., -5., 0., 5., 10.};
  Int_t nVertexBins = sizeof(vertexBins)/sizeof(vertexBins[0]);

  Double_t psiBins[12];
  Int_t nPsiBins = sizeof(psiBins)/sizeof(psiBins[0]);
  for (Int_t iBin = 0; iBin < nPsiBins; ++iBin)
    psiBins[iBin] = iBin * TMath::Pi()/nPsiBins;

  for (Int_t iTrg = 0; iTrg < kTrgLast; ++iTrg) {
    for (Int_t iAss = 0; iAss < kAssLast; ++iAss) {
      GetPoolMgr((Trg_t) iTrg, (Ass_t) iAss) =
	new AliEventPoolManager(10, 100,
				nCentralityBins, centralityBins,
				nVertexBins, vertexBins,
				nPsiBins, psiBins);
      GetPoolMgr((Trg_t) iTrg, (Ass_t) iAss)->SetTargetValues(100, .1, 1);
    }
  }

  fHistCorr = new AliHistCorr*[kEvLast*kCorrLast*kClLast];

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetProtonCorr::~AliAnalysisTaskJetProtonCorr()
{
  // dtor

  // delete [] fHistCorr;
}

void AliAnalysisTaskJetProtonCorr::UserCreateOutputObjects()
{
  // create user output objects

  // setup list
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  // setup histograms
  TH1 *hist;
  TH1 *histStat = AddHistogram(ID(kHistStat), "event statistics;;counts",
                               kStatLast-1, .5, kStatLast-.5);
  histStat->GetXaxis()->SetBinLabel(ID(kStatSeen));
  histStat->GetXaxis()->SetBinLabel(ID(kStatTrg));
  histStat->GetXaxis()->SetBinLabel(ID(kStatEvCuts));
  histStat->GetXaxis()->SetBinLabel(ID(kStatUsed));
  histStat->GetXaxis()->SetBinLabel(ID(kStatCent));
  histStat->GetXaxis()->SetBinLabel(ID(kStatEvPlane));
  histStat->GetXaxis()->SetBinLabel(ID(kStatPID));

  AddHistogram(ID(kHistCentrality), "centrality;C;counts",
	       110, -5., 105.);
  hist = AddHistogram(ID(kHistCentralityUsed), "centrality used;C;event class",
                      110, -5., 105.,
                      kClLast, -.5, kClLast-.5);
  hist->GetYaxis()->SetBinLabel(LAB(kClCentral));
  hist->GetYaxis()->SetBinLabel(LAB(kClSemiCentral));
  hist->GetYaxis()->SetBinLabel(LAB(kClDijet));
  AddHistogram(ID(kHistCentralityCheck), "centrality check;C;counts",
	       110, -5., 105.);
  hist = AddHistogram(ID(kHistCentralityCheckUsed), "centrality check used;C;event class",
                      110, -5., 105.,
                      kClLast, -.5, kClLast-.5);
  hist->GetYaxis()->SetBinLabel(LAB(kClCentral));
  hist->GetYaxis()->SetBinLabel(LAB(kClSemiCentral));
  hist->GetYaxis()->SetBinLabel(LAB(kClDijet));

  AddHistogram(ID(kHistNsigmaTPCTOF), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -5., 5.,
               100, -5., 5.);

  AddHistogram(ID(kHistEvPlane), "event plane angle;#Psi;counts",
               100, -0. * TMath::Pi(), 1. * TMath::Pi());
  AddHistogram(ID(kHistEvPlaneUsed), "event plane angle;#Psi;counts",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);

  AddHistogram(ID(kHistJetPt), "jet spectrum",
               40, 0., 200.);

  AddHistogram(ID(kHistEtaPhiTrgHad), "trg had;#varphi;#eta",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100, -2., 2.);
  AddHistogram(ID(kHistEtaPhiTrgJet), "trg jet;#varphi;#eta",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100, -2., 2.);
  AddHistogram(ID(kHistEtaPhiAssHad), "ass had;#varphi;#eta",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100, -2., 2.);
  AddHistogram(ID(kHistEtaPhiAssProt), "ass proton;#varphi;#eta",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100, -2., 2.);

  for (Int_t iCorr = 0; iCorr < kCorrLast; ++iCorr) {
    for (Int_t iCl = 0; iCl < kClLast; ++iCl) {
      for (Int_t iEv = 0; iEv < kEvLast; ++iEv) {
  	GetHistCorr((CorrType_t) iCorr, (Class_t) iCl, (Ev_t) iEv) =
  	  new AliHistCorr(Form("corr_%s_%s_%s", fkCorrTypeName[iCorr], fkClassName[iCl], fkEvName[iEv]), fOutputList);
      }
    }
  }

  PostData(1, fOutputList);
}

Bool_t AliAnalysisTaskJetProtonCorr::Notify()
{
  // actions to be taken upon notification about input file change

  return AliAnalysisTaskSE::Notify();
}

void AliAnalysisTaskJetProtonCorr::UserExec(Option_t * /* option */)
{
  // actual work

  // setup pointers to input data (null if unavailable)
  // mcEvent:  MC input
  // esdEvent: ESD input
  // outEvent: AOD output
  // aodEvent: AOD input if available, otherwise AOD output

  fMCEvent   = this->MCEvent();
  fESDEvent  = dynamic_cast<AliESDEvent*>(this->InputEvent()); // could also be AOD input
  AliAODEvent* outEvent  = this->AODEvent();
  fAODEvent  = dynamic_cast<AliAODEvent*> (this->InputEvent());
  if (!fAODEvent)
    fAODEvent = outEvent;

  if ((fDebug > 0) && fESDEvent)
    printf("event: %s-%06i\n", CurrentFileName(), fESDEvent->GetEventNumberInFile());

  // record number of sampled events and detect trigger contributions
  FillH1(kHistStat, kStatSeen);
  if (!DetectTriggers()) {
    AliError("Failed to detect the triggers");
    return;
  }

  if (!IsTrigger(kTriggerInt))
    return;

  FillH1(kHistStat, kStatTrg);

  // prepare the event
  // (make sure it is cleaned up in the end)
  if (PrepareEvent()) {
    FillH1(kHistStat, kStatUsed);
    FillH1(kHistCentrality, fCentrality);
    FillH1(kHistCentralityCheck, fCentralityCheck);

    // event cuts
    if (TMath::Abs(fZvtx) > 10.)
      goto stop;
    if (GetCentrality() > 90.)
      goto stop;

    FillH1(kHistStat, kStatEvCuts);

    // event category
    DetectClasses();

    FillH1(kHistEvPlane, fEventPlane);
    for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
      if (IsClass((Class_t) iClass)) {
        FillH2(kHistCentralityUsed, fCentrality, iClass);
        FillH2(kHistCentralityCheckUsed, fCentralityCheck, iClass);
	FillH2(kHistEvPlaneUsed, fEventPlane, iClass);
      }
    }

    // select trigger particles and potential associated particles/protons
    TObjArray trgArray[kTrgLast];
    TObjArray assArray[kAssLast];

    Int_t nPrimTracks = fPrimTrackArray ? fPrimTrackArray->GetEntries() : 0;
    for (Int_t iTrack = 0; iTrack < nPrimTracks; ++iTrack) {
      AliVTrack *trk = (AliVTrack*) fPrimTrackArray->At(iTrack);
      FillH3(kHistNsigmaTPCTOF,
             trk->Pt(),
             fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
             fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
      if (AcceptTrigger(trk)) {
	trgArray[kTrgHad].Add(trk);
	FillH1(kHistEtaPhiTrgHad, trk->Phi(), trk->Eta());
      }
      if (AcceptAssoc(trk)) {
	assArray[kAssHad].Add(trk);
	FillH1(kHistEtaPhiAssHad, trk->Phi(), trk->Eta());
	if (IsProton(trk)) {
	  assArray[kAssProt].Add(trk);
	  FillH1(kHistEtaPhiAssProt, trk->Phi(), trk->Eta());
	}
      }
    }

    // select trigger jet
    Int_t nJets = fJetArray ? fJetArray->GetEntries() : 0;
    for (Int_t iJet = 0; iJet < nJets; ++iJet) {
      AliAODJet *jet = (AliAODJet*) fJetArray->At(iJet);
      FillH1(kHistJetPt, jet->Pt());
      if (AcceptTrigger(jet)) {
	trgArray[kTrgJet].Add(jet);
	FillH1(kHistEtaPhiTrgJet, jet->Phi(), jet->Eta());
      }
    }

    // correlate, both same and mixed event
    for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
      if (IsClass((Class_t) iClass)) {
	for (Int_t iTrg = 0; iTrg < kTrgLast; ++iTrg) {
	  for (Int_t iAss = 0; iAss < kAssLast; ++iAss) {
	    // same event
	    Correlate((Trg_t) iTrg, (Ass_t) iAss, (Class_t) iClass, kEvSame, &trgArray[iTrg], &assArray[iAss]);

	    // mixed event
	    AliEventPool *pool = GetPool((Class_t) iClass, (Trg_t) iTrg, (Ass_t) iAss);
	    if (pool && pool->IsReady()) {
	      // printf("----- using pool: %i %i %i -----\n", iClass, iTrg, iAss);
	      Int_t nEvents = pool->GetCurrentNEvents();
	      for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {
		TObjArray *assTracks = pool->GetEvent(iEvent);
		Correlate((Trg_t) iTrg, (Ass_t) iAss, (Class_t) iClass, kEvMix, &trgArray[iTrg], assTracks, 1./nEvents);
	      }
	    }
	    // if (pool && !pool->IsReady()) {
	    //   printf("----- pool not ready: %i %i %i -----\n", iClass, iTrg, iAss);
	    //   pool->PrintInfo();
	    // }
	  }

	  // fill event pool for mixing
	  if (trgArray[iTrg].GetEntries() > 0) {
	    for (Int_t iAss = 0; iAss < kAssLast; ++iAss) {
	      AliEventPool *pool = GetPool((Class_t) iClass, (Trg_t) iTrg, (Ass_t) iAss);
	      if (pool) {
		pool->UpdatePool(CloneTracks(&assArray[iAss]));
		// printf("----- updating pool: %i %i %i -----\n", iClass, iTrg, iAss);
		// pool->PrintInfo();
	      }
	    }
	  }
	}
      }
    }
  }

 stop:
  CleanUpEvent();

  PostData(1, fOutputList);
}

void AliAnalysisTaskJetProtonCorr::Terminate(const Option_t * /* option */)
{
  // actions at task termination

}

Bool_t AliAnalysisTaskJetProtonCorr::DetectTriggers()
{
  fTriggerMask = 0;

  AliVEvent::EOfflineTriggerTypes physSel = (AliVEvent::EOfflineTriggerTypes) ((AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  TString trgClasses = InputEvent()->GetFiredTriggerClasses();

  // physics selection
  if (physSel & (AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral))
    MarkTrigger(kTriggerInt);

  return kTRUE;
}

Bool_t AliAnalysisTaskJetProtonCorr::DetectClasses()
{
  fClassMask = 0;

  if (IsCentral())
    MarkClass(kClCentral);

  if (IsSemiCentral())
    MarkClass(kClSemiCentral);

  return kTRUE;
}

Bool_t AliAnalysisTaskJetProtonCorr::PrepareEvent()
{
  Bool_t eventGood = kTRUE;

  // retrieve z-vertex position
  const AliVVertex *vtx = InputEvent()->GetPrimaryVertex();
  if (vtx && (vtx->GetNContributors() >= 3.))
    fZvtx = vtx->GetZ();
  else
    fZvtx = 100.;

  // retrieve centrality
  AliCentrality *eventCentrality = InputEvent()->GetCentrality();
  if (eventCentrality) {
    fCentrality = eventCentrality->GetCentralityPercentile("V0M");
    fCentralityCheck = eventCentrality->GetCentralityPercentile("TRK");
    if (fCentrality >= 0.) {
      FillH1(kHistStat, kStatCent);
    } else {
      // centrality estimation not reliable
      eventGood = kFALSE;
      fCentrality = 105.;
    }
  }

  // retrieve event plane
  AliEventplane *eventPlane = InputEvent()->GetEventplane();
  if (eventPlane) {
    FillH1(kHistStat, kStatEvPlane);
    fEventPlane = eventPlane->GetEventplane("Q");
  }
  else
    eventGood = kFALSE;

  // retrieve PID
  fPIDResponse = ((AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  if (fPIDResponse)
    FillH1(kHistStat, kStatPID);
  else
    eventGood = kFALSE;

  // retrieve primary tracks
  if (fESDEvent) {
    fPrimTrackArray = fCutsPrim->GetAcceptedTracks(fESDEvent);
  }
  else if (fAODEvent) {
    fPrimTrackArray = new TObjArray();
    Int_t nTracksAOD = fAODEvent->GetNumberOfTracks();
    for (Int_t iTrack = 0; iTrack < nTracksAOD; ++iTrack) {
      AliAODTrack *trk = fAODEvent->GetTrack(iTrack);
      // track cuts esdTrackCutsH ???
      if (trk->TestFilterMask(1 << 4))
        fPrimTrackArray->Add(trk);
    }
  }
  else
    eventGood = kFALSE;

  // retrieve jet array
  if (fAODEvent) {
    fJetArray = dynamic_cast<TClonesArray*> (fAODEvent->FindListObject(fJetBranchName));
    if (!fJetArray) {
      printf("no jet branch \"%s\" found, in the AODs are:\n", fJetBranchName);
      if (fDebug > 0)
	fAODEvent->GetList()->Print();
    }
  }

  // retrieve event pool for the event category
  if (eventGood) {
    for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
      for (Int_t iTrg = 0; iTrg < kTrgLast; ++iTrg) {
	for (Int_t iAss = 0; iAss < kAssLast; ++iAss) {
	  AliEventPoolManager *mgr = GetPoolMgr((Trg_t) iTrg, (Ass_t) iAss);
	  GetPool((Class_t) iClass, (Trg_t) iTrg, (Ass_t) iAss) =
	    mgr ? mgr->GetEventPool(fCentrality, fZvtx, fEventPlane) : 0x0;
	}
      }
    }
  }

  return eventGood;
}

Bool_t AliAnalysisTaskJetProtonCorr::AcceptTrigger(AliVTrack *trg)
{
  if ((trg->Pt() > fTrgPartPtMin) && (trg->Pt() < fTrgPartPtMax)) {
    if (IsCentral())
      return kTRUE;
    else if (IsSemiCentral()) {
      // check calculation ???
      Float_t deltaPhi = trg->Phi() - GetEventPlane();
      if (deltaPhi > TMath::Pi())
        deltaPhi -= 2. * TMath::Pi();
      else if (deltaPhi < -TMath::Pi())
        deltaPhi += 2. * TMath::Pi();
      if (TMath::Abs(TMath::Pi() - TMath::Abs(deltaPhi)) >= (3./8. * TMath::Pi()))
        return kTRUE;
    }
  }

  return kFALSE;
}

Bool_t AliAnalysisTaskJetProtonCorr::AcceptTrigger(AliAODJet *trg)
{
  if ((trg->Pt() > fTrgJetPtMin) && (trg->Pt() < fTrgJetPtMax)) {
    if (IsCentral())
      return kTRUE;
    else if (IsSemiCentral()) {
      // check calculation ???
      Float_t deltaPhi = trg->Phi() - GetEventPlane();
      if (deltaPhi > TMath::Pi())
        deltaPhi -= 2. * TMath::Pi();
      else if (deltaPhi < -TMath::Pi())
        deltaPhi += 2. * TMath::Pi();
      if (TMath::Abs(TMath::Pi() - TMath::Abs(deltaPhi)) >= (3./8. * TMath::Pi()))
        return kTRUE;
    }
  }

  return kFALSE;
}

Bool_t AliAnalysisTaskJetProtonCorr::AcceptAssoc(AliVTrack *trg)
{
  if ((trg->Pt() > fAssPartPtMin) && (trg->Pt() < fAssPartPtMax))
    return kTRUE;

  return kFALSE;
}

Bool_t AliAnalysisTaskJetProtonCorr::IsProton(AliVTrack *trk)
{
  Double_t nSigmaProtonTPC = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton);
  Double_t nSigmaProtonTOF = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton);

  if ((TMath::Abs(nSigmaProtonTPC) <= 2.) && (TMath::Abs(nSigmaProtonTOF) <= 2.)) {
    return kTRUE;
  }

  return kFALSE;
}

Bool_t AliAnalysisTaskJetProtonCorr::Correlate(CorrType_t corr, Class_t cl, Ev_t ev,
					       TCollection *trgArray, TCollection *assArray, Float_t weight)
{
  AliHistCorr *histCorr = GetHistCorr(corr, cl, ev);

  TIter trgIter(trgArray);

  while (AliVParticle *trgPart = (AliVParticle*) trgIter()) {
    // count the trigger
    histCorr->Trigger(weight);

    // loop over associates
    TIter assIter(assArray);
    while (AliVParticle *assPart = (AliVParticle*) assIter()) {
      histCorr->Fill(trgPart, assPart, weight);
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskJetProtonCorr::Correlate(Trg_t trg, Ass_t ass, Class_t cl, Ev_t ev,
					       TCollection *trgArray, TCollection *assArray, Float_t weight)
{
  CorrType_t corr = (CorrType_t) (2 * trg + ass);
  AliHistCorr *histCorr = GetHistCorr(corr, cl, ev);

  TIter trgIter(trgArray);

  while (AliVParticle *trgPart = (AliVParticle*) trgIter()) {
    // count the trigger
    histCorr->Trigger(weight);

    // loop over associates
    TIter assIter(assArray);
    while (AliVParticle *assPart = (AliVParticle*) assIter()) {
      histCorr->Fill(trgPart, assPart, weight);
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskJetProtonCorr::CleanUpEvent()
{
  if (fAODEvent) {
    delete fPrimTrackArray;
    fPrimTrackArray = 0x0;
  }

  return kTRUE;
}

TObjArray* AliAnalysisTaskJetProtonCorr::CloneTracks(TObjArray* tracks) const
{
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  Int_t nTracks = tracks->GetEntriesFast();
  for (Int_t i = 0; i < nTracks; i++) {
    // tracksClone->Add(new AliDPhiBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge()));

    // WARNING: TObject::Clone() is very!!! expensive, unusable
    // tracksClone->Add(particle->Clone());

    if (AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (tracks->At(i)))
      tracksClone->Add(new AliESDtrack(*esdTrack));
    else if (AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*> (tracks->At(i)))
      tracksClone->Add(new AliAODTrack(*aodTrack));
  }

  return tracksClone;
}

AliAnalysisTaskJetProtonCorr::AliHistCorr::AliHistCorr(TString name, TList *outputList) :
  TNamed(name, name),
  fOutputList(outputList),
  fHistStat(0x0),
  fHistCorrPhi(0x0),
  fHistCorrPhi2(0x0),
  fHistCorrEtaPhi(0x0)
{
  // ctor

  fHistStat = new TH1F(Form("%s_stat", name.Data()), "statistics",
		       1, .5, 1.5);

  fHistCorrPhi = new TH1F(Form("%s_phi", name.Data()), ";#Delta #phi",
			  100, -2.*TMath::Pi(), 2.*TMath::Pi());
  fHistCorrPhi2 = new TH2F(Form("%s_phi2", name.Data()), ";#phi_{trg};#phi_{ass}",
			  100, 0.*TMath::Pi(), 2.*TMath::Pi(),
			  100, 0.*TMath::Pi(), 2.*TMath::Pi());
  fHistCorrEtaPhi = new TH2F(Form("%s_etaphi", name.Data()), ";#Delta#phi;#Delta#eta",
			     100, -1., 2*TMath::Pi()-1.,
			     100, -2., 2.);

  fOutputList->Add(fHistStat);
  fOutputList->Add(fHistCorrPhi);
  fOutputList->Add(fHistCorrPhi2);
  fOutputList->Add(fHistCorrEtaPhi);
}

AliAnalysisTaskJetProtonCorr::AliHistCorr::~AliHistCorr()
{
  // dtor
}

void AliAnalysisTaskJetProtonCorr::AliHistCorr::Fill(AliVParticle *trgPart, AliVParticle *assPart, Float_t weight)
{
  Float_t deltaEta = assPart->Eta() - trgPart->Eta();
  Float_t deltaPhi = assPart->Phi() - trgPart->Phi();
  if (deltaPhi > (2.*TMath::Pi()-1.))
    deltaPhi -= 2. * TMath::Pi();
  else if (deltaPhi < -1.)
    deltaPhi += 2. * TMath::Pi();
  // printf("trg: pt = %5.2f, phi = %5.2f, eta = %5.2f; ass: pt = %5.2f, phi = %5.2f, eta = %5.2f; deltaphi = %5.2f, deltaeta = %5.2f\n",
  // 	 trgPart->Pt(), trgPart->Phi(), trgPart->Eta(), assPart->Pt(), assPart->Phi(), assPart->Eta(), deltaPhi, deltaEta);

  fHistCorrPhi->Fill(deltaPhi);
  fHistCorrPhi2->Fill(trgPart->Phi(), assPart->Phi(), weight);
  fHistCorrEtaPhi->Fill(deltaPhi, deltaEta, weight);
}

// ----- histogram management -----
TH1* AliAnalysisTaskJetProtonCorr::AddHistogram(Hist_t hist, const char *hid, TString title,
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

TH2* AliAnalysisTaskJetProtonCorr::AddHistogram(Hist_t hist, const char *hid, TString title,
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

TH3* AliAnalysisTaskJetProtonCorr::AddHistogram(Hist_t hist, const char *hid, TString title,
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
