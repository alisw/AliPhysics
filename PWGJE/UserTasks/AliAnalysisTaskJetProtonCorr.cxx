// ROOT
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFormula.h"
#include "TRandom.h"

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

#include <iostream>
#include <cmath>

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
  fEventPlaneCheck(5.),
  fPrimTrackArray(0x0),
  fJetArray(0x0),
  fPoolMgr(),
  fPool(),
  fHistCorr(0x0),
  fOutputList(),
  fHist(),
  fShortTaskId("jet_prot_corr"),
  fUseStandardCuts(kTRUE),
  fCutsPrim(0x0),
  fCutsTwoTrackEff(0.02),
  fTrgPartPtMin(6.),
  fTrgPartPtMax(8.),
  fTrgJetPtMin(50.),
  fTrgJetPtMax(80.),
  fTrgJetLeadTrkPtMin(5.),
  fAssPartPtMin(2.),
  fAssPartPtMax(4.),
  fTrgAngleToEvPlane(TMath::Pi() / 4.)
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
  if (fUseStandardCuts) {
    fCutsPrim = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
  } else {
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

    fCutsPrim->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  }

  fCutsPrim->SetEtaRange(-0.9, 0.9);
  fCutsPrim->SetPtRange(0.15, 1E+15);

  // event mixing pool
  Double_t centralityBins[] = {
    0., 2., 4., 6., 8., 10., // central
    30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., // semi-central
    90.
  };
  Int_t nCentralityBins = sizeof(centralityBins)/sizeof(centralityBins[0]);

  Double_t vertexBins[] = {
    -10., -8., -6., -4., -2., 0., 2., 4., 6., 8., 10.
  };
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
				nVertexBins, vertexBins);
				// nPsiBins, psiBins);
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
  histStat->GetXaxis()->SetBinLabel(ID(kStatCentral));
  histStat->GetXaxis()->SetBinLabel(ID(kStatSemiCentral));

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


  AddHistogram(ID(kHistSignalTPC), "TPC dE/dx;p (GeV/c);dE/dx (arb. units)",
	       100, 0., 10., 200, 0., 300.);
  AddHistogram(ID(kHistSignalTOF), "TOF time of flight;p_{T} (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., 50.);
  AddHistogram(ID(kHistBetaTOF), "TOF beta;p (GeV/c); #beta",
	       100, 0., 10.,
	       100, 0., 1.);
  AddHistogram(ID(kHistDeltaTPC), "TPC dE/dx;p (GeV/c);dE/dx (arb. units)",
	       100, 0., 10., 200, -100., 100.);
  AddHistogram(ID(kHistDeltaTPCSemi), "TPC dE/dx;p (GeV/c);dE/dx (arb. units)",
	       100, 0., 10., 200, -100., 100.);
  AddHistogram(ID(kHistDeltaTOF), "TOF time of flight;p_{T} (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFSemi), "TOF time of flight;p_{T} (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);

  // Nsigma templates
  AddHistogram(ID(kHistNsigmaTPCe), "TPC N#sigma - e hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCmu), "TPC N#sigma - #mu hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCpi), "TPC N#sigma - #pi hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCk), "TPC N#sigma - K hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCp), "TPC N#sigma - p hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCd), "TPC N#sigma - d hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCe_e), "TPC N#sigma - e hypothesis (id. e);p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);

  AddHistogram(ID(kHistNsigmaTOFe), "TOF N#sigma - e hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFmu), "TOF N#sigma - #mu hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFpi), "TOF N#sigma - #pi hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFk), "TOF N#sigma - K hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFp), "TOF N#sigma - p hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFd), "TOF N#sigma - d hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFmismatch), "TOF N#sigma - mismatch;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFmismatch2), "TOF N#sigma - mismatch;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);

  // Nsigma templates
  AddHistogram(ID(kHistNsigmaTPCeSemi), "TPC N#sigma - e hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCmuSemi), "TPC N#sigma - #mu hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCpiSemi), "TPC N#sigma - #pi hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCkSemi), "TPC N#sigma - K hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCpSemi), "TPC N#sigma - p hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCdSemi), "TPC N#sigma - d hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);
  AddHistogram(ID(kHistNsigmaTPCe_eSemi), "TPC N#sigma - e hypothesis (id. e);p (GeV/c)",
	       100, 0., 10.,
	       100, -25., 25.);

  AddHistogram(ID(kHistNsigmaTOFeSemi), "TOF N#sigma - e hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFmuSemi), "TOF N#sigma - #mu hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFpiSemi), "TOF N#sigma - #pi hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFkSemi), "TOF N#sigma - K hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFpSemi), "TOF N#sigma - p hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFdSemi), "TOF N#sigma - d hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFmismatchSemi), "TOF N#sigma - mismatch;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTOFmismatch2Semi), "TOF N#sigma - mismatch;p (GeV/c)",
	       100, 0., 10.,
	       200, -50., 50.);

  // delta templates
  AddHistogram(ID(kHistDeltaTOFe), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFmu), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFpi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFk), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFp), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFd), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);

  AddHistogram(ID(kHistDeltaTOFeSemi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFmuSemi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFpiSemi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFkSemi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFpSemi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFdSemi), "TOF #Delta;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);

  // sigma comparisons
  AddHistogram(ID(kHistExpSigmaTOFe), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFmu), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFpi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFk), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFp), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFd), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);

  AddHistogram(ID(kHistExpSigmaTOFeSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFmuSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFpiSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFkSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFpSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);
  AddHistogram(ID(kHistExpSigmaTOFdSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, 0., .25);

  AddHistogram(ID(kHistCmpSigmaTOFe), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFmu), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFpi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFk), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFp), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFd), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);

  AddHistogram(ID(kHistCmpSigmaTOFeSemi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFmuSemi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFpiSemi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFkSemi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFpSemi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);
  AddHistogram(ID(kHistCmpSigmaTOFdSemi), "#sigma comparison;exp #sigma;template #sigma",
	       200, 0., .25, 200, 0., .25);

  // Nsigma distributions
  AddHistogram(ID(kHistNsigmaTPCTOF), "N#sigma TPC-TOF;p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFPt), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsed), "N#sigma TPC-TOF;p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedCentral), "N#sigma TPC-TOF (central);p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedSemiCentral), "N#sigma TPC-TOF (semi-central);p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedPt), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               50, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedPtCentral), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               50, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedPtSemiCentral), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               50, 0., 10.,
               100, -25., 25.,
               200, -50., 50.);

  AddHistogram(ID(kHistEvPlane), "default event plane;#Psi;counts",
               100, -0. * TMath::Pi(), 1. * TMath::Pi());
  AddHistogram(ID(kHistEvPlaneUsed), "default event plane;#Psi;counts",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);
  AddHistogram(ID(kHistEvPlaneCheck), "backup event plane;#Psi;counts",
               100, -1. * TMath::Pi(), 1. * TMath::Pi());
  AddHistogram(ID(kHistEvPlaneCheckUsed), "backup event plane;#Psi;counts",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);
  AddHistogram(ID(kHistEvPlaneCorr), "default - backup event plane;#Psi_{def};#Psi_{bak};counts",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);

  AddHistogram(ID(kHistJetPtCentral), "jet spectrum - central",
               40, 0., 200.);
  AddHistogram(ID(kHistJetPtSemi), "jet spectrum - semi-peripheral",
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
    if (IsClass(kClCentral))
      FillH1(kHistStat, kStatCentral);
    if (IsClass(kClSemiCentral))
      FillH1(kHistStat, kStatSemiCentral);

    FillH1(kHistEvPlane, fEventPlane);
    FillH1(kHistEvPlaneCheck, fEventPlaneCheck);
    for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
      if (IsClass((Class_t) iClass)) {
        FillH2(kHistCentralityUsed, fCentrality, iClass);
        FillH2(kHistCentralityCheckUsed, fCentralityCheck, iClass);
	FillH2(kHistEvPlaneUsed, fEventPlane, iClass);
	FillH2(kHistEvPlaneCheckUsed, fEventPlaneCheck, iClass);
	FillH3(kHistEvPlaneCorr, fEventPlane, fEventPlaneCheck, iClass);
      }
    }

    Bool_t corrEta  = fPIDResponse->UseTPCEtaCorrection();
    Bool_t corrMult = fPIDResponse->UseTPCMultiplicityCorrection();

    // select trigger particles and potential associated particles/protons
    TObjArray trgArray[kTrgLast];
    TObjArray assArray[kAssLast];

    TF1 fTOFsignal("fTOFsignal", &AliAnalysisTaskJetProtonCorr::TOFsignal, -2440., 2440., 4);
    fTOFsignal.SetParameter(0, 1.);
    fTOFsignal.SetParameter(1, 0.);
    fTOFsignal.SetParameter(2, 85.);
    fTOFsignal.SetParameter(3, 80.);

    Int_t nPrimTracks = fPrimTrackArray ? fPrimTrackArray->GetEntries() : 0;
    for (Int_t iTrack = 0; iTrack < nPrimTracks; ++iTrack) {
      AliVTrack *trk = (AliVTrack*) fPrimTrackArray->At(iTrack);
      FillH3(kHistNsigmaTPCTOF,
             trk->P(),
             fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
             fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
      FillH3(kHistNsigmaTPCTOFPt,
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
	FillH3(kHistNsigmaTPCTOFUsed,
	       trk->P(),
	       fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
	       fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
	FillH3(kHistNsigmaTPCTOFUsedPt,
	       trk->Pt(),
	       fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
	       fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
	if (IsClass(kClCentral)) {
	  FillH3(kHistNsigmaTPCTOFUsedCentral,
		 trk->P(),
		 fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
		 fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
	  FillH3(kHistNsigmaTPCTOFUsedPtCentral,
		 trk->Pt(),
		 fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
		 fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
	}
	if (IsClass(kClSemiCentral)) {
	  FillH3(kHistNsigmaTPCTOFUsedSemiCentral,
		 trk->P(),
		 fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
		 fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
	  FillH3(kHistNsigmaTPCTOFUsedPtSemiCentral,
		 trk->Pt(),
		 fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
		 fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
	}
	if (IsProton(trk)) {
	  assArray[kAssProt].Add(trk);
	  FillH1(kHistEtaPhiAssProt, trk->Phi(), trk->Eta());
	}

	// template generation
	Double_t deltaTPC;
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, trk, AliPID::kProton, deltaTPC);
	FillH2(kHistSignalTPC, trk->P(), trk->GetTPCsignal());
	FillH2(kHistDeltaTPC, trk->P(), deltaTPC);
	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trk) == AliPIDResponse::kDetPidOk) {
	  Double_t expTPC = fPIDResponse->GetTPCResponse().GetExpectedSignal(trk, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);
	  Double_t expSigmaTPC = fPIDResponse->GetTPCResponse().GetExpectedSigma(trk, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);

	  // loop over particles
	  for (Int_t iParticle = 0; iParticle <= AliPID::kDeuteron; ++iParticle) {
	    Double_t expTPCx = fPIDResponse->GetTPCResponse().GetExpectedSignal(trk, AliPID::EParticleType(iParticle), AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);
	    Double_t expSigmaTPCx = fPIDResponse->GetTPCResponse().GetExpectedSigma(trk, AliPID::EParticleType(iParticle), AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);
	    Double_t rndTPCx = gRandom->Gaus(expTPCx, expSigmaTPCx);
	    if(IsClass(kClCentral))
	      FillH2(kHistNsigmaTPCe, trk->P(), (rndTPCx-expTPC)/expSigmaTPC, 1., iParticle);
	    if(IsClass(kClSemiCentral))
	      FillH2(kHistNsigmaTPCeSemi, trk->P(), (rndTPCx-expTPC)/expSigmaTPC, 1., iParticle);
	  }
	}

	Double_t deltaTOF;
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, trk, AliPID::kProton, deltaTOF);
	FillH2(kHistSignalTOF, trk->Pt(), trk->GetTOFsignal() * 1.e-3); // ps -> ns
	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk) {
	  AliTOFPIDResponse &tofResponse = fPIDResponse->GetTOFResponse();

	  Float_t p = trk->GetInnerParam()->GetP();

	  Double_t expTOF = fPIDResponse->GetTOFResponse().GetExpectedSignal(trk, AliPID::kProton);
	  Double_t expSigmaTOF = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, expTOF, AliPID::kProton);
	  Double_t length = trk->GetIntegratedLength() * 1.e-2; // cm -> m
	  Double_t tof = trk->GetTOFsignal() * 1.e-12; // ps -> s
	  Double_t beta = length / tof / TMath::C();

	  FillH2(kHistBetaTOF, p, beta);

	  // intrinsic TOF smearing
	  Double_t signalSigma = TMath::Sqrt(fTOFsignal.Variance(-2000., 2000.));
	  Double_t signalSmear = fTOFsignal.GetRandom();
	  // t0 smearing
	  Float_t  timezeroSigma = tofResponse.GetT0binRes(tofResponse.GetMomBin(p));
	  Double_t timezeroSmear  = gRandom->Gaus(0., timezeroSigma);
	  // tracking smearing
	  Double_t fPar[] = { 0.008, 0.008, 0.002, 40.0 };

	  // loop over particles
	  for (Int_t iParticle = 0; iParticle <= AliPID::kDeuteron; ++iParticle) {
	    Double_t expTOFx = fPIDResponse->GetTOFResponse().GetExpectedSignal(trk, AliPID::EParticleType(iParticle));
	    Double_t cmpSigmaTOFx = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, expTOFx, AliPID::EParticleType(iParticle));

	    // tracking smearing
	    Double_t massx = AliPID::ParticleMassZ(AliPID::EParticleType(iParticle));
	    Double_t dppx = fPar[0] + fPar[1] * p + fPar[2] * massx / p;
	    Double_t expSigmaTOFx = dppx * expTOFx / (1.+ p * p / (massx * massx));
	    Double_t texpSmearx = gRandom->Gaus(0., TMath::Sqrt(expSigmaTOFx * expSigmaTOFx + fPar[3]*fPar[3]/p/p));
	    Double_t tmpSigmaTOFx = TMath::Sqrt(signalSigma*signalSigma +
						timezeroSigma*timezeroSigma +
						expSigmaTOFx*expSigmaTOFx + fPar[3]*fPar[3]/p/p);
	    // printf("sigma comparison %i, %f: %f, %f, %f -> %f vs %f\n",
	    // 	   iParticle, expTOFx,
	    // 	   signalSigma, // signal
	    // 	   timezeroSigma, // timezero
	    // 	   expSigmaTOFx, // tracking
	    // 	   tmpSigmaTOFx, // total
	    // 	   cmpSigmaTOFx); // from PID response

	    // TOF signal
	    Double_t rndTOFx = expTOFx + signalSmear + timezeroSmear + texpSmearx;

	    if (IsClass(kClCentral)) {
	      FillH2(kHistNsigmaTOFe, trk->P(), (rndTOFx-expTOF)/expSigmaTOF, 1., iParticle);
	      FillH2(kHistDeltaTOFe, trk->P(), (rndTOFx-expTOF) * 1.e-3, 1., iParticle);
	      FillH2(kHistExpSigmaTOFe, p, cmpSigmaTOFx * 1.e-3, 1., iParticle);
	      FillH2(kHistCmpSigmaTOFe, cmpSigmaTOFx * 1.e-3, tmpSigmaTOFx * 1.e-3, 1., iParticle);
	    }
	    if (IsClass(kClSemiCentral)) {
	      FillH2(kHistNsigmaTOFeSemi, trk->P(), (rndTOFx-expTOF)/expSigmaTOF, 1., iParticle);
	      FillH2(kHistDeltaTOFeSemi, trk->P(), (rndTOFx-expTOF) * 1.e-3, 1., iParticle);
	      FillH2(kHistCmpSigmaTOFeSemi, cmpSigmaTOFx * 1.e-3, tmpSigmaTOFx * 1.e-3, 1., iParticle);
	      FillH2(kHistExpSigmaTOFeSemi, p, cmpSigmaTOFx * 1.e-3, 1., iParticle);
	    }
	  }

	  Double_t rndTOFmismatch = AliTOFPIDResponse::GetMismatchRandomValue(trk->Eta());

	  if(IsClass(kClCentral)) {
	    FillH2(kHistNsigmaTOFmismatch, p, (rndTOFmismatch - expTOF) / expSigmaTOF);
	    FillH2(kHistNsigmaTOFmismatch2, p, (rndTOFmismatch - expTOF) / expSigmaTOF);
	    FillH2(kHistDeltaTOF, trk->P(), deltaTOF * 1.e-3); // ps -> ns
	  }
	  if(IsClass(kClSemiCentral)) {
	    FillH2(kHistNsigmaTOFmismatchSemi, p, (rndTOFmismatch - expTOF) / expSigmaTOF);
	    FillH2(kHistNsigmaTOFmismatch2Semi, p, (rndTOFmismatch - expTOF) / expSigmaTOF);
	    FillH2(kHistDeltaTOFSemi, trk->P(), deltaTOF * 1.e-3); // ps -> ns
	  }
	}
      }
    }

    // select trigger jet
    Int_t nJets = fJetArray ? fJetArray->GetEntries() : 0;
    for (Int_t iJet = 0; iJet < nJets; ++iJet) {
      AliAODJet *jet = (AliAODJet*) fJetArray->At(iJet);
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
	  // >= 0: don't require a trigger in the event
	  // >= 1: require a trigger in the event
	  if (trgArray[iTrg].GetEntries() >= 0) {
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

void AliAnalysisTaskJetProtonCorr::PrintTask(Option_t *option, Int_t indent) const
{
  AliAnalysisTaskSE::PrintTask(option, indent);

  std::cout << std::setw(indent) << " " << "using jet branch: " << fJetBranchName << std::endl;
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
  else
    eventGood = kFALSE;

  // retrieve event plane
  AliEventplane *eventPlane = InputEvent()->GetEventplane();
  if (eventPlane) {
    FillH1(kHistStat, kStatEvPlane);
    fEventPlane = eventPlane->GetEventplane("Q");
    fEventPlaneCheck = eventPlane->GetEventplane("V0", InputEvent());
    if (fEventPlaneCheck < 0)
      fEventPlaneCheck += TMath::Pi();
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
      // 4: track cuts esdTrackCutsH ???
      // 10: R_AA cuts
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
	    // mgr ? mgr->GetEventPool(fCentrality, fZvtx, fEventPlane) : 0x0;
	    mgr ? mgr->GetEventPool(fCentrality, fZvtx) : 0x0;
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
      if (AcceptAngleToEvPlane(trg->Phi(), GetEventPlane())) {
	// printf("track accepted with phi = %f, psi = %f\n", trg->Phi(), GetEventPlane());
        return kTRUE;
      }
    }
  }

  return kFALSE;
}

AliVTrack* AliAnalysisTaskJetProtonCorr::GetLeadingTrack(AliAODJet *jet) const
{
  // return leading track within a jet

  // check contributing tracks
  Int_t nJetTracks = jet->GetRefTracks()->GetEntriesFast();
  Int_t iLeadingTrack = -1;
  Float_t ptLeadingTrack = 0.;
  for (Int_t iTrack=0; iTrack < nJetTracks; ++iTrack) {
    AliAODTrack *track = (AliAODTrack*) jet->GetRefTracks()->At(iTrack);
    // find the leading track
    if (track->Pt() > ptLeadingTrack) {
      ptLeadingTrack = track->Pt();
      iLeadingTrack = iTrack;
    }
  }

  // retrieve the leading track
  return (AliAODTrack*) jet->GetRefTracks()->At(iLeadingTrack);
}

Bool_t AliAnalysisTaskJetProtonCorr::AcceptTrigger(AliAODJet *trg)
{
  // restrict eta
  if (TMath::Abs(trg->Eta()) > .5)
    return kFALSE;

  // require hard leading track
  // (biased jet sample)
  if (GetLeadingTrack(trg)->Pt() < fTrgJetLeadTrkPtMin)
    return kFALSE;

  // check for pt and azimuthal orientation
  if (IsSemiCentral()) {
    if (!AcceptAngleToEvPlane(trg->Phi(), GetEventPlane()))
      return kFALSE;
  }

  // printf("jet accepted with phi = %f, psi = %f\n", trg->Phi(), GetEventPlane());

  // spectrum for cross checks
  if (IsClass(kClCentral))
    FillH1(kHistJetPtCentral, trg->Pt());
  if (IsClass(kClSemiCentral))
    FillH1(kHistJetPtSemi, trg->Pt());

  if ((trg->Pt() < fTrgJetPtMin) || (trg->Pt() > fTrgJetPtMax))
    return kFALSE;

  return kTRUE;
}

Bool_t AliAnalysisTaskJetProtonCorr::AcceptAssoc(AliVTrack *trk)
{
  if ((trk->Pt() > fAssPartPtMin) && (trk->Pt() < fAssPartPtMax) &&
      (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trk) == AliPIDResponse::kDetPidOk) &&
      (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk))
    if ((trk->GetTPCsignalN() >= 60.) &&
	(trk->GetTPCsignal() >= 10) &&
	(TMath::Abs(trk->Eta()) < .9))
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

Bool_t AliAnalysisTaskJetProtonCorr::AcceptAngleToEvPlane(Float_t phi, Float_t psi)
{
  Float_t deltaPhi = phi - psi;

  // map to interval [-pi/2, pi/2)
  deltaPhi = std::fmod(deltaPhi + TMath::Pi()/2., TMath::Pi()) - TMath::Pi()/2.;
  // printf("delta phi = %f with phi = %f, psi = %f\n", deltaPhi, phi, psi);

  if (TMath::Abs(deltaPhi) < fTrgAngleToEvPlane)
    return kTRUE;
  else
    return kFALSE;
}

Float_t AliAnalysisTaskJetProtonCorr::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1,
						  Float_t phi2, Float_t pt2, Float_t charge2,
						  Float_t radius, Float_t bSign)
{
  // calculates dphistar
  // from AliUEHistograms

  Float_t dphistar = phi1 - phi2
    - charge1 * bSign * TMath::ASin(0.075 * radius / pt1)
    + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);

  static const Double_t kPi = TMath::Pi();

  // circularity
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;

  return dphistar;
}


Bool_t AliAnalysisTaskJetProtonCorr::AcceptTwoTracks(AliVParticle *trgPart, AliVParticle *assPart)
{
  // apply two track pair cut

  Float_t phi1 = trgPart->Phi();
  Float_t pt1 = trgPart->Pt();
  Float_t charge1 = trgPart->Charge();

  Float_t phi2 = assPart->Phi();
  Float_t pt2 = assPart->Pt();
  Float_t charge2 = assPart->Charge();

  // check ???
  Float_t bSign = 1.;
  Float_t deta = trgPart->Eta() - assPart->Eta();

  // optimization
  if (TMath::Abs(deta) < fCutsTwoTrackEff * 2.5 * 3) {
    // check first boundaries to see if is worth to loop and find the minimum
    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);

    const Float_t kLimit = fCutsTwoTrackEff * 3;

    Float_t dphistarminabs = 1e5;
    // Float_t dphistarmin = 1e5;
    if ((TMath::Abs(dphistar1) < kLimit) ||
	(TMath::Abs(dphistar2) < kLimit) ||
	(dphistar1 * dphistar2 < 0)) {
      for (Double_t rad=0.8; rad<2.51; rad+=0.01) {
	Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);
	Float_t dphistarabs = TMath::Abs(dphistar);

	if (dphistarabs < dphistarminabs) {
	  // dphistarmin = dphistar;
	  dphistarminabs = dphistarabs;
	}
      }

      if ((dphistarminabs < fCutsTwoTrackEff) &&
	  (TMath::Abs(deta) < fCutsTwoTrackEff))
	return kFALSE;
    }
  }

  return kTRUE;
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
      if (AcceptTwoTracks(trgPart, assPart))
	histCorr->Fill(trgPart, assPart, weight);
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskJetProtonCorr::Correlate(Trg_t trg, Ass_t ass, Class_t cl, Ev_t ev,
					       TCollection *trgArray, TCollection *assArray, Float_t weight)
{
  CorrType_t corr = (CorrType_t) (2 * trg + ass);

  return Correlate(corr, cl, ev, trgArray, assArray, weight);
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
