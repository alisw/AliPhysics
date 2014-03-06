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
#include "AliVTrack.h"
#include "AliVTrdTrack.h"
#include "AliVVertex.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliOADBContainer.h"
#include "AliTOFPIDParams.h"
#include "AliAnalysisTaskVnV0.h"

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
  fRunNumber(-1),
  fOADBContainerTOF(0x0),
  fParamsTOF(0x0),
  fEventplane(0x0),
  fTriggerMask(0),
  fClassMask(0),
  fCentrality(100.),
  fCentralityCheck(100.),
  fZvtx(0.),
  fPIDResponse(0x0),
  fEventPlaneAngle(5.),
  fEventPlaneAngleCheck(5.),
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
  fTrgAngleToEvPlane(TMath::Pi() / 4.),
  fTrgJetPhiModCent(new TF1("jetphimodcent", "1 + 2 * [0] * cos(2*x)", 0., 2 * TMath::Pi())),
  fTrgJetPhiModSemi(new TF1("jetphimodsemi", "1 + 2 * [0] * cos(2*x)", 0., 2 * TMath::Pi())),
  fTrgHadPhiModCent(new TF1("hadphimodcent", "1 + 2 * [0] * cos(2*x)", 0., 2 * TMath::Pi())),
  fTrgHadPhiModSemi(new TF1("hadphimodsemi", "1 + 2 * [0] * cos(2*x)", 0., 2 * TMath::Pi()))
{
  // default ctor

  fkCorrTypeName[kCorrHadHad]  = "hh";
  fkCorrTypeName[kCorrHadProt] = "hp";
  fkCorrTypeName[kCorrJetHad]  = "jh";
  fkCorrTypeName[kCorrJetProt] = "jp";
  fkCorrTypeName[kCorrRndJetHad]  = "rjh";
  fkCorrTypeName[kCorrRndJetProt]  = "rjp";
  fkCorrTypeName[kCorrRndHadHad]  = "rhh";
  fkCorrTypeName[kCorrRndHadProt]  = "rhp";
  fkCorrTypeName[kCorrRndJetExcHad]  = "rjeh";
  fkCorrTypeName[kCorrRndJetExcProt]  = "rjep";
  fkCorrTypeName[kCorrRndHadExcHad]  = "rheh";
  fkCorrTypeName[kCorrRndHadExcProt]  = "rhep";

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

  fTrgJetPhiModCent->SetParameter(0, .05);
  fTrgJetPhiModSemi->SetParameter(0, .20);
  fTrgHadPhiModCent->SetParameter(0, .05);
  fTrgHadPhiModSemi->SetParameter(0, .20);

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

  // open OADB file and retrieve the OADB container with TOF parameters
  TString oadbFileName =
    TString::Format("%s/COMMON/PID/data/TOFPIDParams.root", AliAnalysisManager::GetOADBPath());
  TFile *oadbFile = TFile::Open(oadbFileName); 
  if(!oadbFile->IsOpen())
    AliFatal(Form("Cannot open OADB file %s", oadbFileName.Data()));
  AliOADBContainer *oadbContainer =
    (AliOADBContainer*) oadbFile->Get("TOFoadb");
  if (!oadbContainer)
    AliFatal("Cannot fetch OADB container for VZERO EP selection");
  fOADBContainerTOF = new AliOADBContainer(*oadbContainer);
  oadbFile->Close();
  delete oadbFile;

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
  histStat->GetXaxis()->SetBinLabel(ID(kStatVtx));
  histStat->GetXaxis()->SetBinLabel(ID(kStatEvCuts));
  histStat->GetXaxis()->SetBinLabel(ID(kStatUsed));
  histStat->GetXaxis()->SetBinLabel(ID(kStatCent));
  histStat->GetXaxis()->SetBinLabel(ID(kStatEvPlane));
  histStat->GetXaxis()->SetBinLabel(ID(kStatPID));
  histStat->GetXaxis()->SetBinLabel(ID(kStatCentral));
  histStat->GetXaxis()->SetBinLabel(ID(kStatSemiCentral));

  AddHistogram(ID(kHistVertexNctb), "number of vertex contributors;N_{ctb};counts",
	       100, 0., 2000.);
  AddHistogram(ID(kHistVertexZ), "z-position of primary vertex;z (cm);counts",
	       100, -50., 50.);

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
  AddHistogram(ID(kHistCentralityVsMult), "centrality - multiplicity;centrality percentile (%);N_{prim}",
	       100, 0., 100.,
	       100, 0., 2500.);

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
  AddHistogram(ID(kHistDeltaTOF), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);
  AddHistogram(ID(kHistDeltaTOFSemi), "TOF time of flight;p (GeV/c);t (ns)",
	       100, 0., 10., 200, -2., 2.);

  // Nsigma templates - central
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
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFmu), "TOF N#sigma - #mu hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFpi), "TOF N#sigma - #pi hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFk), "TOF N#sigma - K hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFp), "TOF N#sigma - p hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFd), "TOF N#sigma - d hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFmismatch), "TOF N#sigma - mismatch;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);

  // Nsigma templates - semi-central
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
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFmuSemi), "TOF N#sigma - #mu hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFpiSemi), "TOF N#sigma - #pi hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFkSemi), "TOF N#sigma - K hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFpSemi), "TOF N#sigma - p hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFdSemi), "TOF N#sigma - d hypothesis;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTOFmismatchSemi), "TOF N#sigma - mismatch;p (GeV/c)",
	       100, 0., 10.,
	       200, -100., 100.);

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
  // AddHistogram(ID(kHistExpSigmaTOFe), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFmu), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFpi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFk), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFp), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFd), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);

  // AddHistogram(ID(kHistExpSigmaTOFeSemi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFmuSemi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFpiSemi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFkSemi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFpSemi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);
  // AddHistogram(ID(kHistExpSigmaTOFdSemi), "TOF time of flight;p (GeV/c);t (ns)",
  // 	       100, 0., 10., 200, 0., .25);

  // AddHistogram(ID(kHistCmpSigmaTOFe), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFmu), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFpi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFk), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFp), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFd), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);

  // AddHistogram(ID(kHistCmpSigmaTOFeSemi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFmuSemi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFpiSemi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFkSemi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFpSemi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);
  // AddHistogram(ID(kHistCmpSigmaTOFdSemi), "#sigma comparison;exp #sigma;template #sigma",
  // 	       200, 0., .25, 200, 0., .25);

  // Nsigma distributions
  AddHistogram(ID(kHistNsigmaTPCTOF), "N#sigma TPC-TOF;p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFPt), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsed), "N#sigma TPC-TOF;p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedCentral), "N#sigma TPC-TOF (central);p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedSemiCentral), "N#sigma TPC-TOF (semi-central);p (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               100, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedPt), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               50, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedPtCentral), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               50, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);
  AddHistogram(ID(kHistNsigmaTPCTOFUsedPtSemiCentral), "N#sigma TPC-TOF;p_{T} (GeV/c);N#sigma_{TPC};N#sigma_{TOF}",
               50, 0., 10.,
               100, -25., 25.,
               200, -100., 100.);

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
  AddHistogram(ID(kHistEvPlaneCorr), "default - backup event plane;#Psi_{def};#Psi_{bak};event class",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);
  AddHistogram(ID(kHistEvPlaneCorrNoTrgJets), "event plane w/ and w/o trigger jets;#Psi_{full};#Psi_{notrgjets};event class",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);
  AddHistogram(ID(kHistEvPlaneCorrNoTrgJetsTrgd), "event plane w/ and w/o trigger jets;#Psi_{full};#Psi_{notrgjets};event class",
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
               100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       kClLast, -.5, kClLast-.5);

  AddHistogram(ID(kHistJetPtCentral), "jet spectrum - central;p_{T}^{jet,ch};counts",
               40, 0., 200.);
  AddHistogram(ID(kHistJetPtSemi), "jet spectrum - semi-peripheral;p_{T}^{jet,ch};counts",
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

  AddHistogram(ID(kHistPhiTrgJetEvPlane), "trg jet;#varphi - #Psi_{ev};centrality",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100., 0., 100.);
  AddHistogram(ID(kHistPhiTrgHadEvPlane), "trg had;#varphi - #Psi_{ev};centrality",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100., 0., 100.);
  AddHistogram(ID(kHistPhiAssHadEvPlane), "ass had;#varphi - #Psi_{ev};centrality",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100., 0., 100.);
  AddHistogram(ID(kHistPhiAssProtEvPlane), "ass prot;#varphi - #Psi_{ev};centrality",
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100., 0., 100.);
  AddHistogram(ID(kHistPhiAssHadVsEvPlane), "ass had;#Psi_{ev};#varphi;centrality",
	       100, -0. * TMath::Pi(), 1. * TMath::Pi(),
	       100, -0. * TMath::Pi(), 2. * TMath::Pi(),
	       100., 0., 100.);

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

void AliAnalysisTaskJetProtonCorr::SetParamsTOF()
{
  fParamsTOF =
    dynamic_cast<AliTOFPIDParams*> (fOADBContainerTOF->GetObject(fRunNumber, "TOFparams"));

  if (!fParamsTOF)
    AliError(Form("failed to load TOF parameters for run %i", fRunNumber));
  else if (fDebug > 2) {
    printf("loaded TOF parameters for run %i\n",
	   fRunNumber);
    printf("   intrinsic resolution: %f ps\n",
	   fParamsTOF->GetTOFresolution());
    printf("   tail fraction: %f\n",
	   fParamsTOF->GetTOFtail());
    printf("   start time method: %i\n",
	   fParamsTOF->GetStartTimeMethod());
    printf("   TOF signal parametrization: %f, %f, %f, %f\n",
	   fParamsTOF->GetSigParams(0), fParamsTOF->GetSigParams(1),
	   fParamsTOF->GetSigParams(2), fParamsTOF->GetSigParams(3));
    printf("   matching loss MC: %f%%\n",
	   fParamsTOF->GetTOFmatchingLossMC());
    printf("   additional mismatch for MC: %f%%\n",
	   fParamsTOF->GetTOFadditionalMismForMC());
    printf("   time offset: %f\n",
	   fParamsTOF->GetTOFtimeOffset());
  }
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

    FillH1(kHistVertexZ, fZvtx);

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

    FillH1(kHistEvPlane, fEventPlaneAngle);
    FillH1(kHistEvPlaneCheck, fEventPlaneAngleCheck);
    for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
      if (IsClass((Class_t) iClass)) {
        FillH2(kHistCentralityUsed, fCentrality, iClass);
        FillH2(kHistCentralityCheckUsed, fCentralityCheck, iClass);
	FillH2(kHistEvPlaneUsed, fEventPlaneAngle, iClass);
	FillH2(kHistEvPlaneCheckUsed, fEventPlaneAngleCheck, iClass);
	FillH3(kHistEvPlaneCorr, fEventPlaneAngle, fEventPlaneAngleCheck, iClass);
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
    if (fParamsTOF) {
      Float_t res  = fParamsTOF->GetTOFresolution();
      Float_t tail = fParamsTOF->GetTOFtail() * res;
      fTOFsignal.SetParameter(2, res);
      fTOFsignal.SetParameter(3, tail);
    }
    else {
      fTOFsignal.SetParameter(2, 85.);
      fTOFsignal.SetParameter(3, 80.);
    }

    Int_t nPrimTracks = fPrimTrackArray ? fPrimTrackArray->GetEntries() : 0;
    FillH2(kHistCentralityVsMult, fCentrality, nPrimTracks);
    for (Int_t iTrack = 0; iTrack < nPrimTracks; ++iTrack) {
      AliVTrack *trk = (AliVTrack*) fPrimTrackArray->At(iTrack);
      FillH2(kHistSignalTPC, trk->P(), trk->GetTPCsignal());
      // ??? pt or p?
      FillH2(kHistSignalTOF, trk->P(), trk->GetTOFsignal() * 1.e-3); // ps -> ns
      FillH3(kHistNsigmaTPCTOF,
             trk->P(),
             fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
             fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));
      FillH3(kHistNsigmaTPCTOFPt,
             trk->Pt(),
             fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton),
             fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton));

      FillH3(kHistPhiAssHadVsEvPlane, fEventPlaneAngle, trk->Phi(), fCentrality);
      Float_t phiRel = trk->Phi() - fEventPlaneAngle;
      if (gRandom->Rndm() > .5)
	phiRel -= TMath::Pi();
      if (phiRel < 0.)
	phiRel += 2. * TMath::Pi();
      if (phiRel < 0.)
	AliError(Form("phiRel = %f less than zero, from phi = %f, psi = %f",
		      phiRel, trk->Phi(), fEventPlaneAngle));
      else if (phiRel > 2*TMath::Pi())
	AliError(Form("phiRel = %f greater than 2pi, from phi = %f, psi = %f",
		      phiRel, trk->Phi(), fEventPlaneAngle));

      if (AcceptTrigger(trk)) {
	trgArray[kTrgHad].Add(trk);
	FillH1(kHistEtaPhiTrgHad, trk->Phi(), trk->Eta());
      }
      if (AcceptAssoc(trk)) {
	assArray[kAssHad].Add(trk);
	FillH1(kHistEtaPhiAssHad, trk->Phi(), trk->Eta());
	FillH2(kHistPhiAssHadEvPlane, phiRel, fCentrality);
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
	  if (IsProton(trk))
	    FillH1(kHistEtaPhiAssProt, trk->Phi(), trk->Eta());
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
	  FillH2(kHistPhiAssProtEvPlane, phiRel, fCentrality);
	}

	// template generation
	Double_t deltaTPC;
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, trk, AliPID::kProton, deltaTPC);
	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trk) == AliPIDResponse::kDetPidOk) {
	  Double_t expTPC = fPIDResponse->GetTPCResponse().GetExpectedSignal(trk, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);
	  Double_t expSigmaTPC = fPIDResponse->GetTPCResponse().GetExpectedSigma(trk, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);

	  // loop over particles
	  for (Int_t iParticle = 0; iParticle <= AliPID::kDeuteron; ++iParticle) {
	    Double_t expTPCx = fPIDResponse->GetTPCResponse().GetExpectedSignal(trk, AliPID::EParticleType(iParticle), AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);
	    Double_t expSigmaTPCx = fPIDResponse->GetTPCResponse().GetExpectedSigma(trk, AliPID::EParticleType(iParticle), AliTPCPIDResponse::kdEdxDefault, corrEta, corrMult);
	    Double_t rndTPCx = gRandom->Gaus(expTPCx, expSigmaTPCx);

	    if(IsClass(kClCentral)) {
	      FillH2(kHistNsigmaTPCe, trk->P(), (rndTPCx-expTPC)/expSigmaTPC, 1., iParticle);
	      FillH2(kHistDeltaTPC, trk->P(), deltaTPC);
	    }
	    if(IsClass(kClSemiCentral)) {
	      FillH2(kHistNsigmaTPCeSemi, trk->P(), (rndTPCx-expTPC)/expSigmaTPC, 1., iParticle);
	      FillH2(kHistDeltaTPCSemi, trk->P(), deltaTPC);
	    }
	  }
	}

	Double_t deltaTOF;
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, trk, AliPID::kProton, deltaTOF);
	if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk) {
	  AliTOFPIDResponse &tofResponse = fPIDResponse->GetTOFResponse();

	  Float_t p = trk->P();
	  if (const AliExternalTrackParam *param = trk->GetInnerParam())
	    p = param->GetP();

	  Double_t expTOF = fPIDResponse->GetTOFResponse().GetExpectedSignal(trk, AliPID::kProton);
	  Double_t expSigmaTOF = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, expTOF, AliPID::kProton);
	  Double_t length = trk->GetIntegratedLength() * 1.e-2; // cm -> m
	  Double_t tof = trk->GetTOFsignal() * 1.e-12; // ps -> s
	  Double_t beta = length / tof / TMath::C();

	  FillH2(kHistBetaTOF, p, beta);

	  // intrinsic TOF smearing
	  // Double_t signalSigma = TMath::Sqrt(fTOFsignal.Variance(-2000., 2000.));
	  Double_t signalSmear = fTOFsignal.GetRandom();
	  // t0 smearing
	  Float_t  timezeroSigma = tofResponse.GetT0binRes(tofResponse.GetMomBin(p));
	  Double_t timezeroSmear  = gRandom->Gaus(0., timezeroSigma);
	  // tracking smearing (default parameters, should be overwritten from OADB)
	  Double_t fPar[] = { 0.008, 0.008, 0.002, 40.0 };
	  if (fParamsTOF)
	    for (Int_t i = 0; i < 4; ++i)
	      fPar[i] = fParamsTOF->GetSigParams(i);

	  // loop over particles
	  for (Int_t iParticle = 0; iParticle <= AliPID::kDeuteron; ++iParticle) {
	    Double_t expTOFx = fPIDResponse->GetTOFResponse().GetExpectedSignal(trk, AliPID::EParticleType(iParticle));
	    // Double_t cmpSigmaTOFx = fPIDResponse->GetTOFResponse().GetExpectedSigma(p, expTOFx, AliPID::EParticleType(iParticle));

	    // tracking smearing
	    Double_t massx = AliPID::ParticleMassZ(AliPID::EParticleType(iParticle));
	    Double_t dppx = fPar[0] + fPar[1] * p + fPar[2] * massx / p;
	    Double_t expSigmaTOFx = dppx * expTOFx / (1.+ p * p / (massx * massx));
	    Double_t texpSmearx = gRandom->Gaus(0., TMath::Sqrt(expSigmaTOFx * expSigmaTOFx + fPar[3]*fPar[3]/p/p));
	    // Double_t tmpSigmaTOFx = TMath::Sqrt(signalSigma*signalSigma +
	    // 					timezeroSigma*timezeroSigma +
	    // 					expSigmaTOFx*expSigmaTOFx + fPar[3]*fPar[3]/p/p);
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
	      // FillH2(kHistExpSigmaTOFe, p, cmpSigmaTOFx * 1.e-3, 1., iParticle);
	      // FillH2(kHistCmpSigmaTOFe, cmpSigmaTOFx * 1.e-3, tmpSigmaTOFx * 1.e-3, 1., iParticle);
	    }
	    if (IsClass(kClSemiCentral)) {
	      FillH2(kHistNsigmaTOFeSemi, trk->P(), (rndTOFx-expTOF)/expSigmaTOF, 1., iParticle);
	      FillH2(kHistDeltaTOFeSemi, trk->P(), (rndTOFx-expTOF) * 1.e-3, 1., iParticle);
	      // FillH2(kHistCmpSigmaTOFeSemi, cmpSigmaTOFx * 1.e-3, tmpSigmaTOFx * 1.e-3, 1., iParticle);
	      // FillH2(kHistExpSigmaTOFeSemi, p, cmpSigmaTOFx * 1.e-3, 1., iParticle);
	    }
	  }

	  Double_t rndTOFmismatch = AliTOFPIDResponse::GetMismatchRandomValue(trk->Eta());

	  if(IsClass(kClCentral)) {
	    FillH2(kHistNsigmaTOFmismatch, p, (rndTOFmismatch - expTOF) / expSigmaTOF);
	    FillH2(kHistDeltaTOF, trk->P(), deltaTOF * 1.e-3); // ps -> ns
	  }
	  if(IsClass(kClSemiCentral)) {
	    FillH2(kHistNsigmaTOFmismatchSemi, p, (rndTOFmismatch - expTOF) / expSigmaTOF);
	    FillH2(kHistDeltaTOFSemi, trk->P(), deltaTOF * 1.e-3); // ps -> ns
	  }
	}
      }
    }

    // select trigger jet
    // and remove them from the Q vector
    Int_t nJets = fJetArray ? fJetArray->GetEntries() : 0;
    TVector2 *qVector  = fEventplane->GetQVector();

    for (Int_t iJet = 0; iJet < nJets; ++iJet) {
      AliAODJet *jet = (AliAODJet*) fJetArray->At(iJet);

      if (AcceptTrigger(jet)) {
	trgArray[kTrgJet].Add(jet);
	FillH1(kHistEtaPhiTrgJet, jet->Phi(), jet->Eta());

	if (qVector) {
	  Int_t nRefTracks = jet->GetRefTracks()->GetEntriesFast();
	  for (Int_t iTrack = 0; iTrack < nRefTracks; ++iTrack) {
	    AliVTrack *track = (AliVTrack*) jet->GetRefTracks()->At(iTrack);
	    
	    if (fEventplane) {
	      TVector2 evplaneContrib(fEventplane->GetQContributionX(track),
				      fEventplane->GetQContributionY(track));
	      *qVector -= evplaneContrib;
	    }
	  }
	}
      }
    }
    // printf("event plane angle before/after removal of trigger jets: %f/%f, diff: %f\n",
    // 	   fEventPlaneAngle, qVector->Phi()/2., qVector->Phi()/2. - fEventPlaneAngle);
    if (qVector) {
      for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
	if (IsClass((Class_t) iClass)) {
	  FillH3(kHistEvPlaneCorrNoTrgJets, fEventPlaneAngle, qVector->Phi()/2., iClass);
	  if (trgArray[kTrgJet].GetEntriesFast() > 0)
	    FillH3(kHistEvPlaneCorrNoTrgJetsTrgd, fEventPlaneAngle, qVector->Phi()/2., iClass);
	}
      }
    }

    // invent a trigger jet/hadron and correlated associates
    Float_t pFraction = assArray[kAssHad].GetEntries() > 0 ?
      (Float_t) (assArray[kAssProt].GetEntries()) / assArray[kAssHad].GetEntries() :
      .5;
    GenerateRandom(&trgArray[kTrgJetRnd], &trgArray[kTrgHadRnd],
		   &assArray[kAssHadJetExc], &assArray[kAssProtJetExc],
		   &assArray[kAssHadHadExc], &assArray[kAssProtHadExc],
		   pFraction);

    // correlate, both same and mixed event
    for (Int_t iClass = 0; iClass < kClLast; ++iClass) {
      if (IsClass((Class_t) iClass)) {
	for (Int_t iTrg = 0; iTrg < kTrgLast; ++iTrg) {
	  for (Int_t iAss = 0; iAss <= kAssProt; ++iAss) {
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
	// correlate artificial triggers and associates
	Correlate(kCorrRndJetExcHad,  (Class_t) iClass, kEvSame, &trgArray[kTrgJetRnd], &assArray[kAssHadJetExc]);
	Correlate(kCorrRndJetExcProt, (Class_t) iClass, kEvSame, &trgArray[kTrgJetRnd], &assArray[kAssProtJetExc]);
	Correlate(kCorrRndHadExcHad,  (Class_t) iClass, kEvSame, &trgArray[kTrgHadRnd], &assArray[kAssHadHadExc]);
	Correlate(kCorrRndHadExcProt, (Class_t) iClass, kEvSame, &trgArray[kTrgHadRnd], &assArray[kAssProtHadExc]);
      }
    }

    trgArray[kTrgJetRnd].Delete();
    trgArray[kTrgHadRnd].Clear();
    assArray[kAssHadJetExc].Delete();
    assArray[kAssProtJetExc].Clear();
    assArray[kAssHadHadExc].Delete();
    assArray[kAssProtHadExc].Clear();
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

  AliVEvent::EOfflineTriggerTypes physSel =
    (AliVEvent::EOfflineTriggerTypes) ((AliInputEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
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

  // check for run change
  if (fRunNumber != InputEvent()->GetRunNumber()) {
    fRunNumber = InputEvent()->GetRunNumber();
    SetParamsTOF();
  }

  // retrieve z-vertex position
  fZvtx = 100.;
  const AliVVertex *vtx = InputEvent()->GetPrimaryVertex();
  if (vtx) {
    FillH1(kHistStat, kStatVtx);
    FillH1(kHistVertexNctb, vtx->GetNContributors());
    if (vtx->GetNContributors() >= 3.)
      fZvtx = vtx->GetZ();
  }

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
  fEventplane = InputEvent()->GetEventplane();
  if (fEventplane) {
    fEventPlaneAngle = fEventplane->GetEventplane("Q");
    fEventPlaneAngleCheck = fEventplane->GetEventplane("V0", InputEvent());
    // use V0 event plane angle from flow task:
    // fEventPlaneAngleCheck = AliAnalysisTaskVnV0::GetPsi2V0A();
    // printf("V0A evplane = %f\n", fEventPlaneAngleCheck);
    // fEventPlaneAngleCheck = AliAnalysisTaskVnV0::GetPsi2V0C();
    // printf("V0C evplane = %f\n", fEventPlaneAngleCheck);

    FillH1(kHistStat, kStatEvPlane);
    if (fEventPlaneAngleCheck < 0)
      fEventPlaneAngleCheck += TMath::Pi();
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
	    // mgr ? mgr->GetEventPool(fCentrality, fZvtx, fEventPlaneAngle) : 0x0;
	    mgr ? mgr->GetEventPool(fCentrality, fZvtx) : 0x0;
	}
      }
    }
  }

  return eventGood;
}

Bool_t AliAnalysisTaskJetProtonCorr::AcceptTrigger(AliVTrack *trg)
{
  // check pt interval
  Bool_t acceptPt =
    (trg->Pt() >= fTrgPartPtMin) &&
    (trg->Pt() <= fTrgPartPtMax);

  // for semi-central events check phi w.r.t. event plane
  Bool_t acceptOrientation = IsSemiCentral() ?
    AcceptAngleToEvPlane(trg->Phi(), GetEventPlaneAngle()) :
    kTRUE;

  if (acceptPt) {
    Float_t phiRel = trg->Phi() - fEventPlaneAngle;
    if (gRandom->Rndm() > .5)
      phiRel -= TMath::Pi();
    if (phiRel < 0.)
      phiRel += 2. * TMath::Pi();
    FillH2(kHistPhiTrgHadEvPlane, phiRel, fCentrality);
  }

  return (acceptPt && acceptOrientation);
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
  // (leading track biased jet sample)
  if (GetLeadingTrack(trg)->Pt() < fTrgJetLeadTrkPtMin)
    return kFALSE;

  // check for jet orientation w.r.t. event plane
  // (constrained only for semi-central events)
  Bool_t acceptOrientation = IsSemiCentral() ?
    AcceptAngleToEvPlane(trg->Phi(), GetEventPlaneAngle()) :
    kTRUE;
  // check for jet pt
  Bool_t acceptPt =
    (trg->Pt() >= fTrgJetPtMin) &&
    (trg->Pt() <= fTrgJetPtMax);

  if (acceptPt) {
    // store the azimuthal distribution relative to the event plane
    // for the jets in the pt range of interest
    Float_t phiRel = trg->Phi() - fEventPlaneAngle;
    if (gRandom->Rndm() > .5)
      phiRel -= TMath::Pi();
    if (phiRel < 0.)
      phiRel += 2. * TMath::Pi();
    FillH2(kHistPhiTrgJetEvPlane, phiRel, fCentrality);
  }

  if (acceptOrientation) {
    // spectrum for cross checks
    if (IsClass(kClCentral))
      FillH1(kHistJetPtCentral, trg->Pt());
    if (IsClass(kClSemiCentral))
      FillH1(kHistJetPtSemi, trg->Pt());
  }

  return (acceptPt && acceptOrientation);
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
  if (!histCorr) {
    AliError(Form("no correlation histograms for corr %i, cl %i, ev %i",
		  corr, cl, ev));
    return kFALSE;
  }

  TIter assIter(assArray);
  while (TObject *assObj = assIter()) {
    if (dynamic_cast<AliVParticle*> (assObj)) {
      AliVParticle *assPart = (AliVParticle*) assObj;
      histCorr->Ass(assPart->Phi(), assPart->Eta(), weight);
    }
    else if (dynamic_cast<TLorentzVector*> (assObj)) {
      TLorentzVector *assVec = (TLorentzVector*) assObj;
      histCorr->Ass(assVec->Phi(), assVec->Eta(), weight);
    }
  }

  TIter trgIter(trgArray);
  while (TObject *trgObj = trgIter()) {
    if (AliVParticle *trgPart = dynamic_cast<AliVParticle*> (trgObj)) {
      // count the trigger
      histCorr->Trigger(trgPart->Phi(), trgPart->Eta(), weight);

      // loop over associates
      assIter.Reset();
      while (TObject *assObj = assIter()) {
	if (AliVParticle *assPart = dynamic_cast<AliVParticle*> (assObj)) {
	  if (AcceptTwoTracks(trgPart, assPart))
	    histCorr->Fill(trgPart, assPart, weight);
	}
	else if (TLorentzVector *assVec = dynamic_cast<TLorentzVector*> (assObj)) {
	  AliFatal(Form("got %p, but not implemented", assVec));
	}
      }
    }
    else if (TLorentzVector *trgVec = dynamic_cast<TLorentzVector*> (trgObj)) {
      // count the trigger
      histCorr->Trigger(trgVec->Phi(), trgVec->Eta(), weight);

      // loop over associates
      assIter.Reset();
      while (TObject *assObj = assIter()) {
	if (AliVParticle *assPart = dynamic_cast<AliVParticle*> (assObj)) {
	  histCorr->Fill(trgVec, assPart, weight);
	}
	else if (TLorentzVector *assVec = dynamic_cast<TLorentzVector*> (assObj)) {
	  histCorr->Fill(trgVec, assVec, weight);
	}
      }
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

Bool_t AliAnalysisTaskJetProtonCorr::GenerateRandom(TCollection *trgJetArray,
						    TCollection *trgHadArray,
						    TCollection *assHadJetArray,
						    TCollection *assProtJetArray,
						    TCollection *assHadHadArray,
						    TCollection *assProtHadArray,
						    Float_t pFraction)
{
  // generate random direction

  Float_t meanNoPart = 10;
  Int_t nPart = gRandom->Poisson(meanNoPart);

  Float_t trgJetEta = gRandom->Uniform(-.5, .5);
  Float_t trgJetPhi = 0.;
  if (IsClass(kClSemiCentral)) {
    do {
      trgJetPhi = fTrgJetPhiModSemi->GetRandom();
    } while (!AcceptAngleToEvPlane(trgJetPhi, 0));
    trgJetPhi += fEventPlaneAngle;
    if (trgJetPhi < 0.)
      trgJetPhi += 2. * TMath::Pi();
    else if (trgJetPhi > 2. * TMath::Pi())
      trgJetPhi -= 2. * TMath::Pi();
  }
  else {
    trgJetPhi = fTrgJetPhiModCent->GetRandom();
  }

  // generate one trigger particle
  TLorentzVector *trgJet = new TLorentzVector();
  trgJet->SetPtEtaPhiM(50., trgJetEta, trgJetPhi, .14);
  trgJetArray->Add(trgJet);

  // generate direction for away side
  Float_t trgJetEtaAway = gRandom->Uniform(-.9, .9);
  Float_t trgJetPhiAway = trgJetPhi + TMath::Pi() + gRandom->Gaus(0., .2);

  // generate associated particles
  // with proton/hadron ratio observed in this event
  for (Int_t iPart = 0; iPart < nPart; ++iPart) {
    Float_t deltaEta, deltaPhi;
    const Float_t r = .8;
    do {
      deltaEta = gRandom->Uniform(-r, r);
      deltaPhi = gRandom->Uniform(-r, r);
    } while ((deltaEta * deltaEta + deltaPhi * deltaPhi) > (r * r));

    TLorentzVector *assPart = new TLorentzVector();
    Float_t eta = trgJetEtaAway + deltaEta;
    Float_t phi = trgJetPhiAway + deltaPhi;
    if (eta > .9) {
      delete assPart;
      continue;
    }
    if (phi < 0.)
      phi += 2. * TMath::Pi();
    else if (phi > 2. * TMath::Pi())
      phi -= 2. * TMath::Pi();
    assPart->SetPtEtaPhiM(3., eta, phi, .14);
    assHadJetArray->Add(assPart);
    if (gRandom->Uniform() < pFraction)
      assProtJetArray->Add(assPart);
  }

  // trigger hadron
  Float_t trgHadEta = gRandom->Uniform(-.9, .9);
  Float_t trgHadPhi = 0.;
  if (IsClass(kClSemiCentral)) {
    do {
      trgHadPhi = fTrgHadPhiModSemi->GetRandom();
    } while (!AcceptAngleToEvPlane(trgHadPhi, 0));
    trgHadPhi += fEventPlaneAngle;
    if (trgHadPhi < 0.)
      trgHadPhi += 2. * TMath::Pi();
    else if (trgHadPhi > 2. * TMath::Pi())
      trgHadPhi -= 2. * TMath::Pi();
  }
  else {
    trgHadPhi = fTrgHadPhiModCent->GetRandom();
  }

  // generate one trigger particle
  TLorentzVector *trgHad = new TLorentzVector();
  trgHad->SetPtEtaPhiM(50., trgHadEta, trgHadPhi, .14);
  trgHadArray->Add(trgHad);

  // generate direction for away side
  Float_t trgHadEtaAway = gRandom->Uniform(-.9, .9);
  Float_t trgHadPhiAway = trgHadPhi + TMath::Pi() + gRandom->Gaus(0., .2);

  // generate associated particles
  // with proton/hadron ratio observed in this event
  for (Int_t iPart = 0; iPart < nPart; ++iPart) {
    Float_t deltaEta, deltaPhi;
    const Float_t r = .8;
    do {
      deltaEta = gRandom->Uniform(-r, r);
      deltaPhi = gRandom->Uniform(-r, r);
    } while ((deltaEta * deltaEta + deltaPhi * deltaPhi) > (r * r));

    TLorentzVector *assPart = new TLorentzVector();
    Float_t eta = trgHadEtaAway + deltaEta;
    Float_t phi = trgHadPhiAway + deltaPhi;
    if (eta > .9) {
      delete assPart;
      continue;
    }
    if (phi < 0.)
      phi += 2. * TMath::Pi();
    else if (phi > 2. * TMath::Pi())
      phi -= 2. * TMath::Pi();
    assPart->SetPtEtaPhiM(3., eta, phi, .14);
    assHadHadArray->Add(assPart);
    if (gRandom->Uniform() < pFraction)
      assProtHadArray->Add(assPart);
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
  fHistCorrEtaPhi(0x0),
  fHistCorrAvgEtaPhi(0x0),
  fHistCorrTrgEtaPhi(0x0),
  fHistCorrAssEtaPhi(0x0)
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
  fHistCorrAvgEtaPhi = new TH2F(Form("%s_avgetaphi", name.Data()), ";#Delta#phi;avg #eta",
				100, -1., 2*TMath::Pi()-1.,
				100, -2., 2.);
  fHistCorrTrgEtaPhi = new TH2F(Form("%s_trg_etaphi", name.Data()), ";#varphi;#eta",
				100, 0., 2*TMath::Pi(),
				100, -1., 1.);
  fHistCorrAssEtaPhi = new TH2F(Form("%s_ass_etaphi", name.Data()), ";#varphi;#eta",
				100, 0., 2*TMath::Pi(),
				100, -1., 1.);

  fOutputList->Add(fHistStat);
  fOutputList->Add(fHistCorrPhi);
  fOutputList->Add(fHistCorrPhi2);
  fOutputList->Add(fHistCorrEtaPhi);
  fOutputList->Add(fHistCorrAvgEtaPhi);
  fOutputList->Add(fHistCorrTrgEtaPhi);
  fOutputList->Add(fHistCorrAssEtaPhi);
}

AliAnalysisTaskJetProtonCorr::AliHistCorr::~AliHistCorr()
{
  // dtor
}

void AliAnalysisTaskJetProtonCorr::AliHistCorr::Fill(AliVParticle *trgPart, AliVParticle *assPart, Float_t weight)
{
  Float_t deltaEta = assPart->Eta() - trgPart->Eta();
  Float_t avgEta = (assPart->Eta() + trgPart->Eta()) / 2.;
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
  fHistCorrAvgEtaPhi->Fill(deltaPhi, avgEta, weight);
}

void AliAnalysisTaskJetProtonCorr::AliHistCorr::Fill(TLorentzVector *trgPart, AliVParticle *assPart, Float_t weight)
{
  Float_t deltaEta = assPart->Eta() - trgPart->Eta();
  Float_t avgEta = (assPart->Eta() + trgPart->Eta()) / 2.;
  Float_t deltaPhi = assPart->Phi() - trgPart->Phi();
  if (deltaPhi > (2.*TMath::Pi()-1.))
    deltaPhi -= 2. * TMath::Pi();
  else if (deltaPhi < -1.)
    deltaPhi += 2. * TMath::Pi();
  // printf("trg: pt = %5.2f, phi = %5.2f, eta = %5.2f; ass: pt = %5.2f, phi = %5.2f, eta = %5.2f; deltaphi = %5.2f, deltaeta = %5.2f\n",
  // 	 trgPart->Pt(), trgPart->Phi(), trgPart->Eta(), assPart->Pt(), assPart->Phi(), assPart->Eta(), deltaPhi, deltaEta);

  fHistCorrPhi->Fill(deltaPhi);
  Float_t trgPhi = trgPart->Phi();
  if (trgPhi < 0)
    trgPhi += 2. * TMath::Pi();
  fHistCorrPhi2->Fill(trgPhi, assPart->Phi(), weight);
  fHistCorrEtaPhi->Fill(deltaPhi, deltaEta, weight);
  fHistCorrAvgEtaPhi->Fill(deltaPhi, avgEta, weight);
}

void AliAnalysisTaskJetProtonCorr::AliHistCorr::Fill(TLorentzVector *trgPart, TLorentzVector *assPart, Float_t weight)
{
  Float_t deltaEta = assPart->Eta() - trgPart->Eta();
  Float_t avgEta = (assPart->Eta() + trgPart->Eta()) / 2.;
  Float_t deltaPhi = assPart->Phi() - trgPart->Phi();
  if (deltaPhi > (2.*TMath::Pi()-1.))
    deltaPhi -= 2. * TMath::Pi();
  else if (deltaPhi < -1.)
    deltaPhi += 2. * TMath::Pi();
  // printf("trg: pt = %5.2f, phi = %5.2f, eta = %5.2f; ass: pt = %5.2f, phi = %5.2f, eta = %5.2f; deltaphi = %5.2f, deltaeta = %5.2f\n",
  // 	 trgPart->Pt(), trgPart->Phi(), trgPart->Eta(), assPart->Pt(), assPart->Phi(), assPart->Eta(), deltaPhi, deltaEta);

  fHistCorrPhi->Fill(deltaPhi);
  Float_t trgPhi = trgPart->Phi();
  if (trgPhi < 0)
    trgPhi += 2. * TMath::Pi();
  Float_t assPhi = assPart->Phi();
  if (assPhi < 0)
    assPhi += 2. * TMath::Pi();
  fHistCorrPhi2->Fill(trgPhi, assPhi, weight);
  fHistCorrEtaPhi->Fill(deltaPhi, deltaEta, weight);
  fHistCorrAvgEtaPhi->Fill(deltaPhi, avgEta, weight);
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
