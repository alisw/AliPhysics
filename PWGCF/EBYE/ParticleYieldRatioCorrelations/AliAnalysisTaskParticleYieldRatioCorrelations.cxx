#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TFile.h"
#include <TRandom.h>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskParticleYieldRatioCorrelations.h"
#include "AliAODMCParticle.h"
#include "AliMultSelection.h"
class AliMultSelection;
class AliAnalysisTaskParticleYieldRatioCorrelations;

using namespace std;

ClassImp(AliAnalysisTaskParticleYieldRatioCorrelations)

    AliAnalysisTaskParticleYieldRatioCorrelations::AliAnalysisTaskParticleYieldRatioCorrelations()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fInputList(0),
      NEvents(0),
      fAliEventCuts(0),
      fDeDx(0),
      fTOF(0),
      fHistQAVx(0),
      fHistQAVy(0),
      fHistQAVz(0),
      fHistQAPt(0),
      fHistQAEta(0),
      fHistQAPhi(0),
      fHistQAClustersTPC(0),
      fHistQACrossedRowsTPC(0),
      fHistQAChi2perNDF(0),
      fHistQAClustersITS(0),
      fHistQAEtaPhi(0),
          fHistQAMomPt(0)

{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskParticleYieldRatioCorrelations::AliAnalysisTaskParticleYieldRatioCorrelations(const char *name)
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fInputList(0),
      NEvents(0),
      fAliEventCuts(0),
      fDeDx(0),
      fTOF(0),
      fHistQAVx(0),
      fHistQAVy(0),
      fHistQAVz(0),
      fHistQAPt(0),
      fHistQAEta(0),
      fHistQAPhi(0),
      fHistQAClustersTPC(0),
      fHistQACrossedRowsTPC(0),
      fHistQAChi2perNDF(0),
      fHistQAClustersITS(0),
      fHistQAEtaPhi(0),
          fHistQAMomPt(0)
{
  // constructor
  DefineInput(0, TChain::Class());
  DefineInput(1, TList::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskParticleYieldRatioCorrelations::~AliAnalysisTaskParticleYieldRatioCorrelations()
{
  // destructor
  if (fOutputList)
  {
    delete fOutputList;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelations::UserCreateOutputObjects()
{
  // Get Efficiency maps for futher corrections
  fInputList = (TList *)GetInputData(1);
  for (int iCent = 0; iCent < nCentrClassesUsed; iCent++)
  {
    for (int iEta = 0; iEta < nEtaClasses; iEta++)
    {
      for (int iSort = 0; iSort < 6; ++iSort)
      {
        EfficiencyTracking[iCent][iEta][iSort] = (TH3D *)fInputList->FindObject(Form("EfficiencyTrackingC%dEta%dSort%d", iCent, iEta, iSort));
        ContaminationTracking[iCent][iEta][iSort] = (TH3D *)fInputList->FindObject(Form("ContaminationTrackingC%dEta%dSort%d", iCent, iEta, iSort));
        EfficiencyPID[iCent][iEta][iSort] = (TH3D *)fInputList->FindObject(Form("EfficiencyPIDC%dEta%dSort%d", iCent, iEta, iSort));
        ContaminationPID[iCent][iEta][iSort] = (TH3D *)fInputList->FindObject(Form("ContaminationPIDC%dEta%dSort%d", iCent, iEta, iSort));
      }
    }
  }

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  NEvents = new TH1D("NEvents", "", nCentrClassesUsed, minCent, maxCent);
  for (int iSub = 0; iSub < nSubsamples; iSub++)
  {
    NEventsSub[iSub] = new TH1D(Form("NEventsSub%d", iSub), "", nCentrClassesUsed, minCent, maxCent);
    fOutputList->Add(NEventsSub[iSub]);
  }
  fHistEventsCut = new TH1F("fHistEventsCut", ";NEvents", 6, 0, 6);
  fHistTracksCut = new TH1F("fHistTracksCut", ";Ntracks", 6, 0, 6);
  fHistEventsCut->GetXaxis()->SetBinLabel(1, "minBias");
  fHistEventsCut->GetXaxis()->SetBinLabel(2, "centrality");
  fHistEventsCut->GetXaxis()->SetBinLabel(3, "SPDvsV0M");
  fHistEventsCut->GetXaxis()->SetBinLabel(4, "NContributors");
  fHistEventsCut->GetXaxis()->SetBinLabel(5, "vertex");
  fHistEventsCut->GetXaxis()->SetBinLabel(6, "AliEventCuts");
  fHistTracksCut->GetXaxis()->SetBinLabel(1, "minBias");
  fHistTracksCut->GetXaxis()->SetBinLabel(2, "FilterBit");
  fHistTracksCut->GetXaxis()->SetBinLabel(3, "Eta");
  fHistTracksCut->GetXaxis()->SetBinLabel(4, "pT");
  fHistTracksCut->GetXaxis()->SetBinLabel(4, "TOF");

  fOutputList->Add(NEvents);
  fOutputList->Add(fHistEventsCut);
  fOutputList->Add(fHistTracksCut);

  fDeDx = new TH2D("DeDx", ";p_{TPC};dE/dx (a.u.)", 200, 0.2, 2, 250, 0, 2000);
  fTOF = new TH2D("fTOF", ";p_{TPC};t", 200, 0.2, 2, 200, -4000, 4000);
  fOutputList->Add(fDeDx);
  fOutputList->Add(fTOF);
  for (int iSort = 0; iSort < 4; iSort++)
  {
    fDeDxSorts[iSort] = new TH2D(Form("DeDxSort%d", iSort), ";p_{TPC};dE/dx (a.u.)", 200, 0.2, 2, 250, 0, 2000);
    fTOFSorts[iSort] = new TH2D(Form("fTOFSort%d", iSort), ";p_{TPC};t", 200, 0.2, 2, 200, -4000, 4000);
    fOutputList->Add(fDeDxSorts[iSort]);
    fOutputList->Add(fTOFSorts[iSort]);
  }

  for (int iCent = 0; iCent < nCentrClassesUsed; iCent++)
  {
    fHistMomentsInAllAccRec[iCent] = new TH1D(Form("fHistMomentsInAllAccRec_C%d", iCent), ";;Number Of Tracks", 3 * 2 + 4, 0, 3 * 2 + 4);
    fHistMomentsInAllAccGen[iCent] = new TH1D(Form("fHistMomentsInAllAccGen_C%d", iCent), ";;Number Of Tracks", 3 * 2 + 4, 0, 3 * 2 + 4);
    fHistMomentsInAllAccReg[iCent] = new TH1D(Form("fHistMomentsInAllAccReg_C%d", iCent), ";;Number Of Tracks", 3 * 2 + 4, 0, 3 * 2 + 4);
    fOutputList->Add(fHistMomentsInAllAccRec[iCent]);
    fOutputList->Add(fHistMomentsInAllAccGen[iCent]);
    fOutputList->Add(fHistMomentsInAllAccReg[iCent]);

    for (int iSub = 0; iSub < nSubsamples; iSub++)
    {

      fHistMomentsInAllAccRecSub[iSub][iCent] = new TH1D(Form("fHistMomentsInAllAccRecSub%d_C%d", iSub, iCent), ";;Number Of Tracks", 3 * 2 + 4, 0, 3 * 2 + 4);
      fHistMomentsInAllAccGenSub[iSub][iCent] = new TH1D(Form("fHistMomentsInAllAccGenSub%d_C%d", iSub, iCent), ";;Number Of Tracks", 3 * 2 + 4, 0, 3 * 2 + 4);
      fHistMomentsInAllAccRegSub[iSub][iCent] = new TH1D(Form("fHistMomentsInAllAccRegSub%d_C%d", iSub, iCent), ";;Number Of Tracks", 3 * 2 + 4, 0, 3 * 2 + 4);
      fOutputList->Add(fHistMomentsInAllAccRecSub[iSub][iCent]);
      fOutputList->Add(fHistMomentsInAllAccGenSub[iSub][iCent]);
      fOutputList->Add(fHistMomentsInAllAccRegSub[iSub][iCent]);
    }

    for (int iSort = 0; iSort < 6; iSort++)
    {
      fHistMomentsInEtaRec[iCent][iSort] = new TH1D(Form("fHistMomentsInEtaRec_C%d_S%d", iCent, iSort), ";;Number Of Tracks", nEtaClasses, -0.8, 0.8);
      fHistMomentsInEtaGen[iCent][iSort] = new TH1D(Form("fHistMomentsInEtaGen_C%d_S%d", iCent, iSort), ";;Number Of Tracks", nEtaClasses, -0.8, 0.8);
      fHistMomentsInEtaReg[iCent][iSort] = new TH1D(Form("fHistMomentsInEtaReg_C%d_S%d", iCent, iSort), ";;Number Of Tracks", nEtaClasses, -0.8, 0.8);
      fOutputList->Add(fHistMomentsInEtaRec[iCent][iSort]);
      fOutputList->Add(fHistMomentsInEtaGen[iCent][iSort]);
      fOutputList->Add(fHistMomentsInEtaReg[iCent][iSort]);

      fHistMomentsInPhiRec[iCent][iSort] = new TH1D(Form("fHistMomentsInPhiRec_C%d_S%d", iCent, iSort), ";;Number Of Tracks", nPhiWindows, 0, TMath::TwoPi());
      fHistMomentsInPhiGen[iCent][iSort] = new TH1D(Form("fHistMomentsInPhiGen_C%d_S%d", iCent, iSort), ";;Number Of Tracks", nPhiWindows, 0, TMath::TwoPi());
      fHistMomentsInPhiReg[iCent][iSort] = new TH1D(Form("fHistMomentsInPhiReg_C%d_S%d", iCent, iSort), ";;Number Of Tracks", nPhiWindows, 0, TMath::TwoPi());
      fOutputList->Add(fHistMomentsInPhiRec[iCent][iSort]);
      fOutputList->Add(fHistMomentsInPhiGen[iCent][iSort]);
      fOutputList->Add(fHistMomentsInPhiReg[iCent][iSort]);

      for (int iSub = 0; iSub < nSubsamples; iSub++)
      {
        fHistMomentsInEtaRecSub[iSub][iCent][iSort] = new TH1D(Form("fHistMomentsInEtaRecSub%d_C%d_S%d", iSub, iCent, iSort), ";;Number Of Tracks", nEtaClasses, -0.8, 0.8);
        fHistMomentsInEtaGenSub[iSub][iCent][iSort] = new TH1D(Form("fHistMomentsInEtaGenSub%d_C%d_S%d", iSub, iCent, iSort), ";;Number Of Tracks", nEtaClasses, -0.8, 0.8);
        fHistMomentsInEtaRegSub[iSub][iCent][iSort] = new TH1D(Form("fHistMomentsInEtaRegSub%d_C%d_S%d", iSub, iCent, iSort), ";;Number Of Tracks", nEtaClasses, -0.8, 0.8);
        fOutputList->Add(fHistMomentsInEtaRecSub[iSub][iCent][iSort]);
        fOutputList->Add(fHistMomentsInEtaGenSub[iSub][iCent][iSort]);
        fOutputList->Add(fHistMomentsInEtaRegSub[iSub][iCent][iSort]);

        fHistMomentsInPhiRecSub[iSub][iCent][iSort] = new TH1D(Form("fHistMomentsInPhiRecSub%d_C%d_S%d", iSub, iCent, iSort), ";;Number Of Tracks", nPhiWindows, 0, TMath::TwoPi());
        fHistMomentsInPhiGenSub[iSub][iCent][iSort] = new TH1D(Form("fHistMomentsInPhiGenSub%d_C%d_S%d", iSub, iCent, iSort), ";;Number Of Tracks", nPhiWindows, 0, TMath::TwoPi());
        fHistMomentsInPhiRegSub[iSub][iCent][iSort] = new TH1D(Form("fHistMomentsInPhiRegSub%d_C%d_S%d", iSub, iCent, iSort), ";;Number Of Tracks", nPhiWindows, 0, TMath::TwoPi());
        fOutputList->Add(fHistMomentsInPhiRecSub[iSub][iCent][iSort]);
        fOutputList->Add(fHistMomentsInPhiGenSub[iSub][iCent][iSort]);
        fOutputList->Add(fHistMomentsInPhiRegSub[iSub][iCent][iSort]);
      }
    }

    for (int iSort = 0; iSort < 6; iSort++)
    {
      int skipSort = (2 * (6) - (iSort - 1)) * iSort / 2; // how many hists have filled in cross moments hists array, allow to avoid two same histograms
      for (int jSort = iSort; jSort < 6; jSort++)
      {
        fHistCrossMomentsInEtaRec[iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInEtaRec_C%d_S%dvs%d", iCent, iSort, jSort), ";;Number Of Tracks", nEtaClasses * nEtaClasses, 0, nEtaClasses * nEtaClasses);
        fHistCrossMomentsInEtaGen[iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInEtaGen_C%d_S%dvs%d", iCent, iSort, jSort), ";;Number Of Tracks", nEtaClasses * nEtaClasses, 0, nEtaClasses * nEtaClasses);
        fHistCrossMomentsInEtaReg[iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInEtaReg_C%d_S%dvs%d", iCent, iSort, jSort), ";;Number Of Tracks", nEtaClasses * nEtaClasses, 0, nEtaClasses * nEtaClasses);
        fOutputList->Add(fHistCrossMomentsInEtaRec[iCent][skipSort + jSort - iSort]);
        fOutputList->Add(fHistCrossMomentsInEtaGen[iCent][skipSort + jSort - iSort]);
        fOutputList->Add(fHistCrossMomentsInEtaReg[iCent][skipSort + jSort - iSort]);

        fHistCrossMomentsInPhiRec[iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInPhiRec_C%d_S%dvs%d", iCent, iSort, jSort), ";;Number Of Tracks", nPhiWindows * nPhiWindows, 0, nPhiWindows * nPhiWindows);
        fHistCrossMomentsInPhiGen[iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInPhiGen_C%d_S%dvs%d", iCent, iSort, jSort), ";;Number Of Tracks", nPhiWindows * nPhiWindows, 0, nPhiWindows * nPhiWindows);
        fHistCrossMomentsInPhiReg[iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInPhiReg_C%d_S%dvs%d", iCent, iSort, jSort), ";;Number Of Tracks", nPhiWindows * nPhiWindows, 0, nPhiWindows * nPhiWindows);
        fOutputList->Add(fHistCrossMomentsInPhiRec[iCent][skipSort + jSort - iSort]);
        fOutputList->Add(fHistCrossMomentsInPhiGen[iCent][skipSort + jSort - iSort]);
        fOutputList->Add(fHistCrossMomentsInPhiReg[iCent][skipSort + jSort - iSort]);

        for (int iSub = 0; iSub < nSubsamples; iSub++)
        {
          fHistCrossMomentsInEtaRecSub[iSub][iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInEtaRecSub%d_C%d_S%dvs%d", iSub, iCent, iSort, jSort), ";;Number Of Tracks", nEtaClasses * nEtaClasses, 0, nEtaClasses * nEtaClasses);
          fHistCrossMomentsInEtaGenSub[iSub][iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInEtaGenSub%d_C%d_S%dvs%d", iSub, iCent, iSort, jSort), ";;Number Of Tracks", nEtaClasses * nEtaClasses, 0, nEtaClasses * nEtaClasses);
          fHistCrossMomentsInEtaRegSub[iSub][iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInEtaRegSub%d_C%d_S%dvs%d", iSub, iCent, iSort, jSort), ";;Number Of Tracks", nEtaClasses * nEtaClasses, 0, nEtaClasses * nEtaClasses);
          fOutputList->Add(fHistCrossMomentsInEtaRecSub[iSub][iCent][skipSort + jSort - iSort]);
          fOutputList->Add(fHistCrossMomentsInEtaGenSub[iSub][iCent][skipSort + jSort - iSort]);
          fOutputList->Add(fHistCrossMomentsInEtaRegSub[iSub][iCent][skipSort + jSort - iSort]);

          fHistCrossMomentsInPhiRecSub[iSub][iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInPhiRecSub%d_C%d_S%dvs%d", iSub, iCent, iSort, jSort), ";;Number Of Tracks", nPhiWindows * nPhiWindows, 0, nPhiWindows * nPhiWindows);
          fHistCrossMomentsInPhiGenSub[iSub][iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInPhiGenSub%d_C%d_S%dvs%d", iSub, iCent, iSort, jSort), ";;Number Of Tracks", nPhiWindows * nPhiWindows, 0, nPhiWindows * nPhiWindows);
          fHistCrossMomentsInPhiRegSub[iSub][iCent][skipSort + jSort - iSort] = new TH1D(Form("fHistCrossMomentsInPhiRegSub%d_C%d_S%dvs%d", iSub, iCent, iSort, jSort), ";;Number Of Tracks", nPhiWindows * nPhiWindows, 0, nPhiWindows * nPhiWindows);
          fOutputList->Add(fHistCrossMomentsInPhiRecSub[iSub][iCent][skipSort + jSort - iSort]);
          fOutputList->Add(fHistCrossMomentsInPhiGenSub[iSub][iCent][skipSort + jSort - iSort]);
          fOutputList->Add(fHistCrossMomentsInPhiRegSub[iSub][iCent][skipSort + jSort - iSort]);
        }
      }
    }
  }

  fHistQAPt = new TH1F("fHistQAPt", "p_{T} distribution; p_{T}, GeV/c; dN/dp_{T}", 800, 0.0, 20.0);
  fOutputList->Add(fHistQAPt);

  fHistQAEta = new TH1F("fHistQAEta", "#eta distribution;#eta;dN/d#eta", 500, -10, 10);
  fOutputList->Add(fHistQAEta);

  fHistQAPhi = new TH1F("fHistQAPhi", "#varphi distribution;#varphi;dN/d#varphi", 360, 0, 2 * TMath::Pi());
  fOutputList->Add(fHistQAPhi);

  fHistQAEtaPhi = new TH2D("fHistQAEtaPhi", "N tracks in (#eta, #varphi);#eta;#varphi", 25, -0.8, 0.8, 25, 0, 2 * TMath::Pi());
  fOutputList->Add(fHistQAEtaPhi);

  fHistQAMomPt = new TH2D("fHistQAMomPt", "N tracks;Momentum TPC;Pt", 100, 0, 3, 100, 0, 3);
  fOutputList->Add(fHistQAMomPt);

  for (Int_t i(0); i < 3; i++)
  {
    fHistQASPDTrackletsvsV0MCent[i] = new TH2D(Form("fHistQASPDTrackletsvsV0MCent%d", i), ";V0M Percentile;N Tracklets in SPD", 100, 0, 100, 400, 0, 1e4);
    fOutputList->Add(fHistQASPDTrackletsvsV0MCent[i]);
  }
  for (Int_t i(0); i < 2; i++)
  {
    fHistQAMultTPCvsESD[i] = new TH2D(Form("fHistQAMultTPCvsESD%d", i), ";MultTPC;MultESD", 400, 0, 1e4, 400, 0, 5e4);
    fOutputList->Add(fHistQAMultTPCvsESD[i]);
    fHistQAMultTPCvsV0[i] = new TH2D(Form("fHistQAMultTPCvsV0%d", i), ";MultTPC;V0", 400, 0, 1e4, 400, 0, 5e4);
    fOutputList->Add(fHistQAMultTPCvsV0[i]);
    fHistQAMultTrkvsMultTrkTOF[i] = new TH2D(Form("fHistQAMultTrkvsMultTrkTOF%d", i), ";MultTrk;MultTrkTOF", 400, 0, 5e3, 400, 0, 5e3);
    fOutputList->Add(fHistQAMultTrkvsMultTrkTOF[i]);
  }

  fHistQAVx = new TH1D("fHistQAVx", "Primary vertex distribution - x coordinate;V_{x} (cm);Entries", 400, -1, 1);
  fOutputList->Add(fHistQAVx);
  fHistQAVy = new TH1D("fHistQAVy", "Primary vertex distribution - y coordinate;V_{y} (cm);Entries", 400, -1, 1);
  fOutputList->Add(fHistQAVy);
  fHistQAVz = new TH1D("fHistQAVz", "Primary vertex distribution - z coordinate;V_{z} (cm);Entries", 200, -20., 20.);
  fOutputList->Add(fHistQAVz);

  fHistQAClustersTPC = new TH1D("fHistQAClustersTPC", "N Clusters TPC;N_{TPC clusters};Entries", 161, -0.5, 160.5);
  fOutputList->Add(fHistQAClustersTPC);

  fHistQACrossedRowsTPC = new TH1D("fHistQACrossedRowsTPC", "N Crossed Rows TPC;N_{TPC CrossedRows};Entries", 161, -0.5, 160.5);
  fOutputList->Add(fHistQACrossedRowsTPC);

  fHistQAChi2perNDF = new TH1D("fHistQAChi2perNDF", ";#chi^{2} / ndf;n tracks", 600, -1, 5);
  fOutputList->Add(fHistQAChi2perNDF);

  fHistQAClustersITS = new TH1D("fHistQAClustersITS", "N Clusters ITS;N_ITS_clusters;Entries", 11, -0.5, 10.5);
  fOutputList->Add(fHistQAClustersITS);

  fAliEventCuts = new AliEventCuts();
  if (pbpb)
    fAliEventCuts->SetupPbPb2018();
  else
    fAliEventCuts->SetupRun2pp();
  fAliEventCuts->AddQAplotsToList(fOutputList);

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man)
  {
    AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler *>(man->GetInputEventHandler());
    if (inputHandler)
      fPIDResponse = inputHandler->GetPIDResponse();
    else
      AliFatal("Input handler needed");
  }
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelations::UserExec(Option_t *)
{
  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

  if (!fAOD)
  {
    PostData(1, fOutputList);
    return;
  }
  fHistEventsCut->Fill("minBias", 1);
  if (!fAliEventCuts->AcceptEvent(fAOD))
  {
    PostData(1, fOutputList);
    return;
  }
  fHistEventsCut->Fill("AliEventCuts", 1);
  Float_t centr;
  Int_t NTrackletsSPD;
  if (0)
  {
    AliCentrality *centrality = fAOD->GetCentrality();
    centr = centrality->GetCentralityPercentile("V0M");
  }
  else
  {
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
    if (MultSelection)
    {
      centr = MultSelection->GetMultiplicityPercentile("V0M");
      NTrackletsSPD = MultSelection->GetEstimator("SPDTracklets")->GetValue();
    }
    else
    {
      fHistEventsCut->Fill("AliEventCuts", 1);
      PostData(1, fOutputList);
      return;
    }
  }
  // centrality cut:
  if (centr < minCent || centr > maxCent)
  {
    PostData(1, fOutputList);
    return;
  }
  fHistEventsCut->Fill("centrality", 1);
  // cut fot SPD Tracklets vs V0M Centrality percentle
  fHistQASPDTrackletsvsV0MCent[0]->Fill(centr, NTrackletsSPD);
  if (SPDvsV0MCut)
  {
    TF1 *fSPDvsV0M_DownLimit = new TF1("fSPDvsV0M_DownLimit", "exp(8.1456-0.0354*x-3.8e-04*x*x)", 0, 80);
    TF1 *fSPDvsV0M_UperLimit = new TF1("fSPDvsV0M_UperLimit", "exp(8.58-0.036*x-4.4e-05*x*x)", 0, 80);
    if (NTrackletsSPD < fSPDvsV0M_DownLimit->Eval(centr) || NTrackletsSPD > fSPDvsV0M_UperLimit->Eval(centr))
    {
      PostData(1, fOutputList);
      return;
    }
    fHistEventsCut->Fill("SPDvsV0M", 1);
  }
  fHistQASPDTrackletsvsV0MCent[1]->Fill(centr, NTrackletsSPD);
  // Events with large number of TPC clusters
  Int_t multEsd = ((AliAODHeader *)fAOD->GetHeader())->GetNumberOfESDTracks();
  Int_t multTPC = 0;
  // remaining events cuts
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  if (vtx->GetNContributors() < 1)
  {
    PostData(1, fOutputList);
    return;
  }
  fHistEventsCut->Fill("NContributors", 1);
  if (fabs((Float_t)vtx->GetZ()) > 8)
  {
    PostData(1, fOutputList);
    return;
  }
  fHistEventsCut->Fill("vertex", 1);
  Float_t vertex = (Float_t)vtx->GetZ();
  fHistQAVz->Fill(vtx->GetZ());
  fHistQAVx->Fill(vtx->GetX());
  fHistQAVy->Fill(vtx->GetX());
  Int_t nTracks(fAOD->GetNumberOfTracks());
  Int_t EtaBin, CentrBin, phiBin;
  int sort = 20;
  // CentrBin = (centr - minCent) / ((maxCent - minCent) / nCentrClasses);
  // CentrBin = centr<=5?0:(centr<=10?1: 2+(centr-10)/((maxCent-10)/(nCentrClasses-2)) );
  // CentrBin = 20
  for (Int_t i(0); i < nCentrClassesUsed; i++)
  {
    if (centr >= CentrPercentiles[i] && centr <= CentrPercentiles[i + 1])
    {
      CentrBin = i;
    }
  }
  if (CentrBin < 0 || CentrBin >= nCentrClassesUsed)
  {
    PostData(1, fOutputList);
    return;
  }
  int NAcceptedtracks = 0, NAcceptedtracksinEta[nEtaClasses] = {0};
  Int_t randomSubsample = gRandom->Integer(nSubsamples);

  // declare collector variables
  Double_t NTracksInCutEta[nEtaClasses][nSorts] = {0}, NTracksInCut2Eta[nEtaClasses][nSorts] = {0},
           GenNParticlesEta[nEtaClasses][nSorts] = {0}, NTracksRegEta[nEtaClasses][nSorts] = {0};

  Double_t NTracksInCutPhi[nPhiWindows][nSorts] = {0}, NTracksInCut2Phi[nPhiWindows][nSorts] = {0},
           GenNParticlesPhi[nPhiWindows][nSorts] = {0}, NTracksRegPhi[nPhiWindows][nSorts] = {0};

  fHistTracksCut->Fill("minBias", nTracks);

  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  Float_t PtCut[3] = {0.2, 0.5, 0.5};
  Float_t nSigmaBoundary[3] = {0.5, 0.32, 0.7};
  Int_t multTrk = 0;
  Int_t multTrkTOF = 0;
  // tracks loop
  for (Int_t i(0); i < nTracks; i++)
  {
    AliAODTrack *track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (track->TestFilterBit(32))
    {
      multTrk++;
      if (TMath::Abs(track->GetTOFsignalDz()) <= 10 && track->GetTOFsignal() >= 12000 && track->GetTOFsignal() <= 25000)
        multTrkTOF++;
    }
    if (track->TestFilterBit(128))
      multTPC++;
    if (!track || !track->TestFilterBit(filterBit))
      continue;
    fHistTracksCut->Fill("FilterBit", 1);
    if (fabs(track->Eta()) >= 0.8)
      continue;
    fHistTracksCut->Fill("Eta", 1);
    if (track->Pt() < minP || track->Pt() > maxP)
      continue;
    if (track->GetTPCCrossedRows() <= nCrossedRows)
    {
      continue;
    }
    fHistTracksCut->Fill("pT", 1);
    EtaBin = (0.8 + track->Eta()) / (1.6 / nEtaClasses);
    // phiBin = track->Phi() / (TMath::TwoPi() / nPhiWindows);
    phiBin = (track->Phi() + (movePhi * TMath::TwoPi() / nPhiWindows / 2.0)) / (TMath::TwoPi() / nPhiWindows); // move phi on half of bin for stat errors estimation
    if (phiBin == nPhiWindows)
      phiBin = 0;
    Float_t nOfSigmasTPC_el = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    Float_t nOfSigmasTPC_pi = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    Float_t nOfSigmasTPC_K = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Float_t nOfSigmasTPC_p = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Float_t nOfSigmasTPC_d = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);

    Float_t nOfSigmasTOF_pi = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion, fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nOfSigmasTOF_K = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon, fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nOfSigmasTOF_p = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton, fPIDResponse->GetTOFResponse().GetTimeZero());

    if (fabs(nOfSigmasTOF_pi) > 900 || fabs(nOfSigmasTOF_K) > 900 || fabs(nOfSigmasTOF_p) > 900)
      fHistTracksCut->Fill("TOF", 1);

    Float_t Pt = track->Pt();
    Float_t Px = track->Px();
    Float_t Py = track->Py();
    Float_t Pz = track->Pz();
    Float_t Moment = track->GetTPCmomentum();
    Float_t Phi = track->Phi();
    Float_t Eta = track->Eta();
    Float_t DeDx = track->GetTPCsignal();
    int Charge = track->Charge();
    fDeDx->Fill(Moment, DeDx);
    fTOF->Fill(Moment, track->GetTOFsignal());
    NAcceptedtracksinEta[EtaBin]++;

    fHistQAPt->Fill(Pt);
    fHistQAEta->Fill(Eta);
    fHistQAPhi->Fill(Phi);
    fHistQAEtaPhi->Fill(Eta, Phi);
    fHistQAMomPt->Fill(Moment, Pt);

    fHistQAClustersTPC->Fill(track->GetTPCNcls());
    fHistQACrossedRowsTPC->Fill(track->GetTPCCrossedRows());
    fHistQAClustersITS->Fill(track->GetITSclusters(0));
    fHistQAChi2perNDF->Fill(track->Chi2perNDF());

    // Get correction coefficients from efficiency maps
    Double_t coeff[nSorts] = {0};
    for (int iSort = 0; iSort < 6; iSort++)
    {
      if (!EfficiencyPID[CentrBin][EtaBin][iSort] || !EfficiencyTracking[CentrBin][EtaBin][iSort] || !ContaminationPID[CentrBin][EtaBin][iSort] || !ContaminationTracking[CentrBin][EtaBin][iSort])
        continue;
      if (
          EfficiencyPID[CentrBin][EtaBin][iSort]->GetBinContent(EfficiencyPID[CentrBin][EtaBin][iSort]->FindBin(Pt, Phi, vertex)) > 1e-10 &&
          EfficiencyTracking[CentrBin][EtaBin][iSort]->GetBinContent(EfficiencyTracking[CentrBin][EtaBin][iSort]->FindBin(Pt, Phi, vertex)) > 1e-10)
        coeff[iSort] = 1;
      if (PIDCorr)
        coeff[iSort] *= (1 - ContaminationPID[CentrBin][EtaBin][iSort]->GetBinContent(ContaminationPID[CentrBin][EtaBin][iSort]->FindBin(Pt, Phi, vertex))) /
                        EfficiencyPID[CentrBin][EtaBin][iSort]->GetBinContent(EfficiencyPID[CentrBin][EtaBin][iSort]->FindBin(Pt, Phi, vertex));
      if (TrackingCorr)
        coeff[iSort] *= (1 - ContaminationTracking[CentrBin][EtaBin][iSort]->GetBinContent(ContaminationTracking[CentrBin][EtaBin][iSort]->FindBin(Pt, Phi, vertex))) /
                        EfficiencyTracking[CentrBin][EtaBin][iSort]->GetBinContent(EfficiencyTracking[CentrBin][EtaBin][iSort]->FindBin(Pt, Phi, vertex));
    }

    Float_t nSigma_comb_pi = sqrt(nOfSigmasTPC_pi * nOfSigmasTPC_pi + nOfSigmasTOF_pi * nOfSigmasTOF_pi);
    Float_t nSigma_comb_K = sqrt(nOfSigmasTPC_K * nOfSigmasTPC_K + nOfSigmasTOF_K * nOfSigmasTOF_K);
    Float_t nSigma_comb_p = sqrt(nOfSigmasTPC_p * nOfSigmasTPC_p + nOfSigmasTOF_p * nOfSigmasTOF_p);

    // test if particle fits identification criterias
    // if true, add to collector with correction coefficients as weight factor
    bool selected_pi = false;
    if (Pt < nSigmaBoundary[0] && fabs(nOfSigmasTPC_pi) < nSigma && fabs(nOfSigmasTPC_K) > 3 && fabs(nOfSigmasTPC_p) > 3 && fabs(nOfSigmasTPC_el) > 1)
      selected_pi = true;
    if (Pt > nSigmaBoundary[0] && nSigma_comb_pi < nSigma && nSigma_comb_K > 3 && nSigma_comb_p > 3)
      selected_pi = true;
    if (selected_pi)
    {
      if (Charge < 0)
      {
        NTracksInCutEta[EtaBin][0] += coeff[0];
        NTracksInCut2Eta[EtaBin][0] += coeff[0] * coeff[0];
        NTracksRegEta[EtaBin][0] += 1;
        NTracksInCutPhi[phiBin][0] += coeff[0];
        NTracksInCut2Phi[phiBin][0] += coeff[0] * coeff[0];
        NTracksRegPhi[phiBin][0] += 1;
      }
      if (Charge > 0)
      {
        NTracksInCutEta[EtaBin][1] += coeff[1];
        NTracksInCut2Eta[EtaBin][1] += coeff[1] * coeff[1];
        NTracksRegEta[EtaBin][1] += 1;
        NTracksInCutPhi[phiBin][1] += coeff[1];
        NTracksInCut2Phi[phiBin][1] += coeff[1] * coeff[1];
        NTracksRegPhi[phiBin][1] += 1;
      }
      fDeDxSorts[0]->Fill(Moment, DeDx);
      fTOFSorts[0]->Fill(Moment, track->GetTOFsignal());
    }

    bool selected_K = false;
    if (Pt < nSigmaBoundary[1] && fabs(nOfSigmasTPC_K) < nSigma && fabs(nOfSigmasTPC_pi) > 3 && fabs(nOfSigmasTPC_p) > 3) // && fabs( nOfSigmasTPC_el)>1 )
      selected_K = true;
    if (Pt > nSigmaBoundary[1] && nSigma_comb_K < nSigma && nSigma_comb_pi > 3 && nSigma_comb_p > 3)
      selected_K = true;
    if (selected_K && Pt > PtCut[1])
    {
      if (Charge < 0)
      {
        NTracksInCutEta[EtaBin][2] += coeff[2];
        NTracksInCut2Eta[EtaBin][2] += coeff[2] * coeff[2];
        NTracksRegEta[EtaBin][2] += 1;
        NTracksInCutPhi[phiBin][2] += coeff[2];
        NTracksInCut2Phi[phiBin][2] += coeff[2] * coeff[2];
        NTracksRegPhi[phiBin][2] += 1;
      }
      if (Charge > 0)
      {
        NTracksInCutEta[EtaBin][3] += coeff[3];
        NTracksInCut2Eta[EtaBin][3] += coeff[3] * coeff[3];
        NTracksRegEta[EtaBin][3] += 1;
        NTracksInCutPhi[phiBin][3] += coeff[3];
        NTracksInCut2Phi[phiBin][3] += coeff[3] * coeff[3];
        NTracksRegPhi[phiBin][3] += 1;
      }
      fDeDxSorts[1]->Fill(Moment, DeDx);
      fTOFSorts[1]->Fill(Moment, track->GetTOFsignal());
    }

    bool selected_p = false;
    if (Pt < nSigmaBoundary[2] && fabs(nOfSigmasTPC_p) < nSigma && fabs(nOfSigmasTPC_pi) > 3 && fabs(nOfSigmasTPC_K) > 3 && fabs(nOfSigmasTPC_el) > 1 && fabs(nOfSigmasTPC_d) > 3)
      selected_p = true;
    if (Pt > nSigmaBoundary[2] && nSigma_comb_p < nSigma && nSigma_comb_pi > 3 && nSigma_comb_K > 3)
      selected_p = true;
    if (selected_p && Pt > PtCut[2])
    {
      if (Charge < 0)
      {
        NTracksInCutEta[EtaBin][4] += coeff[4];
        NTracksInCut2Eta[EtaBin][4] += coeff[4] * coeff[4];
        NTracksRegEta[EtaBin][4] += 1;
        NTracksInCutPhi[phiBin][4] += coeff[4];
        NTracksInCut2Phi[phiBin][4] += coeff[4] * coeff[4];
        NTracksRegPhi[phiBin][4] += 1;
      }
      if (Charge > 0)
      {
        NTracksInCutEta[EtaBin][5] += coeff[5];
        NTracksInCut2Eta[EtaBin][5] += coeff[5] * coeff[5];
        NTracksRegEta[EtaBin][5] += 1;
        NTracksInCutPhi[phiBin][5] += coeff[5];
        NTracksInCut2Phi[phiBin][5] += coeff[5] * coeff[5];
        NTracksRegPhi[phiBin][5] += 1;
      }
      fDeDxSorts[2]->Fill(Moment, DeDx);
      fTOFSorts[2]->Fill(Moment, track->GetTOFsignal());
    }
    NAcceptedtracks++;
  }
  // end of track loop

  const AliAODVZERO *vzrData = fAOD->GetVZEROData();
  float sumV0ampl = vzrData->GetMTotV0A() + vzrData->GetMTotV0C();
  fHistQAMultTPCvsV0[0]->Fill(multTPC, sumV0ampl);
  fHistQAMultTPCvsESD[0]->Fill(multTPC, multEsd);
  fHistQAMultTrkvsMultTrkTOF[0]->Fill(multTrk, multTrkTOF);
  if (multEsd - 3.38 * multTPC >= 15000 && LargeTPCCut)
  {
    PostData(1, fOutputList);
    return;
  }
  fHistQASPDTrackletsvsV0MCent[2]->Fill(centr, NTrackletsSPD);
  fHistQAMultTPCvsESD[1]->Fill(multTPC, multEsd);
  fHistQAMultTPCvsV0[1]->Fill(multTPC, multEsd);
  fHistQAMultTrkvsMultTrkTOF[1]->Fill(multTrk, multTrkTOF);

  int nAcceptedGenTracks = 0, NAcceptedGenTracksinEta[nEtaClasses] = {0};
  if (IsMC)
  {
    // MC track loop
    AliMCEvent *fMC = eventHandler->MCEvent();
    Int_t nMCTracks(fMC->GetNumberOfTracks());

    for (Int_t i(0); i < nMCTracks; i++)
    {
      AliAODMCParticle *trackMC = (AliAODMCParticle *)fMC->GetTrack(i);
      if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
        continue;
      if (trackMC && trackMC->IsPhysicalPrimary() && fabs(trackMC->Eta()) < 0.8 && trackMC->Pt() >= minP && trackMC->Pt() <= maxP)
      {
        Float_t GenPt = trackMC->Pt();
        Float_t GenPhi = trackMC->Phi();
        Float_t GenEta = trackMC->Eta();
        int GenCharge = trackMC->Charge();

        EtaBin = (0.8 + GenEta) / (1.6 / nEtaClasses);
        // phiBin = GenPhi / (TMath::TwoPi() / nPhiWindows);
        phiBin = (GenPhi + (movePhi * TMath::TwoPi() / nPhiWindows / 2.0)) / (TMath::TwoPi() / nPhiWindows); // move phi on half of bin for stat errors estimation
        if (phiBin == nPhiWindows)
          phiBin = 0;
        NAcceptedGenTracksinEta[EtaBin]++;

        sort = 20;
        if (fabs(trackMC->GetPdgCode()) == 211)
          sort = 0;
        if (fabs(trackMC->GetPdgCode()) == 321)
          sort = 1;
        if (fabs(trackMC->GetPdgCode()) == 2212)
          sort = 2;
        if (fabs(trackMC->GetPdgCode()) == 11)
          sort = 3;
        if (fabs(trackMC->GetPdgCode()) == 1000010020)
          sort = 4;
        nAcceptedGenTracks++;

        if (sort > 3)
          continue;
        if (GenPt > PtCut[sort] && GenCharge != 0)
        {
          GenNParticlesEta[EtaBin][sort * 2 + (GenCharge < 0 ? 0 : 1)]++;
          GenNParticlesPhi[phiBin][sort * 2 + (GenCharge < 0 ? 0 : 1)]++;
        }
      }
    }
  }

  NEvents->AddBinContent(CentrBin + 1);
  NEventsSub[randomSubsample]->AddBinContent(CentrBin + 1);

  // fill all acceptance hists
  Double_t NTracksInAllAcc[nSorts] = {0}, GenNTracksInAllAcc[nSorts] = {0}, RegNTracksInAllAcc[nSorts] = {0};

  for (int iSort = 0; iSort < 6; iSort++)
  {
    int sortNoSign = iSort / 2;
    int skipSort = (2 * (6) - (iSort - 1)) * iSort / 2;
    for (int iEta = 0; iEta < nEtaClasses; iEta++)
    {
      fHistMomentsInAllAccRec[CentrBin]->AddBinContent(sortNoSign + 3 + 1, NTracksInCutEta[iEta][iSort] - NTracksInCut2Eta[iEta][iSort]);
      fHistMomentsInAllAccRecSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 3 + 1, NTracksInCutEta[iEta][iSort] - NTracksInCut2Eta[iEta][iSort]);
      if (sortNoSign == 0)
      {
        fHistMomentsInAllAccRec[CentrBin]->AddBinContent(3 * 2 + 1 + 2 + iSort, NTracksInCutEta[iEta][iSort] - NTracksInCut2Eta[iEta][iSort]);
        fHistMomentsInAllAccRecSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 1 + 2 + iSort, NTracksInCutEta[iEta][iSort] - NTracksInCut2Eta[iEta][iSort]);
      }

      NTracksInAllAcc[iSort] += NTracksInCutEta[iEta][iSort];
      GenNTracksInAllAcc[iSort] += GenNParticlesEta[iEta][iSort];
      RegNTracksInAllAcc[iSort] += NTracksRegEta[iEta][iSort];
    }

    fHistMomentsInAllAccRec[CentrBin]->AddBinContent(sortNoSign + 1, NTracksInAllAcc[iSort]);
    fHistMomentsInAllAccRec[CentrBin]->AddBinContent(sortNoSign + 3 + 1, pow(NTracksInAllAcc[iSort], 2));
    fHistMomentsInAllAccGen[CentrBin]->AddBinContent(sortNoSign + 1, GenNTracksInAllAcc[iSort]);
    fHistMomentsInAllAccGen[CentrBin]->AddBinContent(sortNoSign + 3 + 1, pow(GenNTracksInAllAcc[iSort], 2));
    fHistMomentsInAllAccReg[CentrBin]->AddBinContent(sortNoSign + 1, RegNTracksInAllAcc[iSort]);
    fHistMomentsInAllAccReg[CentrBin]->AddBinContent(sortNoSign + 3 + 1, pow(RegNTracksInAllAcc[iSort], 2));

    fHistMomentsInAllAccRecSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 1, NTracksInAllAcc[iSort]);
    fHistMomentsInAllAccRecSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 3 + 1, pow(NTracksInAllAcc[iSort], 2));
    fHistMomentsInAllAccGenSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 1, GenNTracksInAllAcc[iSort]);
    fHistMomentsInAllAccGenSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 3 + 1, pow(GenNTracksInAllAcc[iSort], 2));
    fHistMomentsInAllAccRegSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 1, RegNTracksInAllAcc[iSort]);
    fHistMomentsInAllAccRegSub[randomSubsample][CentrBin]->AddBinContent(sortNoSign + 3 + 1, pow(RegNTracksInAllAcc[iSort], 2));

    if (sortNoSign == 0)
    {
      fHistMomentsInAllAccRec[CentrBin]->AddBinContent(3 * 2 + 1 + iSort, NTracksInAllAcc[iSort]);
      fHistMomentsInAllAccRec[CentrBin]->AddBinContent(3 * 2 + 2 + 1 + iSort, pow(NTracksInAllAcc[iSort], 2));
      fHistMomentsInAllAccGen[CentrBin]->AddBinContent(3 * 2 + 1 + iSort, GenNTracksInAllAcc[iSort]);
      fHistMomentsInAllAccGen[CentrBin]->AddBinContent(3 * 2 + 2 + 1 + iSort, pow(GenNTracksInAllAcc[iSort], 2));
      fHistMomentsInAllAccReg[CentrBin]->AddBinContent(3 * 2 + 1 + iSort, RegNTracksInAllAcc[iSort]);
      fHistMomentsInAllAccReg[CentrBin]->AddBinContent(3 * 2 + 2 + 1 + iSort, pow(RegNTracksInAllAcc[iSort], 2));

      fHistMomentsInAllAccRecSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 1 + iSort, NTracksInAllAcc[iSort]);
      fHistMomentsInAllAccRecSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 2 + 1 + iSort, pow(NTracksInAllAcc[iSort], 2));
      fHistMomentsInAllAccGenSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 1 + iSort, GenNTracksInAllAcc[iSort]);
      fHistMomentsInAllAccGenSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 2 + 1 + iSort, pow(GenNTracksInAllAcc[iSort], 2));
      fHistMomentsInAllAccRegSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 1 + iSort, RegNTracksInAllAcc[iSort]);
      fHistMomentsInAllAccRegSub[randomSubsample][CentrBin]->AddBinContent(3 * 2 + 2 + 1 + iSort, pow(RegNTracksInAllAcc[iSort], 2));
    }

    // fill hists for moments in eta
    for (int iEta = 0; iEta < nEtaClasses; iEta++)
    {
      if (IsMC)
      {
        fHistMomentsInEtaGen[CentrBin][iSort]->AddBinContent(iEta + 1, GenNParticlesEta[iEta][iSort]);
        fHistMomentsInEtaGenSub[randomSubsample][CentrBin][iSort]->AddBinContent(iEta + 1, GenNParticlesEta[iEta][iSort]);
      }
      fHistMomentsInEtaRec[CentrBin][iSort]->AddBinContent(iEta + 1, NTracksInCutEta[iEta][iSort]);
      fHistMomentsInEtaRecSub[randomSubsample][CentrBin][iSort]->AddBinContent(iEta + 1, NTracksInCutEta[iEta][iSort]);
      fHistMomentsInEtaReg[CentrBin][iSort]->AddBinContent(iEta + 1, NTracksRegEta[iEta][iSort]);
      fHistMomentsInEtaRegSub[randomSubsample][CentrBin][iSort]->AddBinContent(iEta + 1, NTracksRegEta[iEta][iSort]);

      // fHistMomentsInEtaRec[CentrBin][iSort]->AddBinContent(iEta+1+nEtaClasses,NTracksInCutEta[iEta][iSort]-NTracksInCut2Eta[iEta][iSort]);
      // fHistMomentsInEtaSub[randomSubsample][CentrBin][iSort]->AddBinContent(iEta+1+nEtaClasses,NTracksInCutEta[iEta][iSort]-NTracksInCut2Eta[iEta][iSort]);
      for (int jEta = 0; jEta < nEtaClasses; jEta++)
      {
        for (int jSort = iSort; jSort < 6; jSort++)
        {
          if (IsMC)
          {
            fHistCrossMomentsInEtaGen[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, GenNParticlesEta[iEta][iSort] * GenNParticlesEta[jEta][jSort]);
            fHistCrossMomentsInEtaGenSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, GenNParticlesEta[iEta][iSort] * GenNParticlesEta[jEta][jSort]);
          }
          fHistCrossMomentsInEtaRec[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, NTracksInCutEta[iEta][iSort] * NTracksInCutEta[jEta][jSort]);
          fHistCrossMomentsInEtaReg[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, NTracksRegEta[iEta][iSort] * NTracksRegEta[jEta][jSort]);
          fHistCrossMomentsInEtaRecSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, NTracksInCutEta[iEta][iSort] * NTracksInCutEta[jEta][jSort]);
          fHistCrossMomentsInEtaRegSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, NTracksRegEta[iEta][iSort] * NTracksRegEta[jEta][jSort]);

          if (jEta == iEta && jSort == iSort)
          {
            fHistCrossMomentsInEtaRec[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, NTracksInCutEta[iEta][iSort] - NTracksInCut2Eta[iEta][iSort]);
            fHistCrossMomentsInEtaRecSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iEta * nEtaClasses + jEta, NTracksInCutEta[iEta][iSort] - NTracksInCut2Eta[iEta][iSort]);
          }
        }
      }
    }

    // fill for phi:
    for (int iPhi = 0; iPhi < nPhiWindows; iPhi++)
    {
      fHistMomentsInPhiRec[CentrBin][iSort]->AddBinContent(iPhi + 1, NTracksInCutPhi[iPhi][iSort]);
      fHistMomentsInPhiReg[CentrBin][iSort]->AddBinContent(iPhi + 1, NTracksRegPhi[iPhi][iSort]);
      fHistMomentsInPhiGen[CentrBin][iSort]->AddBinContent(iPhi + 1, GenNParticlesPhi[iPhi][iSort]);

      fHistMomentsInPhiRecSub[randomSubsample][CentrBin][iSort]->AddBinContent(iPhi + 1, NTracksInCutPhi[iPhi][iSort]);
      fHistMomentsInPhiRegSub[randomSubsample][CentrBin][iSort]->AddBinContent(iPhi + 1, NTracksRegPhi[iPhi][iSort]);
      fHistMomentsInPhiGenSub[randomSubsample][CentrBin][iSort]->AddBinContent(iPhi + 1, GenNParticlesPhi[iPhi][iSort]);

      // int skipPhi = (2 * (nPhiWindows) - (iPhi - 1)) * iPhi / 2;
      for (int jPhi = 0; jPhi < nPhiWindows; jPhi++)
      {
        for (int jSort = iSort; jSort < 6; jSort++)
        {
          fHistCrossMomentsInPhiRec[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, NTracksInCutPhi[iPhi][iSort] * NTracksInCutPhi[jPhi][jSort]);
          fHistCrossMomentsInPhiReg[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, NTracksRegPhi[iPhi][iSort] * NTracksRegPhi[jPhi][jSort]);
          fHistCrossMomentsInPhiGen[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, GenNParticlesPhi[iPhi][iSort] * GenNParticlesPhi[jPhi][jSort]);

          fHistCrossMomentsInPhiRecSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, NTracksInCutPhi[iPhi][iSort] * NTracksInCutPhi[jPhi][jSort]);
          fHistCrossMomentsInPhiRegSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, NTracksRegPhi[iPhi][iSort] * NTracksRegPhi[jPhi][jSort]);
          fHistCrossMomentsInPhiGenSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, GenNParticlesPhi[iPhi][iSort] * GenNParticlesPhi[jPhi][jSort]);

          if (jPhi == iPhi && jSort == iSort)
          {
            fHistCrossMomentsInPhiRec[CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, NTracksInCutPhi[iPhi][iSort] - NTracksInCut2Phi[iPhi][iSort]);
            fHistCrossMomentsInPhiRecSub[randomSubsample][CentrBin][skipSort + jSort - iSort]->AddBinContent(1 + iPhi * nPhiWindows + jPhi, NTracksInCutPhi[iPhi][iSort] - NTracksInCut2Phi[iPhi][iSort]);
          }
        }
      }
    }
  }

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskParticleYieldRatioCorrelations::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
