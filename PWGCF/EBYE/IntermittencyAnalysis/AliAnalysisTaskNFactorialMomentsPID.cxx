// Implementation of the Intermittency analysis
// for charged particles, two dimensions (eta and phi)
// Data and MC
// Contact: ramni.gupta@cern.ch
// Contributors: R.Gupta, S.Sharma, S.K.Malik

#include "AliAnalysisTaskNFactorialMomentsPID.h"

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultSelection.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TObjArray.h"
#include "TGrid.h"

class AliAnalysisTaskNFactorialMomentsPID;
using namespace std;

ClassImp(AliAnalysisTaskNFactorialMomentsPID)
  AliAnalysisTaskNFactorialMomentsPID::AliAnalysisTaskNFactorialMomentsPID()
  : AliAnalysisTaskSE(), fAOD(0), fMCEvent(0), fHistList(0), fNtupleListBin1(0), fNtupleListBin2(0), fNtupleListBin3(0), fNtupleListBin4(0), fWCNtupleListBin1(0), fWCNtupleListBin2(0), fWCNtupleListBin3(0), fWCNtupleListBin4(0), zeromult(0), fEventCuts(0), fPIDResponse(0), fPIDCombined(0), mcHeader{ nullptr }, fNtupleList{ nullptr }, fNtupleListCorr{ nullptr }, fNtupleListGen{ nullptr }, fHistQAEta{ nullptr }, fHistQAPhi{ nullptr }, fHistQAVx(0), fHistQAVy(0), counter(0), fHistQAVz(0), fPtBin{ nullptr }, fEtaBin{ nullptr }, fPhiBin{ nullptr }, fMultBin{ nullptr }, fPtBinGen{ nullptr }, fEtaBinGen{ nullptr }, fPhiBinGen{ nullptr }, fMultBinGen{ nullptr }, fHistQAPID{ nullptr }, fHistPDG{ nullptr }, fHistDCAxy{ nullptr }, fHistDCAz{ nullptr }, fHistnITScls{ nullptr }, fHistnTPCcls{ nullptr }, fHistnTPCcrossedrows{ nullptr }, fHistnchi2ITScls{ nullptr }, fHistnchi2TPCcls{ nullptr }, fHistDCAxypT{ nullptr }, fHistDCAzpT{ nullptr }, fEtaPhiBin{ nullptr }, fntpMBin{ nullptr }, fntpMBinCorr{ nullptr }, fntpMBinGen{ nullptr }, fEtaPhiBinGen{ nullptr }, fHistNShCls{ nullptr }, fHistNShClsFra{ nullptr }, fHistNShClsFravspt{ nullptr }, fHistNFoundClsFra{ nullptr }, fHistNFcls{ nullptr }, fEventCounter(0), fTrackCounter(0), fHistQACent(0), ptbin1(0), ptbin2(0), ptbin3(0), ptbin4(0), mfield(0) {}
AliAnalysisTaskNFactorialMomentsPID::AliAnalysisTaskNFactorialMomentsPID(
  const char* name)
  : AliAnalysisTaskSE(name), fAOD(0), fMCEvent(0), fHistList(0), fNtupleListBin1(0), fNtupleListBin2(0), fNtupleListBin3(0), fNtupleListBin4(0), fWCNtupleListBin1(0), fWCNtupleListBin2(0), fWCNtupleListBin3(0), fWCNtupleListBin4(0), zeromult(0), fEventCuts(0), fPIDResponse(0), fPIDCombined(0), mcHeader{ nullptr }, fNtupleList{ nullptr }, fNtupleListCorr{ nullptr }, fNtupleListGen{ nullptr }, fHistQAEta{ nullptr }, fHistQAPhi{ nullptr }, fHistQAVx(0), fHistQAVy(0), counter(0), fHistQAVz(0), fPtBin{ nullptr }, fEtaBin{ nullptr }, fPhiBin{ nullptr }, fMultBin{ nullptr }, fPtBinGen{ nullptr }, fEtaBinGen{ nullptr }, fPhiBinGen{ nullptr }, fMultBinGen{ nullptr }, fHistQAPID{ nullptr }, fHistPDG{ nullptr }, fHistDCAxy{ nullptr }, fHistDCAz{ nullptr }, fHistnITScls{ nullptr }, fHistnTPCcls{ nullptr }, fHistnTPCcrossedrows{ nullptr }, fHistnchi2ITScls{ nullptr }, fHistnchi2TPCcls{ nullptr }, fHistDCAxypT{ nullptr }, fHistDCAzpT{ nullptr }, fEtaPhiBin{ nullptr }, fEtaPhiBinGen{ nullptr }, fntpMBin{ nullptr }, fntpMBinCorr{ nullptr }, fntpMBinGen{ nullptr }, fHistNShCls{ nullptr }, fHistNShClsFra{ nullptr }, fHistNShClsFravspt{ nullptr }, fHistNFoundClsFra{ nullptr }, fHistNFcls{ nullptr }, fEventCounter(0), fTrackCounter(0), fHistQACent(0), ptbin1(0), ptbin2(0), ptbin3(0), ptbin4(0), mfield(0)
{
  // Constructor
  Info("AliAnalysisTaskNFactorialMomentsPID", "Specific Constructor");

  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    nsigmaTPC[i] = -999;
    nsigmaTOF[i] = -999;
  }

  for (Int_t j = 0; j < 40; j++) {
    fhMapEtaPhiBin1M[j] = NULL;
    fhMapEtaPhiBin2M[j] = NULL;
    fhMapEtaPhiBin3M[j] = NULL;
    fhMapEtaPhiBin4M[j] = NULL;

    fWCntpMBin1[j] = NULL;
    fWCntpMBin2[j] = NULL;
    fWCntpMBin3[j] = NULL;
    fWCntpMBin4[j] = NULL;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  if (flagMC) {
    DefineOutput(8, TList::Class());
    DefineOutput(9, TList::Class());
    DefineOutput(10, TList::Class());
    DefineOutput(11, TList::Class());
  }
}

AliAnalysisTaskNFactorialMomentsPID::~AliAnalysisTaskNFactorialMomentsPID()
{
  // destructor
  if (fHistList) {
    delete fHistList;
  }
  if (fPIDResponse) {
    delete fPIDResponse;
    fPIDResponse = 0x0;
  }
  if (fPIDCombined) {
    delete fPIDCombined;
    fPIDCombined = 0x0;
  }
}

void AliAnalysisTaskNFactorialMomentsPID::UserCreateOutputObjects()
{
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  fNtupleListBin1 = new TList(); // NTuple to Store Reconstructed NFM
  fNtupleListBin2 = new TList();
  fNtupleListBin3 = new TList();
  fNtupleListBin4 = new TList();

  fNtupleListBin1->SetOwner(kTRUE);
  fNtupleListBin2->SetOwner(kTRUE);
  fNtupleListBin3->SetOwner(kTRUE);
  fNtupleListBin4->SetOwner(kTRUE);

  TString NtWCnameBin1, NtWCtitleBin1; // Nt--> Ntuple Name, Nttitle-->Ntuple title  ; for without correction reconstructed
  TString NtWCnameBin2, NtWCtitleBin2;
  TString NtWCnameBin3, NtWCtitleBin3;
  TString NtWCnameBin4, NtWCtitleBin4;

  fWCNtupleListBin1 = new TList(); // NTuple to Store Reconstructed NFM
  fWCNtupleListBin2 = new TList();
  fWCNtupleListBin3 = new TList();
  fWCNtupleListBin4 = new TList();

  for (Int_t i = 0; i < mPtBins; ++i) {
    fNtupleList[i] = new TList();
    fNtupleList[i]->SetOwner(kTRUE);

    fNtupleListCorr[i] = new TList();
    fNtupleListCorr[i]->SetOwner(kTRUE);

    if (flagMC) {
      fNtupleListGen[i] = new TList();
      fNtupleListGen[i]->SetOwner(kTRUE);
    }
  }
  // Eventcounter hist defined
  fEventCounter =
    new TH1D("fEventCounter", "histo to keep track", 10, 0.5, 10.5);
  fEventCounter->GetXaxis()->SetBinLabel(1, "Events before cuts");
  fEventCounter->GetXaxis()->SetBinLabel(2, "Pileup cut");
  fEventCounter->GetXaxis()->SetBinLabel(3, "Phy. Sel. cut");
  fEventCounter->GetXaxis()->SetBinLabel(4, "Events with a proper vertex");
  fEventCounter->GetXaxis()->SetBinLabel(5, "Events with 0-5% Centrality");
  fEventCounter->GetXaxis()->SetBinLabel(6, "BField");
  fEventCounter->GetXaxis()->SetBinLabel(7, "Events Analyzed");
  fHistList->Add(fEventCounter);

  // Trackcounter hist defined
  fTrackCounter =
    new TH1D("fTrackCounter", "histo to keep track of tracks", 19, 0.5, 19.5);
  fTrackCounter->GetXaxis()->SetBinLabel(1, "Total tracks");
  fTrackCounter->GetXaxis()->SetBinLabel(2, Form("FilterBit%d", filterBit));
  fTrackCounter->GetXaxis()->SetBinLabel(3, "PIDs");
  fTrackCounter->GetXaxis()->SetBinLabel(4, "Loose cuts");
  fTrackCounter->GetXaxis()->SetBinLabel(5, "Systematic cuts");
  fTrackCounter->GetXaxis()->SetBinLabel(6, "DCAs");
  fTrackCounter->GetXaxis()->SetBinLabel(7, "nSharedCls/nCls");
  fTrackCounter->GetXaxis()->SetBinLabel(8, "nSharedCls/nCrRows");
  fTrackCounter->GetXaxis()->SetBinLabel(9, "nCrRows/nFindableCls");
  for (Int_t i = 0; i < mPtBins; ++i) {
    fTrackCounter->GetXaxis()->SetBinLabel(10 + i, Form("p_{T}_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]));
  }
  fHistList->Add(fTrackCounter);

  fHistQACent = new TH1F("fHistQACent", "Centrality Distribution", 10,
                         minCent - 2, maxCent + 2);

  fHistQAVx =
    new TH1F("fHistQAVx", "Primary vertex distribution -x coordinate;V_{x}",
             1500, -0.3, 0.3);
  fHistQAVy =
    new TH1F("fHistQAVy", "Primary vertex distribution -y coordinate;V_{y}",
             1500, -0.3, 0.3);
  fHistQAVz =
    new TH1F("fHistQAVz", "Primary vertex distribution -z coordinate;V_{z}",
             1500, -20.0, 20.0);
  fHistList->Add(fHistQACent);
  fHistList->Add(fHistQAVx);
  fHistList->Add(fHistQAVy);
  fHistList->Add(fHistQAVz);

  TString name = flagMC ? "Re" : "";

  fHistQAEta[0] = new TH1F(Form("fHistQAEta%s", name.Data()), "Eta distribution of tracks;#eta;Counts", 1000, -1.0, 1.0);
  fHistQAPhi[0] = new TH1F(Form("fHistQAPhi%s", name.Data()), "Phi distribution of tracks;#phi;Counts", 1000, 0.0, 6.5);
  fHistList->Add(fHistQAEta[0]);
  fHistList->Add(fHistQAPhi[0]);
  if (flagMC) {
    fHistQAEta[1] = new TH1F(Form("fHistQAEta%s", name.Data()), "Eta distribution of tracks;#eta;Counts", 1000, -1.0, 1.0);
    fHistQAPhi[1] = new TH1F(Form("fHistQAPhi%s", name.Data()), "Phi distribution of tracks;#phi;Counts", 1000, 0.0, 6.5);
    fHistList->Add(fHistQAEta[1]);
    fHistList->Add(fHistQAPhi[1]);
  }
  fHistDCAxy = new TH1F("fHistDCAxy", "DCAxy distribution;DCAxy;Counts", 900, -4.5, 4.5);
  fHistDCAz = new TH1F("fHistDCAz", "DCAz distribution;DCAz;Counts", 900, -4.5, 4.5);
  fHistDCAxypT = new TH2F("fHistDCAxypT", "DCAxy vs pT;p_{T};dca_{xy}", 50, 0, 5.5, 200, -2.0, 2.0);
  fHistDCAzpT = new TH2F("fHistDCAzpT", "DCAz vs pT;p_{T};dca_{z}", 50, 0, 5.5, 200, -2.0, 2.0);
  fHistnITScls = new TH1F("fHistnITScls", "ITS cluster distribution;ITS cluster;Counts", 10, -0.5, 10.5);
  fHistnTPCcls = new TH1F("fHistnTPCcls", "TPC cluster distribution;TPC cluster;Counts", 200, 0, 200);
  fHistnchi2ITScls = new TH1F("fHistnchi2ITScls", "ITSchi2 cluster distribution;ITSchi2 cluster;Counts", 100, -0.5, 100.5);
  fHistnchi2TPCcls = new TH1F("fHistnchi2TPCcls", "TPCchi2 cluster distribution;TPCchi2 cluster;Counts", 50, 0, 50);
  fHistnTPCcrossedrows = new TH1F("fHistnTPCcrossedrows", "TPC crossed rows distribution;TPC crossed rows;Counts", 200, 0, 200);
  fHistList->Add(fHistDCAxy);
  fHistList->Add(fHistDCAz);
  fHistList->Add(fHistDCAxypT);
  fHistList->Add(fHistDCAzpT);
  fHistList->Add(fHistnITScls);
  fHistList->Add(fHistnTPCcls);
  fHistList->Add(fHistnTPCcrossedrows);
  fHistList->Add(fHistnchi2ITScls);
  fHistList->Add(fHistnchi2TPCcls);

  for (Int_t i = 0; i < 2; ++i) {
    std::array<TString, 2> j = { "before", "after" };
    fHistNShCls[i] = new TH1F(Form("fHistNShCls%s", j[i].Data()), Form("TPC shared cluster distribution %s", j[i].Data()), 200, 0, 200);
    fHistNShClsFra[i] = new TH2F(Form("fHistNShClsFra%s", j[i].Data()), Form("TPC shared cluster fraction %s;sharedcls/ncls;sharedcls/ncrrows", j[i].Data()), 100, 0, 1.0, 100, 0, 1.0);
    fHistNFoundClsFra[i] = new TH2F(Form("fHistNFoundClsFra%s", j[i].Data()), Form("TPC found cluster fraction %s;sharedcls/ncls;ncrrows/findablecls", j[i].Data()), 200, 0, 2, 200, 0, 2);
    fHistNFcls[i] = new TH1F(Form("fHistNFcls%s", j[i].Data()), Form("TPC findable cluster distribution %s;ncrrows/findablecls;counts", j[i].Data()), 40, 0, 4);
    fHistList->Add(fHistNShCls[i]);
    fHistList->Add(fHistNShClsFra[i]);
    fHistList->Add(fHistNFoundClsFra[i]);
    fHistList->Add(fHistNFcls[i]);
  }
  fHistNShClsFravspt = new TH2F("fHistNShClsFraVsPt", "TPC nshared clusters vs #it{p_{T}} before", 500, 0, 5, 400, 0.02, 1);
  fHistList->Add(fHistNShClsFravspt);

  for (Int_t i = 0; i < mPtBins; ++i) {
    fPtBin[i] = new TH1F(Form("fPtBin%s_%.2f-%.2f", name.Data(), ptarray[2 * i], ptarray[2 * i + 1]), Form("Pt distribution of tracks for PtBin%d;#p_{T};Counts", i + 1), 1000, 0.0, 5.0);
    fEtaBin[i] = new TH1F(Form("fEtaBin%s_%.2f-%.2f", name.Data(), ptarray[2 * i], ptarray[2 * i + 1]), Form("Eta distribution of tracks for PtBin%d;#eta;Counts", i + 1), 1000, -1.0, 1.0);
    fPhiBin[i] = new TH1F(Form("fPhiBin%s_%.2f-%.2f", name.Data(), ptarray[2 * i], ptarray[2 * i + 1]), Form("Phi distribution of tracks for PtBin%d;#phi;Counts", i + 1), 1000, 0.0, 6.5);
    fMultBin[i] = new TH1F(Form("fMultBin%s_%.2f-%.2f", name.Data(), ptarray[2 * i], ptarray[2 * i + 1]), Form("Multiplicity distribution of tracks for PtBin%d;Multiplicity;Counts", i + 1), 10000, 0.0, 10000.0);
    fHistList->Add(fPtBin[i]);
    fHistList->Add(fEtaBin[i]);
    fHistList->Add(fPhiBin[i]);
    fHistList->Add(fMultBin[i]);
    if (flagMC) {
      fPtBinGen[i] = new TH1F(Form("fPtBinGen_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]), Form("Pt distribution of Gen tracks for PtBin%d;#p_{T};Counts", i + 1), 1000, 0.0, 5.0);
      fEtaBinGen[i] = new TH1F(Form("fEtaBinGen_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]), Form("Eta distribution of Gen tracks for PtBin%d;#eta;Counts", i + 1), 1000, -1.0, 1.0);
      fPhiBinGen[i] = new TH1F(Form("fPhiBinGen_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]), Form("Phi distribution of Gen tracks for PtBin%d;#phi;Counts", i + 1), 1000, 0.0, 6.5);
      fMultBinGen[i] = new TH1F(Form("fMultBinGen_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]), Form("Multiplicity distribution of Gen tracks for PtBin%d;Multiplicity;Counts", i + 1), 1000, 0.0, 1000.0);
      fHistList->Add(fPtBinGen[i]);
      fHistList->Add(fEtaBinGen[i]);
      fHistList->Add(fPhiBinGen[i]);
      fHistList->Add(fMultBinGen[i]);
    }
    for (Int_t j = 0; j < maxBins; ++j) {
      // Eta-Phi distributions
      fEtaPhiBin[i][j] = new TH2D(Form("fEtaPhiBin%s_%.2f-%.2f_%d", name.Data(), ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mMBin2[j], 0.0, 6.30);
      if (flagSelfAff) {
        fEtaPhiBin[i][j] = new TH2D(Form("fEtaPhiBin%s_%.2f-%.2f_%d", name.Data(), ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mNBin2[j], 0.0, 6.30);
      }
      if (flagMC) {
        fEtaPhiBinGen[i][j] = new TH2D(Form("fEtaPhiBinGen_%.2f-%.2f_%d", ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of Gen tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mMBin2[j], 0.0, 6.30);
        if (flagSelfAff) {
          fEtaPhiBinGen[i][j] = new TH2D(Form("fEtaPhiBinGen_%.2f-%.2f_%d", ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of Gen tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mNBin2[j], 0.0, 6.30);
        }
        fHistList->Add(fEtaPhiBinGen[i][j]);
      }
      fHistList->Add(fEtaPhiBin[i][j]);
      // Tuples for each bin
      fntpMBin[i][j] = new TNtuple(Form("fntpMBin%d_%d%s", i + 1, j + 1, name.Data()), Form("fntpMBin%d_%d", i + 1, j + 1), "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
      fNtupleList[i]->Add(fntpMBin[i][j]);

      fntpMBinCorr[i][j] = new TNtuple(Form("fntpMBin%d_%dCorr", i + 1, j + 1), Form("fntpMBin%d_%dCorr", i + 1, j + 1), "Mbin:WCAv_bincontent:WCFq2e:WCFq3e:WCFq4e:WCFq5e");
      fNtupleList[i]->Add(fntpMBinCorr[i][j]);
      if (flagMC) {
        fntpMBinGen[i][j] = new TNtuple(Form("fntpMBin%d_%dGen", i + 1, j + 1), Form("fntpMBin%d_%dGen", i + 1, j + 1), "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
        fNtupleListGen[i]->Add(fntpMBinGen[i][j]);
      }
    }
  }
  fHistQAPID[0] = new TH2D("TPCdEdx", "TPCdEdx vs momentum;#it{p};TPCdEdx", 100, 0.0, 5.0, 250, 0.0, 500.0);
  fHistQAPID[1] = new TH2D("TPC_Nsigma_Proton", "TPCNsigma vs momentum;#it{p};TPCNsigma", 1000, -10, +10, 500, -5, 5);
  fHistQAPID[2] = new TH2D("TOF_NsigmaProton", "TOFNsigmaProton vs momentum;#it{p};TOFNsigmaProton", 1000, -10.0, 10.0, 500, -5.0, 5.0);
  fHistQAPID[3] = new TH2D("TPCTOFProton", "TPCTOFProton vs momentum;#it{p};TPCTOFProton", 1000, -10, +10, 500, -5, 5);
  fHistQAPID[4] = new TH2D("TOF signal", "momentum vs beta;#it{p};#beta", 1000, 0.5, 5.0, 1000, 0.4, 1.2);
  fHistQAPID[5] = new TH2D("nSigmaTPC  vs nSigmaTOF", "nSigmaTPC  vs nSigmaTOF;nSigmaTPC;nSigmaTOF", 1000, -5.0, 5.0, 1000, -5.0, 5.0);
  fHistQAPID[6] = new TH2D("TPCdEdxAfter", "TPCdEdx vs momentum;#it{p};TPCdEdx", 100, 0.0, 10.0, 250, 0.0, 250.0);
  fHistQAPID[7] = new TH2D("TPC_NsigmaProtonAfter", "TPCNsigmaProton vs momentum;#it{p};TPCNsigmaProton", 1000, -10, +10, 500, -5, 5);
  fHistQAPID[8] = new TH2D("TOF_NsigmaProtonAfter", "TOFNsigmaProton vs momentum;#it{p};TOFNsigmaProton", 1000, -10, +10, 500, -5, 5);
  fHistQAPID[9] = new TH2D("TPCTOF_ProtonAfter", "TPCTOFProton vs momentum;#it{p};TPCTOFProton", 1000, -10, +10, 500, -5, 5);
  fHistQAPID[10] = new TH2D("TOF signalAfter", "momentum vs beta;#it{p};#beta", 1000, 0.5, 5.0, 1000, 0.4, 1.2);
  fHistQAPID[11] = new TH2D("nSigmaTPC  vs nSigmaTOFAfter", "nSigmaTPC  vs nSigmaTOF;nSigmaTPC;nSigmaTOF", 1000, -5.0, 5.0, 1000, -5.0, 5.0);
  fHistQAPID[12] = new TH2D("TPCNsigmaProtonAfter", "TPCNsigmaProton vs momentum;#it{p};TPCNsigmaProton", 100, 0.0, 2.0, 100, -5.0, 5.0);
  fHistPDG[0] = new TH1D("PDGHist", "PDG code of the generated particles;PDG code;Counts", 1000, -500, 500);
  fHistPDG[1] = new TH1D("ElPosPDGHist", "PDG code of the generated electrons and positrons;PDG code;Counts", 100, -100, 100);
  for (Int_t i = 0; i < 13; ++i) {
    fHistList->Add(fHistQAPID[i]);
  }
  for (Int_t i = 0; i < 2; ++i) {
    fHistList->Add(fHistPDG[i]);
  }
  if (!flagMC) {
    fEventCuts.AddQAplotsToList(fHistList, kTRUE);
  }

  DataPosting();
}

// Execution - called per event:
void AliAnalysisTaskNFactorialMomentsPID::UserExec(Option_t*)
{
  fEventCounter->Fill(1);
  AliAODInputHandler* eventHandler = dynamic_cast<AliAODInputHandler*>(
    AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  fPIDResponse =
    dynamic_cast<AliAODInputHandler*>(
      AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())
      ->GetPIDResponse();

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD)
    ::Fatal("AliAnalysisTaskNFactorialMomentsPID::UserExec",
            "No AOD event found!");

  if (flagMC) {
    fMCEvent = MCEvent();
    if (!fMCEvent)
      ::Fatal("AliAnalysisTaskNFactorialMomentsPID::UserExec",
              "No MC event found!");
    mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(
      AliAODMCHeader::StdBranchName());
    if (!mcHeader)
      ::Fatal("AliAnalysisTaskNFactorialMomentsPID::UserExec",
              "No MC header found!");
  } else {
    if (!fPIDResponse)
      ::Fatal("AliAnalysisTaskNFactorialMomentsPID::UserExec",
              "No PIDResponse found!");
  }

  // fEventCuts.SetupLHC15o();
  if (flagPileup && !flagMC) {
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    Bool_t isAcceptedEvent = fEventCuts.AcceptEvent(fAOD); // for QA
    if (!isAcceptedEvent)
      return;
  } else if (flagPileup && flagMC) {
    Bool_t isPileupinGenEvent =
      AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader, fGenName);
    if (isPileupinGenEvent)
      return;
  }

  fEventCounter->Fill(2);

  if (!flagMC) {
    UInt_t fSelectMask = fInputHandler->IsEventSelected();
    if (!(fSelectMask & AliVEvent::kINT7))
      return;
  }
  fEventCounter->Fill(3);

  AliAODVertex* vertex = fAOD->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1)
    return;
  Float_t lvx = (Float_t)vertex->GetX();
  Float_t lvy = (Float_t)vertex->GetY();
  Float_t lvz = (Float_t)vertex->GetZ();

  if ((fabs(lvx) < fVxMax) && (fabs(lvy) < fVyMax) && (fabs(lvz) < fVzMax)) {
    fHistQAVx->Fill(lvx);
    fHistQAVy->Fill(lvy);
    fHistQAVz->Fill(lvz);
  } else
    return;
  fEventCounter->Fill(4);

  // Centrality Selection
  Double_t lMultiPercentile = -1;
  // Centrality slecetion for LHC15o
  AliMultSelection* MultSelection =
    (AliMultSelection*)fAOD->FindListObject("MultSelection");
  lMultiPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  if (lMultiPercentile < minCent || lMultiPercentile >= maxCent)
    return;

  fEventCounter->Fill(5);
  fHistQACent->Fill(lMultiPercentile);

  mfield = fAOD->GetMagneticField();

  if (fBfield == -1 && mfield > 0)
    return;
  if (fBfield == 1 && mfield < 0)
    return;

  fEventCounter->Fill(6);

  // Reset Histos
  ResetHistograms();
  fEventCounter->Fill(7);

  zeromult = kFALSE;

  // Filling Track Information
  FillTrackInfo();
  if (zeromult)
    return;
  CalculateNFMs(fEtaPhiBin, kFALSE);

  if (flagMC) {
    FillMCTrackInfo(); // Fill Generated Information
    CalculateNFMs(fEtaPhiBinGen, kTRUE);
  }

  DataPosting();
}

/*________________________________________________________________________
      Fill Track Information for NFM Calculation (Main Track Loop)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMomentsPID::FillTrackInfo()
{
  // lines added for efficiency maps
  Int_t test = 0;

  if (!gGrid) {
        TGrid::Connect("alien://");
    }

  TFile* f;
  if (test == 0) {
   f =TFile::Open("alien:///alice/cern.ch/user/f/fhaider/TwoDimEtaPhiEffMap_HIJ276_FB128NPCent0to5.root");
    if (!f) {
      AliFatal("Track efficiency histogram file is not open");
      return;
    }
  } else {
    f = TFile::Open("../Efficiency_LHC20e3a/TwoDEtaPhiEfficiencyMapHIJ276_FB768Cent0to5.root");
    if (!f) {
      AliFatal("Track efficiency histogram file is not open");
      return;
    }
  }

  //====================================================================================================
  // READING THE EFFICIENCY MAP FILE
  //====================================================================================================
  for (Int_t ik = 1; ik < 41; ik++) {
    TString RSparseNameBin1 = "f2DimRecBin1M";
    RSparseNameBin1 += ik;
    RSparseNameBin1 += "_proj_0_1";
    fhMapEtaPhiBin1M[ik - 1] = (TH2D*)f->Get(RSparseNameBin1);
    if (!fhMapEtaPhiBin1M[ik - 1])
      cout << "error" << endl;
    // Bin 2 Map
    TString RSparseNameBin2 = "f2DimRecBin2M";
    RSparseNameBin2 += ik;
    RSparseNameBin2 += "_proj_0_1";
    fhMapEtaPhiBin2M[ik - 1] = (TH2D*)f->Get(RSparseNameBin2);
    if (!fhMapEtaPhiBin2M[ik - 1])
      cout << "error" << endl;
    // Bin 3 Map
    TString RSparseNameBin3 = "f2DimRecBin3M";
    RSparseNameBin3 += ik;
    RSparseNameBin3 += "_proj_0_1";
    fhMapEtaPhiBin3M[ik - 1] = (TH2D*)f->Get(RSparseNameBin3);
    if (!fhMapEtaPhiBin3M[ik - 1])
      cout << "error" << endl;
    // Bin 4 Map
    TString RSparseNameBin4 = "f2DimRecBin4M";
    RSparseNameBin4 += ik;
    RSparseNameBin4 += "_proj_0_1";
    fhMapEtaPhiBin4M[ik - 1] = (TH2D*)f->Get(RSparseNameBin4);
    if (!fhMapEtaPhiBin4M[ik - 1])
      cout << "error" << endl;
  }

  Int_t weight = 0;

  Int_t counterBin[maxPtBins] = { 0 };
  Int_t nTracks(fAOD->GetNumberOfTracks());
  Float_t dpstar, deta, dphi;
  Float_t dcaXY, dcaZ;
  Float_t dEdx;

  // Track Loop:
  for (Int_t i = 0; i < nTracks; i++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    fTrackCounter->Fill(1);
    if (!(track) || !(track->TestFilterBit(filterBit)))
      continue;
    counter++;
    fTrackCounter->Fill(2);

    Float_t pt = track->Pt();
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Int_t charge = track->Charge();
    Int_t id = track->GetID();

    if (flagMC) {
      TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(
        fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (AODMCTrackArray == NULL)
        //  return kTRUE;
        continue;

      AliAODMCParticle* particle =
        (AliAODMCParticle*)AODMCTrackArray->At(TMath::Abs(track->GetLabel()));
      if (!particle)
        continue;
    }

    Double_t fTPCnSigmaElectron = -999;
    Double_t fTPCnSigmaProton = -999;
    Double_t fTPCnSigmaPion = -999;
    Double_t fTPCnSigmaKaon = -999;

    Double_t fTOFnSigmaPion = -999;
    Double_t fTOFnSigmaProton = -999;
    Double_t fTOFnSigmaKaon = -999;

    Double_t length = -999., beta = -999., tofTime = -999., tof = -999.;
    Double_t c = TMath::C() * 1.E-9; // m/n

    if (!flagMC) {
      fPIDCombined = new AliPIDCombined();
      fPIDCombined->SetDefaultTPCPriors();
    }

    if (filterBit == 128) {
      dEdx = track->GetTPCsignal();
    } else {
      dEdx = track->GetDetPid()->GetTPCsignal();
    }

    tofTime = track->GetTOFsignal();

    length = track->GetIntegratedLength();
    tof = tofTime * 1E-3; // ns
    if ((tof < 0) || (length < 0))
      continue;
    length = length * 0.01; // in meters
    beta = length / (tof * c);

    // TPC nsigma
    fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    fTPCnSigmaElectron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

    // TOF nsigma
    fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

    if ((fTPCnSigmaPion == -999) || (fTPCnSigmaProton == -999) || (fTPCnSigmaKaon == -999) || (fTPCnSigmaElectron == -999))
      continue;
    if ((fTOFnSigmaPion == -999) || (fTOFnSigmaProton == -999) || (fTOFnSigmaKaon == -999))
      continue;

    Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
    Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
    Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

    fHistQAPID[0]->Fill(track->Pt(), dEdx);

    fHistQAPID[1]->Fill(track->Pt(), fTPCnSigmaProton);
    fHistQAPID[2]->Fill(track->Pt(), fTOFnSigmaProton);

    fHistQAPID[3]->Fill(track->Pt(), fTPCplusTOFnSigmaProton);
    fHistQAPID[4]->Fill(track->Pt(), beta);
    fHistQAPID[5]->Fill(fTPCnSigmaProton, fTOFnSigmaProton);

    if ((TMath::Abs(fTPCnSigmaElectron < 3.0)) && (track->Pt() < 0.6))
      continue;

    if ((track->Pt() < 0.6) && (TMath::Abs(fTPCnSigmaProton) >= 3))
      continue;

    if (track->Pt() < 0.6) {
      fHistQAPID[6]->Fill(track->Pt(), dEdx);
      fHistQAPID[7]->Fill(track->Pt(), fTPCnSigmaProton);
      fHistQAPID[8]->Fill(track->Pt(), fTOFnSigmaProton);
      fHistQAPID[9]->Fill(track->Pt(), fTPCplusTOFnSigmaProton);
      fHistQAPID[10]->Fill(track->Pt(), beta);
      fHistQAPID[11]->Fill(fTPCnSigmaProton, fTOFnSigmaProton);
    }

    Int_t flag = 0;
    if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
      flag += 1;
    if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
      flag += 1;
    if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
      flag += 1;
    if (flag > 1)
      continue;

    if ((track->Pt() > 0.6) && (fTPCplusTOFnSigmaProton > fTPCplusTOFnSigmaKaon))
      continue;
    if ((track->Pt() > 0.6) && (fTPCplusTOFnSigmaProton > fTPCplusTOFnSigmaPion))
      continue;

    if ((track->Pt() > 0.6) && (fTPCplusTOFnSigmaProton >= 3))
      continue;

    if (track->Pt() > 0.6) {
      fHistQAPID[6]->Fill(track->Pt(), dEdx);
      fHistQAPID[7]->Fill(track->Pt(), fTPCnSigmaProton);
      fHistQAPID[8]->Fill(track->Pt(), fTOFnSigmaProton);
      fHistQAPID[9]->Fill(track->Pt(), fTPCplusTOFnSigmaProton);
      fHistQAPID[10]->Fill(track->Pt(), beta);
      fHistQAPID[11]->Fill(fTPCnSigmaProton, fTOFnSigmaProton);
    }

    fTrackCounter->Fill(3);

    if ((eta < minEta) || (eta > maxEta) || (fabs(pt) < 0.2) || charge == 0)
      continue;

    // calculate if pt is in the pt bin range
    GetPtBin(pt);
    fTrackCounter->Fill(4);

    Bool_t mSkipTracks = kFALSE;
    Bool_t mSharedTrack = kFALSE;
    Float_t nSharedCls = track->GetTPCnclsS();
    Float_t nCls = track->GetTPCNcls();
    Float_t nCrossedRows = track->GetTPCCrossedRows();
    Float_t nFindableCls = track->GetTPCNclsF();

    if ((fITSCls > 0.0) && ((track->GetITSchi2() / track->GetITSNcls()) > fITSCls))
      continue;

    if ((fTPCCls > 0.0) && (track->GetTPCNcls() < fTPCCls))
      continue;

    if ((fTPCRows > 0.0) && (track->GetTPCNCrossedRows() < fTPCRows))
      continue;

    fTrackCounter->Fill(5);

    fHistnITScls->Fill(track->GetITSNcls());
    fHistnTPCcls->Fill(track->GetTPCNcls());
    fHistnchi2ITScls->Fill((track->GetITSchi2() / track->GetITSNcls()));
    fHistnchi2TPCcls->Fill((track->GetTPCchi2() / track->GetTPCNcls()));
    fHistnTPCcrossedrows->Fill(track->GetTPCNCrossedRows());

    track->GetImpactParameters(dcaXY, dcaZ);

    if ((fDCAxyMax > 0.0) && (fabs(dcaXY) > (0.0208 + 0.04 / TMath::Power(pt, 1.01))))
      continue;
    if ((fDCAzMax > 0.0) && (fabs(dcaZ) > fDCAzMax))
      continue;
    fHistDCAxy->Fill(dcaXY);
    fHistDCAz->Fill(dcaZ);
    fHistDCAxypT->Fill(pt, dcaXY);
    fHistDCAzpT->Fill(pt, dcaZ);
    fTrackCounter->Fill(6);
    fHistNShCls[0]->Fill(nSharedCls);
    fHistNFcls[0]->Fill(nCrossedRows / nFindableCls);
    fHistNShClsFra[0]->Fill(nSharedCls / nCls, nSharedCls / nCrossedRows);
    fHistNFoundClsFra[0]->Fill(nSharedCls / nCls, nCrossedRows / nFindableCls);

    if ((fSharedClsMax > 0) && ((nSharedCls / nCls) > fSharedClsMax))
      continue;
    fTrackCounter->Fill(7);

    if ((fSharedRowsMax > 0) && ((nSharedCls / nCrossedRows) > fSharedRowsMax))
      continue;
    fTrackCounter->Fill(8);

    if ((fFindableClsMin > 0) && ((nCrossedRows / nFindableCls) < fFindableClsMin))
      continue;
    fTrackCounter->Fill(9);

    fHistNShCls[1]->Fill(nSharedCls);
    fHistNFcls[1]->Fill(nCrossedRows / nFindableCls);
    fHistNShClsFra[1]->Fill(nSharedCls / nCls, nSharedCls / nCrossedRows);
    fHistNFoundClsFra[1]->Fill(nSharedCls / nCls, nCrossedRows / nFindableCls);
    fHistQAEta[0]->Fill(eta);
    fHistQAPhi[0]->Fill(phi);
    fHistNShClsFravspt->Fill(pt, nSharedCls / nCrossedRows);

    for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
      if (ptbin[iPt]) {
        fTrackCounter->Fill(10 + iPt);
        counterBin[iPt]++;
        fPtBin[iPt]->Fill(pt);
        fEtaBin[iPt]->Fill(eta);
        fPhiBin[iPt]->Fill(phi);
        for (Int_t k = 0; k < maxBins; ++k) {
          fEtaPhiBin[iPt][k]->Fill(eta, phi);
        }
      }
    }
  }
  for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
    if (counterBin[iPt] < 1)
      zeromult = kTRUE;
    fMultBin[iPt]->Fill(counterBin[iPt]);
  }
}

/*________________________________________________________________________
          Resetting Histograms after each event
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMomentsPID::ResetHistograms()
{
  for (Int_t i = 0; i < mPtBins; ++i) {
    for (Int_t j = 0; j < maxBins; ++j) {
      if (fEtaPhiBin[i][j])
        fEtaPhiBin[i][j]->Reset();
      if (flagMC) {
        if (fEtaPhiBinGen[i][j])
          fEtaPhiBinGen[i][j]->Reset();
      }
    }
  }
}

/*________________________________________________________________________
          Fill Track Information for NFM Calculation (Generated Level MC)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMomentsPID::FillMCTrackInfo()
{
  Int_t counterBin[maxPtBins] = { 0 };
  for (Int_t i_MCtrk = 0; i_MCtrk < fMCEvent->GetNumberOfTracks(); i_MCtrk++) {
    AliVParticle* lPart = (AliAODMCParticle*)fMCEvent->GetTrack(i_MCtrk);
    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(
      fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL)
      return;

    Bool_t isoobPileup = kFALSE;
    if (flagPileup)
      isoobPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(
        i_MCtrk, mcHeader, AODMCTrackArray);

    Float_t lpt = lPart->Pt();
    Float_t leta = lPart->Eta();
    Float_t lphi = lPart->Phi();
    Int_t lcharge = lPart->Charge();

    if ((TMath::Abs(lPart->PdgCode()) != 2212))
      continue;

    if (!lPart || !lPart->IsPhysicalPrimary() || isoobPileup || lpt < 0.2 ||
        (leta < minEta) || (leta > maxEta) || lcharge == 0)
      continue;
    GetPtBin(lpt);

    fHistQAEta[1]->Fill(leta);
    fHistQAPhi[1]->Fill(lphi);

    for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
      if (ptbin[iPt]) {
        counterBin[iPt]++;
        fPhiBinGen[iPt]->Fill(lphi);
        fEtaBinGen[iPt]->Fill(leta);
        fPtBinGen[iPt]->Fill(lpt);
        for (Int_t k = 0; k < maxBins; ++k) {
          fEtaPhiBinGen[iPt][k]->Fill(leta, lphi);
        }
      }
    }
  }

  for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
    fMultBinGen[iPt]->Fill(counterBin[iPt]);
  }
}

/*________________________________________________________________________
          Return the Pt Bin for the given Pt
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMomentsPID::GetPtBin(Double_t pt)
{
  ptbin.clear();
  for (Int_t i = 0; i < ptarray.GetSize() - 1; i += 2) {
    if (pt >= ptarray[i] && pt <= ptarray[i + 1]) {
      ptbin.push_back(kTRUE);
    } else {
      ptbin.push_back(kFALSE);
    }
  }
}

/*________________________________________________________________________
          Calculation of the Normalized Fq MomentsPID
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMomentsPID::CalculateNFMs(TH2D* h1[maxPtBins][maxBins], Bool_t mcGen)
{
  for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
    for (Int_t iM = 0; iM < maxBins; iM++) {

      Double_t NoOfBins = 0;
      Double_t MSquare = 0;
      if (Mmax == 123)
        NoOfBins = 3 * (iM + 2);
      if (Mmax == 82)
        NoOfBins = 2 * (iM + 2);
      if (flagSelfAff)
        NoOfBins = mMBin2[iM] * mNBin2[iM];
      MSquare = TMath::Power(NoOfBins, maxDim);
      if (flagSelfAff)
        MSquare = mMBin2[iM] * mNBin2[iM];
      Double_t SumOfbincontent = 0;
      Double_t FqEvent[maxqs];
      Double_t sumoff[maxqs];
      Double_t WCsumoff[maxqs];
      Double_t WCFqEvent[maxqs];
      Double_t bincontent = 0.0;
      Double_t WCSumOfbincontent = 0;
      Double_t Mbin = NoOfBins;
      Int_t NofXetabins = 0;
      Int_t NofXphibins = 0;

      for (Int_t index = 0; index < maxqs; index++) {
        FqEvent[index] = 0.0;
        sumoff[index] = 0.0;
        WCFqEvent[index] = 0.0;
        WCsumoff[index] = 0.0;
      }

      NofXetabins = h1[iPt][iM]->GetNbinsX();
      NofXphibins = h1[iPt][iM]->GetNbinsY();

      for (Int_t etabin = 1; etabin <= NofXetabins; etabin++) {
        for (Int_t phibin = 1; phibin <= NofXphibins; phibin++) {
          bincontent = 0.0;
          bincontent = h1[iPt][iM]->GetBinContent(etabin, phibin);

          Float_t CorFactor = 0.0;

          if (iPt == 0)
            CorFactor = (fhMapEtaPhiBin1M[iM]->GetBinContent(etabin, phibin)); //=====Correction factor value for Bin1=====

          if (iPt == 1)
            CorFactor = (fhMapEtaPhiBin2M[iM]->GetBinContent(etabin, phibin)); //====>=====

          if (iPt == 2)
            CorFactor = (fhMapEtaPhiBin3M[iM]->GetBinContent(etabin, phibin)); //====>=====

          if (iPt == 3)
            CorFactor = (fhMapEtaPhiBin4M[iM]->GetBinContent(etabin, phibin)); //====>=====

          SumOfbincontent += bincontent;

          if (TMath::IsNaN(bincontent / CorFactor)) {
            WCSumOfbincontent += 0;
          } else {
            WCSumOfbincontent += (bincontent / CorFactor);
          }

          for (Int_t q = 0; q < maxqs; q++) {
            if (bincontent >= (q + 2)) {

              Double_t Fqeofbin = 0.0;
              Double_t WCFqeofbin = 0.0;

              Fqeofbin = TMath::Factorial(bincontent) /
                         TMath::Factorial(bincontent - (q + 2));

              WCFqeofbin = Fqeofbin / (TMath::Power(CorFactor, (q + 2))); // =====>Efficiency corrected Fq Moment

              if (flagMC) {
                if (TMath::IsNaN(Fqeofbin)) {
                  break;
                }

                if (TMath::IsNaN(WCFqeofbin)) {
                  break;
                }
              }
              sumoff[q] += Fqeofbin;
              WCsumoff[q] += WCFqeofbin;
            }
          }
        }
      }

      Double_t Av_bincontent;
      Double_t WCAv_bincontent;

      Av_bincontent = SumOfbincontent / MSquare;
      WCAv_bincontent = WCSumOfbincontent / MSquare;

      for (Int_t q = 0; q < maxqs; q++) {
        if (sumoff[q] > 0) {
          FqEvent[q] = sumoff[q] / (MSquare);
        }

        if (WCsumoff[q] > 0) {
          WCFqEvent[q] = WCsumoff[q] / (MSquare); //..EFFICIENCY CORRECTED....To be saved in NTu
        }
      }

      Float_t Fq2e = FqEvent[0];
      Float_t Fq3e = FqEvent[1];
      Float_t Fq4e = FqEvent[2];
      Float_t Fq5e = FqEvent[3];
      Float_t Fq6e = FqEvent[4];
      Float_t Fq7e = FqEvent[5];

      Float_t WCFq2e = WCFqEvent[0]; // Bin average of fqe
      Float_t WCFq3e = WCFqEvent[1];
      Float_t WCFq4e = WCFqEvent[2];
      Float_t WCFq5e = WCFqEvent[3];

      if (!mcGen) {
        fntpMBin[iPt][iM]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e,
                                Fq6e, Fq7e);

        fntpMBinCorr[iPt][iM]->Fill(Mbin, WCAv_bincontent, WCFq2e, WCFq3e, WCFq4e, WCFq5e);

      } else {
        fntpMBinGen[iPt][iM]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e,
                                   Fq5e, Fq6e, Fq7e);
      }
    } // end of MBin loop
  }   // end of calculation loop
}

/*________________________________________________________________________
          Post data to the output slot(s)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMomentsPID::DataPosting()
{
  PostData(1, fHistList);
  PostData(2, fNtupleListBin1);
  PostData(3, fNtupleListBin2);
  PostData(4, fNtupleListBin3);
  PostData(5, fNtupleListBin4);
  for (Int_t i = 0; i < mPtBins; i++) {
    PostData(i + 2, fNtupleList[i]);
    if (flagMC)
      PostData(i + 2 + mPtBins, fNtupleListGen[i]);
  }
}

/*________________________________________________________________________
                                Terminate
________________________________________________________________________*/
void AliAnalysisTaskNFactorialMomentsPID::Terminate(Option_t*)
{
  Info("AliAnalysisTaskNFactorialMomentsPID", "Task Successfully finished");
  // terminate
  // called at the END of the analysis (when all events are processed)
}
