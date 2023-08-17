// Implementation of the Intermittency analysis
// for charged particles, two dimensions (eta and phi)
// Data and MC
// Contact: ramni.gupta@cern.ch
// Contributors: R.Gupta, S.Sharma, S.K.Malik

#include "AliAnalysisTaskNFactorialMoments.h"

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

class AliAnalysisTaskNFactorialMoments;
using namespace std;

ClassImp(AliAnalysisTaskNFactorialMoments)
  AliAnalysisTaskNFactorialMoments::AliAnalysisTaskNFactorialMoments()
  : AliAnalysisTaskSE(), fAOD(0), fMCEvent(0), fOutHList(0), fQAList(0), fQAList2(0), fEventCuts(0), fPIDResponse(0), fPIDCombined(0), mcHeader{ nullptr }, fNtupleListBin1(0), fNtupleListBin2(0), fNtupleListBin3(0), fNtupleListBin4(0), fNtupleListBin{ nullptr }, fNtupleListBin1Gen(0), fNtupleListBin2Gen(0), fNtupleListBin3Gen(0), fNtupleListBin4Gen(0), fHistQAEta{ nullptr }, fHistQAPhi{ nullptr }, fHistQAVx(0), fHistdEta{ nullptr }, fHistdPhi{ nullptr }, fHistQAVy(0), counter(0), fHistQAVz(0), fHistPtBin{ nullptr }, fEtaBin{ nullptr }, fPhiBin{ nullptr }, fHistMulBin{ nullptr }, fHistPtBinGen{ nullptr }, fEtaBinGen{ nullptr }, fPhiBinGen{ nullptr }, fHistMulBinGen{ nullptr }, fHistbeforeHBT{ nullptr }, fHistafterHBT{ nullptr }, fHistQAPID{ nullptr }, fHistPDG{ nullptr }, fHistDCAxy{ nullptr }, fHistDCAz{ nullptr }, fHistnITScls{ nullptr }, fHistnTPCcls{ nullptr }, fHistnTPCcrossedrows{ nullptr }, fHistnchi2ITScls{ nullptr }, fHistnchi2TPCcls{ nullptr }, fHistDCAxypT{ nullptr }, fHistDCAzpT{ nullptr }, fHistnsharedcls{ nullptr }, fHistnshclsfra{ nullptr }, fHistnshclsfravspt{ nullptr }, fHistnfoundcls{ nullptr }, fHistnfcls{ nullptr }, fEventCounter(0), fTrackCounter(0), fHistQACent(0), ptbin1(0), ptbin2(0), ptbin3(0), ptbin4(0), mfield(0) {}
AliAnalysisTaskNFactorialMoments::AliAnalysisTaskNFactorialMoments(
  const char* name)
  : AliAnalysisTaskSE(name), fAOD(0), fMCEvent(0), fOutHList(0), fQAList(0), fQAList2(0), fEventCuts(0), fPIDResponse(0), fPIDCombined(0), mcHeader{ nullptr }, fNtupleListBin1(0), fNtupleListBin2(0), fNtupleListBin3(0), fNtupleListBin4(0), fNtupleListBin{ nullptr }, fNtupleListBin1Gen(0), fNtupleListBin2Gen(0), fNtupleListBin3Gen(0), fNtupleListBin4Gen(0), fHistQAEta{ nullptr }, fHistQAPhi{ nullptr }, fHistQAVx(0), fHistdEta{ nullptr }, fHistdPhi{ nullptr }, fHistQAVy(0), counter(0), fHistQAVz(0), fHistPtBin{ nullptr }, fEtaBin{ nullptr }, fPhiBin{ nullptr }, fHistMulBin{ nullptr }, fHistPtBinGen{ nullptr }, fEtaBinGen{ nullptr }, fPhiBinGen{ nullptr }, fHistMulBinGen{ nullptr }, fHistbeforeHBT{ nullptr }, fHistafterHBT{ nullptr }, fHistQAPID{ nullptr }, fHistPDG{ nullptr }, fHistDCAxy{ nullptr }, fHistDCAz{ nullptr }, fHistnITScls{ nullptr }, fHistnTPCcls{ nullptr }, fHistnTPCcrossedrows{ nullptr }, fHistnchi2ITScls{ nullptr }, fHistnchi2TPCcls{ nullptr }, fHistDCAxypT{ nullptr }, fHistDCAzpT{ nullptr }, fHistnsharedcls{ nullptr }, fHistnshclsfra{ nullptr }, fHistnshclsfravspt{ nullptr }, fHistnfoundcls{ nullptr }, fHistnfcls{ nullptr }, fEventCounter(0), fTrackCounter(0), fHistQACent(0), ptbin1(0), ptbin2(0), ptbin3(0), ptbin4(0), mfield(0)
{
  // Constructor
  Info("AliAnalysisTaskNFactorialMoments", "Specific Constructor");
  for (int j = 0; j < M; j++) {
    fHEtaPhiBin1[j] = nullptr;
    fHEtaPhiBin2[j] = nullptr;
    fHEtaPhiBin3[j] = nullptr;
    fHEtaPhiBin4[j] = nullptr;
    fntpMBin1[j] = nullptr;
    fntpMBin2[j] = nullptr;
    fntpMBin3[j] = nullptr;
    fntpMBin4[j] = nullptr;

    if (fismc) {
      fHEtaPhiBin1Gen[j] = nullptr;
      fHEtaPhiBin2Gen[j] = nullptr;
      fHEtaPhiBin3Gen[j] = nullptr;
      fHEtaPhiBin4Gen[j] = nullptr;
      fntpMBin1Gen[j] = nullptr;
      fntpMBin2Gen[j] = nullptr;
      fntpMBin3Gen[j] = nullptr;
      fntpMBin4Gen[j] = nullptr;
    }
  }

  for (int i = 0; i < AliPID::kSPECIES; i++) {
    nsigmaTPC[i] = -999;
    nsigmaTOF[i] = -999;
  }

  DefineInput(0, TChain::Class());

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  if (fismc) {
    DefineOutput(8, TList::Class());
    DefineOutput(9, TList::Class());
    DefineOutput(10, TList::Class());
    DefineOutput(11, TList::Class());
  }
}

AliAnalysisTaskNFactorialMoments::~AliAnalysisTaskNFactorialMoments()
{
  // destructor
  if (fOutHList) {
    delete fOutHList;
  }
  if (fQAList2) {
    delete fQAList2;
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

void AliAnalysisTaskNFactorialMoments::UserCreateOutputObjects()
{
  fOutHList = new TList();
  fQAList2 = new TList();
  fNtupleListBin1 = new TList();
  fNtupleListBin2 = new TList();
  fNtupleListBin3 = new TList();
  fNtupleListBin4 = new TList();

  fOutHList->SetOwner(kTRUE);
  fQAList2->SetOwner(kTRUE);
  fNtupleListBin1->SetOwner(kTRUE);
  fNtupleListBin2->SetOwner(kTRUE);
  fNtupleListBin3->SetOwner(kTRUE);
  fNtupleListBin4->SetOwner(kTRUE);

  if (fismc) {
    fNtupleListBin1Gen = new TList();
    fNtupleListBin2Gen = new TList();
    fNtupleListBin3Gen = new TList();
    fNtupleListBin4Gen = new TList();
    fNtupleListBin1Gen->SetOwner(kTRUE);
    fNtupleListBin2Gen->SetOwner(kTRUE);
    fNtupleListBin3Gen->SetOwner(kTRUE);
    fNtupleListBin4Gen->SetOwner(kTRUE);
  } else {
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAList, kTRUE);
  }

  // QA Histograms for deta-dphi and PID
  fHistdEta = new TH1D("deta", "deta", 400, -0.2, 0.2);
  fQAList2->Add(fHistdEta);
  fHistdPhi = new TH1D("dphi", "dphi", 400, -0.2, 0.2);
  fQAList2->Add(fHistdPhi);
  TString name, name_;
  for (int i = 1; i < 5; i++) {
    name = Form("fdetadphibefore_ptbin%d", i);
    fHistbeforeHBT[i - 1] =
      new TH2D(name, name, 50, -1.6, 1.6, 50, -TMath::Pi(), TMath::Pi());
    fQAList2->Add(fHistbeforeHBT[i - 1]);
  }
  for (int i = 1; i < 5; i++) {
    name_ = Form("fdetadphiafter_ptbin%d", i);
    fHistafterHBT[i - 1] =
      new TH2D(name_, name_, 50, -1.6, 1.6, 50, -TMath::Pi(), TMath::Pi());
    fQAList2->Add(fHistafterHBT[i - 1]);
  }

  if (!fismc) {
    name = "TPCdEdx";
    name_ = "TPCdEdx vs momentum";
    fHistQAPID[0] = new TH2D(name, name_, 100, 0.0, 5.0, 250, 0.0, 500.0);
    fQAList2->Add(fHistQAPID[0]);

    name = "TPCNsigma";
    name_ = "TPCNsigma vs momentum";
    fHistQAPID[12] = new TH2D(name, name_, 100, 0.0, 5.0, 100, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[1]);

    name = "TPCNsigmaProton";
    name_ = "TPCNsigmaProton vs momentum";
    fHistQAPID[1] = new TH2D(name, name_, 100, 0.0, 2.0, 100, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[1]);

    name = "TOFNsigmaProton";
    name_ = "TOFNsigmaProton vs momentum";
    fHistQAPID[2] = new TH2D(name, name_, 100, 0.0, 5.0, 250, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[2]);

    name = "TPCTOFProton";
    name_ = "TPCTOFProton vs momentum";
    fHistQAPID[3] = new TH2D(name, name_, 100, 0.0, 5.0, 250, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[3]);

    name = "TOF signal";
    name_ = "momentum vs beta";
    fHistQAPID[4] = new TH2D(name, name_, 1000, 0.5, 5.0, 1000, 0.4, 1.2);
    fQAList2->Add(fHistQAPID[4]);

    name = "nSigmaTPC  vs nSigmaTOF";
    name_ = "nSigmaTPC  vs nSigmaTOF";
    fHistQAPID[5] = new TH2D(name, name_, 1000, -5.0, 5.0, 1000, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[5]);

    name = "TPCdEdxAfter";
    name_ = "TPCdEdx vs momentum";
    fHistQAPID[6] = new TH2D(name, name_, 100, 0.0, 10.0, 250, 0.0, 250.0);
    fQAList2->Add(fHistQAPID[6]);

    name = "TPCNsigmaProtonAfter";
    name_ = "TPCNsigmaProton vs momentum";
    fHistQAPID[7] = new TH2D(name, name_, 100, 0.0, 2.0, 100, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[7]);

    name = "TOFNsigmaProtonAfter";
    name_ = "TOFNsigmaProton vs momentum";
    fHistQAPID[8] = new TH2D(name, name_, 100, 0.0, 5.0, 250, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[8]);

    name = "TPCTOFProtonAfter";
    name_ = "TPCTOFProton vs momentum";
    fHistQAPID[9] = new TH2D(name, name_, 100, 0.0, 5.0, 250, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[9]);

    name = "TOF signalAfter";
    name_ = "momentum vs beta";
    fHistQAPID[10] = new TH2D(name, name_, 1000, 0.5, 5.0, 1000, 0.4, 1.2);
    fQAList2->Add(fHistQAPID[10]);

    name = "nSigmaTPC  vs nSigmaTOFAfter";
    name_ = "nSigmaTPC  vs nSigmaTOF";
    fHistQAPID[11] = new TH2D(name, name_, 1000, -5.0, 5.0, 1000, -5.0, 5.0);
    fQAList2->Add(fHistQAPID[11]);

  } else {
    name = "PDGHist";
    name_ = "PDG code of the generated particles";
    fHistPDG[0] = new TH1D(name, name_, 1000, -500, 500);
    fQAList2->Add(fHistPDG[0]);
    name = "ElPosPDGHist";
    name_ = "PDG code of the generated electrons and positrons";
    fHistPDG[1] = new TH1D(name, name_, 100, -100, 100);
    fQAList2->Add(fHistPDG[1]);
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
  fOutHList->Add(fEventCounter);

  // Trackcounter hist defined
  fTrackCounter =
    new TH1D("fTrackCounter", "histo to keep track of tracks", 19, 0.5, 19.5);
  fTrackCounter->GetXaxis()->SetBinLabel(1, "Total tracks");
  fTrackCounter->GetXaxis()->SetBinLabel(2, "FilterBit");
  fTrackCounter->GetXaxis()->SetBinLabel(3, "Rej PIDs");
  fTrackCounter->GetXaxis()->SetBinLabel(4, "Loose cuts");
  fTrackCounter->GetXaxis()->SetBinLabel(5, "Systematic cuts");
  fTrackCounter->GetXaxis()->SetBinLabel(6, "DCAs");
  fTrackCounter->GetXaxis()->SetBinLabel(7, "nSharedCls/nCls");
  fTrackCounter->GetXaxis()->SetBinLabel(8, "nSharedCls/nCrRows");
  fTrackCounter->GetXaxis()->SetBinLabel(9, "nCrRows/nFindableCls");
  fTrackCounter->GetXaxis()->SetBinLabel(10, "TwoTracks");
  fTrackCounter->GetXaxis()->SetBinLabel(11, "ptbin1");
  fTrackCounter->GetXaxis()->SetBinLabel(12, "ptbin2");
  fTrackCounter->GetXaxis()->SetBinLabel(13, "ptbin3");
  fTrackCounter->GetXaxis()->SetBinLabel(14, "ptbin4");
  fOutHList->Add(fTrackCounter);

  // Hists defined
  fHistQACent = new TH1F("fHistQACent", "Centrality Distribution", 10,
                         minCent - 2, maxCent + 2);
  fOutHList->Add(fHistQACent);

  fHistQAVx =
    new TH1F("fHistQAVx", "Primary vertex distribution -x coordinate;V_{x}",
             1500, -0.3, 0.3);
  fHistQAVy =
    new TH1F("fHistQAVy", "Primary vertex distribution -y coordinate;V_{y}",
             1500, -0.3, 0.3);
  fHistQAVz =
    new TH1F("fHistQAVz", "Primary vertex distribution -z coordinate;V_{z}",
             1500, -20.0, 20.0);
  fOutHList->Add(fHistQAVx);
  fOutHList->Add(fHistQAVy);
  fOutHList->Add(fHistQAVz);

  name = "fHistQAEta";
  if (fismc)
    name += "Re";
  name_ = "#eta distribution of the track";
  fHistQAEta[0] = new TH1F(name, name_, 1000, minEta, maxEta);
  fOutHList->Add(fHistQAEta[0]);
  if (fismc) {
    name = "fHistQAEtaGen";
    fHistQAEta[1] = new TH1F(name, name_, 1000, minEta, maxEta);
    fOutHList->Add(fHistQAEta[1]);
  }
  name = "fHistQAPhi";
  if (fismc)
    name += "Re";
  name_ = "#phi distribution of the track";
  fHistQAPhi[0] = new TH1F(name, name_, 1000, 0.0, 6.5);
  fOutHList->Add(fHistQAPhi[0]);
  if (fismc) {
    name = "fHistQAPhiGen";
    fHistQAPhi[1] = new TH1F(name, name_, 1000, 0.0, 6.5);
    fOutHList->Add(fHistQAPhi[1]);
  }

  for (int nOfptbins(1); nOfptbins < 5; nOfptbins++) {
    if (fismc) {
      name = Form("fHPtBinRe%d", nOfptbins);
      name_ = Form("fHPtBinGen%d", nOfptbins);
      fHistPtBinGen[nOfptbins - 1] =
        new TH1F(name_, "p_{T}dist", 1000, 0.1, 5.5);
      fOutHList->Add(fHistPtBinGen[nOfptbins - 1]);
    } else
      name = Form("fHPtBin%d", nOfptbins);
    fHistPtBin[nOfptbins - 1] = new TH1F(name, "p_{T}dist", 1000, 0.1, 5.5);
    if (fismc) {
      name = Form("fEtaBinRe%d", nOfptbins);
      name_ = Form("fEtaBinGen%d", nOfptbins);
      fEtaBinGen[nOfptbins - 1] =
        new TH1F(name_, "#eta dist", 1000, minEta, maxEta);
      fOutHList->Add(fEtaBinGen[nOfptbins - 1]);
    } else
      name = Form("fEtaBin%d", nOfptbins);
    fEtaBin[nOfptbins - 1] = new TH1F(name, "#eta dist", 1000, minEta, maxEta);
    if (fismc) {
      name = Form("fPhiBinRe%d", nOfptbins);
      name_ = Form("fPhiBinGen%d", nOfptbins);
      fPhiBinGen[nOfptbins - 1] = new TH1F(name_, "#phi dist", 1000, 0.0, 6.5);
      fOutHList->Add(fPhiBinGen[nOfptbins - 1]);
    } else
      name = Form("fPhiBin%d", nOfptbins);
    fPhiBin[nOfptbins - 1] = new TH1F(name, "#phi dist", 1000, 0., 6.5);
    if (fismc) {
      name = Form("fHMultBinRe%d", nOfptbins);
      name_ = Form("fHMultBinGen%d", nOfptbins);
      fHistMulBinGen[nOfptbins - 1] =
        new TH1F(name_, "Mult dist", 2000, 0.0, 9200.0);
      fOutHList->Add(fHistMulBinGen[nOfptbins - 1]);
    } else
      name = Form("fHMultBin%d", nOfptbins);
    fHistMulBin[nOfptbins - 1] = new TH1F(name, "Mult dist", 2000, 0.0, 9200.0);
    fOutHList->Add(fHistPtBin[nOfptbins - 1]);
    fOutHList->Add(fHistMulBin[nOfptbins - 1]);
    fOutHList->Add(fEtaBin[nOfptbins - 1]);
    fOutHList->Add(fPhiBin[nOfptbins - 1]);
  }

  fHistDCAxy = new TH1F("fHistDCAxy", "DCAxy distribution", 900, -4.5, 4.5);
  fHistDCAz = new TH1F("fHistDCAz", "DCAz distribution", 900, -4.5, 4.5);
  fOutHList->Add(fHistDCAxy);
  fOutHList->Add(fHistDCAz);
  fHistDCAxypT = new TH2F("fHistDCAxypT", "DCAxy vs pT;p_{T};dca_{xy}", 50, 0, 5.5, 40,
                          -0.4, 0.4);
  fHistDCAzpT = new TH2F("fHistDCAzpT", "DCAz vs pT;p_{T};dca_{z}", 50, 0, 5.5, 50, -0.5,
                         0.5);
  fOutHList->Add(fHistDCAxypT);
  fOutHList->Add(fHistDCAzpT);
  fHistnITScls =
    new TH1F("fHistnITScls", "ITS cluster distribution", 10, -0.5, 10.5);
  fOutHList->Add(fHistnITScls);
  fHistnTPCcls =
    new TH1F("fHistnTPCcls", "TPC cluster distribution", 200, 0, 200);
  fOutHList->Add(fHistnTPCcls);
  fHistnchi2ITScls = new TH1F("fHistnchi2ITScls",
                              "ITSchi2 cluster distribution", 100, -0.5, 100.5);
  fOutHList->Add(fHistnchi2ITScls);
  fHistnchi2TPCcls =
    new TH1F("fHistnchi2TPCcls", "TPCchi2 cluster distribution", 50, 0, 50);
  fOutHList->Add(fHistnchi2TPCcls);
  fHistnTPCcrossedrows = new TH1F("fHistnTPCcrossedrows",
                                  "TPC crossed rows distribution", 200, 0, 200);
  fOutHList->Add(fHistnTPCcrossedrows);

  TString name1;
  for (int i = 0; i < 2; ++i) {
    if (i == 0)
      name1 = "before";
    else
      name1 = "after";
    fHistnsharedcls[i] = new TH1F(Form("fHistnsharedcls%s", name1.Data()),
                                  Form("TPC shared cluster distribution %s", name1.Data()), 200, 0, 200);
    fOutHList->Add(fHistnsharedcls[i]);
    fHistnshclsfra[i] = new TH2F(Form("fHistnshclsfra%s", name1.Data()),
                                 Form("TPC shared cluster fraction %s;sharedcls/ncls;sharedcls/ncrrows", name1.Data()), 100, 0,
                                 1.0, 100, 0, 1.0);
    fOutHList->Add(fHistnshclsfra[i]);
    fHistnfoundcls[i] = new TH2F(Form("fHistnfoundcls%s", name1.Data()),
                                 Form("TPC found cluster fraction %s;sharedcls/ncls;ncrrows/findablecls", name1.Data()), 200, 0, 2,
                                 200, 0, 2);
    fOutHList->Add(fHistnfoundcls[i]);
    fHistnfcls[i] = new TH1F(Form("fHistnfcls%s", name1.Data()), Form("TPC findable cluster distribution %s;ncrrows/findablecls;counts", name1.Data()), 40, 0, 4);
    fOutHList->Add(fHistnfcls[i]);
  }

  for (int i = 0; i < 5; ++i) {
    fHistnshclsfravspt[i] = new TH2F(Form("fHistnshclsfravspt%d", i), Form("TPC nshared clusters vs #it{p_{T}} %i", i), 500, 0, 5, 400, 0.02, 1);
    fOutHList->Add(fHistnshclsfravspt[i]);
  }
  for (int nOfHist(1); nOfHist < M + 1; nOfHist++) {
    if (fismc) {
      name = Form("fHistEtaPhiBin1Re_%d", (nOfHist));
      name_ = Form("fHistEtaPhiBin1Gen_%d", (nOfHist));
      fHEtaPhiBin1Gen[nOfHist - 1] = new TH2D(
        name_, Form("PtCut1 (#eta) and (#phi) M%d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
      if (fSelfAffine)
        fHEtaPhiBin1Gen[nOfHist - 1] = new TH2D(
          name_, Form("PtCut1 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
          Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
      fOutHList->Add(fHEtaPhiBin1Gen[nOfHist - 1]);
    } else
      name = Form("fHistEtaPhiBin1_%d", (nOfHist));
    fHEtaPhiBin1[nOfHist - 1] = new TH2D(
      name, Form("PtCut1 (#eta) and (#phi) M%d Binning", (nOfHist)),
      Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
    if (fSelfAffine)
      fHEtaPhiBin1[nOfHist - 1] = new TH2D(
        name, Form("PtCut1 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
    if (fismc) {
      name = Form("fHistEtaPhiBin2Re_%d", (nOfHist));
      name_ = Form("fHistEtaPhiBin2Gen_%d", (nOfHist));
      fHEtaPhiBin2Gen[nOfHist - 1] = new TH2D(
        name_, Form("PtCut2 (#eta) and (#phi) M%d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
      if (fSelfAffine)
        fHEtaPhiBin2Gen[nOfHist - 1] = new TH2D(
          name_, Form("PtCut2 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
          Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
      fOutHList->Add(fHEtaPhiBin2Gen[nOfHist - 1]);
    } else
      name = Form("fHistEtaPhiBin2_%d", (nOfHist));
    fHEtaPhiBin2[nOfHist - 1] = new TH2D(
      name, Form("PtCut2 (#eta) and (#phi) M%d Binning", (nOfHist)),
      Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
    if (fSelfAffine)
      fHEtaPhiBin2[nOfHist - 1] = new TH2D(
        name, Form("PtCut2 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
    if (fismc) {
      name = Form("fHistEtaPhiBin3Re_%d", (nOfHist));
      name_ = Form("fHistEtaPhiBin3Gen_%d", (nOfHist));
      fHEtaPhiBin3Gen[nOfHist - 1] = new TH2D(
        name_, Form("PtCut3 (#eta) and (#phi) M%d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
      if (fSelfAffine)
        fHEtaPhiBin3Gen[nOfHist - 1] = new TH2D(
          name_, Form("PtCut3 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
          Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
      fOutHList->Add(fHEtaPhiBin3Gen[nOfHist - 1]);
    } else
      name = Form("fHistEtaPhiBin3_%d", (nOfHist));
    fHEtaPhiBin3[nOfHist - 1] = new TH2D(
      name, Form("PtCut3 (#eta) and (#phi) M%d Binning", (nOfHist)),
      Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
    if (fSelfAffine)
      fHEtaPhiBin3[nOfHist - 1] = new TH2D(
        name, Form("PtCut3 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
    if (fismc) {
      name = Form("fHistEtaPhiBin4Re_%d", (nOfHist));
      name_ = Form("fHistEtaPhiBin4Gen_%d", (nOfHist));
      fHEtaPhiBin4Gen[nOfHist - 1] = new TH2D(
        name_, Form("PtCut4 (#eta) and (#phi) M%d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
      if (fSelfAffine)
        fHEtaPhiBin4Gen[nOfHist - 1] = new TH2D(
          name_, Form("PtCut4 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
          Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
      fOutHList->Add(fHEtaPhiBin4Gen[nOfHist - 1]);
    } else
      name = Form("fHistEtaPhiBin4_%d", (nOfHist));
    fHEtaPhiBin4[nOfHist - 1] = new TH2D(
      name, Form("PtCut4 (#eta) and (#phi) M%d Binning", (nOfHist)),
      Mbins[nOfHist - 1], minEta, maxEta, Mbins[nOfHist - 1], 0.0, 6.30);
    if (fSelfAffine)
      fHEtaPhiBin4[nOfHist - 1] = new TH2D(
        name, Form("PtCut4 (#eta) and (#phi) MxN %d Binning", (nOfHist)),
        Mbins[nOfHist - 1], minEta, maxEta, Nbins[nOfHist - 1], 0.0, 6.30);
    fOutHList->Add(fHEtaPhiBin1[nOfHist - 1]);
    fOutHList->Add(fHEtaPhiBin2[nOfHist - 1]);
    fOutHList->Add(fHEtaPhiBin3[nOfHist - 1]);
    fOutHList->Add(fHEtaPhiBin4[nOfHist - 1]);
  }
  for (int nOfTuple = 1; nOfTuple < M + 1; nOfTuple++) {
    name_ = Form("Fmsforetaphicut_%d", (nOfTuple));
    if (fismc) {
      name = Form("fntpMBin1Gen_%d", (nOfTuple));
      fntpMBin1Gen[nOfTuple - 1] = new TNtuple(
        name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
      fNtupleListBin1Gen->Add(fntpMBin1Gen[nOfTuple - 1]);
      name = Form("fntpMBin1Re_%d", (nOfTuple));
    } else
      name = Form("fntpMBin1%d", (nOfTuple));
    fntpMBin1[nOfTuple - 1] = new TNtuple(
      name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    if (fismc) {
      name = Form("fntpMBin2Gen_%d", (nOfTuple));
      fntpMBin2Gen[nOfTuple - 1] = new TNtuple(
        name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
      fNtupleListBin2Gen->Add(fntpMBin2Gen[nOfTuple - 1]);
      name = Form("fntpMBin2Re_%d", (nOfTuple));
    } else
      name = Form("fntpMBin2%d", (nOfTuple));
    fntpMBin2[nOfTuple - 1] = new TNtuple(
      name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    if (fismc) {
      name = Form("fntpMBin3Gen_%d", (nOfTuple));
      fntpMBin3Gen[nOfTuple - 1] = new TNtuple(
        name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
      fNtupleListBin3Gen->Add(fntpMBin3Gen[nOfTuple - 1]);
      name = Form("fntpMBin3Re_%d", (nOfTuple));
    } else
      name = Form("fntpMBin3%d", (nOfTuple));
    fntpMBin3[nOfTuple - 1] = new TNtuple(
      name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    if (fismc) {
      name = Form("fntpMBin4Gen_%d", (nOfTuple));
      fntpMBin4Gen[nOfTuple - 1] = new TNtuple(
        name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
      fNtupleListBin4Gen->Add(fntpMBin4Gen[nOfTuple - 1]);
      name = Form("fntpMBin4Re_%d", (nOfTuple));
    } else
      name = Form("fntpMBin4%d", (nOfTuple));
    fntpMBin4[nOfTuple - 1] = new TNtuple(
      name, name_, "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    fNtupleListBin1->Add(fntpMBin1[nOfTuple - 1]);
    fNtupleListBin2->Add(fntpMBin2[nOfTuple - 1]);
    fNtupleListBin3->Add(fntpMBin3[nOfTuple - 1]);
    fNtupleListBin4->Add(fntpMBin4[nOfTuple - 1]);
  }

  DataPosting();
}

// Execution - called per event:
void AliAnalysisTaskNFactorialMoments::UserExec(Option_t*)
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
    ::Fatal("AliAnalysisTaskNFactorialMoments::UserExec",
            "No AOD event found!");

  if (fismc) {
    fMCEvent = MCEvent();
    if (!fMCEvent)
      ::Fatal("AliAnalysisTaskNFactorialMoments::UserExec",
              "No MC event found!");
    mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(
      AliAODMCHeader::StdBranchName());
    if (!mcHeader)
      ::Fatal("AliAnalysisTaskNFactorialMoments::UserExec",
              "No MC header found!");
  } else {
    if (!fPIDResponse)
      ::Fatal("AliAnalysisTaskNFactorialMoments::UserExec",
              "No PIDResponse found!");
  }

  // fEventCuts.SetupLHC15o();
  if (fYear == "2015") {
    if (fpileup && !fismc) {
      fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
      bool isAcceptedEvent = fEventCuts.AcceptEvent(fAOD); // for QA
      if (!isAcceptedEvent)
        return;
    } else if (fpileup && fismc) {
      bool isPileupinGenEvent =
        AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader, fGenName);
      if (isPileupinGenEvent)
        return;
    }
  }

  fEventCounter->Fill(2);

  if (!fismc) {
    UInt_t fSelectMask = fInputHandler->IsEventSelected();
    if (fYear == "2015") {
      if (!(fSelectMask & AliVEvent::kINT7))
        return;
    } else if (fYear == "2010") {
      if (!(fSelectMask & AliVEvent::kMB))
        return;
    }
  }

  fEventCounter->Fill(3);

  AliAODVertex* vertex = fAOD->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1)
    return;
  float lvx = (float)vertex->GetX();
  float lvy = (float)vertex->GetY();
  float lvz = (float)vertex->GetZ();

  if ((fabs(lvx) < fVxMax) && (fabs(lvy) < fVyMax) && (fabs(lvz) < fVzMax)) {
    fHistQAVx->Fill(lvx);
    fHistQAVy->Fill(lvy);
    fHistQAVz->Fill(lvz);
  } else
    return;
  fEventCounter->Fill(4);

  // Centrality Selection
  double lMultiPercentile = -1;
  if (fYear == "2015") {
    // Centrality slecetion for LHC15o
    AliMultSelection* MultSelection =
      (AliMultSelection*)fAOD->FindListObject("MultSelection");
    lMultiPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  } else if (fYear == "2010") {
    AliCentrality* centrality = fAOD->GetCentrality();
    lMultiPercentile = centrality->GetCentralityPercentile("V0M");
  }
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

  // Filling Track Information
  FillTrackInfo();
  CalculateNFMs(fHEtaPhiBin1, fHEtaPhiBin2, fHEtaPhiBin3, fHEtaPhiBin4, kFALSE);

  if (fismc) {
    FillMCTrackInfo(); // Fill Generated Information
    CalculateNFMs(fHEtaPhiBin1Gen, fHEtaPhiBin2Gen, fHEtaPhiBin3Gen,
                  fHEtaPhiBin4Gen, kTRUE);
  }

  DataPosting();
}

/*________________________________________________________________________
      Fill Track Information for NFM Calculation (Main Track Loop)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::FillTrackInfo()
{
  int counterEtacutBin1 = 0, counterEtacutBin2 = 0, counterEtacutBin3 = 0,
      counterEtacutBin4 = 0;
  int nTracks(fAOD->GetNumberOfTracks());
  float dpstar, deta, dphi;
  float dcaXY, dcaZ;
  int trackbefore = 0, trackafter = 0, twotrackscount = 0;

  // Track Loop:
  for (int i(0); i < nTracks; i++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    fTrackCounter->Fill(1);
    if (!(track) || !(track->TestFilterBit(filterBit)))
      continue;
    counter++;
    fTrackCounter->Fill(2);

    float pt = track->Pt();
    float eta = track->Eta();
    float phi = track->Phi();
    int charge = track->Charge();
    int id = track->GetID();

    if (GetParticleID(track, kTRUE))
      continue;
    fTrackCounter->Fill(3);

    if ((fabs(eta) > 0.8) || (fabs(pt) < 0.2) || charge == 0)
      continue;
    // calculate if pt is in the pt bin range
    GetPtBin(pt);
    fTrackCounter->Fill(4);

    bool skiptracks = kFALSE;
    bool sharedtracks = kFALSE;
    float nSharedCls = track->GetTPCnclsS();
    float nCls = track->GetTPCNcls();
    float nCrossedRows = track->GetTPCCrossedRows();
    float nFindableCls = track->GetTPCNclsF();

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

    if ((fDCAxyMax > 0.0) && (fabs(dcaXY) > fDCAxyMax))
      continue;
    if ((fDCAzMax > 0.0) && (fabs(dcaZ) > fDCAzMax))
      continue;
    fHistDCAxy->Fill(dcaXY);
    fHistDCAz->Fill(dcaZ);
    fHistDCAxypT->Fill(pt, dcaXY);
    fHistDCAzpT->Fill(pt, dcaZ);
    fTrackCounter->Fill(6);
    fHistnsharedcls[0]->Fill(nSharedCls);
    fHistnfcls[0]->Fill(nCrossedRows / nFindableCls);
    fHistnshclsfra[0]->Fill(nSharedCls / nCls, nSharedCls / nCrossedRows);
    fHistnfoundcls[0]->Fill(nSharedCls / nCls, nCrossedRows / nFindableCls);

    if ((fSharedClsMax > 0) && ((nSharedCls / nCls) > fSharedClsMax))
      continue;
    fTrackCounter->Fill(7);

    if ((fSharedRowsMax > 0) && ((nSharedCls / nCrossedRows) > fSharedRowsMax))
      continue;
    fTrackCounter->Fill(8);

    if ((fFindableClsMin > 0) && ((nCrossedRows / nFindableCls) < fFindableClsMin))
      continue;
    fTrackCounter->Fill(9);

    fHistnsharedcls[1]->Fill(nSharedCls);
    fHistnfcls[1]->Fill(nCrossedRows / nFindableCls);
    fHistnshclsfra[1]->Fill(nSharedCls / nCls, nSharedCls / nCrossedRows);
    fHistnfoundcls[1]->Fill(nSharedCls / nCls, nCrossedRows / nFindableCls);

    if (fsharity || ftwotrack) {

      TBits clusmap = track->GetTPCClusterMap();
      TBits sharedmap = track->GetTPCSharedMap();

      for (int j(0); j < nTracks; j++) {
        AliAODTrack* track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
        if (!(track2->TestFilterBit(filterBit)))
          continue;

        int charge2 = track2->Charge();
        float pt2 = track2->Pt();
        float eta2 = track2->Eta();
        float phi2 = track2->Phi();
        int id2 = track2->GetID();

        if (GetParticleID(track2, kFALSE))
          continue;

        if (id == id2 || fabs(eta2) > 0.8 || fabs(pt2) < 0.2 || charge2 == 0 ||
            pt < pt2)
          continue;
        if (fsharity) {
          TBits clusmap2 = track2->GetTPCClusterMap();
          TBits sharedmap2 = track2->GetTPCSharedMap();
          float sharity =
            SharedClusterFraction(clusmap, clusmap2, sharedmap, sharedmap2);
          if (sharity > fSharedFraction) {
            sharedtracks = kTRUE;
            twotrackscount++;
          }
        }
        if (sharedtracks)
          break;
        if (ftwotrack) {
          dpstar = dphistarcalculation(phi, eta, pt, charge, phi2, eta2, pt2,
                                       charge2, mfield);
          if (dpstar == 999)
            continue;
          deta = eta2 - eta;
          if (fabs(deta) < fdeta && fabs(dpstar) < fdphi && charge == charge2)
            skiptracks = kTRUE;
        }
        if (skiptracks) {
          if (ptbin[0])
            fHistafterHBT[0]->Fill(deta, dpstar);
          if (ptbin[1])
            fHistafterHBT[1]->Fill(deta, dpstar);
          if (ptbin[2])
            fHistafterHBT[2]->Fill(deta, dpstar);
          if (ptbin[3])
            fHistafterHBT[3]->Fill(deta, dpstar);
          twotrackscount++;
          break;
        }
      }
    }
    if (skiptracks || sharedtracks) {
      continue;
    }
    trackafter++;
    fTrackCounter->Fill(10);

    fHistQAEta[0]->Fill(eta);
    fHistQAPhi[0]->Fill(phi);

    double profileVal[16] = { 0 };

    if (minCent == 0 && maxCent == 5) {
      double profileVal1[16] = { 0.360756, 0.360834, 0.360455, 0.35992, 0.359485, 0.359125, 0.359184, 0.359317, 0.359458, 0.35959, 0.360105, 0.36066, 0.361126, 0.361661, 0.361998, 0.362225 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 5 && maxCent == 10) {
      double profileVal1[16] = { 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379, 0.358205 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 10 && maxCent == 20) {
      double profileVal1[16] = { 0.35707, 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 20 && maxCent == 40) {
      double profileVal1[16] = { 0.35707, 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 40 && maxCent == 80) {
      double profileVal1[16] = { 0.35707, 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else {
      double profileVal1[16] = { 0.360756, 0.360834, 0.360455, 0.35992, 0.359485, 0.359125, 0.359184, 0.359317, 0.359458, 0.35959, 0.360105, 0.36066, 0.361126, 0.361661, 0.361998, 0.362225 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    }

    fHistnshclsfravspt[0]->Fill(pt, nSharedCls / nCrossedRows);

    if (pt > 0.4 && pt < 0.5) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[0]))
        continue;
    }
    if (pt > 0.5 && pt < 0.6) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[1]))
        continue;
    }
    if (pt > 0.6 && pt < 0.7) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[2]))
        continue;
    }
    if (pt > 0.7 && pt < 0.8) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[3]))
        continue;
    }
    if (pt > 0.8 && pt < 0.9) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[4]))
        continue;
    }
    if (pt > 0.9 && pt < 1.0) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[5]))
        continue;
    }
    if (pt > 1.0 && pt < 1.1) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[6]))
        continue;
    }
    if (pt > 1.1 && pt < 1.2) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[7]))
        continue;
    }
    if (pt > 1.2 && pt < 1.3) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[8]))
        continue;
    }
    if (pt > 1.3 && pt < 1.4) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[9]))
        continue;
    }
    if (pt > 1.4 && pt < 1.5) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[10]))
        continue;
    }
    if (pt > 1.5 && pt < 1.6) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[11]))
        continue;
    }
    if (pt > 1.6 && pt < 1.7) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[12]))
        continue;
    }
    if (pt > 1.7 && pt < 1.8) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[13]))
        continue;
    }
    if (pt > 1.8 && pt < 1.9) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[14]))
        continue;
    }
    if (pt > 1.9 && pt < 2.0) {
      if (fSharedFrMean && (nSharedCls / nCrossedRows > profileVal[15]))
        continue;
    }

    if (ptbin[0]) {
      fHistnshclsfravspt[1]->Fill(pt, nSharedCls / nCrossedRows);
      fTrackCounter->Fill(11);
      counterEtacutBin1++;
      fHistPtBin[0]->Fill(pt);
      fEtaBin[0]->Fill(eta);
      fPhiBin[0]->Fill(phi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin1[k]->Fill(eta, phi);
    }
    if (ptbin[1]) {
      fHistnshclsfravspt[2]->Fill(pt, nSharedCls / nCrossedRows);
      fTrackCounter->Fill(12);
      counterEtacutBin2++;
      fHistPtBin[1]->Fill(pt);
      fEtaBin[1]->Fill(eta);
      fPhiBin[1]->Fill(phi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin2[k]->Fill(eta, phi);
    }
    if (ptbin[2]) {
      fHistnshclsfravspt[3]->Fill(pt, nSharedCls / nCrossedRows);
      fTrackCounter->Fill(13);
      counterEtacutBin3++;
      fHistPtBin[2]->Fill(pt);
      fEtaBin[2]->Fill(eta);
      fPhiBin[2]->Fill(phi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin3[k]->Fill(eta, phi);
    }
    if (ptbin[3]) {
      fHistnshclsfravspt[4]->Fill(pt, nSharedCls / nCrossedRows);
      fTrackCounter->Fill(14);
      counterEtacutBin4++;
      fHistPtBin[3]->Fill(pt);
      fEtaBin[3]->Fill(eta);
      fPhiBin[3]->Fill(phi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin4[k]->Fill(eta, phi);
    }
  }

  fHistMulBin[0]->Fill(counterEtacutBin1);
  fHistMulBin[1]->Fill(counterEtacutBin2);
  fHistMulBin[2]->Fill(counterEtacutBin3);
  fHistMulBin[3]->Fill(counterEtacutBin4);
}

/*________________________________________________________________________
          Resetting Histograms after each event
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::ResetHistograms()
{
  for (Int_t i = 0; i < M; i++) {
    if (fHEtaPhiBin1[i])
      fHEtaPhiBin1[i]->Reset();
    if (fHEtaPhiBin2[i])
      fHEtaPhiBin2[i]->Reset();
    if (fHEtaPhiBin3[i])
      fHEtaPhiBin3[i]->Reset();
    if (fHEtaPhiBin4[i])
      fHEtaPhiBin4[i]->Reset();
    if (fismc) {
      if (fHEtaPhiBin1Gen[i])
        fHEtaPhiBin1Gen[i]->Reset();
      if (fHEtaPhiBin2Gen[i])
        fHEtaPhiBin2Gen[i]->Reset();
      if (fHEtaPhiBin3Gen[i])
        fHEtaPhiBin3Gen[i]->Reset();
      if (fHEtaPhiBin4Gen[i])
        fHEtaPhiBin4Gen[i]->Reset();
    }
  }
}

/*________________________________________________________________________
          Fill Track Information for NFM Calculation (Generated Level MC)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::FillMCTrackInfo()
{
  int counterEtacutBin1 = 0, counterEtacutBin2 = 0, counterEtacutBin3 = 0,
      counterEtacutBin4 = 0;
  for (int i_MCtrk = 0; i_MCtrk < fMCEvent->GetNumberOfTracks(); i_MCtrk++) {
    AliVParticle* lPart = (AliAODMCParticle*)fMCEvent->GetTrack(i_MCtrk);
    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(
      fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL)
      return;
    bool isoobPileup = kFALSE;
    if (fpileup)
      isoobPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(
        i_MCtrk, mcHeader, AODMCTrackArray);
    float lpt = lPart->Pt();
    float leta = lPart->Eta();
    float lphi = lPart->Phi();
    int lcharge = lPart->Charge();

    bool isProton = kFALSE;
    bool isElectron = kFALSE;
    if (lPart->PdgCode() == 11 || lPart->PdgCode() == -11)
      isElectron = kTRUE;
    if (lPart->PdgCode() == 2212 || lPart->PdgCode() == -2212)
      isProton = kTRUE;
    if (fRejectEls)
      if (isElectron)
        continue;
    if (fProtonAnalysis)
      if (isProton)
        continue;
    if (!lPart || !lPart->IsPhysicalPrimary() || isoobPileup || lpt < 0.2 ||
        fabs(leta) > 0.8 || lcharge == 0)
      continue;
    GetPtBin(lpt);

    fHistQAEta[1]->Fill(leta);
    fHistQAPhi[1]->Fill(lphi);

    if (ptbin[0]) {
      counterEtacutBin1++;
      fHistPtBinGen[0]->Fill(lpt);
      fEtaBinGen[0]->Fill(leta);
      fPhiBinGen[0]->Fill(lphi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin1Gen[k]->Fill(leta, lphi);
    }
    if (ptbin[1]) {
      counterEtacutBin2++;
      fHistPtBinGen[1]->Fill(lpt);
      fEtaBinGen[1]->Fill(leta);
      fPhiBinGen[1]->Fill(lphi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin2Gen[k]->Fill(leta, lphi);
    }
    if (ptbin[2]) {
      counterEtacutBin3++;
      fHistPtBinGen[2]->Fill(lpt);
      fEtaBinGen[2]->Fill(leta);
      fPhiBinGen[2]->Fill(lphi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin3Gen[k]->Fill(leta, lphi);
    }
    if (ptbin[3]) {
      counterEtacutBin4++;
      fHistPtBinGen[3]->Fill(lpt);
      fEtaBinGen[3]->Fill(leta);
      fPhiBinGen[3]->Fill(lphi);
      for (int k = 0; k < M; k++)
        fHEtaPhiBin4Gen[k]->Fill(leta, lphi);
    }
  }

  fHistMulBinGen[0]->Fill(counterEtacutBin1);
  fHistMulBinGen[1]->Fill(counterEtacutBin2);
  fHistMulBinGen[2]->Fill(counterEtacutBin3);
  fHistMulBinGen[3]->Fill(counterEtacutBin4);
}

/*________________________________________________________________________
          Return the Pt Bin for the given Pt
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::GetPtBin(double pt)
{
  ptbin.clear();
  for (int i = 0; i < ptarray.GetSize() - 1; i += 2) {
    if (pt >= ptarray[i] && pt <= ptarray[i + 1])
      ptbin.push_back(kTRUE);
    else
      ptbin.push_back(kFALSE);
  }
}

/*________________________________________________________________________
                  PID Information (only TPC for now)
________________________________________________________________________*/

bool AliAnalysisTaskNFactorialMoments::GetParticleID(AliAODTrack* trk,
                                                         bool fQA)
{

  if (!fismc) {
    double length = -999., beta = -999., tofTime = -999., tof = -999.;
    double c = TMath::C() * 1.E-9; // m/ns

    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();

    float dEdx = trk->GetDetPid()->GetTPCsignal();
    for (int isp(0); isp < AliPID::kSPECIES; isp++) {
      nsigmaTPC[isp] =
        fPIDResponse->NumberOfSigmasTPC(trk, (AliPID::EParticleType)isp);
      nsigmaTOF[isp] =
        fPIDResponse->NumberOfSigmasTOF(trk, (AliPID::EParticleType)isp);
      if (fQA) {
        fHistQAPID[12]->Fill(trk->Pt(), nsigmaTPC[isp]);
      }
    }
    if (fQA) {
      fHistQAPID[0]->Fill(trk->P() * trk->Charge(), dEdx);
    }

    tofTime = trk->GetTOFsignal();
    length = trk->GetIntegratedLength();

    tof = tofTime * 1E-3; // ns
    if (tof <= 0)
      return kTRUE;
    if (length <= 0)
      return kTRUE;

    length = length * 0.01; // in meters
    beta = length / (tof * c);

    if (fProtonAnalysis) {
      fHistQAPID[1]->Fill(trk->Pt(), nsigmaTPC[4]);
      fHistQAPID[2]->Fill(trk->Pt(), nsigmaTOF[4]);
      double combSquare = nsigmaTPC[4] * nsigmaTPC[4] +
                          nsigmaTOF[4] * nsigmaTOF[4];
      fHistQAPID[3]->Fill(trk->Pt(), combSquare);
      fHistQAPID[4]->Fill(trk->Pt(), beta);
      fHistQAPID[5]->Fill(nsigmaTPC[4], nsigmaTOF[4]);

      if (combSquare > nSigmaPrCut) {
        return kTRUE;
      }

      fHistQAPID[6]->Fill(trk->P() * trk->Charge(), dEdx);
      fHistQAPID[7]->Fill(trk->Pt(), nsigmaTPC[4]);
      fHistQAPID[8]->Fill(trk->Pt(), nsigmaTOF[4]);
      fHistQAPID[9]->Fill(trk->Pt(), combSquare);
      fHistQAPID[10]->Fill(trk->Pt(), beta);
      fHistQAPID[11]->Fill(nsigmaTPC[4], nsigmaTOF[4]);
    }

    if (fRejectEls) {
      if (TMath::Abs(nsigmaTPC[0]) < 1 && (trk->Charge()) > 0)
        return kTRUE;
      if (TMath::Abs(nsigmaTPC[0]) < 1 && (trk->Charge()) < 0)
        return kTRUE;
    }
  } else {
    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(
      fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL)
      return kTRUE;
    AliAODMCParticle* particle =
      (AliAODMCParticle*)AODMCTrackArray->At(TMath::Abs(trk->GetLabel()));
    if (!particle)
      return kTRUE;
    if (!particle->IsPhysicalPrimary())
      return kTRUE;
    if (fProtonAnalysis) {
      if (TMath::Abs(particle->GetPdgCode()) == 2212) {
        fHistPDG[0]->Fill(particle->GetPdgCode());
        return kTRUE;
      }
    }
    if (fRejectEls) {
      bool isElPos = kFALSE;
      if (particle->GetPdgCode() == 11 || particle->GetPdgCode() == -11)
        isElPos = kTRUE;
      if (fQA) {
        fHistPDG[0]->Fill(particle->GetPdgCode());
        if (isElPos)
          fHistPDG[1]->Fill(particle->GetPdgCode());
      }
      if (isElPos)
        return kTRUE;
    }
  }
  return kFALSE;
}

/*________________________________________________________________________
              Two Track Calculation of dphistar
________________________________________________________________________*/

float AliAnalysisTaskNFactorialMoments::dphistarcalculation(
  float phi1, float eta1, float pt1, int charge1, float phi2, float eta2,
  float pt2, int charge2, float bSign)
{
  float radius = 0;
  float dphistartemp = 0;
  float dphistar = 999;
  float deta = eta2 - eta1;
  float kLimit = fdeta * 3;
  double kPi = TMath::Pi();
  float deltaPhi = phi2 - phi1;
  bSign = (bSign > 0) ? 1 : -1;
  // variables and cuts have been taken from
  // https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
  if (abs(eta1 - eta2) < fdeta * 2.5 * 3) {
    // check first boundaries to see if is worth to loop and find the minimum
    float dphistar1 =
      GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
    float dphistar2 =
      GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);

    if (fabs(dphistar1) < kLimit || fabs(dphistar2) < kLimit ||
        (dphistar1 * dphistar2) < 0) {
      // find the smallest dphistar (radius of TPC: 0.8 - 2.5 (m))
      for (int rad(80); rad < 251; rad++) {
        radius = rad * 0.01;
        dphistartemp =
          GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, radius, bSign);
        if (rad == 80)
          dphistar = dphistartemp;
        if (dphistartemp < dphistar)
          dphistar = dphistartemp;
      }
      if (dphistar < kPi)
        dphistar = kPi * 2 - dphistar;
      if (dphistar < -1 * kPi)
        dphistar = -2 * kPi - dphistar;
      if (dphistar > kPi)
        dphistar = kPi * 2 - dphistar;
      if (deltaPhi < -0.5 * TMath::Pi())
        deltaPhi += TMath::TwoPi();
      if (deltaPhi > 1.5 * TMath::Pi())
        deltaPhi -= TMath::TwoPi();
      fHistdPhi->Fill(deltaPhi);
    }
    if (dphistar != 999) {
      fHistdEta->Fill(dphistar);
      GetPtBin(pt1);
      if (ptbin[0])
        fHistbeforeHBT[0]->Fill(deta, dphistar);
      if (ptbin[1])
        fHistbeforeHBT[1]->Fill(deta, dphistar);
      if (ptbin[2])
        fHistbeforeHBT[2]->Fill(deta, dphistar);
      if (ptbin[3])
        fHistbeforeHBT[3]->Fill(deta, dphistar);
    }
  }
  return dphistar;
}

/*________________________________________________________________________
          Shared Cluster Fraction Calculation
________________________________________________________________________*/

float AliAnalysisTaskNFactorialMoments::SharedClusterFraction(TBits& cl1,
                                                                  TBits& cl2,
                                                                  TBits& sh1,
                                                                  TBits& sh2)
{
  int ncl1 = cl1.GetNbits();
  int ncl2 = cl2.GetNbits();
  int sumCls = 0;
  int sumSha = 0;
  int sumQ = 0;
  double shfrac = 0;
  for (int ich = 0; ich < ncl1 && ich < ncl2; ich++) {
    if (cl1.TestBitNumber(ich) && cl2.TestBitNumber(ich)) {   // Both clusters
      if (sh1.TestBitNumber(ich) && sh2.TestBitNumber(ich)) { // Shared
        sumQ++;
        sumCls += 2;
        sumSha += 2;
      } else {
        sumQ--;
        sumCls += 2;
      }
    } else if (cl1.TestBitNumber(ich) || cl2.TestBitNumber(ich)) { // Non shared
      sumQ++;
      sumCls++;
    }
  }
  if (sumCls > 0) {
    shfrac = sumSha * 1.0 / sumCls;
  }

  return shfrac;
}

/*________________________________________________________________________
            Main formula for calculation of dphistar
________________________________________________________________________*/

float AliAnalysisTaskNFactorialMoments::GetDPhiStar(float phi1, float pt1,
                                                        float charge1, float phi2,
                                                        float pt2, float charge2,
                                                        float radius, float bSign)
{
  double kPi = TMath::Pi();
  double deltaPhi = phi2 - phi1;
  if (deltaPhi < -0.5 * TMath::Pi())
    deltaPhi += TMath::TwoPi();
  if (deltaPhi > 1.5 * TMath::Pi())
    deltaPhi -= TMath::TwoPi();
  // calculates dphistar
  float dphistarm = deltaPhi +
                    charge1 * bSign * TMath::ASin(0.075 * radius / (2 * pt1)) -
                    charge2 * bSign * TMath::ASin(0.075 * radius / (2 * pt2));

  // circularity
  if (dphistarm > kPi)
    dphistarm = kPi * 2 - dphistarm;
  if (dphistarm < -kPi)
    dphistarm = -kPi * 2 - dphistarm;
  if (dphistarm > kPi) // might look funny but is needed
    dphistarm = kPi * 2 - dphistarm;

  return dphistarm;
}

/*________________________________________________________________________
          Calculation of the Normalized Fq Moments
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::CalculateNFMs(TH2D* h1[M], TH2D* h2[M],
                                                         TH2D* h3[M], TH2D* h4[M],
                                                         bool mcGen)
{
  for (int fptbins = 0; fptbins < 4; fptbins++) {
    for (int binset = 0; binset < M; binset++) {

      double NoOfBins, MSquare;
      if (Mmax == 123)
        NoOfBins = 3 * (binset + 2);
      if (Mmax == 82)
        NoOfBins = 2 * (binset + 2);
      if (fSelfAffine)
        NoOfBins = Mbins[binset] * Nbins[binset];
      MSquare = TMath::Power(NoOfBins, D);
      if (fSelfAffine)
        MSquare = Mbins[binset] * Nbins[binset];
      double SumOfbincontent = 0;
      double FqEvent[Q];
      double sumoff[Q];
      double bincontent;
      double Mbin = NoOfBins;
      int NofXetabins = 0;
      int NofXphibins = 0;

      for (int index = 0; index < Q; index++) {
        FqEvent[index] = 0.0;
        sumoff[index] = 0.0;
      }

      if (fptbins == 0) {
        NofXetabins = h1[binset]->GetNbinsX();
        NofXphibins = h1[binset]->GetNbinsY();
      }
      if (fptbins == 1) {
        NofXetabins = h2[binset]->GetNbinsX();
        NofXphibins = h2[binset]->GetNbinsY();
      }
      if (fptbins == 2) {
        NofXetabins = h3[binset]->GetNbinsX();
        NofXphibins = h3[binset]->GetNbinsY();
      }
      if (fptbins == 3) {
        NofXetabins = h4[binset]->GetNbinsX();
        NofXphibins = h4[binset]->GetNbinsY();
      }

      for (int etabin = 1; etabin <= NofXetabins; etabin++) {
        for (int phibin = 1; phibin <= NofXphibins; phibin++) {
          bincontent = 0.0;
          if (fptbins == 0)
            bincontent = h1[binset]->GetBinContent(etabin, phibin);
          if (fptbins == 1)
            bincontent = h2[binset]->GetBinContent(etabin, phibin);
          if (fptbins == 2)
            bincontent = h3[binset]->GetBinContent(etabin, phibin);
          if (fptbins == 3)
            bincontent = h4[binset]->GetBinContent(etabin, phibin);
          SumOfbincontent += bincontent;

          for (int q = 0; q < Q; q++) {
            if (bincontent >= (q + 2)) {
              double Fqeofbin = 0.0;
              Fqeofbin = TMath::Factorial(bincontent) /
                         TMath::Factorial(bincontent - (q + 2));
              if (fismc)
                if (TMath::IsNaN(Fqeofbin))
                  break;
              sumoff[q] += Fqeofbin;
            }
          }
        }
      }

      double Av_bincontent = SumOfbincontent / MSquare;
      for (int q = 0; q < Q; q++) {
        if (sumoff[q] > 0.0)
          FqEvent[q] = sumoff[q] / (MSquare);
      }

      float Fq2e = FqEvent[0];
      float Fq3e = FqEvent[1];
      float Fq4e = FqEvent[2];
      float Fq5e = FqEvent[3];
      float Fq6e = FqEvent[4];
      float Fq7e = FqEvent[5];

      if (fptbins == 0) {
        if (!mcGen)
          fntpMBin1[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e,
                                  Fq6e, Fq7e);
        else
          fntpMBin1Gen[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e,
                                     Fq5e, Fq6e, Fq7e);
      }
      if (fptbins == 1) {
        if (!mcGen)
          fntpMBin2[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e,
                                  Fq6e, Fq7e);
        else
          fntpMBin2Gen[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e,
                                     Fq5e, Fq6e, Fq7e);
      }
      if (fptbins == 2) {
        if (!mcGen)
          fntpMBin3[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e,
                                  Fq6e, Fq7e);
        else
          fntpMBin3Gen[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e,
                                     Fq5e, Fq6e, Fq7e);
      }
      if (fptbins == 3) {
        if (!mcGen)
          fntpMBin4[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e,
                                  Fq6e, Fq7e);
        else
          fntpMBin4Gen[binset]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e,
                                     Fq5e, Fq6e, Fq7e);
      }
    }
  } // end of calculation loop
}

/*________________________________________________________________________
          Post data to the output slot(s)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::DataPosting()
{
  PostData(1, fOutHList);
  PostData(2, fNtupleListBin1);
  PostData(3, fNtupleListBin2);
  PostData(4, fNtupleListBin3);
  PostData(5, fNtupleListBin4);

  if (fismc) {
    PostData(6, fNtupleListBin1Gen);
    PostData(7, fNtupleListBin2Gen);
    PostData(8, fNtupleListBin3Gen);
    PostData(9, fNtupleListBin4Gen);
    PostData(10, fQAList2);
  } else {
    PostData(6, fQAList2);
    PostData(7, fQAList);
  }
}

/*________________________________________________________________________
                                Terminate
________________________________________________________________________*/
void AliAnalysisTaskNFactorialMoments::Terminate(Option_t*)
{
  Info("AliAnalysisTaskNFactorialMoments", "Task Successfully finished");
  // terminate
  // called at the END of the analysis (when all events are processed)
}
