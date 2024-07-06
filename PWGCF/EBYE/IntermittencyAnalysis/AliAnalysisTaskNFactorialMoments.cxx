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
  : AliAnalysisTaskSE(), fAOD(0), fMCEvent(0), fHistList(0), fEventCuts(0), fPIDResponse(0), fPIDCombined(0), mcHeader{ nullptr }, fNtupleList{ nullptr }, fNtupleListGen{ nullptr }, fHistQAEta{ nullptr }, fHistQAPhi{ nullptr }, fHistQAVx(0), fHistdEta{ nullptr }, fHistdPhi{ nullptr }, fHistQAVy(0), counter(0), fHistQAVz(0), fPtBin{ nullptr }, fEtaBin{ nullptr }, fPhiBin{ nullptr }, fMultBin{ nullptr }, fPtBinGen{ nullptr }, fEtaBinGen{ nullptr }, fPhiBinGen{ nullptr }, fMultBinGen{ nullptr }, fdEtadPhiBef{ nullptr }, fdEtadPhiAf{ nullptr }, fdEtadPhiChSame{ nullptr }, fdEtadPhiChDiff{ nullptr }, fdEtadPhiPtOrd{ nullptr }, fHistQAPID{ nullptr }, fHistPDG{ nullptr }, fHistDCAxy{ nullptr }, fHistDCAz{ nullptr }, fHistnITScls{ nullptr }, fHistnTPCcls{ nullptr }, fHistnTPCcrossedrows{ nullptr }, fHistnchi2ITScls{ nullptr }, fHistnchi2TPCcls{ nullptr }, fHistDCAxypT{ nullptr }, fHistDCAzpT{ nullptr }, fEtaPhiBin{ nullptr }, fntpMBin{ nullptr }, fntpMBinGen{ nullptr }, fEtaPhiBinGen{ nullptr }, fHistNShCls{ nullptr }, fHistNShClsFra{ nullptr }, fHistNShClsFravspt{ nullptr }, fHistNFoundClsFra{ nullptr }, fHistNFcls{ nullptr }, fEventCounter(0), fTrackCounter(0), fHistQACent(0), ptbin1(0), ptbin2(0), ptbin3(0), ptbin4(0), mfield(0) {}
AliAnalysisTaskNFactorialMoments::AliAnalysisTaskNFactorialMoments(
  const char* name)
  : AliAnalysisTaskSE(name), fAOD(0), fMCEvent(0), fHistList(0), fEventCuts(0), fPIDResponse(0), fPIDCombined(0), mcHeader{ nullptr }, fNtupleList{ nullptr }, fNtupleListGen{ nullptr }, fHistQAEta{ nullptr }, fHistQAPhi{ nullptr }, fHistQAVx(0), fHistdEta{ nullptr }, fHistdPhi{ nullptr }, fHistQAVy(0), counter(0), fHistQAVz(0), fPtBin{ nullptr }, fEtaBin{ nullptr }, fPhiBin{ nullptr }, fMultBin{ nullptr }, fPtBinGen{ nullptr }, fEtaBinGen{ nullptr }, fPhiBinGen{ nullptr }, fMultBinGen{ nullptr }, fdEtadPhiBef{ nullptr }, fdEtadPhiAf{ nullptr }, fdEtadPhiChSame{ nullptr }, fdEtadPhiChDiff{ nullptr }, fdEtadPhiPtOrd{ nullptr }, fHistQAPID{ nullptr }, fHistPDG{ nullptr }, fHistDCAxy{ nullptr }, fHistDCAz{ nullptr }, fHistnITScls{ nullptr }, fHistnTPCcls{ nullptr }, fHistnTPCcrossedrows{ nullptr }, fHistnchi2ITScls{ nullptr }, fHistnchi2TPCcls{ nullptr }, fHistDCAxypT{ nullptr }, fHistDCAzpT{ nullptr }, fEtaPhiBin{ nullptr }, fEtaPhiBinGen{ nullptr }, fntpMBin{ nullptr }, fntpMBinGen{ nullptr }, fHistNShCls{ nullptr }, fHistNShClsFra{ nullptr }, fHistNShClsFravspt{ nullptr }, fHistNFoundClsFra{ nullptr }, fHistNFcls{ nullptr }, fEventCounter(0), fTrackCounter(0), fHistQACent(0), ptbin1(0), ptbin2(0), ptbin3(0), ptbin4(0), mfield(0)
{
  // Constructor
  Info("AliAnalysisTaskNFactorialMoments", "Specific Constructor");

  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
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
  if (flagMC) {
    DefineOutput(8, TList::Class());
    DefineOutput(9, TList::Class());
    DefineOutput(10, TList::Class());
    DefineOutput(11, TList::Class());
  }
}

AliAnalysisTaskNFactorialMoments::~AliAnalysisTaskNFactorialMoments()
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

void AliAnalysisTaskNFactorialMoments::UserCreateOutputObjects()
{
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);

  for (Int_t i = 0; i < mPtBins; ++i) {
    fNtupleList[i] = new TList();
    fNtupleList[i]->SetOwner(kTRUE);
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
  fTrackCounter->GetXaxis()->SetBinLabel(3, "Rej PIDs");
  fTrackCounter->GetXaxis()->SetBinLabel(4, "Loose cuts");
  fTrackCounter->GetXaxis()->SetBinLabel(5, "Systematic cuts");
  fTrackCounter->GetXaxis()->SetBinLabel(6, "DCAs");
  fTrackCounter->GetXaxis()->SetBinLabel(7, "nSharedCls/nCls");
  fTrackCounter->GetXaxis()->SetBinLabel(8, "nSharedCls/nCrRows");
  fTrackCounter->GetXaxis()->SetBinLabel(9, "nCrRows/nFindableCls");
  fTrackCounter->GetXaxis()->SetBinLabel(10, "TwoTracks");
  for (Int_t i = 0; i < mPtBins; ++i) {
    fTrackCounter->GetXaxis()->SetBinLabel(11 + i, Form("p_{T}_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]));
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

  for (Int_t i = 0; i < 2; ++i) {
    std::array<TString, 2> j = { "before", "after" };
    fHistDCAxy[i] = new TH1F(Form("fHistDCAxy%s", j[i].Data()), Form("DCAxy distribution %s;DCAxy;Counts", j[i].Data()), 900, -4.5, 4.5);
    fHistDCAz[i] = new TH1F(Form("fHistDCAz%s", j[i].Data()), Form("DCAz distribution %s;DCAz;Counts", j[i].Data()), 900, -4.5, 4.5);
    fHistnITScls[i] = new TH1F(Form("fHistnITScls%s", j[i].Data()), Form("ITS cluster distribution %s;ITS cluster;Counts", j[i].Data()), 10, -0.5, 10.5);
    fHistnTPCcls[i] = new TH1F(Form("fHistnTPCcls%s", j[i].Data()), Form("TPC cluster distribution %s;TPC cluster;Counts", j[i].Data()), 200, 0, 200);
    fHistnTPCcrossedrows[i] = new TH1F(Form("fHistnTPCcrossedrows%s", j[i].Data()), Form("TPC crossed rows distribution %s;TPC crossed rows;Counts", j[i].Data()), 200, 0, 200);
    fHistnchi2ITScls[i] = new TH1F(Form("fHistnchi2ITScls%s", j[i].Data()), Form("ITSchi2 cluster distribution %s;ITSchi2 cluster;Counts", j[i].Data()), 100, -0.5, 100.5);
    fHistnchi2TPCcls[i] = new TH1F(Form("fHistnchi2TPCcls%s", j[i].Data()), Form("TPCchi2 cluster distribution %s;TPCchi2 cluster;Counts", j[i].Data()), 50, 0, 50);
    fHistDCAxypT[i] = new TH2F(Form("fHistDCAxypT%s", j[i].Data()), Form("DCAxy vs pT %s;p_{T};dca_{xy}", j[i].Data()), 500, 0, 5.5, 400, -4.0, 4.0);
    fHistDCAzpT[i] = new TH2F(Form("fHistDCAzpT%s", j[i].Data()), Form("DCAz vs pT %s;p_{T};dca_{z}", j[i].Data()), 500, 0, 5.5, 400, -4.0, 4.0);
    fHistNShCls[i] = new TH1F(Form("fHistNShCls%s", j[i].Data()), Form("TPC shared cluster distribution %s", j[i].Data()), 200, 0, 200);
    fHistNShClsFra[i] = new TH2F(Form("fHistNShClsFra%s", j[i].Data()), Form("TPC shared cluster fraction %s;sharedcls/ncls;sharedcls/ncrrows", j[i].Data()), 100, 0, 1.0, 100, 0, 1.0);
    fHistNFoundClsFra[i] = new TH2F(Form("fHistNFoundClsFra%s", j[i].Data()), Form("TPC found cluster fraction %s;sharedcls/ncls;ncrrows/findablecls", j[i].Data()), 200, 0, 2, 200, 0, 2);
    fHistNFcls[i] = new TH1F(Form("fHistNFcls%s", j[i].Data()), Form("TPC findable cluster distribution %s;ncrrows/findablecls;counts", j[i].Data()), 40, 0, 4);
    fHistList->Add(fHistDCAxy[i]);
    fHistList->Add(fHistDCAz[i]);
    fHistList->Add(fHistnITScls[i]);
    fHistList->Add(fHistnTPCcls[i]);
    fHistList->Add(fHistnTPCcrossedrows[i]);
    fHistList->Add(fHistnchi2ITScls[i]);
    fHistList->Add(fHistnchi2TPCcls[i]);
    fHistList->Add(fHistDCAxypT[i]);
    fHistList->Add(fHistDCAzpT[i]);
    fHistList->Add(fHistNShCls[i]);
    fHistList->Add(fHistNShClsFra[i]);
    fHistList->Add(fHistNFoundClsFra[i]);
    fHistList->Add(fHistNFcls[i]);
  }
  fHistNShClsFravspt[mPtBins] = new TH2F("fHistNShClsFraVsPtBef", "TPC nshared clusters vs #it{p_{T}} before", 500, 0, 5, 400, 0.02, 1);
  fHistList->Add(fHistNShClsFravspt[4]);
  for (Int_t i = 0; i < mPtBins; ++i) {
    fHistNShClsFravspt[i] = new TH2F(Form("fHistNShClsFraVsPt_%.2f-%.2f", ptarray[i], ptarray[i + 1]), Form("TPC nshared clusters vs #it{p_{T}} %i", i), 500, 0, 5, 400, 0.02, 1);
    fHistList->Add(fHistNShClsFravspt[i]);
  }

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
      fMultBinGen[i] = new TH1F(Form("fMultBinGen_%.2f-%.2f", ptarray[2 * i], ptarray[2 * i + 1]), Form("Multiplicity distribution of Gen tracks for PtBin%d;Multiplicity;Counts", i + 1), 10000, 0.0, 10000.0);
      fHistList->Add(fPtBinGen[i]);
      fHistList->Add(fEtaBinGen[i]);
      fHistList->Add(fPhiBinGen[i]);
      fHistList->Add(fMultBinGen[i]);
    }
    for (Int_t j = 0; j < mMBins; ++j) {
      // Eta-Phi distributions
      fEtaPhiBin[i][j] = new TH2D(Form("fEtaPhiBin%s_%.2f-%.2f_%d", name.Data(), ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mMBin2[j], 0.0, TMath::TwoPi());
      if (flagSelfAff) {
        fEtaPhiBin[i][j] = new TH2D(Form("fEtaPhiBin%s_%.2f-%.2f_%d", name.Data(), ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mNBin2[j], 0.0, TMath::TwoPi());
      }
      if (flagMC) {
        fEtaPhiBinGen[i][j] = new TH2D(Form("fEtaPhiBinGen_%.2f-%.2f_%d", ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of Gen tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mMBin2[j], 0.0, TMath::TwoPi());
        if (flagSelfAff) {
          fEtaPhiBinGen[i][j] = new TH2D(Form("fEtaPhiBinGen_%.2f-%.2f_%d", ptarray[2 * i], ptarray[2 * i + 1], j + 1), Form("Eta-Phi distribution of Gen tracks for PtBin%d and MBin%d;#eta;#phi", i + 1, j + 1), mMBin2[j], minEta, maxEta, mNBin2[j], 0.0, TMath::TwoPi());
        }
        fHistList->Add(fEtaPhiBinGen[i][j]);
      }
      fHistList->Add(fEtaPhiBin[i][j]);
      // Tuples for each bin
      fntpMBin[i][j] = new TNtuple(Form("fntpMBin%d_%d%s", i + 1, j + 1, name.Data()), Form("fntpMBin%d_%d", i + 1, j + 1), "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
      fNtupleList[i]->Add(fntpMBin[i][j]);
      if (flagMC) {
        fntpMBinGen[i][j] = new TNtuple(Form("fntpMBin%d_%dGen", i + 1, j + 1), Form("fntpMBin%d_%dGen", i + 1, j + 1), "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
        fNtupleListGen[i]->Add(fntpMBinGen[i][j]);
      }
    }
  }
  // QA Histograms for detadphi and PID
  fHistdEta = new TH1D("dEta", "dEta;#Delta#eta;Counts", 100, -1.75, 1.75);
  fHistdPhi = new TH1D("dphi", "dPhi;#Delta#phi;Counts", 100, -TMath::Pi(), TMath::Pi());
  fHistList->Add(fHistdEta);
  fHistList->Add(fHistdPhi);

  TString charge[3] = { "Pos", "Neg", "All" };
  TString detadphiRange[2] = { "", "Close" };
  TString ptOrd[4] = { "pt2>pt1_dphi<0", "pt2>pt1_dphi>0", "pt2<pt1_dphi<0", "pt2<pt1_dphi>0" };

  for (Int_t i = 0; i < mPtBins; ++i) {
    fdEtadPhiBef[i] = new TH2D(Form("fdEtadPhiBef%d", i), Form("dEta dPhi for ptbin %d;#Delta#eta;#Delta#phi", i), 35, -1.75, 1.75, 50, -TMath::Pi(), TMath::Pi());
    fdEtadPhiAf[i] = new TH2D(Form("fdEtadPhiAf%d", i), Form("dEta dPhi for ptbin %d;#Delta#eta;#Delta#phi", i), 50, -1.75, 1.75, 50, -TMath::Pi(), TMath::Pi());
    fHistList->Add(fdEtadPhiBef[i]);
    fHistList->Add(fdEtadPhiAf[i]);
    for (Int_t ch = 0; ch < 3; ++ch) {
      fdEtadPhiChSame[i][0][ch] = new TH2D(Form("fdEtadPhiCh_%.2f-%.2f_%s%s", ptarray[2 * i], ptarray[2 * i + 1], charge[ch].Data(), detadphiRange[0].Data()), Form("dEta dPhi for ptbin %d %s;#Delta#eta;#Delta#phi", i, charge[ch].Data()), 35, -1.75, 1.75, 50, -TMath::Pi(), TMath::Pi());
      fdEtadPhiChSame[i][1][ch] = new TH2D(Form("fdEtadPhiCh_%.2f-%.2f_%s%s", ptarray[2 * i], ptarray[2 * i + 1], charge[ch].Data(), detadphiRange[1].Data()), Form("dEta dPhi for ptbin %d %s;#Delta#eta;#Delta#phi", i, charge[ch].Data()), 50, -0.02, 0.02, 50, -0.05, 0.05);
      fHistList->Add(fdEtadPhiChSame[i][0][ch]);
      fHistList->Add(fdEtadPhiChSame[i][1][ch]);
      for (Int_t j = 0; j < 4; ++j) {
        fdEtadPhiPtOrd[i][j][ch][0] = new TH2D(Form("fdEtadPhiPtOrd_%.2f-%.2f_%s_%s%s", ptarray[2 * i], ptarray[2 * i + 1], ptOrd[j].Data(), charge[ch].Data(), detadphiRange[0].Data()), Form("dEta dPhi for ptbin %d %s %s;#Delta#eta;#Delta#phi", i, ptOrd[j].Data(), charge[ch].Data()), 35, -1.75, 1.75, 50, -TMath::Pi(), TMath::Pi());
        fdEtadPhiPtOrd[i][j][ch][1] = new TH2D(Form("fdEtadPhiPtOrd_%.2f-%.2f_%s_%s%s", ptarray[2 * i], ptarray[2 * i + 1], ptOrd[j].Data(), charge[ch].Data(), detadphiRange[1].Data()), Form("dEta dPhi for ptbin %d %s %s;#Delta#eta;#Delta#phi", i, ptOrd[j].Data(), charge[ch].Data()), 50, -0.02, 0.02, 50, -0.05, 0.05);
        fHistList->Add(fdEtadPhiPtOrd[i][j][ch][0]);
        fHistList->Add(fdEtadPhiPtOrd[i][j][ch][1]);
      }
    }
    fdEtadPhiChDiff[i][0] = new TH2D(Form("fdEtadPhiChDiff_%.2f-%.2f%s", ptarray[2 * i], ptarray[2 * i + 1], detadphiRange[0].Data()), Form("dEta dPhi for ptbin %d different ch;#Delta#eta;#Delta#phi", i), 35, -1.75, 1.75, 50, -TMath::Pi(), TMath::Pi());
    fdEtadPhiChDiff[i][1] = new TH2D(Form("fdEtadPhiChDiff_%.2f-%.2f%s", ptarray[2 * i], ptarray[2 * i + 1], detadphiRange[1].Data()), Form("dEta dPhi for ptbin %d different ch;#Delta#eta;#Delta#phi", i), 50, -0.02, 0.02, 50, -0.05, 0.05);
    fHistList->Add(fdEtadPhiChDiff[i][0]);
    fHistList->Add(fdEtadPhiChDiff[i][1]);
  }

  fHistQAPID[0] = new TH2D("TPCdEdx", "TPCdEdx vs momentum;#it{p};TPCdEdx", 100, 0.0, 5.0, 250, 0.0, 500.0);
  fHistQAPID[1] = new TH2D("TPCNsigma", "TPCNsigma vs momentum;#it{p};TPCNsigma", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fHistQAPID[2] = new TH2D("TOFNsigmaProton", "TOFNsigmaProton vs momentum;#it{p};TOFNsigmaProton", 100, 0.0, 5.0, 250, -5.0, 5.0);
  fHistQAPID[3] = new TH2D("TPCTOFProton", "TPCTOFProton vs momentum;#it{p};TPCTOFProton", 100, 0.0, 5.0, 250, -5.0, 5.0);
  fHistQAPID[4] = new TH2D("TOF signal", "momentum vs beta;#it{p};#beta", 1000, 0.5, 5.0, 1000, 0.4, 1.2);
  fHistQAPID[5] = new TH2D("nSigmaTPC  vs nSigmaTOF", "nSigmaTPC  vs nSigmaTOF;nSigmaTPC;nSigmaTOF", 1000, -5.0, 5.0, 1000, -5.0, 5.0);
  fHistQAPID[6] = new TH2D("TPCdEdxAfter", "TPCdEdx vs momentum;#it{p};TPCdEdx", 100, 0.0, 10.0, 250, 0.0, 250.0);
  fHistQAPID[7] = new TH2D("TPCNsigmaProtonAfter", "TPCNsigmaProton vs momentum;#it{p};TPCNsigmaProton", 100, 0.0, 2.0, 100, -5.0, 5.0);
  fHistQAPID[8] = new TH2D("TOFNsigmaProtonAfter", "TOFNsigmaProton vs momentum;#it{p};TOFNsigmaProton", 100, 0.0, 5.0, 250, -5.0, 5.0);
  fHistQAPID[9] = new TH2D("TPCTOFProtonAfter", "TPCTOFProton vs momentum;#it{p};TPCTOFProton", 100, 0.0, 5.0, 250, -5.0, 5.0);
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

  if (!eventHandler || !fAOD) {
    std::cout << "UserExec: "
              << "InputHandler or PIDResponse or AODEvent not found"
              << std::endl;
    return;
  }

  if (flagMC) {
    fMCEvent = MCEvent();
    mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(
      AliAODMCHeader::StdBranchName());
    if (!fMCEvent || !mcHeader) {
      std::cout << "UserExec: "
                << "MCEvent or MCHeader not found"
                << std::endl;
      return;
    }
  }

  // fEventCuts.SetupLHC15o();
  if (fYear == "2015") {
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
  }

  fEventCounter->Fill(2);

  if (!flagMC) {
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

void AliAnalysisTaskNFactorialMoments::FillTrackInfo()
{
  Int_t counterBin[mPtmax] = { 0 };
  Int_t nTracks(fAOD->GetNumberOfTracks());
  Float_t dpstar, deta, dphi;
  Float_t dcaXY, dcaZ;
  Bool_t flagPID = kFALSE;

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
        continue;
      AliAODMCParticle* particle =
        (AliAODMCParticle*)AODMCTrackArray->At(TMath::Abs(track->GetLabel()));
      if (!particle)
        continue;
      if (flagRejEls) {
        Bool_t isElPos = kFALSE;
        if (particle->GetPdgCode() == 11 || particle->GetPdgCode() == -11) {
          isElPos = kTRUE;
        }
        fHistPDG[0]->Fill(particle->GetPdgCode());
        if (isElPos) {
          fHistPDG[1]->Fill(particle->GetPdgCode());
          continue;
        }
      }
    }

    Double_t length = -999., beta = -999., tofTime = -999., tof = -999.;
    Double_t c = TMath::C() * 1.E-9; // m/ns
    Float_t dEdx;

    // TODO: commented out the below statement because 128 seems to not give any values in the PID object

    // NOTE: These are all only PID related QA histograms (mostly for protons). No checks done
    if (flagPID) {
      dEdx = track->GetTPCsignal();
      for (Int_t isp(0); isp < AliPID::kSPECIES; isp++) {
        nsigmaTPC[isp] =
          fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)isp);
        nsigmaTOF[isp] =
          fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)isp);
        fHistQAPID[12]->Fill(pt, nsigmaTPC[isp]);
      }
      fHistQAPID[0]->Fill(track->P() * track->Charge(), dEdx);
      tofTime = track->GetTOFsignal();
      length = track->GetIntegratedLength();
      tof = tofTime * 1E-3; // ns
      if (tof > 0 && length > 0) {
        length = length * 0.01; // in meters
        beta = length / (tof * c);
        fHistQAPID[1]->Fill(track->Pt(), nsigmaTPC[4]);
        fHistQAPID[2]->Fill(track->Pt(), nsigmaTOF[4]);
        Double_t combSquare = TMath::Sqrt(nsigmaTPC[4] * nsigmaTPC[4] +
                                          nsigmaTOF[4] * nsigmaTOF[4]);
        fHistQAPID[3]->Fill(track->Pt(), combSquare);
        fHistQAPID[4]->Fill(track->Pt(), beta);
        fHistQAPID[5]->Fill(nsigmaTPC[4], nsigmaTOF[4]);
        if ((track->Pt() < 0.6) && (TMath::Abs(nsigmaTPC[4]) < 3)) {
          fHistQAPID[6]->Fill(track->P() * track->Charge(), dEdx);
          fHistQAPID[7]->Fill(track->Pt(), nsigmaTPC[4]);
          fHistQAPID[8]->Fill(track->Pt(), nsigmaTOF[4]);
          fHistQAPID[9]->Fill(track->Pt(), combSquare);
          fHistQAPID[10]->Fill(track->Pt(), beta);
          fHistQAPID[11]->Fill(nsigmaTPC[4], nsigmaTOF[4]);
        } else if ((track->Pt() > 0.6) && (combSquare < 3)) {
          fHistQAPID[6]->Fill(track->P() * track->Charge(), dEdx);
          fHistQAPID[7]->Fill(track->Pt(), nsigmaTPC[4]);
          fHistQAPID[8]->Fill(track->Pt(), nsigmaTOF[4]);
          fHistQAPID[9]->Fill(track->Pt(), combSquare);
          fHistQAPID[10]->Fill(track->Pt(), beta);
          fHistQAPID[11]->Fill(nsigmaTPC[4], nsigmaTOF[4]);
        }
      }
      if (flagRejEls) {
        if (TMath::Abs(nsigmaTPC[0]) < 1)
          continue;
      }
    }

    fTrackCounter->Fill(3);

    // NOTE: Some loose cuts on the tracks
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
    Double_t dca[2];

    if (!(GetDCA(track, dca)))
      continue;
    dcaXY = dca[0];
    dcaZ = dca[1];

    if ((fITSCls > 0.0) && ((track->GetITSchi2() / track->GetITSNcls()) > fITSCls))
      continue;

    if ((fTPCCls > 0.0) && (track->GetTPCNcls() < fTPCCls))
      continue;

    if ((fTPCRows > 0.0) && (track->GetTPCNCrossedRows() < fTPCRows))
      continue;

    fTrackCounter->Fill(5);

    fHistnITScls[0]->Fill(track->GetITSNcls());
    fHistnTPCcls[0]->Fill(track->GetTPCNcls());
    fHistnchi2ITScls[0]->Fill((track->GetITSchi2() / track->GetITSNcls()));
    fHistnchi2TPCcls[0]->Fill((track->GetTPCchi2() / track->GetTPCNcls()));
    fHistnTPCcrossedrows[0]->Fill(track->GetTPCNCrossedRows());
    fHistDCAxy[0]->Fill(dcaXY);
    fHistDCAz[0]->Fill(dcaZ);
    fHistDCAxypT[0]->Fill(pt, dcaXY);
    fHistDCAzpT[0]->Fill(pt, dcaZ);

    // NOTE: DCAxy cut ==  0 is loose, 1 is tight
    if ((fDCAxyMax == 0) && (fabs(dcaXY) > (0.0182 + 0.4750 / TMath::Power(pt, 1.95))))
      continue;
    if ((fDCAxyMax == 1) && (fabs(dcaXY) > (0.0198 + 0.3963 / TMath::Power(pt, 1.59))))
      continue;
    if ((fDCAzMax > 0.0) && (fabs(dcaZ) > fDCAzMax))
      continue;

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

    fHistnITScls[1]->Fill(track->GetITSNcls());
    fHistnTPCcls[1]->Fill(track->GetTPCNcls());
    fHistnchi2ITScls[1]->Fill((track->GetITSchi2() / track->GetITSNcls()));
    fHistnchi2TPCcls[1]->Fill((track->GetTPCchi2() / track->GetTPCNcls()));
    fHistnTPCcrossedrows[1]->Fill(track->GetTPCNCrossedRows());
    fHistDCAxy[1]->Fill(dcaXY);
    fHistDCAz[1]->Fill(dcaZ);
    fHistDCAxypT[1]->Fill(pt, dcaXY);
    fHistDCAzpT[1]->Fill(pt, dcaZ);
    fHistNShCls[1]->Fill(nSharedCls);
    fHistNFcls[1]->Fill(nCrossedRows / nFindableCls);
    fHistNShClsFra[1]->Fill(nSharedCls / nCls, nSharedCls / nCrossedRows);
    fHistNFoundClsFra[1]->Fill(nSharedCls / nCls, nCrossedRows / nFindableCls);

    // TODO: These values are calculated for fb768. Have to do these for fb128 as well.
    Double_t profileVal[16] = { 0 };
    if (minCent == 0 && maxCent == 5) {
      Double_t profileVal1[16] = { 0.360756, 0.360834, 0.360455, 0.35992, 0.359485, 0.359125, 0.359184, 0.359317, 0.359458, 0.35959, 0.360105, 0.36066, 0.361126, 0.361661, 0.361998, 0.362225 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 5 && maxCent == 10) {
      Double_t profileVal1[16] = { 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379, 0.358205 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 10 && maxCent == 20) {
      Double_t profileVal1[16] = { 0.35707, 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 20 && maxCent == 40) {
      Double_t profileVal1[16] = { 0.35707, 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else if (minCent == 40 && maxCent == 80) {
      Double_t profileVal1[16] = { 0.35707, 0.357014, 0.356723, 0.356252, 0.355931, 0.355681, 0.35584, 0.356167, 0.356138, 0.356158, 0.356624, 0.357092, 0.357543, 0.357907, 0.35816, 0.358379 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    } else {
      Double_t profileVal1[16] = { 0.360756, 0.360834, 0.360455, 0.35992, 0.359485, 0.359125, 0.359184, 0.359317, 0.359458, 0.35959, 0.360105, 0.36066, 0.361126, 0.361661, 0.361998, 0.362225 };
      std::copy(profileVal1, profileVal1 + 16, profileVal);
    }

    fHistNShClsFravspt[mPtBins]->Fill(pt, nSharedCls / nCrossedRows);

    for (Int_t k = 0; k < 16; ++k) {
      if (pt > 0.4 + (k * 0.1) && pt < 0.5 + (k * 0.1)) {
        if (flagShClPro && (nSharedCls / nCrossedRows > profileVal[k]))
          continue;
      }
    }

    // Two Track Loop:
    // NOTE: No cuts for d\eta and d\phi are applied here. These are only QA histograms to check.
    if (flag2TrackQA) {
      TBits clusmap = track->GetTPCClusterMap();
      TBits sharedmap = track->GetTPCSharedMap();

      for (Int_t j = 0; j < nTracks; j++) {
        AliAODTrack* track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
        if (!(track2->TestFilterBit(filterBit)))
          continue;

        Int_t charge2 = track2->Charge();
        Float_t pt2 = track2->Pt();
        Float_t eta2 = track2->Eta();
        Float_t phi2 = track2->Phi();
        Int_t id2 = track2->GetID();

        if (flagMC) {
          TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(
            fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
          if (AODMCTrackArray == NULL)
            continue;
          AliAODMCParticle* particle =
            (AliAODMCParticle*)AODMCTrackArray->At(TMath::Abs(track2->GetLabel()));
          if (!particle)
            continue;
          if (flagRejEls) {
            Bool_t isElPos = kFALSE;
            if (particle->GetPdgCode() == 11 || particle->GetPdgCode() == -11) {
              isElPos = kTRUE;
            }
            if (isElPos) {
              continue;
            }
          }
        }

        if (id == id2 || (eta2 < minEta) || (eta2 > maxEta) || fabs(pt2) < 0.2 || charge2 == 0)
          continue;
        if (flagSharity) {
          TBits clusmap2 = track2->GetTPCClusterMap();
          TBits sharedmap2 = track2->GetTPCSharedMap();
          Float_t sharity =
            SharedClusterFraction(clusmap, clusmap2, sharedmap, sharedmap2);
          if (sharity > fSharedFraction) {
            mSharedTrack = kTRUE;
          }
        }
        if ((flagSharity) && (mSharedTrack)) {
          break;
        }

        dpstar = CalculateDPhiStar(phi, eta, pt, charge, phi2, eta2, pt2,
                                   charge2, mfield);
        deta = eta2 - eta;
        dphi = phi2 - phi;
        if (dphi > TMath::Pi())
          dphi -= TMath::TwoPi();
        else if (dphi < -TMath::Pi())
          dphi += TMath::TwoPi();
        fHistdEta->Fill(deta);
        fHistdPhi->Fill(dphi);

        Int_t ch = -999;
        if (charge > 0 && charge2 > 0)
          ch = 0;
        if (charge < 0 && charge2 < 0)
          ch = 1;

        for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
          if (ptbin[iPt] && pt > ptarray[2 * iPt] && pt < ptarray[2 * iPt + 1]) {
            fdEtadPhiBef[iPt]->Fill(deta, dphi);
            for (Int_t j = 0; j < 2; ++j) {
              if (charge2 != charge) {
                fdEtadPhiChDiff[iPt][j]->Fill(deta, dphi);
              }
              if (ch != -999) {
                fdEtadPhiChSame[iPt][j][ch]->Fill(deta, dphi);
                fdEtadPhiChSame[iPt][j][2]->Fill(deta, dphi);
                if (pt2 > pt && dphi < 0) {
                  fdEtadPhiPtOrd[iPt][0][ch][j]->Fill(deta, dphi);
                  fdEtadPhiPtOrd[iPt][0][2][j]->Fill(deta, dphi);
                }
                if (pt2 > pt && dphi > 0) {
                  fdEtadPhiPtOrd[iPt][1][ch][j]->Fill(deta, dphi);
                  fdEtadPhiPtOrd[iPt][1][2][j]->Fill(deta, dphi);
                }
                if (pt2 < pt && dphi < 0) {
                  fdEtadPhiPtOrd[iPt][2][ch][j]->Fill(deta, dphi);
                  fdEtadPhiPtOrd[iPt][2][2][j]->Fill(deta, dphi);
                }
                if (pt2 < pt && dphi > 0) {
                  fdEtadPhiPtOrd[iPt][3][ch][j]->Fill(deta, dphi);
                  fdEtadPhiPtOrd[iPt][3][2][j]->Fill(deta, dphi);
                }
              }
            }
          }
        }
        if (dpstar == 999)
          continue;
        if (fabs(deta) < fdeta && fabs(dpstar) < fdphi && charge == charge2)
          mSkipTracks = kTRUE;
        if ((flag2Track) && (mSkipTracks)) {
          break;
        } else {
          for (Int_t k = 0; k < mPtBins; ++k) {
            if (ptbin[k]) {
              fdEtadPhiAf[k]->Fill(deta, dphi);
            }
          }
        }
      } // End of Two Track Loop

      if (((flag2Track) && (mSharedTrack)) || ((flagSharity) && (mSharedTrack))) {
        continue;
      }
    } // End of Two Track QA
    fTrackCounter->Fill(10);
    fHistQAEta[0]->Fill(eta);
    fHistQAPhi[0]->Fill(phi);

    // NOTE: This loop fills all the eta-phi histograms in one go, avoids multiple loops.
    for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
      if (ptbin[iPt]) {
        fHistNShClsFravspt[iPt]->Fill(pt, nSharedCls / nCrossedRows);
        fTrackCounter->Fill(11 + iPt);
        counterBin[iPt]++;
        fPtBin[iPt]->Fill(pt);
        fEtaBin[iPt]->Fill(eta);
        fPhiBin[iPt]->Fill(phi);
        for (Int_t k = 0; k < mMBins; ++k) {
          fEtaPhiBin[iPt][k]->Fill(eta, phi);
        }
      }
    }
  }
  for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
    fMultBin[iPt]->Fill(counterBin[iPt]);
  }
}

/*________________________________________________________________________
          Resetting Histograms after each event
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::ResetHistograms()
{
  for (Int_t i = 0; i < mPtBins; ++i) {
    for (Int_t j = 0; j < mMBins; ++j) {
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

void AliAnalysisTaskNFactorialMoments::FillMCTrackInfo()
{
  Int_t counterBin[mPtmax] = { 0 };
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

    Bool_t isElectron = kFALSE;
    if (lPart->PdgCode() == 11 || lPart->PdgCode() == -11)
      isElectron = kTRUE;

    if (flagRejEls)
      if (isElectron)
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
        for (Int_t k = 0; k < mMBins; ++k) {
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

void AliAnalysisTaskNFactorialMoments::GetPtBin(Double_t pt)
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
              Two Track Calculation of dphistar
________________________________________________________________________*/

Float_t AliAnalysisTaskNFactorialMoments::CalculateDPhiStar(
  Float_t phi1, Float_t eta1, Float_t pt1, Int_t charge1, Float_t phi2, Float_t eta2,
  Float_t pt2, Int_t charge2, Float_t bSign)
{
  Float_t radius = 0;
  Float_t dphistartemp = 0;
  Float_t dphistar = 999;
  Float_t kLimit = fdeta * 3;
  Double_t kPi = TMath::Pi();
  bSign = (bSign > 0) ? 1 : -1;

  // variables and cuts have been taken from
  // https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
  if (abs(eta1 - eta2) < fdeta * 2.5 * 3) {
    // check first boundaries to see if is worth to loop and find the minimum
    Float_t dphistar1 =
      GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
    Float_t dphistar2 =
      GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);

    if (fabs(dphistar1) < kLimit || fabs(dphistar2) < kLimit ||
        (dphistar1 * dphistar2) < 0) {
      // find the smallest dphistar (radius of TPC: 0.8 - 2.5 (m))
      for (Int_t rad(80); rad < 251; rad++) {
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
    }
  }
  return dphistar;
}

/*________________________________________________________________________
          Shared Cluster Fraction Calculation
________________________________________________________________________*/

Float_t AliAnalysisTaskNFactorialMoments::SharedClusterFraction(TBits& cl1,
                                                                TBits& cl2,
                                                                TBits& sh1,
                                                                TBits& sh2)
{
  Int_t ncl1 = cl1.GetNbits();
  Int_t ncl2 = cl2.GetNbits();
  Int_t sumCls = 0;
  Int_t sumSha = 0;
  Int_t sumQ = 0;
  Double_t shfrac = 0;
  for (Int_t ich = 0; ich < ncl1 && ich < ncl2; ich++) {
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

Float_t AliAnalysisTaskNFactorialMoments::GetDPhiStar(Float_t phi1, Float_t pt1,
                                                      Float_t charge1, Float_t phi2,
                                                      Float_t pt2, Float_t charge2,
                                                      Float_t radius, Float_t bSign)
{
  Double_t kPi = TMath::Pi();
  Double_t deltaPhi = phi2 - phi1;
  if (deltaPhi < -0.5 * TMath::Pi())
    deltaPhi += TMath::TwoPi();
  if (deltaPhi > 1.5 * TMath::Pi())
    deltaPhi -= TMath::TwoPi();
  // calculates dphistar
  Float_t dphistarm = deltaPhi +
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
          Calculation of the Normalized Fq MomentsAODs
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::CalculateNFMs(TH2D* h1[mPtmax][mMBins], Bool_t mcGen)
{
  for (Int_t iPt = 0; iPt < mPtBins; ++iPt) {
    for (Int_t iM = 0; iM < mMBins; iM++) {

      Double_t NoOfBins = 0;
      Double_t MSquare = 0;
      if (Mmax == 123)
        NoOfBins = 3 * (iM + 2);
      if (Mmax == 82)
        NoOfBins = 2 * (iM + 2);
      if (flagSelfAff)
        NoOfBins = mMBin2[iM] * mNBin2[iM];
      MSquare = TMath::Power(NoOfBins, mDim);
      if (flagSelfAff)
        MSquare = mMBin2[iM] * mNBin2[iM];
      Double_t SumOfbincontent = 0;
      Double_t FqEvent[mQs];
      Double_t sumoff[mQs];
      Double_t bincontent = 0.0;
      Double_t Mbin = NoOfBins;
      Int_t NofXetabins = 0;
      Int_t NofXphibins = 0;

      for (Int_t index = 0; index < mQs; index++) {
        FqEvent[index] = 0.0;
        sumoff[index] = 0.0;
      }

      NofXetabins = h1[iPt][iM]->GetNbinsX();
      NofXphibins = h1[iPt][iM]->GetNbinsY();

      for (Int_t etabin = 1; etabin <= NofXetabins; etabin++) {
        for (Int_t phibin = 1; phibin <= NofXphibins; phibin++) {
          bincontent = 0.0;
          bincontent = h1[iPt][iM]->GetBinContent(etabin, phibin);
          SumOfbincontent += bincontent;

          for (Int_t q = 0; q < mQs; q++) {
            if (bincontent >= (q + 2)) {
              Double_t Fqeofbin = 0.0;
              Fqeofbin = TMath::Factorial(bincontent) /
                         TMath::Factorial(bincontent - (q + 2));

              if (flagMC) {
                if (TMath::IsNaN(Fqeofbin)) {
                  break;
                }
              }
              sumoff[q] += Fqeofbin;
            }
          }
        }
      }

      Double_t Av_bincontent = SumOfbincontent / MSquare;
      for (Int_t q = 0; q < mQs; q++) {
        if (sumoff[q] > 0.0) {
          FqEvent[q] = sumoff[q] / (MSquare);
        }
      }

      Float_t Fq2e = FqEvent[0];
      Float_t Fq3e = FqEvent[1];
      Float_t Fq4e = FqEvent[2];
      Float_t Fq5e = FqEvent[3];
      Float_t Fq6e = FqEvent[4];
      Float_t Fq7e = FqEvent[5];

      if (!mcGen) {
        fntpMBin[iPt][iM]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e,
                                Fq6e, Fq7e);
      } else {
        fntpMBinGen[iPt][iM]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e,
                                   Fq5e, Fq6e, Fq7e);
      }
    } // end of MBin loop
  }   // end of calculation loop
}

/*________________________________________________________________________
          Get the DCA of the track
________________________________________________________________________*/

Bool_t AliAnalysisTaskNFactorialMoments::GetDCA(AliAODTrack* track, Double_t dca[2])
{
  /// Get track dca aod trck
  if (track->TestBit(AliAODTrack::kIsDCA)) {
    dca[0] = track->DCA();
    dca[1] = track->ZAtDCA();
    return kTRUE;
  }

  Bool_t ok = kFALSE;
  if (fAOD) {
    Double_t covdca[3];
    // AliAODTrack copy(*track);
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(track);

    Float_t xstart = etp.GetX();
    if (xstart > 3.) {
      dca[0] = -999.;
      dca[1] = -999.;
      // printf("This method can be used only for propagation inside the beam pipe \n");
      return kFALSE;
    }

    AliAODVertex* vtx = (AliAODVertex*)(fAOD->GetPrimaryVertex());
    Double_t fBzkG = fAOD->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx, fBzkG, kVeryBig, dca, covdca);
    // ok = copy.PropagateToDCA(vtx,fBzkG,kVeryBig,dca,covdca);
  }
  if (!ok) {
    dca[0] = -999.;
    dca[1] = -999.;
  }
  return ok;
}

/*________________________________________________________________________
          Post data to the output slot(s)
________________________________________________________________________*/

void AliAnalysisTaskNFactorialMoments::DataPosting()
{
  PostData(1, fHistList);
  for (Int_t i = 0; i < mPtBins; i++) {
    PostData(i + 2, fNtupleList[i]);
    if (flagMC)
      PostData(i + 2 + mPtBins, fNtupleListGen[i]);
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
