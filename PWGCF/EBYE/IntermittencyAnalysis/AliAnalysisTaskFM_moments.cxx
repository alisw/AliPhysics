// Implementation of the Intermittency analysis
// for charged particles, two dimensions (eta and phi)
// Data and MC
// Contributors: salman
// Date: 04-May-2025

#include "AliAnalysisTaskFM_moments.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "TMath.h"

class AliAnalysisTaskFM_moments;
using namespace std;

ClassImp(AliAnalysisTaskFM_moments)
    AliAnalysisTaskFM_moments::AliAnalysisTaskFM_moments()
    : AliAnalysisTaskSE(), fAOD(nullptr), fMCEvent(nullptr), _pool(nullptr), fEvPoolMgr(nullptr), fPIDResponse(nullptr), mcHeader(nullptr), fEventCuts(), counter(0), c_field(0), c_percentile(0.0), hEventCounter(nullptr), hTrackCounter(nullptr), hTrackCounterGen(nullptr), hEnvsCounter(nullptr), hQACent{nullptr}, hQAVx(nullptr), hQAVy(nullptr), hQAVz{nullptr}, hDCAxy{nullptr}, hDCAz{nullptr}, hDCAxy_imp{nullptr}, hDCAz_imp{nullptr}, hDCAxypT{nullptr}, hDCAzpT{nullptr}, hnITScls{nullptr}, hnITScls2{nullptr}, hnTPCcls{nullptr}, hnTPCcls2{nullptr}, hnTPCcrossedrows{nullptr}, hmissingCls{nullptr}, htpcTgl{nullptr}, hNShCls{nullptr}, hNShClsFra{nullptr}, hNShClsNcls{nullptr}, hNShClsNcrows{nullptr}, hNFoundClsFra{nullptr}, hNShClsVsPt{nullptr}, hNFcls{nullptr}, hNFindClsFra{nullptr}, hNFindClsVsPt{nullptr}, hTPCsignal{nullptr}, hTPCsignalN{nullptr}, hTPCsignalvsPt{nullptr}, hTPCsignalvsPtN{nullptr}, hTPCsignalvsPtTuned{nullptr}, hTPCsignalvsPtot{nullptr}, hTPCsignalTuned{nullptr}, hITSsignal{nullptr}, hChi2TPC{nullptr}, hChi2ITS{nullptr}, hChi2perclTPC{nullptr}, hPtDis{nullptr}, hPtotDis{nullptr}, hEtaDis{nullptr}, hPhiDis{nullptr}, hPtDisGen{nullptr}, hEtaDisGen{nullptr}, hPhiDisGen{nullptr}, hGoldenChi2{nullptr}, hNumPions{nullptr}, hNumPionsCent{nullptr}, hNumKaons{nullptr}, hNumKaonsCent{nullptr}, hNumProtons{nullptr}, hNumProtonsCent{nullptr}, hPtBin{nullptr}, hEtaBin{nullptr}, hPhiBin{nullptr}, hMultBin{nullptr}, hPtBinGen{nullptr}, hEtaBinGen{nullptr}, hPhiBinGen{nullptr}, hMultBinGen{nullptr}, hEtaPhiBin{nullptr}, hEtaPhiBinGen{nullptr}, fntpMBin{nullptr}, fntpMBinGen{nullptr}, hdEta(nullptr), hdPhi(nullptr), hDEDPraw{nullptr}, hDEDPSel{nullptr}, hDEDPrawChSame{nullptr}, hDEDPrawChDiff{nullptr}, hDEDPrawPtOrd{nullptr}, hList{nullptr}, fNtupleList{nullptr}, fNtupleListGen{nullptr} {}
AliAnalysisTaskFM_moments::AliAnalysisTaskFM_moments(const char *name)
    : AliAnalysisTaskSE(name), fAOD(nullptr), fMCEvent(nullptr), _pool(nullptr), fEvPoolMgr(nullptr), fPIDResponse(nullptr), mcHeader(nullptr), fEventCuts(), counter(0), c_field(0), c_percentile(0.0), hEventCounter(nullptr), hTrackCounter(nullptr), hTrackCounterGen(nullptr), hEnvsCounter(nullptr), hQACent{nullptr}, hQAVx(nullptr), hQAVy(nullptr), hQAVz{nullptr}, hDCAxy{nullptr}, hDCAz{nullptr}, hDCAxy_imp{nullptr}, hDCAz_imp{nullptr}, hDCAxypT{nullptr}, hDCAzpT{nullptr}, hnITScls{nullptr}, hnITScls2{nullptr}, hnTPCcls{nullptr}, hnTPCcls2{nullptr}, hnTPCcrossedrows{nullptr}, hmissingCls{nullptr}, htpcTgl{nullptr}, hNShCls{nullptr}, hNShClsFra{nullptr}, hNShClsNcls{nullptr}, hNShClsNcrows{nullptr}, hNFoundClsFra{nullptr}, hNShClsVsPt{nullptr}, hNFcls{nullptr}, hNFindClsFra{nullptr}, hNFindClsVsPt{nullptr}, hTPCsignal{nullptr}, hTPCsignalN{nullptr}, hTPCsignalvsPt{nullptr}, hTPCsignalvsPtN{nullptr}, hTPCsignalvsPtTuned{nullptr}, hTPCsignalvsPtot{nullptr}, hTPCsignalTuned{nullptr}, hITSsignal{nullptr}, hChi2TPC{nullptr}, hChi2ITS{nullptr}, hChi2perclTPC{nullptr}, hPtDis{nullptr}, hPtotDis{nullptr}, hEtaDis{nullptr}, hPhiDis{nullptr}, hPtDisGen{nullptr}, hEtaDisGen{nullptr}, hPhiDisGen{nullptr}, hGoldenChi2{nullptr}, hNumPions{nullptr}, hNumPionsCent{nullptr}, hNumKaons{nullptr}, hNumKaonsCent{nullptr}, hNumProtons{nullptr}, hNumProtonsCent{nullptr}, hPtBin{nullptr}, hEtaBin{nullptr}, hPhiBin{nullptr}, hMultBin{nullptr}, hPtBinGen{nullptr}, hEtaBinGen{nullptr}, hPhiBinGen{nullptr}, hMultBinGen{nullptr}, hEtaPhiBin{nullptr}, hEtaPhiBinGen{nullptr}, fntpMBin{nullptr}, fntpMBinGen{nullptr}, hdEta(nullptr), hdPhi(nullptr), hDEDPraw{nullptr}, hDEDPSel{nullptr}, hDEDPrawChSame{nullptr}, hDEDPrawChDiff{nullptr}, hDEDPrawPtOrd{nullptr}, hList{nullptr}, fNtupleList{nullptr}, fNtupleListGen{nullptr} {

  Info("AliAnalysisTaskFM_moments", "Specific Constructor");

  DefineInput(0, TChain::Class());
  const Int_t _numOutputs = _FLAG_MC_ ? (_maxPtBin * _numEnvs + _maxPtBin) : (_maxPtBin * _numEnvs);
  DefineOutput(1, TList::Class());
  for (Int_t i = 0; i < _numOutputs; ++i) {
    DefineOutput(i + 2, TList::Class());
  }
  std::fill_n(&foundPart[0][0], _maxPtBin * _numEnvs, 0);
}

AliAnalysisTaskFM_moments::~AliAnalysisTaskFM_moments() {
  if (fPIDResponse) {
    delete fPIDResponse;
    fPIDResponse = 0x0;
  }
}

void AliAnalysisTaskFM_moments::UserCreateOutputObjects() {
  for (Int_t i = 0; i < _numEnvs; ++i) {
    hList[i] = new TList();
    hList[i]->SetOwner(kTRUE);
  }
  for (Int_t i = 0; i < _numPtBin; ++i) {
    for (Int_t j = 0; j < _numEnvs; ++j) {
    fNtupleList[i][j] = new TList();
    fNtupleList[i][j]->SetOwner(kTRUE);
    }
    if (_FLAG_MC_) {
      fNtupleListGen[i] = new TList();
      fNtupleListGen[i]->SetOwner(kTRUE);
    }
  }
  TString dataStr = _FLAG_MC_ ? "rec" : "";
  TString conStr[_numEnvs] = { "fb128", "fb768", "fb768&&rchi2>80", "fb768&&shclsncls<0.3&&shclsncrows<0.25&&nfindableclsncls>0.8", "fb768&&rchi2>80&&shclsncls<0.3&&shclsncrows<0.25&&nfindableclsncls>0.8", "fb768&&pTDCAxy&&rchi2>80&&shclsncls<0.3&&shclsncrows<0.25&&nfindableclsncls>0.8" };
  
  hEventCounter = new TH1F("hEventCounter", ";;no. of events", 9, 0.5, 9.5);
  hEventCounter->SetFillColor(kGreen);
  hList[0]->Add(hEventCounter);
  hTrackCounter = new TH1F("hTrackCounter", ";;no. of tracks", 17, 0.5, 17.5);
  hTrackCounterGen = new TH1F("hTrackCounterGen", ";;no. of tracks", 8, 0.5, 8.5);
  hTrackCounter->SetFillColor(kGreen);
  hTrackCounterGen->SetFillColor(kGreen);
  hList[0]->Add(hTrackCounter);
  if (_FLAG_MC_) hList[0]->Add(hTrackCounterGen);
  hEnvsCounter = new TH1F("hEnvsCounter", ";;no. of tracks", _numEnvs+1, 0.5, _numEnvs + 1.5);
  for (Int_t i = 0; i < _numEnvs; ++i) {
    hEnvsCounter->GetXaxis()->SetBinLabel(i + 1, conStr[i].Data());
  }
  hEnvsCounter->SetFillColor(kGreen);
  hList[0]->Add(hEnvsCounter);

  hQACent[0] = new TH1F("hQACent_raw", "Centrality Distribution", 100, 0, 100);
  hQACent[1] = new TH1F("hQACent_sel", "Centrality Distribution", 100, 0, 100);
  hQAVx = new TH1F("hQAVx", "Primary vertex distribution #minus x coordinate;V_{x}", 1500, -0.3, 0.3);
  hQAVy = new TH1F("hQAVy", "Primary vertex distribution #minus y coordinate;V_{y}", 1500, -0.3, 0.3);
  hQAVz[0] = new TH1F("hQAVz_raw", "Primary vertex distribution #minus z coordinate;V_{z}", 1500, -20.0, 20.0);
  hQAVz[1] = new TH1F("hQAVz_sel", "Primary vertex distribution #minus z coordinate;V_{z}", 500, -20.0, 20.0);
  hList[0]->Add(hQACent[0]);
  hList[0]->Add(hQACent[1]);
  hList[0]->Add(hQAVx);
  hList[0]->Add(hQAVy);
  hList[0]->Add(hQAVz[0]);
  hList[0]->Add(hQAVz[1]);

  TString binCStr[_maxPtBin][_numEnvs];
  TString binStr[_maxPtBin];
  for (Int_t i = 0; i < _numPtBin; ++i) {
    binStr[i] = Form("bin%.2fto%.2f", _ptArray[2 * i], _ptArray[2 * i + 1]);
    for (Int_t j = 0; j < _numEnvs; ++j) binCStr[i][j] = Form("bin%.2fto%.2f_%s", _ptArray[2 * i], _ptArray[2 * i + 1], conStr[j].Data());
    }

  for (Int_t i = 0; i < _numEnvs; ++i) {
    hDCAxy[i] = new TH1F(Form("hDCAxy%s", conStr[i].Data()), "DCAxy distribution;DCAxy;Counts", 4000, -10.5, 10.5);
    hDCAz[i] = new TH1F(Form("hDCAz%s", conStr[i].Data()), "DCAz distribution;DCAz;Counts", 4000, -10.5, 10.5);
    hDCAxy_imp[i] = new TH1F(Form("hDCAxy_imp%s", conStr[i].Data()), "DCAxy calculated with impact parameter;DCAxy;Counts", 4000, -10.5, 10.5);
    hDCAz_imp[i] = new TH1F(Form("hDCAz_imp%s", conStr[i].Data()), "DCAz calculated with impact parameter;DCAz;Counts", 4000, -10.5, 10.5);
    hDCAxypT[i] = new TH2F(Form("hDCAxypT%s", conStr[i].Data()), "DCAxy vs pT;p_{T};dca_{xy}", 2500, 0, 5.5, 2500, -4.0, 4.0);
    hDCAzpT[i] = new TH2F(Form("hDCAzpT%s", conStr[i].Data()), "DCAz vs pT;p_{T};dca_{z}", 2500, 0, 5.5, 2500, -4.0, 4.0);
    hnITScls[i] = new TH1F(Form("hnITScls%s", conStr[i].Data()), "ITS cluster distribution;ITS cluster;Counts", 300, -0.5, 10.5);
    hnITScls2[i] = new TH1F(Form("hnITScls2%s", conStr[i].Data()), "ITS cluster distribution;ITS cluster;Counts", 1e4, 10e3, 80e3);
    hnTPCcls[i] = new TH1F(Form("hnTPCcls%s", conStr[i].Data()), "TPC cluster distribution;TPC cluster;Counts", 300, 0, 200);
    hnTPCcls2[i] = new TH1F(Form("hnTPCcls2%s", conStr[i].Data()), "TPC cluster distribution;TPC cluster;Counts", 1e6, 0, 10e6);
    hnTPCcrossedrows[i] = new TH1F(Form("hnTPCcrossedrows%s", conStr[i].Data()), "TPC crossed rows distribution;TPC crossed rows;Counts", 500, 0, 200);
    hmissingCls[i] = new TH1F(Form("hmissingCls%s", conStr[i].Data()), "TPC missing cluster distribution;missingcls;Counts", 500, 0, 10);
    htpcTgl[i] = new TH1F(Form("htpcTgl%s", conStr[i].Data()), "TPC tgl distribution;TPC tgl;Counts", 1000, -2.0, 2.0);
    hNShCls[i] = new TH1F(Form("hNShCls%s", conStr[i].Data()), "TPC shared cluster distribution", 600, 0, 500);
    hNShClsFra[i] = new TH2F(Form("hNShClsFra%s", conStr[i].Data()), "TPC shared cluster fraction;sharedcls/ncls;sharedcls/ncrrows", 2000, 0, 1.4, 2000, 0, 1.4);
    hNShClsNcls[i] = new TH1F(Form("hNShClsNCls%s", conStr[i].Data()), "TPC shared cluster fraction;sharedcls/ncls", 1000, 0, 1.0);
    hNShClsNcrows[i] = new TH1F(Form("hNShClsNRows%s", conStr[i].Data()), "TPC shared cluster fraction;sharedcls/ncrrows", 1000, 0, 1.0);
    hNFoundClsFra[i] = new TH2F(Form("hNFoundClsFra%s", conStr[i].Data()), "TPC found cluster fraction;sharedcls/ncls;ncrrows/findablecls", 2000, 0, 2, 2000, 0, 10);
    hNShClsVsPt[i] = new TH2F(Form("hNShClsVsPt%s", conStr[i].Data()), "TPC shared cluster fraction;pt;sharedcls/ncls", 1500, 0, 3.0, 1000, 0, 1.2);
    hNFcls[i] = new TH1F(Form("hNFcls%s", conStr[i].Data()), "TPC findable cluster distribution;findablecls;counts", 400, 0, 300);
    hNFindClsFra[i] = new TH1F(Form("hNFoundClsFra%s", conStr[i].Data()), "TPC findable cluster fraction;findablecls/ncls;counts", 400, 0, 60);
    hNFindClsVsPt[i] = new TH2F(Form("hNFoundClsVsPt%s", conStr[i].Data()), "TPC findable cluster fraction;pt;findablecls/ncls", 1500, 0, 3.0, 2000, 0, 10.0);
    hTPCsignal[i] = new TH1F(Form("hTPCsignal%s", conStr[i].Data()), "TPC signal distribution;TPC signal;Counts", 3000, 0, 3500);
    hTPCsignalN[i] = new TH1F(Form("hTPCsignalN%s", conStr[i].Data()), "TPC signal distribution;TPC signal;Counts", 3000, 0, 3500);
    hTPCsignalvsPt[i] = new TH2F(Form("hTPCsignalvsPt%s", conStr[i].Data()), "TPC signal vs pT;pt;TPC signal", 3000, 0, 5.5, 3000, 0, 1500);
    hTPCsignalvsPtN[i] = new TH2F(Form("hTPCsignalvsPtN%s", conStr[i].Data()), "TPC signal vs pT;pt;TPC signal", 3000, 0, 5.5, 3000, 0, 1500);
    hTPCsignalvsPtTuned[i] = new TH2F(Form("hTPCsignalvsPtTuned%s", conStr[i].Data()), "TPC signal vs pT;pt;TPC signal", 3000, 0, 5.5, 3000, 0, 1500);
    hTPCsignalvsPtot[i] = new TH2F(Form("hTPCsignalvsPtot%s", conStr[i].Data()), "TPC signal vs Ptot;Ptot;TPC signal", 3000, 0, 10, 3000, 0, 1500);
    hTPCsignalTuned[i] = new TH1F(Form("hTPCsignalTuned%s", conStr[i].Data()), "TPC signal distribution;TPC signal;Counts", 3000, 0, 3500);
    hITSsignal[i] = new TH1F(Form("hITSsignal%s", conStr[i].Data()), "ITS signal distribution;ITS signal;Counts", 10000, 0, 15000);
    hChi2TPC[i] = new TH1F(Form("hChi2TPC%s", conStr[i].Data()), "TPC chi2 distribution;TPC chi2;Counts", 1000, 0, 2000);
    hChi2ITS[i] = new TH1F(Form("hChi2ITS%s", conStr[i].Data()), "ITS chi2 distribution;ITS chi2;Counts", 1000, 0, 2000);
    hChi2perclTPC[i] = new TH1F(Form("hChi2perclTPC%s", conStr[i].Data()), "TPC chi2 per cluster distribution;TPC chi2 per cluster;Counts", 1000, 0, 100);
    hPtDis[i] = new TH1F(Form("hPtDis%s", conStr[i].Data()), "Pt distribution;Pt;Counts", 1000, 0, 10);
    hPtotDis[i] = new TH1F(Form("hPtotDis%s", conStr[i].Data()), "Ptot distribution;Ptot;Counts", 1000, 0, 10);
    hEtaDis[i] = new TH1F(Form("hEtaDis%s", conStr[i].Data()), "Eta distribution;Eta;Counts", 1000, -1.0, 1.0);
    hPhiDis[i] = new TH1F(Form("hPhiDis%s", conStr[i].Data()), "Phi distribution;Phi;Counts", 1000, 0, 6.5);
    hGoldenChi2[i] = new TH1F(Form("hGoldenChi2%s", conStr[i].Data()), "Golden Chi2 distribution;Golden Chi2;Counts", 1000, 0, 2000);
    hNumPions[i] = new TH1F(Form("hNumPions%s", conStr[i].Data()), "Number of Pions;num of Pions;Counts", 10000, 0, 10000);
    hNumPionsCent[i] = new TH2F(Form("hNumPionsCent%s", conStr[i].Data()), "Number of Pions vs Centrality;cent(%);num of Pions", 500, _minCent, _maxCent, 10000, 0, 10000);
    hNumKaons[i] = new TH1F(Form("hNumKaons%s", conStr[i].Data()), "Number of Kaons;num of Kaons;Counts", 1000, 0, 2000);
    hNumKaonsCent[i] = new TH2F(Form("hNumKaonsCent%s", conStr[i].Data()), "Number of Kaons vs Centrality;cent(%);num of Kaons", 500, _minCent, _maxCent, 1000, 0, 2000);
    hNumProtons[i] = new TH1F(Form("hNumProtons%s", conStr[i].Data()), "Number of Protons;num of Protons;Counts", 1000, 0, 2000);
    hNumProtonsCent[i] = new TH2F(Form("hNumProtonsCent%s", conStr[i].Data()), "Number of Protons vs Centrality;cent(%);num of Protons", 500, _minCent, _maxCent, 1000, 0, 2000);
    hList[i]->Add(hDCAxy[i]);
    hList[i]->Add(hDCAz[i]);
    hList[i]->Add(hDCAxy_imp[i]);
    hList[i]->Add(hDCAz_imp[i]);
    hList[i]->Add(hDCAxypT[i]);
    hList[i]->Add(hDCAzpT[i]);
    hList[i]->Add(hnITScls[i]);
    hList[i]->Add(hnITScls2[i]);
    hList[i]->Add(hnTPCcls[i]);
    hList[i]->Add(hnTPCcls2[i]);
    hList[i]->Add(hnTPCcrossedrows[i]);
    hList[i]->Add(hmissingCls[i]);
    hList[i]->Add(htpcTgl[i]);
    hList[i]->Add(hNShCls[i]);
    hList[i]->Add(hNShClsFra[i]);
    hList[i]->Add(hNShClsNcls[i]);
    hList[i]->Add(hNShClsNcrows[i]);
    hList[i]->Add(hNFoundClsFra[i]);
    hList[i]->Add(hNShClsVsPt[i]);
    hList[i]->Add(hNFcls[i]);
    hList[i]->Add(hNFindClsFra[i]);
    hList[i]->Add(hNFindClsVsPt[i]);
    hList[i]->Add(hTPCsignal[i]);
    hList[i]->Add(hTPCsignalN[i]);
    hList[i]->Add(hTPCsignalvsPt[i]);
    hList[i]->Add(hTPCsignalvsPtN[i]);
    hList[i]->Add(hTPCsignalvsPtTuned[i]);
    hList[i]->Add(hTPCsignalvsPtot[i]);
    hList[i]->Add(hTPCsignalTuned[i]);
    hList[i]->Add(hITSsignal[i]);
    hList[i]->Add(hChi2TPC[i]);
    hList[i]->Add(hChi2ITS[i]);
    hList[i]->Add(hChi2perclTPC[i]);
    hList[i]->Add(hPtDis[i]);
    hList[i]->Add(hPtotDis[i]);
    hList[i]->Add(hEtaDis[i]);
    hList[i]->Add(hPhiDis[i]);
    hList[i]->Add(hGoldenChi2[i]);
    hList[i]->Add(hNumPions[i]);
    hList[i]->Add(hNumPionsCent[i]);
    hList[i]->Add(hNumKaons[i]);
    hList[i]->Add(hNumKaonsCent[i]);
    hList[i]->Add(hNumProtons[i]);
    hList[i]->Add(hNumProtonsCent[i]);

    if (_FLAG_MC_ && i == 0) {
      hPtDisGen = new TH1F("hPtDisGen", "Pt distribution of Gen tracks;#it{p}_{T};Counts", 1000, 0.0, 5.0);
      hEtaDisGen = new TH1F("hEtaDisGen", "Eta distribution of Gen tracks;#eta;Counts", 1000, -1.0, 1.0);hPhiDisGen = new TH1F("hPhiDisGen", "Phi distribution of Gen tracks;#phi;Counts", 1000, 0.0, 6.5);
      hList[i]->Add(hPtDisGen);
      hList[i]->Add(hEtaDisGen);
      hList[i]->Add(hPhiDisGen);
    }
    
    for (Int_t j = 0; j < _numPtBin; ++j) {
      hPtBin[j][i] = new TH1F(Form("hPtBin%s%s", dataStr.Data(), binCStr[j][i].Data()), Form("pt of %s;pt;counts", binCStr[j][i].Data()), 1000, 0.0, 5.0);
      hEtaBin[j][i] = new TH1F(Form("hEtaBin%s%s", dataStr.Data(), binCStr[j][i].Data()), Form("eta of %s;eta;counts", binCStr[j][i].Data()), 1000, -1.0, 1.0);
      hPhiBin[j][i] = new TH1F(Form("hPhiBin%s%s", dataStr.Data(), binCStr[j][i].Data()), Form("phi of %s;phi;counts", binCStr[j][i].Data()), 1000, 0.0, 6.5);
      hMultBin[j][i] = new TH1F(Form("hMultBin%s%s", dataStr.Data(), binCStr[j][i].Data()), Form("multiplicity of %s;Multiplicity;counts", binCStr[j][i].Data()), 10000, 0.0, 10000.0);
      hList[i]->Add(hPtBin[j][i]);
      hList[i]->Add(hEtaBin[j][i]);
      hList[i]->Add(hPhiBin[j][i]);
      hList[i]->Add(hMultBin[j][i]);
      // for MC generated tracks
      if (_FLAG_MC_ && i == 0) {
        hPtBinGen[j] = new TH1F(Form("hPtBinGen%s", binStr[j].Data()), Form("pt of %s for gen;pt;counts", binStr[j].Data()), 1000, 0.0, 5.0);
        hEtaBinGen[j] = new TH1F(Form("hEtaBinGen%s", binStr[j].Data()), Form("eta of %s for gen;eta;counts", binStr[j].Data()), 1000, -1.0, 1.0);
        hPhiBinGen[j] = new TH1F(Form("hPhiBinGen%s", binStr[j].Data()), Form("phi of %s for gen;phi;counts", binStr[j].Data()), 1000, 0.0, 6.5);
        hMultBinGen[j] = new TH1F(Form("hMultBinGen%s", binStr[j].Data()), Form("multiplicity of %s for gen;Multiplicity;counts", binStr[j].Data()), 10000, 0.0, 10000.0);
        hList[i]->Add(hPtBinGen[j]);
        hList[i]->Add(hEtaBinGen[j]);
        hList[i]->Add(hPhiBinGen[j]);
        hList[i]->Add(hMultBinGen[j]);
      } // for MC generated tracks
      // for 2D eta-phi distributions over m phase space bins
      for (Int_t k = 0; k < _numM; ++k) {
        Int_t _tmp = _mBin2[k]*_mBin2[k];
        hEtaPhiBin[j][k][i] = new TH2F(Form("hEtaPhiBin%s%s_M%d", dataStr.Data(), binCStr[j][i].Data(), _tmp), Form("eta#minusphi of %s #minus M^{2} %d;#eta;#phi", binCStr[j][i].Data(), _tmp), _mBin2[k], _minEta, _maxEta, _mBin2[k], 0.0, TMath::TwoPi());
        hList[i]->Add(hEtaPhiBin[j][k][i]);
        if (_FLAG_MC_ && i == 0) {
          hEtaPhiBinGen[j][k][0] = new TH2F(Form("hEtaPhiBinGen%s_M%d", binStr[j].Data(), _tmp), Form("eta#minusphi of %s #minus M^{2} %d;#eta;#phi", binStr[j].Data(), _tmp), _mBin2[k], _minEta, _maxEta, _mBin2[k], 0.0, TMath::TwoPi());
          hList[i]->Add(hEtaPhiBinGen[j][k][0]);
        }
        // for tuples to store values of mult, M2, average bin content and factorial moments
        fntpMBin[j][k][i] = new TNtuple(Form("fntpMBin%s%s_M%d", dataStr.Data(), binCStr[j][i].Data(), k), Form("fntpMBin%s%s_M%d", dataStr.Data(), binCStr[j][i].Data(), k), "Mult:Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
        fNtupleList[j][i]->Add(fntpMBin[j][k][i]);
        if (_FLAG_MC_ && i == 0) {
          fntpMBinGen[j][k] = new TNtuple(Form("fntpMBingen%s_M%d", binStr[j].Data(), k), Form("fntpMBingen%s%s_M%d", dataStr.Data(), binStr[j].Data(), k), "Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
          fNtupleListGen[i]->Add(fntpMBinGen[j][k]);
        }
      } // loop over num of phase space bins
    }// loop over num of pt bins
  } // loop over _numEnvs


  // histograms for detadphi and PID
  TString chStr[3] = {"Pos", "Neg", "All"};
  TString rangeStr[2] = {"", "Close"};
  TString ordStr[4] = {"pt2>pt1_dphi<0", "pt2>pt1_dphi>0", "pt2<pt1_dphi<0", "pt2<pt1_dphi>0"};

  for (Int_t i = 0; i < _numEnvs; ++i) {
    hdEta[i] = new TH1F(Form("hdEta%s", conStr[i].Data()), "dEta;#Delta#eta;Counts", 100, -1.75, 1.75);
    hdPhi[i] = new TH1F(Form("hdPhi%s", conStr[i].Data()), "dPhi;#Delta#phi;Counts", 100, -TMath::Pi(), TMath::Pi());
    hList[i]->Add(hdEta[i]);
    hList[i]->Add(hdPhi[i]);
    for (Int_t j = 0; j < _numPtBin; ++j) {
      hDEDPraw[j][i] = new TH2F(Form("hDEDPraw_%s", binCStr[j][i].Data()), Form("dEta dPhi for %s;#Delta#eta;#Delta#phi", binCStr[j][i].Data()), 100, -1.75, 1.75, 100, -TMath::Pi(), TMath::Pi());
      hDEDPSel[j][i] = new TH2F(Form("hDEDPSel_%s", binCStr[j][i].Data()), Form("dEta dPhi for %s;#Delta#eta;#Delta#phi", binCStr[j][i].Data()), 100, -1.75, 1.75, 100, -TMath::Pi(), TMath::Pi());
      hList[i]->Add(hDEDPraw[j][i]);
      hList[i]->Add(hDEDPSel[j][i]);
      // for all charges
      for (Int_t ch = 0; ch < 3; ++ch) {
        hDEDPrawChSame[j][0][ch][i] = new TH2F(Form("hDEDPChSame_%s_%s%s", binCStr[j][i].Data(), chStr[ch].Data(), rangeStr[0].Data()), Form("dEta dPhi for %s %s;#Delta#eta;#Delta#phi", binCStr[j][i].Data(), chStr[ch].Data()), 100, -1.75, 1.75, 100, -TMath::Pi(), TMath::Pi());
        hDEDPrawChSame[j][1][ch][i] = new TH2F(Form("hDEDPChSame_%s_%s%s", binCStr[j][i].Data(), chStr[ch].Data(), rangeStr[1].Data()), Form("dEta dPhi for %s %s;#Delta#eta;#Delta#phi", binCStr[j][i].Data(), chStr[ch].Data()), 100, -0.02, 0.02, 100, -0.05, 0.05);
        hList[i]->Add(hDEDPrawChSame[j][0][ch][i]);
        hList[i]->Add(hDEDPrawChSame[j][1][ch][i]);
        // for pt ordered
        for (Int_t k = 0; k < 4; ++k) {
          hDEDPrawPtOrd[j][k][ch][0][i] = new TH2F(Form("hDEDPChSame_%s_%s_%s%s", binCStr[j][i].Data(), ordStr[k].Data(), chStr[ch].Data(), rangeStr[0].Data()), Form("dEta dPhi for %s %s %s;#Delta#eta;#Delta#phi", binCStr[j][i].Data(), ordStr[k].Data(), chStr[ch].Data()), 100, -1.75, 1.75, 100, -TMath::Pi(), TMath::Pi());
          hDEDPrawPtOrd[j][k][ch][1][i] = new TH2F(Form("hDEDPChSame_%s_%s_%s%s", binCStr[j][i].Data(), ordStr[k].Data(), chStr[ch].Data(), rangeStr[1].Data()), Form("dEta dPhi for %s %s %s;#Delta#eta;#Delta#phi", binCStr[j][i].Data(), ordStr[k].Data(), chStr[ch].Data()), 100, -0.02, 0.02, 100, -0.05, 0.05);
          hList[i]->Add(hDEDPrawPtOrd[j][k][ch][0][i]);
          hList[i]->Add(hDEDPrawPtOrd[j][k][ch][1][i]);
        } // loop over pt ordered
      } // loop over charges
      // for different charges
      hDEDPrawChDiff[j][0][i] = new TH2F(Form("hDEDPChDiff_%s_%s", binCStr[j][i].Data(), rangeStr[0].Data()), Form("dEta dPhi for %s different ch;#Delta#eta;#Delta#phi", binCStr[j][i].Data()), 100, -1.75, 1.75, 100, -TMath::Pi(), TMath::Pi());
      hDEDPrawChDiff[j][1][i] = new TH2F(Form("hDEDPChDiff_%s_%s", binCStr[j][i].Data(), rangeStr[1].Data()), Form("dEta dPhi for %s different ch;#Delta#eta;#Delta#phi", binCStr[j][i].Data()), 100, -0.02, 0.02, 100, -0.05, 0.05);
      hList[i]->Add(hDEDPrawChDiff[j][0][i]);
      hList[i]->Add(hDEDPrawChDiff[j][1][i]);
    } // loop over pt bins 
  } // loop over _numEnvs 
  
  fEventCuts.AddQAplotsToList(hList[0], kTRUE);
  // initialize mixed event pool manager
  const Int_t _numVtxBins = 10;
  Double_t _binsVertex[_numVtxBins+1] = { -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10 };
  Double_t* _binsVertexP = _binsVertex;
  const Int_t _numBinsCent = 5;
  Double_t _binsCent[_numBinsCent+1] = { 0, 1, 2, 3, 4, 5 };
  Double_t* _binsCentP = _binsCent;
  TAxis* _axis = new TAxis(_numBinsCent, _binsCentP);
  fEvPoolMgr = new AliEventPoolManager(1000, 40000, _numBinsCent, (Double_t*)_axis->GetXbins()->GetArray(), _numVtxBins, _binsVertexP);
  std::cout << "\033[31m" << " == :: == task marooz == :: == " << " Event Pool Manager initialised: " << "\n pool depth: " << 1000 << "\n min tracks: " << 40000 << "\n vtx bins: " << _numVtxBins <<  "\n cent bins: " << _numBinsCent << " \033[0m" << std::endl;

  DataPosting();
}

void AliAnalysisTaskFM_moments::UserExec(Option_t *) {
  hEventCounter->Fill(1);
  hEventCounter->GetXaxis()->SetBinLabel(1, "All events");
  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

  fPIDResponse = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->GetPIDResponse();

  if (!eventHandler || !fAOD || !fPIDResponse)  std::cout << "\033[31m" <<  " UserExec: " << "EventHandlers not found" << " \033[0m" << std::endl;

  if (_FLAG_MC_) {
    fMCEvent = MCEvent();
    mcHeader = (AliAODMCHeader *)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!fMCEvent || !mcHeader) {
      std::cout << " == :: == task marooz == :: == " << " UserExec: " << "MCEvent or MCHeader not found" << std::endl;
      return;
    }
  }

  hEventCounter->Fill(2);
  hEventCounter->GetXaxis()->SetBinLabel(2, "Event Handler");

  bool pileupGen = kFALSE, acceptedEv = kTRUE;
  if (_FLAG_PILEUP_) {
  if (_FLAG_MC_) bool pileupGen = AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader, fGenName);
  else { 
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE); 
    acceptedEv = fEventCuts.AcceptEvent(fAOD); }
  } 
  if (!acceptedEv || pileupGen) return;

  hEventCounter->Fill(3);
  hEventCounter->GetXaxis()->SetBinLabel(3, "Pileup rejection");

/*    if (!_FLAG_MC_) {
      UInt_t fSelectMask = fInputHandler->IsEventSelected();
      if (_YEAR_ == 2015 || _YEAR_ == 2018) {
        if (!(fSelectMask & AliVEvent::kINT7))
          return;
      } else if (_YEAR_ == 2010) {
        if (!(fSelectMask & AliVEvent::kMB))
          return;
      }
    } */

  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  if (!fSelectMask) return;

  hEventCounter->Fill(4);
  hEventCounter->GetXaxis()->SetBinLabel(4, "Event Selection");

  AliAODVertex *vertex = fAOD->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1) return;
  Float_t lvx = (Float_t)vertex->GetX();
  Float_t lvy = (Float_t)vertex->GetY();
  Float_t lvz = (Float_t)vertex->GetZ();
  hQAVz[0]->Fill(lvz);
  if ((fabs(lvx) < 0.3) && (fabs(lvy) < 0.4) && (fabs(lvz) < 10.0)) {
    hQAVx->Fill(lvx);
    hQAVy->Fill(lvy);
    hQAVz[1]->Fill(lvz); _vzMax = lvz;
  } else return;

  hEventCounter->Fill(5);
  hEventCounter->GetXaxis()->SetBinLabel(5, "Vertex cut");

  // cout << "Vertex position: " << lvx << "\t" << lvy << "\t" << lvz << endl;
  //  Centrality Selection
  c_percentile = -1;
  if (_YEAR_ == 2015 || _YEAR_ == 2018) {
    AliMultSelection* MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    c_percentile = MultSelection->GetMultiplicityPercentile("V0M");
    }
  if (_YEAR_ == 2010) {
    AliCentrality* centrality = fAOD->GetCentrality();
    c_percentile = centrality->GetCentralityPercentile("V0M");
    }
 
  hQACent[0]->Fill(c_percentile);
  if (c_percentile < _minCent || c_percentile > _maxCent) return;
  std::cout << "\033[31m" <<  " UserExec: " << " Event passed with centrality " << c_percentile << " \033[0m" << std::endl;

  hEventCounter->Fill(6);
  hEventCounter->GetXaxis()->SetBinLabel(6, "Centrality cut");
  hQACent[1]->Fill(c_percentile);

  c_field = fAOD->GetMagneticField();
  if (fBfield == -1 && c_field > 0) return;
  if (fBfield == 1 && c_field < 0) return;

  hEventCounter->Fill(7);
  hEventCounter->GetXaxis()->SetBinLabel(7, "Magnetic field cut");

  // Reset eta-phi histograms at the beginning of each event
  for (Int_t i = 0; i < _numEnvs; ++i) {
    for (Int_t j = 0; j < _numPtBin; ++j) {
      for (Int_t k = 0; k < _numM; ++k) {
        if (hEtaPhiBin[j][k][i]) hEtaPhiBin[j][k][i]->Reset();
        if (_FLAG_MC_ && i == 0) {
          if (hEtaPhiBinGen[j][k][0]) hEtaPhiBinGen[j][k][0]->Reset();
        }
      }
    }
  }

  hEventCounter->Fill(8);
  hEventCounter->GetXaxis()->SetBinLabel(8, "Final statistics");

  //SetEventMixing(_FLAG_MIXED_);

  FillTrackInfo(); CalculateNFMs(kFALSE);
  if (_FLAG_MC_) {
    FillMCTrackInfo(); CalculateNFMs(kTRUE);
  }

  DataPosting();
}

void AliAnalysisTaskFM_moments::FillTrackInfo() {
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " FillTrackInfo: " << " Filling track information" << std::endl;
  Int_t foundPion[_numEnvs] = {0};
  Int_t foundKaon[_numEnvs] = {0};
  Int_t foundProton[_numEnvs] = {0};
  std::fill_n(&foundPart[0][0], _maxPtBin * _numEnvs, 0);

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  AliAODMCParticle *particle;
  // Track Loop:
  for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
    AliAODTrack *track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
    hTrackCounter->Fill(1);
    hTrackCounter->GetXaxis()->SetBinLabel(1, "All tracks");
    if (!track) continue;

    Int_t label = track->GetLabel();
    if (_FLAG_MC_ && TMath::Abs(label) < AODMCTrackArray->GetEntriesFast())  particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(label)));

    Int_t trackNature = -1;
    if (_FLAG_MC_ && particle->IsPhysicalPrimary()) trackNature = 0;
    if (_FLAG_MC_ && particle->IsSecondaryFromMaterial()) trackNature = 1;
    if (_FLAG_MC_ && particle->IsSecondaryFromWeakDecay()) trackNature = 2;
    if (_FLAG_MC_ && trackNature < 0) continue;
    hTrackCounter->Fill(2);
    hTrackCounter->GetXaxis()->SetBinLabel(2, "track nature");

    const Bool_t tuneOnDataTPC = fPIDResponse->IsTunedOnData() && ((fPIDResponse->GetTunedOnDataMask() & AliPIDResponse::kDetTPC) == AliPIDResponse::kDetTPC);
    Double_t dEdxTPC = tuneOnDataTPC ? fPIDResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();

    Float_t pt = track->Pt();
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Int_t charge = track->Charge();
    Double_t dcaXY, dcaZ;
    Double_t dca[2];
    if (!(GetDCA(track, dca))) continue;
    dcaXY = dca[0];
    dcaZ = dca[1];
    Float_t dcaXY_imp, dcaZ_imp;
    track->GetImpactParameters(dcaXY_imp, dcaZ_imp);
    hTrackCounter->Fill(3);
    hTrackCounter->GetXaxis()->SetBinLabel(3, "dca calc");

    if (eta < _minEta || eta > _maxEta) continue;
    //if (pt < 0.2 || pt > 2.0) continue;
    if (charge == 0) continue;
    hTrackCounter->Fill(4);
    hTrackCounter->GetXaxis()->SetBinLabel(4, "acceptance");

    Float_t nCls = track->GetTPCNcls();
    Int_t id = track->GetID();
    Float_t ptTPC = track->GetTPCmomentum();
    Float_t itsSignal = track->GetITSsignal();
    Float_t tpcSignal = track->GetTPCsignal();
    Float_t tpcSignalN = track->GetTPCsignalN();
    Float_t nSharedCls = track->GetTPCnclsS();
    Float_t nCrossedRows = track->GetTPCCrossedRows();
    Float_t nFindableCls = track->GetTPCNclsF();
    Float_t goldenchi2 = track->GetChi2TPCConstrainedVsGlobal();
    Bool_t isTPConly = track->AliAODTrack::IsTPCOnly();
    Bool_t isTPCconst = track->IsTPCConstrained();
    Bool_t isGlobalconst = track->IsGlobalConstrained();
    Double_t tpcChi2 = track->GetTPCchi2();
    Double_t itsChi2 = track->GetITSchi2();
    Double_t tpcChi2percl = track->GetTPCchi2perCluster();
    Double_t tpcTgl = track->GetTPCTgl();
    Int_t ntracks = fAOD->GetNumberOfTracks();
    Float_t missCls = track->GetTPCClusterInfo(3,0,0,159);
    Float_t nclsITS = track->GetITSNcls();
    Int_t nTPCClusters = fAOD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fAOD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    Bool_t fb128_bug = track->TestFilterBit(128);
    Bool_t fb128 = track->TestFilterBit(128) && (abs(dcaXY) < 2.4 && abs(dcaZ) < 3.2 && nCls > 50 && tpcChi2percl <= 4.0);
    Bool_t fb768 = track->TestFilterBit(768);

    // from note: https://alice-notes.web.cern.ch/node/736
    Bool_t i_tpcRowsClsChi2 = nCrossedRows >= 80 /* && nCls >= 70 */ && tpcChi2percl <= 4.0;
    Bool_t ptDCAxy1 = (abs(dcaXY) < (0.0182 + 0.035 / pow(pt, 1.01))); // run1
    
    // from ilya note (suggested by mesut) https://alice-notes.web.cern.ch/node/1653
    // Bool_t ptDCAxy2 = (abs(dcaXY) < (0.028 + 0.04 * pow(pt, 1.01))); // run2
    // from Ante's email conversation after PAG presentation on 30/Oct/2025
    Bool_t ptDCAxy2 = (abs(dcaXY) < (0.0105 + 0.035 / pow(pt, 1.1))); // run2
    Bool_t i_tpcFracShaClsFindCls = nSharedCls/nCls <= 0.3 && nSharedCls/nCrossedRows <= 0.25 && nFindableCls/nCls >= 0.8;

    Bool_t ptDCAcut = (_YEAR_ == 2010) ? ptDCAxy1 : ptDCAxy2;

    if(fb128) hTrackCounter->Fill(5); hTrackCounter->GetXaxis()->SetBinLabel(5, "fb128)");
    if(fb768) hTrackCounter->Fill(6); hTrackCounter->GetXaxis()->SetBinLabel(6, "fb768");
    if(isTPConly) hTrackCounter->Fill(7); hTrackCounter->GetXaxis()->SetBinLabel(7, "isTPConly");
    if(isTPCconst) hTrackCounter->Fill(8); hTrackCounter->GetXaxis()->SetBinLabel(8, "isTPCconst");
    if(isGlobalconst) hTrackCounter->Fill(9); hTrackCounter->GetXaxis()->SetBinLabel(9, "isGlobalconst");
    if (fb128_bug) hTrackCounter->Fill(10); hTrackCounter->GetXaxis()->SetBinLabel(10, "fb128_bug");
    if (i_tpcRowsClsChi2) hTrackCounter->Fill(11); hTrackCounter->GetXaxis()->SetBinLabel(11, "crows>80 && chi2percl<4");
    if (i_tpcFracShaClsFindCls) hTrackCounter->Fill(12); hTrackCounter->GetXaxis()->SetBinLabel(12, "shcls/ncls<0.3 && shcls/ncrows<0.25 && nfindablecls/ncls>0.8");
    if (ptDCAcut) hTrackCounter->Fill(13); hTrackCounter->GetXaxis()->SetBinLabel(13, Form("ptDCAcut: %s", _YEAR_ == 2010 ? "0.0182 + 0.035 / pow(pt, 1.01)" : "0.0105 + 0.035 / pow(pt, 1.1)"));

    Bool_t _type[_numEnvs] = {kFALSE};
    if (fb128) _type[0] = kTRUE;
    if (fb768) _type[1] = kTRUE;
    if (fb768 && i_tpcRowsClsChi2) _type[2] = kTRUE;
    if (fb768 && i_tpcFracShaClsFindCls) _type[3] = kTRUE;
    if (fb768 && i_tpcRowsClsChi2 && i_tpcFracShaClsFindCls) _type[4] = kTRUE;
    if (fb768 && ptDCAcut && i_tpcRowsClsChi2 && i_tpcFracShaClsFindCls) _type[5] = kTRUE;

    GetPtBin(pt);
    hEnvsCounter->Fill(1);
    hEnvsCounter->GetXaxis()->SetBinLabel(1, "accepted tracks");
    for (Int_t j = 0; j < _numEnvs; ++j) {
    if (_type[j]) {
      hEnvsCounter->Fill(j+2);
      hEnvsCounter->GetXaxis()->SetBinLabel(j+2, Form("type %d", j+1));
      hDCAxy[j]->Fill(dcaXY);
      hDCAz[j]->Fill(dcaZ);
      hDCAxy_imp[j]->Fill(dcaXY_imp);
      hDCAz_imp[j]->Fill(dcaZ_imp);
      hDCAxypT[j]->Fill(pt, dcaXY);
      hDCAzpT[j]->Fill(pt, dcaZ);
      hnITScls[j]->Fill(nclsITS);
      hnITScls2[j]->Fill(nITSClusters);
      hnTPCcls[j]->Fill(nCls);
      hnTPCcls2[j]->Fill(nTPCClusters);
      hnTPCcrossedrows[j]->Fill(nCrossedRows);
      hmissingCls[j]->Fill(missCls);
      htpcTgl[j]->Fill(tpcTgl);
      hNShCls[j]->Fill(nSharedCls);
      hNShClsNcls[j]->Fill(nSharedCls/nCls);
      hNShClsNcrows[j]->Fill(nSharedCls/nCrossedRows);
      hNShClsFra[j]->Fill(nSharedCls/nCls, nSharedCls/nCrossedRows);
      hNFoundClsFra[j]->Fill(nSharedCls/nCls, nCrossedRows/nFindableCls);
      hNShClsVsPt[j]->Fill(pt, nSharedCls/nCls);
      hNFcls[j]->Fill(nFindableCls);
      hNFindClsFra[j]->Fill(nFindableCls/nCls);
      hNFindClsVsPt[j]->Fill(pt, nFindableCls/nCls);
      hTPCsignal[j]->Fill(tpcSignal);
      hTPCsignalN[j]->Fill(tpcSignalN);
      hTPCsignalvsPt[j]->Fill(pt, tpcSignal);
      hTPCsignalvsPtN[j]->Fill(pt, tpcSignalN);
      hTPCsignalvsPtTuned[j]->Fill(pt, dEdxTPC);
      hTPCsignalvsPtot[j]->Fill(ptTPC, dEdxTPC);
      hTPCsignalTuned[j]->Fill(dEdxTPC);
      hITSsignal[j]->Fill(itsSignal);
      hChi2TPC[j]->Fill(tpcChi2);
      hChi2perclTPC[j]->Fill(tpcChi2percl);
      hChi2ITS[j]->Fill(itsChi2/nclsITS);
      hPtDis[j]->Fill(pt);
      hPtotDis[j]->Fill(ptTPC);
      hEtaDis[j]->Fill(eta);
      hPhiDis[j]->Fill(phi);
      hGoldenChi2[j]->Fill(goldenchi2);
      if (_FLAG_MC_) {
	if (TMath::Abs(particle->GetPdgCode()) == 211) foundPion[j]++;
      	if (TMath::Abs(particle->GetPdgCode()) == 321) foundKaon[j]++;
      	if (TMath::Abs(particle->GetPdgCode()) == 2212) foundProton[j]++;
      }

       for (Int_t k = 0; k < _numPtBin; ++k) {
         if (_ptBinFlags[k]) {
           if (j == 0) { hTrackCounter->Fill(14+k); hTrackCounter->GetXaxis()->SetBinLabel(13+k, Form("pt bin %d", k+1)); }
           foundPart[k][j]++;
           hPhiBin[k][j]->Fill(phi);
           hEtaBin[k][j]->Fill(eta);
           hPtBin[k][j]->Fill(pt);
           for (Int_t l = 0; l < _numM; ++l) hEtaPhiBin[k][l][j]->Fill(eta, phi);
          } // end of if ptbin condition
        } // end of ptbin loop
      } // end of if type condition
    } // end of _numEnvs loop

    if (_FLAG_DETADPHI_) {
      for (Int_t ii = i + 1; ii < fAOD->GetNumberOfTracks(); ii++) {
        AliAODTrack *_track = static_cast<AliAODTrack *>(fAOD->GetTrack(ii));
        if (!_track) continue;
        Float_t _pt = _track->Pt();
        Float_t _eta = _track->Eta();
        Float_t _phi = _track->Phi();
        Int_t _charge = _track->Charge();
        Int_t _id = _track->GetID();
        if ((charge == 0) || (_id == id) || (_pt < 0.2) || (_eta < _minEta) || (_eta > _maxEta)) continue;

        Float_t _dPStar = CalculateDPhiStar(phi, eta, pt, charge, _phi, _eta, _pt, _charge, c_field);
        Float_t _dEta = _eta - eta;
        Float_t _dPhi = _phi - phi;
        if (_dPhi > TMath::Pi()) _dPhi -= 2*TMath::Pi();
        if (_dPhi < -TMath::Pi()) _dPhi += 2*TMath::Pi();

        for (Int_t j = 0; j < _numEnvs; ++j) {
          if (_type[j]) {
            hdEta[j]->Fill(_dEta);
            hdPhi[j]->Fill(_dPhi);
            Int_t _ch = -999;
            if (_charge > 0 && _charge > 0) _ch = 0;
            if (_charge < 0 && _charge < 0) _ch = 1;
            if (_ch == -999) continue;

            for (Int_t k = 0; k < _numPtBin; ++k) {
              if (_ptBinFlags[k] && pt > _ptArray[2*k] && pt < _ptArray[2*k+1]) {
                hDEDPraw[k][j]->Fill(_dEta, _dPhi);
                for (Int_t l = 0; l < 2; ++l) {
                  if (_charge != charge) hDEDPrawChDiff[k][l][j]->Fill(_dEta, _dPhi);
                  hDEDPrawChSame[k][l][_ch][j]->Fill(_dEta, _dPhi);
                  hDEDPrawChSame[k][l][2][j]->Fill(_dEta, _dPhi);
                  if (_pt > pt && _dPhi < 0) { hDEDPrawPtOrd[k][0][_ch][l][j]->Fill(_dEta, _dPhi); hDEDPrawPtOrd[k][0][2][l][j]->Fill(_dEta, _dPhi); }
                  if (_pt > pt && _dPhi > 0) { hDEDPrawPtOrd[k][1][_ch][l][j]->Fill(_dEta, _dPhi); hDEDPrawPtOrd[k][1][2][l][j]->Fill(_dEta, _dPhi); }
                  if (_pt < pt && _dPhi < 0) { hDEDPrawPtOrd[k][2][_ch][l][j]->Fill(_dEta, _dPhi); hDEDPrawPtOrd[k][2][2][l][j]->Fill(_dEta, _dPhi); }
                  if (_pt < pt && _dPhi > 0) { hDEDPrawPtOrd[k][3][_ch][l][j]->Fill(_dEta, _dPhi); hDEDPrawPtOrd[k][3][2][l][j]->Fill(_dEta, _dPhi); }
                } // end of loop on two charges
              } //end of loop on ptBin check
            } // end of loop on _numPtBin
          } // end of loop on _numEnvs check
        } // end of loop on _numEnvs
      } // end of loop on ii tracks
    } // end of loop on _FLAG_DETADPHI_

  } // end of track loop

  std::cout << "\033[31m" <<  " ================================= " << std::endl;
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " Filling track information done" << "\033[0m" << std::endl;
  for (Int_t i = 0; i < _numEnvs; ++i) {
    hNumPions[i]->Fill(foundPion[i]);
    hNumKaons[i]->Fill(foundKaon[i]);
    hNumProtons[i]->Fill(foundProton[i]);
    hNumPionsCent[i]->Fill(c_percentile, foundPion[i]);
    hNumKaonsCent[i]->Fill(c_percentile, foundKaon[i]);
    hNumProtonsCent[i]->Fill(c_percentile, foundProton[i]);
    //continue;
    std::cout << "\033[31m" <<  " for condition " << i+1 << " == :: == task marooz == :: == " << " number of pions " << foundPion[i] << " kaons " << foundKaon[i] << " protons " << foundProton[i] << std::endl;
    for (Int_t j = 0; j < _numPtBin; ++j) {
      std::cout << "\033[31m" <<  " ================================= " << std::endl;
      std::cout << "\033[31m" <<  " for pt bin " << j+1 << " == :: == task marooz == :: == " << " number of particles " << foundPart[j][i] << std::endl;
      if (foundPart[j][i] > 0) hMultBin[j][i]->Fill(foundPart[j][i]);
      else std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " No particles found in this bin " << std::endl;
    }
  }

  std::cout << "\033[31m" <<  " ================================= " << std::endl;
  std::cout << "\033[31m" <<  " ================================= " << std::endl;
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " out of FillTrackInfo " << std::endl;
}

void AliAnalysisTaskFM_moments::FillMCTrackInfo() {
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " in FillMCTrackInfo " << std::endl;
  Int_t foundPartGen[_maxPtBin] = {0};
  for (Int_t i_MCtrk = 0; i_MCtrk < fMCEvent->GetNumberOfTracks(); i_MCtrk++) {
    hTrackCounterGen->Fill(1);
    hTrackCounterGen->GetXaxis()->SetBinLabel(1, "All tracks");
    AliVParticle *lPart = (AliAODMCParticle *)fMCEvent->GetTrack(i_MCtrk);
    TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL) return;
    hTrackCounterGen->Fill(2);
    hTrackCounterGen->GetXaxis()->SetBinLabel(2, "AODMCTrackArray");

    Bool_t isoobPileup = kFALSE;
    if (_FLAG_PILEUP_) isoobPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i_MCtrk, mcHeader, AODMCTrackArray);

    Float_t lpt = lPart->Pt();
    Float_t leta = lPart->Eta();
    Float_t lphi = lPart->Phi();
    Int_t lcharge = lPart->Charge();

    if (!lPart || !lPart->IsPhysicalPrimary() || isoobPileup || lpt < 0.2 || (leta < _minEta) || (leta > _maxEta) || lcharge == 0) continue;
    hTrackCounterGen->Fill(3);
    hTrackCounterGen->GetXaxis()->SetBinLabel(3, "acceptance");

    hPtDisGen->Fill(lpt);
    hEtaDisGen->Fill(leta);
    hPhiDisGen->Fill(lphi);

    GetPtBin(lpt);

    for (Int_t iPt = 0; iPt < _numPtBin; ++iPt) {
      if (_ptBinFlags[iPt]) {
        hTrackCounterGen->Fill(4+iPt);
        hTrackCounterGen->GetXaxis()->SetBinLabel(4+iPt, Form("pt bin %d", iPt+1));
        foundPartGen[iPt]++;
        hPhiBinGen[iPt]->Fill(lphi);
        hEtaBinGen[iPt]->Fill(leta);
        hPtBinGen[iPt]->Fill(lpt);
        for (Int_t k = 0; k < _numM; ++k) {
          hEtaPhiBinGen[iPt][k][0]->Fill(leta, lphi);
        } // end of _numM loop
      } // end of pt bin loop
    } // end of pt bin loop
  } // end of MC track loop

  // Fill the multiplicity histograms
  for (Int_t iPt = 0; iPt < _numPtBin; ++iPt) {
    hMultBinGen[iPt]->Fill(foundPartGen[iPt]);
  }
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " out of FillMCTrackInfo " << std::endl;
}

void AliAnalysisTaskFM_moments::GetPtBin(Double_t pt) {
  _ptBinFlags.clear();
  for (Int_t i = 0; i < _ptArray.GetSize() - 1; i += 2) {
    if (pt >= _ptArray[i] && pt <= _ptArray[i + 1])  _ptBinFlags.push_back(kTRUE);
    else  _ptBinFlags.push_back(kFALSE);
  }
}

Float_t AliAnalysisTaskFM_moments::CalculateDPhiStar(Float_t phi1, Float_t eta1, Float_t pt1, Int_t charge1,Float_t phi2, Float_t eta2, Float_t pt2, Int_t charge2, Float_t bSign) {
  Float_t radius = 0;
  Float_t _tmp = 0;
  Float_t dphistar = 999;
  Float_t kLimit = fdeta * 3;
  Double_t kPi = TMath::Pi();
  bSign = (bSign > 0) ? 1 : -1;
  // variables and cuts have been taken from
  // https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
  if (abs(eta1 - eta2) < fdeta * 2.5 * 3) {
    // check first boundaries to see if is worth to loop and find the minimum
    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
    if (fabs(dphistar1) < kLimit || fabs(dphistar2) < kLimit || (dphistar1 * dphistar2) < 0) {
      // find the smallest dphistar (radius of TPC: 0.8 - 2.5 (m))
      for (Int_t rad(80); rad < 251; rad++) {
        radius = rad * 0.01;
        _tmp = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, radius, bSign);
        if (rad == 80) dphistar = _tmp;
        if (_tmp < dphistar) dphistar = _tmp;
      }
      if (dphistar < kPi) dphistar = kPi * 2 - dphistar;
      if (dphistar < -1 * kPi) dphistar = -2 * kPi - dphistar;
      if (dphistar > kPi) dphistar = kPi * 2 - dphistar;
    }
  }
  return dphistar;
}

Float_t AliAnalysisTaskFM_moments::SharedClusterFraction(TBits &cl1, TBits &cl2, TBits &sh1, TBits &sh2) {
  Int_t ncl1 = cl1.GetNbits();
  Int_t ncl2 = cl2.GetNbits();
  Int_t sumCls = 0, sumSha = 0, sumQ = 0;
  Double_t shfrac = 0;
  for (Int_t ich = 0; ich < ncl1 && ich < ncl2; ich++) {
    if (cl1.TestBitNumber(ich) && cl2.TestBitNumber(ich)) {   // Both clusters
      if (sh1.TestBitNumber(ich) && sh2.TestBitNumber(ich)) { // Shared
        sumQ++; sumCls += 2; sumSha += 2;
      } else {
        sumQ--; sumCls += 2;
      }
    } else if (cl1.TestBitNumber(ich) || cl2.TestBitNumber(ich)) { // Non shared
      sumQ++; sumCls++;
    }
  }
  if (sumCls > 0) shfrac = sumSha * 1.0 / sumCls;
  return shfrac;
}

Float_t AliAnalysisTaskFM_moments::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign) {
  Double_t kPi = TMath::Pi();
  Double_t deltaPhi = phi2 - phi1;
  if (deltaPhi < -0.5 * TMath::Pi()) deltaPhi += TMath::TwoPi();
  if (deltaPhi > 1.5 * TMath::Pi()) deltaPhi -= TMath::TwoPi();
  // calculates dphistar
  Float_t dphistarm = deltaPhi + charge1 * bSign * TMath::ASin(0.075 * radius / (2 * pt1)) - charge2 * bSign * TMath::ASin(0.075 * radius / (2 * pt2));
  // circularity
  if (dphistarm > kPi) dphistarm = kPi * 2 - dphistarm;
  if (dphistarm < -kPi) dphistarm = -kPi * 2 - dphistarm;
  if (dphistarm > kPi) // might look funny but is needed
    dphistarm = kPi * 2 - dphistarm;

  return dphistarm;
}

void AliAnalysisTaskFM_moments::CalculateNFMs(Bool_t _isGen) {
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " in CalculateNFMs " << std::endl;

  Int_t _type = _numEnvs;
  if (_isGen) _type = 1;
  for (Int_t i = 0; i < _type; ++i) {
    for (Int_t j = 0; j < _numPtBin; ++j) {
      for (Int_t k = 0; k < _numM; k++) {
      Double_t _nPSBins = 0;
      Double_t _squareM = 0;
      Double_t _sumBinCon = 0;
      Double_t _fqNum[6]; // Fq2e, Fq3e, Fq4e, Fq5e, Fq6e, Fq7e
      Double_t _fqNumBfAvg[6];
      Double_t _fqDenBfAvg= 0.0;
      
      if (_maxPhaseSpaceBins == 123) _nPSBins = 3 * (k + 2);
      if (_maxPhaseSpaceBins == 82) _nPSBins = 2 * (k + 2);
      Double_t Mbin = _nPSBins;
  
      for (Int_t index = 0; index < 6; index++) {
        _fqNum[index] = 0.0;
        _fqNumBfAvg[index] = 0.0;
      }

      TH2F *hist;
      if (_isGen) hist = hEtaPhiBinGen[j][k][0];
      else hist = hEtaPhiBin[j][k][i];
      if (!hist){
        std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " in CalculateNFMs " << " hist not found " << std::endl;
        continue;
      }

      for (Int_t _eb = 1; _eb <= hist->GetNbinsX(); _eb++) {
        for (Int_t _pb = 1; _pb <= hist->GetNbinsY(); _pb++) {
          _fqDenBfAvg= 0.0;
          _fqDenBfAvg= hist->GetBinContent(_eb, _pb);
          _sumBinCon += _fqDenBfAvg;

          for (Int_t q = 0; q < 6; q++) {
            if (_fqDenBfAvg >= (q + 2)) {
              Double_t _fqtmp = 0.0;
              _fqtmp = TMath::Factorial(_fqDenBfAvg) / TMath::Factorial(_fqDenBfAvg- (q + 2));
              if (_FLAG_MC_) {
                if (TMath::IsNaN(_fqtmp)) break;
              }
              _fqNumBfAvg[q] += _fqtmp;
            }
          }
        }
      }

      _squareM = TMath::Power(_nPSBins, 2); // 2 because of two dimensions: eta and phi
      Double_t Av_bincontent = _sumBinCon / _squareM;
      for (Int_t q = 0; q < 6; q++) {
        if (_fqNumBfAvg[q] > 0.0) {
          _fqNum[q] = _fqNumBfAvg[q] / (_squareM);
        }
      }

      Float_t Fq2e = _fqNum[0];
      Float_t Fq3e = _fqNum[1];
      Float_t Fq4e = _fqNum[2];
      Float_t Fq5e = _fqNum[3];
      Float_t Fq6e = _fqNum[4];
      Float_t Fq7e = _fqNum[5];
      Float_t mult = foundPart[j][i];

      if (!_isGen) fntpMBin[j][k][i]->Fill(mult, Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e, Fq6e, Fq7e);
      else fntpMBinGen[j][k]->Fill(Mbin, Av_bincontent, Fq2e, Fq3e, Fq4e, Fq5e, Fq6e, Fq7e);
      } // end of loop on number of phase space bins
    } // end of loop on pt bins
  } // end of Env loop
}

void AliAnalysisTaskFM_moments::SetEventMixing(Bool_t mixed) {
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " in SetEventMixing " << std::endl;
  c_percentile = 0.5;
  _vzMax = 8;
  _pool = fEvPoolMgr->GetEventPool(c_percentile, _vzMax); _pool->PrintInfo();
  TObjArray *fTrack = new TObjArray(); fTrack->SetOwner(kFALSE);

  for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
    AliAODTrack *track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (!track) continue;
    Float_t eta = track->Eta();
    Float_t pt = track->Pt();
    Int_t charge = track->Charge();
    if (eta < _minEta || eta > _maxEta) continue;
    //if (pt < 0.2 || pt > 2.0) continue;
    if (charge == 0) continue;
    fTrack->Add(new AliAODTrack(*track));
  }
  _pool->UpdatePool(fTrack);
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " isready " << _pool->IsReady() << ", " << _pool->NTracksInPool() << ", " << _pool->GetCurrentNEvents() << std::endl;
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " out of SetEventMixing " << std::endl;
}

Bool_t AliAnalysisTaskFM_moments::GetDCA(AliAODTrack *track, Double_t dca[2]) {
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
      return kFALSE;
    }

    AliAODVertex *vtx = (AliAODVertex *)(fAOD->GetPrimaryVertex());
    Double_t fBzkG = fAOD->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx, fBzkG, kVeryBig, dca, covdca);
  }
  if (!ok) {
    dca[0] = -999.;
    dca[1] = -999.;
  }
  return ok;
}

void AliAnalysisTaskFM_moments::DataPosting() {
  for (Int_t i = 0; i < 3; i++) std::cout << "\033[31m =================================================== \033[0m" << std::endl;
  std::cout << "\033[31m" << ":: == task marooz == :: == " << " Data Posting for condition " << _numEnvs << " and ptBins " << _numPtBin << " \033[0m" << std::endl;
  for (Int_t i = 0; i < 3; i++) std::cout << "\033[31m =================================================== \033[0m" << std::endl;
  
  for (Int_t i = 0; i < _numEnvs; i++) PostData(i + 1, hList[i]);
  for (Int_t i = 0; i < _numEnvs; i++) for (Int_t j = 0; j < _numPtBin; j++) PostData(1 + _numEnvs + i * _numPtBin + j, fNtupleList[j][i]); 
  if (_FLAG_MC_) for (Int_t j = 0; j < _numPtBin; j++) PostData(1 + _numEnvs + _numEnvs * _numPtBin + j, fNtupleListGen[j]);
}

void AliAnalysisTaskFM_moments::Terminate(Option_t *) {
  std::cout << "\033[31m" <<  " == :: == task marooz == :: == " << " TASK TERMINATED " << std::endl;
}
