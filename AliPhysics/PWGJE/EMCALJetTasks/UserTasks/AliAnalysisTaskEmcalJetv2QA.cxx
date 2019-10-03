// Jet v2 task using QA method, based on jet sample task (S.Aiola).
//
// Authors: Jason Mueller (CERN summer student 2014) & Alice Ohlson


#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TF1.h>

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskEmcalJetv2QA.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetv2QA)

//________________________________________________________________________
AliAnalysisTaskEmcalJetv2QA::AliAnalysisTaskEmcalJetv2QA() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetv2QA", kTRUE),
  nCentBins(0),
  nCentBins1(1),
  centBins(0x0),
  nJetPtBins(0),
  nJetPtBins1(1),
  jetPtBins(0x0),
  fJetv2(0x0),
  doPtWeight(kFALSE),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistLeadingJetPtCorr(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fDPhiJet(0),
  fDPhiJetPythia(0),
  fDPhiEP(0),
  hGx(0),
  hGy2(0),
  hGxGy2(0),
  hGy4(0),
  hGx2(0),
  hGx2Gy2(0),
  hGxGy4(0),
  hGy6(0),
  hGx2Gy4(0),
  hGxGy6(0),
  hGy8(0),
  hGy(0),
  hN(0),
  htv2std(0),
  htjv2std(0),
  htj2v2std(0),
  hV0jv2std(0),
  htdPsi(0),
  htjdPsi(0),
  htj2dPsi(0),
  hV0jdPsi(0),
  hAx(0),
  hAxDijet(0),
  hQx(0),
  hQy(0),
  hEventData(0),
  hNTracks(0),
  hNTracksCent(0),
  hGxTracks(0),
  hGyTracks(0),
  hGy2Tracks(0),
  hGxGy2Tracks(0),
  hGy4Tracks(0),
  htEPRes(0),
  htjEPRes(0),
  htj2EPRes(0),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetv2QA::AliAnalysisTaskEmcalJetv2QA(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  nCentBins(0),
  nCentBins1(1),
  centBins(0x0),
  nJetPtBins(0),
  nJetPtBins1(1),
  jetPtBins(0x0),
  fJetv2(0x0),
  doPtWeight(kFALSE),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistLeadingJetPtCorr(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fDPhiJet(0),
  fDPhiJetPythia(0),
  fDPhiEP(0),
  hGx(0),
  hGy2(0),
  hGxGy2(0),
  hGy4(0),
  hGx2(0),
  hGx2Gy2(0),
  hGxGy4(0),
  hGy6(0),
  hGx2Gy4(0),
  hGxGy6(0),
  hGy8(0),
  hGy(0),
  hN(0),
  htv2std(0),
  htjv2std(0),
  htj2v2std(0),
  hV0jv2std(0),
  htdPsi(0),
  htjdPsi(0),
  htj2dPsi(0),
  hV0jdPsi(0),
  hAx(0),
  hAxDijet(0),
  hQx(0),
  hQy(0),
  hEventData(0),
  hNTracks(0),
  hNTracksCent(0),
  hGxTracks(0),
  hGyTracks(0),
  hGy2Tracks(0),
  hGxGy2Tracks(0),
  hGy4Tracks(0),
  htEPRes(0),
  htjEPRes(0),
  htj2EPRes(0),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Standard constructor.

  // default binning
  Double_t centBinsTemp[7] = {0,5,10,20,30,40,50};
  SetCentBins(6,centBinsTemp);
  Double_t jetPtBinsTemp[6] = {40.,50.,70.,90.,120.,200.};
  SetJetPtBins(5,jetPtBinsTemp);

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetv2QA::~AliAnalysisTaskEmcalJetv2QA()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetv2QA::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }

  fTracksCont->SetClassName("AliVTrack");
  fCaloClustersCont->SetClassName("AliVCluster");

  fHistTracksPt = new TH1F("fHistTracksPt","fHistTracksPt", fNbins / 2, fMinBinPt, fMaxBinPt / 2);
  fHistTracksPt->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
  fHistTracksPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistTracksPt);

  fHistClustersPt = new TH1F("fHistClustersPt", "fHistClustersPt", fNbins / 2, fMinBinPt, fMaxBinPt / 2);
  fHistClustersPt->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
  fHistClustersPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClustersPt);

  fHistLeadingJetPt = new TH1F("fHistLeadingJetPt", "fHistLeadingJetPt", fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJetPt->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
  fHistLeadingJetPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJetPt);

  fHistLeadingJetPtCorr = new TH1F("fHistLeadingJetPtCorr", "fHistLeadingJetPtCorr", fNbins, fMinBinPt, fMaxBinPt);
  fHistLeadingJetPtCorr->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
  fHistLeadingJetPtCorr->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJetPtCorr);

  fHistJetsPhiEta = new TH2F("fHistJetsPhiEta", "fHistJetsPhiEta", 50, -1, 1, 101, 0, TMath::Pi()*2 + TMath::Pi()/200);
  fHistJetsPhiEta->GetXaxis()->SetTitle("#eta");
  fHistJetsPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJetsPhiEta);

  fHistJetsPtArea = new TH2F("fHistJetsPtArea", "fHistJetsPtArea", fNbins, fMinBinPt, fMaxBinPt, 30, 0, 3);
  fHistJetsPtArea->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
  fHistJetsPtArea->GetYaxis()->SetTitle("area");
  fOutput->Add(fHistJetsPtArea);

  fHistJetsPtLeadHad = new TH2F("fHistJetsPtLeadHad", "fHistJetsPtLeadHad", fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, fMinBinPt, fMaxBinPt / 2);
  fHistJetsPtLeadHad->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
  fHistJetsPtLeadHad->GetYaxis()->SetTitle("p_{T,lead} (GeV/c)");
  fHistJetsPtLeadHad->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPtLeadHad);

  fHistJetsCorrPtArea = new TH2F("fHistJetsCorrPtArea", "fHistJetsCorrPtArea", fNbins*2, -fMaxBinPt, fMaxBinPt, 30, 0, 3);
  fHistJetsCorrPtArea->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
  fHistJetsCorrPtArea->GetYaxis()->SetTitle("area");
  fOutput->Add(fHistJetsCorrPtArea);

  fHistPtDEtaDPhiTrackClus = new TH3F("fHistPtDEtaDPhiTrackClus","fHistPtDEtaDPhiTrackClus;#it{p}_{T}^{track};#Delta#eta;#Delta#varphi",100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);

  fHistPtDEtaDPhiClusTrack = new TH3F("fHistPtDEtaDPhiClusTrack","fHistPtDEtaDPhiClusTrack;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);

  fDPhiJet = new TH1F("fDPhiJet","fDPhiJet",90, -TMath::Pi()/3, 5*TMath::Pi()/3);
  fDPhiJet->GetXaxis()->SetTitle("#Delta#varphi");
  fDPhiJet->GetYaxis()->SetTitle("counts");
  fOutput->Add(fDPhiJet);

  fDPhiJetPythia = new TH1F("fDPhiJetPythia","fDPhiJetPythia",90, -TMath::Pi()/3, 5*TMath::Pi()/3);
  fDPhiJetPythia->GetXaxis()->SetTitle("#Delta#varphi");
  fDPhiJetPythia->GetYaxis()->SetTitle("counts");
  fOutput->Add(fDPhiJetPythia);

  fDPhiEP = new TH1F("fDPhiEP","fDPhiEP",90, 0, 2*TMath::Pi());
  fDPhiEP->GetXaxis()->SetTitle("#Delta#varphi");
  fDPhiEP->GetYaxis()->SetTitle("counts");
  fOutput->Add(fDPhiEP);

  hGx = new TH2D("hGx", "Gx v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGx->GetXaxis()->SetTitle("Centrality (%)");
  hGx->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGx);

  hGy2 = new TH2D("hGy2", "Gy2 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGy2->GetXaxis()->SetTitle("Centrality (%)");
  hGy2->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGy2);

  hGxGy2 = new TH2D("hGxGy2", "GxGy2 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGxGy2->GetXaxis()->SetTitle("Centrality (%)");
  hGxGy2->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGxGy2);

  hGy4 = new TH2D("hGy4", "Gy4 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGy4->GetXaxis()->SetTitle("Centrality (%)");
  hGy4->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGy4);

  hGx2 = new TH2D("hGx2", "Gx2 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGx2->GetXaxis()->SetTitle("Centrality (%)");
  hGx2->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGx2);

  hGx2Gy2 = new TH2D("hGx2Gy2", "Gx2Gy2 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGx2Gy2->GetXaxis()->SetTitle("Centrality (%)");
  hGx2Gy2->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGx2Gy2);

  hGxGy4 = new TH2D("hGxGy4", "GxGy4 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGxGy4->GetXaxis()->SetTitle("Centrality (%)");
  hGxGy4->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGxGy4);

  hGy6 = new TH2D("hGy6", "Gy6 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGy6->GetXaxis()->SetTitle("Centrality (%)");
  hGy6->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGy6);

  hGx2Gy4 = new TH2D("hGx2Gy4", "Gx2Gy4 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGx2Gy4->GetXaxis()->SetTitle("Centrality (%)");
  hGx2Gy4->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGx2Gy4);

  hGxGy6 = new TH2D("hGxGy6", "GxGy6 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGxGy6->GetXaxis()->SetTitle("Centrality (%)");
  hGxGy6->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGxGy6);

  hGy8 = new TH2D("hGy8", "Gy8 v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGy8->GetXaxis()->SetTitle("Centrality (%)");
  hGy8->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGy8);

  hGy = new TH2D("hGy", "Gy v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hGy->GetXaxis()->SetTitle("Centrality (%)");
  hGy->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hGy);

  hN = new TH2D("hN", "N v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hN->GetXaxis()->SetTitle("Centrality (%)");
  hN->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hN);

  htv2std = new TH2D("htv2std", "v2std v Centrality v JetPt w/o jets", nCentBins, centBins, nJetPtBins, jetPtBins);
  htv2std->GetXaxis()->SetTitle("Centrality (%)");
  htv2std->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(htv2std);

  htjv2std = new TH2D("htjv2std", "v2std v Centrality v JetPt w/ jets", nCentBins, centBins, nJetPtBins, jetPtBins);
  htjv2std->GetXaxis()->SetTitle("Centrality (%)");
  htjv2std->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(htjv2std);

  htj2v2std = new TH2D("htj2v2std", "v2std v Centrality v JetPt w/ trackPt < 2 GeV", nCentBins, centBins, nJetPtBins, jetPtBins);
  htj2v2std->GetXaxis()->SetTitle("Centrality (%)");
  htj2v2std->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(htj2v2std);

  hV0jv2std = new TH2D("hV0jv2std", "v2std v Centrality v JetPt", nCentBins, centBins, nJetPtBins, jetPtBins);
  hV0jv2std->GetXaxis()->SetTitle("Centrality (%)");
  hV0jv2std->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hV0jv2std);

  Int_t ndpsibins = 100;
  Double_t dpsibins[101];
  for(Int_t t = 0; t < 101; t++) dpsibins[t] = TMath::Pi()*t/50.;

  htdPsi = new TH3D("htdPsi", "JetAxis - EventPlane w/o jets", nCentBins, centBins, nJetPtBins, jetPtBins, ndpsibins, dpsibins);
  htdPsi->GetZaxis()->SetTitle("#Psi_{jet} - #Psi_{EP}");
  fOutput->Add(htdPsi);

  htjdPsi = new TH3D("htjdPsi", "JetAxis - EventPlane w/ jets", nCentBins, centBins, nJetPtBins, jetPtBins, ndpsibins, dpsibins);
  htjdPsi->GetZaxis()->SetTitle("#Psi_{jet} - #Psi_{EP}");
  fOutput->Add(htjdPsi);

  htj2dPsi = new TH3D("htj2dPsi", "JetAxis - EventPlane w/ trackPt < 2 GeV", nCentBins, centBins, nJetPtBins, jetPtBins, ndpsibins, dpsibins);
  htj2dPsi->GetZaxis()->SetTitle("#Psi_{jet} - #Psi_{EP}");
  fOutput->Add(htj2dPsi);

  hV0jdPsi = new TH3D("hV0jdPsi", "JetAxis - EventPlane", nCentBins, centBins, nJetPtBins, jetPtBins, ndpsibins, dpsibins);
  hV0jdPsi->GetZaxis()->SetTitle("#Psi_{jet} - #Psi_{EP}");
  fOutput->Add(hV0jdPsi);

  hQx = new TH2D("hQx", "Qx Distribution in EP frame", 100, -0.3, 0.3, nCentBins, centBins);
  hQx->GetXaxis()->SetTitle("Qx");
  hQx->GetYaxis()->SetTitle("Centrality (%)");
  fOutput->Add(hQx);

  hQy = new TH2D("hQy", "Qy Distribution in EP frame", 100, -0.3, 0.3, nCentBins, centBins);
  hQy->GetXaxis()->SetTitle("Qy");
  hQy->GetYaxis()->SetTitle("Centrality (%)");
  fOutput->Add(hQy);

  hAx = new TH2D("hAx", "Ax Distribution in EP frame w/o Dijets", 100, -35, 70, nJetPtBins, jetPtBins);
  hAx->GetXaxis()->SetTitle("Ax");
  hAx->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hAx);

  hAxDijet = new TH2D("hAxDijet", "Ax Distribution in EP frame w/ Dijets", 100, -35, 70, nJetPtBins, jetPtBins);
  hAxDijet->GetXaxis()->SetTitle("Ax");
  hAxDijet->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(hAxDijet);

  hEventData = new TH1F("hEventData","Events Kept and Discarded", 9, 0, 9);
  hEventData->GetYaxis()->SetTitle("counts");
  fOutput->Add(hEventData);

  hNTracks = new TH2F("hNTracks","Number of Tracks Per Event", 100, 0, 3000, 3, 0, 3);
  hNTracks->GetXaxis()->SetTitle("# tracks");
  fOutput->Add(hNTracks);

  hNTracksCent = new TH2D("hNTracksCent", "NTracks by centrality", 100, 0, 3000, (Int_t)centBins[nCentBins], centBins[0], centBins[nCentBins]);
  hNTracksCent->GetXaxis()->SetTitle("# tracks");
  hNTracksCent->GetYaxis()->SetTitle("Centrality (%)");
  fOutput->Add(hNTracksCent);

  hGxTracks = new TH2D("hGxTracks", "Gx by NTracks", 200, -200, 200, 100, 0, 3000);
  hGxTracks->GetXaxis()->SetTitle("Gx");
  hGxTracks->GetYaxis()->SetTitle("# tracks");
  fOutput->Add(hGxTracks);

  hGyTracks = new TH2D("hGyTracks", "Gy by NTracks", 200, -200, 200, 100, 0, 3000);
  hGyTracks->GetXaxis()->SetTitle("Gy");
  hGyTracks->GetYaxis()->SetTitle("# tracks");
  fOutput->Add(hGyTracks);

  hGy2Tracks = new TH2D("hGy2Tracks", "Gy2 by NTracks", 100, 0, 20000, 100, 0, 3000);
  hGy2Tracks->GetXaxis()->SetTitle("Gy2");
  hGy2Tracks->GetYaxis()->SetTitle("# tracks");
  fOutput->Add(hGy2Tracks);

  hGxGy2Tracks = new TH2D("hGxGy2Tracks", "GxGy2 by NTracks", 100, -100000, 100000, 100, 0, 3000);
  hGxGy2Tracks->GetXaxis()->SetTitle("GxGy2");
  hGxGy2Tracks->GetYaxis()->SetTitle("# tracks");
  fOutput->Add(hGxGy2Tracks);

  hGy4Tracks = new TH2D("hGy4Tracks", "Gy4 by NTracks", 100, 0, 100000000, 100, 0, 3000);
  hGy4Tracks->GetXaxis()->SetTitle("Gy4");
  hGy4Tracks->GetYaxis()->SetTitle("# tracks");
  fOutput->Add(hGy4Tracks);

  htEPRes = new TH2D("htEPRes", "EP Resolution w/o Jets", nCentBins, centBins, nJetPtBins, jetPtBins);
  htEPRes->GetXaxis()->SetTitle("Centrality (%)");
  htEPRes->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(htEPRes);

  htjEPRes = new TH2D("htjEPRes", "EP Resolution w/ Jets", nCentBins, centBins, nJetPtBins, jetPtBins);
  htjEPRes->GetXaxis()->SetTitle("Centrality (%)");
  htjEPRes->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(htjEPRes);

  htj2EPRes = new TH2D("htj2EPRes", "EP Resolution w/ trackPT < 2 GeV", nCentBins, centBins, nJetPtBins, jetPtBins);
  htj2EPRes->GetXaxis()->SetTitle("Centrality (%)");
  htj2EPRes->GetYaxis()->SetTitle("Leading Jet pT (GeV/c)");
  fOutput->Add(htj2EPRes);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetv2QA::FillHistograms()
{
  // Fill histograms.

  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskEmcalJetv2QA::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetv2QA::Run() // this part loops over each event
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  Double_t jetPhi = -999; // angle of leading jet
  Double_t jetPt = -999; // pt of leading jet
  Double_t jetArea = -999;
  Double_t trackPt = -999;
  Double_t phi = -999; // track phi
  Double_t dPhi = -999; // track phi - jet axis
  Double_t dPhiQA = -999; // track phi - EP
  Double_t tSin = 0; // used for std EP calc
  Double_t tCos = 0; // used for std EP calc
  Double_t jSin = 0; // used for std EP calc
  Double_t jCos = 0; // used for std EP calc
  Double_t tSin2 = 0; // used for std EP calc with trackPt < 2 GeV
  Double_t tCos2 = 0; // used for std EP calc with trackPt < 2 GeV
  Double_t qx = 0; // used for Qx distribution
  Double_t qy = 0; // used for Qy distribution
  Double_t ax = 0; // used for Ax distribution
  Double_t tEP = -999; // EP w/o jets
  Double_t tjEP = -999; // EP w/ jets
  Double_t tjEP2 = -999; // EP w/ jets
  Double_t dPsi = -999; // jet axis - EP
  Double_t gx = 0; // used for G moment calc
  Double_t gy = 0; // used for G moment calc
  Int_t isDijet = 0; // if 0, no dijet. if 1, dijet.
  Int_t nTracksBkgnd = 0; // used to keep track of number of tracks in background
  Int_t nTracksJet =0; // used to keep track of number of tracks in jets

  if(fJetsCont)
    {
      fJetsCont->ResetCurrentID();
      AliEmcalJet *jettest = fJetsCont->GetNextAcceptJet();
      while(jettest) {
      
	fHistJetsPtArea->Fill(jettest->Pt(), jettest->Area());
	fHistJetsPhiEta->Fill(jettest->Eta(), jettest->Phi());
      
	Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jettest);
	fHistJetsPtLeadHad->Fill(jettest->Pt(), ptLeading);
      
	if (fHistJetsCorrPtArea) {
	  Float_t corrPt = jettest->Pt() - fJetsCont->GetRhoVal() * jettest->Area();
	  fHistJetsCorrPtArea->Fill(corrPt, jettest->Area());
	}
	jettest = fJetsCont->GetNextAcceptJet();
      }


      AliEmcalJet *jet = fJetsCont->GetLeadingJet();
      if(jet)
	{
	  jetPhi = jet->Phi(); // get leading jet phi (jet axis)
	  jetPt = jet->Pt(); // get leading jet pT to filter out low pT jet events
	  jetArea = jet->Area();
	}
    }

  // event cuts
  if(!fTracksCont)
    {
      hEventData->Fill("!fTracksCont",1);
      return kFALSE;
    }
  if(!fJetsCont)
    {
      hEventData->Fill("!fJetsCont",1);
      return kFALSE;
    }
  if(jetPt == -999)
    {
      hEventData->Fill("jetPt=-999",1);
      return kFALSE;
    }
  if(jetPt < jetPtBins[0])
    {
      hEventData->Fill("leadingJetPt<jetPtMin",1);
      return kFALSE;
    }
  if(jetPt > jetPtBins[nJetPtBins])
    {
      hEventData->Fill("leadingJetPt>jetPtMax",1);
      return kFALSE;
    }
  if(fCent < centBins[0])
    {
      hEventData->Fill("cent<centMin",1);
      return kFALSE;
    }
  if(fCent > centBins[nCentBins])
    {
      hEventData->Fill("cent>centMax",1);
      return kFALSE;
    }
  hEventData->Fill("good event",1);

  fJetsCont->ResetCurrentID();
  AliEmcalJet *dijet = fJetsCont->GetNextAcceptJet(); // check for dijet events
  while(dijet)
    {
      if(dijet->Pt() > jetPt*2./3. && fabs(jetPhi-dijet->Phi()-TMath::Pi()) < 0.4) // loop over jets with pT>50 and exclude leading jet and check that angular separation is < 0.4
	isDijet = 1;
      dijet = fJetsCont->GetNextAcceptJet();
    }

  if (fCaloClustersCont)
    {
      fCaloClustersCont->ResetCurrentID();
      AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster();
      while(cluster) {
	TLorentzVector nPart;
	cluster->GetMomentum(nPart, fVertex);
	fHistClustersPt->Fill(nPart.Pt());
	cluster = fCaloClustersCont->GetNextAcceptCluster();
      }
    }

  fHistLeadingJetPt->Fill(jetPt);
  fHistLeadingJetPtCorr->Fill(jetPt-fJetsCont->GetRhoVal()*jetArea);

  fTracksCont->ResetCurrentID();
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  while(track)
    { // loop over all particles (including jet tracks)
      trackPt = track->Pt();
      fHistTracksPt->Fill(trackPt);
      phi = track->Phi(); // get track phi

      dPhi = phi-jetPhi; // get track phi - jet axis
      if(dPhi < 0) dPhi += TMath::TwoPi();
      if(dPhi > TMath::TwoPi()) dPhi -= TMath::TwoPi();

      dPhiQA = phi-fEPV0; // get track phi - EP
      if(dPhiQA < 0) dPhiQA += TMath::TwoPi();
      if(dPhiQA > TMath::TwoPi()) dPhiQA -= TMath::TwoPi();
      fDPhiEP->Fill(dPhiQA);

      // fill jet-hadron correlation just to check if track labels make sense...
      if(track->Pt()>1.)
	{
	  Double_t dphiJet = dPhi;
	  if(dphiJet > 5*TMath::Pi()/3) dphiJet -= 2*TMath::Pi();
	  if(dphiJet < -TMath::Pi()/3) dphiJet += 2*TMath::Pi();
	  if(track->GetLabel() == 0)
	    fDPhiJet->Fill(dphiJet);
	  else
	    fDPhiJetPythia->Fill(dphiJet);
	}

      Double_t weight = 1.;
      if(doPtWeight) weight = trackPt;

      gx += weight*cos(2*dPhi);
      gy += weight*sin(2*dPhi);
    
      if(track->GetLabel() == 0)
	{ // sum for std EP method
	  tSin += weight*sin(2*phi); // bkgnd has label = 0
	  tCos += weight*cos(2*phi);
	  qx += weight*cos(2*dPhiQA);
	  qy += weight*sin(2*dPhiQA);
	  nTracksBkgnd += 1;
	}
      else
	{
	  jSin += weight*sin(2*phi); // jets have label =/= 0
	  jCos += weight*cos(2*phi);
	  ax += weight*cos(2*dPhi);
	  nTracksJet += 1;
	}
    
      if(track->Pt() < 2.)
	{ // sum for std EP method w/ trackPt < 2 GeV
	  tSin2 += weight*sin(2*phi);
	  tCos2 += weight*cos(2*phi);
	}
    

      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); // increment to next track
    } // close loop over particles

  hNTracks->Fill(nTracksBkgnd,"Bkgnd Tracks",1);
  hNTracks->Fill(nTracksJet,"Jet Tracks",1);
  hNTracks->Fill(nTracksBkgnd+nTracksJet,"Total Tracks",1);
  hNTracksCent->Fill(nTracksBkgnd+nTracksJet,fCent,1);

  if(nTracksBkgnd == 0)
    {
      hEventData->Fill("no tracks",1);
      hEventData->Fill("good event",-1);
      return kFALSE;
    }

  Double_t v2weight = 1+2*fJetv2*cos(2*(jetPhi-fEPV0)); // set v2 weight for event

  hGx->Fill(fCent, jetPt, v2weight*gx); // fill histograms for QA method
  hGy2->Fill(fCent, jetPt, v2weight*gy*gy);
  hGxGy2->Fill(fCent, jetPt, v2weight*gx*gy*gy);
  hGy4->Fill(fCent, jetPt, v2weight*gy*gy*gy*gy);
  hGx2->Fill(fCent, jetPt, v2weight*gx*gx);
  hGx2Gy2->Fill(fCent, jetPt, v2weight*gx*gx*gy*gy);
  hGxGy4->Fill(fCent, jetPt, v2weight*gx*gy*gy*gy*gy);
  hGy6->Fill(fCent, jetPt, v2weight*gy*gy*gy*gy*gy*gy);
  hGx2Gy4->Fill(fCent, jetPt, v2weight*gx*gx*gy*gy*gy*gy);
  hGxGy6->Fill(fCent, jetPt, v2weight*gx*gy*gy*gy*gy*gy*gy);
  hGy8->Fill(fCent, jetPt, v2weight*gy*gy*gy*gy*gy*gy*gy*gy);
  hGy->Fill(fCent, jetPt, v2weight*gy);
  hN->Fill(fCent, jetPt, v2weight);

  hGxTracks->Fill(gx,nTracksBkgnd+nTracksJet);
  hGyTracks->Fill(gy,nTracksBkgnd+nTracksJet);
  hGy2Tracks->Fill(gy*gy,nTracksBkgnd+nTracksJet);
  hGxGy2Tracks->Fill(gx*gy*gy,nTracksBkgnd+nTracksJet);
  hGy4Tracks->Fill(gy*gy*gy*gy,nTracksBkgnd+nTracksJet);
  hQx->Fill(qx/nTracksBkgnd,fCent);
  hQy->Fill(qy/nTracksBkgnd,fCent);

  if(isDijet == 0)
    hAx->Fill(ax,jetPt);
  if(isDijet == 1)
    hAxDijet->Fill(ax,jetPt);

  tEP = 0.5*atan2(tSin,tCos); // calculate EP w/o jets
  tjEP = 0.5*atan2((tSin+jSin),(tCos+jCos)); // calculate EP w/ jets
  tjEP2 = 0.5*atan2(tSin2,tCos2); // calculate EP w/ trackPt < 2 GeV

  htEPRes->Fill(fCent, jetPt, v2weight*cos(2*(tEP-fEPV0)));
  htjEPRes->Fill(fCent, jetPt, v2weight*cos(2*(tjEP-fEPV0)));
  htj2EPRes->Fill(fCent, jetPt, v2weight*cos(2*(tjEP2-fEPV0)));


  dPsi = jetPhi-tEP;
  if(dPsi < 0) dPsi += TMath::TwoPi();
  if(dPsi > TMath::TwoPi()) dPsi -= TMath::TwoPi();

  htv2std->Fill(fCent, jetPt, v2weight*cos(2*dPsi)); // fill histogram with v2 data w/o jets
  htdPsi->Fill(fCent,jetPt,dPsi); // fill histogram with jet axis - EP w/o jets

  dPsi = jetPhi-tjEP;
  if(dPsi < 0) dPsi += TMath::TwoPi();
  if(dPsi > TMath::TwoPi()) dPsi -= TMath::TwoPi();

  htjv2std->Fill(fCent, jetPt, v2weight*cos(2*dPsi)); // fill histogram with v2 data w/ jets
  htjdPsi->Fill(fCent,jetPt,dPsi); // fill histogram with jet axis - EP w/ jets

  dPsi = jetPhi-tjEP2;
  if(dPsi < 0) dPsi += TMath::TwoPi();
  if(dPsi > TMath::TwoPi()) dPsi -= TMath::TwoPi();

  htj2v2std->Fill(fCent, jetPt, v2weight*cos(2*dPsi)); // fill histogram with v2 data w/ trackPt < 2 GeV
  htj2dPsi->Fill(fCent,jetPt,dPsi); // fill histogram with jet axis - EP w/ trackPt < 2 GeV

  dPsi = jetPhi-fEPV0;
  if(dPsi < 0) dPsi += TMath::TwoPi();
  if(dPsi > TMath::TwoPi()) dPsi -= TMath::TwoPi();

  hV0jv2std->Fill(fCent, jetPt, v2weight*cos(2*dPsi)); // fill histogram with v2 data
  hV0jdPsi->Fill(fCent,jetPt,dPsi); // fill histogram with jet axis - EPV0


  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetv2QA::Terminate(Option_t *)
{
  // Called once at the end of the analysis.
  if(centBins) delete [] centBins;
  if(jetPtBins) delete [] jetPtBins;
}

void AliAnalysisTaskEmcalJetv2QA::SetCentBins(Int_t n, Double_t* bins)
{
  if(centBins) delete [] centBins;
  nCentBins = n;
  nCentBins1 = n+1;
  centBins = new Double_t[nCentBins+1];
  for(Int_t i = 0; i < nCentBins+1; i++)
    centBins[i]=bins[i];
  cout << endl << "Setting " << nCentBins << " centrality bins: " << endl;
  for(Int_t i = 0; i < nCentBins+1; i++)
    cout << centBins[i] << "   ";
  cout << endl << endl;
}

void AliAnalysisTaskEmcalJetv2QA::SetJetPtBins(Int_t n, Double_t* bins)
{
  if(jetPtBins) delete [] jetPtBins;
  nJetPtBins = n;
  nJetPtBins1 = n+1;
  jetPtBins = new Double_t[nJetPtBins+1];
  for(Int_t i = 0; i < nJetPtBins+1; i++)
    jetPtBins[i]=bins[i];
  cout << endl << "Setting " << nJetPtBins << " jet pt bins: " << endl;
  for(Int_t i = 0; i < nJetPtBins+1; i++)
    cout << jetPtBins[i] << "   ";
  cout << endl << endl;
}
