// @(#)root/base:$Id$
// Authors: Alexander Borissov, Sergei Solokhin, Nikita Gladin    03/01/23

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliRun.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSigmaPCMPHOS.h"
#include "AliTriggerAnalysis.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "AliAODMCParticle.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliAODInputHandler.h"
#include "TMath.h"
#include "TPDGCode.h"
#include "AliAODTrack.h"
#include "AliMCEventHandler.h"

// Analysis task to fill histograms with PHOS AOD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskSigmaPCMPHOS)

    ////////////////////////////////////////////////////////////////////////////////
    //Default class constructor
    AliAnalysisTaskSigmaPCMPHOS::AliAnalysisTaskSigmaPCMPHOS() : AliAnalysisTaskSE(),
                                                                 fAODEvent(0), fMCEvent(0x0),
                                                                 fPIDResponse(0),
                                                                 fOutputList(0),
                                                                 fCutsList(0),
                                                                 fMCList(0),
                                                                 fMCPhotonArray(0x0),
                                                                 fMCLambdaArray(0x0),
                                                                 fMCPHOSArray(0x0),
                                                                 fAODv0(0x0), fV0(0x0), ntrack1(0x0), ptrack1(0x0),
                                                                 electronCandidate(0x0),
                                                                 positronCandidate(0x0),
                                                                 fPHOSGeo(0x0),
                                                                 fTriggerAnalysis(new AliTriggerAnalysis),
                                                                 fGamma(0x0),
                                                                 fPCM(0x0),
                                                                 fLambda(0x0),
                                                                 fMixGamma(0x0),
                                                                 fMixLambda(0x0),
                                                                 fAODMCTrackArray(0x0), fConvPhotonArray(0x0)
{
}

////////////////////////////////////////////////////////////////////////////////
// Class constructor for I/O operations
AliAnalysisTaskSigmaPCMPHOS::AliAnalysisTaskSigmaPCMPHOS(const char *name) : AliAnalysisTaskSE(name),
                                                                             fAODEvent(0), fMCEvent(0x0),
                                                                             fPIDResponse(0),
                                                                             fOutputList(0),
                                                                             fCutsList(0),
                                                                             fMCList(0),
                                                                             fMCPhotonArray(0x0),
                                                                             fMCLambdaArray(0x0),
                                                                             fMCPHOSArray(0x0),
                                                                             fAODv0(0x0), fV0(0x0), ntrack1(0x0), ptrack1(0x0),
                                                                             electronCandidate(0x0),
                                                                             positronCandidate(0x0),
                                                                             fPHOSGeo(0x0),
                                                                             fTriggerAnalysis(new AliTriggerAnalysis),
                                                                             fGamma(0x0),
                                                                             fPCM(0x0),
                                                                             fLambda(0x0),
                                                                             fMixGamma(0x0),
                                                                             fMixLambda(0x0),
                                                                             fAODMCTrackArray(0x0), fConvPhotonArray(0x0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}

////////////////////////////////////////////////////////////////////////////////
// Default class destructor
AliAnalysisTaskSigmaPCMPHOS::~AliAnalysisTaskSigmaPCMPHOS()
{
  if (fOutputList)
    delete fOutputList;
  if (fCutsList)
    delete fCutsList;
  if (fMCList)
    delete fMCList;
}

////////////////////////////////////////////////////////////////////////////////
// User-defined output objects, called once
void AliAnalysisTaskSigmaPCMPHOS::UserCreateOutputObjects()
{
  Double_t fLambdaMass = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); //3122 is Lambda, 3212 is Sigma

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fCutsList = new TList();
  fCutsList->SetOwner(kTRUE);

  fMCList = new TList();
  fMCList->SetOwner(kTRUE);

  fMixGamma = new TList();
  fMixGamma->SetOwner(kTRUE);

  fMixLambda = new TList();
  fMixLambda->SetOwner(kTRUE);

  //Gamma-gamma
  const Double_t fMaxEta = 1.2;
  const Int_t nEtaBins = 22;

  const Int_t nMassBins = 100;
  const Double_t fMaxMass = 1.22;
  const Double_t fMinMass = 1.12;
  const Double_t fMaxLambdaMass = 1.22;
  const Double_t fMinLambdaMass = 1.1;
  const Double_t fMaxPhotonMass = 1;
  const Double_t fMinPhotonMass = 0;

  const Int_t nPtBins = 100;
  const Double_t fMaxPt = 10;
  const Double_t fMinPt = 0;
  const Double_t fMaxPtLambda = 10;
  const Double_t fMaxPtPhoton = 4;

  const Int_t nSigmaBins = 100;
  const Double_t fMaxNSigma = 10;

  const Double_t fMaxPHOStime = 1.e-6;

  // 1D histos
  fOutputList->Add(new TH1I("hRunNumbers", "Number of events in runs;run number;N;N_{entries}", 430, 252235, 294925));
  fOutputList->Add(new TH1F("hZvertex", "Z vertex;cm;N_{entries}", 100, -20, 20));
  fOutputList->Add(new TH1F("hNvertexTracks", "N of primary tracks from the primary vertex;cm;N_{entries}", 100, 0, 100));

  fCutsList->Add(new TH1F("hPhoton", "Number of photons after cuts;cut;N_{entries}", 15, 0, 15));
  fCutsList->Add(new TH1F("hLambda", "Number of #Lambda after cuts;cut;N_{entries}", 15, 0, 15));

  // 2D histos
  fOutputList->Add(new TH2F("hTPCResponse", "TPC Response;p [GeV/c];#frac{dE}{dx} [arb.units]", nPtBins, fMinPt, fMaxPt, 250, 0, 250));
  fOutputList->Add(new TH2F("hClusterTOFvsE", "Cluster time vs energy;Time of flight [s];E [GeV]", 200, -fMaxPHOStime, fMaxPHOStime, 50, 0, 2.5));
  fOutputList->Add(new TH2F("hLambdaMass", "#Lambda mass;M_{#Lambda}[GeV/c^{2}];p_{T}[GeV/c]", 50, fLambdaMass - fMaxMassDeviation, fLambdaMass + fMaxMassDeviation, nPtBins, fMinPt, fMaxPtLambda));
  fOutputList->Add(new TH2F("hLambdaBarMass", "#bar{#Lambda} mass;M_{#bar{#Lambda}}[GeV/c^{2}];p_{T}[GeV/c]", 50, fLambdaMass - fMaxMassDeviation, fLambdaMass + fMaxMassDeviation, nPtBins, fMinPt, fMaxPtLambda));

  for (Double_t etaCut : fEtaRangeCuts)
  {
    fOutputList->Add(new TH2F(Form("hDTReGpGp_%g", etaCut), "hDTReGpGp;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPtPhoton));
    fOutputList->Add(new TH2F(Form("hDTReGpGc_%g", etaCut), "hDTReGpGc;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPtPhoton));
    fOutputList->Add(new TH2F(Form("hDTReGcGc_%g", etaCut), "hDTReGcGc;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPtPhoton));
    fOutputList->Add(new TH2F(Form("hDTMiGpGp_%g", etaCut), "hDTMiGpGp;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPtPhoton));

    fOutputList->Add(new TH2F(Form("hDTReLmGp_%g", etaCut), "hDTReLmGp;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fOutputList->Add(new TH2F(Form("hDTMiLmGp_%g", etaCut), "hDTMiLmGp;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fOutputList->Add(new TH2F(Form("hDTM2LmGp_%g", etaCut), "hDTM2LmGp;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fOutputList->Add(new TH2F(Form("hDTReLmGc_%g", etaCut), "hDTReLmGc;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fOutputList->Add(new TH2F(Form("hDTMiLmGc_%g", etaCut), "hDTMiLmGc;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
  }

  // PHOS cuts histos
  for (int i = 0; i < 6; ++i)
  {
    fCutsList->Add(new TH1F(Form("hClusterEnergy_%d", i), "Cluster energy", 100, 0, 5));
    fCutsList->Add(new TH1F(Form("hClusterM02_%d", i), "Cluster assymmetry", 100, -5.0, 5.0));
    fCutsList->Add(new TH1F(Form("hClusterTOF_%d", i), "Cluster Time OF Flight", 100, -fMaxPHOStime, fMaxPHOStime));
    fCutsList->Add(new TH1F(Form("hClusterChi2_%d", i), "Cluster Dispersion", 100, 0, 10));
  }

  // PCM cuts histos
  for (int i = 0; i < 10; ++i)
  {
    fCutsList->Add(new TH1F(Form("hEtaNeg_%d", i), "#eta distribution of conversion electrons", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));
    fCutsList->Add(new TH1F(Form("hEtaPos_%d", i), "#eta distribution of MC conversion positrons", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));
    fCutsList->Add(new TH1F(Form("hnTPCClusNeg_%d", i), "TPC clusters crossed by electron track", 200, 0, 200));
    fCutsList->Add(new TH1F(Form("hnTPCClusPos_%d", i), "TPC clusters crossed by electron track", 200, 0, 200));
    fCutsList->Add(new TH1F(Form("hSigmaElectron_%d", i), "Number of #sigma in TPC for e^{-}", nSigmaBins, -fMaxNSigma, fMaxNSigma));
    fCutsList->Add(new TH1F(Form("hSigmaPositron_%d", i), "Number of #sigma in TPC for e^{+}", nSigmaBins, -fMaxNSigma, fMaxNSigma));
    fCutsList->Add(new TH1F(Form("hV0CPA_%d", i), "", 100, 0.5, 1));
    fCutsList->Add(new TH1F(Form("hV0Radius_%d", i), "", 200, 0, 200));
    fCutsList->Add(new TH1F(Form("hphotonmass_%d", i), "", 100, 0, 1));
    fCutsList->Add(new TH1F(Form("hphotonPt_%d", i), "", 100, 0, fMaxPt));
    fCutsList->Add(new TH2F(Form("hangles_%d", i), "#Delta #theta vs. open angle plot", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi()));
    fCutsList->Add(new TH2F(Form("hArPodPhoton_%d", i), "Armenteros-Podolansky plot for photons", 100, -1, 1, 100, 0, 0.25));
    fCutsList->Add(new TH2F(Form("hPhotonMvsPt_%d", i), "M_{#gamma} vs. p_{T #gamma}", 100, 0, 1, 100, 0, fMaxPt));
  }

  // Lambda cuts histos
  for (int i = 0; i < 9; ++i)
  {
    fCutsList->Add(new TH1F(Form("hLambdaEtaPr_%d", i), "#eta distribution of Proton", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));
    fCutsList->Add(new TH1F(Form("hLambdaEtaPi_%d", i), "#eta distribution of Pion", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));
    fCutsList->Add(new TH1F(Form("hLambdaDCA_%d", i), "#Lambda distance of closest approach", 100, 0, 2));
    fCutsList->Add(new TH1F(Form("hLambdaCPA_%d", i), "#Lambda cosine of pointing angle", 100, 0.9, 1));
    fCutsList->Add(new TH1F(Form("hLambdaV0Raduis_%d", i), "#Lambda vertex radius", 250, 0, 250));
    fCutsList->Add(new TH1F(Form("hLambdaMassCut_%d", i), "#Lambda Mass", 100, 1.0, 1.2));
    fCutsList->Add(new TH1F(Form("hAntiLambdaMass_%d", i), "#bar{#Lambda} Mass", 100, 1.0, 1.2));
    fCutsList->Add(new TH1F(Form("hLambdaPt_%d", i), "#Lambda p_{T}", 100, 0, 10));
    fCutsList->Add(new TH1F(Form("hPosTPCNcls_%d", i), "TPC clusters crossed by positive track", 200, 0, 200));
    fCutsList->Add(new TH1F(Form("hNegTPCNcls_%d", i), "TPC clusters crossed by negative track", 200, 0, 200));

    fCutsList->Add(new TH1F(Form("hSigmaPosProton_%d", i), "Number of #sigma in TPC for p^{+}", nSigmaBins, -fMaxNSigma, fMaxNSigma));
    fCutsList->Add(new TH1F(Form("hSigmaNegProton_%d", i), "Number of #sigma in TPC for p^{-}", nSigmaBins, -fMaxNSigma, fMaxNSigma));
    fCutsList->Add(new TH1F(Form("hSigmaPosPion_%d", i), "Number of #sigma in TPC for #pi^{+}", nSigmaBins, -fMaxNSigma, fMaxNSigma));
    fCutsList->Add(new TH1F(Form("hSigmaNegPion_%d", i), "Number of #sigma in TPC for #pi^{-}", nSigmaBins, -fMaxNSigma, fMaxNSigma));

    fCutsList->Add(new TH2F(Form("hArPodLambda_%d", i), "Armenteros-Podolansky plot for #Lambda", 100, -1, 1, 100, 0, 1));
    fCutsList->Add(new TH1F(Form("hPosDCA_%d", i), "Positive #Lambda daughter DCA", 250, 0, 2));
    fCutsList->Add(new TH1F(Form("hNegDCA_%d", i), "Negative #Lambda daughter DCA", 250, 0, 2));
    fCutsList->Add(new TH1F(Form("hLambdaEta_%d", i), "#Lambda pseudorapidity;#eta;N", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));
  }

  // MC particle cuts histos
  fMCList->Add(new TH1F("hEtaMC", "#eta distribution of MC particles;#eta;N_{entries}", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));

  fMCList->Add(new TH1F("hMCPhotonMass", "Mass distribution of PCM photons;#eta;N_{entries}", nMassBins, fMinPhotonMass, fMaxPhotonMass));
  fMCList->Add(new TH1F("hMCPhotonPt", "p_{T} distribution of PCM photons;#eta;N_{entries}", nPtBins, fMinPt, fMaxPtPhoton));
  fMCList->Add(new TH1F("hMCPhotonEta", "#eta distribution of PCM photons;#eta;N_{entries}", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));

  fMCList->Add(new TH1F("hMCLambdaMass", "Mass distribution of #Lambda;#eta;N_{entries}", nMassBins, fMinLambdaMass, fMaxLambdaMass));
  fMCList->Add(new TH1F("hMCLambdaPt", "p_{T} distribution of #Lambda;#eta;N_{entries}", nPtBins, fMinPt, fMaxPtLambda));
  fMCList->Add(new TH1F("hMCLambdaEta", "#eta distribution of #Lambda;#eta;N_{entries}", 2 * nEtaBins, -2 * fMaxEta, 2 * fMaxEta));

  fMCList->Add(new TH2F("hmc4piLamGam", "hmc4piLamGam;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));

  for (Double_t etaCut : fEtaRangeCuts)
  {
    fMCList->Add(new TH2F(Form("hmc4piLamGam_%g", etaCut), "hmc4piLamGam;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fMCList->Add(new TH2F(Form("hMCReLmGc_%g", etaCut), "hMCReLmGc;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fMCList->Add(new TH2F(Form("hMCReLmGp_%g", etaCut), "hMCReLmGp;M_{#Lambda+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinMass, fMaxMass, nPtBins, fMinPt, fMaxPt));
    fMCList->Add(new TH2F(Form("hMCReGpGc_%g", etaCut), "hMCReGpGc;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPt));
    fMCList->Add(new TH2F(Form("hMCReGcGc_%g", etaCut), "hMCReGcGc;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPt));
    fMCList->Add(new TH2F(Form("hMCReGpGp_%g", etaCut), "hMCReGpGp;M_{#gamma+#gamma}[GeV/c^{2}];p_{T}[GeV/c]", nMassBins, fMinPhotonMass, fMaxPhotonMass, nPtBins, fMinPt, fMaxPt));
  }

  PostData(1, fOutputList);
  PostData(2, fCutsList);
  PostData(3, fMCList);
}

////////////////////////////////////////////////////////////////////////////////
// User-defined event analysis, called for each event
void AliAnalysisTaskSigmaPCMPHOS::UserExec(Option_t *option)
{
  // Main loop, called for each event  Analyze AOD
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr || !mgr->GetInputEventHandler())
    return;

  fAODEvent = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!fAODEvent)
    return;

  Int_t fCurrentRunNumber = InputEvent()->GetRunNumber();
  FillHistogram("hRunNumbers", fCurrentRunNumber);

  Int_t nTracks = fAODEvent->GetNumberOfTracks();
  fPIDResponse = mgr->GetInputEventHandler()->GetPIDResponse();

  // Event selection flags
  Bool_t eventVtxExist = kFALSE;

  // Check the PID response
  if (!fPIDResponse)
    return;

  //Number of Tracks
  if (nTracks == 0)
    return; //No point in continuing if there are no tracks

  fMCEvent = MCEvent(); // Get MC event (called fMCEvent) from the input file

  if (fMCEvent)
  {
    fAODMCTrackArray = dynamic_cast<TClonesArray *>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fAODMCTrackArray)
      return;
    ProcessMC();
  }
  FillV0PhotonArray();

  if (!fPHOSGeo)
    fPHOSGeo = AliPHOSGeometry::GetInstance();

  AliAODVertex *fAODvertex = fAODEvent->GetPrimaryVertex();
  FillHistogram("hNvertexTracks", fAODvertex->GetNContributors());
  if (fAODvertex->GetZ() == 0)
    return;
  FillHistogram("hZvertex", fAODvertex->GetZ());

  for (int i = 0; i < fAODEvent->GetNumberOfTracks(); ++i)
  {
    AliAODTrack *track = dynamic_cast<AliAODTrack *>(fAODEvent->GetTrack(i));
    FillHistogram("hTPCResponse", track->P(), track->GetTPCsignal());
  }

  if (TMath::Abs(fAODvertex->GetZ()) > fMaxZvertex)
    return;

  if (fAODEvent->IsPileupFromSPD())
    return;

  SelectLambda();
  SelectGamma();

  Int_t nGamma = fGamma->GetEntriesFast();
  Int_t nPCM = fPCM->GetEntriesFast();
  Int_t nLambda = fLambda->GetEntriesFast();
  const Double_t fMaxMass = 1.22;

  TLorentzVector pair;

  // fill gamma-gamma  ?repeat without tagged photons
  for (Int_t i = 0; i < nGamma; i++)
  {
    AliCaloPhoton *photonPHOS = dynamic_cast<AliCaloPhoton *>(fGamma->At(i));
    if (!photonPHOS)
      continue;
    for (Int_t j = i + 1; j < nGamma; j++)
    {
      AliCaloPhoton *otherPhotonPHOS = dynamic_cast<AliCaloPhoton *>(fGamma->At(j));
      if (!otherPhotonPHOS)
        continue;

      pair = *photonPHOS + *otherPhotonPHOS;
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(pair.Eta()) < etaCut)
          FillHistogram(Form("hDTReGpGp_%g", etaCut), m, pT);
      }
    }
  }

  //Fill gamma-Lambda
  for (Int_t i = 0; i < nGamma; i++)
  {
    AliCaloPhoton *photonPHOS = dynamic_cast<AliCaloPhoton *>(fGamma->At(i));
    if (!photonPHOS)
      continue;
    for (Int_t j = 0; j < nLambda; j++)
    {
      TLorentzVector *Lambda = dynamic_cast<TLorentzVector *>(fLambda->At(j));
      if (!Lambda)
        continue;

      pair = *photonPHOS + *Lambda;
      Double_t m = pair.M();
      if (m > fMaxMass)
        continue;
      Double_t pT = pair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(pair.Eta()) < etaCut)
          FillHistogram(Form("hDTReLmGp_%g", etaCut), m, pT);
      }
    }
  }

  // Fill gamma_v0-Lambda, idea to separate on Lambdap - particle and Lambdam - antiparticle
  TLorentzVector pvv1, pvv2;
  for (Int_t i = 0; i < nPCM; i++)
  {
    TLorentzVector *photonPCM = dynamic_cast<TLorentzVector *>(fPCM->At(i));
    if (!photonPCM)
      continue;
    for (Int_t j = 0; j < nLambda; j++)
    {
      TLorentzVector *Lambda = dynamic_cast<TLorentzVector *>(fLambda->At(j));
      if (!Lambda)
        continue;
      pair = *photonPCM + *Lambda;
      Double_t pT = pair.Pt();
      Double_t m = pair.M();
      if (m > fMaxMass)
        continue;
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(pair.Eta()) < etaCut)
          FillHistogram(Form("hDTReLmGc_%g", etaCut), m, pT);
      }
    }
    for (Int_t j = 0; j < nGamma; j++)
    {
      AliCaloPhoton *photonPHOS = dynamic_cast<AliCaloPhoton *>(fGamma->At(j));
      if (!photonPHOS)
        continue;

      pair = *photonPCM + *photonPHOS;
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(pair.Eta()) < etaCut)
          FillHistogram(Form("hDTReGpGc_%g", etaCut), m, pT);
      }
    }
    for (Int_t j = i + 1; j < nPCM; j++)
    {

      TLorentzVector *otherPhotonPCM = dynamic_cast<TLorentzVector *>(fPCM->At(i));
      if (!otherPhotonPCM)
        continue;
      pair = *photonPCM + *otherPhotonPCM;
      Double_t m = pair.M();
      Double_t pT = pair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(pair.Eta()) < etaCut)
          FillHistogram(Form("hDTReGcGc_%g", etaCut), m, pT);
      }
    }
  } // end of converted photons

  //================MIXEDs - later, when REAL events will wi checked
  // * Fill gamma-gamma
  for (Int_t m = 0; m < fMixGamma->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixGamma->At(m);
    for (Int_t j = 0; j < tmp->GetEntriesFast(); j++)
    {
      AliCaloPhoton *pv2 = (AliCaloPhoton *)tmp->At(j);
      for (Int_t i = 0; i < nGamma; i++)
      {
        AliCaloPhoton *pv1 = (AliCaloPhoton *)fGamma->At(i);
        pair = *pv1 + *pv2;
        Double_t m = pair.M();
        Double_t pT = pair.Pt();
        for (Double_t etaCut : fEtaRangeCuts)
        {
          if (TMath::Abs(pair.Eta()) < etaCut)
            FillHistogram(Form("hDTMiGpGp_%g", etaCut), m, pT);
        }
      }
    }
  }

  // * Fill gamma-mix.Lambda
  for (Int_t m = 0; m < fMixLambda->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixLambda->At(m);
    for (Int_t j = 0; j < tmp->GetEntriesFast(); j++)
    {
      TLorentzVector *pv2 = (TLorentzVector *)tmp->At(j);
      for (Int_t i = 0; i < nGamma; i++)
      {
        AliCaloPhoton *pv1 = (AliCaloPhoton *)fGamma->At(i);
        pair = *pv1 + *pv2;
        Double_t m = pair.M();
        if (m > fMaxMass)
          continue;
        Double_t pT = pair.Pt();
        for (Double_t etaCut : fEtaRangeCuts)
        {
          if (TMath::Abs(pair.Eta()) < etaCut)
            FillHistogram(Form("hDTM2LmGp_%g", etaCut), m, pT);
        }
      }
    }
  }
  // * Fill mix.gamma-Lambda
  for (Int_t m = 0; m < fMixGamma->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixGamma->At(m);
    for (Int_t i = 0; i < tmp->GetEntriesFast(); i++)
    {
      AliCaloPhoton *pv1 = (AliCaloPhoton *)tmp->At(i);
      for (Int_t j = 0; j < nLambda; j++)
      {
        TLorentzVector *pv2 = (TLorentzVector *)fLambda->At(j);
        pair = *pv1 + *pv2;
        Double_t m = pair.M();
        if (m > fMaxMass)
          continue;
        Double_t pT = pair.Pt();
        for (Double_t etaCut : fEtaRangeCuts)
        {
          if (TMath::Abs(pair.Eta()) < etaCut)
            FillHistogram(Form("hDTMiLmGp_%g", etaCut), m, pT);
        }
      }
    }
  }

  // * Fill gamma.v0-mix.Lambda
  for (Int_t m = 0; m < fMixLambda->GetSize(); m++)
  {
    TClonesArray *tmp = (TClonesArray *)fMixLambda->At(m);
    for (Int_t j = 0; j < tmp->GetEntriesFast(); j++)
    {
      TLorentzVector *pv2 = (TLorentzVector *)tmp->At(j);
      for (Int_t i = 0; i < fConvPhotonArray.size(); i++)
      {
        AliAODv0 *v0_1 = (AliAODv0 *)fAODEvent->GetV0(fConvPhotonArray.at(i));
        if (!v0_1)
          continue;
        pvv1.SetXYZM(v0_1->Px(), v0_1->Py(), v0_1->Pz(), 0);
        pair = pvv1 + *pv2;
        Double_t m = pair.M();
        if (m > fMaxMass)
          continue;
        Double_t pT = pair.Pt();
        for (Double_t etaCut : fEtaRangeCuts)
        {
          if (TMath::Abs(pair.Eta()) < etaCut)
            FillHistogram(Form("hDTMiLmGc_%g", etaCut), m, pT);
        }
      }
    }
  }

  // MC combinations
  for (long unsigned int i = 0; i < fMCPhotonArray.size(); ++i)
  {
    TLorentzVector photon = fMCPhotonArray[i];
    for (long unsigned int j = 0; j < fMCLambdaArray.size(); ++j)
    {
      TLorentzVector Lambda = fMCLambdaArray[j];
      TLorentzVector Sigma0 = Lambda + photon;
      Double_t mass = Sigma0.M();
      if (mass > fMaxMass)
        continue;
      Double_t pT = Sigma0.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(Sigma0.Eta()) < etaCut)
          FillMCHistogram(Form("hMCReLmGc_%g", etaCut), mass, pT);
      }
    }

    for (long unsigned int j = i + 1; j < fMCPhotonArray.size(); ++j)
    {
      TLorentzVector otherPhoton = fMCPhotonArray[j];
      TLorentzVector photonPair = photon + otherPhoton;
      Double_t mass = photonPair.M();
      if (mass > fMaxMass)
        continue;
      Double_t pT = photonPair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(photonPair.Eta()) < etaCut)
          FillMCHistogram(Form("hMCReGcGc_%g", etaCut), mass, pT);
      }
    }

    for (long unsigned int j = 0; j < fMCPHOSArray.size(); ++j)
    {
      TLorentzVector otherPhoton = fMCPHOSArray[i];
      TLorentzVector photonPair = otherPhoton + photon;
      Double_t mass = photonPair.M();
      if (mass > fMaxMass)
        continue;
      Double_t pT = photonPair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(photonPair.Eta()) < etaCut)
          FillMCHistogram(Form("hMCReGpGc_%g", etaCut), mass, pT);
      }
    }
  }

  for (long unsigned int i = 0; i < fMCPHOSArray.size(); ++i)
  {
    TLorentzVector photon = fMCPHOSArray[i];
    for (long unsigned int j = 0; j < fMCLambdaArray.size(); ++j)
    {
      TLorentzVector Lambda = fMCLambdaArray[j];
      TLorentzVector Sigma0 = Lambda + photon;
      Double_t mass = Sigma0.M();
      Double_t pT = Sigma0.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(Sigma0.Eta()) < etaCut)
          FillMCHistogram(Form("hMCReLmGp_%g", etaCut), mass, pT);
      }
    }

    for (long unsigned int j = i + 1; j < fMCPHOSArray.size(); ++j)
    {
      TLorentzVector otherPhoton = fMCPHOSArray[j];
      TLorentzVector photonPair = photon + otherPhoton;
      Double_t mass = photonPair.M();
      Double_t pT = photonPair.Pt();
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(photonPair.Eta()) < etaCut)
          FillMCHistogram(Form("hMCReGpGp_%g", etaCut), mass, pT);
      }
    }
  }

  //Now we either add current events to stack or remove ==> To check, skip for now
  //If no photons in current event - no need to add it to mixed
  const Int_t kMixEvents = 10;
  const Int_t kMixEventsHadr = 100;
  if (fGamma->GetEntriesFast() > 0)
  {
    fMixGamma->AddFirst(fGamma);
    fGamma = 0;
    if (fMixGamma->GetSize() > kMixEvents)
    { //Remove redundant events
      TClonesArray *tmp = dynamic_cast<TClonesArray *>(fMixGamma->Last());
      fMixGamma->RemoveLast();
      delete tmp;
    }
  }

  // make mixed Lambda
  if (fLambda->GetEntriesFast() > 0)
  {
    fMixLambda->AddFirst(fLambda);
    fLambda = 0;
    if (fMixLambda->GetSize() > kMixEventsHadr)
    { //Remove redundant events
      TClonesArray *tmp = dynamic_cast<TClonesArray *>(fMixLambda->Last());
      fMixLambda->RemoveLast();
      delete tmp;
    }
  }

  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fCutsList);
  PostData(3, fMCList);
}

////////////////////////////////////////////////////////////////////////////////
// Process simulated events if Monte-Carlo flag is true
void AliAnalysisTaskSigmaPCMPHOS::ProcessMC()
{
  Int_t nMCTracks = fMCEvent->GetNumberOfTracks();

  for (Int_t iMCtrack = 0; iMCtrack < nMCTracks; iMCtrack++) // index used to be 1
  {
    AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(iMCtrack));
    if (!mcPart)
      continue;
    Int_t MCPartPDGCode = mcPart->GetPdgCode();

    FillMCHistogram("hEtaMC", mcPart->Eta());
    if (TMath::Abs(mcPart->Eta()) > fMaxMCEta)
      continue;

    if (mcPart->GetMother() == -1)
      continue;
    AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(mcPart->GetMother()));
    if (!mother)
      continue;
    Int_t MCMotherPDGCode = mother->GetPdgCode();

    if (MCPartPDGCode == 22 && TMath::Abs(MCMotherPDGCode) == 3212)
    {
      FillMCHistogram("hMCPhotonMass", mcPart->M());
      FillMCHistogram("hMCPhotonPt", mcPart->Pt());
      FillMCHistogram("hMCPhotonEta", mcPart->Eta());
    }

    if (TMath::Abs(MCPartPDGCode) == 3122 && TMath::Abs(MCMotherPDGCode) == 3212)
    {

      FillMCHistogram("hMCLambdaMass", mcPart->M());
      FillMCHistogram("hMCLambdaPt", mcPart->Pt());
      FillMCHistogram("hMCLambdaEta", mcPart->Eta());

      FillMCHistogram("hmc4piLamGam", mother->M(), mother->Pt());
      for (Double_t etaCut : fEtaRangeCuts)
      {
        if (TMath::Abs(mcPart->Eta()) < etaCut)
          FillMCHistogram(Form("hmc4piLamGam_%g", etaCut), mother->M(), mother->Pt());
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Fill array of photons' V0
void AliAnalysisTaskSigmaPCMPHOS::FillV0PhotonArray()
{

  // the same codes as for Sigma+ analysis
  Double_t primaryVtxPos[3] = {0, 0, 0};

  //Clear V0 Photon Array and reset counter
  fConvPhotonArray.clear();
  fMCPhotonArray.clear();

  if (!fPCM)
    fPCM = new TClonesArray("TLorentzVector", 100);
  else
    fPCM->Clear();
  Int_t inPCM = 0;

  Int_t nV0 = fAODEvent->GetNumberOfV0s(); //Number of V0s in the event
  if (nV0 == 0)
    return; //Return if there is no V0 to be processed

  for (Int_t iV0 = 0; iV0 < nV0; iV0++)
  { //Loop over V0s in the event

    //Initialisation of the local Bools
    Bool_t isElectronTPC = kFALSE;
    Bool_t isPositronTPC = kFALSE;
    Bool_t isPhotonTPC = kFALSE;

    TVector3 vecN, vecP, vecM;                 //Momentum Vectors for V0 tracks
    TLorentzVector electron, positron, photon; //Lorentzvectors for invariant mass calculation

    fAODv0 = dynamic_cast<AliAODv0 *>(fAODEvent->GetV0(iV0)); //Get V0 object
    if (!fAODv0)
      continue;
    // Check basic V0 properties: 2 Daughters, opposite charge, total charge = 0
    if (fAODv0->GetNDaughters() != 2)
      continue;
    if (fAODv0->GetNProngs() != 2)
      continue;
    if (fAODv0->GetCharge() != 0)
      continue;
    if (fAODv0->ChargeProng(0) == fAODv0->ChargeProng(1))
      continue;
    // Get daughter tracks
    electronCandidate = dynamic_cast<AliAODTrack *>(fAODv0->GetDaughter(0));
    positronCandidate = dynamic_cast<AliAODTrack *>(fAODv0->GetDaughter(1));
    if (electronCandidate->GetSign() == positronCandidate->GetSign())
      continue;
    if (electronCandidate->Charge() > 0)
    { //Check correct charge
      electronCandidate = dynamic_cast<AliAODTrack *>(fAODv0->GetDaughter(1));
      positronCandidate = dynamic_cast<AliAODTrack *>(fAODv0->GetDaughter(0));
    }

    if (!positronCandidate || !electronCandidate)
      continue;

    if (fMCEvent)
    {
      AliAODMCParticle *MCNegTrack = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(electronCandidate->GetLabel())));
      AliAODMCParticle *MCPosTrack = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(positronCandidate->GetLabel())));
      if (MCPosTrack && MCNegTrack)
      {
        AliAODMCParticle *MCNegCandidate = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCNegTrack->GetMother()));
        AliAODMCParticle *MCPosCandidate = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCPosTrack->GetMother()));
        if (MCPosCandidate && MCNegCandidate && (MCNegCandidate == MCPosCandidate))
        {
          if (MCNegCandidate->GetPdgCode() == 22)
          { // only accept photons as mothers

            AliAODMCParticle *MCNegMother = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCNegCandidate->GetMother()));
            AliAODMCParticle *MCPosMother = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCPosCandidate->GetMother()));
            if (MCNegMother && MCPosMother && (MCNegMother == MCPosMother))
            {
              if (TMath::Abs(MCPosMother->GetPdgCode()) == 3212)
              { //only accept photon from Sigma decay
                TLorentzVector photon(MCNegCandidate->Px(), MCNegCandidate->Py(), MCNegCandidate->Pz(), MCNegCandidate->E());
                fMCPhotonArray.push_back(photon);
                FillCutsHistogram("hPhoton", 14);
              }
            }
          }
        }
      }
    }

    // Check track quality
    Int_t nTPCClustNeg = electronCandidate->GetTPCNcls();
    Int_t nTPCClustPos = positronCandidate->GetTPCNcls();

    // Daughter track PID using TPC
    Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(electronCandidate, AliPID::kElectron);
    Double_t nSigmaTPCpositron = fPIDResponse->NumberOfSigmasTPC(positronCandidate, AliPID::kElectron);
    if (TMath::Abs(nSigmaTPCelectron) < fMaxNsigDaughtTPC)
      isElectronTPC = kTRUE;
    if (TMath::Abs(nSigmaTPCpositron) < fMaxNsigDaughtTPC)
      isPositronTPC = kTRUE;
    if (isElectronTPC && isPositronTPC)
      isPhotonTPC = kTRUE;

    if (!isPhotonTPC)
      continue;
    // Get topological values
    Double_t dcaV0Daughters = TMath::Abs(fAODv0->DcaV0Daughters());
    Double_t dcaPosToPrimVtx = TMath::Abs(fAODv0->DcaPosToPrimVertex());
    Double_t dcaNegToPrimVtx = TMath::Abs(fAODv0->DcaNegToPrimVertex());
    Double_t V0CPA = fAODv0->CosPointingAngle(primaryVtxPos);
    Double_t vtxPosV0[3];
    vtxPosV0[0] = fAODv0->DecayVertexV0X();
    vtxPosV0[1] = fAODv0->DecayVertexV0Y();
    vtxPosV0[2] = fAODv0->DecayVertexV0Z();
    Double_t V0Radius = TMath::Sqrt(vtxPosV0[0] * vtxPosV0[0] + vtxPosV0[1] * vtxPosV0[1]);

    //Get reconstructed cartesian momentum
    vecN.SetXYZ(fAODv0->MomNegX(), fAODv0->MomNegY(), fAODv0->MomNegZ()); //negative daughter
    vecP.SetXYZ(fAODv0->MomPosX(), fAODv0->MomPosY(), fAODv0->MomPosZ()); //positive daughter
    vecM.SetXYZ(fAODv0->MomV0X(), fAODv0->MomV0Y(), fAODv0->MomV0Z());    //mother

    //Custom Armenteros Podolanski calculation since V0 member functions are not reliable!
    Double_t pLNeg = vecN.Dot(vecM) / vecM.Mag(); //Momentum longitudinal
    Double_t pLPos = vecP.Dot(vecM) / vecM.Mag(); //to V0 momentum
    Double_t alpha = (pLPos - pLNeg) / (pLPos + pLNeg);
    Double_t qt = vecN.Perp(vecM);

    // Get kinematic values
    Double_t thetaPos = positronCandidate->Theta();
    Double_t thetaNeg = electronCandidate->Theta();
    Double_t openangle = fAODv0->OpenAngleV0();

    //  Reconstruct photon with TLorentzVector
    Double_t fElectronMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    electron.SetXYZM(vecN(0), vecN(1), vecN(2), fElectronMass);
    positron.SetXYZM(vecP(0), vecP(1), vecP(2), fElectronMass);
    photon = electron + positron;

    // Calculate photon invariant mass with TL
    Double_t photonmass = photon.M();

    // Inv. Mass Cut
    if (photonmass > fMaxphotonMass)
      continue;
    if (photon.Pt() < fMinPhotonPt || photon.Pt() > fMaxPhotonPt)
      continue;

    Double_t photonPt = photon.Pt();
    // Angle calculation
    Double_t deltatheta = thetaPos - thetaNeg;

    // Distributions by cuts parameters

    Double_t params[8] = {openangle, deltatheta, V0CPA, V0Radius, photonmass, alpha, qt, photonPt};
    FillPCMCutHistograms(0, params);
    // Acceptance Cut
    if (TMath::Abs(electronCandidate->Eta()) > fMaxDaughtEta)
      continue;
    if (TMath::Abs(positronCandidate->Eta()) > fMaxDaughtEta)
      continue;
    FillPCMCutHistograms(1, params);
    // Check Track quality and reject poor qualty tracks
    if (nTPCClustNeg < fMinTPCClustDaught)
      continue;
    if (nTPCClustPos < fMinTPCClustDaught)
      continue;
    FillPCMCutHistograms(2, params);
    // Armenteros-Podolanski Cuts
    if (TMath::Abs(alpha) > fMaxAlphaPhoton)
      continue;
    FillPCMCutHistograms(3, params);
    if (TMath::Abs(qt) > fMaxQtPhoton)
      continue;
    FillPCMCutHistograms(4, params);
    // Angle Cut
    if (TMath::Abs(openangle) > fMaxopenangle)
      continue;
    FillPCMCutHistograms(5, params);
    if (TMath::Abs(deltatheta) > fMaxdeltatheta)
      continue;
    FillPCMCutHistograms(6, params);
    // CPA Cut
    if (V0CPA < fMinV0CPA)
      continue;
    FillPCMCutHistograms(7, params);

    //Radius Cut
    if (V0Radius < fMinV0Radius || V0Radius > fMaxV0Radius)
      continue;
    FillPCMCutHistograms(8, params);

    Double_t electronDCA = fAODv0->DcaNegToPrimVertex();
    Double_t positronDCA = fAODv0->DcaPosToPrimVertex();
    if (electronDCA < fMinPhotonDaughterDCA || positronDCA < fMinPhotonDaughterDCA)
      continue;

    Double_t photonDCA = fAODv0->DcaV0Daughters();
    if (photonDCA > fMaxPhotonDaughterTracksDCA)
      continue;

    FillPCMCutHistograms(9, params);
    fConvPhotonArray.push_back(iV0);

    if (inPCM >= fPCM->GetSize())
      fPCM->Expand(inPCM + 50);
    TLorentzVector *photonPCM = new ((*fPCM)[inPCM++]) TLorentzVector(photon.Px(), photon.Py(), photon.Pz(), photon.E());

    //printf("Photon\n");
  } //End of V0 Loop

  //return;
} //End of FillV0PhotonArray()

////////////////////////////////////////////////////////////////////////////////
// Select Lambda hyperons based on cuts specified in the header file
void AliAnalysisTaskSigmaPCMPHOS::SelectLambda()
{
  Double_t fLambdaMass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  fMCLambdaArray.clear();

  if (!fLambda)
    fLambda = new TClonesArray("TLorentzVector", 100);
  else
    fLambda->Clear();

  Int_t nv0 = fAODEvent->GetNumberOfV0s();
  Int_t inLambda = 0;

  const Double_t massLambda = 1.115683;

  while (nv0--)
  {
    fV0 = fAODEvent->GetV0(nv0);
    if (!fV0)
    {
      continue;
    }

    //Use onfly only
    if (!fV0->GetOnFlyStatus())
      continue;

    ntrack1 = (AliAODTrack *)fV0->GetDaughter(1);
    //    if (!AcceptTrack(ntrack1))
    //  continue;

    ptrack1 = (AliAODTrack *)fV0->GetDaughter(0);
    //if (!AcceptTrack(ptrack1))
    //  continue;

    // Remove like-sign
    if (ntrack1->Charge() == ptrack1->Charge())
    {
      continue;
    }

    if (fMCEvent)
    {
      AliAODMCParticle *MCNegTrack = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(ntrack1->GetLabel())));
      AliAODMCParticle *MCPosTrack = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(TMath::Abs(ptrack1->GetLabel())));

      if (MCNegTrack && MCNegTrack)
      {
        AliAODMCParticle *MCNegCandidate = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCNegTrack->GetMother()));
        AliAODMCParticle *MCPosCandidate = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCPosTrack->GetMother()));

        if (MCNegCandidate && MCPosCandidate && (MCNegCandidate == MCPosCandidate) && (TMath::Abs(MCNegCandidate->GetPdgCode()) == 3122))
        {
          AliAODMCParticle *MCNegMother = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCNegCandidate->GetMother()));
          AliAODMCParticle *MCPosMother = dynamic_cast<AliAODMCParticle *>(fAODMCTrackArray->At(MCPosCandidate->GetMother()));

          if (MCNegMother && MCPosMother && (MCNegMother == MCPosMother))
          {
            if (TMath::Abs(MCPosMother->GetPdgCode()) == 3212)
            {
              TLorentzVector Lambda(MCNegCandidate->Px(), MCNegCandidate->Py(), MCNegCandidate->Pz(), MCNegCandidate->M());
              fMCLambdaArray.push_back(Lambda);
              FillCutsHistogram("hLambda", 14);
            }
          }
        }
      }
    }

    if (fV0->Pt() == 0)
    {
      continue;
    }

    if (ntrack1->GetKinkIndex(0) > 0 || ptrack1->GetKinkIndex(0) > 0)
      continue;
    FillLambdaCutHistograms(0);

    // Acceptance Cut
    if (TMath::Abs(ntrack1->Eta()) > fMaxDaughtEta)
      continue;
    if (TMath::Abs(ptrack1->Eta()) > fMaxDaughtEta)
      continue;

    FillLambdaCutHistograms(1);

    Double_t nSigmaPosPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kPion));
    Double_t nSigmaNegPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kPion));
    Double_t nSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kProton));
    Double_t nSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kProton));

    Double_t lV0Position[3];
    fV0->GetXYZ(lV0Position);
    Double_t lV0Radius = TMath::Sqrt(lV0Position[0] * lV0Position[0] + lV0Position[1] * lV0Position[1]);

    if (lV0Radius > fMaxlV0Radius)
      continue;
    FillLambdaCutHistograms(2);

    //DCA V0 daughters
    Double_t dca = fV0->DcaV0Daughters();
    if (dca > fMaxLambdaDaughterTracksDCA)
      continue;
    FillLambdaCutHistograms(3);

    Double_t cpa = fV0->CosPointingAngle(fAODEvent->GetPrimaryVertex());

    if (cpa < fMinCPA)
      continue;
    FillLambdaCutHistograms(4);

    if (fV0->Pt() < fMinLambdaPt || fV0->Pt() > fMaxLambdaPt)
      continue;
    TVector3 vecN, vecP, vecM; //Momentum Vectors for V0 tracks

    //Get reconstructed cartesian momentum
    vecN.SetXYZ(fV0->MomNegX(), fV0->MomNegY(), fV0->MomNegZ()); //negative daughter
    vecP.SetXYZ(fV0->MomPosX(), fV0->MomPosY(), fV0->MomPosZ()); //positive daughter
    vecM.SetXYZ(fV0->MomV0X(), fV0->MomV0Y(), fV0->MomV0Z());    //mother

    //Custom Armenteros Podolanski calculation since V0 member functions are not reliable!
    Double_t pLNeg = vecN.Dot(vecM) / vecM.Mag(); //Momentum longitudinal
    Double_t pLPos = vecP.Dot(vecM) / vecM.Mag(); //to V0 momentum
    Double_t alpha = (pLPos - pLNeg) / (pLPos + pLNeg);
    Double_t qt = vecN.Perp(vecM);

    if (TMath::Abs(alpha) > fMaxAlphaLambda || TMath::Abs(alpha) < fMinAlphaLambda)
      continue;
    if (TMath::Abs(qt) > fMaxQtLambda || TMath::Abs(qt) < fMinQtLambda)
      continue;

    Bool_t isLambda = 0;
    Bool_t isLambdaBar = 0;

    Float_t xyNeg = fV0->DcaNegToPrimVertex();
    Float_t xyPos = fV0->DcaPosToPrimVertex();

    if (TMath::Abs(xyPos) < fMinLambdaDaughterDCA)
      continue;
    if (TMath::Abs(xyNeg) < fMinLambdaDaughterDCA)
      continue;
    FillLambdaCutHistograms(5);

    if (ntrack1->Charge() > 0)
    { //Lambda and proton
      isLambda = (nSigmaNegProton < fMaxSigmaNegProton) && (nSigmaPosPion < fMaxSigmaPosPion);
      isLambdaBar = (nSigmaPosProton < fMaxSigmaPosProton) && (nSigmaNegPion < fMaxSigmaNegPion);
    }
    else
    {
      isLambda = (nSigmaPosProton < fMaxSigmaNegProton) && (nSigmaNegPion < fMaxSigmaPosPion);
      isLambdaBar = (nSigmaNegProton < fMaxSigmaPosProton) && (nSigmaPosPion < fMaxSigmaNegPion);
    }

    if (isLambda && TMath::Abs(fV0->MassLambda() - massLambda) > fMaxMassDeviation)
      continue;
    FillLambdaCutHistograms(6);

    if (isLambdaBar && TMath::Abs(fV0->MassAntiLambda() - massLambda) > fMaxMassDeviation)
      continue;
    FillLambdaCutHistograms(7);

    if (isLambda)
      FillHistogram("hLambdaMass", fV0->MassLambda(), fV0->Pt());
    if (isLambdaBar)
      FillHistogram("hLambdaBarMass", fV0->MassAntiLambda(), fV0->Pt());

    FillLambdaCutHistograms(8);

    //So far combine Lambda and AntiLambda

    if (isLambda || isLambdaBar)
    {
      new ((*fLambda)[inLambda++]) TLorentzVector(fV0->Px(), fV0->Py(), fV0->Pz(), TMath::Sqrt(massLambda * massLambda + fV0->P() * fV0->P()));
      fMCLambdaArray.push_back(TLorentzVector(fV0->Px(), fV0->Py(), fV0->Pz(), TMath::Sqrt(massLambda * massLambda + fV0->P() * fV0->P())));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Reject or accept track
Bool_t AliAnalysisTaskSigmaPCMPHOS::AcceptTrack(const AliAODTrack *t)
{
  if (!t->IsOn(AliAODTrack::kTPCrefit))
    return kFALSE;
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2, 1);
  if (nCrossedRowsTPC < fMinCrossedRowsTPC)
    return kFALSE;
  Int_t findable = t->GetTPCNclsF();
  if (findable <= 0)
    return kFALSE;
  if (nCrossedRowsTPC / findable < fMinTPCCLusterRatio)
    return kFALSE;

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
// Select photons based on cuts specified in the header file
void AliAnalysisTaskSigmaPCMPHOS::SelectGamma()
{

  //Select gamma in PHOS
  Int_t inPHOS = 0;
  if (fGamma)
    fGamma->Clear();
  else
    fGamma = new TClonesArray("AliCaloPhoton", 100);

  const AliAODVertex *fAODv0PHOS = fAODEvent->GetPrimaryVertex();

  Double_t vertexPHOS[3] = {fAODv0PHOS->GetX(), fAODv0PHOS->GetY(), fAODv0PHOS->GetZ()};

  Int_t multClust = fAODEvent->GetNumberOfCaloClusters();

  for (Int_t i = 0; i < multClust; i++)
  {
    AliAODCaloCluster *clu = fAODEvent->GetCaloCluster(i);

    AliCaloPhoton photonCandidate;
    clu->GetMomentum(photonCandidate, vertexPHOS);

    if (fMCEvent)
    {
      TLorentzVector photon(photonCandidate.Px(), photonCandidate.Py(), photonCandidate.Pz(), photonCandidate.E());
      fMCPHOSArray.push_back(photon);
    }

    FillPHOSCutHistograms(clu, 0);
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;
    FillPHOSCutHistograms(clu, 1);
    if (clu->E() > fMaxClusterEnergy || clu->E() < fMinClusterEnergy)
      continue;
    FillPHOSCutHistograms(clu, 2);
    if (clu->GetM02() > fMaxClusterAssymmetry || clu->GetM02() < fMinClusterAssymmetry)
      continue;
    FillPHOSCutHistograms(clu, 3);

    if (TMath::Abs(clu->GetTOF()) > fMaxClusterTOF)
      continue; // TOF cut by D.Peresunko
    FillPHOSCutHistograms(clu, 4);

    //    FillHistogram("hClusterTOFvsE", clu->GetTOF(), clu->E());

    if (!(clu->GetDispersion() > 0.6 && clu->GetDispersion() < 4))
      continue;

    if (TMath::Abs(clu->GetTOF()) > fMaxClusterTOF)
      continue; // TOF cut by D.Peresunko
    FillPHOSCutHistograms(clu, 5);

    FillHistogram("hClusterTOFvsE", clu->GetTOF(), clu->E());

    if (inPHOS >= fGamma->GetSize())
      fGamma->Expand(inPHOS + 50);
    AliCaloPhoton *p = new ((*fGamma)[inPHOS++]) AliCaloPhoton(photonCandidate.Px(), photonCandidate.Py(), photonCandidate.Pz(), photonCandidate.E());
  }
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 1D Cuts histogram
void AliAnalysisTaskSigmaPCMPHOS::FillCutsHistogram(const char *key, Double_t x)
{
  //FillHistogram
  TH1 *hist = dynamic_cast<TH1 *>(fCutsList->FindObject(key));
  if (hist)
    hist->Fill(x);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 2D Cuts histogram
void AliAnalysisTaskSigmaPCMPHOS::FillCutsHistogram(const char *key, Double_t x, Double_t y)
{
  //FillHistogram
  TH2 *th2 = dynamic_cast<TH2 *>(fCutsList->FindObject(key));
  if (th2)
    th2->Fill(x, y);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 1D MC histogram
void AliAnalysisTaskSigmaPCMPHOS::FillMCHistogram(const char *key, Double_t x)
{
  TH1 *hist = dynamic_cast<TH1 *>(fMCList->FindObject(key));
  if (hist)
    hist->Fill(x);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 2D MC histogram
void AliAnalysisTaskSigmaPCMPHOS::FillMCHistogram(const char *key, Double_t x, Double_t y)
{
  TH2 *th2 = dynamic_cast<TH2 *>(fMCList->FindObject(key));
  if (th2)
    th2->Fill(x, y);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 1D histogram
void AliAnalysisTaskSigmaPCMPHOS::FillHistogram(const char *key, Double_t x)
{
  TH1 *hist = dynamic_cast<TH1 *>(fOutputList->FindObject(key));
  if (hist)
    hist->Fill(x);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 2D histogram
void AliAnalysisTaskSigmaPCMPHOS::FillHistogram(const char *key, Double_t x, Double_t y)
{
  TH2 *th2 = dynamic_cast<TH2 *>(fOutputList->FindObject(key));
  if (th2)
    th2->Fill(x, y);
}

////////////////////////////////////////////////////////////////////////////////
// Template to fill 3D histogram
void AliAnalysisTaskSigmaPCMPHOS::FillHistogram(const char *key, Double_t x, Double_t y, Double_t z)
{
  TH3 *th3 = dynamic_cast<TH3 *>(fOutputList->FindObject(key));
  if (th3)
    th3->Fill(x, y, z);
}

////////////////////////////////////////////////////////////////////////////////
// Class Termination
void AliAnalysisTaskSigmaPCMPHOS::Terminate(Option_t *option){};

////////////////////////////////////////////////////////////////////////////////
// Fill PHOS cuts histograms
void AliAnalysisTaskSigmaPCMPHOS::FillPHOSCutHistograms(const AliAODCaloCluster *clu, const Int_t i)
{
  FillCutsHistogram(Form("hClusterEnergy_%d", i), clu->E());
  FillCutsHistogram(Form("hClusterM02_%d", i), clu->GetM02());
  FillCutsHistogram(Form("hClusterTOF_%d", i), clu->GetTOF());
  FillCutsHistogram(Form("hClusterChi2_%d", i), clu->GetDispersion());
}

////////////////////////////////////////////////////////////////////////////////
// Fill PCM cuts histograms
void AliAnalysisTaskSigmaPCMPHOS::FillPCMCutHistograms(const Int_t i, const Double_t params[7])
{
  // Daughter track PID using TPC
  Double_t nSigmaTPCelectron = fPIDResponse->NumberOfSigmasTPC(electronCandidate, AliPID::kElectron);
  Double_t nSigmaTPCpositron = fPIDResponse->NumberOfSigmasTPC(positronCandidate, AliPID::kElectron);
  if (TMath::Abs(nSigmaTPCelectron) < fMaxNsigDaughtTPC)
    FillCutsHistogram(Form("hSigmaElectron_%d", i), nSigmaTPCelectron);
  if (TMath::Abs(nSigmaTPCpositron) < fMaxNsigDaughtTPC)
    FillCutsHistogram(Form("hSigmaPositron_%d", i), nSigmaTPCpositron);

  FillCutsHistogram(Form("hEtaNeg_%d", i), electronCandidate->Eta());
  FillCutsHistogram(Form("hEtaPos_%d", i), positronCandidate->Eta());
  FillCutsHistogram(Form("hnTPCClusNeg_%d", i), electronCandidate->GetTPCNcls());
  FillCutsHistogram(Form("hnTPCClusPos_%d", i), positronCandidate->GetTPCNcls());
  FillCutsHistogram(Form("hangles_%d", i), params[1], params[0]);
  FillCutsHistogram(Form("hV0CPA_%d", i), params[2]);
  FillCutsHistogram(Form("hV0Radius_%d", i), params[3]);
  FillCutsHistogram(Form("hphotonmass_%d", i), params[4]);
  FillCutsHistogram(Form("hArPodPhoton_%d", i), params[5], params[6]);
  FillCutsHistogram(Form("hphotonPt_%d", i), params[7]);
  FillCutsHistogram(Form("hPhotonMvsPt_%d", i), params[4], params[7]);
  FillCutsHistogram("hPhoton", i);
}

////////////////////////////////////////////////////////////////////////////////
// Fill Lambda cuts histograms
void AliAnalysisTaskSigmaPCMPHOS::FillLambdaCutHistograms(const Int_t i)
{

  Bool_t isLambda = 0;
  Bool_t isLambdaBar = 0;

  Double_t lV0Position[3];
  fV0->GetXYZ(lV0Position);
  Double_t lV0Radius = TMath::Sqrt(lV0Position[0] * lV0Position[0] + lV0Position[1] * lV0Position[1]);

  TVector3 vecN, vecP, vecM; //Momentum Vectors for V0 tracks

  //Get reconstructed cartesian momentum
  vecN.SetXYZ(fV0->MomNegX(), fV0->MomNegY(), fV0->MomNegZ()); //negative daughter
  vecP.SetXYZ(fV0->MomPosX(), fV0->MomPosY(), fV0->MomPosZ()); //positive daughter
  vecM.SetXYZ(fV0->MomV0X(), fV0->MomV0Y(), fV0->MomV0Z());    //mother

  //Custom Armenteros Podolanski calculation since V0 member functions are not reliable!
  Double_t pLNeg = vecN.Dot(vecM) / vecM.Mag(); //Momentum longitudinal
  Double_t pLPos = vecP.Dot(vecM) / vecM.Mag(); //to V0 momentum
  Double_t alpha = (pLPos - pLNeg) / (pLPos + pLNeg);
  Double_t qt = vecN.Perp(vecM);

  Double_t nSigmaPosPion = fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kPion);
  Double_t nSigmaNegPion = fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kPion);
  Double_t nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(ptrack1, AliPID::kProton);
  Double_t nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(ntrack1, AliPID::kProton);

  Bool_t isProton = (TMath::Abs(nSigmaPosProton) < fMaxSigmaPosProton);
  Bool_t isAntiProton = (TMath::Abs(nSigmaNegProton) < fMaxSigmaNegProton);
  Bool_t isPosPion = (TMath::Abs(nSigmaPosPion) < fMaxSigmaPosPion);
  Bool_t isNegPion = (TMath::Abs(nSigmaNegPion) < fMaxSigmaNegPion);

  isLambda = isProton && isNegPion;
  isLambdaBar = isAntiProton && isPosPion;

  TLorentzVector positiveVector(ptrack1->Px(), ptrack1->Py(), ptrack1->Pz(), ptrack1->E());
  TLorentzVector negativeVector(ntrack1->Px(), ntrack1->Py(), ntrack1->Pz(), ntrack1->E());
  TLorentzVector Lambda = positiveVector + negativeVector;

  if (isLambda)
  {
    FillCutsHistogram(Form("hLambdaEtaPr_%d", i), ptrack1->Eta());
    FillCutsHistogram(Form("hLambdaEtaPi_%d", i), ntrack1->Eta());
    FillCutsHistogram(Form("hSigmaPosProton_%d", i), nSigmaPosProton);
    FillCutsHistogram(Form("hSigmaNegPion_%d", i), nSigmaNegPion);
    FillCutsHistogram(Form("hLambdaMassCut_%d", i), Lambda.M());
  }

  if (isLambdaBar)
  {
    FillCutsHistogram(Form("hLambdaEtaPr_%d", i), ntrack1->Eta());
    FillCutsHistogram(Form("hLambdaEtaPi_%d", i), ptrack1->Eta());
    FillCutsHistogram(Form("hSigmaNegProton_%d", i), nSigmaNegProton);
    FillCutsHistogram(Form("hSigmaPosPion_%d", i), nSigmaPosPion);
    FillCutsHistogram(Form("hAntiLambdaMass_%d", i), Lambda.M());
  }
  FillCutsHistogram("hLambda", i);
  FillCutsHistogram(Form("hPosTPCNcls_%d", i), ptrack1->GetTPCNcls());
  FillCutsHistogram(Form("hNegTPCNcls_%d", i), ntrack1->GetTPCNcls());
  FillCutsHistogram(Form("hLambdaDCA_%d", i), fV0->DcaV0Daughters());
  FillCutsHistogram(Form("hLambdaCPA_%d", i), fV0->CosPointingAngle(fAODEvent->GetPrimaryVertex()));
  FillCutsHistogram(Form("hLambdaV0Raduis_%d", i), lV0Radius);
  FillCutsHistogram(Form("hLambdaPt_%d", i), Lambda.Pt());
  FillCutsHistogram(Form("hPosDCA_%d", i), fV0->DcaPosToPrimVertex());
  FillCutsHistogram(Form("hNegDCA_%d", i), fV0->DcaNegToPrimVertex());
  FillCutsHistogram(Form("hLambdaEta_%d", i), fV0->Eta());
  FillCutsHistogram(Form("hArPodLambda_%d", i), alpha, qt);
}
