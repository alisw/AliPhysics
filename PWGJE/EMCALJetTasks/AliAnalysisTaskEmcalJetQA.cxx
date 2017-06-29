/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <cstring>

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <THashList.h>

#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliTrackContainer.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"

#include "AliAnalysisTaskEmcalJetQA.h"

ClassImp(AliAnalysisTaskEmcalJetQA)

//________________________________________________________________________
AliAnalysisTaskEmcalJetQA::AliAnalysisTaskEmcalJetQA() :
  AliAnalysisTaskEmcalLight("AliAnalysisTaskEmcalJetQA", kTRUE),
  fCellEnergyCut(0.05),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fDoLeadingObjectPosition(0),
  fMaxCellsInCluster(50),
  fPtBinWidth(0.5),
  fMaxPt(150),
  fSeparateEMCalDCal(kTRUE),
  fIsEmbedded(kFALSE),
  fCent2(0),
  fCent3(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0),
  fNTotTracks(0),
  fLeadingTrack(),
  fHistManager("AliAnalysisTaskEmcalJetQA")
{
  // Default constructor.

  memset(fNTotClusters, 0, sizeof(Int_t)*3);

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetQA::AliAnalysisTaskEmcalJetQA(const char *name) :
  AliAnalysisTaskEmcalLight(name, kTRUE),
  fCellEnergyCut(0.05),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fDoLeadingObjectPosition(0),
  fMaxCellsInCluster(50),
  fPtBinWidth(0.5),
  fMaxPt(150),
  fSeparateEMCalDCal(kTRUE),
  fIsEmbedded(kFALSE),
  fCent2(0),
  fCent3(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0),
  fNTotTracks(0),
  fLeadingTrack(),
  fHistManager(name)
{
  // Standard

  memset(fNTotClusters, 0, sizeof(Int_t)*3);

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetQA::~AliAnalysisTaskEmcalJetQA()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetQA::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  TString histname;
  TString title;

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  for (auto cont_it : fParticleCollArray) {
    AliParticleContainer* cont = cont_it.second;
    if (!fParticleLevel && fIsMC) {
      for (Int_t i = 0; i < GetNCentBins(); i++) {
        histname = TString::Format("%s/fHistTrNegativeLabels_%d", cont->GetArrayName().Data(), i);
        title = histname + ";% of negative labels;counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), 500, 0, 1);

        histname = TString::Format("%s/fHistTrZeroLabels_%d", cont->GetArrayName().Data(), i);
        title = histname + ";% of zero labels;counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), 500, 0, 1);
      }
    }

    Int_t nlabels = 4;
    if (fParticleLevel) nlabels = 1;

    for (Int_t i = 0; i < GetNCentBins(); i++) {
      histname = TString::Format("%s/fHistRejectionReason_%d", cont->GetArrayName().Data(), i);
      title = histname + ";Rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 40, 0, 100);
      SetRejectionReasonLabels(hist->GetXaxis());

      for (Int_t j = 0; j < nlabels; j++) {
        histname = TString::Format("%s/fHistTrPhiEtaPt_%d_%d", cont->GetArrayName().Data(), i, j);
        title = histname + ";#eta;#phi;#it{p}_{T} (GeV/#it{c})";
        fHistManager.CreateTH3(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi(), nPtBins, 0, fMaxPt);
      }

      if (!fParticleLevel) {
        if (fIsMC) {
          histname = TString::Format("%s/fHistTrPhiEtaZeroLab_%d", cont->GetArrayName().Data(), i);
          title = histname + ";#eta;#phi;counts";
          fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi());

          histname = TString::Format("%s/fHistTrPtZeroLab_%d", cont->GetArrayName().Data(), i);
          title = histname + ";#it{p}_{T} (GeV/#it{c});counts";
          fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
        }

        histname = TString::Format("%s/fHistTrEmcPhiEta_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#eta;#phi;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi());

        histname = TString::Format("%s/fHistTrEmcPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

        histname = TString::Format("%s/fHistTrPhiEtaNonProp_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#eta;#phi;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi());

        histname = TString::Format("%s/fHistTrPtNonProp_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

        histname = TString::Format("%s/fHistDeltaEtaPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});#delta#eta;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});#delta#phi;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});#delta#it{p}_{T} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, -fMaxPt/2, fMaxPt/2);
      }
    }
  }

  for (auto cont_it : fClusterCollArray) {
    AliClusterContainer* cont = cont_it.second;
    for (Int_t i = 0; i < GetNCentBins(); i++) {
      const Int_t nSM = 20;

      for (Int_t sm = 0; sm < nSM; sm++) {
        histname = TString::Format("%s/BySM/fHistClusEnergy_SM%d_%d", cont->GetArrayName().Data(), sm, i);
        title = histname + ";#it{E}_{cluster} (GeV);counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
      }

      histname = TString::Format("%s/fHistRejectionReason_%d", cont->GetArrayName().Data(), i);
      title = histname + ";Rejection reason;#it{E}_{cluster} (GeV);counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 40, 0, 100);
      SetRejectionReasonLabels(hist->GetXaxis());

      histname = TString::Format("%s/fHistClusPosition_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{x} (cm);#it{y} (cm);#it{z} (cm)";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 50, -500, 500, 50, -500, 500, 50, -500, 500);

      histname = TString::Format("%s/fHistClusPhiEtaEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#eta;#phi;#it{E}_{cluster} (GeV)";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi(), nPtBins, 0, fMaxPt);

      if (fForceBeamType != kpp) {
        histname = TString::Format("%s/fHistClusDeltaPhiEPEnergy_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{E}_{cluster} (GeV);#phi_{cluster} - #psi_{EP};counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 100, 0, TMath::Pi());
      }

      if (fIsEmbedded) {
        histname = TString::Format("%s/fHistClusMCEnergyFraction_%d", cont->GetArrayName().Data(), i);
        title = histname + ";MC fraction;counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, 1.2);
      }

      histname = TString::Format("%s/fHistClusTimeEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);time (s);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, -5e-6, 5e-6);

      Int_t nbins = fMaxCellsInCluster;
      while (nbins > nPtBins) nbins /= 2;
      histname = TString::Format("%s/fHistNCellsEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);#it{N}_{cells};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nbins, -0.5, fMaxCellsInCluster-0.5);

      histname = TString::Format("%s/fHistFcrossEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);#it{F}_{cross};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 200, -3.5, 1.5);

      histname = TString::Format("%s/fHistClusEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

      histname = TString::Format("%s/fHistClusNonLinCorrEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

      histname = TString::Format("%s/fHistClusHadCorrEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
    }
  }

  if (!fCaloCellsName.IsNull()) {
    for (Int_t i = 0; i < GetNCentBins(); i++) {
      histname = TString::Format("%s/fHistCellsAbsIdEnergy_%d", fCaloCellsName.Data(), i);
      title = histname + ";cell abs. Id;#it{E}_{cell} (GeV);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 20000,0,20000,(Int_t)(nPtBins / 2), 0, fMaxPt / 2);

      histname = TString::Format("%s/fHistCellsAbsIdTime_%d", fCaloCellsName.Data(), i);
      title = histname + ";cell abs. Id;#it{time}_{cell} (s);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 20000,0,20000,(Int_t)(nPtBins / 2), -5e-6, 5e-6);
    }
  }

  Int_t dim = 0;
  TString axistitle[40];
  Int_t nbins[40] = {0};
  Double_t min[40] = {0};
  Double_t max[40] = {0};

  if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
    axistitle[dim] = "Centrality %";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 100;
    dim++;

    if (!fCentMethod2.IsNull()) {
      axistitle[dim] = Form("Centrality %s %%", fCentMethod2.Data());
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = 100;
      dim++;
    }

    if (!fCentMethod3.IsNull()) {
      axistitle[dim] = Form("Centrality %s %%", fCentMethod3.Data());
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = 100;
      dim++;
    }

    if (fDoV0QA==1) {
      axistitle[dim] = "V0A total multiplicity";
      nbins[dim] = 200;
      min[dim] = 0;
      max[dim] = 20000;
      dim++;

      axistitle[dim] = "V0C total multiplicity";
      nbins[dim] = 200;
      min[dim] = 0;
      max[dim] = 20000;
      dim++;
    }
    else if (fDoV0QA==2) {
      axistitle[dim] = "V0A+V0C total multiplicity";
      nbins[dim] = 300;
      min[dim] = 0;
      max[dim] = 30000;
      dim++;
    }

    if (fDoEPQA) {
      axistitle[dim] = "#psi_{EP}";
      nbins[dim] = 200;
      min[dim] = -TMath::Pi();
      max[dim] = TMath::Pi();
      dim++;
    }
  }

  if (fParticleCollArray.size() > 0) {
    axistitle[dim] = "No. of tracks";
    if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
      nbins[dim] = 3000;
      min[dim] = -0.5;
      max[dim] = 6000-0.5;
    }
    else {
      nbins[dim] = 200;
      min[dim] = 0;
      max[dim] = 200;
    }
    dim++;

    axistitle[dim] = "#it{p}_{T,track}^{leading} (GeV/c)";
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;

    if (fDoLeadingObjectPosition) {
      axistitle[dim] = "#eta_{track}^{leading}";
      nbins[dim] = 100;
      min[dim] = -1;
      max[dim] = 1;
      dim++;

      axistitle[dim] = "#phi_{track}^{leading}";
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = TMath::TwoPi();
      dim++;
    }
  }

  if (fClusterCollArray.size() > 0) {
    axistitle[dim] = "No. of clusters";

    if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
      nbins[dim] = 2000;
      min[dim] = -0.5;
      max[dim] = 4000-0.5;
    }
    else {
      nbins[dim] = 200;
      min[dim] = 0;
      max[dim] = 200;
    }
    dim++;

    if (fSeparateEMCalDCal) {
      axistitle[dim] = "#it{E}_{EMCal cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;

      axistitle[dim] = "#it{E}_{DCal cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;

      axistitle[dim] = "#it{E}_{PHOS cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;

      if (fDoLeadingObjectPosition) {
        axistitle[dim] = "#eta_{EMCal cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = -1;
        max[dim] = 1;
        dim++;

        axistitle[dim] = "#phi_{EMCal cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = 0;
        max[dim] = TMath::TwoPi();
        dim++;

        axistitle[dim] = "#eta_{DCal cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = -1;
        max[dim] = 1;
        dim++;

        axistitle[dim] = "#phi_{DCal cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = 0;
        max[dim] = TMath::TwoPi();
        dim++;

        axistitle[dim] = "#eta_{PHOS cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = -1;
        max[dim] = 1;
        dim++;

        axistitle[dim] = "#phi_{PHOS cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = 0;
        max[dim] = TMath::TwoPi();
        dim++;
      }
    }
    else {
      axistitle[dim] = "#it{E}_{cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;

      if (fDoLeadingObjectPosition) {
        axistitle[dim] = "#eta_{cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = -1;
        max[dim] = 1;
        dim++;

        axistitle[dim] = "#phi_{cluster}^{leading}";
        nbins[dim] = 100;
        min[dim] = 0;
        max[dim] = TMath::TwoPi();
        dim++;
      }
    }
  }

  if (!fCaloCellsName.IsNull()) {
    axistitle[dim] = "No. of cells";

    if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
      nbins[dim] = 5000;
      min[dim] = -0.5;
      max[dim] = 10000-0.5;
    }
    else {
      nbins[dim] = 500;
      min[dim] = -0.5;
      max[dim] = 500-0.5;
    }

    dim++;
  }

  THnSparse* hn = fHistManager.CreateTHnSparse("fHistEventQA","fHistEventQA",dim,nbins,min,max);
  for (Int_t i = 0; i < dim; i++)
    hn->GetAxis(i)->SetTitle(axistitle[i]);

  fOutput->Add(fHistManager.GetListOfHistograms());

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetQA::ExecOnce()
{
  if (fClusterCollArray.size() == 0  && fCaloCellsName.IsNull()) {
    fNeedEmcalGeom = kFALSE;
  }
  else {
    fNeedEmcalGeom = kTRUE;
  }

  AliAnalysisTaskEmcalLight::ExecOnce();

  if (fDoV0QA) {
    fVZERO = InputEvent()->GetVZEROData();
    if (!fVZERO) {
      AliError("AliVVZERO not available");
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalLight::RetrieveEventObjects()) return kFALSE;

  if (!fCentMethod2.IsNull() || !fCentMethod3.IsNull()) {
    if (fBeamType == kAA || fBeamType == kpA ) {
      AliCentrality *aliCent = InputEvent()->GetCentrality();
      if (aliCent) {
        if (!fCentMethod2.IsNull())
          fCent2 = aliCent->GetCentralityPercentile(fCentMethod2);
        if (!fCentMethod3.IsNull())
          fCent3 = aliCent->GetCentralityPercentile(fCentMethod3);
      }
    }
  }

  if (fVZERO) {
    fV0ATotMult = AliESDUtils::GetCorrV0A(fVZERO->GetMTotV0A(),fVertex[2]);
    fV0CTotMult = AliESDUtils::GetCorrV0C(fVZERO->GetMTotV0C(),fVertex[2]);
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetQA::FillHistograms()
{
  // Fill histograms.

  if (fCentBin < 0) {
    AliError(Form("fCentBin is %d! fCent = %.3f. Fix the centrality bins to include all possible values of centrality.", fCentBin, fCent));
    return kFALSE;
  }

  EventQA_t eventQA;

  DoTrackLoop();
  AliDebug(2,Form("%d tracks found in the event", fNTotTracks));
  eventQA.fMaxTrack = fLeadingTrack;

  DoClusterLoop();
  AliDebug(2,Form("%d clusters found in EMCal, %d in DCal and %d in PHOS", fNTotClusters[0], fNTotClusters[1], fNTotClusters[2]));
  for (Int_t i = 0; i < 3; i++) {
    eventQA.fMaxCluster[i] = fLeadingCluster[i];
  }

  if (fCaloCells) {
    eventQA.fNCells = DoCellLoop();
    AliDebug(2,Form("%d cells found in the event", eventQA.fNCells));
  }

  eventQA.fCent = fCent;
  eventQA.fCent2 = fCent2;
  eventQA.fCent3 = fCent3;
  eventQA.fV0A = fV0ATotMult;
  eventQA.fV0C = fV0CTotMult;
  eventQA.fEP = fEPV0;
  eventQA.fNTracks = fNTotTracks;
  eventQA.fNClusters[0] = fNTotClusters[0];
  eventQA.fNClusters[1] = fNTotClusters[1];
  eventQA.fNClusters[2] = fNTotClusters[2];

  FillEventQAHisto(eventQA);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetQA::FillEventQAHisto(const EventQA_t& eventQA)
{
  Double_t contents[40]={0};

  Int_t globalNclusters = eventQA.fNClusters[0] + eventQA.fNClusters[1] + eventQA.fNClusters[2];

  AliTLorentzVector globalMaxCluster;
  for (Int_t i = 0; i < 3; i++) {
    if (globalMaxCluster.E() < eventQA.fMaxCluster[i].E())  globalMaxCluster = eventQA.fMaxCluster[i];
  }

  THnSparse* histEventQA = static_cast<THnSparse*>(fHistManager.FindObject("fHistEventQA"));

  for (Int_t i = 0; i < histEventQA->GetNdimensions(); i++) {
    TString title(histEventQA->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = eventQA.fCent;
    else if (title==Form("Centrality %s %%", fCentMethod2.Data()))
      contents[i] = eventQA.fCent2;
    else if (title==Form("Centrality %s %%", fCentMethod3.Data()))
      contents[i] = eventQA.fCent3;
    else if (title=="V0A total multiplicity")
      contents[i] = eventQA.fV0A;
    else if (title=="V0C total multiplicity")
      contents[i] = eventQA.fV0C;
    else if (title=="V0A+V0C total multiplicity")
      contents[i] = eventQA.fV0A+eventQA.fV0C;
    else if (title=="#psi_{RP}")
      contents[i] = eventQA.fEP;
    else if (title=="No. of tracks")
      contents[i] = eventQA.fNTracks;
    else if (title=="No. of clusters")
      contents[i] = globalNclusters;
    else if (title=="No. of cells")
      contents[i] = eventQA.fNCells;
    else if (title=="#it{p}_{T,track}^{leading} (GeV/c)")
      contents[i] = eventQA.fMaxTrack.Pt();
    else if (title=="#eta_{track}^{leading}")
      contents[i] = eventQA.fMaxTrack.Eta();
    else if (title=="#phi_{track}^{leading}")
      contents[i] = eventQA.fMaxTrack.Phi_0_2pi();
    else if (title=="#it{E}_{cluster}^{leading} (GeV)")
      contents[i] = globalMaxCluster.E();
    else if (title=="#eta_{cluster}^{leading}")
      contents[i] = globalMaxCluster.Eta();
    else if (title=="#phi_{cluster}^{leading}")
      contents[i] = globalMaxCluster.Phi();
    else if (title=="#it{E}_{EMCal cluster}^{leading} (GeV)")
      contents[i] = eventQA.fMaxCluster[0].E();
    else if (title=="#eta_{EMCal cluster}^{leading}")
      contents[i] = eventQA.fMaxCluster[0].Phi_0_2pi();
    else if (title=="#phi_{EMCal cluster}^{leading}")
      contents[i] = eventQA.fMaxCluster[0].Eta();
    else if (title=="#it{E}_{DCal cluster}^{leading} (GeV)")
      contents[i] = eventQA.fMaxCluster[1].E();
    else if (title=="#phi_{DCal cluster}^{leading}")
      contents[i] = eventQA.fMaxCluster[1].Phi_0_2pi();
    else if (title=="#eta_{DCal cluster}^{leading}")
      contents[i] = eventQA.fMaxCluster[1].Eta();
    else if (title=="#it{E}_{PHOS cluster}^{leading} (GeV)")
      contents[i] = eventQA.fMaxCluster[2].E();
    else if (title=="#phi_{PHOS cluster}^{leading}")
      contents[i] = eventQA.fMaxCluster[2].Phi_0_2pi();
    else if (title=="#eta_{PHOS cluster}^{leading}")
      contents[i] = eventQA.fMaxCluster[2].Eta();
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  histEventQA->Fill(contents);
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetQA::DoCellLoop()
{
  // Do cell loop.

  if (!fCaloCells) return 0;

  TString histname_en = TString::Format("%s/fHistCellsAbsIdEnergy_%d", fCaloCellsName.Data(), fCentBin);
  TString histname_tm = TString::Format("%s/fHistCellsAbsIdTime_%d", fCaloCellsName.Data(), fCentBin);

  const Int_t ncells = fCaloCells->GetNumberOfCells();
  Int_t nAccCells = 0;

  for (Int_t pos = 0; pos < ncells; pos++) {
    Float_t amp   = fCaloCells->GetAmplitude(pos);
    Float_t time   = fCaloCells->GetTime(pos);
    Int_t   absId = fCaloCells->GetCellNumber(pos);

    if (amp < fCellEnergyCut) continue;

    fHistManager.FillTH2(histname_en, absId, amp);
    fHistManager.FillTH2(histname_tm, absId, time);
    nAccCells++;
  }

  return nAccCells;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetQA::GetFcross(AliVCluster *cluster, AliVCaloCells *cells)
{
  Int_t    AbsIdseed  = -1;
  Double_t Eseed      = 0;
  for (Int_t i = 0; i < cluster->GetNCells(); i++) {
    if (cells->GetCellAmplitude(cluster->GetCellAbsId(i)) > Eseed) {
      Eseed     = cells->GetCellAmplitude(cluster->GetCellAbsId(i));
      AbsIdseed = cluster->GetCellAbsId(i);
    }
  }

  if (Eseed < 1e-9)
    return 100;

  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
  fGeom->GetCellIndex(AbsIdseed,imod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,iphi,ieta);

  //Get close cells index and energy, not in corners

  Int_t absID1 = -1;
  Int_t absID2 = -1;

  if (iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if (iphi > 0)                                 absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

  // In case of cell in eta = 0 border, depending on SM shift the cross cell index

  Int_t absID3 = -1;
  Int_t absID4 = -1;

  if (ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2)) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
  }
  else if (ieta == 0 && imod%2) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
  }
  else {
    if (ieta < AliEMCALGeoParams::fgkEMCALCols-1)
      absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if (ieta > 0)
      absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
  }

  Double_t  ecell1 = cells->GetCellAmplitude(absID1);
  Double_t  ecell2 = cells->GetCellAmplitude(absID2);
  Double_t  ecell3 = cells->GetCellAmplitude(absID3);
  Double_t  ecell4 = cells->GetCellAmplitude(absID4);

  Double_t Ecross = ecell1 + ecell2 + ecell3 + ecell4;

  Double_t Fcross = 1 - Ecross/Eseed;

  return Fcross;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetQA::DoClusterLoop()
{
  // Do cluster loop.

  TString histname;

  memset(fNTotClusters, 0, sizeof(Int_t)*3);
  for (Int_t i = 0; i < 3; i++) fLeadingCluster[i].SetPxPyPzE(0,0,0,0);

  for (auto cont_it : fClusterCollArray) {
    AliClusterContainer* clusters = cont_it.second;
    // Cluster loop (EMCal, DCal, PHOS)
    AliClusterIterableMomentumContainer itcont = clusters->all_momentum();
    for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      UInt_t rejectionReason = 0;
      if (!clusters->AcceptCluster(it.current_index(), rejectionReason)) {
        histname = TString::Format("%s/fHistRejectionReason_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
        continue;
      }

      Float_t pos[3]={0};
      it->second->GetPosition(pos);
      histname = TString::Format("%s/fHistClusPosition_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH3(histname, pos[0], pos[1], pos[2]);

      histname = TString::Format("%s/fHistClusPhiEtaEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH3(histname, it->first.Eta(), it->first.Phi_0_2pi(), it->first.E());

      histname = TString::Format("%s/fHistClusDeltaPhiEPEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      if (fHistManager.FindObject(histname)) {
        Double_t ep = it->first.Phi_0_2pi() - fEPV0;
        while (ep < 0) ep += TMath::Pi();
        while (ep >= TMath::Pi()) ep -= TMath::Pi();
        fHistManager.FillTH2(histname, it->first.E(), ep);
      }

      histname = TString::Format("%s/fHistNCellsEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH2(histname, it->first.E(), it->second->GetNCells());

      histname = TString::Format("%s/fHistClusTimeEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH2(histname, it->first.E(), it->second->GetTOF());

      histname = TString::Format("%s/fHistClusMCEnergyFraction_%d", clusters->GetArrayName().Data(), fCentBin);
      if (fHistManager.FindObject(histname)) {
        fHistManager.FillTH1(histname, it->second->GetMCEnergyFraction());
      }

      histname = TString::Format("%s/fHistClusEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH1(histname, it->second->E());

      if (it->second->GetNonLinCorrEnergy() > 0.) {
        histname = TString::Format("%s/fHistClusNonLinCorrEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH1(histname, it->second->GetNonLinCorrEnergy());
      }

      if (it->second->GetHadCorrEnergy() > 0.) {
        histname = TString::Format("%s/fHistClusHadCorrEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH1(histname, it->second->GetHadCorrEnergy());
      }

      // The following histograms are filled only for EMCal/DCal
      if (it->second->IsEMCAL()) {
        Double_t phi = it->first.Phi_0_2pi();
        Int_t isDcal = Int_t(phi > fgkEMCalDCalPhiDivide);
        if (fLeadingCluster[isDcal].E() < it->first.E()) fLeadingCluster[isDcal] = it->first;
        fNTotClusters[isDcal]++;

        if (fCaloCells) {
          histname = TString::Format("%s/fHistFcrossEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
          fHistManager.FillTH2(histname, it->first.E(), GetFcross(it->second, fCaloCells));
        }

        Int_t sm = fGeom->GetSuperModuleNumber(it->second->GetCellAbsId(0));
        if (sm >=0 && sm < 20) {
          histname = TString::Format("%s/BySM/fHistClusEnergy_SM%d_%d", clusters->GetArrayName().Data(), sm, fCentBin);
          fHistManager.FillTH1(histname, it->second->E());
        }
        else {
          AliError(Form("Supermodule %d does not exist!", sm));
        }
      }
      else if (it->second->IsPHOS()) {
        fNTotClusters[2]++;
        if (fLeadingCluster[2].E() < it->first.E()) fLeadingCluster[2] = it->first;
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetQA::DoTrackLoop()
{
  // Do track loop.

  Int_t neg = 0;
  Int_t zero = 0;

  fNTotTracks = 0;
  fLeadingTrack.SetPxPyPzE(0,0,0,0);

  TString histname;

  for (auto cont_it : fParticleCollArray) {
    AliParticleContainer* particles = cont_it.second;
    AliParticleIterableMomentumContainer itcont = particles->all_momentum();
    for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {

      UInt_t rejectionReason = 0;
      if (!particles->AcceptParticle(it.current_index(), rejectionReason)) {
        histname = TString::Format("%s/fHistRejectionReason_%d", particles->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname, particles->GetRejectionReasonBitPosition(rejectionReason), it->first.Pt());
        continue;
      }

      fNTotTracks++;

      if (fLeadingTrack.Pt() < it->first.Pt()) fLeadingTrack = it->first;

      if (fParticleLevel) {
        histname = TString::Format("%s/fHistTrPhiEtaPt_%d_0", particles->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH3(histname, it->first.Eta(), it->first.Phi_0_2pi(), it->first.Pt());
      }
      else {
        if (it->second->GetLabel() == 0) {
          zero++;
          histname = TString::Format("%s/fHistTrPhiEtaZeroLab_%d", particles->GetArrayName().Data(), fCentBin);
          if (fHistManager.FindObject(histname)) {

            fHistManager.FillTH2(histname, it->first.Eta(), it->first.Phi_0_2pi());

          }
          histname = TString::Format("%s/fHistTrPtZeroLab_%d", particles->GetArrayName().Data(), fCentBin);
          if (fHistManager.FindObject(histname)) {
            fHistManager.FillTH1(histname, it->first.Pt());
          }
        }

        if (it->second->GetLabel() < 0) neg++;

        Int_t type = 0;
        AliTrackContainer* tracks = dynamic_cast<AliTrackContainer*>(particles);
        if (tracks) {
          // Track type (hybrid)
          if (tracks->GetTrackFilterType() == AliEmcalTrackSelection::kHybridTracks) {
            type = tracks->GetTrackType(it.current_index());
          }
        }

        // Track propagation to EMCal surface
        AliVTrack* vtrack = dynamic_cast<AliVTrack*>(it->second);
        if (vtrack) {
          if (vtrack->GetTrackEtaOnEMCal() == -999 || vtrack->GetTrackPhiOnEMCal() == -999) {
            histname = TString::Format("%s/fHistTrPhiEtaNonProp_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname)) {
              fHistManager.FillTH2(histname, it->first.Eta(), it->first.Phi_0_2pi());
            }
            histname = TString::Format("%s/fHistTrPtNonProp_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname)) {
              fHistManager.FillTH1(histname, it->first.Pt());
            }
          }
          else {
            histname = TString::Format("%s/fHistTrEmcPhiEta_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH2(histname, vtrack->GetTrackEtaOnEMCal(), vtrack->GetTrackPhiOnEMCal());

            histname = TString::Format("%s/fHistTrEmcPt_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH1(histname, vtrack->GetTrackPtOnEMCal());

            histname = TString::Format("%s/fHistDeltaEtaPt_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH2(histname, it->first.Pt(), it->first.Eta() - vtrack->GetTrackEtaOnEMCal());

            histname = TString::Format("%s/fHistDeltaPhiPt_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH2(histname, it->first.Pt(), it->first.Phi_0_2pi() - vtrack->GetTrackPhiOnEMCal());

            histname = TString::Format("%s/fHistDeltaPtvsPt_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH2(histname, it->first.Pt(), it->first.Pt() - vtrack->GetTrackPtOnEMCal());
          }
        }

        if (type >= 0 && type <= 3) {
          histname = TString::Format("%s/fHistTrPhiEtaPt_%d_%d", particles->GetArrayName().Data(), fCentBin, type);
          fHistManager.FillTH3(histname, it->first.Eta(), it->first.Phi_0_2pi(), it->first.Pt());
        }
        else {
          AliWarning(Form("%s: track type %d not recognized!", GetName(), type));
        }
      }
    }

    histname = TString::Format("%s/fHistTrNegativeLabels_%d", particles->GetArrayName().Data(), fCentBin);
    if (fHistManager.FindObject(histname))
      fHistManager.FillTH1(histname, 1. * neg / fNTotTracks);

    histname = TString::Format("%s/fHistTrZeroLabels_%d", particles->GetArrayName().Data(), fCentBin);
    if (fHistManager.FindObject(histname))
      fHistManager.FillTH1(histname, 1. * zero / fNTotTracks);
  }
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetQA* AliAnalysisTaskEmcalJetQA::AddTaskEmcalJetQA(TString ntracks, TString nclusters, TString ncells, TString subdir, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalJetQA", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskEmcalJetQA", "This task requires an input event handler");
    return NULL;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings

  if (ntracks == "usedefault") {
    if (dataType == kESD) {
      ntracks = "Tracks";
    }
    else if (dataType == kAOD) {
      ntracks = "tracks";
    }
    else {
      ntracks = "";
    }
  }

  if (nclusters == "usedefault") {
    if (dataType == kESD) {
      nclusters = "CaloClusters";
    }
    else if (dataType == kAOD) {
      nclusters = "caloClusters";
    }
    else {
      nclusters = "";
    }
  }

  if (ncells == "usedefault") {
    if (dataType == kESD) {
      ncells = "EMCALCells";
    }
    else if (dataType == kAOD) {
      ncells = "emcalCells";
    }
    else {
      ncells = "";
    }
  }

  TString name("AliAnalysisTaskEmcalJetQA");
  if (!ntracks.IsNull()) {
    name += "_";
    name += ntracks;
  }
  if (!nclusters.IsNull()) {
    name += "_";
    name += nclusters;
  }
  if (!ncells.IsNull()) {
    name += "_";
    name += ncells;
  }
  if (!suffix.IsNull() != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmcalJetQA* qaTask = new AliAnalysisTaskEmcalJetQA(name);
  qaTask->SetCaloCellsName(ncells);
  qaTask->SetVzRange(-10,10);
  qaTask->AddParticleContainer(ntracks.Data());
  qaTask->AddClusterContainer(nclusters.Data());

  // Final settings, pass to manager and set the containers
  mgr->AddTask(qaTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  TString contName = TString::Format("%s_histos", name.Data());
  TString commonoutput;
  if (subdir.IsNull()) {
    commonoutput = mgr->GetCommonFileName();
  }
  else {
    commonoutput = TString::Format("%s:%s", mgr->GetCommonFileName(), subdir.Data());
  }

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),
                  TList::Class(),AliAnalysisManager::kOutputContainer,
                  commonoutput);
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );

  return qaTask;
}

/**
 * Add this task to the QA train
 * \param runnumber Run number
 */
void AliAnalysisTaskEmcalJetQA::AddTaskEmcalJetQA_QAtrain(Int_t runnumber)
{
  EBeamType_t beam = BeamTypeFromRunNumber(runnumber);
  Int_t nCentBins = 0;
  if (beam == kpA || beam == kAA) nCentBins = 4;
  std::vector<std::string> triggerClasses = {"CINT7", "CEMC7", "CDMC7", "EG1", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2" };
  for (auto triggerClass : triggerClasses) {
    TString suffix(triggerClass.c_str());
    suffix.ReplaceAll("-", "_");
    AliAnalysisTaskEmcalJetQA* task = AddTaskEmcalJetQA("", "usedefault", "usedefault", "CaloQA_default", suffix);
    task->AddAcceptedTriggerClass(triggerClass.c_str());
    task->SetForceBeamType(beam);
    if (runnumber == 0 || (runnumber >= 265077 && runnumber <= 999999)) { // Run-2 p-Pb (LHC16q): disabling vertex cut
      task->SetVzRange(-999,999);
    }
  }
}
