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
  AliAnalysisTaskEmcal("AliAnalysisTaskEmcalJetQA", kTRUE),
  fCellEnergyCut(0.05),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fDoLeadingObjectPosition(0),
  fMaxCellsInCluster(50),
  fSeparateEMCalDCal(kTRUE),
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

  memset(fNTotClusters, 0, sizeof(Int_t)*2);
  memset(fLeadingCluster, 0, sizeof(AliVCluster*)*2);

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetQA::AliAnalysisTaskEmcalJetQA(const char *name) :
  AliAnalysisTaskEmcal(name, kTRUE),
  fCellEnergyCut(0.05),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fDoLeadingObjectPosition(0),
  fMaxCellsInCluster(50),
  fSeparateEMCalDCal(kTRUE),
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

  memset(fNTotClusters, 0, sizeof(Int_t)*2);
  memset(fLeadingCluster, 0, sizeof(AliVCluster*)*2);

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

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  AliEmcalContainer* cont = 0;

  TString histname;
  TString title;

  TIter nextPartColl(&fParticleCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextPartColl()))) {
    fHistManager.CreateHistoGroup(cont->GetArrayName());
    if (!fParticleLevel && fIsMC) {
      for (Int_t i = 0; i < fNcentBins; i++) {
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

    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = TString::Format("%s/fHistRejectionReason_%d", cont->GetArrayName().Data(), i);
      title = histname + ";Rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 40, 0, 100);
      SetRejectionReasonLabels(hist->GetXaxis());

      for (Int_t j = 0; j < nlabels; j++) {
        histname = TString::Format("%s/fHistTrPhiEtaPt_%d_%d", cont->GetArrayName().Data(), i, j);
        title = histname + ";#eta;#phi;#it{p}_{T} (GeV/#it{c})";
        fHistManager.CreateTH3(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi(), fNbins, fMinBinPt, fMaxBinPt);
      }

      if (!fParticleLevel) {
        if (fIsMC) {
          histname = TString::Format("%s/fHistTrPhiEtaZeroLab_%d", cont->GetArrayName().Data(), i);
          title = histname + ";#eta;#phi;counts";
          fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi());

          histname = TString::Format("%s/fHistTrPtZeroLab_%d", cont->GetArrayName().Data(), i);
          title = histname + ";#it{p}_{T} (GeV/#it{c});counts";
          fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt);
        }

        histname = TString::Format("%s/fHistTrEmcPhiEta_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#eta;#phi;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi());

        histname = TString::Format("%s/fHistTrEmcPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt);

        histname = TString::Format("%s/fHistTrPhiEtaNonProp_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#eta;#phi;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1, 1, 100, 0, TMath::TwoPi());

        histname = TString::Format("%s/fHistTrPtNonProp_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt);

        histname = TString::Format("%s/fHistDeltaEtaPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});#delta#eta;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});#delta#phi;counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T} (GeV/#it{c});#delta#it{p}_{T} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, -fMaxBinPt/2, fMaxBinPt/2);
      }
    }
  }

  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) {
    fHistManager.CreateHistoGroup(cont->GetArrayName());
    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = TString::Format("%s/fHistRejectionReason_%d", cont->GetArrayName().Data(), i);
      title = histname + ";Rejection reason;#it{E}_{cluster} (GeV);counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 40, 0, 100);
      SetRejectionReasonLabels(hist->GetXaxis());

      histname = TString::Format("%s/fHistClusPosition_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{x} (cm);#it{y} (cm);#it{z} (cm)";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 50, -500, 500, 50, -500, 500, 50, -500, 500);

      histname = TString::Format("%s/fHistClusPhiEtaEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#eta;#phi;#it{E}_{cluster} (GeV)";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 50, -1, 1, 50, 0, TMath::TwoPi(), fNbins, fMinBinPt, fMaxBinPt);

      if (fForceBeamType != kpp) {
        histname = TString::Format("%s/fHistClusDeltaPhiEPEnergy_%d", cont->GetArrayName().Data(), i);
        title = histname + ";#it{E}_{cluster} (GeV);#phi_{cluster} - #psi_{EP};counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 100, 0, TMath::Pi());
      }

      if (fIsEmbedded) {
        histname = TString::Format("%s/fHistClusMCEnergyFraction_%d", cont->GetArrayName().Data(), i);
        title = histname + ";MC fraction;counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, 0, 1.2);
      }

      histname = TString::Format("%s/fHistClusTimeEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);time (s);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, -1e-6, 1e-6);

      Int_t nbins = fMaxCellsInCluster;
      while (nbins > fNbins) nbins /= 2;
      histname = TString::Format("%s/fHistNCellsEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);#it{N}_{cells};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, nbins, -0.5, fMaxCellsInCluster-0.5);

      histname = TString::Format("%s/fHistFcrossEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);#it{F}_{cross};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 200, -3.5, 1.5);

      histname = TString::Format("%s/fHistClusEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/fHistClusNonLinCorrEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/fHistClusHadCorrEnergy_%d", cont->GetArrayName().Data(), i);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt);
    }
  }

  if (!fCaloCellsName.IsNull()) {
    fHistManager.CreateHistoGroup(fCaloCellsName);
    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = TString::Format("%s/fHistCellsAbsIdEnergy_%d", fCaloCellsName.Data(), i);
      title = histname + ";cell abs. Id;#it{E}_{cell} (GeV);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 20000,0,20000,(Int_t)(fNbins / 2), fMinBinPt, fMaxBinPt / 2);
    }
  }

  Int_t dim = 0;
  TString axistitle[40];
  Int_t nbins[40] = {0};
  Double_t min[40] = {0};
  Double_t max[40] = {0};

  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    axistitle[dim] = "Centrality %";
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 101;
    dim++;

    if (!fCentMethod2.IsNull()) {
      axistitle[dim] = Form("Centrality %s %%", fCentMethod2.Data());
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = 101;
      dim++;
    }

    if (!fCentMethod3.IsNull()) {
      axistitle[dim] = Form("Centrality %s %%", fCentMethod3.Data());
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = 101;
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

  if (fParticleCollArray.GetEntriesFast()>0) {
    axistitle[dim] = "No. of tracks";
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
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
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
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

  if (fClusterCollArray.GetEntriesFast()>0) {
    axistitle[dim] = "No. of clusters";

    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
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
      nbins[dim] = fNbins;
      min[dim] = fMinBinPt;
      max[dim] = fMaxBinPt;
      dim++;

      axistitle[dim] = "#it{E}_{DCal cluster}^{leading} (GeV)";
      nbins[dim] = fNbins;
      min[dim] = fMinBinPt;
      max[dim] = fMaxBinPt;
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
      }
    }
    else {
      axistitle[dim] = "#it{E}_{cluster}^{leading} (GeV)";
      nbins[dim] = fNbins;
      min[dim] = fMinBinPt;
      max[dim] = fMaxBinPt;
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

    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
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
  if (fClusterCollArray.GetEntriesFast() == 0  && fCaloCellsName.IsNull()) {
    fNeedEmcalGeom = kFALSE;
  }
  else {
    fNeedEmcalGeom = kTRUE;
  }

  AliAnalysisTaskEmcal::ExecOnce();

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

  if (!AliAnalysisTaskEmcal::RetrieveEventObjects()) return kFALSE;

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

  EventQA_t eventQA;

  if (fTracks) {
    DoTrackLoop();
    AliDebug(2,Form("%d tracks found in the event", fNTotTracks));

    eventQA.fMaxTrack = fLeadingTrack;
  } 

  if (fCaloClusters) {
    DoClusterLoop();
    AliDebug(2,Form("%d clusters found in EMCal and %d in DCal", fNTotClusters[0], fNTotClusters[1]));


    for (Int_t i = 0; i < 2; i++) {
      if (!fLeadingCluster[i]) continue;
      fLeadingCluster[i]->GetMomentum(eventQA.fMaxCluster[i], fVertex);
    }
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

  FillEventQAHisto(eventQA);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetQA::FillEventQAHisto(const EventQA_t& eventQA)
{
  Double_t contents[40]={0};

  Int_t globalNclusters = eventQA.fNClusters[0] + eventQA.fNClusters[1];

  AliTLorentzVector globalMaxCluster = eventQA.fMaxCluster[0].E() > eventQA.fMaxCluster[1].E() ?
      eventQA.fMaxCluster[0] : eventQA.fMaxCluster[1];

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

  TString histname = TString::Format("%s/fHistCellsAbsIdEnergy_%d", fCaloCellsName.Data(), fCentBin);

  const Int_t ncells = fCaloCells->GetNumberOfCells();
  Int_t nAccCells = 0;

  for (Int_t pos = 0; pos < ncells; pos++) {
    Float_t amp   = fCaloCells->GetAmplitude(pos);
    Int_t   absId = fCaloCells->GetCellNumber(pos);

    if (amp < fCellEnergyCut) continue;

    fHistManager.FillTH2(histname, absId,amp);
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

  memset(fNTotClusters, 0, sizeof(Int_t)*2);
  memset(fLeadingCluster, 0, sizeof(AliVCluster*)*2);

  AliClusterContainer* clusters = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    // Cluster loop
    AliVCluster* cluster = 0;
    clusters->ResetCurrentID();
    while ((cluster = clusters->GetNextCluster())) {
      AliTLorentzVector nPart;
      Double_t energy = 0;

      if (clusters->GetDefaultClusterEnergy() >= 0 && clusters->GetDefaultClusterEnergy() <= AliVCluster::kLastUserDefEnergy) {
        energy = cluster->GetUserDefEnergy(clusters->GetDefaultClusterEnergy());
        cluster->GetMomentum(nPart, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      }
      else {
        energy = cluster->E();
        cluster->GetMomentum(nPart, fVertex);
      }

      if (energy <= 0) continue;

      UInt_t rejectionReason = 0;
      if (!clusters->AcceptCluster(clusters->GetCurrentID(), rejectionReason)) {
        histname = TString::Format("%s/fHistRejectionReason_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), energy);
        continue;
      }

      Float_t pos[3]={0};
      cluster->GetPosition(pos);
      histname = TString::Format("%s/fHistClusPosition_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH3(histname, pos[0], pos[1], pos[2]);

      Double_t phi = nPart.Phi_0_2pi();

      Int_t isDcal = Int_t(phi > fgkEMCalDCalPhiDivide);

      histname = TString::Format("%s/fHistClusPhiEtaEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH3(histname, nPart.Eta(), phi, energy);

      if (!fLeadingCluster[isDcal] || fLeadingCluster[isDcal]->E() < cluster->E()) fLeadingCluster[isDcal] = cluster;

      histname = TString::Format("%s/fHistClusDeltaPhiEPEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      if (fHistManager.FindObject(histname)) {
        Double_t ep = nPart.Phi_0_2pi() - fEPV0;
        while (ep < 0) ep += TMath::Pi();
        while (ep >= TMath::Pi()) ep -= TMath::Pi();
        fHistManager.FillTH2(histname, cluster->E(), ep);
      }

      histname = TString::Format("%s/fHistNCellsEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH2(histname, cluster->E(), cluster->GetNCells());

      histname = TString::Format("%s/fHistClusTimeEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH2(histname, cluster->E(), cluster->GetTOF());

      if (fCaloCells) {
        histname = TString::Format("%s/fHistFcrossEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname, cluster->E(), GetFcross(cluster, fCaloCells));
      }

      histname = TString::Format("%s/fHistClusMCEnergyFraction_%d", clusters->GetArrayName().Data(), fCentBin);
      if (fHistManager.FindObject(histname)) {
        fHistManager.FillTH1(histname, cluster->GetMCEnergyFraction());
      }

      histname = TString::Format("%s/fHistClusEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->E());

      if (cluster->GetNonLinCorrEnergy() > 0.) {
        histname = TString::Format("%s/fHistClusNonLinCorrEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());
      }

      if (cluster->GetHadCorrEnergy() > 0.) {
        histname = TString::Format("%s/fHistClusHadCorrEnergy_%d", clusters->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());
      }

      fNTotClusters[isDcal]++;
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

  AliParticleContainer* particles = 0;
  TIter nextPartColl(&fParticleCollArray);
  while ((particles = static_cast<AliParticleContainer*>(nextPartColl()))) {
    particles->ResetCurrentID();
    AliVParticle* track = 0;
    AliTLorentzVector mom;
    while ((track = particles->GetNextParticle())) {
      particles->GetMomentum(mom, particles->GetCurrentID());

      UInt_t rejectionReason = 0;
      if (!particles->AcceptParticle(particles->GetCurrentID(), rejectionReason)) {
        histname = TString::Format("%s/fHistRejectionReason_%d", particles->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname, particles->GetRejectionReasonBitPosition(rejectionReason), mom.Pt());
        continue;
      }

      fNTotTracks++;

      if (fLeadingTrack.Pt() < mom.Pt()) fLeadingTrack = mom;

      if (fParticleLevel) {
        histname = TString::Format("%s/fHistTrPhiEtaPt_%d", particles->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname, mom.Eta(), mom.Phi_0_2pi(), mom.Pt());
      }
      else {
        if (track->GetLabel() == 0) {
          zero++;
          histname = TString::Format("%s/fHistTrPhiEtaZeroLab_%d", particles->GetArrayName().Data(), fCentBin);
          if (fHistManager.FindObject(histname)) {

            fHistManager.FillTH2(histname, mom.Eta(), mom.Phi_0_2pi());

          }
          histname = TString::Format("%s/fHistTrPtZeroLab_%d", particles->GetArrayName().Data(), fCentBin);
          if (fHistManager.FindObject(histname)) {
            fHistManager.FillTH1(histname, mom.Pt());
          }
        }

        if (track->GetLabel() < 0) neg++;

        Int_t type = 0;
        AliTrackContainer* tracks = dynamic_cast<AliTrackContainer*>(particles);

        if (tracks) {
          // Track type (hybrid)
          if (tracks->GetTrackFilterType() == AliEmcalTrackSelection::kHybridTracks) {
            type = tracks->GetTrackType(tracks->GetCurrentID());
          }

          // Track propagation to EMCal surface
          AliVTrack* vtrack = static_cast<AliVTrack*>(track);

          if (vtrack->GetTrackEtaOnEMCal() == -999 || vtrack->GetTrackPhiOnEMCal() == -999) {
            histname = TString::Format("%s/fHistTrPhiEtaNonProp_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname)) {
              fHistManager.FillTH2(histname, mom.Eta(), mom.Phi_0_2pi());
            }
            histname = TString::Format("%s/fHistTrPtNonProp_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname)) {
              fHistManager.FillTH1(histname, mom.Pt());
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
              fHistManager.FillTH2(histname, mom.Pt(), mom.Eta() - vtrack->GetTrackEtaOnEMCal());

            histname = TString::Format("%s/fHistDeltaPhiPt_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH2(histname, mom.Pt(), mom.Phi_0_2pi() - vtrack->GetTrackPhiOnEMCal());

            histname = TString::Format("%s/fHistDeltaPtvsPt_%d", particles->GetArrayName().Data(), fCentBin);
            if (fHistManager.FindObject(histname))
              fHistManager.FillTH2(histname, mom.Pt(), mom.Pt() - vtrack->GetTrackPtOnEMCal());
          }
        }

        if (type >= 0 && type <= 3) {
          histname = TString::Format("%s/fHistTrPhiEtaPt_%d_%d", particles->GetArrayName().Data(), fCentBin, type);
          fHistManager.FillTH3(histname, mom.Eta(), mom.Phi_0_2pi(), mom.Pt());
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
