/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

#include <array>
#include <memory>

#include <TRandom3.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>

#include "AliEmcalJet.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliClusterContainer.h"
#include "AliParticleContainer.h"

#include "AliAnalysisTaskJetUEStudies.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetUEStudies);
/// \endcond

/// Default constructor for ROOT I/O purposes
AliAnalysisTaskJetUEStudies::AliAnalysisTaskJetUEStudies() :
  AliAnalysisTaskJetUE(),
  fAlternativeRho(),
  fHistManager(),
  fRandom(nullptr)

{
}

/// Standard named constructor
///
/// \param name Name of the task
AliAnalysisTaskJetUEStudies::AliAnalysisTaskJetUEStudies(const char *name) :
  AliAnalysisTaskJetUE(name, kTRUE),
  fAlternativeRho(),
  fHistManager(name),
  fRandom(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);
}

/// Overloads base class method. Creates output objects
void AliAnalysisTaskJetUEStudies::UserCreateOutputObjects()
{
  AliInfo(Form("CreateOutputObjects of task %s", GetName()));

  AliAnalysisTaskEmcalJetLight::UserCreateOutputObjects();

  Int_t maxTracks = 6000;
  Int_t constituentsNbins = 250;
  Double_t constituentsMax = 249.5;
  Double_t maxRho = 500;
  Double_t minCorrPt = -250;

  if (fForceBeamType == kpp) {
    constituentsNbins = 50;
    constituentsMax = 49.5;
    maxRho = 50;
    maxTracks = 200;
    minCorrPt = -50;
  }
  else if (fForceBeamType == kpA) {
    constituentsNbins = 100;
    constituentsMax = 99.5;
    maxRho = 200;
    maxTracks = 500;
    minCorrPt = -100;
  }

  // Avoid half bins
  minCorrPt -= minCorrPt - TMath::FloorNint(minCorrPt / fPtBinWidth) * fPtBinWidth;

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);
  Int_t nCorrPtBins = TMath::CeilNint((fMaxPt - minCorrPt) / fPtBinWidth);

  TString histname;
  TString title;

  for (auto cont_it : fJetCollArray) {
    AliJetContainer* jets = cont_it.second;
    if (jets->GetRhoName().IsNull()) continue;
    fDefaultRhoNames.insert(jets->GetRhoName());
    fAlternativeRho.insert(std::make_pair(jets->GetRhoName(), nullptr));
  }

  for (auto rho1 : fAlternativeRho) {
    for (auto rho2 : fAlternativeRho) {
      if (rho1.first == rho2.first) continue;
      histname = TString::Format("%s/fHistDelta%sOverRhoVsCent", rho1.first.Data(), rho2.first.Data());
      title = TString::Format("%s;Centrality (%%);2#frac{%s - %s}{%s + %s};counts", histname.Data(), rho1.first.Data(), rho2.first.Data(), rho1.first.Data(), rho2.first.Data());
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, 500, -2.5, 2.5);

      histname = TString::Format("%s/fHistDelta%sVsRho", rho1.first.Data(), rho2.first.Data());
      title = TString::Format("%s;#frac{%s + %s}{2} (GeV/#it{c} #times rad^{-1});%s - %s (GeV/#it{c} #times rad^{-1});counts", histname.Data(), rho1.first.Data(), rho2.first.Data(), rho1.first.Data(), rho2.first.Data());
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*4, minCorrPt*2, -minCorrPt*2);

      histname = TString::Format("%s/fHist%sVs%s", rho1.first.Data(), rho2.first.Data(), rho1.first.Data());
      title = TString::Format("%s;%s (GeV/#it{c} #times rad^{-1}); %s (GeV/#it{c} #times rad^{-1});counts", histname.Data(), rho1.first.Data(), rho2.first.Data());
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, 500, 0, maxRho);
    }
  }

  for (auto cont_it : fJetCollArray) {
    AliJetContainer* jets = cont_it.second;

    histname = TString::Format("%s/fHistLeadingJetPtVsCent", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T,jet}^{lead} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, nPtBins, 0, fMaxPt);

    histname = TString::Format("%s/fHistLeadingJetPhiVsEta", jets->GetArrayName().Data());
    title = histname + ";#eta_{jet};#phi_{jet} (rad);counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 50, -1, 1, 150, 0, TMath::TwoPi());

    histname = TString::Format("%s/fHistB2BLeadingJetPtVsCent", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T,jet}^{lead} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, nPtBins, 0, fMaxPt);

    for (auto rho : fAlternativeRho) {
      for (auto rho2 : fAlternativeRho) {
        if (rho.first == rho2.first) continue;
        histname = TString::Format("%s/%s/fHistB2BDelta%sOverRhoVsCent", jets->GetArrayName().Data(), rho.first.Data(), rho2.first.Data());
        title = TString::Format("%s;Centrality (%%);2#frac{%s - %s}{%s + %s};counts", histname.Data(), rho.first.Data(), rho2.first.Data(), rho.first.Data(), rho2.first.Data());
        fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, 500, -2.5, 2.5);

        histname = TString::Format("%s/%s/fHistB2BDelta%sVsRho", jets->GetArrayName().Data(), rho.first.Data(), rho2.first.Data());
        title = TString::Format("%s;#frac{%s + %s}{2} (GeV/#it{c} #times rad^{-1});%s - %s (GeV/#it{c} #times rad^{-1});counts", histname.Data(), rho.first.Data(), rho2.first.Data(), rho.first.Data(), rho2.first.Data());
        fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*4, minCorrPt*2, -minCorrPt*2);

        histname = TString::Format("%s/fHistB2B%sVs%s", rho.first.Data(), rho2.first.Data(), rho.first.Data());
        title = TString::Format("%s;%s (GeV/#it{c} #times rad^{-1}); %s (GeV/#it{c} #times rad^{-1});counts", histname.Data(), rho.first.Data(), rho2.first.Data());
        fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, 500, 0, maxRho);
      }

      histname = TString::Format("%s/%s/fHistLeadingJetCorrPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#it{p}_{T,jet}^{lead} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, nCorrPtBins, minCorrPt, fMaxPt);

      histname = TString::Format("%s/%s/fHistB2BLeadingJetCorrPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#it{p}_{T,jet}^{lead} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, nCorrPtBins, minCorrPt, fMaxPt);

      histname = TString::Format("%s/%s/fHistRCPhiVsEta", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#eta_{RC};#phi_{RC} (rad);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, -1, 1, 150, 0, TMath::TwoPi());

      histname = TString::Format("%s/%s/fHistRCPerpPhiVsEta", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#eta_{RC};#phi_{RC} (rad);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, -1, 1, 150, 0, TMath::TwoPi());

      histname = TString::Format("%s/%s/fHistRCExclLeadJetPhiVsEta", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#eta_{RC};#phi_{RC} (rad);counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, -1, 1, 150, 0, TMath::TwoPi());

      histname = TString::Format("%s/%s/fHistRCDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";Centrality (%);#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#rho (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#rho (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#rho (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#rho (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#rho (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#rho (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, maxRho, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);

      histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});#Delta#it{p}_{T}^{RC} = #it{p}_{T}^{RC} - #it{A}_{jet}#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, (nCorrPtBins - nPtBins)*2, minCorrPt, -minCorrPt);
    }

    if (fCentBins.size() > 1) {
      for (Int_t i = 0; i < fCentBins.size()-1; i++) {
        histname = TString::Format("%s/fHistJetPt_Cent%d_%d", jets->GetArrayName().Data(), TMath::CeilNint(fCentBins[i]), TMath::CeilNint(fCentBins[i+1]));
        title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

        histname = TString::Format("%s/fHistB2BJetPt_Cent%d_%d", jets->GetArrayName().Data(), TMath::CeilNint(fCentBins[i]), TMath::CeilNint(fCentBins[i+1]));
        title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

        for (auto rho : fAlternativeRho) {
          histname = TString::Format("%s/%s/fHistJetCorrPt_Cent%d_%d", jets->GetArrayName().Data(), rho.first.Data(), TMath::CeilNint(fCentBins[i]), TMath::CeilNint(fCentBins[i+1]));
          title = histname + ";#it{p}_{T,jet}^{lead} - #it{A}_{jet}#rho (GeV/#it{c});counts";
          fHistManager.CreateTH1(histname.Data(), title.Data(), nCorrPtBins, minCorrPt, fMaxPt);

          histname = TString::Format("%s/%s/fHistB2BJetCorrPt_Cent%d_%d", jets->GetArrayName().Data(), rho.first.Data(), TMath::CeilNint(fCentBins[i]), TMath::CeilNint(fCentBins[i+1]));
          title = histname + ";#it{p}_{T,jet}^{lead} - #it{A}_{jet}#rho (GeV/#it{c});counts";
          fHistManager.CreateTH1(histname.Data(), title.Data(), nCorrPtBins, minCorrPt, fMaxPt);
        }
      }
    }
    else {
      histname = TString::Format("%s/fHistJetPt", jets->GetArrayName().Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

      histname = TString::Format("%s/fHistB2BJetPt", jets->GetArrayName().Data());
      title = histname + ";#it{p}_{T,jet}^{lead} (GeV/#it{c});counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);

      for (auto rho : fAlternativeRho) {
        histname = TString::Format("%s/%s/fHistJetCorrPt", jets->GetArrayName().Data(), rho.first.Data());
        title = histname + ";#it{p}_{T,jet}^{lead} - #it{A}_{jet}#rho (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nCorrPtBins, minCorrPt, fMaxPt);

        histname = TString::Format("%s/%s/fHistB2BJetCorrPt", jets->GetArrayName().Data(), rho.first.Data());
        title = histname + ";#it{p}_{T,jet}^{lead} - #it{A}_{jet}#rho (GeV/#it{c});counts";
        fHistManager.CreateTH1(histname.Data(), title.Data(), nCorrPtBins, minCorrPt, fMaxPt);
      }
    }
  }

  TIter nextElement(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = nextElement())) fOutput->Add(obj);
  PostData(1, fOutput);
}

/**
 * Perform steps needed to initialize the analysis.
 */
void AliAnalysisTaskJetUEStudies::ExecOnce()
{
  AliAnalysisTaskEmcalJetLight::ExecOnce();

  for (auto rhoIt = fAlternativeRho.begin(); rhoIt != fAlternativeRho.end();) {
    if (rhoIt->second) continue;
    rhoIt->second = dynamic_cast<AliRhoParameter*>(fInputEvent->FindListObject(rhoIt->first));
    if (rhoIt->second) {
      rhoIt++;
    }
    else {
      AliError(Form("%s: Could not retrieve rho %s! This rho name will be ignored", GetName(), rhoIt->first.Data()));
      rhoIt = fAlternativeRho.erase(rhoIt);
    }
  }

  fRandom = new TRandom3(0);
}

/**
 * Run the analysis.
 */
Bool_t AliAnalysisTaskJetUEStudies::Run()
{
  CalculateEventProperties();

  return kTRUE;
}

/**
 * Overloads base class method. Fills the output histograms
 * @return kTRUE if successful
 */
Bool_t AliAnalysisTaskJetUEStudies::FillHistograms()
{
  TString histname;

  if (fCentBin < 0) {
    AliError(Form("fCentBin is %d! fCent = %.3f. Fix the centrality bins to include all possible values of centrality.", fCentBin, fCent));
    return kFALSE;
  }

  for (auto rho1 : fAlternativeRho) {
    for (auto rho2 : fAlternativeRho) {
      if (rho1 == rho2) continue;
      histname = TString::Format("%s/fHistDelta%sOverRhoVsCent", rho1.first.Data(), rho2.first.Data());
      fHistManager.FillTH2(histname, fCent, 2.0 * (rho1.second->GetVal() - rho2.second->GetVal()) / (rho1.second->GetVal() + rho2.second->GetVal()));

      histname = TString::Format("%s/fHistDelta%sVsRho", rho1.first.Data(), rho2.first.Data());
      fHistManager.FillTH2(histname, (rho1.second->GetVal() + rho2.second->GetVal()) / 2, rho1.second->GetVal() - rho2.second->GetVal());

      histname = TString::Format("%s/fHist%sVs%s", rho1.first.Data(), rho2.first.Data(), rho1.first.Data());
      fHistManager.FillTH2(histname, rho1.second->GetVal(), rho2.second->GetVal());
    }
  }

  for (auto cont_it : fJetCollArray) {
    Bool_t isB2B = IsB2BEvent(cont_it.first);

    AliJetContainer* jets = cont_it.second;

    std::shared_ptr<AliEmcalJet> random_cone(GetRandomCone(jets));
    std::shared_ptr<AliEmcalJet> random_cone_perp(nullptr);
    std::shared_ptr<AliEmcalJet> random_cone_excl_lead(nullptr);
    AliEmcalJet* leadingJet = fLeadingJet[cont_it.first];

    if (leadingJet) {
      random_cone_perp.reset(GetRandomConePerp(jets, leadingJet));
      random_cone_excl_lead.reset(GetRandomConeExclLead(jets, leadingJet));
    }
    else {
      random_cone_perp = random_cone;
      random_cone_excl_lead = random_cone;
    }

    if (leadingJet) {
      histname = TString::Format("%s/fHistLeadingJetPtVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname, fCent, leadingJet->Pt());

      histname = TString::Format("%s/fHistLeadingJetPhiVsEta", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname, leadingJet->Eta(), leadingJet->Phi());

      if (isB2B){
        histname = TString::Format("%s/fHistB2BLeadingJetPtVsCent", jets->GetArrayName().Data());
        fHistManager.FillTH2(histname, fCent, leadingJet->Pt());
      }
    }

    for (auto rho : fAlternativeRho) {
      if (leadingJet) {
        histname = TString::Format("%s/%s/fHistLeadingJetCorrPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, fCent, leadingJet->Pt() - rho.second->GetVal() * leadingJet->Area());
      }

      histname = TString::Format("%s/%s/fHistRCPhiVsEta", jets->GetArrayName().Data(), rho.first.Data());
      fHistManager.FillTH2(histname, random_cone->Eta(), random_cone->Phi());

      histname = TString::Format("%s/%s/fHistRCDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      fHistManager.FillTH2(histname, fCent, random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

      histname = TString::Format("%s/%s/fHistRCDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      fHistManager.FillTH2(histname, rho.second->GetVal(), random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

      histname = TString::Format("%s/%s/fHistRCPerpPhiVsEta", jets->GetArrayName().Data(), rho.first.Data());
      fHistManager.FillTH2(histname, random_cone_perp->Eta(), random_cone_perp->Phi());

      histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
      fHistManager.FillTH2(histname, fCent, random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

      histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
      fHistManager.FillTH2(histname, rho.second->GetVal(), random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

      if (random_cone_excl_lead) {
        histname = TString::Format("%s/%s/fHistRCExclLeadJetPhiVsEta", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, random_cone_excl_lead->Eta(), random_cone_excl_lead->Phi());

        histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, fCent, random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());

        histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, rho.second->GetVal(), random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());
      }

      if (leadingJet) {
        histname = TString::Format("%s/%s/fHistRCDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, leadingJet->Pt(), random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

        histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, leadingJet->Pt(), random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

        if (random_cone_excl_lead) {
          histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, leadingJet->Pt(), random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());
        }
      }
      else {
        histname = TString::Format("%s/%s/fHistRCDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, -1, random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

        histname = TString::Format("%s/%s/fHistRCPerpDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, -1, random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

        histname = TString::Format("%s/%s/fHistRCExclLeadJetDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, -1, random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());
      }

      if (isB2B) {
        if (leadingJet) {
          histname = TString::Format("%s/%s/fHistB2BLeadingJetCorrPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, fCent, leadingJet->Pt() - rho.second->GetVal() * leadingJet->Area());
        }

        histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, fCent, random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

        histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, rho.second->GetVal(), random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

        histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, fCent, random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

        histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
        fHistManager.FillTH2(histname, rho.second->GetVal(), random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

        if (random_cone_excl_lead) {
          histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsCent", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, fCent, random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());

          histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsRho", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, rho.second->GetVal(), random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());
        }

        if (leadingJet) {
          histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, leadingJet->Pt(), random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

          histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, leadingJet->Pt(), random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

          if (random_cone_excl_lead) {
            histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
            fHistManager.FillTH2(histname, leadingJet->Pt(), random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());
          }
        }
        else {
          histname = TString::Format("%s/%s/fHistB2BRCDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, -1, random_cone->Pt() - rho.second->GetVal() * random_cone->Area());

          histname = TString::Format("%s/%s/fHistB2BRCPerpDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, -1, random_cone_perp->Pt() - rho.second->GetVal() * random_cone_perp->Area());

          histname = TString::Format("%s/%s/fHistB2BRCExclLeadJetDeltaPtVsLeadJetPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH2(histname, -1, random_cone_excl_lead->Pt() - rho.second->GetVal() * random_cone_excl_lead->Area());
        }

        for (auto rho2 : fAlternativeRho) {
          if (rho == rho2) continue;
          histname = TString::Format("%s/%s/fHistB2BDelta%sOverRhoVsCent", jets->GetArrayName().Data(), rho.first.Data(), rho2.first.Data());
          fHistManager.FillTH2(histname, fCent, 2.0 * (rho.second->GetVal() - rho2.second->GetVal()) / (rho.second->GetVal() + rho2.second->GetVal()));

          histname = TString::Format("%s/%s/fHistB2BDelta%sVsRho", jets->GetArrayName().Data(), rho.first.Data(), rho2.first.Data());
          fHistManager.FillTH2(histname, (rho.second->GetVal() + rho2.second->GetVal()) / 2, rho.second->GetVal() - rho2.second->GetVal());

          histname = TString::Format("%s/fHistB2B%sVs%s", rho.first.Data(), rho2.first.Data(), rho.first.Data());
          fHistManager.FillTH2(histname, rho.second->GetVal(), rho2.second->GetVal());
        }
      }
    }

    for (auto jet : jets->accepted()) {
      if (fCentBins.size() > 1) {
        histname = TString::Format("%s/fHistJetPt_Cent%d_%d", jets->GetArrayName().Data(), TMath::CeilNint(fCentBins[fCentBin]), TMath::CeilNint(fCentBins[fCentBin+1]));
        fHistManager.FillTH1(histname, jet->Pt());

        if (isB2B) {
          histname = TString::Format("%s/fHistB2BJetPt_Cent%d_%d", jets->GetArrayName().Data(), TMath::CeilNint(fCentBins[fCentBin]), TMath::CeilNint(fCentBins[fCentBin+1]));
          fHistManager.FillTH1(histname, jet->Pt());
        }

        for (auto rho : fAlternativeRho) {
          histname = TString::Format("%s/%s/fHistJetCorrPt_Cent%d_%d", jets->GetArrayName().Data(), rho.first.Data(), TMath::CeilNint(fCentBins[fCentBin]), TMath::CeilNint(fCentBins[fCentBin+1]));
          fHistManager.FillTH1(histname, jet->Pt() - jet->Area() * rho.second->GetVal());

          if (isB2B) {
            histname = TString::Format("%s/%s/fHistB2BJetCorrPt_Cent%d_%d", jets->GetArrayName().Data(), rho.first.Data(), TMath::CeilNint(fCentBins[fCentBin]), TMath::CeilNint(fCentBins[fCentBin+1]));
            fHistManager.FillTH1(histname, jet->Pt() - jet->Area() * rho.second->GetVal());
          }
        }
      }
      else {
        histname = TString::Format("%s/fHistJetPt", jets->GetArrayName().Data());
        fHistManager.FillTH1(histname, jet->Pt());

        if (isB2B) {
          histname = TString::Format("%s/fHistB2BJetPt", jets->GetArrayName().Data());
          fHistManager.FillTH1(histname, jet->Pt());
        }

        for (auto rho : fAlternativeRho) {
          histname = TString::Format("%s/%s/fHistJetCorrPt", jets->GetArrayName().Data(), rho.first.Data());
          fHistManager.FillTH1(histname, jet->Pt() - jet->Area() * rho.second->GetVal());

          if (isB2B) {
            histname = TString::Format("%s/%s/fHistB2BJetCorrPt", jets->GetArrayName().Data(), rho.first.Data());
            fHistManager.FillTH1(histname, jet->Pt() - jet->Area() * rho.second->GetVal());
          }
        }
      }
    }
  }

  return kTRUE;
}

/**
 * This function sums all the particles in a rigid cone. The particle are extracted from
 * a collection of AliEmcalContainer, i.e. multiple AliEmcalContainer are possible however they
 * must be all of the same type (AliParticleContainer, AliClusterContainer etc.).
 * @param[out] pt Sum of the transverse momenta of the particles
 * @param[in] eta Eta direction of the cone axis
 * @param[in] phi Phi direction of the cone axis
 * @param[in] maxD2 Radius squared of the con
 * @param[in] CollArray Map containing the array of AliEmcalContainer objects
 * @param[out] ConstList List of the particles that are included in the rigid cone
 * @return The number of particles included in the rigid cone
 */
template <class T, Int_t MAX_CONSTITUENTS>
Int_t AliAnalysisTaskJetUEStudies::SumParticles(Double_t& pt, Double_t eta, Double_t phi, Double_t maxD2, std::map<std::string, T*>& CollArray, std::array<Int_t, MAX_CONSTITUENTS>& ConstList)
{
  auto IndexMap = T::GetEmcalContainerIndexMap();
  Int_t N = 0;
  for (auto coll : CollArray) {
    auto itcont = coll.second->accepted_momentum();
    for (auto itmom = itcont.begin(); itmom != itcont.end(); itmom++) {
      Double_t etaDiff = eta - itmom->first.Eta();
      Double_t phiDiff = AliEmcalContainer::RelativePhi(phi, itmom->first.Phi());
      Double_t d2 = etaDiff*etaDiff + phiDiff*phiDiff;
      if (d2 > maxD2) continue;
      ConstList[N] = IndexMap.GlobalIndexFromLocalIndex(coll.second->GetArray(), itmom.current_index());
      pt += itmom->first.Pt();
      N++;
      if (N >= MAX_CONSTITUENTS) {
        AliError(Form("Reached the maximum number of constituents = %d", MAX_CONSTITUENTS));
        return N;
      }
    }
  }
  return N;
}

/**
 * This function generates a rigid cone jet, given a radius and a direction in the eta/phi plane
 * @param radius Radius of the rigid cone
 * @param eta Eta direction
 * @param phi Phi direction
 */
AliEmcalJet* AliAnalysisTaskJetUEStudies::GetJetCone(Double_t radius, Double_t eta, Double_t phi)
{
  Double_t pt = 0;
  Double_t maxD2 = radius * radius;

  const Int_t MAX_CONSTITUENTS = 100;

  static std::array<Int_t, MAX_CONSTITUENTS> tracks;
  static std::array<Int_t, MAX_CONSTITUENTS> clusters;

  Int_t nParticles = SumParticles<AliParticleContainer, MAX_CONSTITUENTS>(pt, eta, phi, maxD2, fParticleCollArray, tracks);
  Int_t nClusters = SumParticles<AliClusterContainer, MAX_CONSTITUENTS>(pt, eta, phi, maxD2, fClusterCollArray, clusters);

  AliEmcalJet* jet = new AliEmcalJet(pt, eta, phi, 0);
  jet->SetArea(maxD2*TMath::Pi());

  jet->SetNumberOfTracks(nParticles);
  jet->SetNumberOfClusters(nClusters);

  for (Int_t i = 0; i < nParticles; i++) jet->AddTrackAt(tracks[i], i);
  for (Int_t i = 0; i < nClusters; i++) jet->AddClusterAt(clusters[i], i);

  return jet;
}

/**
 * This function generates a random rigid cone jet using acceptance and radius from a given jet container
 * @param jetCont A jet container
 */
AliEmcalJet* AliAnalysisTaskJetUEStudies::GetRandomCone(AliJetContainer* jetCont)
{
  Double_t eta = fRandom->Uniform(jetCont->GetMinEta() + jetCont->GetJetRadius(), jetCont->GetMaxEta() - jetCont->GetJetRadius());
  Double_t phi = 0;
  if (jetCont->GetMinPhi() == -10 && jetCont->GetMaxPhi() == 10) {
    phi = fRandom->Uniform(0, TMath::TwoPi());
  }
  else {
    phi = fRandom->Uniform(jetCont->GetMinPhi() + jetCont->GetJetRadius(), jetCont->GetMaxPhi() - jetCont->GetJetRadius());
  }

  return GetJetCone(jetCont->GetJetRadius(), eta, phi);
}

/**
 * This function generates a random rigid cone jet using acceptance and radius from a given jet container.
 * The random direction is chosen in the azimuthal direction perpendicular to the jet provided.
 * @param jetCont A jet container
 * @param leadJet Leading jet
 */
AliEmcalJet* AliAnalysisTaskJetUEStudies::GetRandomConePerp(AliJetContainer* jetCont, AliEmcalJet* leadJet)
{
  Double_t eta = fRandom->Uniform(jetCont->GetMinEta() + jetCont->GetJetRadius(), jetCont->GetMaxEta() - jetCont->GetJetRadius());

  Double_t double_radius = jetCont->GetJetRadius()*2;

  Double_t minPhi = leadJet->Phi() + double_radius;
  Double_t maxPhi = leadJet->Phi() + TMath::Pi() - double_radius;

  Double_t sign = fRandom->Rndm();
  if (sign > 0.5) {
    minPhi += TMath::Pi();
    maxPhi += TMath::Pi();
  }

  Double_t phi = TVector2::Phi_0_2pi(fRandom->Uniform(minPhi, maxPhi));

  return GetJetCone(jetCont->GetJetRadius(), eta, phi);
}

/**
 * This function generates a random rigid cone jet using acceptance and radius from a given jet container.
 * The randon rigid cone is force to have no overlap with the jet provided.
 * @param jetCont A jet container
 * @param leadJet Leading jet
 */
AliEmcalJet* AliAnalysisTaskJetUEStudies::GetRandomConeExclLead(AliJetContainer* jetCont, AliEmcalJet* leadJet)
{
  Double_t eta = fRandom->Uniform(jetCont->GetMinEta() + jetCont->GetJetRadius(), jetCont->GetMaxEta() - jetCont->GetJetRadius());
  Double_t phi = 0;
  if (jetCont->GetMinPhi() == -10 && jetCont->GetMaxPhi() == 10) {
    phi = fRandom->Uniform(0, TMath::TwoPi());
  }
  else {
    phi = fRandom->Uniform(jetCont->GetMinPhi() + jetCont->GetJetRadius(), jetCont->GetMaxPhi() - jetCont->GetJetRadius());
  }

  const Int_t MAX_ITER = 20;
  AliEmcalJet* jet = nullptr;
  for (Int_t i = 0; i < MAX_ITER; i++) {
    jet = GetJetCone(jetCont->GetJetRadius(), eta, phi);
    if (!AreJetsOverlapping(leadJet, jet)) {
      break;
    }
    else {
      delete jet;
      jet = nullptr;
    }
  }

  if (!jet) AliError(Form("Could not find a random cone that does not overlap with the leading jet after %d iterations!", MAX_ITER));

  return jet;
}

/// Create an instance of this class and add it to the analysis manager
///
/// \param ntracks name of the track collection
/// \param nclusters name of the calorimeter cluster collection
/// \param trackPtCut minimum transverse momentum of tracks
/// \param clusECut minimum energy of calorimeter clusters
/// \param suffix additional suffix that can be added at the end of the task name
/// \return pointer to the new AddTaskEmcalJetSpectraQA task
AliAnalysisTaskJetUEStudies* AliAnalysisTaskJetUEStudies::AddTaskJetUEStudies(TString trackName, TString clusName,
    Double_t trackPtCut, Double_t clusECut, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskJetUEStudies", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AliAnalysisTaskJetUEStudies", "This task requires an input event handler");
    return nullptr;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString name("AliAnalysisTaskJetUEStudies");
  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskJetUEStudies* jetTask = new AliAnalysisTaskJetUEStudies(name);
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);
  AliParticleContainer *partCont = jetTask->AddParticleContainer(trackName.Data());
  if (partCont) partCont->SetParticlePtCut(trackPtCut);

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(clusName.Data());
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  // Final settings, pass to manager and set the containers
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                  TList::Class(),AliAnalysisManager::kOutputContainer,
                  Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );

  return jetTask;
}
