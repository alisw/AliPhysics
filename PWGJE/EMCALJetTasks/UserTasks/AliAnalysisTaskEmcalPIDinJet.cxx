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

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliEmcalCorrectionTask.h"
#include "THistManager.h"

#include "AliAnalysisTaskEmcalPIDinJet.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalPIDinJet);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalPIDinJet::AliAnalysisTaskEmcalPIDinJet() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager(),
  fPIDResponse(0)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalPIDinJet::AliAnalysisTaskEmcalPIDinJet(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),
  fPIDResponse(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalPIDinJet::~AliAnalysisTaskEmcalPIDinJet()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalPIDinJet::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateTrackHistograms();
  AllocateJetHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskEmcalPIDinJet::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    // for (Int_t cent = 0; cent < fNcentBins; cent++) {
  
      if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
        // CreateTHnSparse array
        Int_t nbin3[3] = {150, 36, 1200};
        Double_t max3[3] = {15, 0.9, 15};
        Double_t min3[3] = {0, -0.9, -15};

        // ITS, TPC and TOF sigmal
        histname = TString::Format("%s/histTrackTPC", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});dE/dx (a.u.)", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 4000, 0.1, 20, 4000, 20, 1000);
        histname = TString::Format("%s/histTrackITS", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});dE/dx (a.u.)", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 1000, 0.1, 20, 1000, 20, 1000);
        histname = TString::Format("%s/histTrackTOF", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});TOF #beta", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 1000, 0.3, 20, 800, 0.4, 1.2);
        // ITS, TPC and TOF signal after ITS+TPC+TOF Cut
        histname = TString::Format("%s/histTrackTPC2", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});dE/dx (a.u.)", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 4000, 0.1, 20, 4000, 20, 1000);
        histname = TString::Format("%s/histTrackITS2", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});dE/dx (a.u.)", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 1000, 0.1, 20, 1000, 20, 1000);
        histname = TString::Format("%s/histTrackTOF2", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});dE/dx (a.u.)", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 1000, 0.3, 20, 800, 0.4, 1.2);
        // TPC signal after TPC Cut
        histname = TString::Format("%s/histTrackTPC3", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{track} (GeV/#it{c});dE/dx (a.u.)", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 4000, 0.1, 20, 4000, 20, 1000);

        // TPC vs TOF (with ITS cut)
        histname = TString::Format("%s/histTrack2DPionNsigma2426", groupname.Data());
        histtitle = TString::Format("%s;N#sigma_{TPC};N#sigma_{TOF}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 100, -10, 10, 100, -10, 10);
        histname = TString::Format("%s/histTrack2DPionNsigma2628", groupname.Data());
        histtitle = TString::Format("%s;N#sigma_{TPC};N#sigma_{TOF}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 100, -10, 10, 100, -10, 10);
        histname = TString::Format("%s/histTrack2DPionNsigma2830", groupname.Data());
        histtitle = TString::Format("%s;N#sigma_{TPC};N#sigma_{TOF}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 100, -10, 10, 100, -10, 10);
        histname = TString::Format("%s/histTrack2DProtonNsigma2426", groupname.Data());
        histtitle = TString::Format("%s;N#sigma_{TPC};N#sigma_{TOF}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 100, -10, 10, 100, -10, 10);
        histname = TString::Format("%s/histTrack2DProtonNsigma2628", groupname.Data());
        histtitle = TString::Format("%s;N#sigma_{TPC};N#sigma_{TOF}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 100, -10, 10, 100, -10, 10);
        histname = TString::Format("%s/histTrack2DProtonNsigma2830", groupname.Data());
        histtitle = TString::Format("%s;N#sigma_{TPC};N#sigma_{TOF}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 100, -10, 10, 100, -10, 10);
        
        // Nsigma vs Pt about Pion and Proton with ITS+TPC+TOF Cut
        histname = TString::Format("%s/histTrackPionNsigPt", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);
        histname = TString::Format("%s/histTrackProtonNsigPt", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);

        histname = TString::Format("%s/histTrackPionNsigPtTOF", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);
        histname = TString::Format("%s/histTrackProtonNsigPtTOF", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);

        histname = TString::Format("%s/histTrackPionNsigPtITS", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);
        histname = TString::Format("%s/histTrackProtonNsigPtITS", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);

        // Pion and Proton Pt distribution with ITS+TPC+TOF Cut
        histname = TString::Format("%s/histTrackPionPt", groupname.Data());
        histtitle = TString::Format("%s;#eta ;#it{p}_{T,track} (GeV/#it{c})", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 36, -0.9, 0.9, 200, 0, 20);
        histname = TString::Format("%s/histTrackProtonPt", groupname.Data());
        histtitle = TString::Format("%s;#eta ;#it{p}_{T,track} (GeV/#it{c})", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 36, -0.9, 0.9, 200, 0, 20);

        // Pion and Proton Pt distribution with TPC Cut
        histname = TString::Format("%s/histTrack2PionPt", groupname.Data());
        histtitle = TString::Format("%s;#it{#eta};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 36, -0.9, 0.9, 200, 0, 20);
        histname = TString::Format("%s/histTrack2ProtonPt", groupname.Data());
        histtitle = TString::Format("%s;#it{#eta};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 36, -0.9, 0.9, 200, 0, 20);

        // Nsigma vs Eta vs Pt about pion and Proton with TPC
        histname = TString::Format("%s/histTrackPionNsigEta", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);
        histname = TString::Format("%s/histTrackProtonNsigEta", groupname.Data());
        histtitle = TString::Format("%s", histname.Data());
        fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);

        // track QA
        histname = TString::Format("%s/histTrackPt", groupname.Data());
        histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

        histname = TString::Format("%s/histTrackPhi", groupname.Data());
        histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

        histname = TString::Format("%s/histTrackEta", groupname.Data());
        histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

        histname = TString::Format("%s/histTrackEtaPhi", groupname.Data());
        histtitle = TString::Format("%s;#it{#eta}_{track};{#phi}_{track}", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 6, -1, 1, fNbins / 2, 0, TMath::TwoPi());
      }

      histname = TString::Format("%s/histNTracks", groupname.Data());
      histtitle = TString::Format("%s;number of tracks;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    // }
  }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalPIDinJet::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    jetCont->SetJetPtCut(1);
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    // for (Int_t cent = 0; cent < fNcentBins; cent++) {

      // Jet QA  *********************************************************************************************      
      histname = TString::Format("%s/histJetRhoValue", groupname.Data());
      histtitle = TString::Format("%s;Centrality ; #rho", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 50, 0, 50, 500, 0, 500);

      histname = TString::Format("%s/histJetContainer", groupname.Data());
      histtitle = TString::Format("%s;cotainer;count", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 10, 0, 10);

      histname = TString::Format("%s/histJetPt", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histJetCorrPt", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 250, -50, 200);

      histname = TString::Format("%s/histJetvsJetCorr", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,jet}^{corr} (GeV/#it{c})", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0, 200, 250, -50, 200);

      histname = TString::Format("%s/histJetPhi", groupname.Data());
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta", groupname.Data());
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histJetEtaPhi", groupname.Data());
      histtitle = TString::Format("%s;#it{#eta}_{jet};#phi_{jet}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, fNbins / 6, -1, 1, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetPtvsArea", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});#it{A}_{jet}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 1);

      // Get Leading and Subleading Jet  **********************************************************************
      histname = TString::Format("%s/histLeadingJetPt", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,leadingjet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histSubLeadingJetPt", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,subleadingjet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);
    
      histname = TString::Format("%s/histLeSubPhivsLeadPt", groupname.Data());
      histtitle = TString::Format("%s;#Delta#it{#phi}_{jet}; #it{p}_{T,jet}^{leading}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 100, 0, TMath::TwoPi(), 300, 0, 300);

      histname = TString::Format("%s/histLeSubPhivsdE", groupname.Data());
      histtitle = TString::Format("%s;#Delta#it{#phi}_{jet}; LeSub #Delta#it{p}_{T,jet}^{corr}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 100, 0, TMath::TwoPi(), 250, 0, 250);

      histname = TString::Format("%s/histNumberOfJet", groupname.Data());
      histtitle = TString::Format("%s;Number of Jet ;counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);

      // Do Track Loop inside(outside) Jet  **********************************************************************
      // jet track Pt ***************************************************************************
      histname = TString::Format("%s/histJetTrackCorrPt", groupname.Data());
      histtitle = TString::Format("%s;#Delta#eta; #Delta#phi", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 100, -0.7, 0.7, 100, -TMath::Pi(), TMath::Pi());

      histname = TString::Format("%s/histJetTrackCorr", groupname.Data());
      histtitle = TString::Format("%s;#Delta#eta; #Delta#phi", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 100, -0.7, 0.7, 100, -TMath::Pi(), TMath::Pi());

      // PPiRatio at each ditance R ********************************************************************************
      // CreateTHnSparse array
      Int_t nbin[4] = {19, 70, 100, 36};
      Double_t max[4] = {200, 3.5, 20, 0.9};
      Double_t min[4] = {10, 0, 0, -0.9};

      Int_t nbin4[5] = {19, 14, 25, 36, 150};
      Double_t max4[5] = {200, 0.7, 8, 0.9, 15};
      Double_t min4[5] = {10, 0, 3, -0.9, -15};

      histname = TString::Format("%s/histPionPt", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 4, nbin, min, max);

      histname = TString::Format("%s/histProtonPt", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 4, nbin, min, max);
      
      histname = TString::Format("%s/histPionPtBG1", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 4, nbin, min, max);

      histname = TString::Format("%s/histProtonPtBG1", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 4, nbin, min, max);

      histname = TString::Format("%s/histPionPtBG2", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 4, nbin, min, max);

      histname = TString::Format("%s/histProtonPtBG2", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 4, nbin, min, max);

      histname = TString::Format("%s/histJetNsigmaPion", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 5, nbin4, min4, max4);

      histname = TString::Format("%s/histJetNsigmaProton", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 5, nbin4, min4, max4);

      // Jet Rho *****************************************************************************************
      // CreateTHnSparse array
      Int_t nbin2[3] = {19, 14, 200};
      Double_t max2[3] = {200, 0.7, 100};
      Double_t min2[3] = {10, 0, 0};

      Int_t nbin3[3] = {190, 14, 300};
      Double_t max3[3] = {200, 0.7, 300};
      Double_t min3[3] = {10, 0, 0};

      histname = TString::Format("%s/histJetRhoEach", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});Distance #it{r}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 300, 0, 300, 14, 0, 0.7);

      histname = TString::Format("%s/histRhoEachBG1", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});Distance #it{r}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 300, 0, 300, 14, 0, 0.7);

      histname = TString::Format("%s/histRhoEachBG2", groupname.Data());
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});Distance #it{r}", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 300, 0, 300, 14, 0, 0.7);

      histname = TString::Format("%s/histRho", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin2, min2, max2);

      histname = TString::Format("%s/histRhoBG1", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin2, min2, max2);

      histname = TString::Format("%s/histRhoBG2", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin2, min2, max2);
      // ============================================================
      histname = TString::Format("%s/histJetTrkR", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);
      histname = TString::Format("%s/histJetTrkRBG1", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);
      histname = TString::Format("%s/histJetTrkRBG2", groupname.Data());
      histtitle = TString::Format("%s", histname.Data());
      fHistManager.CreateTHnSparse(histname, histtitle, 3, nbin3, min3, max3);

    // }
  }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalPIDinJet::FillHistograms()
{
  
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
  }

  DoJetLoop();
  DoTrackLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalPIDinJet::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  Double_t container = 0;
  
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    container++;
    groupname = jetCont->GetName();
    jetCont->SetJetPtCut(10);
    UInt_t count2 = 0;
    Double_t pt1 = 0, eta1 = 0, phi1 = 0,  pt2 = 0, ptr2 = 0, eta2 = 0, phi2 = 0;
    Double_t Rho = 0;

    if(fCent >= 10) continue;
    
    if (jetCont->GetRhoParameter()) {
      Rho = jetCont->GetRhoVal();
      histname = TString::Format("%s/histJetRhoValue", groupname.Data());
      fHistManager.FillTH2(histname, fCent, Rho);
    }
    if (Rho == 0) continue;

    histname = TString::Format("%s/histJetContainer", groupname.Data());
    fHistManager.FillTH1(histname, container);

    for(auto jet : jetCont->accepted()) {
      if (jet->GetNumberOfTracks() < 2) continue;
      if (jet->GetLeadingTrack()->Pt() < 5) continue;
      Double_t PtCorr = 0;

      // Rho Value in Jet Area ****************************
      PtCorr = jet->Pt() - Rho*(jet->Area());

      // Jet QA *******************************************************************************************
      histname = TString::Format("%s/histJetPt", groupname.Data());
      fHistManager.FillTH1(histname, jet->Pt());
      histname = TString::Format("%s/histJetCorrPt", groupname.Data());
      fHistManager.FillTH1(histname, PtCorr);
      histname = TString::Format("%s/histJetvsJetCorr", groupname.Data());
      fHistManager.FillTH2(histname, jet->Pt(), PtCorr);
      histname = TString::Format("%s/histJetPhi", groupname.Data());
      fHistManager.FillTH1(histname, jet->Phi());
      histname = TString::Format("%s/histJetEta", groupname.Data());
      fHistManager.FillTH1(histname, jet->Eta());
      histname = TString::Format("%s/histJetEtaPhi", groupname.Data());
      fHistManager.FillTH2(histname, jet->Eta() , jet->Phi());
      histname = TString::Format("%s/histJetPtvsArea", groupname.Data());
      fHistManager.FillTH2(histname, PtCorr, jet->Area());

      // Get Leading and Subleading Jet  **********************************************************************
      if (PtCorr < 10) continue;
      count2++;
      if (PtCorr > pt1) {
        pt2 = pt1;
        eta2 = eta1;
        phi2 = phi1;

        pt1 = PtCorr;
        eta1 = jet->Eta();
        phi1 = jet->Phi();
      } else {
        if (PtCorr > pt2) {
          pt2 = PtCorr;
          eta2 = jet->Eta();
          phi2 = jet->Phi();
        }
      }
    }

    Double_t DPhi = TMath::Pi()-std::abs(std::abs(phi2-phi1)-TMath::Pi());
    if (pt1 < 10) continue;
    if (pt2 == 0) {
      histname = TString::Format("%s/histLeadingJetPt", groupname.Data());
      fHistManager.FillTH1(histname, pt1);
    } else if (DPhi >= TMath::Pi()*5.0/6.0) {
      histname = TString::Format("%s/histLeadingJetPt", groupname.Data());
      fHistManager.FillTH1(histname, pt1);
      histname = TString::Format("%s/histSubLeadingJetPt", groupname.Data());
      fHistManager.FillTH1(histname, pt2);
      histname = TString::Format("%s/histLeSubPhivsLeadPt", groupname.Data());
      fHistManager.FillTH2(histname, DPhi, pt1);
      histname = TString::Format("%s/histLeSubPhivsdE", groupname.Data());
      fHistManager.FillTH2(histname, DPhi, pt1-pt2);
    }
    
    // Number of Jet ****************************
    histname = TString::Format("%s/histNumberOfJet", groupname.Data());
    fHistManager.FillTH1(histname, count2);

    // Do Track Loop inside(outside) Jet  **********************************************************************
    AliParticleContainer* partCont2 = 0;
    UInt_t count3 = 0;
    Double_t Ptsum[14] ={0}, PtsumBG1[14] = {0}, PtsumBG2[14] = {0};
    Double_t phi3 = phi1+(TMath::Pi()/2);
    Double_t phi4 = phi1-(TMath::Pi()/2);

    TIter next(&fParticleCollArray);
    while ((partCont2 = static_cast<AliParticleContainer*>(next()))) {
      for(auto part2 : partCont2->accepted()) {
        if (!part2) continue;

        Double_t preDeltaPhi1 = 0, DeltaPhi1 = 0, DeltaEta1 = 0, DeltaR1 = 1;
        Double_t preDeltaPhi3 = 0, DeltaPhi3 = 0, DeltaEta3 = 0, DeltaR3 = 1;
        Double_t preDeltaPhi4 = 0, DeltaPhi4 = 0, DeltaEta4 = 0, DeltaR4 = 1;
        Double_t nSigmaPion_TPC = 100, nSigmaProton_TPC = 100, nSigmaPion_TOF = 100, nSigmaProton_TOF = 100, nSigmaPion_ITS = 100, nSigmaProton_ITS = 100;
        Int_t Pion = 0, Proton = 0;

        const Double_t PiP1 = 146.4426245, PiP2 = 0.044101702, PiP3 = 0.063117611, PiP4 = 17.67452851, PiP5 = 3.049185967, PiP6 = -130.7208566;
        const Double_t PiB1 = 176.3921758, PiB2 = 0.039143447, PiB3 = 0.5110546, PiB4 = 11.87755816, PiB5 = 3.125664511, PiB6 = -159.3522847;

        if (partCont2->GetLoadedClass()->InheritsFrom("AliVTrack")) {
          const AliVTrack* track = static_cast<const AliVTrack*>(part2);

          // calculation Delta R from Jet Axis  *****************************************************************
          preDeltaPhi1 = track->Phi() - phi1;
          DeltaEta1 = track->Eta() - eta1;
          if (preDeltaPhi1 > TMath::Pi()) {
            DeltaPhi1 = preDeltaPhi1 - TMath::TwoPi();
          } else if (preDeltaPhi1 < -TMath::Pi()) {
            DeltaPhi1 = preDeltaPhi1 + TMath::TwoPi();
          } else {
            DeltaPhi1 = preDeltaPhi1;
          }
          DeltaR1 = std::sqrt(std::pow(DeltaEta1, 2.) + std::pow(DeltaPhi1, 2.));
          
          histname = TString::Format("%s/histJetTrackCorrPt", groupname.Data());
          fHistManager.FillTH2(histname, DeltaEta1, DeltaPhi1, track->Pt());
          histname = TString::Format("%s/histJetTrackCorr", groupname.Data());
          fHistManager.FillTH2(histname, DeltaEta1, DeltaPhi1);

          preDeltaPhi3 = track->Phi() - phi3;
          DeltaEta3 = track->Eta() - eta1;
          if (preDeltaPhi3 > TMath::Pi()) {
            DeltaPhi3 = preDeltaPhi3 - TMath::TwoPi();
          } else if (preDeltaPhi3 < -TMath::Pi()) {
            DeltaPhi3 = preDeltaPhi3 + TMath::TwoPi();
          } else {
            DeltaPhi3 = preDeltaPhi3;
          }
          preDeltaPhi4 = track->Phi() - phi4;
          DeltaEta4 = track->Eta() - eta1;
          if (preDeltaPhi4 > TMath::Pi()) {
            DeltaPhi4 = preDeltaPhi4 - TMath::TwoPi();
          } else if (preDeltaPhi4 < -TMath::Pi()) {
            DeltaPhi4 = preDeltaPhi4 + TMath::TwoPi();
          } else {
            DeltaPhi4 = preDeltaPhi4;
          }

          DeltaR3 = std::sqrt(std::pow(DeltaEta3, 2.) + std::pow(DeltaPhi3, 2.));
          DeltaR4 = std::sqrt(std::pow(DeltaEta4, 2.) + std::pow(DeltaPhi4, 2.));

          // *****************************************************************************************************
          // **** R dependency of track **************************************************************************
          // *****************************************************************************************************

          // PID nsigma *************************************************** 
          Double_t theta = 2*atan(exp(-(track->Eta())));
          Double_t Path1 = 1/sin(theta);
          Double_t Pimain = PiP1/(pow(track->P(),PiP2)) + PiP3*(1-exp(-(track->P()-PiP4)/PiP5)) + PiP6;
          Double_t Pibase = PiB1/(pow(track->P(),PiB2)) + PiB3*(1-exp(-(track->P()-PiB4)/PiB5)) + PiB6;
          Double_t Pmain = 0.788063667245128/pow(track->P(), 0.827848392474029);
          Double_t Pbase = 0.948451434/pow(track->P(), 0.714745037812307);
          Double_t PiITSmain = 0.31324467;
          Double_t PiITSbase = -0.021591449*(track->P())*(track->P()) + 0.312802826*(track->P()) + 0.752130515;
          Double_t PITSmain = -0.028610726*(track->P())*(track->P()) + 0.181788188*(track->P()) + 0.056034205;
          Double_t PITSbase = -0.039332093*(track->P())*(track->P()) + 0.30964846*(track->P()) + 0.352195793;

          nSigmaPion_TPC = (fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion))-(Pimain*Path1-Pibase);
          nSigmaProton_TPC = (fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton))-(Pmain*Path1-Pbase);
          nSigmaPion_TOF = (fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion))-(0.072720624*(track->P())-0.072833382);
          nSigmaProton_TOF = (fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton))-(0.020667827*(track->P())+0.157612143);
          nSigmaPion_ITS = (fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion))-(PiITSmain*Path1-PiITSbase);
          nSigmaProton_ITS = (fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton))-(PITSmain*Path1-PITSbase);

          if (track->Pt() < 3.0) {
            Pion = ((std::abs(nSigmaPion_ITS) < 2.5) && (std::abs(nSigmaPion_TPC) < 2.5) && (std::abs(nSigmaPion_TOF) < 2.5));
            Proton = ((std::abs(nSigmaProton_ITS) < 2.5) && (std::abs(nSigmaProton_TPC) < 2.5) && (std::abs(nSigmaProton_TOF) < 2.5));
          } else {
            Pion = (std::abs(nSigmaPion_TPC) < 2.5);
            Proton = (std::abs(nSigmaProton_TPC) < 2.5);
          }

          // cheak Nsigma distribution inside jet ****
          Double_t JetNsigPion[5] = {pt1, DeltaR1, track->Pt(), track->Eta(), nSigmaPion_TPC};
          histname = TString::Format("%s/histJetNsigmaPion", groupname.Data());
          fHistManager.FillTHnSparse(histname, JetNsigPion);

          Double_t JetNsigProton[5] = {pt1, DeltaR1, track->Pt(), track->Eta(), nSigmaProton_TPC};
          histname = TString::Format("%s/histJetNsigmaProton", groupname.Data());
          fHistManager.FillTHnSparse(histname, JetNsigProton);

          // Pt_Sum and pion/proton Pt in difference R *********************************************
          if(DeltaR1 < 0.05) Ptsum[0] += track->Pt();
          else if (DeltaR1 < 0.1) Ptsum[1] += track->Pt();
          else if (DeltaR1 < 0.15) Ptsum[2] += track->Pt();
          else if (DeltaR1 < 0.2) Ptsum[3] += track->Pt();
          else if (DeltaR1 < 0.25) Ptsum[4] += track->Pt();
          else if (DeltaR1 < 0.3) Ptsum[5] += track->Pt();
          else if (DeltaR1 < 0.35) Ptsum[6] += track->Pt();
          else if (DeltaR1 < 0.4) Ptsum[7] += track->Pt();
          else if (DeltaR1 < 0.45) Ptsum[8] += track->Pt();
          else if (DeltaR1 < 0.5) Ptsum[9] += track->Pt();
          else if (DeltaR1 < 0.55) Ptsum[10] += track->Pt();
          else if (DeltaR1 < 0.6) Ptsum[11] += track->Pt();
          else if (DeltaR1 < 0.65) Ptsum[12] += track->Pt();
          else if (DeltaR1 < 0.7) Ptsum[13] += track->Pt();

          histname = TString::Format("%s/histJetRhoEach", groupname.Data());
          fHistManager.FillTH2(histname, track->Pt(), DeltaR1, track->Pt()/(pt1*0.05));

          Double_t rhoCorr[3] = {pt1, DeltaR1, track->Pt()};
          histname = TString::Format("%s/histJetTrkR", groupname.Data());
          fHistManager.FillTHnSparse(histname, rhoCorr, track->Pt()/(pt1*0.05));
          
          if (Pion) {
            Double_t PionN[4] = {pt1, DeltaR1, track->Pt(), track->Eta()};
            histname = TString::Format("%s/histPionPt", groupname.Data());
            fHistManager.FillTHnSparse(histname, PionN);
          }
          if (Proton) {
            Double_t ProtonN[4] = {pt1, DeltaR1, track->Pt(), track->Eta()};
            histname = TString::Format("%s/histProtonPt", groupname.Data());
            fHistManager.FillTHnSparse(histname, ProtonN);
          }
        
          // Rho BackGround1 in difference R ******************************************************************
          if(DeltaR3 < 0.05) PtsumBG1[0] += track->Pt();
          else if (DeltaR3 < 0.1) PtsumBG1[1] += track->Pt();
          else if (DeltaR3 < 0.15) PtsumBG1[2] += track->Pt();
          else if (DeltaR3 < 0.2) PtsumBG1[3] += track->Pt();
          else if (DeltaR3 < 0.25) PtsumBG1[4] += track->Pt();
          else if (DeltaR3 < 0.3) PtsumBG1[5] += track->Pt();
          else if (DeltaR3 < 0.35) PtsumBG1[6] += track->Pt();
          else if (DeltaR3 < 0.4) PtsumBG1[7] += track->Pt();
          else if (DeltaR3 < 0.45) PtsumBG1[8] += track->Pt();
          else if (DeltaR3 < 0.5) PtsumBG1[9] += track->Pt();
          else if (DeltaR3 < 0.55) PtsumBG1[10] += track->Pt();
          else if (DeltaR3 < 0.6) PtsumBG1[11] += track->Pt();
          else if (DeltaR3 < 0.65) PtsumBG1[12] += track->Pt();
          else if (DeltaR3 < 0.7) PtsumBG1[13] += track->Pt();

          histname = TString::Format("%s/histRhoEachBG1", groupname.Data());
          fHistManager.FillTH2(histname, track->Pt(), DeltaR3, track->Pt()/(pt1*0.05));

          Double_t rhoCorrBG1[3] = {pt1, DeltaR3, track->Pt()};
          histname = TString::Format("%s/histJetTrkRBG1", groupname.Data());
          fHistManager.FillTHnSparse(histname, rhoCorrBG1, track->Pt()/(pt1*0.05));
          
          if (Pion) {
            Double_t PionN[4] = {pt1, DeltaR3, track->Pt(), track->Eta()};
            histname = TString::Format("%s/histPionPtBG1", groupname.Data());
            fHistManager.FillTHnSparse(histname, PionN);
          }
          if (Proton) {
            Double_t ProtonN[4] = {pt1, DeltaR3, track->Pt(), track->Eta()};
            histname = TString::Format("%s/histProtonPtBG1", groupname.Data());
            fHistManager.FillTHnSparse(histname, ProtonN);
          }
          
          // Rho BackGround2 in difference R ******************************************************************
          if(DeltaR4 < 0.05) PtsumBG2[0] += track->Pt();
          else if (DeltaR4 < 0.1) PtsumBG2[1] += track->Pt();
          else if (DeltaR4 < 0.15) PtsumBG2[2] += track->Pt();
          else if (DeltaR4 < 0.2) PtsumBG2[3] += track->Pt();
          else if (DeltaR4 < 0.25) PtsumBG2[4] += track->Pt();
          else if (DeltaR4 < 0.3) PtsumBG2[5] += track->Pt();
          else if (DeltaR4 < 0.35) PtsumBG2[6] += track->Pt();
          else if (DeltaR4 < 0.4) PtsumBG2[7] += track->Pt();
          else if (DeltaR4 < 0.45) PtsumBG2[8] += track->Pt();
          else if (DeltaR4 < 0.5) PtsumBG2[9] += track->Pt();
          else if (DeltaR4 < 0.55) PtsumBG2[10] += track->Pt();
          else if (DeltaR4 < 0.6) PtsumBG2[11] += track->Pt();
          else if (DeltaR4 < 0.65) PtsumBG2[12] += track->Pt();
          else if (DeltaR4 < 0.7) PtsumBG2[13] += track->Pt();

          histname = TString::Format("%s/histRhoEachBG2", groupname.Data());
          fHistManager.FillTH2(histname, track->Pt(), DeltaR4, track->Pt()/(pt1*0.05));

          Double_t rhoCorrBG2[3] = {pt1, DeltaR4, track->Pt()};
          histname = TString::Format("%s/histJetTrkRBG2", groupname.Data());
          fHistManager.FillTHnSparse(histname, rhoCorrBG2, track->Pt()/(pt1*0.05));
          
          if (Pion) {
            Double_t PionN[4] = {pt1, DeltaR4, track->Pt(), track->Eta()};
            histname = TString::Format("%s/histPionPtBG2", groupname.Data());
            fHistManager.FillTHnSparse(histname, PionN);
          }
          if (Proton) {
            Double_t ProtonN[4] = {pt1, DeltaR4, track->Pt(), track->Eta()};
            histname = TString::Format("%s/histProtonPtBG2", groupname.Data());
            fHistManager.FillTHnSparse(histname, ProtonN);
          }
        }
      }
    }

    for (Int_t i=0; i<14; i++){
      Double_t JetX[3] = {pt1, 0.05*i+0.025, Ptsum[i]/(0.05*pt1)};
      histname = TString::Format("%s/histRho", groupname.Data());
      fHistManager.FillTHnSparse(histname, JetX);
    }
    for (Int_t i=0; i<14; i++){
      Double_t JetXBG1[3] = {pt1, 0.05*i+0.025, PtsumBG1[i]/(0.05*pt1)};
      histname = TString::Format("%s/histRhoBG1", groupname.Data());
      fHistManager.FillTHnSparse(histname, JetXBG1);
    }
    for (Int_t i=0; i<14; i++){
      Double_t JetXBG2[3] = {pt1, 0.05*i+0.025, PtsumBG2[i]/(0.05*pt1)};
      histname = TString::Format("%s/histRhoBG2", groupname.Data());
      fHistManager.FillTHnSparse(histname, JetXBG2);
    }
  }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalPIDinJet::DoTrackLoop()
{
  AliClusterContainer* clusCont = GetClusterContainer(0);

  TString histname;
  TString groupname;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    if(fCent >= 10) continue;
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;
      Double_t nSigmaPion1_TPC = 100, nSigmaProton1_TPC = 100, nSigmaPion1_TOF = 100, nSigmaProton1_TOF = 100, nSigmaPion1_ITS = 100, nSigmaProton1_ITS = 100;
      Int_t Pion1 = 0, Proton1 = 0, Pion2 = 0, Proton2 = 0;

      // Pion Nsigma correction factor
      const Double_t PiP1 = 146.4426245, PiP2 = 0.044101702, PiP3 = 0.063117611, PiP4 = 17.67452851, PiP5 = 3.049185967, PiP6 = -130.7208566;
      const Double_t PiB1 = 176.3921758, PiB2 = 0.039143447, PiB3 = 0.5110546, PiB4 = 11.87755816, PiB5 = 3.125664511, PiB6 = -159.3522847;

      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track2 = static_cast<const AliVTrack*>(part);

        histname = TString::Format("%s/histTrackTPC", groupname.Data());
        fHistManager.FillTH1(histname, track2->P(), track2->GetTPCsignal());
        histname = TString::Format("%s/histTrackITS", groupname.Data());
        fHistManager.FillTH1(histname, track2->P(), track2->GetITSsignal());
        histname = TString::Format("%s/histTrackTOF", groupname.Data());
        fHistManager.FillTH1(histname, track2->P(), (100/3)*(track2->GetIntegratedLength())/(track2->GetTOFsignal()));

        Double_t theta = 2*atan(exp(-(track2->Eta())));
        Double_t Path1 = 1/sin(theta);
        Double_t Pimain = PiP1/(pow(track2->P(),PiP2)) + PiP3*(1-exp(-(track2->P()-PiP4)/PiP5)) + PiP6;
        Double_t Pibase = PiB1/(pow(track2->P(),PiB2)) + PiB3*(1-exp(-(track2->P()-PiB4)/PiB5)) + PiB6;
        Double_t Pmain = 0.788063667245128/pow(track2->P(), 0.827848392474029);
        Double_t Pbase = 0.948451434/pow(track2->P(), 0.714745037812307);
        Double_t PiITSmain = 0.31324467;
        Double_t PiITSbase = -0.021591449*(track2->P())*(track2->P()) + 0.312802826*(track2->P()) + 0.752130515;
        Double_t PITSmain = -0.028610726*(track2->P())*(track2->P()) + 0.181788188*(track2->P()) + 0.056034205;
        Double_t PITSbase = -0.039332093*(track2->P())*(track2->P()) + 0.30964846*(track2->P()) + 0.352195793;

        nSigmaPion1_TPC = (fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kPion))-(Pimain*Path1-Pibase);
        nSigmaProton1_TPC = (fPIDResponse->NumberOfSigmasTPC(track2,AliPID::kProton))-(Pmain*Path1-Pbase);
        nSigmaPion1_TOF = (fPIDResponse->NumberOfSigmasTOF(track2,AliPID::kPion))-(0.072720624*(track2->P())-0.072833382);
        nSigmaProton1_TOF = (fPIDResponse->NumberOfSigmasTOF(track2,AliPID::kProton))-(0.020667827*(track2->P())+0.107612143);
        nSigmaPion1_ITS = (fPIDResponse->NumberOfSigmasITS(track2,AliPID::kPion))-(PiITSmain*Path1-PiITSbase);
        nSigmaProton1_ITS = (fPIDResponse->NumberOfSigmasITS(track2,AliPID::kProton))-(PITSmain*Path1-PITSbase);

        if ( std::abs(nSigmaPion1_ITS) ) {
          if( track2->Pt() >= 2.4 && track2->Pt() < 2.6 ){
            histname = TString::Format("%s/histTrack2DPionNsigma2426", groupname.Data());
            fHistManager.FillTH2(histname, nSigmaPion1_TPC, nSigmaPion1_TOF);
          } else if( track2->Pt() >= 2.6 && track2->Pt() < 2.8 ){
            histname = TString::Format("%s/histTrack2DPionNsigma2628", groupname.Data());
            fHistManager.FillTH2(histname, nSigmaPion1_TPC, nSigmaPion1_TOF);
          } else if( track2->Pt() >= 2.8 && track2->Pt() < 3.0 ){
            histname = TString::Format("%s/histTrack2DPionNsigma2830", groupname.Data());
            fHistManager.FillTH2(histname, nSigmaPion1_TPC, nSigmaPion1_TOF);
          }
        }
        if ( std::abs(nSigmaProton1_ITS) ) {
          if( track2->Pt() >= 2.4 && track2->Pt() < 2.6 ){
            histname = TString::Format("%s/histTrack2DProtonNsigma2426", groupname.Data());
            fHistManager.FillTH2(histname, nSigmaProton1_TPC, nSigmaProton1_TOF);
          } else if( track2->Pt() >= 2.6 && track2->Pt() < 2.8 ){
            histname = TString::Format("%s/histTrack2DProtonNsigma2628", groupname.Data());
            fHistManager.FillTH2(histname, nSigmaProton1_TPC, nSigmaProton1_TOF);
          } else if( track2->Pt() >= 2.8 && track2->Pt() < 3.0 ){
            histname = TString::Format("%s/histTrack2DProtonNsigma2830", groupname.Data());
            fHistManager.FillTH2(histname, nSigmaProton1_TPC, nSigmaProton1_TOF);
          }
        }
        
        // Nsigma Cut ITS+TPC+TOF 
        Pion1 = (std::abs(nSigmaPion1_TPC) < 2.5 && std::abs(nSigmaPion1_TOF) < 2.5 && std::abs(nSigmaPion1_ITS) < 2.5);
        Proton1 = (std::abs(nSigmaProton1_TPC) < 2.5 && std::abs(nSigmaProton1_TOF) < 2.5 && std::abs(nSigmaProton1_ITS) < 2.5);
        // Nsigma Cut TPC 
        Pion2 = (std::abs(nSigmaPion1_TPC) < 2.5);
        Proton2 = (std::abs(nSigmaProton1_TPC) < 2.5);

        Double_t PionX[3] = {track2->Pt(), track2->Eta(), nSigmaPion1_TPC};
        Double_t ProtonX[3] = {track2->Pt(), track2->Eta(), nSigmaProton1_TPC};
        Double_t PionXTOF[3] = {track2->Pt(), track2->Eta(), nSigmaPion1_TOF};
        Double_t ProtonXTOF[3] = {track2->Pt(), track2->Eta(), nSigmaProton1_TOF};
        Double_t PionXITS[3] = {track2->Pt(), track2->Eta(), nSigmaPion1_ITS};
        Double_t ProtonXITS[3] = {track2->Pt(), track2->Eta(), nSigmaProton1_ITS};

        // ITS+TPC+TOF Cut
        if (Pion1 || Proton1) {
          histname = TString::Format("%s/histTrackTPC2", groupname.Data());
          fHistManager.FillTH2(histname, track2->P(), track2->GetTPCsignal());
          histname = TString::Format("%s/histTrackITS2", groupname.Data());
          fHistManager.FillTH2(histname, track2->P(), track2->GetITSsignal());
          histname = TString::Format("%s/histTrackTOF2", groupname.Data());
          fHistManager.FillTH2(histname, track2->P(), (100/3)*(track2->GetIntegratedLength())/(track2->GetTOFsignal()));

          if(Pion1){
            histname = TString::Format("%s/histTrackPionPt", groupname.Data());
            fHistManager.FillTH2(histname, track2->Eta(), track2->Pt());

            histname = TString::Format("%s/histTrackPionNsigPt", groupname.Data());
            fHistManager.FillTHnSparse(histname, PionX);
            histname = TString::Format("%s/histTrackPionNsigPtTOF", groupname.Data());
            fHistManager.FillTHnSparse(histname, PionXTOF);
            histname = TString::Format("%s/histTrackPionNsigPtITS", groupname.Data());
            fHistManager.FillTHnSparse(histname, PionXITS);
          }
          if(Proton1) {
            histname = TString::Format("%s/histTrackProtonPt", groupname.Data());
            fHistManager.FillTH2(histname, track2->Eta(), track2->Pt());
            
            histname = TString::Format("%s/histTrackProtonNsigPt", groupname.Data());
            fHistManager.FillTHnSparse(histname, ProtonX);
            histname = TString::Format("%s/histTrackProtonNsigPtTOF", groupname.Data());
            fHistManager.FillTHnSparse(histname, ProtonXTOF);
            histname = TString::Format("%s/histTrackProtonNsigPtITS", groupname.Data());
            fHistManager.FillTHnSparse(histname, ProtonXITS);
          }
        }

        if (Pion2 || Proton2) {
          histname = TString::Format("%s/histTrackTPC3", groupname.Data());
          fHistManager.FillTH2(histname, track2->P(), track2->GetTPCsignal());

          if(Pion2) {
            histname = TString::Format("%s/histTrack2PionPt", groupname.Data());
            fHistManager.FillTH2(histname, track2->Eta(), track2->Pt());
          }
          if(Proton2) {
            histname = TString::Format("%s/histTrack2ProtonPt", groupname.Data());
            fHistManager.FillTH2(histname, track2->Eta(), track2->Pt());
          }
        }

        // Nsigma vs Eta vs Pt
        histname = TString::Format("%s/histTrackPionNsigEta", groupname.Data());
        fHistManager.FillTHnSparse(histname, PionX);
        histname = TString::Format("%s/histTrackProtonNsigEta", groupname.Data());
        fHistManager.FillTHnSparse(histname, ProtonX);

        // track QA 
        histname = TString::Format("%s/histTrackPt", groupname.Data());
        fHistManager.FillTH1(histname, track2->Pt());

        histname = TString::Format("%s/histTrackPhi", groupname.Data());
        fHistManager.FillTH1(histname, track2->Phi());

        histname = TString::Format("%s/histTrackEta", groupname.Data());
        fHistManager.FillTH1(histname, track2->Eta());

        histname = TString::Format("%s/histTrackEtaPhi", groupname.Data());
        fHistManager.FillTH2(histname, track2->Eta() , track2->Phi());
      }
    }
    histname = TString::Format("%s/histNTracks", groupname.Data());
    fHistManager.FillTH1(histname, count);
  }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalPIDinJet::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalPIDinJet::Run()
{
 return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalPIDinJet::Terminate(Option_t *) 
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskEmcalPIDinJet * AliAnalysisTaskEmcalPIDinJet::AddTaskEmcalPIDinJet(
  const char *ntracks,
  const char *nclusters,
  const char* ncells,
  const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalPIDinJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalPIDinJet", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);

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

  if (cellName == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellName = "emcalCells";
    }
    else {
      cellName = "";
    }
  }

  TString name("AliAnalysisTaskEmcalPIDinJet");
  if (!trackName.IsNull()) {
    name += "_";
    name += trackName;
  }
  if (!clusName.IsNull()) {
    name += "_";
    name += clusName;
  }
  if (!cellName.IsNull()) {
    name += "_";
    name += cellName;
  }
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmcalPIDinJet* sampleTask = new AliAnalysisTaskEmcalPIDinJet(name);
  sampleTask->SetCaloCellsName(cellName);
  sampleTask->SetVzRange(-10,10);

  if (trackName == "mcparticles") {
    sampleTask->AddMCParticleContainer(trackName);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    sampleTask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    sampleTask->AddParticleContainer(trackName);
  }
  sampleTask->AddClusterContainer(clusName);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(sampleTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (sampleTask, 0,  cinput1 );
  mgr->ConnectOutput (sampleTask, 1, coutput1 );

  return sampleTask;
}
