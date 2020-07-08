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
#include <TF1.h>

#include <AliAnalysisManager.h>
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

#include "AliAnalysisTaskJetVn.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetVn);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskJetVn::AliAnalysisTaskJetVn() :
  AliAnalysisTaskEmcalJet(),
  fEventCuts(),
  fFlowQnVectorMgr(0x0),
  fHistManager()
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskJetVn::AliAnalysisTaskJetVn(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fEventCuts(),
  fFlowQnVectorMgr(0x0),
  fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskJetVn::~AliAnalysisTaskJetVn()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskJetVn::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>
    (AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (flowQnVectorTask != 0x0) {
      fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
  else {
    AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
  }

  AllocateEventHistograms();
  AllocateTrackHistograms();
  AllocateJetHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

void AliAnalysisTaskJetVn::AllocateEventHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  groupname="Event";
  fHistManager.CreateHistoGroup(groupname);

  histname  = TString::Format("%s/histCentrality", groupname.Data());
  histtitle = TString::Format("%s;Centrality;counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);

  histname  = TString::Format("%s/histEventPlanePsi1_V0A", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,V0A};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi1_V0C", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,V0C};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi1_TPC", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,TPC};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi2_V0A", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,V0A};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi2_V0C", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,V0C};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi2_TPC", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,TPC};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi3_V0A", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,V0A};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi3_V0C", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,V0C};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

  histname  = TString::Format("%s/histEventPlanePsi3_TPC", groupname.Data());
  histtitle = TString::Format("%s;#psi_{2,TPC};counts", histname.Data());
  fHistManager.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskJetVn::AllocateTrackHistograms()
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

        // adding histo for counting events
        histname = TString::Format("Hist nEvents");
        histtitle = TString::Format("Number of Events");
        fHistManager.CreateTH1(histname, histtitle, 1, 0.0, 1.0);

        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

            histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

            histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;number of tracks;events", histname.Data());
            if (fForceBeamType != kpp) {
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
            }

            // event plane angle histograms
            histname = TString::Format("%s/histPsi2VZERO_%d", groupname.Data(), cent);
            histtitle = "Psi2 from VZERO";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi3VZERO_%d", groupname.Data(), cent);
            histtitle = "Psi3 from VZERO";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi2VZEROA_%d", groupname.Data(), cent);
            histtitle = "Psi2 from VZEROA";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi3VZEROA_%d", groupname.Data(), cent);
            histtitle = "Psi3 from VZEROA";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi2VZEROC_%d", groupname.Data(), cent);
            histtitle = "Psi2 from VZEROC";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi3VZEROC_%d", groupname.Data(), cent);
            histtitle = "Psi3 from VZEROC";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi2TPC_%d", groupname.Data(), cent);
            histtitle = "Psi2 from TPC";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi3TPC_%d", groupname.Data(), cent);
            histtitle = "Psi3 from TPC";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi2TPCPositiveEta_%d", groupname.Data(), cent);
            histtitle = "Psi2 from TPC positive eta side";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi3TPCPositiveEta_%d", groupname.Data(), cent);
            histtitle = "Psi3 from TPC positive eta side";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi2TPCNegativeEta_%d", groupname.Data(), cent);
            histtitle = "Psi2 from TPC negative eta";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histPsi3TPCNegativeEta_%d", groupname.Data(), cent);
            histtitle = "Psi3 from TPC negative eta side";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            // histograms for particle angle relative to the event plane
            histname = TString::Format("%s/histTrackPhiMinusPsi2_%d", groupname.Data(), cent);
            histtitle = "Track phi minus Psi2";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histTrackPhiMinusPsi3_%d", groupname.Data(), cent);
            histtitle = "Track phi minus Psi3";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            // sub-event plane resolution histograms for calculating event plane resolution
            histname = TString::Format("%s/histSubEventPlaneResolution1_%d", groupname.Data(), cent);
            histtitle = "Sub event plane resolution 1";
            fHistManager.CreateTH1(histname, histtitle, 100, -1, 1);

            histname = TString::Format("%s/histSubEventPlaneResolution2_%d", groupname.Data(), cent);
            histtitle = "Sub event plane resolution 2";
            fHistManager.CreateTH1(histname, histtitle, 100, -1, 1);

            histname = TString::Format("%s/histSubEventPlaneResolution3_%d", groupname.Data(), cent);
            histtitle = "Sub event plane resolution 3";
            fHistManager.CreateTH1(histname, histtitle, 100, -1, 1);

            // event plane histos from Qn vector corrections framework
            histname = TString::Format("%s/histQnCorrPsi2TPCPosEta_%d", groupname.Data(), cent);
            histtitle = "Qn Vector Corrected Event Plane Psi2 in TPC Positive Eta";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histQnCorrPsi3TPCPosEta_%d", groupname.Data(), cent);
            histtitle = "Qn Vector Corrected Event Plane Psi3 in TPC Positive Eta";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histQnCorrPsi2TPCNegEta_%d", groupname.Data(), cent);
            histtitle = "Qn Vector Corrected Event Plane Psi2 in TPC Negative Eta";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histQnCorrPsi3TPCNegEta_%d", groupname.Data(), cent);
            histtitle = "Qn Vector Corrected Event Plane Psi3 in TPC Negative Eta";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histQnCorrPsi2VZERO_%d", groupname.Data(), cent);
            histtitle = "Qn Vector Corrected Event Plane Psi2 in VZERO";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());

            histname = TString::Format("%s/histQnCorrPsi3VZERO_%d", groupname.Data(), cent);
            histtitle = "Qn Vector Corrected Event Plane Psi3 in VZERO";
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, -TMath::TwoPi(), 2*TMath::TwoPi());
        }
    }

    histname = "fHistSumNTracks";
    histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
    if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
    }
    else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
    }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskJetVn::AllocateJetHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        //cout << "Groupname = " << groupname << "\n";
        // Protect against creating the histograms twice
        if (fHistManager.FindObject(groupname)) {
            AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
            continue;
        }
        fHistManager.CreateHistoGroup(groupname);
        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

            histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

            histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

            histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;number of jets;events", histname.Data());
            if (fForceBeamType != kpp) {
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
            }

            // histograms for jet angle relative to the event plane
            histname = TString::Format("%s/histJetPhiMinusPsi2_%d", groupname.Data(), cent);
            histtitle = "Jet phi minus psi2";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());
            histname = TString::Format("%s/histJetPhiMinusPsi3_%d", groupname.Data(), cent);
            histtitle = "Jet phi minus psi3";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            // histograms for in-plane vs out-of-plane jets
            histname = TString::Format("%s/histJetsInPlaneOutOfPlaneV2_%d", groupname.Data(), cent);
            histtitle = "In-plane vs out-of-plane jets (v2)";
            fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);
            histname = TString::Format("%s/histJetsInPlaneOutOfPlaneV3_%d", groupname.Data(), cent);
            histtitle = "In-plane vs out-of-plane jets (v2)";
            fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);

            if (!jetCont->GetRhoName().IsNull()) {
                // Rho histograms
                histname = TString::Format("%s/histJetRho_%d", groupname.Data(), cent);
                histtitle = "Rho";
                fHistManager.CreateTH1(histname, histtitle, fNbins, 0.0, 300.0);
                histname = TString::Format("%s/histJetRhoLocal_%d", groupname.Data(), cent);
                histtitle = "Rho Local";
                fHistManager.CreateTH1(histname, histtitle, fNbins, 0.0, 300.0);
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
                fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
                histname = TString::Format("%s/histJetCorrPtLocal_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} local (GeV/#it{c});counts", histname.Data());
                fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);

                // histograms for corrected in-plane vs out-of-plane jets
                histname = TString::Format("%s/histJetsCorrInPlaneOutOfPlaneV2_%d", groupname.Data(), cent);
                histtitle = "Corrected in-plane vs out-of-plane jets (v2)";
                fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);
                histname = TString::Format("%s/histJetsCorrInPlaneOutOfPlaneV3_%d", groupname.Data(), cent);
                histtitle = "Corrected in-plane vs out-of-plane jets (v2)";
                fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);
                histname = TString::Format("%s/histJetsCorrLocalInPlaneOutOfPlaneV2_%d", groupname.Data(), cent);
                histtitle = "Local rho corrected in-plane vs out-of-plane jets (v2)";
                fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);
                histname = TString::Format("%s/histJetsCorrLocalInPlaneOutOfPlaneV3_%d", groupname.Data(), cent);
                histtitle = "Local rho corrected in-plane vs out-of-plane jets (v2)";
                fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);

                histname = TString::Format("%s/histJetLocalRhoVsAverageRho_%d", groupname.Data(), cent);
                histtitle = "Local rho versus average rho";
                fHistManager.CreateTH2(histname, histtitle, fNbins, 0.0, 300.0, fNbins, 0.0, 300.0);
                histname = TString::Format("%s/histJetCorrPtLocalVsJetCorrPt_%d", groupname.Data(), cent);
                histtitle = "Local rho adjusted jet pT versus average rho adjusted jet pT";
                fHistManager.CreateTH2(histname, histtitle, fNbins, 0.0, 100.0, fNbins, 0.0, 100.0);

                // histo local rho vs delta phi
                histname = TString::Format("%s/histJetRhoVsDeltaPhi_%d", groupname.Data(), cent);
                histtitle = "Local rho versus angle relative to event plane";
                fHistManager.CreateTH2(histname, histtitle, fNbins, 0.0, TMath::TwoPi(), fNbins, 0.0, 300.0);
            }
        }
    }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskJetVn::FillHistograms()
{
  FillEventHistos();
  DoJetLoop();
  DoTrackLoop();
  //DoClusterLoop();
  //DoCellLoop();

  return kTRUE;
}

void AliAnalysisTaskJetVn::FillEventHistos()
{
  TString histname;
  TString groupname;
  groupname="Event";
  const AliQnCorrectionsQnVector *pQnVectorV0A  = 0x0;
  const AliQnCorrectionsQnVector *pQnVectorV0C  = 0x0;
  const AliQnCorrectionsQnVector *pQnVectorTPC  = 0x0;
  pQnVectorV0A  = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","align","align");
  pQnVectorV0C  = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","align","align");
  pQnVectorTPC  = fFlowQnVectorMgr->GetDetectorQnVector("TPC", "latest", "latest");

  histname = TString::Format("%s/histCentrality", groupname.Data());
  fHistManager.FillTH1(histname, fCent);

  if(pQnVectorV0A) {
    histname = TString::Format("%s/histEventPlanePsi1_V0A", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorV0A->EventPlane(1));

    histname = TString::Format("%s/histEventPlanePsi2_V0A", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorV0A->EventPlane(2));

    histname = TString::Format("%s/histEventPlanePsi3_V0A", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorV0A->EventPlane(3));
  }

  if(pQnVectorV0C) {
    histname = TString::Format("%s/histEventPlanePsi1_V0C", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorV0C->EventPlane(1));

    histname = TString::Format("%s/histEventPlanePsi2_V0C", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorV0C->EventPlane(2));

    histname = TString::Format("%s/histEventPlanePsi3_V0C", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorV0C->EventPlane(3));
  }


  if(pQnVectorTPC) {
    histname = TString::Format("%s/histEventPlanePsi1_TPC", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorTPC->EventPlane(1));

    histname = TString::Format("%s/histEventPlanePsi2_TPC", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorTPC->EventPlane(2));

    histname = TString::Format("%s/histEventPlanePsi3_TPC", groupname.Data());
    fHistManager.FillTH1(histname, pQnVectorTPC->EventPlane(3));
  }
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetVn::DoJetLoop()
{
    Double_t psi2, phiMinusPsi2; // angles for jet vn calculation
    Double_t psi3, phiMinusPsi3;
    //Double_t totalRhoVal = CalculateTotalRho(param0, param1);

    TString histname;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        UInt_t count = 0;

        // Qn Corrected Psi angles
        Double_t corrPsi2TPCPosEta = 0.0;
        Double_t corrPsi3TPCPosEta = 0.0;
        const AliQnCorrectionsQnVector *corrQnVectorTPCPosEta;
        corrQnVectorTPCPosEta = fFlowQnVectorMgr->GetDetectorQnVector("TPCPosEtaQoverM");
        if (corrQnVectorTPCPosEta != NULL) {
          corrPsi2TPCPosEta = corrQnVectorTPCPosEta->EventPlane(2);
          if (corrPsi2TPCPosEta < 0.0) corrPsi2TPCPosEta += TMath::TwoPi()/2;
          corrPsi3TPCPosEta = corrQnVectorTPCPosEta->EventPlane(3);
          if (corrPsi3TPCPosEta < 0.0) corrPsi3TPCPosEta += TMath::TwoPi()/3;
        }
        Double_t corrPsi2TPCNegEta = 0.0;
        Double_t corrPsi3TPCNegEta = 0.0;
        const AliQnCorrectionsQnVector *corrQnVectorTPCNegEta;
        corrQnVectorTPCNegEta = fFlowQnVectorMgr->GetDetectorQnVector("TPCNegEtaQoverM");
        if (corrQnVectorTPCNegEta != NULL) {
          corrPsi2TPCNegEta = corrQnVectorTPCNegEta->EventPlane(2);
          if (corrPsi2TPCNegEta < 0.0) corrPsi2TPCNegEta += TMath::TwoPi()/2;
          corrPsi3TPCNegEta = corrQnVectorTPCNegEta->EventPlane(3);
          if (corrPsi3TPCNegEta < 0.0) corrPsi3TPCNegEta += TMath::TwoPi()/3;
        }
        Double_t corrPsi2VZERO = 0.0;
        Double_t corrPsi3VZERO = 0.0;
        const AliQnCorrectionsQnVector *corrQnVectorVZERO;
        corrQnVectorVZERO = fFlowQnVectorMgr->GetDetectorQnVector("VZEROQoverM");
        if (corrQnVectorVZERO != NULL) {
          corrPsi2VZERO = corrQnVectorVZERO->EventPlane(2);
          if (corrPsi2VZERO < 0.0) corrPsi2VZERO += TMath::TwoPi()/2;
          corrPsi3VZERO = corrQnVectorVZERO->EventPlane(3);
          if (corrPsi3VZERO < 0.0) corrPsi3VZERO += TMath::TwoPi()/3;
        }

        for(auto jet : jetCont->accepted()) {
            if (!jet) continue;
            count++;

            histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Pt());
            histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Area());
            histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Phi());
            histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Eta());

            // Filling histos for angle relative to event plane
            histname = TString::Format("%s/histJetPhiMinusPsi2_%d", groupname.Data(), fCentBin);
            phiMinusPsi2 = jet->Phi() - corrPsi2VZERO;
            if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
            fHistManager.FillTH1(histname, phiMinusPsi2);
            histname = TString::Format("%s/histJetPhiMinusPsi3_%d", groupname.Data(), fCentBin);
            phiMinusPsi3 = jet->Phi() - corrPsi3VZERO;
            if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
            fHistManager.FillTH1(histname, phiMinusPsi3);

            //if (jetCont->GetRhoParameter()) {
                Double_t localRhoVal = 0.0;
                Double_t localRhoValScaled = 0.0;
                Double_t jetPtCorr = 0.0;
                Double_t jetPtCorrLocal = 0.0;
                if (jet->Phi() - corrPsi2VZERO >= 0.0) localRhoVal = CalculateLocalRho(jet->Phi() - corrPsi2VZERO, param0, param1);
                if (jet->Phi() - corrPsi2VZERO < 0.0) localRhoVal = CalculateLocalRho(jet->Phi() - corrPsi2VZERO + TMath::TwoPi(), param0, param1);
                //if (jet->Phi() - psi2 >= 0.0) localRhoVal = CalculateLocalRho(jet->Phi() - psi2, 1, param1);
                //if (jet->Phi() - psi2 < 0.0) localRhoVal = CalculateLocalRho(jet->Phi() - psi2 + TMath::TwoPi(), 1, param1);
                //localRhoValScaled = localRhoVal * jetCont->GetRhoVal() / param0;
                //localRhoValScaled = (jetCont->GetRhoVal() / (2*0.2*totalRhoVal)) * localRhoVal;
                localRhoValScaled = localRhoVal * jetCont->GetRhoVal() / param0;
                jetPtCorr = jet->Pt() - jetCont->GetRhoVal() * jet->Area();
                jetPtCorrLocal = jet->Pt() - localRhoValScaled * jet->Area();

                histname = TString::Format("%s/histJetRho_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, jetCont->GetRhoVal());
                histname = TString::Format("%s/histJetRhoLocal_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, localRhoVal); // trying out local rho val
                //fHistManager.FillTH1(histname, InputEvent()->FindListObject(fLocalRho)->GetLocalVal(jet->Phi(), 0.2)); // trying out local rho val
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, jetPtCorr);
                histname = TString::Format("%s/histJetCorrPtLocal_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, jetPtCorrLocal);
                //fHistManager.FillTH1(histname, jet->Pt() - localRhoValScaled * 0.4);

                //cout << "Jet radius = " << GetJetContainer()->GetJetRadius() << "; Jet area = " << jet->Area() << endl;
                //cout << "Average rho = " << jetCont->GetRhoVal() << "; Total rho = " << totalRhoVal << endl;
                //cout << "Average rho = " << jetCont->GetRhoVal() << endl;
                //cout << "Local rho = " << localRhoVal << "; Local rho scaled = " << localRhoValScaled << endl;
                //cout << "Jet pT = " << jet->Pt() << "; rho adj jet pT = " << jetPtCorr << "; Local rho adj jet pT = " << jetPtCorrLocal << endl;

                histname = TString::Format("%s/histJetLocalRhoVsAverageRho_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH2(histname, jetCont->GetRhoVal(), localRhoValScaled);
                histname = TString::Format("%s/histJetCorrPtLocalVsJetCorrPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH2(histname, jetPtCorr, jetPtCorrLocal);
                histname = TString::Format("%s/histJetRhoVsDeltaPhi_%d", groupname.Data(), fCentBin);
                if (jet->Phi() - corrPsi2VZERO >= 0.0) fHistManager.FillTH2(histname, jet->Phi() - corrPsi2VZERO, localRhoValScaled);
                if (jet->Phi() - corrPsi2VZERO < 0.0) fHistManager.FillTH2(histname, jet->Phi() - corrPsi2VZERO + TMath::TwoPi(), localRhoValScaled);
            //}

            if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)) {
                histname = TString::Format("%s/histJetsInPlaneOutOfPlaneV2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jet->Pt());
                histname = TString::Format("%s/histJetsCorrInPlaneOutOfPlaneV2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jetPtCorr);
                histname = TString::Format("%s/histJetsCorrLocalInPlaneOutOfPlaneV2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jetPtCorrLocal);
            }
            if ((phiMinusPsi2 >= TMath::Pi()/4 && phiMinusPsi2 < 3*TMath::Pi()/4) || (phiMinusPsi2 >= 5*TMath::Pi()/4 && phiMinusPsi2 < 7*TMath::Pi()/4)) {
                histname = TString::Format("%s/histJetsInPlaneOutOfPlaneV2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jet->Pt());
                histname = TString::Format("%s/histJetsCorrInPlaneOutOfPlaneV2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jetPtCorr);
                histname = TString::Format("%s/histJetsCorrLocalInPlaneOutOfPlaneV2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jetPtCorrLocal);
            }
            if ((phiMinusPsi3 < TMath::Pi()/6) || (phiMinusPsi3 >= TMath::Pi()/2 && phiMinusPsi3 < 5*TMath::Pi()/6) ||
                (phiMinusPsi3 >= 7*TMath::Pi()/6 && phiMinusPsi3 < 3*TMath::Pi()/2) || (phiMinusPsi3 >= 11*TMath::Pi()/6)) {
                histname = TString::Format("%s/histJetsInPlaneOutOfPlaneV3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jet->Pt());
                histname = TString::Format("%s/histJetsCorrInPlaneOutOfPlaneV3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jet->Pt()- jetCont->GetRhoVal() * jet->Area());
                histname = TString::Format("%s/histJetsCorrLocalInPlaneOutOfPlaneV3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jetPtCorrLocal);
            }
            if ((phiMinusPsi3 >= TMath::Pi()/6 && phiMinusPsi3 < TMath::Pi()/2) || (phiMinusPsi3 >= 5*TMath::Pi()/6 && phiMinusPsi3 < 7*TMath::Pi()/6) ||
                (phiMinusPsi3 >= 3*TMath::Pi()/2 && phiMinusPsi3 < 11*TMath::Pi()/6)) {
                histname = TString::Format("%s/histJetsInPlaneOutOfPlaneV3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jet->Pt());
                histname = TString::Format("%s/histJetsCorrInPlaneOutOfPlaneV3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jetPtCorr);
                histname = TString::Format("%s/histJetsCorrLocalInPlaneOutOfPlaneV3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jetPtCorrLocal);
            }
        }
        histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

void AliAnalysisTaskJetVn::DoTrackLoop()
{
    //AliClusterContainer* clusCont = GetClusterContainer(0);
    Double_t minEta = -0.7, maxEta = 0.7;
    Double_t psi2, phiMinusPsi2; // angles for jet vn calculation
    Double_t psi3, phiMinusPsi3;
    Double_t psi2VZERO, psi2VZEROA, psi2VZEROC;
    Double_t psi3VZERO, psi3VZEROA, psi3VZEROC;
    Double_t psi2TPC = -999.0, psi2TPCPositiveEta = -999.0, psi2TPCNegativeEta = -999.0;
    Double_t psi3TPC = -999.0, psi3TPCPositiveEta = -999.0, psi3TPCNegativeEta = -999.0;
    Double_t QxTPC = 0, QyTPC = 0, QxTPCPositiveEta = 0, QyTPCPositiveEta = 0, QxTPCNegativeEta = 0, QyTPCNegativeEta = 0;

    TString histname;
    TString groupname;
    UInt_t sumAcceptedTracks = 0;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = partCont->GetName();

        UInt_t count = 0;

        // counting events
        histname = TString::Format("Hist nEvents");
        fHistManager.FillTH1(histname, 0.5);

        // filling Qn Corrected Psi histograms
        Double_t corrPsi2TPCPosEta = 0.0;
        Double_t corrPsi3TPCPosEta = 0.0;
        const AliQnCorrectionsQnVector *corrQnVectorTPCPosEta;
        corrQnVectorTPCPosEta = fFlowQnVectorMgr->GetDetectorQnVector("TPCPosEtaQoverM");
        if (corrQnVectorTPCPosEta != NULL) {
          corrPsi2TPCPosEta = corrQnVectorTPCPosEta->EventPlane(2);
          if (corrPsi2TPCPosEta < 0.0) corrPsi2TPCPosEta += TMath::TwoPi()/2;
          corrPsi3TPCPosEta = corrQnVectorTPCPosEta->EventPlane(3);
          if (corrPsi3TPCPosEta < 0.0) corrPsi3TPCPosEta += TMath::TwoPi()/3;
          histname = TString::Format("%s/histQnCorrPsi2TPCPosEta_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, corrPsi2TPCPosEta);
          histname = TString::Format("%s/histQnCorrPsi3TPCPosEta_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, corrPsi3TPCPosEta);
        }
        Double_t corrPsi2TPCNegEta = 0.0;
        Double_t corrPsi3TPCNegEta = 0.0;
        const AliQnCorrectionsQnVector *corrQnVectorTPCNegEta;
        corrQnVectorTPCNegEta = fFlowQnVectorMgr->GetDetectorQnVector("TPCNegEtaQoverM");
        if (corrQnVectorTPCNegEta != NULL) {
          corrPsi2TPCNegEta = corrQnVectorTPCNegEta->EventPlane(2);
          if (corrPsi2TPCNegEta < 0.0) corrPsi2TPCNegEta += TMath::TwoPi()/2;
          corrPsi3TPCNegEta = corrQnVectorTPCNegEta->EventPlane(3);
          if (corrPsi3TPCNegEta < 0.0) corrPsi3TPCNegEta += TMath::TwoPi()/3;
          histname = TString::Format("%s/histQnCorrPsi2TPCNegEta_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, corrPsi2TPCNegEta);
          histname = TString::Format("%s/histQnCorrPsi3TPCNegEta_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, corrPsi3TPCNegEta);
        }
        Double_t corrPsi2VZERO = 0.0;
        Double_t corrPsi3VZERO = 0.0;
        const AliQnCorrectionsQnVector *corrQnVectorVZERO;
        corrQnVectorVZERO = fFlowQnVectorMgr->GetDetectorQnVector("VZEROQoverM");
        if (corrQnVectorVZERO != NULL) {
          corrPsi2VZERO = corrQnVectorVZERO->EventPlane(2);
          if (corrPsi2VZERO < 0.0) corrPsi2VZERO += TMath::TwoPi()/2;
          corrPsi3VZERO = corrQnVectorVZERO->EventPlane(3);
          if (corrPsi3VZERO < 0.0) corrPsi3VZERO += TMath::TwoPi()/3;
          histname = TString::Format("%s/histQnCorrPsi2VZERO_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, corrPsi2VZERO);
          histname = TString::Format("%s/histQnCorrPsi3VZERO_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, corrPsi3VZERO);
        }

        //psi2 = CalculateEventPlaneVZERO(2);
        //psi3 = CalculateEventPlaneVZERO(3);

        for(auto part : partCont->accepted()) {
            if (!part) continue;
            count++;

            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, part->Pt());

            histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, part->Phi());

            histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, part->Eta());

            if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
                const AliVTrack* track = static_cast<const AliVTrack*>(part);
                // Filling histos for angle relative to event plane
                phiMinusPsi2 = track->Phi() - corrPsi2VZERO;
                phiMinusPsi3 = track->Phi() - corrPsi3VZERO;
                if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
                if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
                histname = TString::Format("%s/histTrackPhiMinusPsi2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, phiMinusPsi2);
                histname = TString::Format("%s/histTrackPhiMinusPsi3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, phiMinusPsi3);

                // Calculating event plane using TPC
                //cout << "Particle eta = " << track->Eta() << ", phi = " << track->Phi() << endl;
                if (track->Eta() >= minEta && track->Eta() <= maxEta) {
                    QxTPC += TMath::Cos(2.0*track->Phi());
                    QyTPC += TMath::Sin(2.0*track->Phi());
                    //cout << "Particle inside TPC eta range. QxTPC = " << QxTPC << ", QyTPC = " << QyTPC << "." << endl;
                }
                if (track->Eta() >= 0.0 && track->Eta() <= maxEta) {
                    QxTPCPositiveEta += TMath::Cos(2.0*track->Phi());
                    QyTPCPositiveEta += TMath::Sin(2.0*track->Phi());
                    //cout << "Particle inside TPC positive eta range. QxTPCPositiveEta = " << QxTPCPositiveEta << ", QyTPCPositiveEta = " << QyTPCPositiveEta << "." << endl;
                }
                if (track->Eta() >= minEta && track->Eta() <= 0.0) {
                    QxTPCNegativeEta += TMath::Cos(2.0*track->Phi());
                    QyTPCNegativeEta += TMath::Sin(2.0*track->Phi());
                    //cout << "Particle inside TPC negative eta range. QxTPCNegativeEta = " << QxTPCNegativeEta << ", QyTPCNegativeEta = " << QyTPCNegativeEta << "." << endl;
                }
                if (track->Eta() > maxEta || track->Eta() < minEta) {
                    //cout << "Particle outside of TPC eta range." << endl;
                }
            }
        }
        sumAcceptedTracks += count;

        histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);

        // Filling VZERO event plane histograms
        psi2VZERO = CalculateEventPlaneVZERO(2);
        psi3VZERO = CalculateEventPlaneVZERO(3);
        histname = TString::Format("%s/histPsi2VZERO_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi2VZERO);
        histname = TString::Format("%s/histPsi3VZERO_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi3VZERO);
        psi2VZEROA = CalculateEventPlaneVZEROA(2);
        psi3VZEROA = CalculateEventPlaneVZEROA(3);
        histname = TString::Format("%s/histPsi2VZEROA_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi2VZEROA);
        histname = TString::Format("%s/histPsi3VZEROA_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi3VZEROA);
        psi2VZEROC = CalculateEventPlaneVZEROC(2);
        psi3VZEROC = CalculateEventPlaneVZEROC(3);
        histname = TString::Format("%s/histPsi2VZEROC_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi2VZEROC);
        histname = TString::Format("%s/histPsi3VZEROC_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi3VZEROC);

        // Calculating event plane using TPC and filling TPC event plane histograms
        //cout << "TPC Q Vector = (" << QxTPC << "," << QyTPC << ")" << endl;
        //cout << "TPC Positive Eta Q Vector = (" << QxTPCPositiveEta << "," << QyTPCPositiveEta << ")" << endl;
        //cout << "TPC Negative Eta Q Vector = (" << QxTPCNegativeEta << "," << QyTPCNegativeEta << ")" << endl;
        if (QxTPC != 0.0 && QyTPC != 0.0 && QxTPCPositiveEta != 0.0 && QyTPCPositiveEta != 0.0 && QxTPCNegativeEta != 0.0 && QyTPCNegativeEta != 0.0) {
            psi2TPC = TMath::ATan2(QyTPC, QxTPC)/2;
            if (psi2TPC < 0.0) psi2TPC += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi2TPC_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psi2TPC);
            //cout << "TPC event plane = " << psi2TPC << endl;
            psi2TPCPositiveEta = TMath::ATan2(QyTPCPositiveEta, QxTPCPositiveEta)/2;
            if (psi2TPCPositiveEta < 0.0) psi2TPCPositiveEta += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi2TPCPositiveEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psi2TPCPositiveEta);
            //cout << "TPC positive eta event plane = " << psi2TPCPositiveEta << endl;
            psi2TPCNegativeEta = TMath::ATan2(QyTPCNegativeEta, QxTPCNegativeEta)/2;
            if (psi2TPCNegativeEta < 0.0) psi2TPCNegativeEta += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi2TPCNegativeEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psi2TPCNegativeEta);
            //cout << "TPC negative eta event plane = " << psi2TPCNegativeEta << endl;

            psi3TPC = TMath::ATan2(QyTPC, QxTPC)/2;
            if (psi3TPC < 0.0) psi3TPC += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi3TPC_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psi3TPC);
            psi3TPCPositiveEta = TMath::ATan2(QyTPCPositiveEta, QxTPCPositiveEta)/2;
            if (psi3TPCPositiveEta < 0.0) psi3TPCPositiveEta += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi3TPCPositiveEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psi3TPCPositiveEta);
            psi3TPCNegativeEta = TMath::ATan2(QyTPCNegativeEta, QxTPCNegativeEta)/2;
            if (psi3TPCNegativeEta < 0.0) psi3TPCNegativeEta += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi3TPCNegativeEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psi3TPCNegativeEta);

            // Calculating sub-event plane resolutionfor different eta sections
            histname = TString::Format("%s/histSubEventPlaneResolution1_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, CalculateSubEventPlaneResolution(2, psi2VZERO, psi2TPCPositiveEta));
            histname = TString::Format("%s/histSubEventPlaneResolution2_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, CalculateSubEventPlaneResolution(2, psi2VZERO, psi2TPCNegativeEta));
            histname = TString::Format("%s/histSubEventPlaneResolution3_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, CalculateSubEventPlaneResolution(2, psi2TPCPositiveEta, psi2TPCNegativeEta));
        }
    }
    //histname = "fHistSumNTracks";
    //fHistManager.FillTH1(histname, sumAcceptedTracks);
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskJetVn::ExecOnce()
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
Bool_t AliAnalysisTaskJetVn::Run()
{
  if(fEventCuts.AcceptEvent(InputEvent()))
    return kTRUE;
  else
    return kFALSE;

  Float_t corrPsi2VZERO = 0.0;
  const AliQnCorrectionsQnVector *corrQnVectorVZERO;
  corrQnVectorVZERO = fFlowQnVectorMgr->GetDetectorQnVector("VZEROQoverM");
  if (corrQnVectorVZERO != NULL) {
    corrPsi2VZERO = (Float_t)corrQnVectorVZERO->EventPlane(2);
  }

  TF1* rhoFit = new TF1("rhoFit", "[0] * (1 + 2 * ([1] * cos(2 * x)))", 0.0, TMath::TwoPi());

  TH1* histBackgroundTracks = new TH1F("histBackgroundTracks", "histBackgroundTracks", 100, 0.0, TMath::TwoPi());

  TString histname;
  TString groupname;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        Float_t trackDeltaPhi = track->Phi() - corrPsi2VZERO;
        if (trackDeltaPhi < 0.0) trackDeltaPhi += TMath::TwoPi();
        histBackgroundTracks->Fill(trackDeltaPhi);
      }
    }
  }

  histBackgroundTracks->Fit("rhoFit");
  Float_t param0, param1;
  param0 = rhoFit->GetParameter(0);
  param1 = rhoFit->GetParameter(1);
  delete rhoFit;

  //DoTrackLoop();
  //DoJetLoop(param0, param1);
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskJetVn::Terminate(Option_t *)
{
}

// Function to calculate psi2 and psi3 from vzero a and c data
Double_t AliAnalysisTaskJetVn::CalculateEventPlaneVZERO(Int_t n)
{
    Double_t psi = -1.0;
    //Double_t Q2[] = {-999., -999.};
    //Double_t Q3[] = {-999., -999.};
    Double_t Qxan = 0, Qyan = 0;
    Double_t Qxcn = 0, Qycn = 0;
    Double_t Qxa3 = 0, Qya3 = 0;
    Double_t Qxc3 = 0, Qyc3 = 0;
    Double_t sumMa = 0, sumMc = 0;

    //AliVVZERO* aodV0 =(InputEvent())->GetVZEROData();

    for (Int_t iV0 = 0; iV0 < 64; iV0++) {                                      // iterate over 64 vzero segments
        Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);                      // azimuthal angle of VZERO segment
        //Float_t multV0 = aodV0->GetMultiplicity(iV0);                           // multiplicity of VZERO segment
        Float_t multV0 = InputEvent()->GetVZEROEqMultiplicity(iV0);             // gain adjusted multiplicity of VZERO segment
        if (iV0 < 32){                                                          // c side
            Qxcn += TMath::Cos(2.*phiV0) * multV0;                              // x axis, c side, v2
            Qycn += TMath::Sin(2.*phiV0) * multV0;                              // y axis, c side, v2
            Qxc3 += TMath::Cos(3.*phiV0) * multV0;                              // x axis, c side, v3
            Qyc3 += TMath::Sin(3.*phiV0) * multV0;                              // y axis, c side. v3
            sumMc += multV0;                                                    // running total multiplicity c side
        }
        else {                                                                  // a side
            Qxan += TMath::Cos(2.*phiV0) * multV0;                              // x axis, a side, v2
            Qyan += TMath::Sin(2.*phiV0) * multV0;                              // y axis, a side, v2
            Qxa3 += TMath::Cos(3.*phiV0) * multV0;                              // x axis, a side, v3
            Qya3 += TMath::Sin(3.*phiV0) * multV0;                              // y axis, a side, v3
            sumMa += multV0;                                                    // running total multiplicity a side
        }
    }
    if (sumMa <=0 || sumMc <= 0) return psi;
    if (n==2) {
        psi = .5*TMath::ATan2(Qyan+Qycn,Qxan+Qxcn);
        if (psi < 0.) psi += TMath::TwoPi()/2;
    }
    if (n==3) {
        psi = (1./3.)*TMath::ATan2(Qyan+Qycn,Qxan+Qxcn);
        if (psi < 0.) psi += TMath::TwoPi()/3;
    }

    return psi;
}

Double_t AliAnalysisTaskJetVn::CalculateEventPlaneVZEROA(Int_t n)
{
    Double_t psi = -999.0, phiV0 = -999.0;
    Double_t Qxan = 0.0, Qyan = 0.0;
    Double_t sumMc = 0.0;
    Float_t multV0 = 0.0;
    //AliVVZERO* aodV0 = (InputEvent())->GetVZEROData();
    for (Int_t iV0 = 0; iV0 < 32; iV0++) {                                      // iterate over 32 VZEROA segments
        phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);                               // azimuthal angle of VZERO segment
        //multV0 = aodV0->GetMultiplicity(iV0);                                   // multiplicity of vzero segment
        multV0 = InputEvent()->GetVZEROEqMultiplicity(iV0);                     // gain adjusted multiplicity of VZERO segment
        Qxan += TMath::Cos(n*phiV0) * multV0;                                   // x axis, c side, vn
        Qyan += TMath::Sin(n*phiV0) * multV0;                                   // y axis, c side, vn                                             // running total multiplicity c side
        sumMc += multV0;
    }
    if (sumMc <= 0) return psi;
    else {
        psi = TMath::ATan2(Qyan, Qxan)/n;
        if (psi < 0.) psi += TMath::TwoPi()/n;
    }
    return psi;
}

Double_t AliAnalysisTaskJetVn::CalculateEventPlaneVZEROC(Int_t n)
{
    Double_t psi = -999.0, phiV0 = -999.0;
    Double_t Qxcn = 0.0, Qycn = 0.0;
    Double_t sumMc = 0.0;
    Float_t multV0 = 0.0;
    for (Int_t iV0 = 32; iV0 < 64; iV0++) {                                     // iterate over 32 VZEROC segments
        phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);                               // azimuthal angle of VZERO segment
        multV0 = InputEvent()->GetVZEROEqMultiplicity(iV0);                     // gain adjusted multiplicity of VZERO segment
        Qxcn += TMath::Cos(n*phiV0) * multV0;                                   // x axis, c side, vn
        Qycn += TMath::Sin(n*phiV0) * multV0;                                   // y axis, c side, vn
        sumMc += multV0;                                                        // running total multiplicity c side
    }
    if (sumMc <= 0) return psi;
    else {
        psi = TMath::ATan2(Qycn, Qxcn)/n;
        if (psi < 0.) psi += TMath::TwoPi()/n;
    }
    return psi;
}

Double_t AliAnalysisTaskJetVn::CalculateSubEventPlaneResolution(Int_t n, Double_t psiI, Double_t psiJ)
{
    return TMath::Abs(TMath::Cos(n*(psiI - psiJ)));
}

Double_t AliAnalysisTaskJetVn::CalculateLocalRho(Double_t jetAngle, Double_t param0, Double_t param1)
{
    TF1* rhoFit = new TF1("rhoFit", "[0] * (1 + 2 * ([1] * cos(2 * x)))", 0.0, TMath::TwoPi());
    rhoFit->SetParameter(0, param0);
    rhoFit->SetParameter(1, param1);
    Double_t localRho;
    localRho = param0 * (1 + 2 * (param1 * cos(2 * jetAngle)));
    delete rhoFit;
    return localRho;
}

Double_t AliAnalysisTaskJetVn::CalculateTotalRho(Double_t param0, Double_t param1)
{
    TF1* rhoFit = new TF1("rhoFit", "[0] * (1 + 2 * ([1] * cos(2 * x)))", 0.0, TMath::TwoPi());
    rhoFit->SetParameter(0, param0);
    rhoFit->SetParameter(1, param1);

    Double_t totalRho;
    totalRho = rhoFit->Integral(0.0, TMath::TwoPi());

    delete rhoFit;

    return totalRho;
}


/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskJetVn * AliAnalysisTaskJetVn::AddTaskEmcalJetSample(
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
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
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

  TString name("AliAnalysisTaskJetVn");
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

  AliAnalysisTaskJetVn* sampleTask = new AliAnalysisTaskJetVn(name);
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
