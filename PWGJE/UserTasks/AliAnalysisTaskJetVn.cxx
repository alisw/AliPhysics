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

//#include "AliAnalysisTaskEmcalJetSample.h"
#include "AliAnalysisTaskJetVn.h"

// aliroot includes - these are stored in /home/william/alice/sw/ubuntu1604_x86-64/AliRoot/0-1/include
#include <AliCentrality.h>
#include <AliVVertex.h>
#include <AliVTrack.h>
#include <AliVVZERO.h>
#include <AliESDEvent.h>
#include <AliAODTrack.h>
#include <AliOADBContainer.h>
#include <AliMultSelection.h>
#include <AliInputEventHandler.h>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetVn);
/// \endcond

//Default constructor. Needed by ROOT I/O
AliAnalysisTaskJetVn::AliAnalysisTaskJetVn() :
    AliAnalysisTaskEmcalJet(),
    fHistManager(0x0)
{
}

// Standard constructor. Should be used by the user.
// @param[in] name Name of the task
AliAnalysisTaskJetVn::AliAnalysisTaskJetVn(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fHistManager(name)
{
    SetMakeGeneralHistograms(kTRUE);
}

// Destructor
AliAnalysisTaskJetVn::~AliAnalysisTaskJetVn()
{
}

// Performing run-independent initialization.
// Here the histograms should be instantiated.
void AliAnalysisTaskJetVn::UserCreateOutputObjects()
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

// This function allocates the histograms for basic tracking QA.
// A set of histograms (pT, eta, phi, difference between kinematic properties
// at the vertex and at the EMCal surface, number of tracks) is allocated
// per each particle container and per each centrality bin.
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
        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

            histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

            /*if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
                histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
                fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

                histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
                fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 200, -2, 2);

                histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
                fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, -fMaxBinPt/2, fMaxBinPt/2);

                histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
                fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 4);
            }*/

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
            //histtitle = TString::Format("%s/histPsi2VZERO_%d", groupname.Data(), cent);
            histtitle = "Psi2 from VZERO";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histPsi3VZERO_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histPsi3VZERO_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            histname = TString::Format("%s/histPsi2VZEROA_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histPsi2VZEROA_%d", groupname.Data(), cent);
            histtitle = "Psi2 from VZEROA";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histPsi3VZEROA_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histPsi3VZEROA_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            histname = TString::Format("%s/histPsi2VZEROC_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histPsi2VZEROC_%d", groupname.Data(), cent);
            histtitle = "Psi2 from VZEROC";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histPsi3VZEROC_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histPsi3VZEROC_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            histname = TString::Format("%s/histPsi2TPC_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histPsi2TPC_%d", groupname.Data(), cent);
            histtitle = "Psi2 from TPC";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histPsi3TPC_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histPsi3TPC_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            histname = TString::Format("%s/histPsi2TPCPositiveEta_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histPsi2TPCPositiveEta_%d", groupname.Data(), cent);
            histtitle = "Psi2 from TPC positive eta side";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histPsi3TPCPositiveEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histPsi3TPCPositiveEta_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            histname = TString::Format("%s/histPsi2TPCNegativeEta_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histPsi2TPCNegativeEta_%d", groupname.Data(), cent);
            histtitle = "Psi2 from TPC negative eta";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histPsi3TPCNegativeEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histPsi3TPCNegativeEta_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            // histograms for particle angle relative to the event plane
            histname = TString::Format("%s/histTrackPhiMinusPsi2_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histTrackPhiMinusPsi2_%d", groupname.Data(), cent);
            histtitle = "Track phi minus Psi2";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histTrackPhiMinusPsi3_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s/histTrackPhiMinusPsi3_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            // sub-event plane resolution histograms for calculating event plane resolution
            histname = TString::Format("%s/histSubEventPlaneResolution1_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histSubEventPlaneResolution1_%d", groupname.Data(), cent);
            histtitle = "Sub event plane resolution 1";
            fHistManager.CreateTH1(histname, histtitle, 100, -1, 1);

            histname = TString::Format("%s/histSubEventPlaneResolution2_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histSubEventPlaneResolution2_%d", groupname.Data(), cent);
            histtitle = "Sub event plane resolution 2";
            fHistManager.CreateTH1(histname, histtitle, 100, -1, 1);

            histname = TString::Format("%s/histSubEventPlaneResolution3_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histSubEventPlaneResolution3_%d", groupname.Data(), cent);
            histtitle = "Sub event plane resolution 3";
            fHistManager.CreateTH1(histname, histtitle, 100, -1, 1);
        }
    }

    /*histname = "fHistSumNTracks";
    histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
    if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
    }
    else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
    }*/
}

// This function allocates the histograms for basic jet QA.
// A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
// per each jet container and per each centrality bin.
void AliAnalysisTaskJetVn::AllocateJetHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        cout << "Groupname = " << groupname << "\n";
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

            //histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
            //fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

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

            /*if (!jetCont->GetRhoName().IsNull()) {
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
                fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
            }*/

            // histograms for jet angle relative to the event plane
            histname = TString::Format("%s/histJetPhiMinusPsi2_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histJetPhiMinusPsi2_%d", groupname.Data(), cent);
            histtitle = "Jet phi minus psi2";
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            /*histname = TString::Format("%s/histJetPhiMinusPsi3_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s#i/histJetPhiMinusPsi3_%d", groupname.Data(), cent);
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());*/

            // histograms for in-plane vs out-of-plane jets
            histname = TString::Format("%s/histJetsInPlaneOutOfPlane_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s/histJetsInPlaneOutOfPlane_%d", groupname.Data(), cent);
            histtitle = "Jet phi minus psi2";
            fHistManager.CreateTH2(histname, histtitle, 2, 0, 2, 10, 0, 100);
        }
    }
}

// The body of this function should contain instructions to fill the output histograms.
// This function is called inside the event loop, after the function Run() has been
// executed successfully (i.e. it returned kTRUE).
// @return Always kTRUE
Bool_t AliAnalysisTaskJetVn::FillHistograms()
{
    DoJetLoop();
    DoTrackLoop();

    return kTRUE;
}

// This function performs a loop over the reconstructed jets
// in the current event and fills the relevant histograms.
void AliAnalysisTaskJetVn::DoJetLoop()
{
    //Int_t nJetIn[10][fNcentBins], nJetOut[10][fNcentBins];
    Double_t psi2, psi3, phiMinusPsi2, phiMinusPsi3; // angles for jet vn calculation

    TString histname;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        UInt_t count = 0;

        psi2 = CalculateEventPlaneVZERO(2);
        //psi3 = CalculateEventPlaneVZERO(3);

        for(auto jet : jetCont->accepted()) {
            if (!jet) continue;
            count++;

            histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Pt());

            //histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, jet->Area());

            histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Phi());

            histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Eta());

            // Filling histos for angle relative to event plane
            histname = TString::Format("%s/histJetPhiMinusPsi2_%d", groupname.Data(), fCentBin);
            phiMinusPsi2 = jet->Phi() - psi2;
            if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
            fHistManager.FillTH1(histname, phiMinusPsi2);
            /*histname = TString::Format("%s/histJetPhiMinusPsi3_%d", groupname.Data(), fCentBin);
            phiMinusPsi3 = jet->Phi() - psi3;
            if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
            fHistManager.FillTH1(histname, phiMinusPsi3);*/

            /*if (jetCont->GetRhoParameter()) {
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
            }*/

            if ((phiMinusPsi2 < TMath::Pi()/4) || (phiMinusPsi2 >= 3*TMath::Pi()/4 && phiMinusPsi2 < 5*TMath::Pi()/4) || (phiMinusPsi2 >= 7*TMath::Pi()/4)) {
                histname = TString::Format("%s/histJetsInPlaneOutOfPlane_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 0.5, jet->Pt());
            }
            if ((phiMinusPsi2 >= TMath::Pi()/4 && phiMinusPsi2 < 3*TMath::Pi()/4) || (phiMinusPsi2 >= 5*TMath::Pi()/4 && phiMinusPsi2 < 7*TMath::Pi()/4)) {
                histname = TString::Format("%s/histJetsInPlaneOutOfPlane_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, 1.5, jet->Pt());
            }
        }
        histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

// This function performs a loop over the reconstructed tracks
// in the current event and fills the relevant histograms.
void AliAnalysisTaskJetVn::DoTrackLoop()
{
    AliClusterContainer* clusCont = GetClusterContainer(0);

    Double_t minEta = -0.7, maxEta = 0.7;
    Double_t psi2, psi3, phiMinusPsi2, phiMinusPsi3; // angles for jet vn calculation
    Double_t psi2VZERO, psi3VZERO, psi2VZEROA, psi3VZEROA, psi2VZEROC, psi3VZEROC;
    Double_t psiTPC = -999.0, psiTPCPositiveEta = -999.0, psiTPCNegativeEta = -999.0;
    Double_t QxTPC = 0, QyTPC = 0, QxTPCPositiveEta = 0, QyTPCPositiveEta = 0, QxTPCNegativeEta = 0, QyTPCNegativeEta = 0;

    TString histname;
    TString groupname;
    UInt_t sumAcceptedTracks = 0;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = partCont->GetName();

        UInt_t count = 0;

        psi2 = CalculateEventPlaneVZERO(2);
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

                /*histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

                histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

                histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());*/

                // Filling histos for angle relative to event plane
                phiMinusPsi2 = track->Phi() - psi2;
                //phiMinusPsi3 = track->Phi() - psi3;
                if (phiMinusPsi2 < 0.0) phiMinusPsi2 += TMath::TwoPi();
                //if (phiMinusPsi3 < 0.0) phiMinusPsi3 += TMath::TwoPi();
                histname = TString::Format("%s/histTrackPhiMinusPsi2_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, phiMinusPsi2);
                /*histname = TString::Format("%s/histTrackPhiMinusPsi3_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, phiMinusPsi3);*/

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

                if (clusCont) {
                    Int_t iCluster = track->GetEMCALcluster();
                    if (iCluster >= 0) {
                        AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
                        /*if (cluster) {
                            histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
                            fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
                        }*/
                    }
                }
            }
        }
        sumAcceptedTracks += count;

        histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);

        // Filling VZERO event plane histograms
        psi2VZERO = CalculateEventPlaneVZERO(2);
        //psi3VZERO = CalculateEventPlaneVZERO(3);
        histname = TString::Format("%s/histPsi2VZERO_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi2VZERO);
        //histname = TString::Format("%s/histPsi3VZERO_%d", groupname.Data(), fCentBin);
        //fHistManager.FillTH1(histname, psi3VZERO);
        psi2VZEROA = CalculateEventPlaneVZEROA(2);
        //psi3VZEROA = CalculateEventPlaneVZEROA(3);
        histname = TString::Format("%s/histPsi2VZEROA_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi2VZEROA);
        //histname = TString::Format("%s/histPsi3VZEROA_%d", groupname.Data(), fCentBin);
        //fHistManager.FillTH1(histname, psi3VZEROA);
        psi2VZEROC = CalculateEventPlaneVZEROC(2);
        //psi3VZEROC = CalculateEventPlaneVZEROC(3);
        histname = TString::Format("%s/histPsi2VZEROC_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, psi2VZEROC);
        //histname = TString::Format("%s/histPsi3VZEROC_%d", groupname.Data(), fCentBin);
        //fHistManager.FillTH1(histname, psi3VZEROC);

        // Calculating event plane using TPC and filling TPC event plane histograms
        //cout << "TPC Q Vector = (" << QxTPC << "," << QyTPC << ")" << endl;
        //cout << "TPC Positive Eta Q Vector = (" << QxTPCPositiveEta << "," << QyTPCPositiveEta << ")" << endl;
        //cout << "TPC Negative Eta Q Vector = (" << QxTPCNegativeEta << "," << QyTPCNegativeEta << ")" << endl;
        if (QxTPC != 0.0 && QyTPC != 0.0 && QxTPCPositiveEta != 0.0 && QyTPCPositiveEta != 0.0 && QxTPCNegativeEta != 0.0 && QyTPCNegativeEta != 0.0) {
            psiTPC = TMath::ATan2(QyTPC, QxTPC)/2;
            if (psiTPC < 0.0) psiTPC += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi2TPC_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psiTPC);
            //cout << "TPC event plane = " << psiTPC << endl;
            psiTPCPositiveEta = TMath::ATan2(QyTPCPositiveEta, QxTPCPositiveEta)/2;
            if (psiTPCPositiveEta < 0.0) psiTPCPositiveEta += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi2TPCPositiveEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psiTPCPositiveEta);
            //cout << "TPC positive eta event plane = " << psiTPCPositiveEta << endl;
            psiTPCNegativeEta = TMath::ATan2(QyTPCNegativeEta, QxTPCNegativeEta)/2;
            if (psiTPCNegativeEta < 0.0) psiTPCNegativeEta += TMath::TwoPi()/2;
            histname = TString::Format("%s/histPsi2TPCNegativeEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, psiTPCNegativeEta);
            //cout << "TPC negative eta event plane = " << psiTPCNegativeEta << endl;

            // Calculating sub-event plane resolutionfor different eta sections
            histname = TString::Format("%s/histSubEventPlaneResolution1_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, CalculateSubEventPlaneResolution(2, psi2VZERO, psiTPCPositiveEta));
            histname = TString::Format("%s/histSubEventPlaneResolution2_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, CalculateSubEventPlaneResolution(2, psi2VZERO, psiTPCNegativeEta));
            histname = TString::Format("%s/histSubEventPlaneResolution3_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, CalculateSubEventPlaneResolution(2, psiTPCPositiveEta, psiTPCNegativeEta));
        }
    }

    //histname = "fHistSumNTracks";
    //fHistManager.FillTH1(histname, sumAcceptedTracks);
}

// This function is executed automatically for the first event.
// Some extra initialization can be performed here.
void AliAnalysisTaskJetVn::ExecOnce()
{
    AliAnalysisTaskEmcalJet::ExecOnce();
}

// Run analysis code here, if needed.
// It will be executed before FillHistograms().
// If this function return kFALSE, FillHistograms() will *not*
// be executed for the current event
// @return Always kTRUE
Bool_t AliAnalysisTaskJetVn::Run()
{
    return kTRUE;
}

// This function is called once at the end of the analysis.
void AliAnalysisTaskJetVn::Terminate(Option_t *) {}

// Function to calculate psi2 and psi3 from vzero a and c data
Double_t AliAnalysisTaskJetVn::CalculateEventPlaneVZERO(Int_t n)
{
    Double_t psi = -1.0;
    Double_t Q2[] = {-999., -999.};
    Double_t Q3[] = {-999., -999.};
    Double_t Qxan = 0, Qyan = 0;
    Double_t Qxcn = 0, Qycn = 0;
    Double_t Qxa3 = 0, Qya3 = 0;
    Double_t Qxc3 = 0, Qyc3 = 0;
    Double_t sumMa = 0, sumMc = 0;

    AliVVZERO* aodV0 =(InputEvent())->GetVZEROData();

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
    /*Double_t iCentSPD = GetCentrality("CL1");
    Double_t QyanCor = (Qyan - fMQ[1][0][0]->GetBinContent(iCentSPD+1))/fWQ[1][0][0]->GetBinContent(iCentSPD+1); // used for flattening q vector histo
    Double_t QycnCor = (Qycn - fMQ[1][1][0]->GetBinContent(iCentSPD+1))/fWQ[1][1][0]->GetBinContent(iCentSPD+1);
    Double_t QxanCor = (Qxan - fMQ[0][0][0]->GetBinContent(iCentSPD+1))/fWQ[0][0][0]->GetBinContent(iCentSPD+1);
    Double_t QxcnCor = (Qxcn - fMQ[0][1][0]->GetBinContent(iCentSPD+1))/fWQ[0][1][0]->GetBinContent(iCentSPD+1);*/
    //comb[0] = .5*TMath::ATan2(QyanCor+QycnCor,QxanCor+QxcnCor);                 // combined psi2 angle
    if (n==2) {
        psi = .5*TMath::ATan2(Qyan+Qycn,Qxan+Qxcn);
        if (psi < 0.) psi += TMath::TwoPi()/2;
    }
    /*QyanCor = (Qya3 - fMQ[1][0][1]->GetBinContent(iCentSPD+1))/fWQ[1][0][1]->GetBinContent(iCentSPD+1);
    QycnCor = (Qyc3 - fMQ[1][1][1]->GetBinContent(iCentSPD+1))/fWQ[1][1][1]->GetBinContent(iCentSPD+1);
    QxanCor = (Qxa3 - fMQ[0][0][1]->GetBinContent(iCentSPD+1))/fWQ[0][0][1]->GetBinContent(iCentSPD+1);
    QxcnCor = (Qxc3 - fMQ[0][1][1]->GetBinContent(iCentSPD+1))/fWQ[0][1][1]->GetBinContent(iCentSPD+1);*/
    //comb[1] = (1./3.)*TMath::ATan2(QyanCor+QycnCor,QxanCor+QxcnCor);            // combined psi3 angle
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
    AliVVZERO* aodV0 = (InputEvent())->GetVZEROData();
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
    AliVVZERO* aodV0 = (InputEvent())->GetVZEROData();
    for (Int_t iV0 = 32; iV0 < 64; iV0++) {                                     // iterate over 32 VZEROC segments
        phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);                               // azimuthal angle of VZERO segment
        //multV0 = aodV0->GetMultiplicity(iV0);                                   // multiplicity of vzero segment
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

/*Double_t AliAnalysisTaskJetVn::CalculateEventPlaneTPC(Int_t n)
{
    Double_t psi = -1.0;
    Double_t Qx = 0, Qy = 0;
    Double_t sumMc = 0;
    Double_t phiTrack, etaTrack;
    Double_t minEta = -0.7, maxEta = 0.7;
    Int_t nTracks = 0;
    AliVParticle* partTPC;

    //AliVVZERO* aodV0 =(InputEvent())->GetVZEROData();
    nTracks = InputEvent()->GetNumberOfTracks();
    for(Int_t i = 0; i <= nTracks; i++) {
        partTPC = InputEvent()->GetTrack(i);
        //phiTrack = (InputEvent())->GetTrack(i)->Phi();
        phiTrack = partTPC->Phi();
        //etaTrack = (InputEvent())->GetTrack(i)->Eta();
        etaTrack = partTPC->Eta();
        if (etaTrack >= minEta && etaTrack <= maxEta) {
            Qx += TMath::Cos(2.0*phiTrack);
            Qy += TMath::Sin(2.0*phiTrack);
        }
    }

    //if (sumMc <= 0) return psi;
    if (n==2) {
        psi = .5*TMath::ATan2(Qy, Qx);
        if (psi < 0.) psi += TMath::TwoPi()/2;
    }
    if (n==3) {
        psi = (1./3.)*TMath::ATan2(Qy, Qx);
        if (psi < 0.) psi += TMath::TwoPi()/3;
    }

    return psi;
}*/

/*Double_t AliAnalysisTaskJetVn::CalculateEventPlaneTPCPositiveEta(Int_t n)
{
    Double_t psi = -1.0;
    return psi;
}*/

/*Double_t AliAnalysisTaskJetVn::CalculateEventPlaneTPCNegativeEta(Int_t n)
{
    Double_t psi = -1.0;
    return psi;
}*/

/*Double_t AliAnalysisTaskJetVn::CalculateEventPlaneResolution(Double_t psiVZERO, Double_t psiTPCPositiveEta, Double_t psiTPCNegativeEta)
{
    Double_t R = TMath::Sqrt(TMath::Cos(2*(psiVZERO - psiTPCPositiveEta)) * TMath::Cos(2*(psiVZERO - psiTPCNegativeEta)) / TMath::Cos(2*(psiTPCPositiveEta - psiTPCNegativeEta)));
    return R;
}*/

Double_t AliAnalysisTaskJetVn::CalculateSubEventPlaneResolution(Int_t n, Double_t psiI, Double_t psiJ)
{
  return TMath::Abs(TMath::Cos(n*(psiI - psiJ)));
}

// This function adds the task to the analysis manager. Often, this function is called
// by an AddTask C macro. However, by compiling the code, it ensures that we do not
// have to deal with difficulties caused by CINT.
AliAnalysisTaskJetVn * AliAnalysisTaskJetVn::AddTaskJetVn(
    const char* nTracks,
    const char* nClusters,
    const char* nCells,
    const char* suffix)
{
    //==========================================================================
    // Get the pointer to the existing analysis manager via the static access method.
    //==========================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
        return 0;
    }

    //==========================================================================
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==========================================================================
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

    //==========================================================================
    // Init the task and do settings
    //==========================================================================
    TString trackName(nTracks);
    TString clusName(nClusters);
    TString cellName(nCells);

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

    //AliAnalysisTaskJetVn* sampleTask = new AliAnalysisTaskEmcalJetSample(name);
    AliAnalysisTaskJetVn* jetVnTask = new AliAnalysisTaskJetVn(name);
    jetVnTask->SetCaloCellsName(cellName);
    jetVnTask->SetVzRange(-10,10);

    if (trackName == "mcparticles") {
        jetVnTask->AddMCParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
        jetVnTask->AddTrackContainer(trackName);
    }
    else if (!trackName.IsNull()) {
        jetVnTask->AddParticleContainer(trackName);
    }
    jetVnTask->AddClusterContainer(clusName);

    //==========================================================================
    // Final settings, pass to manager and set the containers
    //==========================================================================
    mgr->AddTask(jetVnTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
        TList::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput  (jetVnTask, 0,  cinput1 );
    mgr->ConnectOutput (jetVnTask, 1, coutput1 );

    return jetVnTask;
}
