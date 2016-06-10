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

// general ROOT includes
#include <TCanvas.h>
#include <TChain.h>
#include <TMath.h>
#include <TProfile.h>
#include <TAxis.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TObjArray.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliAnalysisTaskEmcalJetHF.h"

using std::cout;
using std::endl;


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetHF);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetHF::AliAnalysisTaskEmcalJetHF() :
AliAnalysisTaskEmcalJet(),
fVevent(),
fESD(),
fAOD(),
fpidResponse(),
fEventCounter(),
fInvmassCut(),
fHistManager(),
fdEdx(),  //Histo's
fM20(),
fM02(),
fM20EovP(),
fM02EovP(),
fInvmassLS(),
fInvmassULS(),
fEMCTrketa(),
fEMCTrkphi(),
fHistJetEovP(),
fHistJetEovPvPt(),
fHistClusEovP(),
fHistClusEovPnonlin(),
fHistClusEovPHadCorr()
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetHF::AliAnalysisTaskEmcalJetHF(const char *name) :
AliAnalysisTaskEmcalJet(name, kTRUE),
fVevent(0),
fESD(0),
fAOD(0),
fpidResponse(0),
fEventCounter(0),
fInvmassCut(0),
fHistManager(name),
fdEdx(0),  //Histo's
fM20(0),
fM02(0),
fM20EovP(0),
fM02EovP(0),
fInvmassLS(0),
fInvmassULS(0),
fEMCTrketa(0),
fEMCTrkphi(0),
fHistJetEovP(0),
fHistJetEovPvPt(0),
fHistClusEovP(0),
fHistClusEovPnonlin(0),
fHistClusEovPHadCorr(0)
{
    SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetHF::~AliAnalysisTaskEmcalJetHF()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetHF::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    /////////////////////////////////////////////////
    //Automatic determination of the analysis mode//
    ////////////////////////////////////////////////
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");


    AllocateClusterHistograms();
    AllocateTrackHistograms();
    AllocateJetHistograms();
    AllocateCellHistograms();
    
    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
        fOutput->Add(obj);
    }
    fOutput->Add(fdEdx);
    fOutput->Add(fM20);
    fOutput->Add(fM02);
    fOutput->Add(fM20EovP);
    fOutput->Add(fM02EovP);
    fOutput->Add(fInvmassLS);
    fOutput->Add(fInvmassULS);
    fOutput->Add(fEMCTrketa);
    fOutput->Add(fEMCTrkphi);
    fOutput->Add(fHistJetEovP);
    fOutput->Add(fHistJetEovPvPt);
    fOutput->Add(fHistClusEovP);
    fOutput->Add(fHistClusEovPnonlin);
    fOutput->Add(fHistClusEovPHadCorr);
    
    
    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic EMCal cluster QA.
 * A set of histograms (energy, eta, phi, number of cluster) is allocated
 * per each cluster container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetHF::AllocateClusterHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliClusterContainer* clusCont = 0;
    TIter next(&fClusterCollArray);
    while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
        groupname = clusCont->GetName();
        fHistManager.CreateHistoGroup(groupname);

        fHistClusEovP = new TH1F("ClusEovP","Cluster EovP",130,0,1.3);
        fHistClusEovPnonlin = new TH1F("ClusEovPnonlin","nonlin Cluster EovP",130,0,1.3);
        fHistClusEovPHadCorr = new TH1F("ClusEovPHaddCorr","HadCorr Cluster EovP",130,0,1.3);
        fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
        fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);

        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

            histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster}^{exotic} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

            histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster}^{non-lin.corr.} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

            histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster}^{had.corr.} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

            histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{custer};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

            histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{custer};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

            histname = TString::Format("%s/histNClusters_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;number of clusters;events", histname.Data());
            if (fForceBeamType != kpp) {
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
            }
        }
    }
}

/*
 * This function allocates the histograms for basic EMCal QA.
 * One 2D histogram with the cell energy spectra and the number of cells
 * per event is allocated per each centrality bin.
 */
void AliAnalysisTaskEmcalJetHF::AllocateCellHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname(fCaloCellsName);
    
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
        histname = TString::Format("%s/histCellEnergyvsAbsId_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;cell abs. ID;#it{E}_{cell} (GeV);counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 20000, 0, 20000, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

        histname = TString::Format("%s/histNCells_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;number of cells;events", histname.Data());
        if (fForceBeamType != kpp) {
            fHistManager.CreateTH1(histname, histtitle, 500, 0, 6000);
        }
        else {
            fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
        }
    }
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetHF::AllocateTrackHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = partCont->GetName();
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

            if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
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
            }
            
            histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;number of tracks;events", histname.Data());
            if (fForceBeamType != kpp) {
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
            }
        }
    }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetHF::AllocateJetHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        fHistManager.CreateHistoGroup(groupname);

        fHistJetEovP    = new TH1F("JetEovP","HFE Jet EovP",200,0,5);
        fHistJetEovPvPt = new TH2F("JetEovPvPt","HFE Jet EovP vs Pt",100,0,1,100,0,100);

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

            if (!jetCont->GetRhoName().IsNull()) {
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
                fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
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
Bool_t AliAnalysisTaskEmcalJetHF::FillHistograms()
{
    DoJetLoop();
    DoTrackLoop();
    DoClusterLoop();
    DoCellLoop();
    
    return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHF::DoJetLoop()
{
    TString histname;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);

    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        UInt_t count = 0;
        Double_t rhoVal = 0;
        rhoVal = jetCont->GetRhoVal();

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
            
                if (jetCont->GetRhoParameter()) {
                    histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
                }
            histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, count);



            Float_t ptLeading = jetCont->GetLeadingHadronPt(jet);
            Float_t corrPt = jet->Pt() - rhoVal * jet->Area();
            TLorentzVector leadPart;
            jetCont->GetLeadingHadronMomentum(leadPart, jet);

            //cout<<"JET Corr Pt: "<<corrPt<<"  Pt of Leading Hadron: "<< ptLeading<<endl;

            AliParticleContainer* tracks = jetCont->GetParticleContainer();
            if (tracks) {
                for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
                    AliVParticle *track = jet->TrackAt(it, tracks->GetArray());

                    if (track) {
                        Double_t JetTrackP = track->P();
                        Double_t JetTrackPt = track->Pt();


                        Double_t dphi = TVector2::Phi_0_2pi(track->Phi() - jet->Phi());
                        Double_t deta = track->Eta() - jet->Eta();
                        Double_t dist = TMath::Sqrt(deta * deta + dphi * dphi);

                        //+++++----Track matching to EMCAL
                        AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(track);
                        Int_t EMCalIndex = -1;
                        EMCalIndex = Vtrack->GetEMCALcluster();
                        if(EMCalIndex < 0) continue;

                        AliVCluster *clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
                        if(!(clustMatch && clustMatch->IsEMCAL())) continue;

                        Double_t clustMatchE = clustMatch->E();
                        Double_t clustMatchNonLinE = clustMatch->GetNonLinCorrEnergy();
                        Double_t clustMatchHadCorr = clustMatch->GetHadCorrEnergy();

                        //cout<<"Other Track Pt: "<< JetTrackPt << " Track Matched to Cluster Energy: "<< clustMatchE<<endl;
                        fHistJetEovP->Fill(clustMatchNonLinE / Vtrack->P());


                    }// Accepted Jet Track
                }//Track loop
            }//Track Container
        } //jet loop
    }//Jet Container
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHF::DoTrackLoop()
{
    AliClusterContainer* clusCont = GetClusterContainer(0);

    TString histname;
    TString groupname;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = partCont->GetName();
        UInt_t count = 0;
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

                histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

                histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

                histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());

                Double_t TrackP = track->P();

                //////////////////////
                //Photonic Electrons//
                //////////////////////

                Bool_t fFlagPhotonicElec = kFALSE;

                //----+++Identify Non-HFE
                //SelectPhotonicElectron(iTracks,track,fFlagPhotonicElec);

                if (clusCont) {
                    Int_t iCluster = track->GetEMCALcluster();
                    if (iCluster >= 0) {
                        AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
                        if (cluster) {

                            Double_t ClusterE = cluster->E();
                            Double_t ClusterNonLinE = cluster->GetNonLinCorrEnergy();
                            Double_t ClusterHadCorr = cluster->GetHadCorrEnergy();

                            Double_t sigma = 0.0;


                            sigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
                            //cout<<"Sigma: " << sigma << endl;
                            fHistClusEovP->Fill(ClusterE / TrackP);
                            fHistClusEovPnonlin->Fill( ClusterNonLinE / TrackP);
                            fHistClusEovPHadCorr->Fill( ClusterHadCorr / TrackP);

                            histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
                            fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
                        }
                    }
                }
            }
        }

        histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHF::DoClusterLoop()
{
    TString histname;
    TString groupname;
    AliClusterContainer* clusCont = 0;
    TIter next(&fClusterCollArray);
    while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
        groupname = clusCont->GetName();

        for(auto cluster : clusCont->all()) {
            if (!cluster) continue;
            
            if (cluster->GetIsExotic()) {
                histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, cluster->E());
            }
        }

        UInt_t count = 0;
        for(auto cluster : clusCont->accepted()) {
            if (!cluster) continue;
            count++;
            
            AliTLorentzVector nPart;
            cluster->GetMomentum(nPart, fVertex);
            
            histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, cluster->E());
            
            histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());
            
            histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());
            
            histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, nPart.Phi_0_2pi());
            
            histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, nPart.Eta());
        }

        histname = TString::Format("%s/histNClusters_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

/**
 * This function performs a loop over the reconstructed EMCal cells
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHF::DoCellLoop()
{
    if (!fCaloCells) return;

    TString histname;

    const Short_t ncells = fCaloCells->GetNumberOfCells();

    histname = TString::Format("%s/histNCells_%d", fCaloCellsName.Data(), fCentBin);
    fHistManager.FillTH1(histname, ncells);

    histname = TString::Format("%s/histCellEnergyvsAbsId_%d", fCaloCellsName.Data(), fCentBin);
    for (Short_t pos = 0; pos < ncells; pos++) {
        Double_t amp   = fCaloCells->GetAmplitude(pos);
        Short_t absId  = fCaloCells->GetCellNumber(pos);

        fHistManager.FillTH2(histname, absId, amp);
    }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetHF::ExecOnce()
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
Bool_t AliAnalysisTaskEmcalJetHF::Run()
{
    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return kFALSE;
    }

    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (fESD) {
        //   printf("fESD available\n");
        //return;
    }

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (fAOD) {
        // printf("fAOD available\n");
        //return;
    }

    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();

    fpidResponse = fInputHandler->GetPIDResponse();
    if (!fpidResponse) {
        printf("ERROR: fpidResponse not available\n");
        return kFALSE;
    }

    Int_t ntracks;
    ntracks = fVevent->GetNumberOfTracks();
    //printf("There are %d tracks in this event\n",ntracks);

    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    Double_t NcontV = pVtx->GetNContributors();

   // if(NcontV<2)return kTRUE;  //Events with 2 Tracks

    //////////////////////
    //EMcal cluster info//
    //////////////////////
    EMCalClusterInfo();


    Int_t numberofvertices = 100;
    if(fAOD) numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    if(IsAODanalysis())
    {
        for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
            AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
            if(!aodvertex) continue;
            if(aodvertex->GetType()==AliAODVertex::kKink) {
                AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
                if(!mother) continue;
                Int_t idmother = mother->GetID();
                listofmotherkink[numberofmotherkink] = idmother;
                numberofmotherkink++;
            }
        }
    } //+++

    fEventCounter++;

    //cout<<"Event: "<<fEventCounter<<"  Number of Vertices: "<<numberofvertices<<"  Number of Tracks:  "<<ntracks<< endl;

    return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::EMCalClusterInfo()
{
    /////////////////////////////
    //EMCAL cluster information//
    ////////////////////////////

    Int_t Nclust = -999;
    TVector3 clustpos;
    Float_t  emcx[3]; // cluster pos
    Double_t clustE=-999, emcphi = -999, emceta=-999;
    Nclust = fVevent->GetNumberOfCaloClusters();
    for(Int_t icl=0; icl<Nclust; icl++)
    {
        AliVCluster *clust = 0x0;
        clust = fVevent->GetCaloCluster(icl);
        if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);

        if(clust && clust->IsEMCAL())
        {
            clustE = clust->E();
            clust->GetPosition(emcx);
            clustpos.SetXYZ(emcx[0],emcx[1],emcx[2]);
            emcphi = clustpos.Phi();
            emceta = clustpos.Eta();
            //fHistClustE->Fill(clustE);
            //fEMCClsEtaPhi->Fill(emceta,emcphi);
        }
    }
}
//________________________________________________________________________

//___________________________________________________________________________

void AliAnalysisTaskEmcalJetHF::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec)
{
    ///////////////////////////////
    //Photonic electron selection//
    ///////////////////////////////

    Bool_t flagPhotonicElec = kFALSE;
    Double_t ptAsso=-999., nsigma=-999.0;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;

    for(Int_t jTracks = 0; jTracks<fVevent->GetNumberOfTracks(); jTracks++){
        if(jTracks==itrack) continue;

        AliVParticle* VtrackAsso = fVevent->GetTrack(jTracks);
        if (!VtrackAsso) {
            printf("ERROR: Could not receive track %d\n", jTracks);
            continue;
        }

        AliVTrack *trackAsso = dynamic_cast<AliVTrack*>(VtrackAsso);
        if(!trackAsso) continue;

        if(IsAODanalysis()) {
            AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
            if(!atrackAsso) continue;
            if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(atrackAsso->GetTPCNcls() < 70) continue;
            if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
            if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
        }

        nsigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        ptAsso = trackAsso->Pt();
        Int_t chargeAsso = trackAsso->Charge();
        Int_t charge = track->Charge();

        if(ptAsso <0.2) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;

        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;

        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;

        AliKFParticle::SetField(fVevent->GetMagneticField());

        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);

        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

        Double_t mass=-999., width = -999;
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);

        if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
        if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);

        if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
            flagPhotonicElec = kTRUE;
        }
    }
    fFlagPhotonicElec = flagPhotonicElec;
}


/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetHF::Terminate(Option_t *)
{
    cout<<"#######################"<<endl;
    cout<<"#### Task Finished ####"<<endl;
    cout<<"#######################"<<endl;
}
