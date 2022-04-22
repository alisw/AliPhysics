/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskAlphaPiAOD class
// analysis task for the Hyper Hydrogen 4 Analysis
//
// Copied from AliAnalysisTaskAlphapiAOD class
//
// Author:
// Bong-Hwi Lim bong-hwi.lim@cern.ch
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskAlphaPiAOD.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "AliAnalysisDataContainer.h"
#include "AliLog.h"
using std::string;

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>
#include <TRandom3.h>
#include <TTree.h>

// ALIROOT includes
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODcascade.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"
#include "Math/Vector4D.h"
#include "THistManager.h"
#include "AliKFVertex.h"
#include "AliKFParticle.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskAlphaPiAOD);
///\endcond

namespace {

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

double Sq(double x) {
    return x * x;
}

constexpr int kHyperPdg{1010010040};
constexpr double kHyperMass{3.931};  /// from AliPDG.cxx

}  // namespace

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskAlphaPiAOD::AliAnalysisTaskAlphaPiAOD(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                                    fEventCut{false},
                                                                                    fHistos{nullptr},
                                                                                    fMC{isMC} {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskAlphaPiAOD::~AliAnalysisTaskAlphaPiAOD() {
    if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
        return;
    if (fHistos)
        delete fHistos;
    if (fTree)
        delete fTree;
    if (!fMC) {
        delete fRecHyper;
    }
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskAlphaPiAOD::UserCreateOutputObjects() {
    fHistos = new THistManager("hyperhists");
    fEventCut.AddQAplotsToList(fHistos->GetListOfHistograms());  // EventCuts QA Histograms

    // QA Histograms
    fHistos->CreateTH2("QA/hTPCPIDAlpha", "TPC PID #alpha;#it{p}_{T} (GeV/#it{c});n_{#sigma} #alpha", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDAntiAlpha", "TPC PID #bar{#alpha};#it{p}_{T} (GeV/#it{c});n_{#sigma} #bar{#alpha}", 200, 0, 20, 200, 0, 200);

    // MC QA Histograms
    // if (fMC) {
    //     fHistos->CreateTH2("QA_MC/hTPCPIDAlpha", "TPC PID #alpha;#it{p}_{T} (GeV/#it{c});n_{#sigma} #alpha", 200, 0, 20, 200, 0, 200);
    //     fHistos->CreateTH2("QA_MC/hTPCPIDAntiAlpha", "TPC PID #bar{#alpha};#it{p}_{T} (GeV/#it{c});n_{#sigma} #bar{#alpha}", 200, 0, 20, 200, 0, 200);
    // }

    fRecHyper = fMC ? &fGenHyper : new StructHyper;
    OpenFile(2);
    fTree = new TTree("HyperTree", "HyperTree");
    if (fMC) {
        fTree->Branch("HyperMC", &fGenHyper);
        fMCEvent = MCEvent();
    } else {
        fTree->Branch("StructHyper", fRecHyper);
    }

    PostAllData();
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskAlphaPiAOD::UserExec(Option_t *) {
    AliAODEvent *ev = (AliAODEvent *)InputEvent();
    if (!fEventCut.AcceptEvent(ev)) {
        PostAllData();
        return;
    }

    AliAODMCHeader *header{nullptr};
    TClonesArray *MCTrackArray{nullptr};
    if (fMC) {
        //OOB pileup
        header = static_cast<AliAODMCHeader *>(ev->FindListObject(AliAODMCHeader::StdBranchName()));
        if (!header) {
            AliWarning("No header found.");
            PostAllData();
            return;
        }
        MCTrackArray = dynamic_cast<TClonesArray *>(ev->FindListObject(AliAODMCParticle::StdBranchName()));
        if (!MCTrackArray) {
            AliWarning("No MC track array found.");
            PostAllData();
            return;
        }
    }

    auto pvObj = fEventCut.GetPrimaryVertex();
    double pv[3];
    pvObj->GetXYZ(pv);
    fRecHyper->fZ = pv[2];

    /// To perform the majority of the analysis - and also this one - the standard PID handler is
    /// required.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();
    fPID = handl->GetPIDResponse();

    fRecHyper->centrality = fEventCut.GetCentrality();

    unsigned char tgr{0u};

    if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)
        tgr |= kINT7;
    if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)
        tgr |= kCentral;
    if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
        tgr |= kSemiCentral;
    if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)
        tgr |= kHighMultV0;
    int magField = ev->GetMagneticField() > 0 ? kPositiveB : 0;

    fRecHyper->trigger = tgr + magField;

    std::vector<int> checkedHyperLabel;
    fGenHyper.isReconstructed = true;
    // V0 method
    if(fUseV0Method) {
        for (int iV0{0}; iV0 < ev->GetNumberOfV0s(); ++iV0) {
            AliAODv0 *v0{ev->GetV0(iV0)};
            if (!v0)
                continue;
            if (v0->GetOnFlyStatus() != fUseOnTheFly)
                continue;
            //get daughter tracks (positive, negative and bachelor)
            AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
            AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
            if (!pTrack || !nTrack) {
                AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the lambdas ...\n");
                continue;
            }

            if (!(pTrack->GetStatus() & AliVTrack::kTPCrefit) || !(nTrack->GetStatus() & AliVTrack::kTPCrefit) ||
                pTrack->GetTPCsignalN() < 50 || nTrack->GetTPCsignalN() < 50 ||
                std::abs(pTrack->Eta()) > 0.8 || std::abs(nTrack->Eta()) > 0.8) {
                continue;
            }

            int hyperLabel{-1};
            if (fMC) {
                fGenHyper.pdg = 0;
                auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(pTrack->GetLabel()));
                auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(nTrack->GetLabel()));
                if ((posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kAlpha) &&
                    negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kPion)) ||
                    (negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kAlpha) &&
                    posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kPion))) {
                    // Check hyper
                    int labMothPos = posPart->GetMother();
                    int labMothNeg = negPart->GetMother();
                    if (labMothNeg >= 0 && labMothNeg == labMothPos && !AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(labMothNeg, header, MCTrackArray)) {
                        auto hyper = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
                        if (hyper && std::abs(hyper->GetPdgCode()) == kHyperPdg) {
                            hyperLabel = labMothNeg;
                            fGenHyper.pdg = hyper->GetPdgCode();
                            fGenHyper.ptMC = hyper->Pt();
                            fGenHyper.etaMC = hyper->Eta();
                            fGenHyper.yMC = hyper->Y();
                            double ov[3], dv[3];
                            hyper->XvYvZv(ov);
                            posPart->XvYvZv(dv);
                            fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * hyper->M() / hyper->P();
                        }
                    }
                }

                if (fOnlyTrueCandidates && hyperLabel < 0)
                    continue;
            }

            double pNsigma{!fUseCustomPID ? fPID->NumberOfSigmasTPC(pTrack, AliPID::kAlpha)
                        : fMC          ? fPID->NumberOfSigmasTPC(pTrack, AliPID::kAlpha)
                                        : customNsigma(pTrack->GetTPCmomentum(), pTrack->GetTPCsignal())};
            double nNsigma{!fUseCustomPID ? fPID->NumberOfSigmasTPC(pTrack, AliPID::kAlpha)
                        : fMC          ? fPID->NumberOfSigmasTPC(nTrack, AliPID::kAlpha)
                                        : customNsigma(nTrack->GetTPCmomentum(), nTrack->GetTPCsignal())};
            if (std::abs(pNsigma) > 5 && std::abs(nNsigma) > 5) {
                continue;
            }
            if (std::abs(pNsigma) < 5 && std::abs(nNsigma) < 5) {
                continue;
            }
            fHistos->FillTH2("QA/hTPCPIDAlpha", pTrack->GetTPCmomentum(), pTrack->GetTPCsignal());
            fHistos->FillTH2("QA/hTPCPIDAntiAlpha", nTrack->GetTPCmomentum(), nTrack->GetTPCsignal());

            fRecHyper->Matter = std::abs(pNsigma) < 5;
            auto alpha = fRecHyper->Matter ? pTrack : nTrack;
            auto pion = fRecHyper->Matter ? nTrack : pTrack;

            double sv[3]{v0->GetSecVtxX(), v0->GetSecVtxY(), v0->GetSecVtxZ()};
            double deltaPos[3]{sv[0] - pv[0], sv[1] - pv[1], sv[2] - pv[2]};

            LVector_t alphaVector, piVector, hyperVector;
            double alphaP[3]{fRecHyper->Matter ? v0->MomPosX() : v0->MomNegX(), fRecHyper->Matter ? v0->MomPosY() : v0->MomNegY(), fRecHyper->Matter ? v0->MomPosZ() : v0->MomNegZ()};
            double piP[3]{!fRecHyper->Matter ? v0->MomPosX() : v0->MomNegX(), !fRecHyper->Matter ? v0->MomPosY() : v0->MomNegY(), !fRecHyper->Matter ? v0->MomPosZ() : v0->MomNegZ()};

            alphaVector.SetCoordinates(alphaP[0] * 2, alphaP[1] * 2, alphaP[2] * 2, AliPID::ParticleMass(AliPID::kAlpha));
            piVector.SetCoordinates(piP[0], piP[1], piP[2], AliPID::ParticleMass(AliPID::kPion));
            hyperVector = piVector + alphaVector;
            if (hyperVector.mass() > fMassRange[1] || hyperVector.mass() < fMassRange[0]) {
                continue;
            }

            double cpa = (deltaPos[0] * hyperVector.px() +
                        deltaPos[1] * hyperVector.py() +
                        deltaPos[2] * hyperVector.pz()) /
                        std::sqrt(hyperVector.P2() * (Sq(deltaPos[0]) + Sq(deltaPos[1]) + Sq(deltaPos[2])));

            fRecHyper->pt = hyperVector.pt();
            fRecHyper->m = hyperVector.mass();
            fRecHyper->V0CosPA = cpa;
            fRecHyper->Rapidity = Eta2y(fRecHyper->pt, kHyperMass, hyperVector.eta());

            fRecHyper->V0radius = v0->RadiusSecVtx();
            fRecHyper->Lrec = v0->DecayLengthV0(pv);
            fRecHyper->ct = fRecHyper->Lrec * kHyperMass / hyperVector.P();

            fRecHyper->alphaProngPvDCA = fRecHyper->Matter ? v0->DcaPosToPrimVertex() : v0->DcaNegToPrimVertex();
            fRecHyper->PiProngPvDCA = fRecHyper->Matter ? v0->DcaNegToPrimVertex() : v0->DcaPosToPrimVertex();
            float _dummy, xy;
            alpha->GetImpactParameters(xy, _dummy);
            fRecHyper->alphaProngPvDCAXY = xy;
            pion->GetImpactParameters(xy, _dummy);
            fRecHyper->PiProngPvDCAXY = xy;
            fRecHyper->ProngsDCA = v0->DcaV0Daughters();
            fRecHyper->TPCmomalpha = alpha->GetTPCmomentum();
            fRecHyper->TPCsignalalpha = alpha->GetTPCsignal();
            fRecHyper->NitsClustersalpha = alpha->GetITSNcls();
            fRecHyper->TPCnSigmaPi = fPID->NumberOfSigmasTPC(pion, AliPID::kPion);
            fRecHyper->TPCnSigmaalpha = fRecHyper->Matter ? pNsigma : nNsigma;
            fRecHyper->NpidClustersPion = pion->GetTPCsignalN();
            fRecHyper->NpidClustersalpha = alpha->GetTPCsignalN();

            if (hyperLabel != -1) {
                if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), hyperLabel) != checkedHyperLabel.end()) {
                    fGenHyper.isDuplicated = true;
                } else {
                    fGenHyper.isDuplicated = false;
                    checkedHyperLabel.push_back(hyperLabel);
                }
            }
            fTree->Fill();
        }
    }
    else {
        // Trackloop method (for trial)
        // iTrack loop -> alpha
        for (int iTrack{0}; iTrack < ev->GetNumberOfTracks(); ++iTrack) {
            AliAODTrack *alphaTrack = dynamic_cast<AliAODTrack *>(ev->GetTrack(iTrack));
            if (alphaTrack) {
                AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the lambdas ...\n");
                continue;
            }
            if (!(alphaTrack->GetStatus() & AliVTrack::kTPCrefit) || alphaTrack->GetTPCsignalN() < 50 
                || std::abs(alphaTrack->Eta()) > 0.8) {
                continue;
            }

            double pNsigma{fPID->NumberOfSigmasTPC(alphaTrack, AliPID::kAlpha)};
            if (std::abs(pNsigma) > 5 || std::abs(pNsigma) < 5) {
                continue;
            }
            if(alphaTrack->GetSign() > 0)
                fHistos->FillTH2("QA/hTPCPIDAlpha", alphaTrack->GetTPCmomentum(), alphaTrack->GetTPCsignal());
            else
                fHistos->FillTH2("QA/hTPCPIDAntiAlpha", alphaTrack->GetTPCmomentum(), alphaTrack->GetTPCsignal());

            // jTrack loop - pion
            for (int jTrack{0}; jTrack < ev->GetNumberOfTracks(); ++jTrack) {
                if (jTrack == iTrack)
                    continue;
                AliAODTrack *pionTrack = dynamic_cast<AliAODTrack *>(ev->GetTrack(jTrack));
                if (pionTrack) {
                AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the lambdas ...\n");
                    continue;
                }
                if (!(pionTrack->GetStatus() & AliVTrack::kTPCrefit) || pionTrack->GetTPCsignalN() < 50 
                    || std::abs(pionTrack->Eta()) > 0.8) {
                    continue;
                }
                double nNsigma{fPID->NumberOfSigmasTPC(pionTrack, AliPID::kPion)};
                if (std::abs(nNsigma) > 5 && std::abs(nNsigma) < 5) {
                    continue;
                }

                int hyperLabel{-1};
                if (fMC) {
                    fGenHyper.pdg = 0;
                    auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(alphaTrack->GetLabel()));
                    auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(pionTrack->GetLabel()));
                    if ((posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kAlpha) &&
                        negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kPion)) ||
                        (negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kAlpha) &&
                        posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kPion))) {
                        // Check hyper
                        int labMothPos = posPart->GetMother();
                        int labMothNeg = negPart->GetMother();
                        if (labMothNeg >= 0 && labMothNeg == labMothPos && !AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(labMothNeg, header, MCTrackArray)) {
                            auto hyper = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
                            if (hyper && std::abs(hyper->GetPdgCode()) == kHyperPdg) {
                                hyperLabel = labMothNeg;
                                fGenHyper.pdg = hyper->GetPdgCode();
                                fGenHyper.ptMC = hyper->Pt();
                                fGenHyper.etaMC = hyper->Eta();
                                fGenHyper.yMC = hyper->Y();
                                double ov[3], dv[3];
                                hyper->XvYvZv(ov);
                                posPart->XvYvZv(dv);
                                fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * hyper->M() / hyper->P();
                            }
                        }
                    }

                    if (fOnlyTrueCandidates && hyperLabel < 0)
                        continue;
                }


                fRecHyper->Matter = alphaTrack->GetSign() > 0;

                // Vertexing
                Double_t cv[21], dztemp[2], covartemp[3], pos[3], cov[6], chi2perNDF, dispersion;
                for (int i = 0; i < 21; i++)
                    cv[i] = 0;
                const AliESDVertex *vtxT3D = dynamic_cast<AliESDEvent *>(InputEvent())->GetPrimaryVertex();
                AliExternalTrackParam pionTrk(*(AliESDtrack *)ev->GetTrack(jTrack));
                AliExternalTrackParam alphaTrk(*(AliESDtrack *)ev->GetTrack(iTrack));
                AliExternalTrackParam *pPionTrk = &pionTrk, *pAlphaTrk = &alphaTrk;

                pPionTrk->PropagateToDCA(vtxT3D, magField, 250, dztemp, covartemp);
                pAlphaTrk->PropagateToDCA(vtxT3D, magField, 250, dztemp, covartemp);

                AliESDVertex *vertexESD = 0;
                AliAODVertex *vertexAOD = 0;

                AliKFParticle::SetField(magField);
                AliKFVertex vertexKF;

                AliESDtrack *esdTrackPion = (AliESDtrack *)pPionTrk;
                AliKFParticle daughterKFPion(*esdTrackPion, 211);
                vertexKF.AddDaughter(daughterKFPion);

                AliESDtrack *esdTrackAlpha = (AliESDtrack *)pAlphaTrk;
                AliKFParticle daughterKFXi(*esdTrackAlpha, 211);
                vertexKF.AddDaughter(daughterKFXi);

                vertexESD = new AliESDVertex(vertexKF.Parameters(),
                                            vertexKF.CovarianceMatrix(),
                                            vertexKF.GetChi2(),
                                            vertexKF.GetNContributors());

                vertexESD->GetXYZ(pos);        // position
                vertexESD->GetCovMatrix(cov);  //covariance matrix
                chi2perNDF = vertexESD->GetChi2toNDF();
                dispersion = vertexESD->GetDispersion();
                delete vertexESD;
                vertexESD = NULL;
                vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, 2);  // Hyper Vertex
                
                double sv[3];
                vertexAOD->GetXYZ(sv);  // Secondary vertex
                double deltaPos[3]{sv[0] - pv[0], sv[1] - pv[1], sv[2] - pv[2]};

                LVector_t alphaVector, piVector, hyperVector;
                double alphaP[3];
                alphaTrack->GetPxPyPz(alphaP);
                double piP[3];
                pionTrack->GetPxPyPz(piP);

                alphaVector.SetCoordinates(alphaP[0] * 2, alphaP[1] * 2, alphaP[2] * 2, AliPID::ParticleMass(AliPID::kAlpha));
                piVector.SetCoordinates(piP[0], piP[1], piP[2], AliPID::ParticleMass(AliPID::kPion));
                hyperVector = piVector + alphaVector;
                if (hyperVector.mass() > fMassRange[1] || hyperVector.mass() < fMassRange[0]) {
                    continue;
                }

                fRecHyper->pt = hyperVector.pt();
                fRecHyper->m = hyperVector.mass();
                fRecHyper->V0CosPA = (deltaPos[0] * hyperVector.px() +
                                      deltaPos[1] * hyperVector.py() +
                                      deltaPos[2] * hyperVector.pz()) /
                                      std::sqrt(hyperVector.P2() * (Sq(deltaPos[0]) + Sq(deltaPos[1]) + Sq(deltaPos[2])));
                fRecHyper->Rapidity = Eta2y(fRecHyper->pt, kHyperMass, hyperVector.eta());

                fRecHyper->V0radius = TMath::Hypot(sv[0] - pv[0], sv[1] - pv[1]);
                fRecHyper->Lrec = TMath::Sqrt(TMath::Power(sv[0] - pv[0],2) +
                                TMath::Power(sv[1] - pv[1],2) +
                                TMath::Power(sv[2] - pv[2],2 ));
                fRecHyper->ct = fRecHyper->Lrec * kHyperMass / hyperVector.P();
                Float_t b[2];     // Float due to the function input
                Float_t bCov[3];  // Float due to the function input
                alphaTrack->GetImpactParameters(b, bCov);
                fRecHyper->alphaProngPvDCA = b[0];
                pionTrack->GetImpactParameters(b, bCov);
                fRecHyper->PiProngPvDCA = b[0];
                float _dummy, xy;
                alphaTrack->GetImpactParameters(xy, _dummy);
                fRecHyper->alphaProngPvDCAXY = xy;
                pionTrack->GetImpactParameters(xy, _dummy);
                fRecHyper->PiProngPvDCAXY = xy;
                Double_t xdummy, ydummy;
                fRecHyper->ProngsDCA = pPionTrk->GetDCA(pAlphaTrk, magField, xdummy, ydummy);
                fRecHyper->TPCmomalpha = alphaTrack->GetTPCmomentum();
                fRecHyper->TPCsignalalpha = alphaTrack->GetTPCsignal();
                fRecHyper->NitsClustersalpha = alphaTrack->GetITSNcls();
                fRecHyper->TPCnSigmaPi = fPID->NumberOfSigmasTPC(pionTrack, AliPID::kPion);
                fRecHyper->TPCnSigmaalpha = fPID->NumberOfSigmasTPC(alphaTrack, AliPID::kAlpha);
                fRecHyper->NpidClustersPion = pionTrack->GetTPCsignalN();
                fRecHyper->NpidClustersalpha = alphaTrack->GetTPCsignalN();

                if (hyperLabel != -1) {
                    if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), hyperLabel) != checkedHyperLabel.end()) {
                        fGenHyper.isDuplicated = true;
                    } else {
                        fGenHyper.isDuplicated = false;
                        checkedHyperLabel.push_back(hyperLabel);
                    }
                }
                fTree->Fill();
            }
            
        }
    }
    if (fMC) {
        fGenHyper.isReconstructed = false;
        //loop on generated
        for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT) {
            auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
            int pdg = std::abs(track->GetPdgCode());
            if (pdg != kHyperPdg) {
                continue;
            }
            if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), iT) != checkedHyperLabel.end()) {
                continue;
            }

            if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray)) {
                continue;
            }
            fGenHyper.ptMC = track->Pt();
            fGenHyper.etaMC = track->Eta();
            fGenHyper.yMC = track->Y();
            fGenHyper.pdg = track->GetPdgCode();
            double ov[3], dv[3];
            track->XvYvZv(ov);
            bool otherDecayChannel{true};
            for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++) {
                auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
                if (!daugh) {
                    continue;
                }
                if (std::abs(daugh->GetPdgCode()) == AliPID::ParticleCode(AliPID::kAlpha)) {
                    otherDecayChannel = false;
                    daugh->XvYvZv(dv);
                    break;
                }
            }
            if (otherDecayChannel)
                continue;
            fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * track->M() / track->P();
            fTree->Fill();
        }
    }

    PostAllData();
}

AliAnalysisTaskAlphaPiAOD *AliAnalysisTaskAlphaPiAOD::AddTask(bool isMC, TString tskname, TString suffix) {
    // Get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskHe4piAOD", "No analysis manager found.");
        return nullptr;
    }

    // Check the analysis type using the event handlers connected to the analysis
    // manager.
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskHe4piAOD", "This task requires an input event handler");
        return nullptr;
    }

    tskname.Append(suffix.Data());
    AliAnalysisTaskAlphaPiAOD *task = new AliAnalysisTaskAlphaPiAOD(isMC, tskname.Data());

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
        Form("%s_output", tskname.Data()), TList::Class(),
        AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

    AliAnalysisDataContainer *coutput2 =
        mgr->CreateContainer(Form("%s_treeHyper", tskname.Data()), TTree::Class(),
                             AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    coutput2->SetSpecialOutput();

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1);
    mgr->ConnectOutput(task, 2, coutput2);
    return task;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskAlphaPiAOD::Eta2y(float pt, float m, float eta) const {
    return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}

void AliAnalysisTaskAlphaPiAOD::PostAllData() {
    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fTree);
}

void AliAnalysisTaskAlphaPiAOD::SetCustomBetheBloch(float res, const float *bethe) {
    fCustomResolution = res;
    std::copy(bethe, bethe + 5, fCustomBethe);
}

double AliAnalysisTaskAlphaPiAOD::customNsigma(double mom, double sig) {
    const float bg = mom / AliPID::ParticleMass(AliPID::kAlpha);
    const float *p = fCustomBethe;
    const float expS = AliExternalTrackParam::BetheBlochAleph(bg, p[0], p[1], p[2], p[3], p[4]);
    return (sig - expS) / (fCustomResolution * expS);
}