/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnalysisTaskSigma1385PM
 *
 *  Test code for the reconstructing Sigma(1385)^{+-}
 *  Output will be saved to nTuple -> can be used for TMVA input
 *
 *  Author: Bong-Hwi Lim
 *
 */

#include <TDatabasePDG.h>
#include <math.h>
#include <iostream>
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "TChain.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAODv0.h"
#include "AliAnalysisTaskSigma1385PM.h"
#include "AliESDv0.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "THistManager.h"

const Double_t pi = TMath::Pi();
const Double_t pionMass = AliPID::ParticleMass(AliPID::kPion);
const Double_t v0Mass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

enum {
    kSigmaStarPCode = 3224,  // Sigma(1385)+
    kSigmaStarNCode = 3114,  // Sigma(1385)-
    kLambdaCode = 3122,      // Lambda
    kProtonCode = 2212,      // Proton+
    kPionCode = 211,         // Pion+
};
enum {
    kSigmaStarP = 1,
    kSigmaStarN,
    kAntiSigmaStarP,
    kAntiSigmaStarN,
    kSigmaStarP_MIX,  // 5
    kSigmaStarN_MIX,
    kAntiSigmaStarP_MIX,
    kAntiSigmaStarN_MIX,
    kSigmaStarP_GEN,  // 9
    kSigmaStarN_GEN,
    kAntiSigmaStarP_GEN,
    kAntiSigmaStarN_GEN,
    kSigmaStarP_REC,  // 13
    kSigmaStarN_REC,
    kAntiSigmaStarP_REC,
    kAntiSigmaStarN_REC,
    kAllType
};

class AliAnalysisTaskSigma1385PM;

ClassImp(AliAnalysisTaskSigma1385PM)  // classimp: necessary for root

    AliAnalysisTaskSigma1385PM::AliAnalysisTaskSigma1385PM()
    : AliAnalysisTaskSE(), fEvt(0), fNtupleSigma1385(0) {}
//_____________________________________________________________________________
AliAnalysisTaskSigma1385PM::AliAnalysisTaskSigma1385PM(const char* name,
                                                           Bool_t MCcase)
    : AliAnalysisTaskSE(name), fEvt(0), IsMC(MCcase), fNtupleSigma1385(0) {
    // constructor
    DefineInput(
        0, TChain::Class());  // define the input of the analysis: in this case
                              // we take a 'chain' of events this chain is
                              // created by the analysis manager, so no need to
                              // worry about it, it does its work automatically
    DefineOutput(1, TList::Class());  // define the ouptut of the analysis: in
                                      // this case it's a list of histograms
    DefineOutput(
        2, TNtupleD::Class());  // you can add more output objects by calling
                                // DefineOutput(2, classname::Class())
}
//_____________________________________________________________________________
AliAnalysisTaskSigma1385PM::~AliAnalysisTaskSigma1385PM() {}
//___________________________________________________________________
void AliAnalysisTaskSigma1385PM::SetCutOpen() {
    // Pion cuts
    SetFilterbitSigmaStarPion(1);
    SetMaxNsigSigmaStarPion(5);
    SetMaxEtaSigmaStarPion(0.8);
    SetMaxVertexZSigmaStarPion(99);
    SetMaxVertexXYsigSigmaStarPion(99);

    // Lambda cuts
    SetMaxNsigV0Proton(5);
    SetMaxNsigV0Pion(5);
    SetMaxDCAV0daughters(999);
    SetMaxDCAPVV0(999);
    SetMinCPAV0(0.9);
    SetMaxMassWindowV0(999);

    // Sigma Star cut
    SetSigmaStarRapidityCutHigh(1);
    SetSigmaStarRapidityCutLow(-1);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma1385PM::UserCreateOutputObjects() {
    fTrackCuts = new AliESDtrackCuts();
    fTrackCuts->GetStandardITSTPCTrackCuts2011(kTRUE, kTRUE);
    fTrackCuts->SetEtaRange(-0.8, 0.8);
    fTrackCuts->SetPtRange(0.15, 1e20);

    fHistos = new THistManager("Sigma1385hists");
    fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());

    auto binType = AxisStr(
        "Type", {"SigmaStarP", "SigmaStarN", "AntiSigmaStarP", "AntiSigmaStarN",
                 "SigmaStarP_mix", "SigmaStarN_mix", "AntiSigmaStarP_mix",
                 "AntiSigmaStarN_mix", "SigmaStarP_gen", "SigmaStarN_gen",
                 "AntiSigmaStarP_gen", "AntiSigmaStarN_gen", "SigmaStarP_rec",
                 "SigmaStarN_rec", "AntiSigmaStarP_rec", "AntiSigmaStarN_rec"});

    std::vector<double> centaxisbin = {
        0,  1,  5,  10, 15, 20, 30,
        40, 50, 60, 70, 80, 90, 100};  // can be use from pp to PbPb
    binCent = AxisVar("Cent", centaxisbin);
    auto binPt = AxisFix("Pt", 200, 0, 20);
    auto binMass = AxisFix("Mass", 2000, 1.0, 3.0);
    binZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10});

    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 100, "s");
    CreateTHnSparse("Sigma1385", "Sigma1385", 4,
                    {binType, binCent, binPt, binMass}, "s");

    fEMpool.resize(binCent.GetNbins() + 1,
                   std::vector<eventpool>(binZ.GetNbins() + 1));

    fNtupleSigma1385 = new TNtupleD(
        "fNtupleSigma1385", "Sigma1385",
        "PIDSigmaStarPion:DCASigmaStarPionToPrimVertexZ:"
        "DCASigmaStarPionToPrimVertexR:EtaSigmaStarPion:PhiSigmaStarPion:"
        "PIDV0pTrackProton:PIDV0pTrackPion:PIDV0nTrackProton:PIDV0pTrackPion:"
        "DCAV0Daughters:DCAV0ToPrimVertex:CosPointingAngleV0:V0Mass:EtaV0:"
        "PhiV0:MCflag");

    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fNtupleSigma1385);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma1385PM::UserExec(Option_t*) {
    AliVEvent* event = InputEvent();
    if (!event) {
        PostData(1, fHistos->GetListOfHistograms());
        PostData(2, fNtupleSigma1385);
        AliInfo("Could not retrieve event");
        return;
    }
    AliNanoAODHeader* nanoHeader =
        dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());

    event->IsA() == AliESDEvent::Class()
        ? fEvt = dynamic_cast<AliESDEvent*>(event)
        : fEvt = dynamic_cast<AliAODEvent*>(event);
    if (!fEvt) {
        PostData(1, fHistos->GetListOfHistograms());
        PostData(2, fNtupleSigma1385);
        return;
    }
    AliInputEventHandler* inputHandler =
        (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();

    bool IsEvtSelected{false}, IsINEL0True{false};
    if (!nanoHeader) {
        IsEvtSelected = fEventCuts.AcceptEvent(event);
        if (IsMC) {
            if (fEvt->IsA() != AliESDEvent::Class())
                fMCArray = (TClonesArray*)fEvt->FindListObject(
                    "mcparticles");  // AOD Case
            fMCEvent = MCEvent();
            IsINEL0True = fEventCuts.IsTrueINELgtZero(fEvt, true);
        }
        fCent = fEventCuts.GetCentrality(0);
        fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
        if (!fPIDResponse)
            AliInfo("No PIDd");
    } else {
        IsEvtSelected = true;
        fCent = nanoHeader->GetCentr("V0M");
    }

    if (!IsEvtSelected) {
        PostData(1, fHistos->GetListOfHistograms());
        PostData(2, fNtupleSigma1385);
        return;  // event cut
    }

    fHistos->FillTH1("hMultiplicity", (double)fCent);

    if (IsMC)
        FillMCinput(fMCEvent);
    if (fEvt->IsA() == AliAODEvent::Class())
        vertex = ((AliAODEvent*)fEvt)->GetPrimaryVertex();
    const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
    const AliVVertex* spdVtx = fEvt->GetPrimaryVertexSPD();
    PVx = pVtx->GetX();
    PVy = pVtx->GetY();
    PVz = pVtx->GetZ();
    fZ = spdVtx->GetZ();

    // Event Mixing pool -----------------------------------------------------
    zbin = binZ.FindBin(fZ) - 1;           // Event mixing z-bin
    centbin = binCent.FindBin(fCent) - 1;  // Event mixing cent bin

    bool checkPion = GoodTracksSelection();
    bool checkV0 = GoodV0Selection();

    if (checkPion && checkV0) {
        FillTracks();  // Fill the histogram
        if (fFillnTuple)
            FillNtuples();
    }
    if (fsetmixing && goodtrackindices.size()) {
        FillTrackToEventPool();  // use only pion track pool.
    }
    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fNtupleSigma1385);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma1385PM::Terminate(Option_t*) {}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigma1385PM::GoodTracksSelection() {
    const UInt_t nTracks = fEvt->GetNumberOfTracks();
    goodtrackindices.clear();
    AliVTrack* track;

    for (UInt_t it = 0; it < nTracks; it++) {
        if (fEvt->IsA() == AliESDEvent::Class()) {
            track = (AliESDtrack*)fEvt->GetTrack(it);
            if (!track)
                continue;
            if (!fTrackCuts->AcceptTrack((AliESDtrack*)track))
                continue;
        }  // ESD Case
        else {
            track = (AliAODTrack*)fEvt->GetTrack(it);
            if (!track)
                continue;
            if (!((AliAODTrack*)track)->TestFilterBit(fFilterBit))
                continue;
        }  // AOD Case

        Double_t fTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
        Double_t pionZ = abs(track->GetZ() - fZ);
        Double_t pionPt = track->Pt();
        Double_t pionSigmaDCA_r = (0.0026 + 0.0050 / pionPt);
        Double_t pionDCA_r =
            sqrt(pow(track->GetX(), 2) + pow(track->GetY(), 2));

        if (abs(fTPCNSigPion) > fTPCNsigSigmaStarPionCut)
            continue;
        if (abs(track->Eta()) > fSigmaStarPionEtaCut)
            continue;
        if (pionPt < 0.15)
            continue;
        if (pionZ > fSigmaStarPionZVertexCut)
            continue;
        if (pionDCA_r > pionSigmaDCA_r * fSigmaStarPionXYVertexSigmaCut)
            continue;

        goodtrackindices.push_back(it);
    }
    return goodtrackindices.size();
}
Bool_t AliAnalysisTaskSigma1385PM::GoodV0Selection() {
    goodv0indices.clear();
    const UInt_t nV0 = fEvt->GetNumberOfV0s();

    AliESDv0* v0ESD;
    AliAODv0* v0AOD;
    Double_t LambdaX, LambdaY, LambdaZ;
    Double_t fDCADist_LambdaProton_PV, fDCADist_LambdaPion_PV;

    Bool_t AcceptedV0 = kTRUE;
    if (fEvt->IsA() == AliESDEvent::Class()) {  // ESD case
        for (UInt_t it = 0; it < nV0; it++) {
            AcceptedV0 = kTRUE;
            v0ESD = ((AliESDEvent*)fEvt)->GetV0(it);
            if (!v0ESD)
                continue;

            if (TMath::Abs(v0ESD->GetPindex()) ==
                TMath::Abs(v0ESD->GetNindex()))
                continue;

            AliESDtrack* pTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
            AliESDtrack* nTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
            if (TMath::Abs(((pTrackV0->GetSign()) - (nTrackV0->GetSign()))) <
                0.1)
                AcceptedV0 = kFALSE;

            // PID cuts
            Double_t fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            Double_t fTPCNSigAntiProton =
                GetTPCnSigma(nTrackV0, AliPID::kProton);
            Double_t fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            Double_t fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            Bool_t fPIDLambda =
                (abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigPion) < fTPCNsigLambdaPionCut);
            Bool_t fPIDAntiLambda =
                (abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut);

            if (fPIDLambda)
                v0ESD->ChangeMassHypothesis(kLambda0);
            if (fPIDAntiLambda)
                v0ESD->ChangeMassHypothesis(kLambda0Bar);

            if (!fPIDLambda && !fPIDAntiLambda)
                AcceptedV0 = kFALSE;

            // DCA cut
            // DCA between Dautgher particles
            Double_t fDCADistLambda = fabs(v0ESD->GetDcaV0Daughters());

            if (fDCADistLambda > fDCADistLambdaDaughtersCut)
                AcceptedV0 = kFALSE;  // DCA proton-pion

            // DCA to PV
            Double_t fDCADistLambda_PV = fabs(v0ESD->GetD(PVx, PVy, PVz));

            if (fDCADistLambda_PV > fDCArDistLambdaPVCut)
                AcceptedV0 = kFALSE;

            // CPA cut
            Double_t fLambdaCPA =
                v0ESD->GetV0CosineOfPointingAngle(PVx, PVy, PVz);

            if (fLambdaCPA < fV0CosineOfPointingAngleCut)
                AcceptedV0 = kFALSE;

            // Mass window cut
            Double_t fMassV0 = v0ESD->GetEffMass();

            if (fabs(fMassV0 - v0Mass) > fV0MassWindowCut)
                AcceptedV0 = kFALSE;

            // After selection above
            if (AcceptedV0)
                goodv0indices.push_back(it);  // for standard V0
        }                                     // All V0 loop
    }                                         // ESD case
    else {
        for (UInt_t it = 0; it < nV0; it++) {
            AcceptedV0 = kTRUE;
            v0AOD = ((AliAODEvent*)fEvt)->GetV0(it);
            if (!v0AOD)
                continue;

            if (TMath::Abs(v0AOD->GetPosID()) == TMath::Abs(v0AOD->GetNegID()))
                continue;

            AliAODTrack* pTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack* nTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

            // PID cuts
            Double_t fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            Double_t fTPCNSigAntiProton =
                GetTPCnSigma(nTrackV0, AliPID::kProton);
            Double_t fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            Double_t fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            Bool_t fPIDLambda =
                (abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigPion) < fTPCNsigLambdaPionCut);
            Bool_t fPIDAntiLambda =
                (abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut);

            if (!fPIDLambda && !fPIDAntiLambda)
                AcceptedV0 = kFALSE;

            // DCA cut
            Double_t fDCADistLambda = fabs(v0AOD->DcaV0Daughters());
            if (fDCADistLambda > fDCADistLambdaDaughtersCut)
                AcceptedV0 = kFALSE;  // DCA proton-pion

            // DCA to PV
            Double_t fDCADistLambda_PV = fabs(v0AOD->DcaV0ToPrimVertex());

            if (fDCADistLambda_PV > fDCArDistLambdaPVCut)
                AcceptedV0 = kFALSE;

            Double_t fLambdaCPA = v0AOD->CosPointingAngle(vertex);

            if (fLambdaCPA < fV0CosineOfPointingAngleCut)
                AcceptedV0 = kFALSE;

            // Mass window cut
            Double_t fMassV0 = v0AOD->MassLambda();
            if (fabs(fMassV0 - v0Mass) > fV0MassWindowCut)
                AcceptedV0 = kFALSE;

            // After selection above
            if (AcceptedV0) {
                goodv0indices.push_back(it);
            }
        }  // All v0 loop
    }      // AOD case

    return goodv0indices.size();
}
void AliAnalysisTaskSigma1385PM::FillTracks() {
    AliVTrack* track1;
    AliESDv0* v0ESD;
    AliAODv0* v0AOD;
    Double_t fTPCNSigProton, fTPCNSigAntiProton, fTPCNSigPion, fTPCNSigAntiPion;
    Bool_t fPIDLambda, fPIDAntiLambda;
    Bool_t SkipMixing = kFALSE;
    Int_t pID, nID;

    TLorentzVector temp1, temp2;
    TLorentzVector vecsum;
    const UInt_t nV0 = goodv0indices.size();
    const UInt_t nTracks = goodtrackindices.size();

    tracklist trackpool;
    if (fsetmixing) {
        eventpool& ep = fEMpool[centbin][zbin];
        if ((int)ep.size() < (int)fnMix)
            SkipMixing = kTRUE;
        if (!SkipMixing) {
            for (auto pool : ep) {
                for (auto track : pool)
                    trackpool.push_back((AliVTrack*)track);
            }
        }
    }
    for (UInt_t i = 0; i < nV0; i++) {
        if (fEvt->IsA() == AliESDEvent::Class()) {
            v0ESD = ((AliESDEvent*)fEvt)->GetV0(goodv0indices[i]);
            if (!v0ESD)
                continue;
            AliESDtrack* pTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
            AliESDtrack* nTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();
            // Lambda check
            fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            fTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
            fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            fPIDLambda = (abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) &&
                         (abs(fTPCNSigPion) < fTPCNsigLambdaPionCut);
            fPIDAntiLambda =
                (abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut);
            if (fPIDLambda)
                v0ESD->ChangeMassHypothesis(kLambda0);
            if (fPIDAntiLambda)
                v0ESD->ChangeMassHypothesis(kLambda0Bar);

            temp1.SetXYZM(v0ESD->Px(), v0ESD->Py(), v0ESD->Pz(),
                          v0ESD->GetEffMass());
        } else {
            v0AOD = ((AliAODEvent*)fEvt)->GetV0(goodv0indices[i]);
            if (!v0AOD)
                continue;
            AliAODTrack* pTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack* nTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();
            // Lambda check
            fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            fTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
            fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            fPIDLambda = (abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) &&
                         (abs(fTPCNSigPion) < fTPCNsigLambdaPionCut);
            fPIDAntiLambda =
                (abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut);

            if (fPIDLambda)
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassLambda());
            if (fPIDAntiLambda)
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassAntiLambda());
        }

        for (UInt_t j = 0; j < nTracks; j++) {
            track1 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[j]);
            if (!track1)
                continue;

            if (track1->GetID() == pID || track1->GetID() == nID)
                continue;

            temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), pionMass);

            vecsum = temp1 + temp2;  // temp1 = cascade, temp2=pion
            // Y cut
            if ((vecsum.Rapidity() > fSigmaStarYCutHigh) ||
                (vecsum.Rapidity() < fSigmaStarYCutLow))
                continue;

            auto sign = kAllType;

            Bool_t pionN = track1->Charge() == -1;
            Bool_t pionP = track1->Charge() == +1;

            if (fPIDLambda && pionP)
                sign = kSigmaStarP;
            if (fPIDLambda && pionN)
                sign = kSigmaStarN;
            if (fPIDAntiLambda && pionP)
                sign = kAntiSigmaStarN;
            if (fPIDAntiLambda && pionN)
                sign = kAntiSigmaStarP;

            FillTHnSparse("Sigma1385", {(double)sign, (double)fCent,
                                        vecsum.Pt(), vecsum.M()});

            if (IsMC &&
                IsTrueSigmaStar(goodv0indices[i], goodtrackindices[j])) {
                if (fPIDLambda && pionP)
                    sign = kSigmaStarP_REC;
                if (fPIDLambda && pionN)
                    sign = kSigmaStarN_REC;
                if (fPIDAntiLambda && pionP)
                    sign = kAntiSigmaStarN_REC;
                if (fPIDAntiLambda && pionN)
                    sign = kAntiSigmaStarP_REC;

                FillTHnSparse("Sigma1385", {(double)sign, (double)fCent,
                                            vecsum.Pt(), vecsum.M()});
            }
        }  // pion loop

        if ((centbin >= 0) && (zbin >= 0) && fsetmixing && !SkipMixing) {
            auto sign = kAllType;
            for (UInt_t jt = 0; jt < trackpool.size(); jt++) {
                track1 = trackpool.at(jt);
                if (track1->GetID() == pID || track1->GetID() == nID)
                    continue;
                temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(),
                              pionMass);
                vecsum = temp1 + temp2;
                // Y cut
                if ((vecsum.Rapidity() > fSigmaStarYCutHigh) ||
                    (vecsum.Rapidity() < fSigmaStarYCutLow))
                    continue;

                Bool_t pionN = track1->Charge() == -1;
                Bool_t pionP = track1->Charge() == +1;

                if (fPIDLambda && pionP)
                    sign = kSigmaStarP_MIX;
                if (fPIDLambda && pionN)
                    sign = kSigmaStarN_MIX;
                if (fPIDAntiLambda && pionP)
                    sign = kAntiSigmaStarN_MIX;
                if (fPIDAntiLambda && pionN)
                    sign = kAntiSigmaStarP_MIX;

                FillTHnSparse("Sigma1385", {(double)sign, (double)fCent,
                                            vecsum.Pt(), vecsum.M()});
            }
        }
    }
}
void AliAnalysisTaskSigma1385PM::FillNtuples() {
    AliVTrack* track1 = nullptr;
    AliESDv0* v0ESD = nullptr;
    AliAODv0* v0AOD = nullptr;
    Double_t fTPCNSigProton, fTPCNSigAntiProton, fTPCNSigPion, fTPCNSigAntiPion;
    Bool_t fPIDLambda, fPIDAntiLambda;
    Int_t pID, nID;
    Double_t tmp[16];
    for (UInt_t i = 0; i < 16; i++)
        tmp[i] = -999;  // initial value

    TLorentzVector temp1, temp2;
    TLorentzVector vecsum;
    const UInt_t nV0 = goodv0indices.size();
    const UInt_t nTracks = goodtrackindices.size();
    for (UInt_t i = 0; i < nV0; i++) {
        if (fEvt->IsA() == AliESDEvent::Class()) {
            v0ESD = ((AliESDEvent*)fEvt)->GetV0(goodv0indices[i]);
            if (!v0ESD)
                continue;
            AliESDtrack* pTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
            AliESDtrack* nTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();
            // Lambda check
            fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            fTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
            fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            fPIDLambda = (abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) &&
                         (abs(fTPCNSigPion) < fTPCNsigLambdaPionCut);
            fPIDAntiLambda =
                (abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut);

            if (fPIDLambda)
                v0ESD->ChangeMassHypothesis(kLambda0);
            if (fPIDAntiLambda)
                v0ESD->ChangeMassHypothesis(kLambda0Bar);

            temp1.SetXYZM(v0ESD->Px(), v0ESD->Py(), v0ESD->Pz(),
                          v0ESD->GetEffMass());

            // nTuple
            tmp[5] = fTPCNSigProton;                     // PIDV0pTrackProton
            tmp[6] = fTPCNSigAntiProton;                 // PIDV0pTrackPion
            tmp[7] = fTPCNSigPion;                       // PIDV0nTrackProton
            tmp[8] = fTPCNSigAntiPion;                   // PIDV0nTrackPion
            tmp[9] = fabs(v0ESD->GetDcaV0Daughters());   // DCAV0Daughters
            tmp[10] = fabs(v0ESD->GetD(PVx, PVy, PVz));  // DCAV0ToPrimVertex
            tmp[11] = v0ESD->GetV0CosineOfPointingAngle(
                PVx, PVy, PVz);             // CosPointingAngleV0
            tmp[12] = v0ESD->GetEffMass();  // V0Mass
            tmp[13] = v0ESD->Eta();         // EtaV0
            tmp[14] = v0ESD->Phi();         // PhiV0
        } else {
            v0AOD = ((AliAODEvent*)fEvt)->GetV0(goodv0indices[i]);
            if (!v0AOD)
                continue;
            AliAODTrack* pTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack* nTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();
            // Lambda check
            fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            fTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
            fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            fPIDLambda = (abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) &&
                         (abs(fTPCNSigPion) < fTPCNsigLambdaPionCut);
            fPIDAntiLambda =
                (abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) &&
                (abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut);

            if (fPIDLambda)
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassLambda());
            if (fPIDAntiLambda)
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassAntiLambda());
            // nTuple
            tmp[5] = fTPCNSigProton;                     // PIDV0pTrackProton
            tmp[6] = fTPCNSigAntiProton;                 // PIDV0pTrackPion
            tmp[7] = fTPCNSigPion;                       // PIDV0nTrackProton
            tmp[8] = fTPCNSigAntiPion;                   // PIDV0nTrackPion
            tmp[9] = fabs(v0AOD->DcaV0Daughters());      // DCAV0Daughters
            tmp[10] = fabs(v0AOD->DcaV0ToPrimVertex());  // DCAV0ToPrimVertex
            tmp[11] = v0AOD->CosPointingAngle(vertex);   // CosPointingAngleV0
            if (fPIDLambda)
                tmp[12] = v0AOD->MassLambda();  // V0Mass
            if (fPIDAntiLambda)
                tmp[12] = v0AOD->MassAntiLambda();  // V0Mass
            tmp[13] = v0AOD->Eta();                 // EtaV0
            tmp[14] = v0AOD->Phi();                 // PhiV0
        }

        for (UInt_t j = 0; j < nTracks; j++) {
            track1 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[j]);
            if (!track1)
                continue;

            if (track1->GetID() == pID || track1->GetID() == nID)
                continue;

            temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), pionMass);

            vecsum = temp1 + temp2;  // temp1 = cascade, temp2=pion
            // Y cut
            if ((vecsum.Rapidity() > fSigmaStarYCutHigh) ||
                (vecsum.Rapidity() < fSigmaStarYCutLow))
                continue;

            auto sign = kAllType;

            Bool_t pionN = track1->Charge() == -1;
            Bool_t pionP = track1->Charge() == +1;

            if (fPIDLambda && pionP)
                sign = kSigmaStarP;
            if (fPIDLambda && pionN)
                sign = kSigmaStarN;
            if (fPIDAntiLambda && pionP)
                sign = kAntiSigmaStarN;
            if (fPIDAntiLambda && pionN)
                sign = kAntiSigmaStarP;

            tmp[0] = GetTPCnSigma(track1, AliPID::kPion);  // PIDSigmaStarPion
            tmp[1] = abs(track1->GetZ() - fZ);  // DCASigmaStarPionToPrimVertexZ
            tmp[2] =
                sqrt(pow(track1->GetX(), 2) +
                     pow(track1->GetY(), 2));  // DCASigmaStarPionToPrimVertexR
            tmp[3] = track1->Eta();            // EtaSigmaStarPion
            tmp[4] = track1->Phi();            // PhiSigmaStarPion
            // nTuple -> See above
            // tmp[5] = fTPCNSigProton;//PIDV0pTrackProton
            // tmp[6] = fTPCNSigAntiProton;//PIDV0pTrackPion
            // tmp[7] = fTPCNSigPion;//PIDV0nTrackProton
            // tmp[8] = fTPCNSigAntiPion;//PIDV0nTrackPion
            // tmp[9] = fabs(v0ESD->GetDcaV0Daughters()); //DCAV0Daughters
            // tmp[10] = fabs(v0ESD->GetD(PVx, PVy, PVz)); //DCAV0ToPrimVertex
            // tmp[11] = v0ESD->GetV0CosineOfPointingAngle(PVx, PVy,
            // PVz);//CosPointingAngleV0 tmp[12] = v0ESD->GetEffMass(); //V0Mass
            // tmp[13] = v0ESD->Eta(); //EtaV0
            // tmp[14] = v0ESD->Phi(); //PhiV0

            if (IsMC) {
                if (IsTrueSigmaStar(goodv0indices[i], goodtrackindices[j]))
                    tmp[15] = (int)sign;  // MCflag
                else
                    tmp[15] = 5;  // MCflag -> not true
            } else
                tmp[15] = 0;  // MCflag -> data
        }                     // pion loop
        fNtupleSigma1385->Fill(tmp);
    }
}
void AliAnalysisTaskSigma1385PM::FillMCinput(AliMCEvent* fMCEvent) {
    auto sign = kAllType;
    if (fEvt->IsA() == AliESDEvent::Class()) {
        for (Int_t it = 0; it < fMCEvent->GetNumberOfPrimaries(); it++) {
            TParticle* mcInputTrack =
                (TParticle*)fMCEvent->GetTrack(it)->Particle();
            if (!mcInputTrack) {
                Error("UserExec", "Could not receive MC track %d", it);
                continue;
            }
            Int_t v0PdgCode = mcInputTrack->GetPdgCode();
            if ((abs(v0PdgCode) != kSigmaStarPCode) &&
                (abs(v0PdgCode) != kSigmaStarNCode))
                continue;
            if (IsPrimaryMC && !mcInputTrack->IsPrimary())
                continue;
            // Y cut
            if ((mcInputTrack->Y() > fSigmaStarYCutHigh) ||
                (mcInputTrack->Y() < fSigmaStarYCutLow))
                continue;
            if (v0PdgCode == kSigmaStarPCode)
                sign = kSigmaStarP_GEN;
            if (v0PdgCode == -kSigmaStarPCode)
                sign = kAntiSigmaStarP_GEN;
            if (v0PdgCode == kSigmaStarNCode)
                sign = kSigmaStarN_GEN;
            if (v0PdgCode == -kSigmaStarNCode)
                sign = kAntiSigmaStarN_GEN;

            FillTHnSparse("Sigma1385",
                          {(double)sign, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        }
    } else {
        for (Int_t it = 0; it < fMCArray->GetEntriesFast(); it++) {
            AliAODMCParticle* mcInputTrack =
                (AliAODMCParticle*)fMCArray->At(it);
            if (!mcInputTrack) {
                Error("UserExec", "Could not receive MC track %d", it);
                continue;
            }

            Int_t v0PdgCode = mcInputTrack->GetPdgCode();

            if ((abs(v0PdgCode) != kSigmaStarPCode) &&
                (abs(v0PdgCode) != kSigmaStarNCode))
                continue;
            if (IsPrimaryMC && !mcInputTrack->IsPrimary())
                continue;

            // Y cut
            if ((mcInputTrack->Y() > fSigmaStarYCutHigh) ||
                (mcInputTrack->Y() < fSigmaStarYCutLow))
                continue;

            if (v0PdgCode == kSigmaStarPCode)
                sign = kSigmaStarP_GEN;
            if (v0PdgCode == -kSigmaStarPCode)
                sign = kAntiSigmaStarP_GEN;
            if (v0PdgCode == kSigmaStarNCode)
                sign = kSigmaStarN_GEN;
            if (v0PdgCode == -kSigmaStarNCode)
                sign = kAntiSigmaStarN_GEN;

            FillTHnSparse("Sigma1385",
                          {(double)sign, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        }
    }
}
Bool_t AliAnalysisTaskSigma1385PM::IsTrueSigmaStar(UInt_t v0Index,
                                                     UInt_t pionIndex) {
    Bool_t trueSigmaStar = kFALSE;
    AliVTrack* track1;
    AliESDv0* v0ESD;
    AliAODv0* v0AOD;

    track1 = (AliVTrack*)fEvt->GetTrack(pionIndex);

    if (fEvt->IsA() == AliESDEvent::Class()) {
        v0ESD = ((AliESDEvent*)fEvt)->GetV0(v0Index);
        if (!v0ESD)
            return kFALSE;
        AliESDtrack* pTrackV0 =
            ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
        AliESDtrack* nTrackV0 =
            ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));

        TParticle* MCLamD1 =
            (TParticle*)fMCEvent->GetTrack(abs(pTrackV0->GetLabel()))
                ->Particle();
        TParticle* MCLamD2 =
            (TParticle*)fMCEvent->GetTrack(abs(nTrackV0->GetLabel()))
                ->Particle();
        TParticle* MCLam;
        TParticle* MCPi;
        TParticle* MCSigmaStar;
        if ((abs(MCLamD1->GetPdgCode()) == kProtonCode &&
             abs(MCLamD2->GetPdgCode()) == kPionCode) ||
            (abs(MCLamD1->GetPdgCode()) == kPionCode &&
             abs(MCLamD2->GetPdgCode()) == kProtonCode)) {
            if (MCLamD1->GetMother(0) == MCLamD2->GetMother(0)) {
                MCLam =
                    (TParticle*)fMCEvent->GetTrack(abs(MCLamD1->GetMother(0)))
                        ->Particle();
                if (abs(MCLam->GetPdgCode()) == kLambdaCode) {
                    MCPi =
                        (TParticle*)fMCEvent->GetTrack(abs(track1->GetLabel()))
                            ->Particle();
                    if (abs(MCPi->GetPdgCode()) == kPionCode) {
                        if (MCPi->GetMother(0) == MCLam->GetMother(0)) {
                            MCSigmaStar =
                                (TParticle*)fMCEvent
                                    ->GetTrack(abs(MCLam->GetMother(0)))
                                    ->Particle();
                            if ((abs(MCSigmaStar->GetPdgCode()) ==
                                 kSigmaStarPCode) ||
                                (abs(MCSigmaStar->GetPdgCode()) ==
                                 kSigmaStarNCode)) {
                                if (IsPrimaryMC) {
                                    if (MCSigmaStar->IsPrimary()) {
                                        trueSigmaStar = kTRUE;
                                    }  // Primary(input) SigmaStar check
                                } else {
                                    trueSigmaStar = kTRUE;
                                }
                            }  // Sigma Star check
                        }      // Pion mother = Lambda mother
                    }          // pion check
                }              // Lambda check
            }                  // Lambda duather's mother check
        }                      // Lambda daughter (pion, proton) check
    } else {
        v0AOD = ((AliAODEvent*)fEvt)->GetV0(v0Index);
        if (!v0AOD)
            return kFALSE;
        AliAODTrack* pTrackV0 =
            (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
        AliAODTrack* nTrackV0 =
            (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

        AliAODMCParticle* MCLamD1 =
            (AliAODMCParticle*)fMCArray->At(abs(pTrackV0->GetLabel()));
        AliAODMCParticle* MCLamD2 =
            (AliAODMCParticle*)fMCArray->At(abs(nTrackV0->GetLabel()));
        AliAODMCParticle* MCLam;
        AliAODMCParticle* MCPi;
        AliAODMCParticle* MCSigmaStar;

        if ((abs(MCLamD1->GetPdgCode()) == kProtonCode &&
             abs(MCLamD2->GetPdgCode()) == kPionCode) ||
            (abs(MCLamD1->GetPdgCode()) == kPionCode &&
             abs(MCLamD2->GetPdgCode()) == kProtonCode)) {
            if (MCLamD1->GetMother() == MCLamD2->GetMother()) {
                MCLam =
                    (AliAODMCParticle*)fMCArray->At(abs(MCLamD1->GetMother()));
                if (abs(MCLam->GetPdgCode()) == kLambdaCode) {
                    MCPi = (AliAODMCParticle*)fMCArray->At(
                        abs(track1->GetLabel()));
                    if (abs(MCPi->GetPdgCode()) == kPionCode) {
                        if (MCLam->GetMother() == MCPi->GetMother()) {
                            MCSigmaStar = (AliAODMCParticle*)fMCArray->At(
                                abs(MCLam->GetMother()));
                            if ((abs(MCSigmaStar->GetPdgCode()) ==
                                 kSigmaStarPCode) ||
                                (abs(MCSigmaStar->GetPdgCode()) ==
                                 kSigmaStarNCode)) {
                                if (IsPrimaryMC) {
                                    if (MCSigmaStar->IsPrimary()) {
                                        trueSigmaStar = kTRUE;
                                    }  // Primary(input) SigmaStar check
                                } else {
                                    trueSigmaStar = kTRUE;
                                }
                            }  // Sigma Star check
                        }      // Pion mother = Lambda mother
                    }          // pion check
                }              // Lambda check
            }                  // Lambda duather's mother check
        }                      // Lambda daughter (pion, proton) check
    }
    return trueSigmaStar;
}

THnSparse* AliAnalysisTaskSigma1385PM::CreateTHnSparse(
    TString name,
    TString title,
    Int_t ndim,
    std::vector<TAxis> bins,
    Option_t* opt) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    const TAxis* axises[bins.size()];
    for (UInt_t i = 0; i < bins.size(); i++)
        axises[i] = &bins[i];
    THnSparse* h = fHistos->CreateTHnSparse(name, title, ndim, axises, opt);
    return h;
}
Long64_t AliAnalysisTaskSigma1385PM::FillTHnSparse(TString name,
                                                     std::vector<Double_t> x,
                                                     Double_t w) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    auto hsparse = dynamic_cast<THnSparse*>(fHistos->FindObject(name));
    if (!hsparse) {
        std::cout << "ERROR : no " << name << std::endl;
        exit(1);
    }
    return FillTHnSparse(hsparse, x, w);
}

Long64_t AliAnalysisTaskSigma1385PM::FillTHnSparse(THnSparse* h,
                                                     std::vector<Double_t> x,
                                                     Double_t w) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    if (int(x.size()) != h->GetNdimensions()) {
        std::cout << "ERROR : wrong sized of array while Fill " << h->GetName()
                  << std::endl;
        exit(1);
    }
    return h->Fill(&x.front(), w);
}
TAxis AliAnalysisTaskSigma1385PM::AxisFix(TString name,
                                            int nbin,
                                            Double_t xmin,
                                            Double_t xmax) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(nbin, xmin, xmax);
    axis.SetName(name);
    return axis;
}
TAxis AliAnalysisTaskSigma1385PM::AxisStr(TString name,
                                            std::vector<TString> bin) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
    UInt_t i = 1;
    for (auto blabel : bin)
        ax.SetBinLabel(i++, blabel);
    return ax;
}

TAxis AliAnalysisTaskSigma1385PM::AxisVar(TString name,
                                            std::vector<Double_t> bin) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(bin.size() - 1, &bin.front());
    axis.SetName(name);
    return axis;
}
double AliAnalysisTaskSigma1385PM::GetTPCnSigma(AliVTrack* track,
                                                  AliPID::EParticleType type) {
    AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
    if (nanoT) {
        static bool used = false;
        if (!used) {
            AliNanoAODTrack::InitPIDIndex();
            used = true;
        }
        return nanoT->GetVar(
            AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, type));
    } else
        return fPIDResponse->NumberOfSigmasTPC(track, type);
}
void AliAnalysisTaskSigma1385PM::FillTrackToEventPool() {
    // Fill Selected tracks to event mixing pool
    if ((centbin < 0) || (zbin < 0))
        return;
    AliVTrack* goodtrack;

    tracklist* etl;
    eventpool* ep;
    // Event mixing pool

    ep = &fEMpool[centbin][zbin];
    ep->push_back(tracklist());
    etl = &(ep->back());
    // Fill selected tracks
    for (UInt_t i = 0; i < goodtrackindices.size(); i++) {
        goodtrack = (AliVTrack*)fEvt->GetTrack(goodtrackindices[i]);
        if (!goodtrack)
            continue;
        etl->push_back((AliVTrack*)goodtrack->Clone());
    }
    if (!goodtrackindices.size())
        ep->pop_back();
    if ((int)ep->size() > (int)fnMix) {
        for (auto it : ep->front())
            delete it;
        ep->pop_front();
    }
}