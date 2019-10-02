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
 *  Output could be saved to nTuple by using SetFillnTuple(kTRUE) 
 *    -> can be used for TMVA input
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

ClassImp(AliAnalysisTaskSigma1385PM)

    AliAnalysisTaskSigma1385PM::AliAnalysisTaskSigma1385PM()
    : AliAnalysisTaskSE(), fEvt(0), fNtupleSigma1385(0) {}
//_____________________________________________________________________________
AliAnalysisTaskSigma1385PM::AliAnalysisTaskSigma1385PM(const char* name,
                                                           Bool_t MCcase)
    : AliAnalysisTaskSE(name), fEvt(0), IsMC(MCcase), fNtupleSigma1385(0) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TNtupleD::Class());
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
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

    fHistos = new THistManager("Sigma1385hists");
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

    CreateTHnSparse("Sigma1385", "Sigma1385", 4,
                    {binType, binCent, binPt, binMass}, "s");
    fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 100, "s");

    fHistos->CreateTH2("QA/hTPCPIDPion", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH1("QA/hEtaPion", "", 40, -2, 2);
    fHistos->CreateTH1("QA/hDCAPVPion", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCArPVPion", "", 300, 0, 3, "s");
    fHistos->CreateTH2("QA/hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDLambdaPion", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 2000, 0,
                       200);
    fHistos->CreateTH2("QA/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 2000, 0,
                       200);
    fHistos->CreateTH1("QA/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCAlambdaPV", "", 500, 0, 5, "s");
    fHistos->CreateTH1("QA/hCosPAlambda", "", 150, 0.85, 1.0, "s");
    fHistos->CreateTH1("QA/hTotMomQA", "", 100, 0, 10, "s");
    fHistos->CreateTH1("QA/hYlambda", "", 400, -2, 2, "s");
    fHistos->CreateTH1("QA/hDecayLengthQA", "", 100, 0, 100, "s");
    fHistos->CreateTH1("QA/hMassLambda", "", 500, 1.0, 1.5, "s");
    fHistos->CreateTH1("QA/hLambdaRxy", "", 400, 0, 400);

    fHistos->CreateTH2("QAcut/hTPCPIDPion", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH1("QAcut/hEtaPion", "", 40, -2, 2);
    fHistos->CreateTH1("QAcut/hDCAPVPion", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVPion", "", 300, 0, 3, "s");
    fHistos->CreateTH2("QAcut/hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDLambdaPion", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 2000, 0,
                       200);
    fHistos->CreateTH2("QAcut/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 2000, 0,
                       200);
    fHistos->CreateTH1("QAcut/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCAlambdaPV", "", 500, 0, 5, "s");
    fHistos->CreateTH1("QAcut/hCosPAlambda", "", 150, 0.85, 1.0, "s");
    fHistos->CreateTH1("QAcut/hTotMomQA", "", 100, 0, 10, "s");
    fHistos->CreateTH1("QAcut/hYlambda", "", 400, -2, 2, "s");
    fHistos->CreateTH1("QAcut/hDecayLengthQA", "", 100, 0, 100, "s");
    fHistos->CreateTH1("QAcut/hMassLambda", "", 500, 1.0, 1.5, "s");
    fHistos->CreateTH1("QAcut/hLambdaRxy", "", 400, 0, 400);

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
        if(!IsNano)
            IsNano = kTRUE;
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
    lPosPV[0] = pVtx->GetX();
    lPosPV[1] = pVtx->GetY();
    lPosPV[2] = pVtx->GetZ();

    // Event Mixing pool -----------------------------------------------------
    zbin = binZ.FindBin(lPosPV[2]) - 1;           // Event mixing z-bin
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
    Float_t b[2];
    Float_t bCov[3];
    Double_t fTPCNSigPion, pionZ, pionPt, pionSigmaDCA_r, pionDCA_r, fEta;

    for (UInt_t it = 0; it < nTracks; it++) {
        track = (AliVTrack*)fEvt->GetTrack(it);
        if (!track)
            continue;
        fHistos->FillTH1("hNofTracks", 0.5);

        GetImpactParam(track, b, bCov);

        // ---------- Track selection begin ----------
        if (fEvt->IsA() == AliESDEvent::Class()) {
            if (!fTrackCuts->AcceptTrack((AliESDtrack*)track))
                continue;
        }  // ESD Case
        else {
            if (!IsNano && !((AliAODTrack*)track)->TestFilterBit(32))
                continue;
        }  // AOD Case

        pionZ = b[1];
        fTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
        pionPt = track->Pt();
        pionSigmaDCA_r = (0.0026 + 0.0050 / pionPt);
        pionDCA_r = b[0];
        fEta = TMath::Abs(track->Eta());

        fHistos->FillTH1("QA/hDCAPVPion", pionZ);
        fHistos->FillTH1("QA/hDCArPVPion", pionDCA_r);
        fHistos->FillTH1("QA/hEtaPion", fEta);
        fHistos->FillTH2("QA/hTPCPIDPion", track->GetTPCmomentum(), track->GetTPCsignal());

        if (TMath::Abs(fTPCNSigPion) > fTPCNsigSigmaStarPionCut)
            continue;
        if (fEta > fSigmaStarPionEtaCut)
            continue;
        if (pionPt < 0.15)
            continue;
        if (pionZ > fSigmaStarPionZVertexCut)
            continue;
        if (pionDCA_r > pionSigmaDCA_r * fSigmaStarPionXYVertexSigmaCut)
            continue;
        
        fHistos->FillTH1("QAcut/hDCAPVPion", pionZ);
        fHistos->FillTH1("QAcut/hDCArPVPion", pionDCA_r);
        fHistos->FillTH1("QAcut/hEtaPion", fEta);
        fHistos->FillTH2("QAcut/hTPCPIDPion", track->GetTPCmomentum(), track->GetTPCsignal());

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
    Bool_t fPIDLambda, fPIDAntiLambda;
    Double_t fDCADist_LambdaProton_PV, fDCADist_LambdaPion_PV;
    Double_t v0Position[3];
    Double_t fMassV0;
    Double_t radius, lV0TotalMomentum, fLength, lLifetime;
    UInt_t isAnti = 0;

    Bool_t AcceptedV0 = kTRUE;
    if (fEvt->IsA() == AliESDEvent::Class()) {  // ESD case
        for (UInt_t it = 0; it < nV0; it++) {
            fPIDLambda = kFALSE;
            fPIDAntiLambda = kFALSE;
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
            
            // filter like-sign V0
            if (pTrackV0->GetSign()*nTrackV0->GetSign() > 0)
                AcceptedV0 = kFALSE;

            // PID cuts
            Double_t fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            Double_t fTPCNSigAntiProton =
                GetTPCnSigma(nTrackV0, AliPID::kProton);
            Double_t fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            Double_t fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            if ((TMath::Abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) 
                && (TMath::Abs(fTPCNSigPion) < fTPCNsigLambdaPionCut) 
                && (nTrackV0->GetSign() < 0) ){
                
                fPIDLambda = kTRUE;
                v0ESD->ChangeMassHypothesis(kLambda0);
                isAnti = 0;

                fHistos->FillTH2("QA/hTPCPIDLambdaProton",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA/hTPCPIDLambdaPion",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
            }
            else if ((TMath::Abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) 
                    && (TMath::Abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut) 
                    && (nTrackV0->GetSign() > 0)){

                v0ESD->ChangeMassHypothesis(kLambda0Bar);
                fPIDAntiLambda = kTRUE;
                isAnti = 1;

                fHistos->FillTH2("QA/hTPCPIDAntiLambdaProton",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA/hTPCPIDAntiLambdaPion",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
            }
            else
                AcceptedV0 = kFALSE;

            // DCA cut
            // DCA between Dautgher particles
            Double_t fDCADistLambda = TMath::Abs(v0ESD->GetDcaV0Daughters());
            fHistos->FillTH1("QA/hDCA_lambdaDaughters", fDCADistLambda);
            if (fDCADistLambda > fDCADistLambdaDaughtersCut)
                AcceptedV0 = kFALSE;  // DCA proton-pion

            // DCA to PV
            Double_t fDCADistLambda_PV = TMath::Abs(v0ESD->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));
            fHistos->FillTH1("QA/hDCAlambdaPV", fDCADistLambda_PV);
            if (fDCADistLambda_PV > fDCArDistLambdaPVCut)
                AcceptedV0 = kFALSE;

            // CPA cut
            Double_t fLambdaCPA =
                TMath::Abs(v0ESD->GetV0CosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]));
            fHistos->FillTH1("QA/hCosPAlambda", fLambdaCPA);
            if (fLambdaCPA < fV0CosineOfPointingAngleCut || fLambdaCPA >= 1)
                AcceptedV0 = kFALSE;
            
            // Rapidity cut
            fHistos->FillTH1("QA/hYlambda", v0ESD->RapLambda());
            if (TMath::Abs(v0ESD->RapLambda()) > fMaxLambdaRapidity)
                return kFALSE;

            // Radius cut
            v0ESD->GetXYZ(v0Position[0], v0Position[1], v0Position[2]);
            radius = TMath::Hypot(v0Position[0], v0Position[1]);
            fHistos->FillTH1("QA/hLambdaRxy", radius);
            if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
                AcceptedV0 = kFALSE;

            // Life time cut
            lV0TotalMomentum = v0ESD->P();
            fHistos->FillTH1("QA/hTotMomQA", lV0TotalMomentum);
            fLength = TMath::Sqrt(TMath::Power(v0Position[0] - lPosPV[0], 2) +
                                  TMath::Power(v0Position[1] - lPosPV[1], 2) +
                                  TMath::Power(v0Position[2] - lPosPV[2], 2));
            
            lLifetime = TMath::Abs(v0Mass * fLength / lV0TotalMomentum);
            if (lLifetime > fLambdaLifetime)
                AcceptedV0 = kFALSE;

            // Mass window cut
            fMassV0 = v0ESD->GetEffMass();
            fHistos->FillTH1("QA/hMassLambda", fMassV0);
            if (TMath::Abs(fMassV0 - v0Mass) > fV0MassWindowCut)
                AcceptedV0 = kFALSE;

            // After selection above
            if (AcceptedV0){
                goodv0indices.push_back({it, isAnti});  // for standard V0
                if (fPIDLambda) {
                    fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                }
                if (fPIDAntiLambda) {
                    fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaProton",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaPion",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                }
                fHistos->FillTH1("QAcut/hDCA_lambdaDaughters", fDCADistLambda);
                fHistos->FillTH1("QAcut/hDCAlambdaPV", fDCADistLambda_PV);
                fHistos->FillTH1("QAcut/hCosPAlambda", fLambdaCPA);
                fHistos->FillTH1("QAcut/hYlambda", v0ESD->RapLambda());
                fHistos->FillTH1("QAcut/hLambdaRxy", radius);
                fHistos->FillTH1("QAcut/hTotMomQA", lV0TotalMomentum);
                fHistos->FillTH1("QAcut/hMassLambda", fMassV0);
                fHistos->FillTH1("QAcut/hDecayLengthQA", fLength);
            }
        }                                     // All V0 loop
    }                                         // ESD case
    else {
        for (UInt_t it = 0; it < nV0; it++) {
            fPIDLambda = kFALSE;
            fPIDAntiLambda = kFALSE;
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

            // filter like-sign V0
            if (pTrackV0->GetSign()*nTrackV0->GetSign() > 0)
                AcceptedV0 = kFALSE;

            // PID cuts
            Double_t fTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            Double_t fTPCNSigAntiProton =
                GetTPCnSigma(nTrackV0, AliPID::kProton);
            Double_t fTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            Double_t fTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            if ((TMath::Abs(fTPCNSigProton) < fTPCNsigLambdaProtonCut) 
                && (TMath::Abs(fTPCNSigPion) < fTPCNsigLambdaPionCut) 
                && (nTrackV0->GetSign() < 0) ){
                
                fPIDLambda = kTRUE;
                isAnti = 0;

                fHistos->FillTH2("QA/hTPCPIDLambdaProton",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA/hTPCPIDLambdaPion",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
            }
            else if ((TMath::Abs(fTPCNSigAntiProton) < fTPCNsigLambdaProtonCut) 
                    && (TMath::Abs(fTPCNSigAntiPion) < fTPCNsigLambdaPionCut) 
                    && (nTrackV0->GetSign() > 0)){

                fPIDAntiLambda = kTRUE;
                isAnti = 1;

                fHistos->FillTH2("QA/hTPCPIDAntiLambdaProton",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA/hTPCPIDAntiLambdaPion",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
            }
            else
                AcceptedV0 = kFALSE;

            // DCA cut
            Double_t fDCADistLambda = TMath::Abs(v0AOD->DcaV0Daughters());
            fHistos->FillTH1("QA/hDCA_lambdaDaughters", fDCADistLambda);
            if (fDCADistLambda > fDCADistLambdaDaughtersCut)
                AcceptedV0 = kFALSE;  // DCA proton-pion

            // DCA to PV
            Double_t fDCADistLambda_PV = TMath::Abs(v0AOD->DcaV0ToPrimVertex());
            fHistos->FillTH1("QA/hDCAlambdaPV", fDCADistLambda_PV);
            if (fDCADistLambda_PV > fDCArDistLambdaPVCut)
                AcceptedV0 = kFALSE;

            // CPA cut
            Double_t fLambdaCPA = TMath::Abs(v0AOD->CosPointingAngle(vertex));
            fHistos->FillTH1("QA/hCosPAlambda", fLambdaCPA);
            if (fLambdaCPA < fV0CosineOfPointingAngleCut || fLambdaCPA >= 1)
                AcceptedV0 = kFALSE;

            // Rapidity cut
            fHistos->FillTH1("QA/hYlambda", v0AOD->RapLambda());
            if (TMath::Abs(v0AOD->RapLambda()) > fMaxLambdaRapidity)
                return kFALSE;

            // Radius cut
            radius = v0AOD->RadiusV0();
            fHistos->FillTH1("QA/hLambdaRxy", radius);
            if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
                AcceptedV0 = kFALSE;

            // Life time cut
            lV0TotalMomentum = TMath::Sqrt(v0AOD->Ptot2V0());
            fHistos->FillTH1("QA/hTotMomQA", lV0TotalMomentum);
            fLength = v0AOD->DecayLength(lPosPV);
            fHistos->FillTH1("QA/hDecayLengthQA", fLength);
            lLifetime = TMath::Abs(v0Mass * fLength / lV0TotalMomentum);
            if (lLifetime > fLambdaLifetime)
                AcceptedV0 = kFALSE;

            // Mass window cut
            fMassV0 = -999;
            if (fPIDLambda)
                fMassV0 = v0AOD->MassLambda();
            else if (fPIDAntiLambda)
                fMassV0 = v0AOD->MassAntiLambda();
            fHistos->FillTH1("QA/hMassLambda", fMassV0);
            if (TMath::Abs(fMassV0 - v0Mass) > fV0MassWindowCut)
                AcceptedV0 = kFALSE;

            // After selection above
            if (AcceptedV0) {
                goodv0indices.push_back({it, isAnti});
                if (fPIDLambda) {
                    fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                }
                if (fPIDAntiLambda) {
                    fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaProton",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaPion",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                }
                fHistos->FillTH1("QAcut/hDCA_lambdaDaughters", fDCADistLambda);
                fHistos->FillTH1("QAcut/hDCAlambdaPV", fDCADistLambda_PV);
                fHistos->FillTH1("QAcut/hCosPAlambda", fLambdaCPA);
                fHistos->FillTH1("QAcut/hYlambda", v0AOD->RapLambda());
                fHistos->FillTH1("QAcut/hLambdaRxy", radius);
                fHistos->FillTH1("QAcut/hTotMomQA", lV0TotalMomentum);
                fHistos->FillTH1("QAcut/hMassLambda", fMassV0);
                fHistos->FillTH1("QAcut/hDecayLengthQA", fLength);
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
    Bool_t isAnti, isPionPlus;
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
            v0ESD = ((AliESDEvent*)fEvt)->GetV0(goodv0indices[i][0]);
            if (!v0ESD)
                continue;
            AliESDtrack* pTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
            AliESDtrack* nTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();

            if (goodv0indices[i][1] > 0)
                isAnti = true;
            else
                isAnti = false;

            if (!isAnti)
                v0ESD->ChangeMassHypothesis(kLambda0);
            else
                v0ESD->ChangeMassHypothesis(kLambda0Bar);

            temp1.SetXYZM(v0ESD->Px(), v0ESD->Py(), v0ESD->Pz(),
                          v0ESD->GetEffMass());
        } else {
            v0AOD = ((AliAODEvent*)fEvt)->GetV0(goodv0indices[i][0]);
            if (!v0AOD)
                continue;
            AliAODTrack* pTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack* nTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();

            if (goodv0indices[i][1] > 0)
                isAnti = true;
            else
                isAnti = false;

            if (!isAnti)
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassLambda());
            else
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

            if (track1->Charge() > 0)
                isPionPlus = true;
            else
                isPionPlus = false;

            if (!isAnti && isPionPlus)
                sign = kSigmaStarP;
            if (!isAnti && !isPionPlus)
                sign = kSigmaStarN;
            if (isAnti && isPionPlus)
                sign = kAntiSigmaStarN;
            if (isAnti && !isPionPlus)
                sign = kAntiSigmaStarP;

            FillTHnSparse("Sigma1385", {(double)sign, (double)fCent,
                                        vecsum.Pt(), vecsum.M()});

            if (IsMC &&
                IsTrueSigmaStar(goodv0indices[i][0], goodtrackindices[j])) {
                if (!isAnti && isPionPlus)
                    sign = kSigmaStarP_REC;
                if (!isAnti && !isPionPlus)
                    sign = kSigmaStarN_REC;
                if (isAnti && isPionPlus)
                    sign = kAntiSigmaStarN_REC;
                if (isAnti && !isPionPlus)
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

                if (goodv0indices[i][1] > 0)
                    isAnti = true;
                else
                    isAnti = false;

                if (!isAnti && isPionPlus)
                    sign = kSigmaStarP_MIX;
                if (!isAnti && !isPionPlus)
                    sign = kSigmaStarN_MIX;
                if (isAnti && isPionPlus)
                    sign = kAntiSigmaStarN_MIX;
                if (isAnti && !isPionPlus)
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
    Bool_t isAnti, isPionPlus;
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
            v0ESD = ((AliESDEvent*)fEvt)->GetV0(goodv0indices[i][0]);
            if (!v0ESD)
                continue;
            AliESDtrack* pTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
            AliESDtrack* nTrackV0 =
                ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();

            if (goodv0indices[i][1] > 0)
                isAnti = true;
            else
                isAnti = false;

            if (!isAnti)
                v0ESD->ChangeMassHypothesis(kLambda0);
            else
                v0ESD->ChangeMassHypothesis(kLambda0Bar);

            temp1.SetXYZM(v0ESD->Px(), v0ESD->Py(), v0ESD->Pz(),
                          v0ESD->GetEffMass());

            // nTuple
            tmp[5] = fTPCNSigProton;                     // PIDV0pTrackProton
            tmp[6] = fTPCNSigAntiProton;                 // PIDV0pTrackPion
            tmp[7] = fTPCNSigPion;                       // PIDV0nTrackProton
            tmp[8] = fTPCNSigAntiPion;                   // PIDV0nTrackPion
            tmp[9] = TMath::Abs(v0ESD->GetDcaV0Daughters());   // DCAV0Daughters
            tmp[10] = TMath::Abs(v0ESD->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));  // DCAV0ToPrimVertex
            tmp[11] = v0ESD->GetV0CosineOfPointingAngle(
                lPosPV[0], lPosPV[1], lPosPV[2]);             // CosPointingAngleV0
            tmp[12] = v0ESD->GetEffMass();  // V0Mass
            tmp[13] = v0ESD->Eta();         // EtaV0
            tmp[14] = v0ESD->Phi();         // PhiV0
        } else {
            v0AOD = ((AliAODEvent*)fEvt)->GetV0(goodv0indices[i][0]);
            if (!v0AOD)
                continue;
            AliAODTrack* pTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack* nTrackV0 =
                (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
            pID = pTrackV0->GetID();
            nID = nTrackV0->GetID();

            if (goodv0indices[i][1] > 0)
                isAnti = true;
            else
                isAnti = false;

            if (!isAnti)
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassLambda());
            else
                temp1.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                              v0AOD->MassAntiLambda());
            // nTuple
            tmp[5] = fTPCNSigProton;                     // PIDV0pTrackProton
            tmp[6] = fTPCNSigAntiProton;                 // PIDV0pTrackPion
            tmp[7] = fTPCNSigPion;                       // PIDV0nTrackProton
            tmp[8] = fTPCNSigAntiPion;                   // PIDV0nTrackPion
            tmp[9] = TMath::Abs(v0AOD->DcaV0Daughters());      // DCAV0Daughters
            tmp[10] = TMath::Abs(v0AOD->DcaV0ToPrimVertex());  // DCAV0ToPrimVertex
            tmp[11] = v0AOD->CosPointingAngle(vertex);   // CosPointingAngleV0
            if (!isAnti)
                tmp[12] = v0AOD->MassLambda();  // V0Mass
            else
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

            if (goodv0indices[i][1] > 0)
                isAnti = true;
            else
                isAnti = false;

            if (!isAnti && isPionPlus)
                sign = kSigmaStarP;
            if (!isAnti && !isPionPlus)
                sign = kSigmaStarN;
            if (isAnti && isPionPlus)
                sign = kAntiSigmaStarN;
            if (isAnti && !isPionPlus)
                sign = kAntiSigmaStarP;

            tmp[0] = GetTPCnSigma(track1, AliPID::kPion);  // PIDSigmaStarPion
            tmp[1] = TMath::Abs(track1->GetZ() - lPosPV[2]);  // DCASigmaStarPionToPrimVertexZ
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
            // tmp[9] = TMath::Abs(v0ESD->GetDcaV0Daughters()); //DCAV0Daughters
            // tmp[10] = TMath::Abs(v0ESD->GetD(lPosPV[0], lPosPV[1], lPosPV[2])); //DCAV0ToPrimVertex
            // tmp[11] = v0ESD->GetV0CosineOfPointingAngle(lPosPV[0], lPosPV[1],
            // lPosPV[2]);//CosPointingAngleV0 tmp[12] = v0ESD->GetEffMass(); //V0Mass
            // tmp[13] = v0ESD->Eta(); //EtaV0
            // tmp[14] = v0ESD->Phi(); //PhiV0

            if (IsMC) {
                if (IsTrueSigmaStar(goodv0indices[i][0], goodtrackindices[j]))
                    tmp[15] = (int)sign;  // MCflag
                else
                    tmp[15] = 5;  // MCflag -> not true
            } else
                tmp[15] = 0;  // MCflag -> data
            fNtupleSigma1385->Fill(tmp);
        }                     // pion loop
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
            if ((TMath::Abs(v0PdgCode) != kSigmaStarPCode) &&
                (TMath::Abs(v0PdgCode) != kSigmaStarNCode))
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

            if ((TMath::Abs(v0PdgCode) != kSigmaStarPCode) &&
                (TMath::Abs(v0PdgCode) != kSigmaStarNCode))
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
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(pTrackV0->GetLabel()))
                ->Particle();
        TParticle* MCLamD2 =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(nTrackV0->GetLabel()))
                ->Particle();
        TParticle* MCLam;
        TParticle* MCPi;
        TParticle* MCSigmaStar;
        if ((TMath::Abs(MCLamD1->GetPdgCode()) == kProtonCode &&
             TMath::Abs(MCLamD2->GetPdgCode()) == kPionCode) ||
            (TMath::Abs(MCLamD1->GetPdgCode()) == kPionCode &&
             TMath::Abs(MCLamD2->GetPdgCode()) == kProtonCode)) {
            if (MCLamD1->GetMother(0) == MCLamD2->GetMother(0)) {
                MCLam =
                    (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLamD1->GetMother(0)))
                        ->Particle();
                if (TMath::Abs(MCLam->GetPdgCode()) == kLambdaCode) {
                    MCPi =
                        (TParticle*)fMCEvent->GetTrack(TMath::Abs(track1->GetLabel()))
                            ->Particle();
                    if (TMath::Abs(MCPi->GetPdgCode()) == kPionCode) {
                        if (MCPi->GetMother(0) == MCLam->GetMother(0)) {
                            MCSigmaStar =
                                (TParticle*)fMCEvent
                                    ->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                                    ->Particle();
                            if ((TMath::Abs(MCSigmaStar->GetPdgCode()) ==
                                 kSigmaStarPCode) ||
                                (TMath::Abs(MCSigmaStar->GetPdgCode()) ==
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
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(pTrackV0->GetLabel()));
        AliAODMCParticle* MCLamD2 =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(nTrackV0->GetLabel()));
        AliAODMCParticle* MCLam;
        AliAODMCParticle* MCPi;
        AliAODMCParticle* MCSigmaStar;

        if ((TMath::Abs(MCLamD1->GetPdgCode()) == kProtonCode &&
             TMath::Abs(MCLamD2->GetPdgCode()) == kPionCode) ||
            (TMath::Abs(MCLamD1->GetPdgCode()) == kPionCode &&
             TMath::Abs(MCLamD2->GetPdgCode()) == kProtonCode)) {
            if (MCLamD1->GetMother() == MCLamD2->GetMother()) {
                MCLam =
                    (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLamD1->GetMother()));
                if (TMath::Abs(MCLam->GetPdgCode()) == kLambdaCode) {
                    MCPi = (AliAODMCParticle*)fMCArray->At(
                        TMath::Abs(track1->GetLabel()));
                    if (TMath::Abs(MCPi->GetPdgCode()) == kPionCode) {
                        if (MCLam->GetMother() == MCPi->GetMother()) {
                            MCSigmaStar = (AliAODMCParticle*)fMCArray->At(
                                TMath::Abs(MCLam->GetMother()));
                            if ((TMath::Abs(MCSigmaStar->GetPdgCode()) ==
                                 kSigmaStarPCode) ||
                                (TMath::Abs(MCSigmaStar->GetPdgCode()) ==
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
void AliAnalysisTaskSigma1385PM::GetImpactParam(AliVTrack* track,
                                                Float_t p[2],
                                                Float_t cov[3]) {
    AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
    if (nanoT)
        nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
    else
        track->GetImpactParameters(p, cov);
}