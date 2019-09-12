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

/* AliAnalysisTaskNanoCheck
 *
 *  Test code for the nano check
 *
 *  Author: Bong-Hwi Lim
 *          Maximiliano Puccio
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
#include "AliAnalysisTaskNanoCheck.h"
#include "AliESDv0.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "THistManager.h"

const Double_t pi = TMath::Pi();
const Double_t pionMass = AliPID::ParticleMass(AliPID::kPion);
const Double_t v0Mass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

enum {
    kBefore = 1, // 1
    kTrackCut,
    kEtaCut,
    kpTcut,
    kPIDcut,
    kPVDCAcut,
    krDCACut,
    kAll
};

class AliAnalysisTaskNanoCheck;

ClassImp(AliAnalysisTaskNanoCheck)  // classimp: necessary for root

    AliAnalysisTaskNanoCheck::AliAnalysisTaskNanoCheck()
    : AliAnalysisTaskSE(), fEvt(0) {}
//_____________________________________________________________________________
AliAnalysisTaskNanoCheck::AliAnalysisTaskNanoCheck(const char* name,
                                                           Bool_t MCcase)
    : AliAnalysisTaskSE(name), fEvt(0), IsMC(MCcase) {
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
AliAnalysisTaskNanoCheck::~AliAnalysisTaskNanoCheck() {}
//_____________________________________________________________________________
void AliAnalysisTaskNanoCheck::UserCreateOutputObjects() {
    // Standard Track cut for ESD
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

    fHistos = new THistManager("NanoCheckhists");
    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 100, "s");
    fHistos->CreateTH1("hNofTracks", "", 2, 0, 2, "s");
    fHistos->CreateTH1("hNofV0s", "", 2, 0, 2, "s");
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
    
    std::vector<double> centaxisbin = {
        0,  1,  5,  10, 15, 20, 30,
        40, 50, 60, 70, 80, 90, 100};  // can be use from pp to PbPb
    binCent = AxisVar("Cent", centaxisbin);
    auto binPt = AxisFix("Pt", 200, 0, 20);
    auto binMass = AxisFix("Mass", 2000, 1.0, 3.0);
    binZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10});
    auto hNumofTracks = AxisStr("Type", {"before", "trackCut", "etaCut", "pTcut", "nPIDcut",
                                    "dPVDCAcut", "rDCACut", "all"});

    fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());

    PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoCheck::UserExec(Option_t*) {
    AliVEvent* event = InputEvent();
    if (!event) {
        PostData(1, fHistos->GetListOfHistograms());
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
        return;  // event cut
    }

    fHistos->FillTH1("hMultiplicity", (double)fCent);

    if (fEvt->IsA() == AliAODEvent::Class())
        vertex = ((AliAODEvent*)fEvt)->GetPrimaryVertex();
    const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
    lPosPV[0] = pVtx->GetX();
    lPosPV[1] = pVtx->GetY();
    lPosPV[2] = pVtx->GetZ();

    bool checkPion = GoodTracksSelection();
    bool checkV0 = GoodV0Selection();

    PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoCheck::Terminate(Option_t*) {}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskNanoCheck::GoodTracksSelection() {
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

        if(IsNano)
            ((AliNanoAODTrack*)track)->GetImpactParameters(b[0], b[1]);
        else
            track->GetImpactParameters(b, bCov);

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
        
        if (TMath::Abs(fTPCNSigPion) > fTPCNsigNanoCheckerPionCut)
            continue;
        if (fEta > fNanoCheckerPionEtaCut)
            continue;
        if (pionPt < 0.15)
            continue;
        if (pionZ > fNanoCheckerPionZVertexCut)
            continue;
        if (pionDCA_r > pionSigmaDCA_r * fNanoCheckerPionXYVertexSigmaCut)
            continue;
        // ---------- Track selection done ----------
        fHistos->FillTH1("hNofTracks", 1.5);

        fHistos->FillTH1("QAcut/hDCAPVPion", pionZ);
        fHistos->FillTH1("QAcut/hDCArPVPion", pionDCA_r);
        fHistos->FillTH1("QAcut/hEtaPion", fEta);
        fHistos->FillTH2("QAcut/hTPCPIDPion", track->GetTPCmomentum(), track->GetTPCsignal());

        goodtrackindices.push_back(it);
    }
    return goodtrackindices.size();
}
Bool_t AliAnalysisTaskNanoCheck::GoodV0Selection() {
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
            fHistos->FillTH1("hNofV0s", 0.5);
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
                AcceptedV0 = kFALSE;

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
            fHistos->FillTH1("QA/hDecayLengthQA", fLength);
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
                fHistos->FillTH1("hNofV0s", 1.5);
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
            fHistos->FillTH1("hNofV0s", 0.5);
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
                fHistos->FillTH1("hNofV0s", 1.5);
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

THnSparse* AliAnalysisTaskNanoCheck::CreateTHnSparse(
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
Long64_t AliAnalysisTaskNanoCheck::FillTHnSparse(TString name,
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

Long64_t AliAnalysisTaskNanoCheck::FillTHnSparse(THnSparse* h,
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
TAxis AliAnalysisTaskNanoCheck::AxisFix(TString name,
                                            int nbin,
                                            Double_t xmin,
                                            Double_t xmax) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(nbin, xmin, xmax);
    axis.SetName(name);
    return axis;
}
TAxis AliAnalysisTaskNanoCheck::AxisStr(TString name,
                                            std::vector<TString> bin) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
    UInt_t i = 1;
    for (auto blabel : bin)
        ax.SetBinLabel(i++, blabel);
    return ax;
}

TAxis AliAnalysisTaskNanoCheck::AxisVar(TString name,
                                            std::vector<Double_t> bin) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(bin.size() - 1, &bin.front());
    axis.SetName(name);
    return axis;
}
double AliAnalysisTaskNanoCheck::GetTPCnSigma(AliVTrack* track,
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