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
#include "AliMultSelectionTask.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAODcascade.h"
#include "AliAODv0.h"
#include "AliAnalysisTaskNanoCheck.h"
#include "AliESDcascade.h"
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
    fHistos->CreateTH1("hMultiplicity", "", 101, -1, 100, "s");
    if (checkTracks)
        fHistos->CreateTH1("hNofTracks", "", 2, 0, 2, "s");
    if (checkV0s) {
        fHistos->CreateTH1("hNofV0s", "", 2, 0, 2, "s");
        fHistos->CreateTH1("hNofV0sDistribution", "", 500, 0, 500, "s");
    }
    if (checkCascades) {
        fHistos->CreateTH1("hNofCascades", "", 2, 0, 2, "s");
        fHistos->CreateTH1("hNofCascadesDistribution", "", 50, 0, 50, "s");
    }
    if (checkTracks) {
        fHistos->CreateTH2("QA_pion/hTPCPIDPion", "", 200, 0, 20, 2000, 0, 200);
        fHistos->CreateTH1("QA_pion/hEtaPion", "", 40, -2, 2);
        fHistos->CreateTH1("QA_pion/hDCAPVPion", "", 300, 0, 3, "s");
        fHistos->CreateTH1("QA_pion/hDCArPVPion", "", 300, 0, 3, "s");
    }
    if (checkV0s) {
        fHistos->CreateTH2("QA_V0/hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0, 200);
        fHistos->CreateTH2("QA_V0/hTPCPIDLambdaPion", "", 200, 0, 20, 2000, 0, 200);
        fHistos->CreateTH2("QA_V0/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 2000, 0,
                        200);
        fHistos->CreateTH2("QA_V0/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 2000, 0,
                        200);
        fHistos->CreateTH1("QA_V0/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
        fHistos->CreateTH1("QA_V0/hDCAlambdaPV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_V0/hDCAlambdaPionPV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_V0/hDCAlambdaProtonPV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_V0/hCosPAlambda", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QA_V0/hTotMomQA", "", 100, 0, 10, "s");
        fHistos->CreateTH1("QA_V0/hYlambda", "", 400, -2, 2, "s");
        fHistos->CreateTH1("QA_V0/hDecayLengthQA", "", 100, 0, 100, "s");
        fHistos->CreateTH1("QA_V0/hMassLambda", "", 500, 1.0, 1.5, "s");
        fHistos->CreateTH1("QA_V0/hLambdaRxy", "", 400, 0, 400);
    }
    if (checkCascades) {
        fHistos->CreateTH2("QA_Xi/hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0,
                           200);
        fHistos->CreateTH2("QA_Xi/hTPCPIDLambdaPion", "", 200, 0, 20, 2000, 0,
                           200);
        fHistos->CreateTH2("QA_Xi/hTPCPIDBachelorPion", "", 200, 0, 20, 2000, 0,
                           200);
        fHistos->CreateTH1("QA_Xi/hTPCPIDsignalLambdaProton", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QA_Xi/hTPCPIDsignalLambdaPion", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QA_Xi/hTPCPIDsignalBachelorPion", "", 100, -5, 5,
                           "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_Lambda_BTW_Daughters", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_Xi_BTW_Daughters", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_lambda_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_Xi_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_LambdaProton_to_PV", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_LambdaPion_to_PV", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("QA_Xi/hDCADist_BachelorPion_to_PV", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("QA_Xi/hCosPA_lambda", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QA_Xi/hCosPA_Xi", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QA_Xi/hMass_Xi", "", 500, 1.0, 1.5, "s");
        fHistos->CreateTH2("QA_Xi/hPhiEta_Xi", "", 180, 0, 2 * pi, 40, -2, 2);
        fHistos->CreateTH2("QA_Xi/hLambda_Rxy", "", 400, -200, 200, 400, -200,
                           200);
        fHistos->CreateTH2("QA_Xi/hXi_Rxy", "", 400, -200, 200, 400, -200,
                           200);
    }

    if (checkTracks) {
        fHistos->CreateTH2("QA_pionCut/hTPCPIDPion", "", 200, 0, 20, 2000, 0, 200);
        fHistos->CreateTH1("QA_pionCut/hEtaPion", "", 40, -2, 2);
        fHistos->CreateTH1("QA_pionCut/hDCAPVPion", "", 300, 0, 3, "s");
        fHistos->CreateTH1("QA_pionCut/hDCArPVPion", "", 300, 0, 3, "s");
    }
    if (checkV0s) {
        fHistos->CreateTH2("QA_V0Cut/hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0, 200);
        fHistos->CreateTH2("QA_V0Cut/hTPCPIDLambdaPion", "", 200, 0, 20, 2000, 0, 200);
        fHistos->CreateTH2("QA_V0Cut/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 2000, 0,
                        200);
        fHistos->CreateTH2("QA_V0Cut/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 2000, 0,
                        200);
        fHistos->CreateTH1("QA_V0Cut/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
        fHistos->CreateTH1("QA_V0Cut/hDCAlambdaPV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_V0Cut/hDCAlambdaPionPV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_V0Cut/hDCAlambdaProtonPV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_V0Cut/hCosPAlambda", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QA_V0Cut/hTotMomQA", "", 100, 0, 10, "s");
        fHistos->CreateTH1("QA_V0Cut/hYlambda", "", 400, -2, 2, "s");
        fHistos->CreateTH1("QA_V0Cut/hDecayLengthQA", "", 100, 0, 100, "s");
        fHistos->CreateTH1("QA_V0Cut/hMassLambda", "", 500, 1.0, 1.5, "s");
        fHistos->CreateTH1("QA_V0Cut/hLambdaRxy", "", 400, 0, 400);
    }
    if (checkCascades) {
        fHistos->CreateTH2("QA_XiCut/hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0,
                           200);
        fHistos->CreateTH2("QA_XiCut/hTPCPIDLambdaPion", "", 200, 0, 20, 2000,
                           0, 200);
        fHistos->CreateTH2("QA_XiCut/hTPCPIDBachelorPion", "", 200, 0, 20, 2000,
                           0, 200);
        fHistos->CreateTH1("QA_XiCut/hTPCPIDsignalLambdaProton", "", 100, -5, 5,
                           "s");
        fHistos->CreateTH1("QA_XiCut/hTPCPIDsignalLambdaPion", "", 100, -5, 5,
                           "s");
        fHistos->CreateTH1("QA_XiCut/hTPCPIDsignalBachelorPion", "", 100, -5, 5,
                           "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_Lambda_BTW_Daughters", "", 300, 0,
                           3, "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_Xi_BTW_Daughters", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_lambda_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_Xi_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_LambdaProton_to_PV", "", 500, 0,
                           5, "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_LambdaPion_to_PV", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("QA_XiCut/hDCADist_BachelorPion_to_PV", "", 500, 0,
                           5, "s");
        fHistos->CreateTH1("QA_XiCut/hCosPA_lambda", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QA_XiCut/hCosPA_Xi", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QA_XiCut/hMass_Xi", "", 500, 1.0, 1.5, "s");
        fHistos->CreateTH2("QA_XiCut/hPhiEta_Xi", "", 180, 0, 2 * pi, 40, -2, 2);
        fHistos->CreateTH2("QA_XiCut/hLambda_Rxy", "", 400, -200, 200, 400,
                           -200, 200);
        fHistos->CreateTH2("QA_XiCut/hXi_Rxy", "", 400, -200, 200, 400, -200,
                           200);
    }
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
    bool IsEvtSelected{false}, IsINEL0True{false};

    // Input Events
    AliVEvent* event = InputEvent();
    if (!event) {
        PostData(1, fHistos->GetListOfHistograms());
        AliInfo("Could not retrieve event");
        return;
    }
    // Nano?
    AliNanoAODHeader* nanoHeader =
        dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());

    // AOD? ESD?
    event->IsA() == AliESDEvent::Class()
        ? fEvt = dynamic_cast<AliESDEvent*>(event)
        : fEvt = dynamic_cast<AliAODEvent*>(event);
    if (!fEvt) {
        PostData(1, fHistos->GetListOfHistograms());
        return;
    }
    if (!IsAOD && (event->IsA() != AliESDEvent::Class()))
        IsAOD = true;

    // Input hander for PID
    AliInputEventHandler* inputHandler =
        (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();

    if (!nanoHeader) { // ESD and AOD case
        IsEvtSelected = fEventCuts.AcceptEvent(event); // Initialize AliEventCuts
        if (IsMC) { // MC
            if (fEvt->IsA() != AliESDEvent::Class())
                fMCArray = (TClonesArray*)fEvt->FindListObject(
                    "mcparticles");  // AOD Case
            fMCEvent = MCEvent();
            IsINEL0True = fEventCuts.IsTrueINELgtZero(fEvt, true);
        }
        // Get default values from AliEventCuts
        fCent = AliMultSelectionTask::IsINELgtZERO(event) 
                    ? fEventCuts.GetCentrality() 
                    : -0.5;
        // PID response
        fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
        if (!fPIDResponse)
            AliInfo("No PIDd");
    } else { // Nano case
        if (!IsNano)
            IsNano = kTRUE;
        IsEvtSelected = true; // all events are already filtered
        // Get default values from NanoHeader
        fCent = nanoHeader->GetCentr("V0M");
        static int inel_index = -1;
        if (inel_index < 0) inel_index = nanoHeader->GetVarIndex("cstINELgt0");
        if ((inel_index > 0) && (nanoHeader->GetVar(inel_index) < 0.5))
            fCent = -0.5;
    }

    if (!IsEvtSelected) {
        PostData(1, fHistos->GetListOfHistograms());
        return;  // event cut
    }
    // Event cuts passed!
    // Basic QA plots
    fHistos->FillTH1("hMultiplicity", (double)fCent);

    // Get Vertex informations
    if (IsAOD)
        vertex = ((AliAODEvent*)fEvt)->GetPrimaryVertex();
    const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
    lPosPV[0] = pVtx->GetX();
    lPosPV[1] = pVtx->GetY();
    lPosPV[2] = pVtx->GetZ();

    // Do the jobs
    bool checkPion{false}, checkV0{false}, checkCascade{false};
    if (checkTracks)
        checkPion = GoodTracksSelection();
    if (checkV0s)
        checkV0 = GoodV0Selection();
    if (checkCascades)
        checkCascade = GoodCascadeSelection();

    PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoCheck::Terminate(Option_t*) { }
//_____________________________________________________________________________
Bool_t AliAnalysisTaskNanoCheck::GoodTracksSelection() {
    const UInt_t nTracks = fEvt->GetNumberOfTracks();
    goodtrackindices.clear();
    AliVTrack* track;
    Float_t b[2];
    Float_t bCov[3];
    Double_t TPCNSigPion, pionZ, pionPt, pionSigmaDCAr, pionDCAr, pionEta;

    for (UInt_t it = 0; it < nTracks; it++) {
        track = (AliVTrack*)fEvt->GetTrack(it);
        if (!track)
            continue;
        fHistos->FillTH1("hNofTracks", 0.5);

        GetImpactParam(track, b, bCov);

        // ---------- Track selection begin ----------
        if (!IsAOD) {
            if (!fTrackCuts->AcceptTrack((AliESDtrack*)track))
                continue;
        }  // ESD Case
        else {
            if (!IsNano) {
                if (!((AliAODTrack*)track)->TestFilterBit(32))
                    continue;
            } else {
                if (!(static_cast<AliNanoAODTrack*>(track)->TestFilterBit(32)))
                    continue;
            }
        }

        pionZ = b[1];
        TPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
        pionPt = track->Pt();
        pionSigmaDCAr = (0.0026 + 0.0050 / pionPt);
        pionDCAr = b[0];
        pionEta = TMath::Abs(track->Eta());

        fHistos->FillTH1("QA_pion/hDCAPVPion", pionZ);
        fHistos->FillTH1("QA_pion/hDCArPVPion", pionDCAr);
        fHistos->FillTH1("QA_pion/hEtaPion", pionEta);
        fHistos->FillTH2("QA_pion/hTPCPIDPion", track->GetTPCmomentum(),
                         track->GetTPCsignal());

        if (TMath::Abs(TPCNSigPion) > fTPCNsigNanoCheckerPionCut)
            continue;
        if (pionEta > fNanoCheckerPionEtaCut)
            continue;
        if (pionPt < 0.15)
            continue;
        if (pionZ > fNanoCheckerPionZVertexCut)
            continue;
        if (pionDCAr > pionSigmaDCAr * fNanoCheckerPionXYVertexSigmaCut)
            continue;
        // ---------- Track selection done ----------
        fHistos->FillTH1("hNofTracks", 1.5);

        fHistos->FillTH1("QA_pionCut/hDCAPVPion", pionZ);
        fHistos->FillTH1("QA_pionCut/hDCArPVPion", pionDCAr);
        fHistos->FillTH1("QA_pionCut/hEtaPion", pionEta);
        fHistos->FillTH2("QA_pionCut/hTPCPIDPion", track->GetTPCmomentum(),
                         track->GetTPCsignal());

        goodtrackindices.push_back(it);
    }
    return goodtrackindices.size();
}
Bool_t AliAnalysisTaskNanoCheck::GoodV0Selection() {
    goodv0indices.clear();
    const UInt_t nV0 = fEvt->GetNumberOfV0s();
    fHistos->FillTH1("hNofV0sDistribution", nV0);

    AliESDv0* v0ESD;
    AliAODv0* v0AOD;
    Double_t TPCNSigProton, TPCNSigAntiProton, TPCNSigPion, TPCNSigAntiPion;
    Double_t DCADistLambda_PV, DCADistLambda, LambdaCPA;
    Double_t v0Position[3];
    Double_t MassV0;
    Double_t radius, lV0TotalMomentum, Length, lLifetime;
    UInt_t isAnti = 0;
    Float_t b[2];
    Float_t bCov[3];
    Double_t DCADist_LambdaPionPV, DCADist_LambdaProtonPV;

    Bool_t AcceptedV0 = kTRUE;
    if (!IsAOD) {  // ESD case
        for (UInt_t it = 0; it < nV0; it++) {
            AcceptedV0 = kTRUE;
            v0ESD = ((AliESDEvent*)fEvt)->GetV0(it);
            
            if (fOnlyUseOnTheFlyV0 && !v0ESD->GetOnFlyStatus()) 
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
            TPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            TPCNSigAntiProton =
                GetTPCnSigma(nTrackV0, AliPID::kProton);
            TPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            TPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            if ((TMath::Abs(TPCNSigProton) < fTPCNsigLambdaProtonCut) 
                && (TMath::Abs(TPCNSigPion) < fTPCNsigLambdaPionCut) ){
                
                v0ESD->ChangeMassHypothesis(kLambda0);
                isAnti = 0;

                fHistos->FillTH2("QA_V0/hTPCPIDLambdaProton",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA_V0/hTPCPIDLambdaPion",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
            }
            else if ((TMath::Abs(TPCNSigAntiProton) < fTPCNsigLambdaProtonCut) 
                    && (TMath::Abs(TPCNSigAntiPion) < fTPCNsigLambdaPionCut) ){

                v0ESD->ChangeMassHypothesis(kLambda0Bar);
                isAnti = 1;

                fHistos->FillTH2("QA_V0/hTPCPIDAntiLambdaProton",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA_V0/hTPCPIDAntiLambdaPion",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
            }
            else
                AcceptedV0 = kFALSE;

            // DCA cut
            // DCA between Dautgher particles
            DCADistLambda = TMath::Abs(v0ESD->GetDcaV0Daughters());
            fHistos->FillTH1("QA_V0/hDCA_lambdaDaughters", DCADistLambda);
            if (DCADistLambda > fDCADistLambdaDaughtersCut)
                AcceptedV0 = kFALSE;

            // DCA to PV
            DCADistLambda_PV = TMath::Abs(v0ESD->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));
            fHistos->FillTH1("QA_V0/hDCAlambdaPV", DCADistLambda_PV);
            if (DCADistLambda_PV > fDCArDistLambdaPVCut)
                AcceptedV0 = kFALSE;
            if (!isAnti)
                GetImpactParam(nTrackV0, b, bCov);
            else
                GetImpactParam(pTrackV0, b, bCov);
            DCADist_LambdaPionPV = b[0];
            fHistos->FillTH1("QA_V0/hDCAlambdaPionPV", DCADist_LambdaPionPV);
            if (!isAnti)
                GetImpactParam(pTrackV0, b, bCov);
            else
                GetImpactParam(nTrackV0, b, bCov);
            DCADist_LambdaProtonPV = b[0];
            fHistos->FillTH1("QA_V0/hDCAlambdaProtonPV", DCADist_LambdaProtonPV);

            // CPA cut
            LambdaCPA =
                TMath::Abs(v0ESD->GetV0CosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]));
            fHistos->FillTH1("QA_V0/hCosPAlambda", LambdaCPA);
            if (LambdaCPA < fV0CosineOfPointingAngleCut || LambdaCPA >= 1)
                AcceptedV0 = kFALSE;
            
            // Rapidity cut
            fHistos->FillTH1("QA_V0/hYlambda", v0ESD->RapLambda());
            if (TMath::Abs(v0ESD->RapLambda()) > fMaxLambdaRapidity)
                AcceptedV0 = kFALSE;

            // Radius cut
            v0ESD->GetXYZ(v0Position[0], v0Position[1], v0Position[2]);
            radius = TMath::Hypot(v0Position[0], v0Position[1]);
            fHistos->FillTH1("QA_V0/hLambdaRxy", radius);
            if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
                AcceptedV0 = kFALSE;

            // Life time cut
            lV0TotalMomentum = v0ESD->P();
            fHistos->FillTH1("QA_V0/hTotMomQA", lV0TotalMomentum);
            Length = TMath::Sqrt(TMath::Power(v0Position[0] - lPosPV[0], 2) +
                                  TMath::Power(v0Position[1] - lPosPV[1], 2) +
                                  TMath::Power(v0Position[2] - lPosPV[2], 2));
            fHistos->FillTH1("QA_V0/hDecayLengthQA", Length);
            lLifetime = TMath::Abs(v0Mass * Length / lV0TotalMomentum);
            if (lLifetime > fLambdaLifetime)
                AcceptedV0 = kFALSE;

            // Mass window cut
            MassV0 = v0ESD->GetEffMass();
            fHistos->FillTH1("QA_V0/hMassLambda", MassV0);
            if (TMath::Abs(MassV0 - v0Mass) > fV0MassWindowCut)
                AcceptedV0 = kFALSE;

            // After selection above
            if (AcceptedV0){
                fHistos->FillTH1("hNofV0s", 1.5);
                goodv0indices.push_back({it, isAnti});  // for standard V0
                if (!isAnti) {
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDLambdaProton",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDLambdaPion",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                }
                else {
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDAntiLambdaProton",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDAntiLambdaPion",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                }
                fHistos->FillTH1("QA_V0Cut/hDCA_lambdaDaughters", DCADistLambda);
                fHistos->FillTH1("QA_V0Cut/hDCAlambdaPV", DCADistLambda_PV);
                fHistos->FillTH1("QA_V0Cut/hDCAlambdaPionPV", DCADist_LambdaPionPV);
                fHistos->FillTH1("QA_V0Cut/hDCAlambdaProtonPV", DCADist_LambdaProtonPV);
                fHistos->FillTH1("QA_V0Cut/hCosPAlambda", LambdaCPA);
                fHistos->FillTH1("QA_V0Cut/hYlambda", v0ESD->RapLambda());
                fHistos->FillTH1("QA_V0Cut/hLambdaRxy", radius);
                fHistos->FillTH1("QA_V0Cut/hTotMomQA", lV0TotalMomentum);
                fHistos->FillTH1("QA_V0Cut/hMassLambda", MassV0);
                fHistos->FillTH1("QA_V0Cut/hDecayLengthQA", Length);
            }
        }                                     // All V0 loop
    }                                         // ESD case
    else {
        for (UInt_t it = 0; it < nV0; it++) {
            AcceptedV0 = kTRUE;
            v0AOD = ((AliAODEvent*)fEvt)->GetV0(it);

            if (fOnlyUseOnTheFlyV0 && !v0AOD->GetOnFlyStatus()) 
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
            TPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
            TPCNSigAntiProton =
                GetTPCnSigma(nTrackV0, AliPID::kProton);
            TPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
            TPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

            if ((TMath::Abs(TPCNSigProton) <= fTPCNsigLambdaProtonCut) 
                && (TMath::Abs(TPCNSigPion) <= fTPCNsigLambdaPionCut) ){
                
                isAnti = 0;

                fHistos->FillTH2("QA_V0/hTPCPIDLambdaProton",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA_V0/hTPCPIDLambdaPion",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
            }
            else if ((TMath::Abs(TPCNSigAntiProton) <= fTPCNsigLambdaProtonCut) 
                    && (TMath::Abs(TPCNSigAntiPion) <= fTPCNsigLambdaPionCut) ){

                isAnti = 1;

                fHistos->FillTH2("QA_V0/hTPCPIDAntiLambdaProton",
                                 nTrackV0->GetTPCmomentum(),
                                 nTrackV0->GetTPCsignal());
                fHistos->FillTH2("QA_V0/hTPCPIDAntiLambdaPion",
                                 pTrackV0->GetTPCmomentum(),
                                 pTrackV0->GetTPCsignal());
            }
            else
                AcceptedV0 = kFALSE;

            // DCA cut
            DCADistLambda = TMath::Abs(v0AOD->DcaV0Daughters());
            fHistos->FillTH1("QA_V0/hDCA_lambdaDaughters", DCADistLambda);
            if (DCADistLambda > fDCADistLambdaDaughtersCut)
                AcceptedV0 = kFALSE;  // DCA proton-pion

            // DCA to PV
            DCADistLambda_PV = TMath::Abs(v0AOD->DcaV0ToPrimVertex());
            fHistos->FillTH1("QA_V0/hDCAlambdaPV", DCADistLambda_PV);
            if (DCADistLambda_PV > fDCArDistLambdaPVCut)
                AcceptedV0 = kFALSE;
            if (!isAnti){
                DCADist_LambdaPionPV = TMath::Abs(v0AOD->DcaNegToPrimVertex());
                DCADist_LambdaProtonPV = TMath::Abs(v0AOD->DcaPosToPrimVertex());
            }
            else {
                DCADist_LambdaPionPV = TMath::Abs(v0AOD->DcaPosToPrimVertex());
                DCADist_LambdaProtonPV = TMath::Abs(v0AOD->DcaNegToPrimVertex());
            }
            fHistos->FillTH1("QA_V0/hDCAlambdaPionPV", DCADist_LambdaPionPV);
            fHistos->FillTH1("QA_V0/hDCAlambdaProtonPV", DCADist_LambdaProtonPV);

            // CPA cut
            LambdaCPA = TMath::Abs(v0AOD->CosPointingAngle(vertex));
            fHistos->FillTH1("QA_V0/hCosPAlambda", LambdaCPA);
            if (LambdaCPA < fV0CosineOfPointingAngleCut || LambdaCPA >= 1)
                AcceptedV0 = kFALSE;

            // Rapidity cut
            fHistos->FillTH1("QA_V0/hYlambda", v0AOD->RapLambda());
            if (TMath::Abs(v0AOD->RapLambda()) > fMaxLambdaRapidity)
                AcceptedV0 = kFALSE;

            // Radius cut
            radius = v0AOD->RadiusV0();
            fHistos->FillTH1("QA_V0/hLambdaRxy", radius);
            if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
                AcceptedV0 = kFALSE;

            // Life time cut
            lV0TotalMomentum = TMath::Sqrt(v0AOD->Ptot2V0());
            fHistos->FillTH1("QA_V0/hTotMomQA", lV0TotalMomentum);
            Length = v0AOD->DecayLength(lPosPV);
            fHistos->FillTH1("QA_V0/hDecayLengthQA", Length);
            lLifetime = TMath::Abs(v0Mass * Length / lV0TotalMomentum);
            if (lLifetime > fLambdaLifetime)
                AcceptedV0 = kFALSE;

            // Mass window cut
            MassV0 = -999;
            if (!isAnti)
                MassV0 = v0AOD->MassLambda();
            else
                MassV0 = v0AOD->MassAntiLambda();
            fHistos->FillTH1("QA_V0/hMassLambda", MassV0);
            if (TMath::Abs(MassV0 - v0Mass) > fV0MassWindowCut)
                AcceptedV0 = kFALSE;

            // After selection above
            if (AcceptedV0) {
                fHistos->FillTH1("hNofV0s", 1.5);
                goodv0indices.push_back({it, isAnti});
                if (!isAnti) {
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDLambdaProton",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDLambdaPion",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                }
                else {
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDAntiLambdaProton",
                                     nTrackV0->GetTPCmomentum(),
                                     nTrackV0->GetTPCsignal());
                    fHistos->FillTH2("QA_V0Cut/hTPCPIDAntiLambdaPion",
                                     pTrackV0->GetTPCmomentum(),
                                     pTrackV0->GetTPCsignal());
                }
                fHistos->FillTH1("QA_V0Cut/hDCA_lambdaDaughters", DCADistLambda);
                fHistos->FillTH1("QA_V0Cut/hDCAlambdaPV", DCADistLambda_PV);
                fHistos->FillTH1("QA_V0Cut/hDCAlambdaPionPV", DCADist_LambdaPionPV);
                fHistos->FillTH1("QA_V0Cut/hDCAlambdaProtonPV", DCADist_LambdaProtonPV);
                fHistos->FillTH1("QA_V0Cut/hCosPAlambda", LambdaCPA);
                fHistos->FillTH1("QA_V0Cut/hYlambda", v0AOD->RapLambda());
                fHistos->FillTH1("QA_V0Cut/hLambdaRxy", radius);
                fHistos->FillTH1("QA_V0Cut/hTotMomQA", lV0TotalMomentum);
                fHistos->FillTH1("QA_V0Cut/hMassLambda", MassV0);
                fHistos->FillTH1("QA_V0Cut/hDecayLengthQA", Length);
            }
        }  // All v0 loop
    }      // AOD case

    return goodv0indices.size();
}
Bool_t AliAnalysisTaskNanoCheck::GoodCascadeSelection(){
    goodcascadeindices.clear();
    const UInt_t nXi = fEvt->GetNumberOfCascades();
    fHistos->FillTH1("hNofCascadesDistribution", nXi);

    AliESDcascade* Xicandidate;
    AliAODcascade* Xicandidate_aod;
    Double_t TPCNSigProton, TPCNSigLambdaPion, TPCNSigBachelorPion;
    Double_t pTPCmom, nTPCmom, bTPCmom;
    Double_t pTPCsig, nTPCsig, bTPCsig;
    Double_t DCADist_BachelorPion_PV, DCADist_Lambda_PV, DCADist_Xi_PV;
    Double_t DCADist_LambdaProton_PV, DCADist_LambdaPion_PV;
    Double_t DCADist_Lambda, DCADist_Xi;
    Double_t LambdaCPA, XiCPA;
    Double_t XiPosition[3];
    Double_t XiMomentum[3];
    Double_t LambdaPosition[3];
    Double_t Mass_Xi, Xi_momsum;
    Double_t Xi_eta, Xi_phi;
    UInt_t isXiMinus = 0;
    Float_t b[2];
    Float_t bCov[3];

    Bool_t StandardXi = kTRUE;
    if (!IsAOD) {  // ESD case
        for (UInt_t it = 0; it < nXi; it++) {
            StandardXi = kTRUE;
            Xicandidate = ((AliESDEvent*)fEvt)->GetCascade(it);
            fHistos->FillTH1("hNofCascades", 0.5);
            if (!Xicandidate)
                continue;

            if (TMath::Abs(Xicandidate->GetPindex()) ==
                TMath::Abs(Xicandidate->GetNindex()))
                continue;

            if (TMath::Abs(Xicandidate->GetPindex()) ==
                TMath::Abs(Xicandidate->GetBindex()))
                continue;

            if (TMath::Abs(Xicandidate->GetNindex()) ==
                TMath::Abs(Xicandidate->GetBindex()))
                continue;

            AliESDtrack* pTrackXi = (AliESDtrack*)fEvt->GetTrack(
                TMath::Abs(Xicandidate->GetPindex()));
            AliESDtrack* nTrackXi = (AliESDtrack*)fEvt->GetTrack(
                TMath::Abs(Xicandidate->GetNindex()));
            AliESDtrack* bTrackXi = (AliESDtrack*)fEvt->GetTrack(
                TMath::Abs(Xicandidate->GetBindex()));

            // PID cuts for Xi daughters
            if (Xicandidate->Charge() == -1)  
                isXiMinus = 1;
            else
                isXiMinus = 0;

            if (isXiMinus) {  // Xi- has +proton, -pion
                pTPCmom = pTrackXi->GetTPCmomentum();
                pTPCsig = pTrackXi->GetTPCsignal();
                nTPCmom = nTrackXi->GetTPCmomentum();
                nTPCsig = nTrackXi->GetTPCsignal();
                TPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                TPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            } else {  // Xi+ has -proton, +pion
                nTPCmom = pTrackXi->GetTPCmomentum();
                nTPCsig = pTrackXi->GetTPCsignal();
                pTPCmom = nTrackXi->GetTPCmomentum();
                pTPCsig = nTrackXi->GetTPCsignal();
                TPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                TPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }
            bTPCmom = bTrackXi->GetTPCmomentum();
            bTPCsig = bTrackXi->GetTPCsignal();
            TPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion
            fHistos->FillTH2("QA_Xi/hTPCPIDLambdaProton", pTPCmom, pTPCsig);
            fHistos->FillTH2("QA_Xi/hTPCPIDLambdaPion", nTPCmom, nTPCsig);
            fHistos->FillTH2("QA_Xi/hTPCPIDBachelorPion", bTPCmom, bTPCsig);
            fHistos->FillTH1("QA_Xi/hTPCPIDsignalLambdaProton", TPCNSigProton);
            fHistos->FillTH1("QA_Xi/hTPCPIDsignalLambdaPion",
                             TPCNSigLambdaPion);
            fHistos->FillTH1("QA_Xi/hTPCPIDsignalBachelorPion",
                             TPCNSigBachelorPion);

            if (TMath::Abs(TPCNSigProton) > fTPCNsigLambdaProtonCut)
                StandardXi = kFALSE;  // PID for proton
            if (TMath::Abs(TPCNSigLambdaPion) > fTPCNsigLambdaPionCut)
                StandardXi = kFALSE;  // PID for 1st pion
            if (TMath::Abs(TPCNSigBachelorPion) > fTPCNsigBachelorPionCut)
                StandardXi = kFALSE;  // PID for 2nd pion

            // DCA cut
            // DCA between Dautgher particles
            DCADist_Lambda = TMath::Abs(Xicandidate->GetDcaV0Daughters());
            DCADist_Xi = TMath::Abs(Xicandidate->GetDcaXiDaughters());
            fHistos->FillTH1("QA_Xi/hDCADist_Lambda_BTW_Daughters",
                             DCADist_Lambda);
            fHistos->FillTH1("QA_Xi/hDCADist_Xi_BTW_Daughters", DCADist_Xi);

            if (DCADist_Lambda > fDCADist_LambdaDaughtersCut)
                StandardXi = kFALSE;  // DCA proton-pion
            if (DCADist_Xi > fDCADist_XiDaughtersCut)
                StandardXi = kFALSE;  // DCA Lambda-pion

            // DCA to PV
            DCADist_Lambda_PV =
                TMath::Abs(Xicandidate->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));
            DCADist_Xi_PV = TMath::Abs(
                Xicandidate->GetDcascade(lPosPV[0], lPosPV[1], lPosPV[2]));
            if (isXiMinus) {  // Xi- has +proton, -pion
                GetImpactParam(pTrackXi, b, bCov);
                DCADist_LambdaProton_PV = b[0];
                GetImpactParam(nTrackXi, b, bCov);
                DCADist_LambdaPion_PV = b[0];
            } else {
                GetImpactParam(nTrackXi, b, bCov);
                DCADist_LambdaProton_PV = b[0];
                GetImpactParam(pTrackXi, b, bCov);
                DCADist_LambdaPion_PV = b[0];
            }
            GetImpactParam(bTrackXi, b, bCov);
            DCADist_BachelorPion_PV = b[0];
            fHistos->FillTH1("QA_Xi/hDCADist_lambda_to_PV",
                                DCADist_Lambda_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_Xi_to_PV", DCADist_Xi_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_LambdaProton_to_PV",
                                DCADist_LambdaProton_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_LambdaPion_to_PV",
                                DCADist_LambdaPion_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_BachelorPion_to_PV",
                                DCADist_BachelorPion_PV);

            if (DCADist_Lambda_PV < fDCADist_Lambda_PVCut)
                StandardXi = kFALSE;  // DCA Lambda-vertex

            // CPA cut
            Xicandidate->GetXYZcascade(XiPosition[0], XiPosition[1],
                                       XiPosition[2]);
            LambdaCPA = Xicandidate->GetV0CosineOfPointingAngle(
                XiPosition[0], XiPosition[1], XiPosition[2]);
            // LambdaCPA = Xicandidate->GetV0CosineOfPointingAngle(lPosPV[0],
            // lPosPV[1], lPosPV[2]); // PV check
            XiCPA = Xicandidate->GetCascadeCosineOfPointingAngle(
                lPosPV[0], lPosPV[1], lPosPV[2]);
            fHistos->FillTH1("QA_Xi/hCosPA_lambda", LambdaCPA);
            fHistos->FillTH1("QA_Xi/hCosPA_Xi", XiCPA);

            if (LambdaCPA < fV0CosineOfPointingAngleCut)
                StandardXi = kFALSE;
            if (XiCPA < fCascadeCosineOfPointingAngleCut)
                StandardXi = kFALSE;

            // Mass window cut
            Mass_Xi = Xicandidate->GetEffMassXi();
            fHistos->FillTH1("QA_Xi/hMass_Xi", Mass_Xi);
            /*
            // Disable for Xi study
            if (TMath::Abs(Mass_Xi - Ximass) > fXiMassWindowCut)
                StandardXi = kFALSE;
            */
            // Eta cut
            Xi_eta = Xicandidate->Eta();
            Xi_phi = Xicandidate->Phi();
            if (TMath::Abs(Xi_eta) > fXiEtaCut)
                StandardXi = kFALSE;
            fHistos->FillTH2("QA_Xi/hPhiEta_Xi", Xi_phi, Xi_eta);

            // XY Raidus cut(experiemntal)
            Xicandidate->GetXYZ(LambdaPosition[0], LambdaPosition[1],
                                LambdaPosition[2]);
            fHistos->FillTH2("QA_Xi/hLambda_Rxy", LambdaPosition[0],
                             LambdaPosition[1]);
            // if(sqrt( pow(LambdaX,2) + pow(LambdaY,2) ) > 100)
            // StandardXi=kFALSE; // NOT USING
            fHistos->FillTH2("QA_Xi/hXi_Rxy", XiPosition[0], XiPosition[1]);

            // After selection above
            if (StandardXi) {  // Save only the Xi is good candidate
                fHistos->FillTH2("QA_XiCut/hTPCPIDLambdaProton", pTPCmom,
                                 pTPCsig);
                fHistos->FillTH2("QA_XiCut/hTPCPIDLambdaPion", nTPCmom,
                                 nTPCsig);
                fHistos->FillTH2("QA_XiCut/hTPCPIDBachelorPion", bTPCmom,
                                 bTPCsig);
                fHistos->FillTH1("QA_XiCut/hTPCPIDsignalLambdaProton",
                                 TPCNSigProton);
                fHistos->FillTH1("QA_XiCut/hTPCPIDsignalLambdaPion",
                                 TPCNSigLambdaPion);
                fHistos->FillTH1("QA_XiCut/hTPCPIDsignalBachelorPion",
                                 TPCNSigBachelorPion);
                fHistos->FillTH1("QA_XiCut/hDCADist_Lambda_BTW_Daughters",
                                 DCADist_Lambda);
                fHistos->FillTH1("QA_XiCut/hDCADist_Xi_BTW_Daughters",
                                 DCADist_Xi);
                fHistos->FillTH1("QA_XiCut/hDCADist_lambda_to_PV",
                                 DCADist_Lambda_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_Xi_to_PV", DCADist_Xi_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_LambdaProton_to_PV",
                                 DCADist_LambdaProton_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_LambdaPion_to_PV",
                                 DCADist_LambdaPion_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_BachelorPion_to_PV",
                                 DCADist_BachelorPion_PV);
                fHistos->FillTH1("QA_XiCut/hCosPA_lambda", LambdaCPA);
                fHistos->FillTH1("QA_XiCut/hCosPA_Xi", XiCPA);
                fHistos->FillTH1("QA_XiCut/hMass_Xi", Mass_Xi);
                fHistos->FillTH2("QA_XiCut/hPhiEta_Xi", Xi_phi, Xi_eta);
                fHistos->FillTH2("QA_XiCut/hLambda_Rxy", LambdaPosition[0],
                                 LambdaPosition[1]);
                fHistos->FillTH2("QA_XiCut/hXi_Rxy", XiPosition[0],
                                 XiPosition[1]);

                goodcascadeindices.push_back({it, isXiMinus});
                fHistos->FillTH1("hNofCascades", 1.5);
            }  // for standard Xi
        }      // All Xi loop
    }          // ESD case
    else {
        for (UInt_t it = 0; it < nXi; it++) {
            StandardXi = kTRUE;
            Xicandidate_aod = ((AliAODEvent*)fEvt)->GetCascade(it);
            fHistos->FillTH1("hNofCascades", 0.5);
            if (!Xicandidate_aod)
                continue;

            if (TMath::Abs(Xicandidate_aod->GetPosID()) ==
                TMath::Abs(Xicandidate_aod->GetNegID()))
                continue;
            if (TMath::Abs(Xicandidate_aod->GetPosID()) ==
                TMath::Abs(Xicandidate_aod->GetBachID()))
                continue;
            if (TMath::Abs(Xicandidate_aod->GetNegID()) ==
                TMath::Abs(Xicandidate_aod->GetBachID()))
                continue;

            AliAODTrack* pTrackXi =
                (AliAODTrack*)(Xicandidate_aod->GetDaughter(0));
            AliAODTrack* nTrackXi =
                (AliAODTrack*)(Xicandidate_aod->GetDaughter(1));
            AliAODTrack* bTrackXi =
                (AliAODTrack*)(Xicandidate_aod->GetDecayVertexXi()->GetDaughter(
                    0));

            // PID cuts for Xi daughters
            if (Xicandidate_aod->ChargeXi() == -1)
                isXiMinus = 1;
            else
                isXiMinus = 0;

            if (isXiMinus) {  // Xi- has +proton, -pion
                pTPCmom = pTrackXi->GetTPCmomentum();
                pTPCsig = pTrackXi->GetTPCsignal();
                nTPCmom = nTrackXi->GetTPCmomentum();
                nTPCsig = nTrackXi->GetTPCsignal();
                TPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                TPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            } else {  // Xi+ has -proton, +pion
                nTPCmom = pTrackXi->GetTPCmomentum();
                nTPCsig = pTrackXi->GetTPCsignal();
                pTPCmom = nTrackXi->GetTPCmomentum();
                pTPCsig = nTrackXi->GetTPCsignal();
                TPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                TPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }
            bTPCmom = bTrackXi->GetTPCmomentum();
            bTPCsig = bTrackXi->GetTPCsignal();
            TPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion
            fHistos->FillTH2("QA_Xi/hTPCPIDLambdaProton", pTPCmom, pTPCsig);
            fHistos->FillTH2("QA_Xi/hTPCPIDLambdaPion", nTPCmom, nTPCsig);
            fHistos->FillTH2("QA_Xi/hTPCPIDBachelorPion", bTPCmom, bTPCsig);
            fHistos->FillTH1("QA_Xi/hTPCPIDsignalLambdaProton", TPCNSigProton);
            fHistos->FillTH1("QA_Xi/hTPCPIDsignalLambdaPion",
                             TPCNSigLambdaPion);
            fHistos->FillTH1("QA_Xi/hTPCPIDsignalBachelorPion",
                             TPCNSigBachelorPion);

            if (TMath::Abs(TPCNSigProton) > fTPCNsigLambdaProtonCut) 
                StandardXi = kFALSE;  // PID for proton
            if (TMath::Abs(TPCNSigLambdaPion) > fTPCNsigLambdaPionCut)
                StandardXi = kFALSE;  // PID for 1st pion
            if (TMath::Abs(TPCNSigBachelorPion) > fTPCNsigBachelorPionCut)
                StandardXi = kFALSE;  // PID for 2nd pion

            // DCA cut
            // DCA between Dautgher particles
            DCADist_Lambda = TMath::Abs(Xicandidate_aod->DcaV0Daughters());
            DCADist_Xi = TMath::Abs(Xicandidate_aod->DcaXiDaughters());
            fHistos->FillTH1("QA_Xi/hDCADist_Lambda_BTW_Daughters", DCADist_Lambda);
            fHistos->FillTH1("QA_Xi/hDCADist_Xi_BTW_Daughters", DCADist_Xi);

            if (DCADist_Lambda > fDCADist_LambdaDaughtersCut)
                StandardXi = kFALSE;  // DCA proton-pion
            if (DCADist_Xi > fDCADist_XiDaughtersCut)
                StandardXi = kFALSE;  // DCA Lambda-pion

            // DCA to PV
            DCADist_Lambda_PV =
                TMath::Abs(Xicandidate_aod->DcaV0ToPrimVertex());
            DCADist_Xi_PV = TMath::Abs(Xicandidate_aod->DcaXiToPrimVertex(
                lPosPV[0], lPosPV[1], lPosPV[2]));
            if (isXiMinus) {  // Xi- has +proton, -pion
                GetImpactParam(pTrackXi, b, bCov);
                DCADist_LambdaProton_PV = b[0];
                GetImpactParam(nTrackXi, b, bCov);
                DCADist_LambdaPion_PV = b[0];
            } else {
                GetImpactParam(nTrackXi, b, bCov);
                DCADist_LambdaProton_PV = b[0];
                GetImpactParam(pTrackXi, b, bCov);
                DCADist_LambdaPion_PV = b[0];
            }
            GetImpactParam(bTrackXi, b, bCov);
            DCADist_BachelorPion_PV = b[0];
            /*
            // typical AOD method
            if (isXiMinus) {
                DCADist_LambdaProton_PV =
                    TMath::Abs(Xicandidate_aod->DcaPosToPrimVertex());
                DCADist_LambdaPion_PV =
                    TMath::Abs(Xicandidate_aod->DcaNegToPrimVertex());
            } else {
                DCADist_LambdaProton_PV =
                    TMath::Abs(Xicandidate_aod->DcaNegToPrimVertex());
                DCADist_LambdaPion_PV =
                    TMath::Abs(Xicandidate_aod->DcaPosToPrimVertex());
            }
            
            DCADist_BachelorPion_PV =
                TMath::Abs(Xicandidate_aod->DcaBachToPrimVertex());
            */
            fHistos->FillTH1("QA_Xi/hDCADist_lambda_to_PV", DCADist_Lambda_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_Xi_to_PV", DCADist_Xi_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_LambdaProton_to_PV",
                             DCADist_LambdaProton_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_LambdaPion_to_PV",
                             DCADist_LambdaPion_PV);
            fHistos->FillTH1("QA_Xi/hDCADist_BachelorPion_to_PV",
                             DCADist_BachelorPion_PV);

            if (DCADist_Lambda_PV < fDCADist_Lambda_PVCut)
                StandardXi = kFALSE;  // DCA Lambda-vertex

            // CPA cut
            XiPosition[0] = Xicandidate_aod->DecayVertexXiX();
            XiPosition[1] = Xicandidate_aod->DecayVertexXiY();
            XiPosition[2] = Xicandidate_aod->DecayVertexXiZ();
            LambdaCPA = Xicandidate_aod->CosPointingAngle(XiPosition);
            // LambdaCPA = Xicandidate_aod->CosPointingAngle(lPosPV);
            XiCPA = Xicandidate_aod->CosPointingAngleXi(lPosPV[0], lPosPV[1],
                                                         lPosPV[2]);
            fHistos->FillTH1("QA_Xi/hCosPA_lambda", LambdaCPA);
            fHistos->FillTH1("QA_Xi/hCosPA_Xi", XiCPA);

            if (LambdaCPA < fV0CosineOfPointingAngleCut)
                StandardXi = kFALSE;
            if (XiCPA < fCascadeCosineOfPointingAngleCut)
                StandardXi = kFALSE;

            // Mass window cut
            Mass_Xi = Xicandidate_aod->MassXi();
            fHistos->FillTH1("QA_Xi/hMass_Xi", Mass_Xi);
            /*
            // Disable for Xi study
            if (TMath::Abs(Mass_Xi - Ximass) > fXiMassWindowCut)
                StandardXi = kFALSE;
            */

            // Eta cut
            // Eta: 0.5*TMath::Log((P()+Pz())/(P()-Pz()+1.e-13))
            // Phi: TMath::Pi()+TMath::ATan2(-Py(),-Px())
            XiMomentum[0] = Xicandidate_aod->MomXiX();
            XiMomentum[1] = Xicandidate_aod->MomXiY();
            XiMomentum[2] = Xicandidate_aod->MomXiZ();
            Xi_momsum =
                TMath::Sqrt(XiMomentum[0] * XiMomentum[0] + XiMomentum[1] * XiMomentum[1] +
                            XiMomentum[2] * XiMomentum[2]);
            Xi_eta = 0.5 * TMath::Log((Xi_momsum + XiMomentum[2]) /
                                      (Xi_momsum - XiMomentum[2] + 1.e-13));
            Xi_phi = TMath::Pi() + TMath::ATan2(-XiMomentum[1], -XiMomentum[0]);
            if (TMath::Abs(Xi_eta) > fXiEtaCut)
                StandardXi = kFALSE;
            fHistos->FillTH2("QA_Xi/hPhiEta_Xi", Xi_phi, Xi_eta);

            // XY Raidus cut(experiemntal)
            // Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);
            fHistos->FillTH2("QA_Xi/hLambda_Rxy", Xicandidate_aod->DecayVertexV0X(),
                             Xicandidate_aod->DecayVertexV0Y());
            // if(sqrt( Xicandidate->RadiusV0() ) > 100)
            // StandardXi=kFALSE; // NOT USING

            fHistos->FillTH2("QA_Xi/hXi_Rxy", XiPosition[0], XiPosition[1]);

            // After selection above
            if (StandardXi) {  // Save only the Xi is good candidate
                fHistos->FillTH2("QA_XiCut/hTPCPIDLambdaProton", pTPCmom,
                                 pTPCsig);
                fHistos->FillTH2("QA_XiCut/hTPCPIDLambdaPion", nTPCmom,
                                 nTPCsig);
                fHistos->FillTH2("QA_XiCut/hTPCPIDBachelorPion", bTPCmom,
                                 bTPCsig);
                fHistos->FillTH1("QA_XiCut/hTPCPIDsignalLambdaProton",
                                 TPCNSigProton);
                fHistos->FillTH1("QA_XiCut/hTPCPIDsignalLambdaPion",
                                 TPCNSigLambdaPion);
                fHistos->FillTH1("QA_XiCut/hTPCPIDsignalBachelorPion",
                                 TPCNSigBachelorPion);
                fHistos->FillTH1("QA_XiCut/hDCADist_Lambda_BTW_Daughters",
                                 DCADist_Lambda);
                fHistos->FillTH1("QA_XiCut/hDCADist_Xi_BTW_Daughters",
                                 DCADist_Xi);
                fHistos->FillTH1("QA_XiCut/hDCADist_lambda_to_PV",
                                 DCADist_Lambda_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_Xi_to_PV", DCADist_Xi_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_LambdaProton_to_PV",
                                 DCADist_LambdaProton_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_LambdaPion_to_PV",
                                 DCADist_LambdaPion_PV);
                fHistos->FillTH1("QA_XiCut/hDCADist_BachelorPion_to_PV",
                                 DCADist_BachelorPion_PV);
                fHistos->FillTH1("QA_XiCut/hCosPA_lambda", LambdaCPA);
                fHistos->FillTH1("QA_XiCut/hCosPA_Xi", XiCPA);
                fHistos->FillTH1("QA_XiCut/hMass_Xi", Mass_Xi);
                fHistos->FillTH2("QA_XiCut/hPhiEta_Xi", Xi_phi, Xi_eta);
                fHistos->FillTH2("QA_XiCut/hLambda_Rxy", LambdaPosition[0],
                                 LambdaPosition[1]);
                fHistos->FillTH2("QA_XiCut/hXi_Rxy", XiPosition[0],
                                 XiPosition[1]);
                goodcascadeindices.push_back({it,isXiMinus});
                fHistos->FillTH1("hNofCascades", 1.5);
            }  // for standard Xi
        }      // All Xi loop
    }          // AOD case

    return goodcascadeindices.size();
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
void AliAnalysisTaskNanoCheck::GetImpactParam(AliVTrack* track, Float_t p[2], Float_t cov[3]) {
    AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
    if (nanoT)
        nanoT->GetImpactParameters(p[0], p[1]);
    else
        track->GetImpactParameters(p, cov);
}