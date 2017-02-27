/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron Flow with TPC plus TOF            //
//  Non-Photonic Electron identified with Invariant mass              //
//  analysis methos in function  SelectPhotonicElectron               //
//                                                                    //
//                                                                    //
//  Author: Andrea Dubla (Utrecht University)                         //
//                                                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskFlowITSTPCTOFQCSP.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAODInputHandler.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliCentralitySelectionTask.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloTrigger.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"
#include "AliEMCALTrack.h"
//#include "AliEMCALTracker.h"
#include "AliMagF.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliCentrality.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliFlowTrack.h"
#include "AliAnalysisTaskVnV0.h"
#include "AliSelectNonHFE.h"


class AliFlowTrackCuts;

using namespace std;

ClassImp(AliAnalysisTaskFlowITSTPCTOFQCSP)
//________________________________________________________________________
AliAnalysisTaskFlowITSTPCTOFQCSP::AliAnalysisTaskFlowITSTPCTOFQCSP(const char *name)
: AliAnalysisTaskSE(name)
,fDebug(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fCuts(0)
,fIdentifiedAsOutInz(kFALSE)
,fPassTheEventCut(kFALSE)
,fCFM(0)
,fPID(0)
,ftpcpid(0)
,fPIDqa(0)
,fCutsRP(0)     // track cuts for reference particles
,fNullCuts(0) // dummy cuts for flow event tracks
,fFlowEvent(0) //! flow events (one for each inv mass band)
,fkCentralityMethod(0)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(0)
,fInvmassCut(0)
,fpTCutmin(0)
,fpTCutmax(5)
,fTrigger(0)
,fPhi(0)
,fEta(0)
,fVZEROA(0)
,fVZEROC(0)
,fTPCM(0)
,fNoEvents(0)
,fInclusiveElecPt(0)
,fTPCnsigma(0)
,fTPCnsigmaAft(0)
,fITSnsigmaAft(0)
,fTPCnsigmaVSptAft(0)
,fTOFns(0)
,fTOFnsAft(0)
,fTOFBetaAft(0)
,fCentralityPass(0)
,fCentralityNoPass(0)
,fInvmassLS1(0)
,fInvmassULS1(0)
,fPhotoElecPt(0)
,fSemiInclElecPt(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fMultCorAfterCuts(0)
,fMultvsCentr(0)
,fSubEventDPhiv2(0)
,EPVzA(0)
,EPVzC(0)
,EPTPC(0)
,fV2Phi(0)
,fvertex(0)
,fMultCorBeforeCuts(0)
,fQAPid(0)
,fminITSnsigmaLowpT(-1)
,fmaxITSnsigmaLowpT(1)
,fminITSnsigmaHighpT(-2)
,fmaxITSnsigmaHighpT(2)
,fminTPCnsigmaLowpT(-1)
,fmaxTPCnsigmaLowpT(3)
,fminTPCnsigmaHighpT(0)
,fmaxTPCnsigmaHighpT(3)
//,fQAPIDSparse(kFALSE)
,fminTOFnSigma(-2)
,fmaxTOFnSigma(2)
,fQAPidSparse(0)
,fTPCS(0)
,fVz(0)
,fOpeningAngleCut(0.1)
,fOP_angle(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fNonHFE(new AliSelectNonHFE)
,fDCA(0)
,fITSnsigma(0)
,fITSnsigmaElect(0)
,fTPCnsigmaAftITSTOF(0)
,fTPCnsigmaAftTOF(0)
,fITSnsigmaAftTOF(0)
,fITSvsTOF(0)
,fTPCvsITS(0)
,fTPCvsITSafterTOF(0)
,fTPCvsTOF(0)
,fMultCut(0)
,fMultCorAfterCentrBeforeCuts(0)
,fMultCorAfterVZTRKComp(0)
,fCentralityBeforePileup(0)
,fCentralityAfterVZTRK(0)
,fCentralityAfterCorrCut(0)
,fMultCorAfterCorrCut(0)
,EPVz(0)
,EPTPCp(0)
,EPTPCn(0)
,fSubEventDPhiv2new(0)
,fV2Phivzerotot(0)
,fHistCentrDistr(0x0)
,fCentralityNoPassForFlattening(0)
,fptminAsso(0)
,fSparsephipsiULS(0)
,fSparsephipsiLS(0)
,fSparseMassULS(0)
,fSparseMassLS(0)
,fAssoTPCCluster(0)
,fAssoITSRefit(0)
,fPhiCut(0)
,fHistEPDistrWeight(0)
,EPweights(0)
,EPVzAftW(0)
,multCorrection(0)
,fEtaMinimumPositive(0)
,fEtaMinimumNegative(0)
,fCentralityAll(0)
,fCentFlatMine(0)
,fmineta(-0.8)
,fmaxeta(0.8)
{
    //Named constructor
    
    fPID = new AliHFEpid("hfePid");
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, AliFlowEventSimple::Class());
}

//________________________________________________________________________
AliAnalysisTaskFlowITSTPCTOFQCSP::AliAnalysisTaskFlowITSTPCTOFQCSP()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElectFlow")
,fDebug(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fCuts(0)
,fIdentifiedAsOutInz(kFALSE)
,fPassTheEventCut(kFALSE)
,fCFM(0)
,fPID(0)
,ftpcpid(0)
,fPIDqa(0)
,fCutsRP(0)     // track cuts for reference particles
,fNullCuts(0) // dummy cuts for flow event tracks
,fFlowEvent(0) //! flow events (one for each inv mass band)
,fkCentralityMethod(0)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(0)
,fInvmassCut(0)
,fpTCutmin(0)
,fpTCutmax(5)
,fTrigger(0)
,fPhi(0)
,fEta(0)
,fVZEROA(0)
,fVZEROC(0)
,fTPCM(0)
,fNoEvents(0)
,fInclusiveElecPt(0)
,fTPCnsigma(0)
,fTPCnsigmaAft(0)
,fITSnsigmaAft(0)
,fTPCnsigmaVSptAft(0)
,fTOFns(0)
,fTOFnsAft(0)
,fTOFBetaAft(0)
,fCentralityPass(0)
,fCentralityNoPass(0)
,fInvmassLS1(0)
,fInvmassULS1(0)
,fPhotoElecPt(0)
,fSemiInclElecPt(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fMultCorAfterCuts(0)
,fMultvsCentr(0)
,fSubEventDPhiv2(0)
,EPVzA(0)
,EPVzC(0)
,EPTPC(0)
,fV2Phi(0)
,fvertex(0)
,fMultCorBeforeCuts(0)
,fQAPid(0)
,fminITSnsigmaLowpT(-1)
,fmaxITSnsigmaLowpT(1)
,fminITSnsigmaHighpT(-2)
,fmaxITSnsigmaHighpT(2)
,fminTPCnsigmaLowpT(-1)
,fmaxTPCnsigmaLowpT(3)
,fminTPCnsigmaHighpT(0)
,fmaxTPCnsigmaHighpT(3)
//,fQAPIDSparse(kFALSE)
,fminTOFnSigma(-2)
,fmaxTOFnSigma(2)
,fQAPidSparse(0)
,fTPCS(0)
,fVz(0)
,fOpeningAngleCut(0.1)
,fOP_angle(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fNonHFE(new AliSelectNonHFE)
,fDCA(0)
,fITSnsigma(0)
,fITSnsigmaElect(0)
,fTPCnsigmaAftITSTOF(0)
,fTPCnsigmaAftTOF(0)
,fITSnsigmaAftTOF(0)
,fITSvsTOF(0)
,fTPCvsITS(0)
,fTPCvsITSafterTOF(0)
,fTPCvsTOF(0)
,fMultCut(0)
,fMultCorAfterCentrBeforeCuts(0)
,fMultCorAfterVZTRKComp(0)
,fCentralityBeforePileup(0)
,fCentralityAfterVZTRK(0)
,fCentralityAfterCorrCut(0)
,fMultCorAfterCorrCut(0)
,EPVz(0)
,EPTPCp(0)
,EPTPCn(0)
,fSubEventDPhiv2new(0)
,fV2Phivzerotot(0)
,fHistCentrDistr(0x0)
,fCentralityNoPassForFlattening(0)
,fptminAsso(0)
,fSparsephipsiULS(0)
,fSparsephipsiLS(0)
,fSparseMassULS(0)
,fSparseMassLS(0)
,fAssoTPCCluster(0)
,fAssoITSRefit(0)
,fPhiCut(0)
,fHistEPDistrWeight(0)
,EPweights(0)
,EPVzAftW(0)
,multCorrection(0)
,fEtaMinimumPositive(0)
,fEtaMinimumNegative(0)
,fCentralityAll(0)
,fCentFlatMine(0)
,fmineta(-0.8)
,fmaxeta(0.8)
{
    //Default constructor
    fPID = new AliHFEpid("hfePid");
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, AliFlowEventSimple::Class());
    //DefineOutput(3, TTree::Class());
}
//_________________________________________

AliAnalysisTaskFlowITSTPCTOFQCSP::~AliAnalysisTaskFlowITSTPCTOFQCSP()
{
    //Destructor
    
    delete fOutputList;
    delete fPID;
    //  delete fPIDResponse;
    delete fCFM;
    delete fPIDqa;
    if (fOutputList) delete fOutputList;
    if (fFlowEvent) delete fFlowEvent;
    delete fNonHFE;
}
//_________________________________________

void AliAnalysisTaskFlowITSTPCTOFQCSP::UserExec(Option_t*)
{
    //Main loop
    //Called for each event
    
    // create pointer to event
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    
    
    if (!fAOD)
    {
        printf("ERROR: fAOD not available\n");
        return;
    }
    
    if(!fCuts)
    {
        AliError("HFE cuts not available");
        return;
    }
    
    if(!fPID->IsInitialized())
    {
        // Initialize PID with the given run number
        AliWarning("PID not initialised, get from Run no");
        fPID->InitializePID(fAOD->GetRunNumber());
    }
    
    // cout << "kTrigger   ==   " << fTrigger <<endl;
    
    if(fTrigger==0){
        if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral)) return;
    }
    if(fTrigger==1){
        
        if ( !(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny) ) return;
        
        TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
        
        if ( ! ( firedTriggerClasses.Contains("CVLN_B2-B-NOPF-ALLNOTRD") || firedTriggerClasses.Contains("CVLN_R1-B-NOPF-ALLNOTRD") || firedTriggerClasses.Contains("CSEMI_R1-B-NOPF-ALLNOTRD") ) ) return;
    }
    if(fTrigger==2){
        if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMCEGA)) return;
    }
    if(fTrigger==3){
        if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;
    }
    if(fTrigger==4){
        if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kCentral | AliVEvent::kSemiCentral))) return;
    }
    if(fTrigger==5){
        if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kSemiCentral))) return;
    }
    
    
    
    
    //---------------CENTRALITY AND EVENT SELECTION-----------------------
    
    
    
    Int_t fNOtrks =  fAOD->GetNumberOfTracks();
    Float_t vtxz = -999;
    const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
    if (!trkVtx || trkVtx->GetNContributors()<=0)return;
    TString vtxTtl = trkVtx->GetTitle();
    if (!vtxTtl.Contains("VertexerTracks"))return;
    const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
    if (!spdVtx || spdVtx->GetNContributors()<=0)return;
    if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5)return;
    vtxz = trkVtx->GetZ();
    if(TMath::Abs(vtxz)>fVz)return;
    
    // Event cut
    if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fAOD)) return;
    if(fNOtrks<2) return;
    
    
    Bool_t pass = kFALSE; //to select centrality
    CheckCentrality(fAOD,pass);
    if(!pass)return;
    
    fvertex->Fill(vtxz);
    
    
    fNoEvents->Fill(0);
    PlotVZeroMultiplcities(fAOD);
    
    SetNullCuts(fAOD);
    PrepareFlowEvent(fAOD->GetNumberOfTracks(),fFlowEvent);    //Calculate event plane Qvector and EP resolution for inclusive
    
    AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
    if(!pidResponse)
    {
        AliDebug(1, "Using default PID Response");
        pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
    }
    
    fPID->SetPIDResponse(pidResponse);
    
    fCFM->SetRecEventInfo(fAOD);
    
    // Look for kink mother
    Int_t numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
        AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
        if(!aodvertex) continue;
        if(aodvertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[numberofmotherkink] = idmother;
            //printf("ID %d\n",idmother);
            numberofmotherkink++;
        }
    }
    
    //=============================================V0EP from Alex======================================================================
    Double_t qxEPa = 0, qyEPa = 0;
    Double_t qxEPc = 0, qyEPc = 0;
    Double_t qxEP = 0, qyEP = 0;
    
    Double_t evPlAngV0A = fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD, 8, 2, qxEPa, qyEPa);
    Double_t evPlAngV0C = fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD, 9, 2, qxEPc, qyEPc);
    Double_t evPlAngV0 = fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD, 10, 2, qxEP, qyEP);
    
    
    Double_t Qx2 = 0, Qy2 = 0;
    Double_t Qx2p = 0, Qy2p = 0;
    Double_t Qx2n = 0, Qy2n = 0;
    
    for (Int_t iT = 0; iT < fAOD->GetNumberOfTracks(); iT++){
        
        AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iT));
        if(!aodTrack) AliFatal("Not a standard AOD");
        
        if (!aodTrack)
            continue;
        
        if ((TMath::Abs(aodTrack->Eta()) > 0.8) || (aodTrack->Pt() < 0.2) || (aodTrack->GetTPCNcls() < 70) || (aodTrack->Pt() >= 20.0))
            continue;
        
        if (!aodTrack->TestFilterBit(128))
            continue;
        
        
        if(aodTrack->Eta()>fEtaMinimumPositive && aodTrack->Eta()<0.8){
            
            Qx2p += TMath::Cos(2*aodTrack->Phi());
            Qy2p += TMath::Sin(2*aodTrack->Phi());
        }
        if(aodTrack->Eta()<fEtaMinimumNegative && aodTrack->Eta()> -0.8){
            
            Qx2n += TMath::Cos(2*aodTrack->Phi());
            Qy2n += TMath::Sin(2*aodTrack->Phi());
        }
        
        
        Qx2 += TMath::Cos(2*aodTrack->Phi());
        Qy2 += TMath::Sin(2*aodTrack->Phi());
        
        
        
        
    }
    
    Double_t evPlAngTPC = TMath::ATan2(Qy2, Qx2)/2.;
    Double_t evPlAngTPCn = TMath::ATan2(Qy2n, Qx2n)/2.;
    Double_t evPlAngTPCp = TMath::ATan2(Qy2p, Qx2p)/2.;
    
    EPVzA->Fill(evPlAngV0A);
    EPVzC->Fill(evPlAngV0C);
    EPTPC->Fill(evPlAngTPC);
    
    EPTPCn->Fill(evPlAngTPCn);
    EPTPCp->Fill(evPlAngTPCp);
    EPVz->Fill(evPlAngV0);
    
    
    
    Double_t weightEP =1;
    if(EPweights){
        weightEP = GiveMeWeight(evPlAngV0);
        EPVzAftW->Fill(evPlAngV0,weightEP);
        
    }
    
    
    
    fSubEventDPhiv2->Fill(0.5, TMath::Cos(2.*(evPlAngV0A-evPlAngTPC))); // vzeroa - tpc
    fSubEventDPhiv2->Fill(1.5, TMath::Cos(2.*(evPlAngV0A-evPlAngV0C))); // vzeroa - vzeroc
    fSubEventDPhiv2->Fill(2.5, TMath::Cos(2.*(evPlAngV0C-evPlAngTPC))); // tpc - vzeroc
    
    
    if(EPweights){
        fSubEventDPhiv2new->Fill(0.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCp)),weightEP); // vzero - tpcp
        fSubEventDPhiv2new->Fill(1.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCn)),weightEP); // vzero - tpcn
        fSubEventDPhiv2new->Fill(2.5, TMath::Cos(2.*(evPlAngTPCp-evPlAngTPCn))); // tpcp - tpcn
    }
    if(!EPweights){
        fSubEventDPhiv2new->Fill(0.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCp))); // vzero - tpcp
        fSubEventDPhiv2new->Fill(1.5, TMath::Cos(2.*(evPlAngV0-evPlAngTPCn))); // vzero - tpcn
        fSubEventDPhiv2new->Fill(2.5, TMath::Cos(2.*(evPlAngTPCp-evPlAngTPCn))); // tpcp - tpcn
    }
    //====================================================================================================================
    
    AliAODTrack *track = NULL;
    
    // Track loop
    for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++)
    {
        track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
        if(!track) AliFatal("Not a standard AOD");
        if (!track)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;  // TESTBIT FOR AOD double Counting
        
        //--------------------------------------hfe begin-----------------------------------------------------------
        //==========================================================================================================
        //======================================track cuts==========================================================
        if(track->Eta()<fmineta || track->Eta()>fmaxeta)    continue;    //eta cuts on candidates
        
        // RecKine: ITSTPC cuts
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
        // Reject kink mother
        Bool_t kinkmotherpass = kTRUE;
        for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
            if(track->GetID() == listofmotherkink[kinkmother]) {
                kinkmotherpass = kFALSE;
                continue;
            }
        }
        if(!kinkmotherpass) continue;
        
        // RecPrim
        //  if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;  //deleted for DCA absence
        // HFEcuts: ITS layers cuts
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        // HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        //==========================================================================================================
        Double_t eta = track->Eta();
        Double_t phi = track->Phi();
        
        if(fPhiCut){
            if(phi<1.4 || phi >3.14)continue; //to have same EMCal phi acceptance
        }
        
        
        
        Double_t pt = track->Pt();         //pt track after cuts
        if(pt<fpTCutmin || pt>fpTCutmax) continue;
        //==========================================================================================================
        //=========================================PID==============================================================
        if(track->GetTPCsignalN() < fTPCS) continue;
        Float_t fITSnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasITS(track, AliPID::kElectron) : 1000;
        Float_t fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
        Float_t fTOFnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTOF(track, AliPID::kElectron) : 1000;
        //   Float_t eDEDX = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, kTRUE);
        
        Float_t fTPCnSigmaPI = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kPion) : 1000;
        
        
        Double_t CorrectTPCNSigma;
        Double_t mult = fVevent->GetNumberOfESDTracks()/8;
        
        if(multCorrection){
            CorrectTPCNSigma = ftpcpid->GetCorrectedTPCnSigma(eta, mult, fTPCnSigma);
            // cout <<fTPCnSigma << "   ====  " <<COrrectTPCNSigma<<endl;
            fTPCnSigma = CorrectTPCNSigma;
            // cout <<fTPCnSigma << "   ====  " <<COrrectTPCNSigma<<endl;
        }
        
        //       if(TMath::Abs(fTOFnSigma) < fmaxTOFnSigma && TMath::Abs(fTPCnSigma) < fmaxTPCnsigmaLowpT && fTPCnSigmaPI > 4)
        // {
        //     fITSnsigmaElect->Fill(track->P(),fITSnSigma);
        // }
        
        
        fITSnsigma->Fill(track->P(),fITSnSigma);
        fTPCnsigma->Fill(track->P(),fTPCnSigma);
        fTOFns->Fill(track->P(),fTOFnSigma);
        fITSvsTOF->Fill(fTOFnSigma,fITSnSigma);
        fTPCvsITS->Fill(fTPCnSigma,fITSnSigma);
        fTPCvsTOF->Fill(fTPCnSigma,fTOFnSigma);
        
        // if( pt >= 0.3){
        if(fTOFnSigma < fminTOFnSigma || fTOFnSigma > fmaxTOFnSigma) continue;
        // }//cuts on nsigma tof full pt range
        
        fITSnsigmaAftTOF->Fill(track->P(),fITSnSigma);
        fTPCnsigmaAftTOF->Fill(track->P(),fTPCnSigma);
        fTPCvsITSafterTOF->Fill(fTPCnSigma,fITSnSigma);
        
        Double_t valPidSparse[3] = {
            pt,
            fITSnSigma,
            fTPCnSigma,
        };
        fQAPidSparse->Fill(valPidSparse);
        
        
        if( pt < 1.5){
            if(fITSnSigma < fminITSnsigmaLowpT || fITSnSigma > fmaxITSnsigmaLowpT)continue;
        }//cuts on nsigma its low pt
        if( pt >= 1.5){
            if(fITSnSigma < fminITSnsigmaHighpT || fITSnSigma > fmaxITSnsigmaHighpT)continue;
        }//cuts on nsigma its high pt
        fTPCnsigmaAftITSTOF->Fill(track->P(),fTPCnSigma);
        if(pt < 1.5){
            if(fTPCnSigma < fminTPCnsigmaLowpT || fTPCnSigma > fmaxTPCnsigmaLowpT) continue;
        }//cuts on nsigma tpc lowpt
        if(pt >= 1.5){
            if(fTPCnSigma < fminTPCnsigmaHighpT || fTPCnSigma > fmaxTPCnsigmaHighpT) continue;
        }//cuts on nsigma tpc high pt
        fTPCnsigmaAft->Fill(track->P(),fTPCnSigma);
        fTPCnsigmaVSptAft->Fill(pt,fTPCnSigma);
        
        //==========================================================================================================
        //=========================================QA PID SPARSE====================================================
        Float_t timeTOF = track->GetTOFsignal();
        Double_t intTime[5] = {-99., -99., -99., -99., -99.};
        track->GetIntegratedTimes(intTime);
        Float_t timeElec = intTime[0];
        Float_t intLength = 2.99792458e-2* timeElec;
        Double_t beta = 0.1;
        if ((intLength > 0) && (timeTOF > 0))
            beta = intLength/2.99792458e-2/timeTOF;
        
        //     if(fQAPIDSparse){
        //         Double_t valPid[4] = {
        //             track->P(),
        //             track->GetTPCsignal(),
        //             beta,
        //             track->Charge()
        //         };
        //         fQAPid->Fill(valPid);
        //     }
        
        
        fITSnsigmaAft->Fill(track->P(),fITSnSigma);
        fTPCnsigmaAft->Fill(track->P(),fTPCnSigma);
        fTOFnsAft->Fill(track->P(),fTOFnSigma);
        fTOFBetaAft->Fill(track->P(),beta);
        fInclusiveElecPt->Fill(pt);
        fPhi->Fill(phi);
        fEta->Fill(eta);
        //=========================================================================================================
        //----------------------Flow of Inclusive Electrons--------------------------------------------------------
        AliFlowTrack *sTrack = new AliFlowTrack();
        sTrack->Set(track);
        sTrack->SetID(track->GetID());
        sTrack->SetForRPSelection(kTRUE);
        sTrack->SetForPOISelection(kTRUE);
        sTrack->SetMass(263732);
        for(int iRPs=0; iRPs!=fFlowEvent->NumberOfTracks(); ++iRPs)
        {
            //   cout << " no of rps " << iRPs << endl;
            AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack( iRPs ));
            if (!iRP) continue;
            if (!iRP->InRPSelection()) continue;
            if( sTrack->GetID() == iRP->GetID())
            {
                if(fDebug) printf(" was in RP set");
                //  cout << sTrack->GetID() <<"   ==  " << iRP->GetID() << " was in RP set====REMOVED" <<endl;
                iRP->SetForRPSelection(kFALSE);
                // fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
            }
        } //end of for loop on RPs
        fFlowEvent->InsertTrack(((AliFlowTrack*) sTrack));
        fFlowEvent->SetNumberOfPOIs(fFlowEvent->GetNumberOfPOIs()+1);
        //============================Event Plane Method with V0====================================================
        Double_t v2PhiV0A = TMath::Cos(2*(phi - evPlAngV0A));
        Double_t v2PhiV0C = TMath::Cos(2*(phi - evPlAngV0C));
        Double_t v2Phi[3] = {
            v2PhiV0A,
            v2PhiV0C,
            pt};
        fV2Phi->Fill(v2Phi);
        
        Double_t v2PhiVz = TMath::Cos(2*(phi - evPlAngV0));
        Double_t v2PhiV0tot[2] = {
            v2PhiVz,
            pt};
        
        if(EPweights) fV2Phivzerotot->Fill(v2PhiV0tot,weightEP);
        if(!EPweights) fV2Phivzerotot->Fill(v2PhiV0tot);
        
        
        //==========================================================================================================
        //=========================================================================================================
        
        
        if(fDCA){
            fNonHFE = new AliSelectNonHFE();
            fNonHFE->SetAODanalysis(kTRUE);
            fNonHFE->SetInvariantMassCut(fInvmassCut);
            if(fOP_angle) fNonHFE->SetOpeningAngleCut(fOpeningAngleCut);
            //fNonHFE->SetChi2OverNDFCut(fChi2Cut);
            //if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
            fNonHFE->SetAlgorithm("DCA"); //KF
            fNonHFE->SetPIDresponse(pidResponse);
            fNonHFE->SetTrackCuts(-3,3);
            
            fNonHFE->SetHistAngleBack(fOpeningAngleLS);
            fNonHFE->SetHistAngle(fOpeningAngleULS);
            //fNonHFE->SetHistDCABack(fDCABack);
            //fNonHFE->SetHistDCA(fDCA);
            fNonHFE->SetHistMassBack(fInvmassLS1);
            fNonHFE->SetHistMass(fInvmassULS1);
            
            fNonHFE->FindNonHFE(iTracks,track,fAOD);
            
            // Int_t *fUlsPartner = fNonHFE->GetPartnersULS();
            // Int_t *fLsPartner = fNonHFE->GetPartnersLS();
            // Bool_t fUlsIsPartner = kFALSE;
            // Bool_t fLsIsPartner = kFALSE;
            if(fNonHFE->IsULS()){
                for(Int_t kULS =0; kULS < fNonHFE->GetNULS(); kULS++){
                    fULSElecPt->Fill(track->Pt());
                }
            }
            
            if(fNonHFE->IsLS()){
                for(Int_t kLS =0; kLS < fNonHFE->GetNLS(); kLS++){
                    fLSElecPt->Fill(track->Pt());
                }
            }
        }
        if(!fDCA){
            //=========================================================================================================
            //----------------------Selection and Flow of Photonic Electrons-----------------------------
            Bool_t fFlagPhotonicElec = kFALSE;
            SelectPhotonicElectron(iTracks,track,fTPCnSigma,evPlAngV0,fFlagPhotonicElec,weightEP,mult);
            if(fFlagPhotonicElec){
                fPhotoElecPt->Fill(pt);
                fITSnsigmaElect->Fill(track->P(),fITSnSigma);
            }
            // Semi inclusive electron
            if(!fFlagPhotonicElec){fSemiInclElecPt->Fill(pt);}
        }
        //-------------------------------------------------------------------------------------------
        
    }//end loop on track
    PostData(1, fOutputList);
    PostData(2, fFlowEvent);
    
    //----------hfe end---------
}
//_________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SelectPhotonicElectron(Int_t itrack,const AliAODTrack *track,Float_t fTPCnSigma,Double_t evPlAngV0, Bool_t &fFlagPhotonicElec, Double_t weightEPflat, Double_t multev)
{
    
    //Identify non-heavy flavour electrons using Invariant mass method
    Bool_t flagPhotonicElec = kFALSE;
    
    for(Int_t jTracks = 0; jTracks<fAOD->GetNumberOfTracks(); jTracks++){
        AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(jTracks));
        if(!trackAsso) AliFatal("Not a standard AOD");
        if (!trackAsso) {
            printf("ERROR: Could not receive track %d\n", jTracks);
            continue;
        }
        //  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;  // TESTBIT FOR AOD double Counting
        if(!trackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        
        if(fAssoITSRefit){
            if(!(trackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
        }
        
        if(!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
        
        //	if((!(trackAsso->GetStatus()&AliESDtrack::kITSrefit)|| (!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
        
        if(jTracks == itrack) continue;
        Double_t ptAsso=-999., nsigma=-999.0;
        Double_t mass=-999., width = -999;
        Double_t openingAngle = -999.;
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        
        
        ptAsso = trackAsso->Pt();
        Short_t chargeAsso = trackAsso->Charge();
        Short_t charge = track->Charge();
        // nsigma = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(ptAsso < fptminAsso) continue;
        
        nsigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
        
        
        Double_t CorrectTPCNSigma;
        if(multCorrection){
            CorrectTPCNSigma = ftpcpid->GetCorrectedTPCnSigma(trackAsso->Eta(), multev, nsigma);
            nsigma = CorrectTPCNSigma;
        }
        
        
        
        //if(trackAsso->GetTPCNcls() < 80) continue;
        if(trackAsso->GetTPCNcls() < fAssoTPCCluster) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fAOD->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        recg.GetMass(mass,width);
        
        openingAngle = ge1.GetAngle(ge2);
        if(fFlagLS) fOpeningAngleLS->Fill(openingAngle);
        if(fFlagULS) fOpeningAngleULS->Fill(openingAngle);
        if(fOP_angle){if(openingAngle > fOpeningAngleCut) continue;}
        
        
        if(fFlagLS) fInvmassLS1->Fill(mass);
        if(fFlagULS) fInvmassULS1->Fill(mass);
        
        if(fFlagULS){
            Double_t MassSparseULS[3] = {
                track->Pt(),
                mass
            };
            fSparseMassULS->Fill(MassSparseULS);
        }
        if(fFlagLS){
            Double_t MassSparseLS[3] = {
                track->Pt(),
                mass
            };
            fSparseMassLS->Fill(MassSparseLS);
        }
        
        
        if(mass<fInvmassCut){
            if(fFlagULS){fULSElecPt->Fill(track->Pt());}
            if(fFlagLS){fLSElecPt->Fill(track->Pt());}
        }
        
        
        Double_t phi = track->Phi();
        Float_t DeltaPhi_eEP = TVector2::Phi_0_2pi(phi - evPlAngV0);
        if(DeltaPhi_eEP > TMath::Pi()) {DeltaPhi_eEP = DeltaPhi_eEP - TMath::Pi();}
        
        
        if(mass<fInvmassCut){
            if(fFlagULS){
                Double_t ulsSparse[3] = {
                    track->Pt(),
                    fTPCnSigma,
                    DeltaPhi_eEP
                };
                if(EPweights) fSparsephipsiULS->Fill(ulsSparse,weightEPflat);
                if(!EPweights) fSparsephipsiULS->Fill(ulsSparse);
            }
            if(fFlagLS){
                Double_t lsSparse[3] = {
                    track->Pt(),
                    fTPCnSigma,
                    DeltaPhi_eEP
                };
                if(EPweights) fSparsephipsiLS->Fill(lsSparse,weightEPflat);
                if(!EPweights)fSparsephipsiLS->Fill(lsSparse);
            }
        }
        
        
        
        if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
            flagPhotonicElec = kTRUE;
        }
    }//track loop
    
    fFlagPhotonicElec = flagPhotonicElec;
}
//___________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::UserCreateOutputObjects()
{
    
    //Create histograms
    //----------hfe initialising begin---------
    fNullCuts = new AliFlowTrackCuts("null_cuts");
    
    AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
    cc->SetNbinsMult(10000);
    cc->SetMultMin(0);
    cc->SetMultMax(10000);
    
    cc->SetNbinsPt(100);
    cc->SetPtMin(0);
    cc->SetPtMax(5);
    
    cc->SetNbinsPhi(180);
    cc->SetPhiMin(0.0);
    cc->SetPhiMax(TMath::TwoPi());
    
    cc->SetNbinsEta(30);
    cc->SetEtaMin(-8.0);
    cc->SetEtaMax(+8.0);
    
    cc->SetNbinsQ(500);
    cc->SetQMin(0.0);
    cc->SetQMax(3.0);
    
    
    //   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    //   AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    //    if (!inputHandler){
    //        AliFatal("Input handler needed");
    //   }
    //   else{
    //       fPIDResponse=inputHandler->GetPIDResponse();
    //   }
    //pid response object
    //  if (!fPIDResponse)AliError("PIDResponse object was not created");
    
    
    //--------Initialize PID
    fPID->SetHasMCData(kFALSE);
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("ITS", 0);
        fPID->AddDetector("TOF", 1);
        fPID->AddDetector("TPC", 2);
        
    }
    
    fPID->SortDetectors();
    fPIDqa = new AliHFEpidQAmanager();
    fPIDqa->Initialize(fPID);
    
    
    
    //--------Initialize correction Framework and Cuts
    fCFM = new AliCFManager;
    const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
    fCFM->SetNStepParticle(kNcutSteps);
    for(Int_t istep = 0; istep < kNcutSteps; istep++)
        fCFM->SetParticleCutsList(istep, NULL);
    
    if(!fCuts){
        AliWarning("Cuts not available. Default cuts will be used");
        fCuts = new AliHFEcuts;
        fCuts->CreateStandardCuts();
    }
    
    fCuts->SetAOD();
    fCuts->Initialize(fCFM);
    //----------hfe initialising end--------
    //---------Output Tlist
    fOutputList = new TList();
    fOutputList->SetOwner();
    fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    fNoEvents = new TH1F("fNoEvents","",1,0,1) ;
    fOutputList->Add(fNoEvents);
    
    fITSnsigma = new TH2F("fITSnsigma", "ITS - n sigma before HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fITSnsigma);
    
    
    fITSnsigmaElect = new TH2F("fITSnsigmaElect", "fITSnsigmaElect - n sigma before HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fITSnsigmaElect);
    
    
    fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma before HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigma);
    
    fITSnsigmaAft = new TH2F("fITSnsigmaAft", "ITS - n sigma after HFE pid",1000,0,10,300,-10,20);
    fOutputList->Add(fITSnsigmaAft);
    fITSvsTOF = new TH2F("fITSvsTOF", "ITS tof",400,-20,20,400,-20,20);
    fOutputList->Add(fITSvsTOF);
    fTPCvsITS = new TH2F("TPCvsITS", "TPC ITS",400,-20,20,400,-20,20);
    fOutputList->Add(fTPCvsITS);
    fTPCvsTOF = new TH2F("TPCvsTOF", "TPC TOF",400,-20,20,400,-20,20);
    fOutputList->Add(fTPCvsTOF);
    fTPCvsITSafterTOF = new TH2F("TPCvsITSafterTOF", "TPC ITS",400,-20,20,400,-20,20);
    fOutputList->Add(fTPCvsITSafterTOF);
    
    
    fITSnsigmaAftTOF = new TH2F("fITSnsigmaAftTOF", "ITS - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fITSnsigmaAftTOF);
    
    fTPCnsigmaAft = new TH2F("fTPCnsigmaAft", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaAft);
    
    fTPCnsigmaVSptAft = new TH2F("fTPCnsigmaVSptAft", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaVSptAft);
    
    fTPCnsigmaAftITSTOF = new TH2F("fTPCnsigmaAftITSTOF", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaAftITSTOF);
    
    fTPCnsigmaAftTOF = new TH2F("fTPCnsigmaAftTOF", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaAftTOF);
    
    fTOFns = new TH2F("fTOFns","track TOFnSigma",600,0,6,400,-20,20);
    fOutputList->Add(fTOFns);
    
    fTOFnsAft = new TH2F("fTOFnsAft","track TOFnSigma",600,0,6,400,-20,20);
    fOutputList->Add(fTOFnsAft);
    
    fTOFBetaAft = new TH2F("fTOFBetaAft","track TOFBeta",600,0,6,120,0,1.2);
    fOutputList->Add(fTOFBetaAft);
    
    fInclusiveElecPt = new TH1F("fInclElecPt", "Inclusive electron pt",100,0,5);
    fOutputList->Add(fInclusiveElecPt);
    
    fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",100,0,5);
    fOutputList->Add(fPhotoElecPt);
    
    fSemiInclElecPt = new TH1F("fSemiInclElecPt", "Semi-inclusive electron pt",100,0,5);
    fOutputList->Add(fSemiInclElecPt);
    
    fULSElecPt = new TH1F("fULSElecPt", "ULS electron pt",100,0,5);
    fOutputList->Add(fULSElecPt);
    
    fLSElecPt = new TH1F("fLSElecPt", "LS electron pt",100,0,5);
    fOutputList->Add(fLSElecPt);
    
    fInvmassLS1 = new TH1F("fInvmassLS1", "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassLS1);
    
    fInvmassULS1 = new TH1F("fInvmassULS1", "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassULS1);
    
    fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass", 101, -1, 100);
    fOutputList->Add(fCentralityPass);
    
    fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass", 101, -1, 100);
    fOutputList->Add(fCentralityNoPass);
    
    fCentralityNoPassForFlattening = new TH1F("fCentralityNoPassForFlattening", "Centrality No Pass for flattening", 101, -1, 100);
    fOutputList->Add(fCentralityNoPassForFlattening);
    
    fCentralityBeforePileup = new TH1F("fCentralityBeforePileup", "fCentralityBeforePileup Pass", 101, -1, 100);
    fOutputList->Add(fCentralityBeforePileup);
    
    fCentralityAfterVZTRK = new TH1F("fCentralityAfterVZTRK", "fCentralityAfterVZTRK Pass", 101, -1, 100);
    fOutputList->Add(fCentralityAfterVZTRK);
    
    fCentralityAfterCorrCut = new TH1F("fCentralityAfterCorrCut", "fCentralityAfterCorrCut Pass", 101, -1, 100);
    fOutputList->Add(fCentralityAfterCorrCut);
    
    fPhi = new TH1F("fPhi", "#phi distribution", 100, -.5, 7);
    fOutputList->Add(fPhi);
    
    fEta = new TH1F("fEta", "#eta distribution", 100, -1.1, 1.1);
    fOutputList->Add(fEta);
    
    fVZEROA = new TH1F("fVZEROA", "VZERO A Multiplicity", 1000, 0, 10000);
    fOutputList->Add(fVZEROA);
    
    fVZEROC = new TH1F("fVZEROC", "VZERO C Multiplicity", 1000, 0, 10000);
    fOutputList->Add(fVZEROC);
    
    fTPCM = new TH1F("fTPCM", "TPC multiplicity", 1000, 0, 10000);
    fOutputList->Add(fTPCM);
    
    fvertex = new TH1D("fvertex", "vertex distribution", 300, -15,15);
    fOutputList->Add(fvertex);
    
    fMultCorBeforeCuts = new TH2F("fMultCorBeforeCuts", "TPC vs Global multiplicity (Before cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fOutputList->Add(fMultCorBeforeCuts);
    
    fMultCorAfterCuts = new TH2F("fMultCorAfterCuts", "TPC vs Global multiplicity (After cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fOutputList->Add(fMultCorAfterCuts);
    
    fMultCorAfterCentrBeforeCuts = new TH2F("fMultCorAfterCentrBeforeCuts", "TPC vs Global multiplicity (After CC before cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fOutputList->Add(fMultCorAfterCentrBeforeCuts);
    
    fMultCorAfterVZTRKComp = new TH2F("fMultCorAfterVZTRKComp", "TPC vs Global multiplicity (After V0-TRK); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fOutputList->Add(fMultCorAfterVZTRKComp);
    
    fMultCorAfterCorrCut = new TH2F("fMultCorAfterCorrCut", "TPC vs Global multiplicity (After CorrCut); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
    fOutputList->Add(fMultCorAfterCorrCut);
    
    fMultvsCentr = new TH2F("fMultvsCentr", "Multiplicity vs centrality; centrality; Multiplicity", 100, 0., 100, 100, 0, 3000);
    fOutputList->Add(fMultvsCentr);
    
    fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
    fOutputList->Add(fOpeningAngleLS);
    
    fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
    fOutputList->Add(fOpeningAngleULS);
    
    
    
    //----------------------------------------------------------------------------
    EPVzA = new TH1D("EPVzA", "EPVzA", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPVzA);
    EPVzC = new TH1D("EPVzC", "EPVzC", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPVzC);
    EPTPC = new TH1D("EPTPC", "EPTPC", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPTPC);
    EPVz = new TH1D("EPVz", "EPVz", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPVz);
    EPTPCp = new TH1D("EPTPCp", "EPTPCp", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPTPCp);
    EPTPCn = new TH1D("EPTPCn", "EPTPCn", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPTPCn);
    
    
    //----------------------------------------------------------------------------
    fSubEventDPhiv2 = new TProfile("fSubEventDPhiv2", "fSubEventDPhiv2", 3, 0, 3);
    fSubEventDPhiv2->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
    fSubEventDPhiv2->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
    fSubEventDPhiv2->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
    fSubEventDPhiv2->Sumw2();
    fOutputList->Add(fSubEventDPhiv2);
    
    fSubEventDPhiv2new = new TProfile("fSubEventDPhiv2new", "fSubEventDPhiv2new", 3, 0, 3);
    fSubEventDPhiv2new->GetXaxis()->SetBinLabel(1, "<cos(2(#Psi_{a} - #Psi_{b}))>");
    fSubEventDPhiv2new->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{a} - #Psi_{c}>))");
    fSubEventDPhiv2new->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{b} - #Psi_{c}>))");
    fSubEventDPhiv2new->Sumw2();
    fOutputList->Add(fSubEventDPhiv2new);
    
    //================================Event Plane with VZERO A and C=====================
    const Int_t nPtBins = 12;
    Double_t binsPt[nPtBins+1] = {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5};
    // v2A, v2C, pt
    Int_t    bins[3] = {  50,  50, nPtBins};
    Double_t xmin[3] = { -1., -1.,   0};
    Double_t xmax[3] = {  1.,  1.,   5};
    fV2Phi = new THnSparseF("fV2Phi", "v2A:v2C:pt", 3, bins, xmin, xmax);
    // Set bin limits for axes which are not standard binned
    fV2Phi->SetBinEdges(2, binsPt);
    // set axes titles
    fV2Phi->GetAxis(0)->SetTitle("v_{2} (V0A)");
    fV2Phi->GetAxis(1)->SetTitle("v_{2} (V0C)");
    fV2Phi->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
    fV2Phi->Sumw2();
    fOutputList->Add(fV2Phi);
    
    
    //================================Event Plane with VZERO=====================
    Int_t    binsV[2] = {  50,  100};
    Double_t xminV[2] = { -1.,   0};
    Double_t xmaxV[2] = {  1.,   5};
    fV2Phivzerotot = new THnSparseF("fV2Phivzerotot", "v2:pt", 2, binsV, xminV, xmaxV);
    // Set bin limits for axes which are not standard binned
    //fV2Phivzerotot->SetBinEdges(1, binsPt);
    // set axes titles
    fV2Phivzerotot->GetAxis(0)->SetTitle("v_{2} (V0)");
    fV2Phivzerotot->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    fV2Phivzerotot->Sumw2();
    fOutputList->Add(fV2Phivzerotot);
    
    
    
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //   if(fQAPIDSparse){
    //       Int_t    binsQA[4] = { 150,  100,  120,    3};
    //       Double_t xminQA[4] = { 0.,    50,   0, -1.5};
    //       Double_t xmaxQA[4] = { 15.,  150, 1.2,  1.5};
    //       fQAPid = new THnSparseF("fQAPid", "p:dEdx:beta:ch", 4, binsQA, xminQA, xmaxQA);
    //       fQAPid->GetAxis(0)->SetTitle("p (Gev/c");
    //       fQAPid->GetAxis(1)->SetTitle("dE/dx");
    //        fQAPid->GetAxis(2)->SetTitle("#beta (TOF)");
    //        fQAPid->GetAxis(3)->SetTitle("charge");
    //        fOutputList->Add(fQAPid);
    //    }
    //===========================================================================
    Int_t    binsQA2[3] = { 100,  40, 150/*,  60*/};
    Double_t xminQA2[3] = { 0.,   -2, -15/*,  -3*/};
    Double_t xmaxQA2[3] = { 5.,    2,  15/*,   3*/};
    fQAPidSparse = new THnSparseF("fQAPidSparse", "pt:itsnsigma:tpcnsigma", 3, binsQA2, xminQA2, xmaxQA2);
    fQAPidSparse->GetAxis(0)->SetTitle("pt (Gev/c)");
    fQAPidSparse->GetAxis(1)->SetTitle("itsnsigma");
    fQAPidSparse->GetAxis(2)->SetTitle("tpcnsigma");
    fOutputList->Add(fQAPidSparse);
    //===========================================================================
    
    
    
    Int_t    binsphipsi[3] = { 100, 150,           6};
    Double_t xminphipsi[3] = { 0.,  -15,           0};
    Double_t xmaxphipsi[3] = { 5.,   15, TMath::Pi()};
    fSparsephipsiULS = new THnSparseF("fSparsephipsiULS", "pt:tpcnsigma:DeltaPhiULS", 3, binsphipsi, xminphipsi, xmaxphipsi);
    fSparsephipsiULS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparsephipsiULS->GetAxis(1)->SetTitle("tpcnsigma");
    fSparsephipsiULS->GetAxis(2)->SetTitle("DeltaPhiULS");
    fSparsephipsiULS->Sumw2();
    fOutputList->Add(fSparsephipsiULS);
    
    fSparsephipsiLS = new THnSparseF("fSparsephipsiLS", "pt:tpcnsigma:DeltaPhiLS", 3, binsphipsi, xminphipsi, xmaxphipsi);
    fSparsephipsiLS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparsephipsiLS->GetAxis(1)->SetTitle("tpcnsigma");
    fSparsephipsiLS->GetAxis(2)->SetTitle("DeltaPhiLS");
    fSparsephipsiLS->Sumw2();
    fOutputList->Add(fSparsephipsiLS);
    
    
    Int_t    binsmass[2] = { 100, 200};
    Double_t xminmass[2] = { 0.,  0};
    Double_t xmaxmass[2] = { 5., 1.};
    fSparseMassULS = new THnSparseF("fSparseMassULS", "pt:mass (GeV/c^{2})", 2, binsmass, xminmass, xmaxmass);
    fSparseMassULS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparseMassULS->GetAxis(1)->SetTitle("mass");
    fOutputList->Add(fSparseMassULS);
    
    fSparseMassLS = new THnSparseF("fSparseMassLS", "pt:mass (GeV/c^{2})", 2, binsmass, xminmass, xmaxmass);
    fSparseMassLS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparseMassLS->GetAxis(1)->SetTitle("mass");
    fOutputList->Add(fSparseMassLS);
    
    
    EPVzAftW = new TH1D("EPVzAftW", "EPVzAftW",60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputList->Add(EPVzAftW);
    
    
    fOutputList->Add(fHistEPDistrWeight);
    
    fCentralityAll = new TH1F("fCentralityAll", "Centrality Pass All", 101, -1, 100);
    fOutputList->Add(fCentralityAll);
    
    
    
    PostData(1,fOutputList);
    // create and post flowevent
    fFlowEvent = new AliFlowEvent(10000);
    PostData(2, fFlowEvent);
    
}
//________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::Terminate(Option_t *)
{
    // Info("Terminate");
    AliAnalysisTaskSE::Terminate();
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskFlowITSTPCTOFQCSP::PlotVZeroMultiplcities(const T* event) const
{
    // QA multiplicity plots
    fVZEROA->Fill(event->GetVZEROData()->GetMTotV0A());
    fVZEROC->Fill(event->GetVZEROData()->GetMTotV0C());
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskFlowITSTPCTOFQCSP::SetNullCuts(T* event)
{
    //Set null cuts
    if (fDebug) cout << " fCutsRP " << fCutsRP << endl;
    fCutsRP->SetEvent(event, MCEvent());
    fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
    fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
    fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
    fNullCuts->SetEvent(event, MCEvent());
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const
{
    //Prepare flow events
    FlowEv->ClearFast();
    FlowEv->Fill(fCutsRP, fNullCuts);
    FlowEv->SetReferenceMultiplicity(iMulti);
    FlowEv->DefineDeadZone(0, 0, 0, 0);
    //  FlowEv->TagSubeventsInEta(-0.7, 0, 0, 0.7);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowITSTPCTOFQCSP::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    // Check single track cuts for a given cut step
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}
//_________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::CheckCentrality(AliAODEvent* event, Bool_t &centralitypass)
{
    //============================Multiplicity TPV vs Global===============================================================================
    const Int_t nGoodTracks = event->GetNumberOfTracks();
    Float_t multTPC(0.); // tpc mult estimate
    Float_t multGlob(0.); // global multiplicity
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill tpc mult
        AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(iTracks));
        if(!trackAOD) AliFatal("Not a standard AOD");
        if (!trackAOD) continue;
        if (!(trackAOD->TestFilterBit(1))) continue;
        if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2)) continue;
        multTPC++;
    }
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill global mult
        AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(iTracks));
        if(!trackAOD) AliFatal("Not a standard AOD");
        if (!trackAOD) continue;
        if (!(trackAOD->TestFilterBit(16))) continue;
        if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1)) continue;
        Double_t b[2] = {-99., -99.};
        Double_t bCov[3] = {-99., -99., -99.};
        if (!(trackAOD->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov))) continue;
        if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
        multGlob++;
    } //track loop
    fMultCorBeforeCuts->Fill(multGlob, multTPC);//before all cuts...even before centrality selectrion
    //============================================================================================================================
    // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
    if (!fkCentralityMethod) AliFatal("No centrality method set! FATAL ERROR!");
    fCentrality = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethod);
    fCentralityAll->Fill(fCentrality);
    
    //   cout << "--------------Centrality evaluated-------------------------"<<endl;
    if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
    {
        fCentralityNoPass->Fill(fCentrality);
        //    cout << "--------------Fill no pass-----"<< fCentrality <<"--------------------"<<endl;
        centralitypass = kFALSE;
    }else
    {
        //    cout << "--------------Fill pass----"<< fCentrality <<"---------------------"<<endl;
        centralitypass = kTRUE;
    }
    if (centralitypass){
        fMultCorAfterCentrBeforeCuts->Fill(multGlob, multTPC);
        fCentralityBeforePileup->Fill(fCentrality);
    }//...after centrality selectrion
    //============================================================================================================================
    //to remove the bias introduced by multeplicity outliers---------------------
    Float_t centTrk = event->GetCentrality()->GetCentralityPercentile("TRK");
    Float_t centv0 = event->GetCentrality()->GetCentralityPercentile("V0M");
    if (TMath::Abs(centv0 - centTrk) > 5.0){
        centralitypass = kFALSE;
        fCentralityNoPass->Fill(fCentrality);
    }
    if (centralitypass){
        fMultCorAfterVZTRKComp->Fill(multGlob, multTPC);
        fCentralityAfterVZTRK->Fill(fCentrality);
    }//...after centrality selectrion
    //============================================================================================================================
    if(fMultCut){
        if(fTrigger==1 || fTrigger==4 || fTrigger==5){
            if(! (multTPC > (-36.73 + 1.48*multGlob) && multTPC < (62.87 + 1.78*multGlob))){
                //   cout <<" Trigger ==" <<fTrigger<< endl;
                centralitypass = kFALSE;
                fCentralityNoPass->Fill(fCentrality);
            }//2011 Semicentral
        }
        if(fTrigger==0){
            if(! (multTPC > (77.9 + 1.395*multGlob) && multTPC < (187.3 + 1.665*multGlob))){
                //     cout <<" Trigger ==" <<fTrigger<< endl;
                centralitypass = kFALSE;
                fCentralityNoPass->Fill(fCentrality);
            }//2011
        }//2011 Central
    }
    if (centralitypass){
        fMultCorAfterCorrCut->Fill(multGlob, multTPC);
        fCentralityAfterCorrCut->Fill(fCentrality);
    }//...after CORR CUT
    //=================================All cuts are passed==================++++==================================================
    //=================================Now Centrality flattening for central trigger==================++++==================================================
    if(!fCentFlatMine){
        if(fTrigger==0 || fTrigger==4){
            if(!IsEventSelectedForCentrFlattening(fCentrality)){
                centralitypass = kFALSE;
                fCentralityNoPassForFlattening->Fill(fCentrality);
            }
        }
    }
    if(fCentFlatMine){
        if(fTrigger==0 || fTrigger==4){
            if(!IsEventSelectedForCentrFlattening_Bis(fCentrality)){
                centralitypass = kFALSE;
                fCentralityNoPassForFlattening->Fill(fCentrality);
            }
        }
    }
    //==============================fill histo after all cuts==============================++++==================================================
    if(centralitypass){
        fCentralityPass->Fill(fCentrality);
        fMultCorAfterCuts->Fill(multGlob, multTPC);
        fMultvsCentr->Fill(fCentrality, multTPC);
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod)
{
    // Set a centrality range ]min, max] and define the method to use for centrality selection
    fCentralityMin = CentralityMin;
    fCentralityMax = CentralityMax;
    fkCentralityMethod = CentralityMethod;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetIDCuts(Double_t minTOFnSigma, Double_t maxTOFnSigma, Double_t minITSnsigmalowpt, Double_t maxITSnsigmalowpt,  Double_t minITSnsigmahighpt, Double_t maxITSnsigmahighpt, Double_t minTPCnsigmalowpt, Double_t maxTPCnsigmalowpt,  Double_t minTPCnsigmahighpt, Double_t maxTPCnsigmahighpt)
{
    //Set PID cuts
    fminTOFnSigma = minTOFnSigma;
    fmaxTOFnSigma = maxTOFnSigma;
    fminITSnsigmaLowpT = minITSnsigmalowpt;
    fmaxITSnsigmaLowpT = maxITSnsigmalowpt;
    fminITSnsigmaHighpT = minITSnsigmahighpt;
    fmaxITSnsigmaHighpT = maxITSnsigmahighpt;
    fminTPCnsigmaLowpT = minTPCnsigmalowpt;
    fmaxTPCnsigmaLowpT = maxTPCnsigmalowpt;
    fminTPCnsigmaHighpT = minTPCnsigmahighpt;
    fmaxTPCnsigmaHighpT = maxTPCnsigmahighpt;
    
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetpTCuttrack(Double_t ptmin, Double_t ptmax)
{
    //Set pt cuts
    fpTCutmin = ptmin;
    fpTCutmax = ptmax;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetHistoForCentralityFlattening(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef,Int_t switchTRand){
    // set the histo for centrality flattening
    // the centrality is flatten in the range minCentr,maxCentr
    // if centrRef is zero, the minimum in h within (minCentr,maxCentr) defines the reference
    //                positive, the value of h(centrRef) defines the reference (-> the centrality distribution might be not flat in the whole desired range)
    //                negative, h(bin with max in range)*centrRef is used to define the reference (-> defines the maximum loss of events, also in this case the distribution might be not flat)
    // switchTRand is used to set the unerflow bin of the histo: if it is < -1 in the analysis the random event selection will be done on using TRandom
    
    if(maxCentr<minCentr){
        AliWarning("AliAnalysisCheckCorrdist::Wrong centralities values while setting the histogram for centrality flattening");
    }
    
    if(fHistCentrDistr)delete fHistCentrDistr;
    fHistCentrDistr=(TH1F*)h->Clone("hCentralityFlat");
    fHistCentrDistr->SetTitle("Reference histo for centrality flattening");
    Int_t minbin=fHistCentrDistr->FindBin(minCentr*1.00001); // fast if fix bin width
    Int_t maxbin=fHistCentrDistr->FindBin(maxCentr*0.9999);
    fHistCentrDistr->GetXaxis()->SetRange(minbin,maxbin);
    Double_t ref=0.,bincont=0.,binrefwidth=1.;
    Int_t binref=0;
    if(TMath::Abs(centrRef)<0.0001){
        binref=fHistCentrDistr->GetMinimumBin();
        binrefwidth=fHistCentrDistr->GetBinWidth(binref);
        ref=fHistCentrDistr->GetBinContent(binref)/binrefwidth;
    }
    else if(centrRef>0.){
        binref=h->FindBin(centrRef);
        if(binref<1||binref>h->GetNbinsX()){
            AliWarning("AliRDHFCuts::Wrong centrality reference value while setting the histogram for centrality flattening");
        }
        binrefwidth=fHistCentrDistr->GetBinWidth(binref);
        ref=fHistCentrDistr->GetBinContent(binref)/binrefwidth;
    }
    else{
        if(centrRef<-1) AliWarning("AliRDHFCuts: with this centrality reference no flattening will be applied");
        binref=fHistCentrDistr->GetMaximumBin();
        binrefwidth=fHistCentrDistr->GetBinWidth(binref);
        ref=fHistCentrDistr->GetMaximum()*TMath::Abs(centrRef)/binrefwidth;
    }
    
    for(Int_t j=1;j<=h->GetNbinsX();j++){// Now set the "probabilities"
        if(h->GetBinLowEdge(j)*1.0001>=minCentr&&h->GetBinLowEdge(j+1)*0.9999<=maxCentr){
            bincont=h->GetBinContent(j);
            fHistCentrDistr->SetBinContent(j,ref/bincont*h->GetBinWidth(j));
            fHistCentrDistr->SetBinError(j,h->GetBinError(j)*ref/bincont);
        }
        else{
            h->SetBinContent(j,1.1);// prob > 1 to assure that events will not be rejected
        }
    }
    
    fHistCentrDistr->SetBinContent(0,switchTRand);
    return;
    
}

//-------------------------------------------------
Bool_t AliAnalysisTaskFlowITSTPCTOFQCSP::IsEventSelectedForCentrFlattening(Float_t centvalue){
    //
    //  Random event selection, based on fHistCentrDistr, to flatten the centrality distribution
    //  Can be faster if it was required that fHistCentrDistr covers
    //  exactly the desired centrality range (e.g. part of the lines below should be done during the
    // setting of the histo) and TH1::SetMinimum called
    //
    
    if(!fHistCentrDistr) return kTRUE;
    // Int_t maxbin=fHistCentrDistr->FindBin(fMaxCentrality*0.9999);
    //   if(maxbin>fHistCentrDistr->GetNbinsX()){
    //     AliWarning("AliRDHFCuts: The maximum centrality exceeds the x-axis limit of the histogram for centrality flattening");
    //   }
    
    Int_t bin=fHistCentrDistr->FindBin(centvalue); // Fast if the histo has a fix bin
    Double_t bincont=fHistCentrDistr->GetBinContent(bin);
    Double_t centDigits=centvalue-(Int_t)(centvalue*100.)/100.;// this is to extract a random number between 0 and 0.01
    
    if(fHistCentrDistr->GetBinContent(0)<-0.9999){
        if(gRandom->Uniform(1.)<bincont)return kTRUE;
        return kFALSE;
    }
    
    if(centDigits*100.<bincont)return kTRUE;
    return kFALSE;
    
}
//---------------------------------------------------------------------------
//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetHistoForEPFlattWeights(TH1D *h){
    
    if(fHistEPDistrWeight)delete fHistEPDistrWeight;
    fHistEPDistrWeight=(TH1D*)h->Clone("fHistEPDistrWeight");
    Double_t Inte = fHistEPDistrWeight->Integral()/fHistEPDistrWeight->GetNbinsX();
    
    
    
    for(Int_t j=1;j<=h->GetNbinsX();j++){// Now set the "probabilities"
        Double_t w = Inte/fHistEPDistrWeight->GetBinContent(j);
        fHistEPDistrWeight->SetBinError(j,0./*h->GetBinError(j)*ref/bincont*/);
        
        fHistEPDistrWeight->SetBinContent(j,w);
    }
    return;
    
}
//-------------------------------------------------
Double_t AliAnalysisTaskFlowITSTPCTOFQCSP::GiveMeWeight(Double_t EP){
    
    Int_t Bin = fHistEPDistrWeight->FindBin(EP);
    Double_t ww = fHistEPDistrWeight->GetBinContent(Bin);
    return ww;
    
}
//-------------------------------------------------


//_____________________________________________________________________________
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetHistoForCentralityFlattening_Bis(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef){
    // set the histo for centrality flattening
    // the centrality is flatten in the range minCentr,maxCentr
    // if centrRef is zero, the minimum in h within (minCentr,maxCentr) defines the reference
    //                positive, the value of h(centrRef) defines the reference (-> the centrality distribution might be not flat in the whole desired range)
    //                negative, h(bin with max in range)*centrRef is used to define the reference (-> defines the maximum loss of events, also in this case the distribution might be not flat)
    // switchTRand is used to set the unerflow bin of the histo: if it is < -1 in the analysis the random event selection will be done on using TRandom
    
    if(maxCentr<minCentr){
        AliWarning("AliAnalysisCheckCorrdist::Wrong centralities values while setting the histogram for centrality flattening");
    }
    
    if(fHistCentrDistr)delete fHistCentrDistr;
    fHistCentrDistr=(TH1F*)h->Clone("fHistCentrDistr");
    Int_t minbin=fHistCentrDistr->FindBin(minCentr*1.00001); // fast if fix bin width
    Int_t maxbin=fHistCentrDistr->FindBin(maxCentr*0.9999);
    fHistCentrDistr->GetXaxis()->SetRange(minbin,maxbin);
    Double_t ref=0.,bincont=0.,binrefwidth=1.;
    Int_t binref=0;
    if(TMath::Abs(centrRef)<0.0001){
        binref=fHistCentrDistr->GetMinimumBin();
        binrefwidth=fHistCentrDistr->GetBinWidth(binref);
        ref=fHistCentrDistr->GetBinContent(binref)/binrefwidth;
    }
    
    for(Int_t j=minbin;j<=maxbin;j++){// Now set the "probabilities"
        bincont=h->GetBinContent(j);
        fHistCentrDistr->SetBinContent(j,ref/bincont*h->GetBinWidth(j));
        fHistCentrDistr->SetBinError(j,0./*h->GetBinError(j)*ref/bincont*/);
    }
    
    return;
    
}

//-------------------------------------------------
Bool_t AliAnalysisTaskFlowITSTPCTOFQCSP::IsEventSelectedForCentrFlattening_Bis(Double_t centvalue){
    //
    //  Random event selection, based on fHistCentrDistr, to flatten the centrality distribution
    //  Can be faster if it was required that fHistCentrDistr covers
    //  exactly the desired centrality range (e.g. part of the lines below should be done during the
    // setting of the histo) and TH1::SetMinimum called
    //
    
    if(!fHistCentrDistr) return kTRUE;
    
    Int_t bin=fHistCentrDistr->FindBin(centvalue); // Fast if the histo has a fix bin
    
    Double_t bincont=fHistCentrDistr->GetBinContent(bin);
    Double_t centDigits=centvalue-(Int_t)(centvalue*100.)/100.;// this is to extract a random number between 0 and 0.01
    
    if(fHistCentrDistr->GetBinContent(0)<-0.9999){
        if(gRandom->Uniform(1.)<bincont)return kTRUE;
        return kFALSE;
    }
    
    //if(TMath::Abs(centDigits*100.)<bincont)return kTRUE;
    if(centDigits*100.<bincont)return kTRUE;  // original from centrality flatt
    
    return kFALSE;
    
}
//---------------------------------------------------------------------------
void AliAnalysisTaskFlowITSTPCTOFQCSP::SetEtaRange(Double_t etaminimum, Double_t etamaximum)
{
    //Set PID cuts
    fmineta = etaminimum;
    fmaxeta = etamaximum;
}
//_____________________________________________________________________________
//_____________________________________________________________________________

