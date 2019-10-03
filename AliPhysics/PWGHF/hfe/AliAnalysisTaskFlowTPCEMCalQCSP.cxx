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
//  Task for Heavy Flavour Electron Flow with TPC plus EMCal          //
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
#include "AliAnalysisTaskFlowTPCEMCalQCSP.h"
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
#include "AliEMCALGeometry.h"
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
#include "AliPID.h"
#include "AliPIDResponse.h"
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

ClassImp(AliAnalysisTaskFlowTPCEMCalQCSP)
//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalQCSP::AliAnalysisTaskFlowTPCEMCalQCSP(const char *name)
: AliAnalysisTaskSE(name)
,fDebug(0)
,fAOD(0)
,fVevent(0)
,fGeom(0)
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
,fpTCut(0)
,fTrigger(0)
,fPhi(0)
,fEta(0)
,fVZEROA(0)
,fVZEROC(0)
,fTPCM(0)
,fNoEvents(0)
,fTrkEovPBef(0)
//,fdEdxBef(0)
,fInclusiveElecPt(0)
,fTPCnsigma(0)
,fTPCnsigmaAft(0)
,fCentralityPass(0)
,fCentralityNoPass(0)
,fInvmassLS1(0)
,fInvmassULS1(0)
,fPhotoElecPt(0)
,fSemiInclElecPt(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fminTPC(-1)
,fmaxTPC(3)
,fminEovP(0.8)
,fmaxEovP(1.2)
,fminM20(0.03)
,fmaxM20(0.3)
,fminM02(0.03)
,fmaxM02(0.5)
,fDispersion(1)
,fMultCorAfterCuts(0)
,fMultvsCentr(0)
,fSubEventDPhiv2(0)
,EPVzA(0)
,EPVzC(0)
,EPTPC(0)
,fV2Phi(0)
,fSparseElectronHadron(0)
,fvertex(0)
,fMultCorBeforeCuts(0)
,fSideBandsFlow(kFALSE)
,fPhiminusPsi(kFALSE)
,fFlowEventCont(0) //! flow events (one for each inv mass band)
,fpurity(kFALSE)
,fSparseElectronpurity(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fNonHFE(new AliSelectNonHFE)
,fDCA(0)
,fOpeningAngleCut(0)
,fOP_angle(0)
,fAssoTPCCluster(0)
,fAssoITSRefit(0)
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
,fInvmassLS1highpt(0)
,fInvmassULS1highpt(0)
,fSparsephipsiULS(0)
,fSparsephipsiLS(0)
,fSparseMassULS(0)
,fSparseMassLS(0)
,fHistEPDistrWeight(0)
,EPweights(0)
,EPVzAftW(0)
,multCorrection(0)
,fptminAsso(0)
,fEtaMinimumPositive(0)
,fEtaMinimumNegative(0)
,fCentralityAll(0)
,fCentFlatMine(0)
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
    if(fSideBandsFlow){
        DefineOutput(3, AliFlowEventSimple::Class());
    }
    //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalQCSP::AliAnalysisTaskFlowTPCEMCalQCSP()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskFlowTPCEMCalQCSP")
,fDebug(0)
,fAOD(0)
,fVevent(0)
,fGeom(0)
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
,fpTCut(0)
,fTrigger(0)
,fPhi(0)
,fEta(0)
,fVZEROA(0)
,fVZEROC(0)
,fTPCM(0)
,fNoEvents(0)
,fTrkEovPBef(0)
//,fdEdxBef(0)
,fInclusiveElecPt(0)
,fTPCnsigma(0)
,fTPCnsigmaAft(0)
,fCentralityPass(0)
,fCentralityNoPass(0)
,fInvmassLS1(0)
,fInvmassULS1(0)
,fPhotoElecPt(0)
,fSemiInclElecPt(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fminTPC(-1)
,fmaxTPC(3)
,fminEovP(0.8)
,fmaxEovP(1.2)
,fminM20(0.03)
,fmaxM20(0.3)
,fminM02(0.03)
,fmaxM02(0.5)
,fDispersion(1)
,fMultCorAfterCuts(0)
,fMultvsCentr(0)
,fSubEventDPhiv2(0)
,EPVzA(0)
,EPVzC(0)
,EPTPC(0)
,fV2Phi(0)
,fSparseElectronHadron(0)
,fvertex(0)
,fMultCorBeforeCuts(0)
,fSideBandsFlow(kFALSE)
,fPhiminusPsi(kFALSE)
,fFlowEventCont(0) //! flow events (one for each inv mass band)
,fpurity(kFALSE)
,fSparseElectronpurity(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fNonHFE(new AliSelectNonHFE)
,fDCA(0)
,fOpeningAngleCut(0)
,fOP_angle(0)
,fAssoTPCCluster(0)
,fAssoITSRefit(0)
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
,fInvmassLS1highpt(0)
,fInvmassULS1highpt(0)
,fSparsephipsiULS(0)
,fSparsephipsiLS(0)
,fSparseMassULS(0)
,fSparseMassLS(0)
,fHistEPDistrWeight(0)
,EPweights(0)
,EPVzAftW(0)
,multCorrection(0)
,fptminAsso(0)
,fEtaMinimumPositive(0)
,fEtaMinimumNegative(0)
,fCentralityAll(0)
,fCentFlatMine(0)
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
    //  DefineOutput(3, TTree::Class());
    if(fSideBandsFlow){
        DefineOutput(3, AliFlowEventSimple::Class());
    }
    //DefineOutput(3, TTree::Class());
}
//_________________________________________

AliAnalysisTaskFlowTPCEMCalQCSP::~AliAnalysisTaskFlowTPCEMCalQCSP()
{
    //Destructor
    
    delete fOutputList;
    delete fGeom;
    delete fPID;
    delete fCFM;
    delete fPIDqa;
    if (fDCA)  delete fNonHFE;
    if (fOutputList) delete fOutputList;
    if (fFlowEvent) delete fFlowEvent;
    if (fFlowEventCont) delete fFlowEventCont;
    
}
//_________________________________________

void AliAnalysisTaskFlowTPCEMCalQCSP::UserExec(Option_t*)
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
    
    //  cout << "kTrigger   ==   " << fTrigger <<endl;
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
    if(TMath::Abs(vtxz)>10)return;
    
    // Event cut
    if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fAOD)) return;
    if(fNOtrks<2) return;
    
    Bool_t pass = kFALSE; //to select centrality and pile up protection
    CheckCentrality(fAOD,pass);
    if(!pass)return;
    fvertex->Fill(vtxz);
    
    fNoEvents->Fill(0);
    PlotVZeroMultiplcities(fAOD);
    
    SetNullCuts(fAOD);
    PrepareFlowEvent(fAOD->GetNumberOfTracks(),fFlowEvent);    //Calculate event plane Qvector and EP resolution for inclusive
    
    if(fSideBandsFlow){
        PrepareFlowEvent(fAOD->GetNumberOfTracks(),fFlowEventCont);    //Calculate event plane Qvector and EP resolution for inclusive
    }
    
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
        //----------hfe begin---------
        if(track->Eta()<-0.7 || track->Eta()>0.7)    continue;    //eta cuts on candidates
        
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
        
        Double_t fClsE = -999, p = -999, fEovP=-999, pt = -999, fTPCnSigma=0;
        // Track extrapolation
        Int_t fClsId = track->GetEMCALcluster();
        if(fClsId < 0) continue;
        AliAODCaloCluster *cluster = fAOD->GetCaloCluster(fClsId);
        if(TMath::Abs(cluster->GetTrackDx()) > 0.05 || TMath::Abs(cluster->GetTrackDz()) > 0.05) continue;
        
        pt = track->Pt();         //pt track after cuts
        if(pt<fpTCut) continue;
        fClsE = cluster->E();
        p = track->P();
        // dEdx = track->GetTPCsignal();
        fEovP = fClsE/p;
        fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
        
        
        Double_t CorrectTPCNSigma;
        Double_t mult = fVevent->GetNumberOfESDTracks()/8;
        
        if(multCorrection){
            CorrectTPCNSigma = ftpcpid->GetCorrectedTPCnSigma(track->Eta(), mult, fTPCnSigma);
            // cout <<fTPCnSigma << "   ====  " <<COrrectTPCNSigma<<endl;
            fTPCnSigma = CorrectTPCNSigma;
            // cout <<fTPCnSigma << "   ====  " <<COrrectTPCNSigma<<endl;
        }
        
        
        Double_t m20 =cluster->GetM20();
        Double_t m02 =cluster->GetM02();
        Double_t disp=cluster->GetDispersion();
        if(fTPCnSigma >= -1 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,fEovP);
        fTPCnsigma->Fill(p,fTPCnSigma);
        //   fdEdxBef->Fill(p,dEdx);
        Double_t eta = track->Eta();
        Double_t phi = track->Phi();
        //-----------------------Phiminupsi method to remove the contamination-----------------------------------------------
        //-----------------------fTPCnSigma < -3.5 hadrons will be selected from this region--------------------------
        Float_t dPhi_aeh = TVector2::Phi_0_2pi(phi - evPlAngV0A);
        if(dPhi_aeh > TMath::Pi()) dPhi_aeh = dPhi_aeh - TMath::Pi();
        Float_t dPhi_ceh = TVector2::Phi_0_2pi(phi - evPlAngV0C);
        if(dPhi_ceh > TMath::Pi()) dPhi_ceh = dPhi_ceh - TMath::Pi();
        
        if(fPhiminusPsi){
            Double_t valueElh[8] = {
                pt,
                fEovP,
                fTPCnSigma,
                m20,
                m02,
                disp,
                dPhi_aeh,
                dPhi_ceh};
            fSparseElectronHadron->Fill(valueElh);
        }
        //----------------------------------------------------------------------------------------------------------
        //---------------------------From here usual electron selection---------------------------------------------
        //----------------------------------------------------------------------------------------------------------
        if(m20 < fminM20 || m20 > fmaxM20) continue;
        if(m02 < fminM02 || m02 > fmaxM02) continue;
        if(disp > fDispersion ) continue;
        //---------------------------------for purity---------------------------------------------------------------
        if(fpurity){
            Double_t valuepurity[3] = {
                pt,
                fEovP,
                fTPCnSigma};
            fSparseElectronpurity->Fill(valuepurity);
        }
        //----------------------------------------------------------------------------------------------------------
        //----------------------------------------------------------------------------------------------------------
        if(fTPCnSigma < fminTPC || fTPCnSigma > fmaxTPC) continue; //cuts on nsigma tpc and EoP
        //===============================Flow Event for Contamination=============================================
        if(fSideBandsFlow){
            if(fEovP>0 && fEovP<0.6){
                AliFlowTrack *sTrackCont = new AliFlowTrack();
                sTrackCont->Set(track);
                sTrackCont->SetID(track->GetID());
                sTrackCont->SetForRPSelection(kTRUE);
                sTrackCont->SetForPOISelection(kTRUE);
                sTrackCont->SetMass(2637);
                for(int iRPs=0; iRPs!=fFlowEventCont->NumberOfTracks(); ++iRPs)
                {
                    //   cout << " no of rps " << iRPs << endl;
                    AliFlowTrack *iRPCont = dynamic_cast<AliFlowTrack*>(fFlowEventCont->GetTrack( iRPs ));
                    if (!iRPCont) continue;
                    if (!iRPCont->InRPSelection()) continue;
                    if( sTrackCont->GetID() == iRPCont->GetID())
                    {
                        if(fDebug) printf(" was in RP set");
                        //       cout << sTrack->GetID() <<"   ==  " << iRP->GetID() << " was in RP set" <<endl;
                        iRPCont->SetForRPSelection(kFALSE);
                        //  fFlowEventCont->SetNumberOfRPs(fFlowEventCont->GetNumberOfRPs() - 1);
                    }
                } //end of for loop on RPs
                fFlowEventCont->InsertTrack(((AliFlowTrack*) sTrackCont));
                fFlowEventCont->SetNumberOfPOIs(fFlowEventCont->GetNumberOfPOIs()+1);
                
            }
        }
        //==========================================================================================================
        //===============================From here eovP cut is used fro QC, SP and EPV0=============================
        if(fEovP < fminEovP || fEovP >fmaxEovP) continue;
        //==========================================================================================================
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
        
        
        
        //=========================================================================================================
        fTPCnsigmaAft->Fill(p,fTPCnSigma);
        fInclusiveElecPt->Fill(pt);
        fPhi->Fill(phi);
        fEta->Fill(eta);
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
                //       cout << sTrack->GetID() <<"   ==  " << iRP->GetID() << " was in RP set" <<endl;
                iRP->SetForRPSelection(kFALSE);
                // fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
            }
        } //end of for loop on RPs
        fFlowEvent->InsertTrack(((AliFlowTrack*) sTrack));
        fFlowEvent->SetNumberOfPOIs(fFlowEvent->GetNumberOfPOIs()+1);
        
        
        if(fDCA){
            //----------------------Selection of Photonic Electrons DCA-----------------------------
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
            //----------------------Selection of Photonic Electrons KFParticle-----------------------------
            Bool_t fFlagPhotonicElec = kFALSE;
            SelectPhotonicElectron(iTracks,track,fEovP, evPlAngV0, fFlagPhotonicElec,weightEP,mult);
            if(fFlagPhotonicElec){fPhotoElecPt->Fill(pt);}
            // Semi inclusive electron
            if(!fFlagPhotonicElec){fSemiInclElecPt->Fill(pt);}
        }
        
        
        
    }//end loop on track
    
    PostData(1, fOutputList);
    PostData(2, fFlowEvent);
    if(fSideBandsFlow){
        PostData(3, fFlowEventCont);
    }
    
    //----------hfe end---------
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::SelectPhotonicElectron(Int_t itrack,const AliAODTrack *track,Double_t fEovP,Double_t evPlAngV0, Bool_t &fFlagPhotonicElec, Double_t weightEPflat, Double_t multev)
{
    //Identify non-heavy flavour electrons using Invariant mass method KF
    
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
        //    if((!(trackAsso->GetStatus()&AliESDtrack::kITSrefit) || (!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
        if(fAssoITSRefit){
            if(!(trackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
        }
        
        if(!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
        
        if(jTracks == itrack) continue;
        Double_t ptAsso=-999., nsigma=-999.0;
        Double_t mass=-999., width = -999;
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t openingAngle = -999.;
        Double_t ptcutonmasshighpt = track->Pt();
        
        ptAsso = trackAsso->Pt();
        Short_t chargeAsso = trackAsso->Charge();
        Short_t charge = track->Charge();
        
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(ptAsso <fptminAsso) continue;
        nsigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
        
        Double_t CorrectTPCNSigma;
        if(multCorrection){
            CorrectTPCNSigma = ftpcpid->GetCorrectedTPCnSigma(trackAsso->Eta(), multev, nsigma);
            nsigma = CorrectTPCNSigma;
        }
        
        //80
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
        if(fOP_angle)if(openingAngle > fOpeningAngleCut) continue;
        
        
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
        
        
        if(ptcutonmasshighpt >= 8.){
            if(fFlagLS) fInvmassLS1highpt->Fill(mass);
            if(fFlagULS) fInvmassULS1highpt->Fill(mass);
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
                    fEovP,
                    DeltaPhi_eEP
                };
                if(EPweights) fSparsephipsiULS->Fill(ulsSparse,weightEPflat);
                if(!EPweights) fSparsephipsiULS->Fill(ulsSparse);
            }
            if(fFlagLS){
                Double_t lsSparse[3] = {
                    track->Pt(),
                    fEovP,
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
//__________________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::UserCreateOutputObjects()
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
    cc->SetPtMax(50);
    
    cc->SetNbinsPhi(180);
    cc->SetPhiMin(0.0);
    cc->SetPhiMax(TMath::TwoPi());
    
    cc->SetNbinsEta(30);
    cc->SetEtaMin(-7.0);
    cc->SetEtaMax(+7.0);
    
    cc->SetNbinsQ(500);
    cc->SetQMin(0.0);
    cc->SetQMax(3.0);
    
    //--------Initialize PID
    fPID->SetHasMCData(kFALSE);
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("TPC", 0);
        fPID->AddDetector("EMCAL", 1);
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
    
    fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma before HFE pid",1000,0,50,200,-10,10);
    fOutputList->Add(fTPCnsigma);
    
    fTPCnsigmaAft = new TH2F("fTPCnsigmaAft", "TPC - n sigma after HFE pid",1000,0,50,200,-10,10);
    fOutputList->Add(fTPCnsigmaAft);
    
    fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",1000,0,50,100,0,2);
    fOutputList->Add(fTrkEovPBef);
    
    //    fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",1000,0,50,150,0,150);
    //     fOutputList->Add(fdEdxBef);
    
    fInclusiveElecPt = new TH1F("fInclElecPt", "Inclusive electron pt",1000,0,100);
    fOutputList->Add(fInclusiveElecPt);
    
    fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",1000,0,100);
    fOutputList->Add(fPhotoElecPt);
    
    fSemiInclElecPt = new TH1F("fSemiInclElecPt", "Semi-inclusive electron pt",1000,0,100);
    fOutputList->Add(fSemiInclElecPt);
    
    fULSElecPt = new TH1F("fULSElecPt", "ULS electron pt",1000,0,100);
    fOutputList->Add(fULSElecPt);
    
    fLSElecPt = new TH1F("fLSElecPt", "LS electron pt",1000,0,100);
    fOutputList->Add(fLSElecPt);
    
    fInvmassLS1 = new TH1F("fInvmassLS1", "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassLS1);
    
    fInvmassULS1 = new TH1F("fInvmassULS1", "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassULS1);
    
    fInvmassLS1highpt = new TH1F("fInvmassLS1highpt", "Inv mass of LS (e,e); mass(GeV/c^2) highpt; counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassLS1highpt);
    
    fInvmassULS1highpt = new TH1F("fInvmassULS1highpt", "Inv mass of ULS (e,e); mass(GeV/c^2) highpt; counts;", 1000,0,1.0);
    fOutputList->Add(fInvmassULS1highpt);
    
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
    
    
    //================================Event Plane with VZERO A & C=====================
    const Int_t nPtBins = 12;
    Double_t binsPt[nPtBins+1] = {0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10, 13};
    // v2A, v2C, pt
    Int_t    bins[3] = {  50,  50, nPtBins};
    Double_t xmin[3] = { -1., -1.,   0};
    Double_t xmax[3] = {  1.,  1.,   13};
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
    // const Int_t nPtBins = 10;
    // Double_t binsPt[nPtBins+1] = {0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    // v2, pt
    Int_t    binsV[2] = {  50,  nPtBins};
    Double_t xminV[2] = { -1.,   0};
    Double_t xmaxV[2] = {  1.,   13};
    fV2Phivzerotot = new THnSparseF("fV2Phivzerotot", "v2:pt", 2, binsV, xminV, xmaxV);
    // Set bin limits for axes which are not standard binned
    fV2Phivzerotot->SetBinEdges(1, binsPt);
    // set axes titles
    fV2Phivzerotot->GetAxis(0)->SetTitle("v_{2} (V0)");
    fV2Phivzerotot->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    fV2Phivzerotot->Sumw2();
    
    fOutputList->Add(fV2Phivzerotot);
    
    
    
    //----------------------------------------------------------------------------
    
    if(fPhiminusPsi){
        Int_t binsvElectH[8]={ 600,  200, 200 ,100,  100,  100,   10,          10}; //pt, E/p,TPCnSigma,M20,M02,Disp Phi-psiV0A ,Phi-PsiV0C,eta (commented)
        Double_t xminvElectH[8]={0,    0, -10 ,  0,    0,    0,    0,           0};
        Double_t xmaxvElectH[8]={20,   2,  10 ,  2,    2,    2,  TMath::Pi(), TMath::Pi()};
        fSparseElectronHadron = new THnSparseD("ElectronHadron","ElectronHadron",8,binsvElectH,xminvElectH,xmaxvElectH);
        fSparseElectronHadron->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fSparseElectronHadron->GetAxis(1)->SetTitle("EovP");
        fSparseElectronHadron->GetAxis(2)->SetTitle("TPCnSigma");
        fSparseElectronHadron->GetAxis(3)->SetTitle("M20");
        fSparseElectronHadron->GetAxis(4)->SetTitle("M02");
        fSparseElectronHadron->GetAxis(5)->SetTitle("Disp");
        fSparseElectronHadron->GetAxis(6)->SetTitle("phiminuspsi V0A");
        fSparseElectronHadron->GetAxis(7)->SetTitle("phiminuspsi V0C");
        fOutputList->Add(fSparseElectronHadron);
    }
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    if(fpurity){
        Int_t binsvpurity[3]={   600,200, 200}; //pt, E/p,TPCnSigma
        Double_t xminvpurity[3]={0,    0, -10};
        Double_t xmaxvpurity[3]={30,   2,  10};
        fSparseElectronpurity = new THnSparseD("Electronpurity","Electronpurity",3,binsvpurity,xminvpurity,xmaxvpurity);
        fSparseElectronpurity->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fSparseElectronpurity->GetAxis(1)->SetTitle("EovP");
        fSparseElectronpurity->GetAxis(2)->SetTitle("TPCnSigma");
        fOutputList->Add(fSparseElectronpurity);
    }
    //----------------------------------------------------------------------------
    
    
    Int_t    binsphipsi[3] = { 100,   200,           6};
    Double_t xminphipsi[3] = { 0.,      0,           0};
    Double_t xmaxphipsi[3] = { 10.,      2, TMath::Pi()};
    fSparsephipsiULS = new THnSparseF("fSparsephipsiULS", "pt:eop:DeltaPhiULS", 3, binsphipsi, xminphipsi, xmaxphipsi);
    fSparsephipsiULS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparsephipsiULS->GetAxis(1)->SetTitle("eop");
    fSparsephipsiULS->GetAxis(2)->SetTitle("DeltaPhiULS");
    fSparsephipsiULS->Sumw2();
    fOutputList->Add(fSparsephipsiULS);
    
    fSparsephipsiLS = new THnSparseF("fSparsephipsiLS", "pt:eop:DeltaPhiLS", 3, binsphipsi, xminphipsi, xmaxphipsi);
    fSparsephipsiLS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparsephipsiLS->GetAxis(1)->SetTitle("eop");
    fSparsephipsiLS->GetAxis(2)->SetTitle("DeltaPhiLS");
    fSparsephipsiLS->Sumw2();
    fOutputList->Add(fSparsephipsiLS);
    
    Int_t    binsmass[2] = { 100, 200};
    Double_t xminmass[2] = { 0.,  0};
    Double_t xmaxmass[2] = { 10., 1.};
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
    
    if(fSideBandsFlow){
        fFlowEventCont = new AliFlowEvent(10000);
        PostData(3, fFlowEventCont);
    }
}

//________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::Terminate(Option_t *)
{
    // Info("Terminate");
    AliAnalysisTaskSE::Terminate();
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskFlowTPCEMCalQCSP::PlotVZeroMultiplcities(const T* event) const
{
    // QA multiplicity plots
    fVZEROA->Fill(event->GetVZEROData()->GetMTotV0A());
    fVZEROC->Fill(event->GetVZEROData()->GetMTotV0C());
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskFlowTPCEMCalQCSP::SetNullCuts(T* event)
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
void AliAnalysisTaskFlowTPCEMCalQCSP::PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const
{
    //Prepare flow events
    FlowEv->ClearFast();
    FlowEv->Fill(fCutsRP, fNullCuts);
    FlowEv->SetReferenceMultiplicity(iMulti);
    FlowEv->DefineDeadZone(0, 0, 0, 0);
    //  FlowEv->TagSubeventsInEta(-0.7, 0, 0, 0.7);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalQCSP::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    // Check single track cuts for a given cut step
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::CheckCentrality(AliAODEvent* event, Bool_t &centralitypass)
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
void AliAnalysisTaskFlowTPCEMCalQCSP::SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod)
{
    // Set a centrality range ]min, max] and define the method to use for centrality selection
    fCentralityMin = CentralityMin;
    fCentralityMax = CentralityMax;
    fkCentralityMethod = CentralityMethod;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::SetIDCuts(Double_t minTPC, Double_t maxTPC, Double_t minEovP, Double_t maxEovP, Double_t minM20, Double_t maxM20, Double_t minM02, Double_t maxM02, Double_t Dispersion)
{
    //Set ID cuts
    fminTPC = minTPC;
    fmaxTPC = maxTPC;
    fminEovP = minEovP;
    fmaxEovP = maxEovP;
    fminM20 = minM20;
    fmaxM20 = maxM20;
    fminM02 = minM02;
    fmaxM02 = maxM02;
    fDispersion = Dispersion;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::SetHistoForCentralityFlattening(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef,Int_t switchTRand){
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
Bool_t AliAnalysisTaskFlowTPCEMCalQCSP::IsEventSelectedForCentrFlattening(Float_t centvalue){
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
void AliAnalysisTaskFlowTPCEMCalQCSP::SetHistoForEPFlattWeights(TH1D *h){
    
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
Double_t AliAnalysisTaskFlowTPCEMCalQCSP::GiveMeWeight(Double_t EP){
    
    Int_t Bin = fHistEPDistrWeight->FindBin(EP);
    Double_t ww = fHistEPDistrWeight->GetBinContent(Bin);
    return ww;
    
}
//-------------------------------------------------

//_____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalQCSP::SetHistoForCentralityFlattening_Bis(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef){
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
Bool_t AliAnalysisTaskFlowTPCEMCalQCSP::IsEventSelectedForCentrFlattening_Bis(Double_t centvalue){
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





//_____________________________________________________________________________




