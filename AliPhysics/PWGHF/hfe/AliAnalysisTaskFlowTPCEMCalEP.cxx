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

// Class for heavy-flavour electron v2 with EMCal triggered events
// Author: Denise Godoy


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
#include "AliVEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisTaskFlowTPCEMCalEP.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliEMCALTrack.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliQnCorrectionsManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"


ClassImp(AliAnalysisTaskFlowTPCEMCalEP)

using std::cout;
using std::endl;

//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalEP::AliAnalysisTaskFlowTPCEMCalEP(const char *name)
: AliAnalysisTaskSE(name)
,fWhichPeriod(2015)
,fAssPtCut(0.5)
,fITSncut(3)
,fAssTPCnCut(80)
,fTPCnCut(100)
,fAssITSrefitCut(kTRUE)
,fUseNewEP(kTRUE)
,fUseTender(kTRUE)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(0.)
,fTimeCut(kFALSE)
,fWeightSyst(kFALSE)
,fSystTOFcut(kFALSE)
,fCutM02(2.)
,fCutM20(2.)
,fSScut(kFALSE)
,fEnablePileupRejVZEROTPCout(kFALSE)
//,fESD(0)
,fAOD(0)
,fVevent(0)
,fpidResponse(0)
,fMultSelection(0)
,fCentrality(0)
,fMC(0)
,fStack(0)
,fMCparticle(0)
,fMCarray(0)
,fMCheader(0)
,fTracks_tender(0)
,fCaloClusters_tender(0)
,fOutputList(0)
//,fTrackCuts(0)
//,fAssTrackCuts(0)
,fCuts(0)
,fIdentifiedAsOutInz(kFALSE)
,fPassTheEventCut(kFALSE)
,fRejectKinkMother(kFALSE)
,fIsMC(kFALSE)
,fIsAOD(kFALSE)
,fSetMassConstraint(kFALSE)
,fVz(0.0)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,flowQnVectorTask(0)
,fFlowQnVectorMgr(0)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fWhichDecay(0)
,fPi0EtaWeight(1.)
,fCentAftThr(0)
,fTrigger(0)
,fNoEvents(0)
,fTrkpt(0)
,fTrackPtBefTrkCuts(0)
,fTrackPtAftTrkCuts(0)
,fChargedParticlePhi(0)
,fElectronPhi(0)
,fCent(0)
,fCentAftFlt(0)
,fTPCsubEPres(0)
,fCorr(0)
,fElecMC(0)
,fcpV2_3040(0)
,fcpV2_4050(0)
,fpassV0("align")
,fpassTPC("twist")

{
    //Named constructor
    
    for(Int_t k = 0; k < 3; k++) {
        fevPlaneV0[k] = NULL;
        fevPlaneV0AftThr[k] = NULL;
        feTPCV2[k] = NULL;
        feV2[k] = NULL;
        fChargPartV2[k] = NULL;
        fMtcPartV2[k] = NULL;
        fEPres[k] = NULL;
        
        fPi0Pt[k] = NULL;
        fEtaPt[k] = NULL;
        fElecPtULSInvmassCut[k] = NULL;
        fElecPtLSInvmassCut[k] = NULL;
        fElecPtInvmassCut[k] = NULL;
        fInclElec[k] = NULL;
        fInvmassLS[k] = NULL;
        fInvmassULS[k] = NULL;
        fOpeningAngleLS[k] = NULL;
        fOpeningAngleULS[k] = NULL;
        
        fHistITSnSig[k] = NULL;
        fHistTOFnSig[k] = NULL;
        fHistTPCnSig[k] = NULL;
        fHistTPCnSigITScut[k] = NULL;
        fHistTPCnSigTOFcut[k] = NULL;
        fHistTPCnSigITSTOFcut[k] = NULL;
        fHistTPCnSigEop[k] = NULL;
        fHistTPCnSigEMCalnSig[k] = NULL;
        
        fHistM02sig[k] = NULL;
        fHistM20sig[k] = NULL;
        fHistM02backg[k] = NULL;
        fHistM20backg[k] = NULL;
        fHistM02EoverP[k] = NULL;
        fHistM20EoverP[k] = NULL;
        
        fEoverPsignalTPC[k] = NULL;
        fEoverPsignalTPCM02[k] = NULL;
        fEoverPbackg[k] = NULL;
    }
    
    
    fPID = new AliHFEpid("hfePid");
    //fTrackCuts = new AliAODTrackCuts();
    //fAssTrackCuts = new AliAODTrackCuts();
    
    InitParameters();
    
    DefineInput(0, TChain::Class());
    //DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalEP::AliAnalysisTaskFlowTPCEMCalEP()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
,fWhichPeriod(2015)
,fAssPtCut(0.5)
,fITSncut(3)
,fAssTPCnCut(80)
,fTPCnCut(100)
,fAssITSrefitCut(kTRUE)
,fUseNewEP(kTRUE)
,fUseTender(kTRUE)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(0.)
,fTimeCut(kFALSE)
,fWeightSyst(kFALSE)
,fSystTOFcut(kFALSE)
,fCutM02(2.)
,fCutM20(2.)
,fSScut(kFALSE)
,fEnablePileupRejVZEROTPCout(kFALSE)
//,fESD(0)
,fAOD(0)
,fVevent(0)
,fpidResponse(0)
,fMultSelection(0)
,fCentrality(0)
,fMC(0)
,fStack(0)
,fMCparticle(0)
,fMCarray(0)
,fMCheader(0)
,fTracks_tender(0)
,fCaloClusters_tender(0)
,fOutputList(0)
//,fTrackCuts(0)
//,fAssTrackCuts(0)
,fCuts(0)
,fIdentifiedAsOutInz(kFALSE)
,fPassTheEventCut(kFALSE)
,fRejectKinkMother(kFALSE)
,fIsMC(kFALSE)
,fIsAOD(kFALSE)
,fSetMassConstraint(kFALSE)
,fVz(0.0)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,flowQnVectorTask(0)
,fFlowQnVectorMgr(0)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fWhichDecay(0)
,fPi0EtaWeight(1.)
,fCentAftThr(0)
,fTrigger(0)
,fNoEvents(0)
,fTrkpt(0)
,fTrackPtBefTrkCuts(0)
,fTrackPtAftTrkCuts(0)
,fChargedParticlePhi(0)
,fElectronPhi(0)
,fCent(0)
,fCentAftFlt(0)
,fTPCsubEPres(0)
,fCorr(0)
,fElecMC(0)
,fcpV2_3040(0)
,fcpV2_4050(0)
,fpassV0("align")
,fpassTPC("twist")

{
    
    //Default constructor
    
    for(Int_t k = 0; k < 3; k++) {
        fevPlaneV0[k] = NULL;
        fevPlaneV0AftThr[k] = NULL;
        feTPCV2[k] = NULL;
        feV2[k] = NULL;
        fChargPartV2[k] = NULL;
        fMtcPartV2[k] = NULL;
        fEPres[k] = NULL;
        
        fPi0Pt[k] = NULL;
        fEtaPt[k] = NULL;
        fElecPtULSInvmassCut[k] = NULL;
        fElecPtLSInvmassCut[k] = NULL;
        fElecPtInvmassCut[k] = NULL;
        fInclElec[k] = NULL;
        fInvmassLS[k] = NULL;
        fInvmassULS[k] = NULL;
        fOpeningAngleLS[k] = NULL;
        fOpeningAngleULS[k] = NULL;
        
        fHistITSnSig[k] = NULL;
        fHistTOFnSig[k] = NULL;
        fHistTPCnSig[k] = NULL;
        fHistTPCnSigITScut[k] = NULL;
        fHistTPCnSigTOFcut[k] = NULL;
        fHistTPCnSigITSTOFcut[k] = NULL;
        fHistTPCnSigEop[k] = NULL;
        fHistTPCnSigEMCalnSig[k] = NULL;
        
        fHistM02sig[k] = NULL;
        fHistM20sig[k] = NULL;
        fHistM02backg[k] = NULL;
        fHistM20backg[k] = NULL;
        fHistM02EoverP[k] = NULL;
        fHistM20EoverP[k] = NULL;
        
        fEoverPsignalTPC[k] = NULL;
        fEoverPsignalTPCM02[k] = NULL;
        fEoverPbackg[k] = NULL;
    }
    
    
    
    fPID = new AliHFEpid("hfePid");
    //fTrackCuts = new AliAODTrackCuts();
    //fAssTrackCuts = new AliAODTrackCuts();
    
    InitParameters();
    
    DefineInput(0, TChain::Class());
    //DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
}
//_________________________________________

AliAnalysisTaskFlowTPCEMCalEP::~AliAnalysisTaskFlowTPCEMCalEP()
{
    //Destructor
    
    delete fOutputList;
    delete fPID;
    delete fCFM;
    delete fPIDqa;
    delete fTracks_tender;
    delete fCaloClusters_tender;
    //delete fTrackCuts;
    //delete fAssTrackCuts;
}
//_________________________________________

void AliAnalysisTaskFlowTPCEMCalEP::UserExec(Option_t*)
{
    //Main loop
    //Called for each event
    
    // create pointer to event
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if(!fVevent){
        printf("ERROR: fVEvent not available\n");
        return;
    }
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD){
        printf("ERROR: fAOD not available\n");
        return;
    }
    
    
    if(!fCuts){
        AliError("HFE cuts not available");
        return;
    }
    
    if(!fPID->IsInitialized()){
        // Initialize PID with the given run number
        AliWarning("PID not initialised, get from Run no");
        fPID->InitializePID(fAOD->GetRunNumber());
    }
    
    if(fUseTender){
        fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
        fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
        
        if (!fTracks_tender || !fCaloClusters_tender) return;
    }
    
    
    fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    Int_t NpureMC = -1;
    
    if (fMCheader){
        TList *lh=fMCheader->GetCocktailHeaders();
        if(lh){
            AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(0); //  0 for HIJING
            NpureMC = gh->NProduced();
        }
    }
    
    
    Int_t ntracks = -999;
    if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender)  ntracks = fTracks_tender->GetEntries();
    
    
    if(fIsMC)fMC = MCEvent();
    if(fIsMC && fMC) fStack = fMC->Stack();
    
    Int_t fNOtrks = fAOD->GetNumberOfTracks();
    
    
    // suggested by DPG to remove outliers
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    const AliAODVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    
    Double_t pVtxZ = -999;
    pVtxZ = pVtx->GetZ();
    
    if (pVtx->GetNContributors()<2 || spdVtx->GetNContributors()<1) return; // one of vertices is missing
    
    double covTrc[6],covSPD[6];
    pVtx->GetCovarianceMatrix(covTrc);
    spdVtx->GetCovarianceMatrix(covSPD);
    double dz = pVtx->GetZ()-spdVtx->GetZ();
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
    if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing
    
    if(TMath::Abs(pVtxZ)>10) return;
    
    // suggested by Ionut to remove pile-up
    
    Int_t nTPCout=0;
    Float_t mTotV0=0;
    
    AliAODVZERO* v0data=(AliAODVZERO*) fAOD->GetVZEROData();
    Float_t mTotV0A=v0data->GetMTotV0A();
    Float_t mTotV0C=v0data->GetMTotV0C();
    
    mTotV0=mTotV0A+mTotV0C;
    
    for(Int_t itrack=0; itrack<fNOtrks; itrack++) { // loop on tracks
        AliAODTrack * track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
        if(!track) {AliFatal("Not a standard AOD");}
        if(track->GetID()<0)continue;
        if((track->GetFlags())&(AliESDtrack::kTPCout)) nTPCout++;
        else continue;
    }
    Float_t mV0Cut=-2200.+(2.5*nTPCout)+(0.000012*nTPCout*nTPCout); //function to apply to pile-up rejection
    
    if(fEnablePileupRejVZEROTPCout && (mTotV0<mV0Cut)) return;

    
    fNoEvents->Fill(0);
    
    if(fNOtrks<2) return;
    
    fpidResponse = fInputHandler->GetPIDResponse();
    if(!fpidResponse){
        AliDebug(1, "Using default PID Response");
        fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
    }
    
    fPID->SetPIDResponse(fpidResponse);
    
    fCFM->SetRecEventInfo(fAOD);
    
    Float_t cent = -1.;
    
    fCentrality = fAOD->GetCentrality();
    
    fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    
    if(fMultSelection)
        cent = fMultSelection->GetMultiplicityPercentile("V0M", kFALSE);
    
    else
        cent = fCentrality->GetCentralityPercentile("V0M");
    
    
    
    Int_t iPt=8, iCent=3, iDeltaphi=4;
    
    fCent->Fill(cent);
    
    
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[nVertices];
    Int_t nMotherKink = 0;
    for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
        AliAODVertex *vertex = fAOD->GetVertex(ivertex);
        if(!vertex) continue;
        if(vertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
            if(!mother) continue;
            Int_t idmother = mother->GetID();
            listofmotherkink[nMotherKink] = idmother;
            nMotherKink++;
        }
    }
    
    
    if (cent>=0  && cent<10) iCent=0;
    if (cent>=10 && cent<30) iCent=1;
    if (cent>=30 && cent<50) iCent=2;
    if (!fUseNewEP && (cent<0 || cent>=50)) return;
    if (fUseNewEP && (cent<30 || cent>=50)) return; // only calibration for 30-50% at the moment
    
    
    // Centrality flattening
    
    Bool_t rejectToFlattenCent = kFALSE;
    Bool_t weightToFlattenCent = kTRUE;
    
    
    Double_t wEvent = 1.;
    Double_t EPweight = 1.; // only for old EP framework
    Double_t centWeight = 1.;
    
    
    Bool_t rejectEvent = kFALSE;
    Int_t centBin = fCent->FindBin(cent);
    
    if (weightToFlattenCent) centWeight =GetCentWeight(centBin);
    if (rejectToFlattenCent){// not implemented for 2015 yet
        rejectEvent = RejectEvent(cent,centBin);
        //if (iCent==0 && rejectEvent) return;
    }
    
    
    // Event plane
    Bool_t rejectToFlattenEP = kFALSE, weightToFlattenEP = kFALSE, rejectEventPlane = kFALSE;
    
    Double_t evPlaneV0A = -999., evPlaneV0C = -999., evPlaneV0 = -999., evPlaneTPC = -999., evPlaneTPCneg = -999., evPlaneTPCpos = -999.;
    
    TList *qnlist = 0x0;
    
    if (fUseNewEP){
        
        wEvent = EPweight*centWeight; // only centrality is flattened
        
        flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *> (AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
        
        if (flowQnVectorTask != NULL) {
            fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
            //fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);
        }
        else {
            AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
            return;
        }
        
        qnlist = fFlowQnVectorMgr->GetQnVectorList();
        if (!qnlist) {
            AliWarning("qnlist not found!!");
            return;
        }
        
        const AliQnCorrectionsQnVector *qnTPCpos;
        const AliQnCorrectionsQnVector *qnTPCneg;
        const AliQnCorrectionsQnVector *qnV0;
        
        qnV0      = GetQnVectorFromList(qnlist, "VZEROQoverM", fpassV0, fpassV0);
        qnTPCpos  = GetQnVectorFromList(qnlist, "TPCPosEtaQoverM", fpassTPC, fpassTPC);
        qnTPCneg  = GetQnVectorFromList(qnlist, "TPCNegEtaQoverM", fpassTPC, fpassTPC);
        
        // new way to get the EP, test it later
        //         qnV0 = fFlowQnVectorMgr->GetDetectorQnVector("VZERO", "latest", "raw");
        //         qnTPCpos = fFlowQnVectorMgr->GetDetectorQnVector("TPCPosEta", "latest", "plain");
        //         qnTPCneg = fFlowQnVectorMgr->GetDetectorQnVector("TPCNegEta", "latest", "plain");
        
        
        if (qnTPCpos != NULL)  evPlaneTPCpos = qnTPCpos->EventPlane(2);
        if (qnTPCneg != NULL)  evPlaneTPCneg = qnTPCneg->EventPlane(2);
        if (qnV0 != NULL)      evPlaneV0 = qnV0->EventPlane(2);
        
        if(evPlaneTPCpos <0) evPlaneTPCpos += TMath::Pi();
        if(evPlaneTPCneg <0) evPlaneTPCneg += TMath::Pi();
        if(evPlaneV0 <0)     evPlaneV0 += TMath::Pi();
        
    }
    
    
    
    
    if (!fUseNewEP){
        
        Double_t qxV0A = 0, qyV0A = 0, qxV0C = 0, qyV0C = 0, qxV0 = 0, qyV0 = 0,  qxTPC = 0, qyTPC = 0;
        TVector2 *qTPC = 0x0;
        TVector2 qVectorfortrack;
        
        evPlaneV0 = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,10,2,qxV0,qyV0));
        evPlaneV0A = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,8,2,qxV0A,qyV0A));
        evPlaneV0C = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,9,2,qxV0C,qyV0C));
        
        if(evPlaneV0 > TMath::Pi())  evPlaneV0 = evPlaneV0 - TMath::Pi();
        if(evPlaneV0A > TMath::Pi()) evPlaneV0A = evPlaneV0A - TMath::Pi();
        if(evPlaneV0C > TMath::Pi()) evPlaneV0C = evPlaneV0C - TMath::Pi();
        
        AliEventplane* TPCep =  fAOD->GetEventplane();
        if (TPCep && TPCep->GetQVector()) {
            qTPC = TPCep->GetQVector();
            qxTPC = qTPC->X();
            qyTPC = qTPC->Y();
            
            qVectorfortrack.Set(qxTPC,qyTPC);
            evPlaneTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;
        }
        
        
        // Event plane flattening
        
        weightToFlattenEP = kTRUE;
        Int_t epBin = fevPlaneV0[0]->FindBin(evPlaneV0);
        if (weightToFlattenEP) EPweight = GetEPweight(epBin,iCent);
        wEvent = EPweight*centWeight;
        
        
        if (rejectToFlattenCent){
            rejectEventPlane = RejectEventPlane(evPlaneV0,epBin);
            //if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA && rejectEventPlane) return;
        }
        
        
        // Event plane resolutions
        
        // >> 2 subevent method (only for TPC EP)
        
        TVector2 *qsub1a = TPCep->GetQsub1();
        TVector2 *qsub2a = TPCep->GetQsub2();
        Double_t evPlaneResTPC = -999.;
        if(qsub1a && qsub2a){
            evPlaneResTPC = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
        }
        
        fTPCsubEPres->Fill(evPlaneResTPC);
        
        // >> 3 event method (V0, V0A, and V0C EP)
        
        Double_t Qx2pos = 0., Qy2pos = 0., Qx2neg = 0., Qy2neg = 0., Qweight = 1;
        
        for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
            
            
            AliVParticle* Vtrack = 0x0;
            if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
            if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
            
            if (!Vtrack) {
                printf("ERROR: Could not receive track for EP %d\n", iTracks);
                continue;
            }
            AliAODTrack *trackEP = dynamic_cast<AliAODTrack*>(Vtrack);
            
            if (!trackEP) {
                printf("ERROR: Could not receive track %d\n", iTracks);
                continue;
            }
            
            if(!trackEP->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
            
            
            if(TMath::Abs(trackEP->Eta())>0.7) continue;
            
            //if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, trackEP)) continue;
            
            if (trackEP->Pt() < 2) Qweight = trackEP->Pt()/2;
            if (trackEP->Pt() >= 2) Qweight = 1;
            
            
            if(trackEP->Eta()>0 && trackEP->Eta()<0.8){
                Qx2pos += Qweight*TMath::Cos(2*trackEP->Phi());
                Qy2pos += Qweight*TMath::Sin(2*trackEP->Phi());
            }
            if(trackEP->Eta()<0 && trackEP->Eta()>-0.8){
                Qx2neg += Qweight*TMath::Cos(2*trackEP->Phi());
                Qy2neg += Qweight*TMath::Sin(2*trackEP->Phi());
            }
        }//track loop only for EP
        
        evPlaneTPCneg = TMath::ATan2(Qy2neg, Qx2neg)/2;
        evPlaneTPCpos = TMath::ATan2(Qy2pos, Qx2pos)/2;
        
    }// old EP
    
    
    // centrality after flattening
    if (weightToFlattenCent) fCentAftFlt->Fill(cent,wEvent);
    
    if (weightToFlattenEP){
        fEPres[0]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCpos),wEvent);
        fEPres[1]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCneg),wEvent);
        fEPres[2]->Fill(cent,GetCos2DeltaPhi(evPlaneTPCpos,evPlaneTPCneg),wEvent);
        fevPlaneV0[iCent]->Fill(evPlaneV0,wEvent);
    }
    else{
        fEPres[0]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCpos));
        fEPres[1]->Fill(cent,GetCos2DeltaPhi(evPlaneV0,evPlaneTPCneg));
        fEPres[2]->Fill(cent,GetCos2DeltaPhi(evPlaneTPCpos,evPlaneTPCneg));
        fevPlaneV0[iCent]->Fill(evPlaneV0);
    }
    
    // Selection of pi0 and eta in MC to compute the weight
    
    Bool_t isPrimary = kFALSE, isFromLMdecay = kTRUE, isFromHFdecay=kTRUE;
    
    if(fMCarray){
        Int_t nParticles = fMCarray->GetEntries();
        for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
            AliAODMCParticle* particle = (AliAODMCParticle*) fMCarray->At(iParticle);
            int fPDG = particle->GetPdgCode();
            double pTMC = particle->Pt();
            
            Bool_t iEnhance = kFALSE;
            if(iParticle>=NpureMC)iEnhance = kTRUE;
            
            Double_t etaMC = particle->Eta();
            if (TMath::Abs(etaMC)>1.2)continue;
            
            isPrimary = IsPrimary(particle);
            isFromLMdecay = IsFromLMdecay(particle);
            isFromHFdecay = IsFromHFdecay(particle);
            
            if (isPrimary){
                if(fPDG==111) fPi0Pt[iCent]->Fill(iEnhance,pTMC); //pi0
                if(fPDG==221) fEtaPt[iCent]->Fill(iEnhance,pTMC); //eta
            }
            if (!isFromHFdecay && !isFromLMdecay){
                if(fPDG==111) fPi0Pt[iCent]->Fill(iEnhance+2,pTMC); //pi0
                if(fPDG==221) fEtaPt[iCent]->Fill(iEnhance+2,pTMC); //eta
            }
        }
    }//MC
    
    
    Double_t ptRange[7] = {1.5,2,2.5,3,4,6,8};
    Double_t deltaPhiRange[4];
    for(Int_t j=0;j<4;j++){
        deltaPhiRange[j] = j*(TMath::Pi()/4);
    }
    
    Bool_t IsSameEvent = kFALSE;
    
    // Track loop
    for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        
        AliVParticle* Vtrack = 0x0;
        if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
        if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
        
        if (!Vtrack) {
            printf("ERROR: Could not receive tagged  track %d\n", iTracks);
            continue;
        }
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(Vtrack);
        if(!track) continue;
        
        if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
        
        fTrackPtBefTrkCuts->Fill(track->Pt());
        
        // HFE cuts
        
        //if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
        if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
            if(track->GetKinkIndex(0) != 0) continue;
        }
        
        //if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
        
        //if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        //if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        
        // track cuts (same as HFE)
        
        //        Bool_t isGoodIEtrack = InclElecTrackCuts(track);
        //        if(!isGoodIEtrack) continue;
        
        if(TMath::Abs(track->Eta())>0.7) continue;
        if(track->GetTPCNcls() < fTPCnCut) continue;
        
        if (track->Pt()<3 && track->GetITSNcls() < fITSncut) continue; // for ITS+TOF+TPC analysis, ITS<5
        if (track->Pt()>=3 && track->GetITSNcls() < 3) continue; // for TPC+EMcal analysis, ITS<3
        
        
        if(!track->IsOn(AliAODTrack::kITSrefit)) continue;
        if(!track->IsOn(AliAODTrack::kTPCrefit)) continue;
        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
        if(track->GetTPCFoundFraction() < 0.6) continue;
        
        
        if (track->Pt()<0.5) continue;
        
        Bool_t kinkmotherpass = kTRUE;
        for(Int_t kinkmother = 0; kinkmother < nMotherKink; kinkmother++) {
            if(track->GetID() == listofmotherkink[kinkmother]) {
                kinkmotherpass = kFALSE;
                continue;
            }
        }
        if(!kinkmotherpass) continue;
        
        Double_t d0z0[2]={-999,-999}, cov[3];
        
        if(track->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
            if(TMath::Abs(d0z0[0]) > 2.4 || TMath::Abs(d0z0[1]) > 3.2) continue;
        
        // analysis
        

        
        fTrackPtAftTrkCuts->Fill(track->Pt());
        
        Double_t clsE=-9.,p=-99.,EovP=-99.,pt=-99.,dEdx=-99.,fTPCnSigma=9.,phi=-9.,m02=-9.,m20=-9.,fEMCalnSigma=9.,dphi=9.,cosdphi=9.,fITSnSigma=9,fTOFnSigma=9, clsTime = -9;
        
        pt = track->Pt();
        fTrkpt->Fill(pt);
        for(Int_t i=0;i<6;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                iPt=i;
                continue;
            }
        }
        
        p = track->P();
        phi = track->Phi();
        dEdx = track->GetTPCsignal();
        
        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        fITSnSigma = fpidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
        fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        
        Double_t emcphimim = 1.39;
        Double_t emcphimax = 3.265;
        
        Int_t clsId = track->GetEMCALcluster();
        if (clsId>0){
            AliVCluster *cluster=0x0;
            if(!fUseTender) cluster = (AliVCluster*)fVevent->GetCaloCluster(clsId);
            if(fUseTender) cluster = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(clsId));
            
            if(cluster && cluster->IsEMCAL() && phi>emcphimim && phi<emcphimax){
                clsE = cluster->E();
                m20 = cluster->GetM20();
                m02 = cluster->GetM02();
                clsTime = cluster->GetTOF()*1e+9; // ns
            }
        }
        
        
        if(fUseTender && !fMCarray && fTimeCut && pt>=3 && (TMath::Abs(clsTime)>50)) continue;
        
        EovP = clsE/p;
        
        dphi = GetDeltaPhi(phi,evPlaneV0);
        cosdphi = GetCos2DeltaPhi(phi,evPlaneV0);
        
        
        Double_t fTPCnSigma_Pion=9., fTOFnSigma_Pion=9.;
        fTPCnSigma_Pion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        fTOFnSigma_Pion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        
        
        if (cent>=30 && cent <40 && TMath::Abs(fTPCnSigma_Pion)<2 && TMath::Abs(fTOFnSigma_Pion)<2) fcpV2_3040->Fill(pt,cosdphi,wEvent);
        if (cent>=40 && cent <50 && TMath::Abs(fTPCnSigma_Pion)<2 && TMath::Abs(fTOFnSigma_Pion)<2) fcpV2_4050->Fill(pt,cosdphi,wEvent);
        
        
        // not implemented yet for 2015
        //if(fIsMC) fEMCalnSigma = GetSigmaEMCalMC(EovP, pt, iCent);
        //else fEMCalnSigma = GetSigmaEMCal(EovP, pt, iCent);
        
        fChargedParticlePhi->Fill(phi); // phi of charged particles
        
        
        fHistITSnSig[iCent]->Fill(p,fITSnSigma);
        fHistTOFnSig[iCent]->Fill(p,fTOFnSigma);
        fHistTPCnSig[iCent]->Fill(p,fTPCnSigma);
        if (fITSnSigma>-2 && fITSnSigma<2) fHistTPCnSigITScut[iCent]->Fill(p,fTPCnSigma);
        if (fTOFnSigma>-2 && fTOFnSigma<2)fHistTPCnSigTOFcut[iCent]->Fill(p,fTPCnSigma);
        if (fITSnSigma>-2 && fITSnSigma<2 && fTOFnSigma>-2 && fTOFnSigma<2) fHistTPCnSigITSTOFcut[iCent]->Fill(p,fTPCnSigma);
        if (pt>2) fHistTPCnSigEop[iCent]->Fill(fTPCnSigma,EovP);
        if (pt>2) fHistTPCnSigEMCalnSig[iCent]->Fill(fTPCnSigma,fEMCalnSigma);

        
        if (fTPCnSigma>0  && fTPCnSigma<3 && fTOFnSigma>-2 && fTOFnSigma<2 && fITSnSigma>-2 && fITSnSigma<2)
            fElectronPhi->Fill(phi); // phi of electron candidates
        
        if(!fMCarray){
            if (fTPCnSigma>-1  && fTPCnSigma<3 && EovP>0.8 && EovP<1.2){
                fHistM02sig[iCent]->Fill(pt,m02);
                fHistM20sig[iCent]->Fill(pt,m20);
            }
            
            if (fTPCnSigma<-4  && EovP<0.5){
                fHistM02backg[iCent]->Fill(pt,m02);
                fHistM20backg[iCent]->Fill(pt,m20);
            }
        }

        if (pt>3) fHistM02EoverP[iCent]->Fill(m02,EovP);
        if (pt>3) fHistM20EoverP[iCent]->Fill(m20,EovP);
        
        
        if (pt<3 && fTPCnSigma>0  && fTPCnSigma<3 && fTOFnSigma>-2 && fTOFnSigma<2 && fITSnSigma>-2 && fITSnSigma<2)
            fEoverPsignalTPC[iCent]->Fill(pt,EovP);
        
        if (pt>=3 && fTPCnSigma>-1  && fTPCnSigma<3) fEoverPsignalTPC[iCent]->Fill(pt,EovP);
        
        if (fTPCnSigma>-1  && fTPCnSigma<3 && m02>0.02 && m20>0.02 && m02<0.3) fEoverPsignalTPCM02[iCent]->Fill(pt,EovP);
        if (fTPCnSigma>-5  && fTPCnSigma<-3.5) fEoverPbackg[iCent]->Fill(pt,EovP);
        
        //--- track accepted
        //        AliHFEpidObject hfetrack;
        //        hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
        //        hfetrack.SetRecTrack(track);
        //        hfetrack.SetPbPb();
        //        if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) continue;
        
        
        // checking centrality and event plane distributions for events with electron above the trigger threshold
        if (fTPCnSigma>=-1 && fTPCnSigma<3 && fEMCalnSigma>0 && fEMCalnSigma<3 && pt>=8 && GetCollisionCandidates()==AliVEvent::kEMCEGA && !IsSameEvent){
            fevPlaneV0AftThr[iCent]->Fill(evPlaneV0);
            fCentAftThr->Fill(cent);
            IsSameEvent=kTRUE;
        }
        
        Bool_t fFlagPhotonicElec = kFALSE;
        Bool_t fFlagPhotonicElecBCG = kFALSE;
        Double_t MCweight = 1.;
        Int_t iDecay = 0;
        
        Int_t partPDG = -99;
        Double_t partPt = -99.;
        Bool_t MChijing;
        
        if (pt<3){
            if ((fITSnSigma < (-1.*fSigmaITScut) ) || (fITSnSigma > fSigmaITScut)) continue;
            if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) continue;
            if ((fTPCnSigma<fSigmaTPCcut)  || (fTPCnSigma>3)) continue;
        }

        if (pt>=3 && fSystTOFcut){
            if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) continue;
        }
            
        
        for(Int_t i=0;i<3;i++) {
            if (dphi>=deltaPhiRange[i] && dphi<deltaPhiRange[i+1]){
                iDeltaphi=i;
                continue;
            }
        }
        
        if (pt>=3 && fSScut && (m02<0.02 || m02>fCutM02)) continue;
        if (pt>=3 && fSScut && (m20<0.02 || m20>fCutM20)) continue;
        
        if(!fMCarray){
            SelectPhotonicElectron(iTracks,track, fFlagPhotonicElec, fFlagPhotonicElecBCG,1,iCent,0,0,EovP,fTPCnSigma,evPlaneV0);
            
            if (pt<3 || (pt>=3 && fTPCnSigma>=-1 && fTPCnSigma<3 && EovP>0.8 && EovP<1.2))
                fInclElec[iCent]->Fill(pt,0.);
            
            Double_t corr[7]={static_cast<Double_t>(iCent),pt,fTPCnSigma,EovP,dphi,cosdphi,0.};
            fCorr->Fill(corr,wEvent);
        }
        
        // MC
        if(fMCarray){
            Int_t label = track->GetLabel();
            if(label!=0){
                fMCparticle = (AliAODMCParticle*) fMCarray->At(label);
                if(fMCparticle){
                    partPDG = fMCparticle->GetPdgCode();
                    partPt = fMCparticle->Pt();
                    
                    if (TMath::Abs(partPDG)==11){
                        fHistM02sig[iCent]->Fill(pt,m02);
                        fHistM20sig[iCent]->Fill(pt,m20);
                    }
                    
                    if (TMath::Abs(partPDG)!=11){
                        fHistM02backg[iCent]->Fill(pt,m02);
                        fHistM20backg[iCent]->Fill(pt,m20);
                    }
                    
                    
                    Bool_t iEnhance = kFALSE;
                    
                    if(label>=NpureMC)iEnhance = kTRUE;
                    
                    GetWeightAndDecay(fMCparticle,iCent,iDecay,MCweight);
                    
                    Double_t corr[7]={static_cast<Double_t>(iCent),pt,fTPCnSigma,EovP,dphi,cosdphi,static_cast<Double_t>(iDecay)};
                    fCorr->Fill(corr,MCweight);
                    
                    if (TMath::Abs(partPDG)!=11) continue;
                    
                    fInclElec[iCent]->Fill(pt,(Double_t)iDecay,MCweight);
                    if (pt<3 || (pt>=3 && fTPCnSigma>=-1 && fTPCnSigma<3 && EovP>0.8 && EovP<1.2))
                        fInclElec[iCent]->Fill(pt,(Double_t)iDecay);
                    
                    SelectPhotonicElectron(iTracks,track, fFlagPhotonicElec, fFlagPhotonicElecBCG,MCweight,iCent,iEnhance,iDecay,EovP,fTPCnSigma,evPlaneV0);
                    
                }// end particle
            }// end label
        }//end MC
        
        if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA ) fChargPartV2[iCent]->Fill(iPt,cosdphi,wEvent);
        else fChargPartV2[iCent]->Fill(iPt,cosdphi);
        
        if (clsE>0){
            if (iCent==0 && GetCollisionCandidates()!=AliVEvent::kEMCEGA )  fMtcPartV2[iCent]->Fill(iPt,cosdphi,wEvent);
            else  fMtcPartV2[iCent]->Fill(iPt,cosdphi);
        }
        
        if (fTPCnSigma>=-0.5 && fTPCnSigma<3) feTPCV2[iCent]->Fill(iPt,cosdphi);
        if (fTPCnSigma>=-0.5 && fTPCnSigma<3 && fEMCalnSigma>-1 && fEMCalnSigma<3) feV2[iCent]->Fill(iPt,cosdphi);
        
    }//end of track loop
    PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::UserCreateOutputObjects()
{
    //--- Check MC
    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
        fIsMC = kTRUE;
        printf("+++++ MC Data available");
    }
    //--------Initialize PID
    
    fPID->SetHasMCData(fIsMC);
    
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
        AliWarning("Cuts not available. Default cuts will be used");//same as in the config file (to be be removed)
        fCuts = new AliHFEcuts;
        fCuts->CreateStandardCuts();
        fCuts->SetMinNClustersTPC(fTPCnCut);
        fCuts->SetMinRatioTPCclusters(0.6);
        fCuts->SetMaxChi2perClusterTPC(3.5);
        fCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
        //hfecuts->SetMinNClustersITS(ITSncut); // it depends on pt
        fCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
        fCuts->SetCheckITSLayerStatus(kFALSE);
        fCuts->SetVertexRange(10.);
        fCuts->SetPtRange(1.5, 50);
        fCuts->SetMaxImpactParam(2.4,3.2); // radial, z
    }
    
    fCuts->SetAOD();
    fCuts->Initialize(fCFM);
    
    //---------Output Tlist
    fOutputList = new TList();
    fOutputList->SetOwner();
    fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    int nbin_v2 = 6;
    double bin_v2[7] = {1.5,2,2.5,3,4,6,8};
    
    fCentAftThr = new TH1F("fCentAftThr","Centrality for electron pt > 8 GeV/c",100,0,100) ;
    fOutputList->Add(fCentAftThr);
    
    fTrigger = new TH2F("fTrigger","",5,0,5,3,0,3) ;
    fOutputList->Add(fTrigger);
    
    fNoEvents = new TH1F("fNoEvents","",1,0,1) ;
    fOutputList->Add(fNoEvents);
    
    fTrkpt = new TH1F("fTrkpt","track pt",100,0,50);
    fOutputList->Add(fTrkpt);
    
    fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",20,0,10);
    fOutputList->Add(fTrackPtBefTrkCuts);
    
    fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",20,0,10);
    fOutputList->Add(fTrackPtAftTrkCuts);
    
    fChargedParticlePhi = new TH1F("fChargedParticlePhi","track phi",100,0,TMath::TwoPi());
    fOutputList->Add(fChargedParticlePhi);
    
    fElectronPhi = new TH1F("fElectronPhi","electron phi",100,0,TMath::TwoPi());
    fOutputList->Add(fElectronPhi);
    

    fCent = new TH1F("fCent","Centrality",100,0,100) ;
    fOutputList->Add(fCent);
    
    fCentAftFlt = new TH1F("fCentAftFlt","Centrality",100,0,100) ;
    fOutputList->Add(fCentAftFlt);
    
    fTPCsubEPres = new TH1F("fTPCsubEPres","TPC subevent plane resolution",100,-1,1);
    fOutputList->Add(fTPCsubEPres);
    
    
    int nbin_cpv2 = 28;
    double bin_cpv2[29] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,4.5,5.0,5.5,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,20.0};
    
    fcpV2_3040 = new TH2F("fcpV2_3040","",nbin_cpv2,bin_cpv2,20,-1,1) ;
    fOutputList->Add(fcpV2_3040);
    
    fcpV2_4050 = new TH2F("fcpV2_4050","",nbin_cpv2,bin_cpv2,20,-1,1) ;
    fOutputList->Add(fcpV2_4050);
    
    //iCent,pt,fTPCnSigma,EovP,dphi,cosphi,iDecay
    Int_t binsv2[7]=   {3,20,20,40,4        ,20,7};
    Double_t xminv2[7]={0,0 ,-5, 0,0          , -1,0};
    Double_t xmaxv2[7]={3,10, 5, 2,TMath::Pi(),  1,7};
    fCorr = new THnSparseD ("fCorr","Correlations",7,binsv2,xminv2,xmaxv2);
    fCorr->Sumw2();
    fOutputList->Add(fCorr);
    
    //iCent,pt,mass,fFlagLS,iEnhance,iDecay,EovP,fTPCnSigma,dphiPhotElec
    Int_t binsv3[9]=   {3,20,30 ,2,2,7,40,20,4};
    Double_t xminv3[9]={0,0 ,0  ,0,0,0,0 ,-5,0};
    Double_t xmaxv3[9]={3,10,0.3,2,2,7,2 , 5,TMath::Pi()};
    fElecMC = new THnSparseD ("fElecMC","MC",9,binsv3,xminv3,xmaxv3);
    fElecMC->Sumw2();
    fOutputList->Add(fElecMC);
    
    for(Int_t i=0; i<3; i++) {
        fevPlaneV0[i] = new TH1F(Form("fevPlaneV0%d",i),"V0 EP",100,0,TMath::Pi());
        fevPlaneV0[i]->Sumw2();
        fOutputList->Add(fevPlaneV0[i]);
        
        fevPlaneV0AftThr[i] = new TH1F(Form("fevPlaneV0AftThr%d",i),"V0 EP for electron pt > 8 GeV/c",100,0,TMath::Pi());
        fevPlaneV0AftThr[i]->Sumw2();
        fOutputList->Add(fevPlaneV0AftThr[i]);
        
        feTPCV2[i] = new TH2F(Form("feTPCV2%d",i), "", 8,0,8,100,-1,1);
        feTPCV2[i]->Sumw2();
        //fOutputList->Add(feTPCV2[i]);
        
        feV2[i] = new TH2F(Form("feV2%d",i), "", 8,0,8,100,-1,1);
        feV2[i]->Sumw2();
        //fOutputList->Add(feV2[i]);
        
        fChargPartV2[i] = new TH2F(Form("fChargPartV2%d",i), "", 8,0,8,100,-1,1);
        fChargPartV2[i]->Sumw2();
        fOutputList->Add(fChargPartV2[i]);
        
        fMtcPartV2[i] = new TH2F(Form("fMtcPartV2%d",i), "", 8,0,8,100,-1,1);
        fMtcPartV2[i]->Sumw2();
        fOutputList->Add(fMtcPartV2[i]);
        
        fEPres[i] = new TH2F(Form("fEPres%d",i), "", 100,0,100,100,-1,1);
        fEPres[i]->Sumw2();
        fOutputList->Add(fEPres[i]);
        
        fElecPtULSInvmassCut[i] = new TH2F(Form("fElecPtULSInvmassCut%d",i), "electron pt, ULS, invariant mass cut",nbin_v2,bin_v2,7,0,7);
        fElecPtULSInvmassCut[i]->Sumw2();
        fOutputList->Add(fElecPtULSInvmassCut[i]);
        
        fElecPtLSInvmassCut[i] = new TH2F(Form("fElecPtLSInvmassCut%d",i), "electron pt, LS, invariant mass cut",nbin_v2,bin_v2,7,0,7);
        fElecPtLSInvmassCut[i]->Sumw2();
        fOutputList->Add(fElecPtLSInvmassCut[i]);
        
        fElecPtInvmassCut[i] = new TH2F(Form("fElecPtInvmassCut%d",i), "electron pt, invariant mass cut",nbin_v2,bin_v2,7,0,7);
        fElecPtInvmassCut[i]->Sumw2();
        fOutputList->Add(fElecPtInvmassCut[i]);
        
        fInclElec[i] = new TH2F(Form("fInclElec%d",i), "inclusive electron pt", nbin_v2,bin_v2,7,0,7);
        fInclElec[i]->Sumw2();
        fOutputList->Add(fInclElec[i]);
        
        fInvmassLS[i] = new TH2F(Form("fInvmassLS%d",i), "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 500,0,0.5,100,0,50);
        fInvmassLS[i]->Sumw2();
        fOutputList->Add(fInvmassLS[i]);
        
        fInvmassULS[i] = new TH2F(Form("fInvmassULS%d",i), "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 500,0,0.5,100,0,50);
        fInvmassULS[i]->Sumw2();
        fOutputList->Add(fInvmassULS[i]);
        
        fOpeningAngleLS[i] = new TH2F(Form("fOpeningAngleLS%d",i),"Opening angle for LS pairs",100,0,1,100,0,50);
        fOpeningAngleLS[i]->Sumw2();
        fOutputList->Add(fOpeningAngleLS[i]);
        
        fOpeningAngleULS[i] = new TH2F(Form("fOpeningAngleULS%d",i),"Opening angle for ULS pairs",100,0,1,100,0,50);
        fOpeningAngleULS[i]->Sumw2();
        fOutputList->Add(fOpeningAngleULS[i]);
        
        fPi0Pt[i] = new TH2F(Form("fPi0Pt%d",i), "Pi0 weight",4,0,4,100,0,50);
        fPi0Pt[i]->Sumw2();
        fOutputList->Add(fPi0Pt[i]);
        
        fEtaPt[i] = new TH2F(Form("fEtaPt%d",i), "Eta weight",4,0,4,100,0,50);
        fEtaPt[i]->Sumw2();
        fOutputList->Add(fEtaPt[i]);
        
        fHistITSnSig[i] = new TH2F(Form("fHistITSnSig%d",i),Form("fHistITSnSig%d",i),50,0,5,200,-10,10);
        fOutputList->Add(fHistITSnSig[i]);
        
        fHistTOFnSig[i] = new TH2F(Form("fHistTOFnSig%d",i),Form("fHistTOFnSig%d",i),100,0,10,200,-10,10);
        fOutputList->Add(fHistTOFnSig[i]);
        
        fHistTPCnSig[i] = new TH2F(Form("fHistTPCnSig%d",i),Form("fHistTPCnSig%d",i),150,0,15,200,-10,10);
        fOutputList->Add(fHistTPCnSig[i]);
        
        fHistTPCnSigITScut[i] = new TH2F(Form("fHistTPCnSigITScut%d",i),Form("fHistTPCnSigITScut%d",i),150,0,15,200,-10,10);
        fOutputList->Add(fHistTPCnSigITScut[i]);
        
        fHistTPCnSigTOFcut[i] = new TH2F(Form("fHistTPCnSigTOFcut%d",i),Form("fHistTPCnSigTOFcut%d",i),150,0,15,200,-10,10);
        fOutputList->Add(fHistTPCnSigTOFcut[i]);
        
        fHistTPCnSigITSTOFcut[i] = new TH2F(Form("fHistTPCnSigITSTOFcut%d",i),Form("fHistTPCnSigITSTOFcut%d",i),150,0,15,200,-10,10);
        fOutputList->Add(fHistTPCnSigITSTOFcut[i]);
        
        fHistTPCnSigEop[i]  = new TH2F(Form("fHistTPCnSigEop%d",i),Form("fHistTPCnSigEop%d",i),100,-5,5,100,0,2);
        fOutputList->Add(fHistTPCnSigEop[i]);
        
        fHistTPCnSigEMCalnSig[i]  = new TH2F(Form("fHistTPCnSigEMCalnSig%d",i),Form("fHistTPCnSigEMCalnSig%d",i),100,-5,5,100,-5,5);
        fOutputList->Add(fHistTPCnSigEMCalnSig[i]);
        
        fHistM02sig[i]  = new TH2F(Form("fHistM02sig%d",i),Form("fHistM02sig%d",i),nbin_v2,bin_v2,200,0,2);
        //fOutputList->Add(fHistM02sig[i]);
        
        fHistM20sig[i]  = new TH2F(Form("fHistM20sig%d",i),Form("fHistM20sig%d",i),nbin_v2,bin_v2,200,0,2);
        //fOutputList->Add(fHistM20sig[i]);
        
        fHistM02backg[i]  = new TH2F(Form("fHistM02backg%d",i),Form("fHistM02backg%d",i),nbin_v2,bin_v2,200,0,2);
        //fOutputList->Add(fHistM02backg[i]);
        
        fHistM20backg[i]  = new TH2F(Form("fHistM20backg%d",i),Form("fHistM20backg%d",i),nbin_v2,bin_v2,200,0,2);
        //fOutputList->Add(fHistM20backg[i]);
        
        fHistM02EoverP[i]  = new TH2F(Form("fHistM02EoverP%d",i),Form("fHistM02EoverP%d",i),200,0,2,200,0,2);
        //fOutputList->Add(fHistM02EoverP[i]);
        
        fHistM20EoverP[i]  = new TH2F(Form("fHistM20EoverP%d",i),Form("fHistM20EoverP%d",i),200,0,2,200,0,2);
        //fOutputList->Add(fHistM20EoverP[i]);
        
        fEoverPsignalTPC[i]  = new TH2F(Form("fEoverPsignalTPC%d",i),Form("fEoverPsignalTPC%d",i),nbin_v2,bin_v2,40,0,2);
        fOutputList->Add(fEoverPsignalTPC[i]);
        
        fEoverPsignalTPCM02[i]  = new TH2F(Form("fEoverPsignalTPCM02%d",i),Form("fEoverPsignalTPCM02%d",i),nbin_v2,bin_v2,40,0,2);
        fOutputList->Add(fEoverPsignalTPCM02[i]);
        
        fEoverPbackg[i]  = new TH2F(Form("fEoverPbackg%d",i),Form("fEoverPbackg%d",i),nbin_v2,bin_v2,40,0,2);
        fOutputList->Add(fEoverPbackg[i]);
        
    }
    
    PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::Terminate(Option_t *)
{
    // Info("Terminate");
    AliAnalysisTaskSE::Terminate();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    // Check single track cuts for a given cut step
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetCos2DeltaPhi(Double_t phiA,Double_t phiB) const
{
    //Get cos[2(phi-psi_EP)] or cos[2(psi_subEP1 - psi_subEP2)]
    Double_t dPhi = TVector2::Phi_0_2pi(phiA - phiB);
    if(dPhi > TMath::Pi()) dPhi = dPhi - TMath::Pi();
    Double_t cos2DeltaPhi = TMath::Cos(2*dPhi);
    
    return cos2DeltaPhi;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetDeltaPhi(Double_t phiA,Double_t phiB) const
{
    //Get phi-psi_EP
    Double_t dPhi = TVector2::Phi_0_2pi(phiA - phiB);
    if(dPhi > TMath::Pi()) dPhi = dPhi - TMath::Pi();
    
    return dPhi;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetPi0weight(Double_t mcPi0pT, Int_t iCent) const
{
    //Get Pi0 weight
    double weight = 1.0;
    
    if (fWhichPeriod==2011){
        
        if (iCent==0){
            double parLowPt[4] = {0.00386062,0.913163,54.9096,84.0142};
            double parHighPt[4] = {0.02824,0.1246,3.56903,5.60296};
            
            if(mcPi0pT>0.0 && mcPi0pT<5.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=5.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==1){
            double parLowPt[4] = {0.000197581,0.960398,174.617,269.436};
            double parHighPt[4] = {0.0435973,0.0732613,3.43644,5.56708};
            
            if(mcPi0pT>0.0 && mcPi0pT<5.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=5.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==2){
            double parLowPt[4] = {0.00395183,0.905839,52.3325,78.9736};
            double parHighPt[4] = {0.0639772,0.0954623,3.21665,5.67225};
            
            if(mcPi0pT>0.0 && mcPi0pT<5.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=5.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }
    }
    
    
    if (fWhichPeriod==2015){
        
        if (iCent==0){
            double parLowPt[4] = {1.77623,0.204811,3.3,5.45959};
            double parHighPt[4] = {1.84871,0.509792,5.07536,7.53539};
            
            if(mcPi0pT>0.0 && mcPi0pT<4.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=4.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==1){
            double parLowPt[4] = {1.29732,0.208539,3.41144,5.58408};
            double parHighPt[4] = {1.36774,0.455898,4.75742,7.18348};
            
            if(mcPi0pT>0.0 && mcPi0pT<4.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=4.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==2 && !fWeightSyst){
            double parLowPt[4] = {1.8991,0.15029,2.60735,5.07192};
            double parHighPt[4] = {2.15436,0.600959,4.80365,7.87579};
            
            if(mcPi0pT>0.0 && mcPi0pT<4.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=4.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==2 && fWeightSyst){
            double parLowPt[4] = {0.937028,0.674846,9.02659,10.};
            double parHighPt[4] = {2.7883,0.,2.5684,5.63827};
            
            if(mcPi0pT>0.0 && mcPi0pT<4.0) weight = (parLowPt[0]*mcPi0pT)/TMath::Power(parLowPt[1]+mcPi0pT/parLowPt[2],parLowPt[3]);
            if(mcPi0pT>=4.0) weight = (parHighPt[0]*mcPi0pT)/TMath::Power(parHighPt[1]+mcPi0pT/parHighPt[2],parHighPt[3]);
        }

    }
    
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetEtaweight(Double_t mcEtapT, Int_t iCent) const
{
    //Get eta weight
    double weight = 1.0;
    
    if (fWhichPeriod==2011){
        if (iCent==0){
            double parLowPt[4] = {0.00218816,0.903496,52.9872,73.6404};
            double parHighPt[4] = {0.0742314,0.296077,3.33914,5.86723};
            
            if(mcEtapT>0.0 && mcEtapT<5.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=5.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==1){
            double parLowPt[4] = {0.00218739,0.889904,49.5944,66.8576};
            double parHighPt[4] = {0.124957,0.216147,3.09109,5.76838};
            
            if(mcEtapT>0.0 && mcEtapT<5.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=5.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==2){
            double parLowPt[4] = {0.00326269,0.911628,57.3255,78.6103};
            double parHighPt[4] = {0.134015,0.207723,3.00919,5.83206};
            
            if(mcEtapT>0.0 && mcEtapT<5.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=5.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }
    }
    
    if (fWhichPeriod==2015){
        if (iCent==0){
            double parLowPt[4] = {5.39949,0.294842,2.68408,5.51972};
            double parHighPt[4] = {7.042,0.639658,4.18394,7.48064};
            
            if(mcEtapT>0.0 && mcEtapT<4.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=4.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==1){
            double parLowPt[4] = {4.95436,0.335366,2.7807,5.79148};
            double parHighPt[4] = {5.33649,0.609717,4.08239,7.35917};
            
            if(mcEtapT>0.0 && mcEtapT<4.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=4.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==2 && !fWeightSyst){
            double parLowPt[4] = {3.76961,0.324609,2.60977,5.70101};
            double parHighPt[4] = {5.45548,0.668754,3.9863,7.63197};
            
            if(mcEtapT>0.0 && mcEtapT<4.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=4.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }
        if (iCent==2 && fWeightSyst){
            double parLowPt[4] = {2.26982,0.75242,7.12772,10.};
            double parHighPt[4] = {2.57403,0.,2.28527,5.659};
            
            if(mcEtapT>0.0 && mcEtapT<4.0) weight = (parLowPt[0]*mcEtapT)/TMath::Power(parLowPt[1]+mcEtapT/parLowPt[2],parLowPt[3]);
            if(mcEtapT>=4.0) weight = (parHighPt[0]*mcEtapT)/TMath::Power(parHighPt[1]+mcEtapT/parHighPt[2],parHighPt[3]);
        }

    }
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetSigmaEMCal(Double_t EoverP, Double_t pt, Int_t iCent) const
{
    //Get sigma for EMCal PID
    Double_t NumberOfSigmasEMCal = 99.;
    Double_t ptRange[9] = {1.5,2,2.5,3,4,6,8,10,13};
    
    if (iCent==0){
        Double_t mean[8]={1.04892,1.04471,1.04397,1.04715,1.04617,1.04147,1.05363,1.04902};
        Double_t sigma[8]={0.157249,0.143196,0.130118,0.118269,0.105798,0.100816,0.0910207,0.0953318};
        for(Int_t i=0;i<8;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
                continue;
            }
        }
    }
    if (iCent==1){
        Double_t mean[8]={1.01201,1.01064,1.01248,1.01728,1.02346,1.02177,1.04038,1.03314};
        Double_t sigma[8]={0.144614,0.126229,0.120568,0.107897,0.0919854,0.0920917,0.0859356,0.085302};
        for(Int_t i=0;i<8;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
                continue;
            }
        }
    }
    if (iCent==2){
        Double_t mean[8]={0.975778,0.975963,0.983835,0.988513,0.999726,1.00552,1.01144,1.00319};
        Double_t sigma[8]={0.130389,0.117007,0.10375,0.0971151,0.0893869,0.0873147,0.083138,0.0874688};
        for(Int_t i=0;i<8;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
                continue;
            }
        }
    }
    return NumberOfSigmasEMCal;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetSigmaEMCalMC(Double_t EoverP, Double_t pt, Int_t iCent) const
{
    //Get sigma for EMCal PID
    Double_t NumberOfSigmasEMCal = 99.;
    Double_t ptRange[9] = {1.5,2,2.5,3,4,6,8,10,13};
    
    if (iCent==0){
        Double_t mean[8]={1.01076,1.00735,1.00386,1.00281,1.00114,0.998282,0.995936,0.998286};
        Double_t sigma[8]={0.153704,0.137907,0.127886,0.115947,0.102482,0.0921989,0.0896079,0.0944837};
        for(Int_t i=0;i<8;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
                continue;
            }
        }
    }
    if (iCent==1){
        Double_t mean[8]={0.97531,0.973007,0.971888,0.972424,0.97437,0.976057,0.977703,0.984494};
        Double_t sigma[8]={0.132568,0.119308,0.107527,0.099176,0.0873851,0.0779302,0.0779114,0.0834648};
        for(Int_t i=0;i<8;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
                continue;
            }
        }
    }
    if (iCent==2){
        Double_t mean[8]={0.954379,0.952449,0.952901,0.955364,0.961415,0.965205,0.968959,0.976448};
        Double_t sigma[8]={0.120315,0.106597,0.0968691,0.0879189,0.0784124,0.0719245,0.0704888,0.080023};
        for(Int_t i=0;i<8;i++) {
            if (pt>=ptRange[i] && pt<ptRange[i+1]){
                NumberOfSigmasEMCal = (EoverP-mean[i])/sigma[i];
                continue;
            }
        }
    }
    return NumberOfSigmasEMCal;
}
//________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::GetWeightAndDecay(AliAODMCParticle *particle, Int_t iCent, Int_t &decay, Double_t &weight)
{
    //Get pi0/eta weight for MC with enchanced signal and decay channel
    Double_t w = 1.;
    Int_t d = 0;
    Int_t partPDG = particle->GetPdgCode();
    
    if (TMath::Abs(partPDG)==11){
        Int_t idMother = particle->GetMother();
        
        if (idMother>0){
            AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
            Int_t motherPDG = mother->GetPdgCode();
            Double_t motherPt = mother->Pt();
            
            Bool_t isMotherPrimary = IsPrimary(mother);
            Bool_t isMotherFromHF = IsFromHFdecay(mother);
            Bool_t isMotherFromLM = IsFromLMdecay(mother);
            
            if (motherPDG==111 && (isMotherPrimary || (!isMotherFromHF && !isMotherFromLM))){ // pi0 -> e
                d = 1;
                w = GetPi0weight(motherPt,iCent);
            }
            
            if (motherPDG==221  && (isMotherPrimary || (!isMotherFromHF && !isMotherFromLM))){ // eta -> e
                d = 2;
                w = GetEtaweight(motherPt,iCent);
            }
            
            //Int_t idSecondMother = particle->GetSecondMother();
            Int_t idSecondMother = mother->GetMother();
            
            if (idSecondMother>0){
                AliAODMCParticle* secondMother = (AliAODMCParticle*) fMCarray->At(idSecondMother);
                Int_t secondMotherPDG = secondMother->GetPdgCode();
                Double_t secondMotherPt = secondMother->Pt();
                
                Bool_t isSecondMotherPrimary = IsPrimary(secondMother);
                Bool_t isSecondMotherFromHF = IsFromHFdecay(secondMother);
                Bool_t isSecondMotherFromLM = IsFromLMdecay(secondMother);
                
                if (motherPDG==22 && secondMotherPDG==111 && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //pi0 -> g -> e
                    d = 3;
                    w = GetPi0weight(secondMotherPt,iCent);
                }
                
                if (motherPDG==22 && secondMotherPDG==221  && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //eta -> g -> e
                    d = 4;
                    w = GetEtaweight(secondMotherPt,iCent);
                }
                
                if (motherPDG==111 && secondMotherPDG==221  && (isSecondMotherPrimary || (!isSecondMotherFromHF && !isSecondMotherFromLM))){ //eta -> pi0 -> e
                    d = 5;
                    w = GetEtaweight(secondMotherPt,iCent);
                }
                
                Int_t idThirdMother = secondMother->GetMother();
                if (idThirdMother>0){
                    AliAODMCParticle* thirdMother = (AliAODMCParticle*) fMCarray->At(idThirdMother);
                    Int_t thirdMotherPDG = thirdMother->GetPdgCode();
                    Double_t thirdMotherPt = thirdMother->Pt();
                    
                    Bool_t isThirdMotherPrimary = IsPrimary(thirdMother);
                    Bool_t isThirdMotherFromHF = IsFromHFdecay(thirdMother);
                    Bool_t isThirdMotherFromLM = IsFromLMdecay(thirdMother);
                    
                    if (motherPDG==22 && secondMotherPDG==111 && thirdMotherPDG==221 && (isThirdMotherPrimary || (!isThirdMotherFromHF && !isThirdMotherFromLM))){//eta->pi0->g-> e
                        d = 6;
                        w = GetEtaweight(thirdMotherPt,iCent);
                    }
                }//third mother
            }//second mother
        }//mother
    }// if electron
    decay = d;
    weight = w;
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetCentWeight(Int_t centbin){
    // Get cebtrality weight for flattening (0-10%)
    Int_t wBin = centbin-1;
    if (wBin<0 || wBin>9) return 1;
    
    if (fWhichPeriod==2011){
        Double_t weightcent[] = {0.996425,0.987564,0.966774,0.966422,0.967739,0.991296,0.983746,0.991528,1.01627,1.15803};
        return weightcent[wBin];
    }
    
    if (fWhichPeriod==2015){
        Double_t weightcent[] = {1.0043,0.999405,0.999784,0.999727,0.999566,0.999633,0.999406,0.999317,0.99948,0.999407};
        return weightcent[wBin];
    }
    
    else
        return 1;
    
}
//_________________________________________
Double_t AliAnalysisTaskFlowTPCEMCalEP::GetEPweight(Int_t bin, Int_t iCent)
{
    //Get event plane weight for flattening
    Int_t wBin = bin-1;
    if (wBin<0 || wBin>99) return 1;
    
    if (fWhichPeriod==2011 && iCent==0){
        Double_t weightEP[] = {
            0.982991,0.988171,0.9899237,0.9914497,0.9906325,0.9956888,0.9972689,1.000973,1.002418,1.006948,
            1.007226,1.008336,1.01335,1.011154,1.018333,1.019898,1.026543,1.023092,1.028325,1.026844,
            1.031437,1.03014,1.031728,1.030307,1.037547,1.03471,1.03722,1.039466,1.037632,1.041682,
            1.042824,1.037494,1.046057,1.046622,1.042124,1.043161,1.040339,1.040997,1.043782,
            1.039092,1.039026,1.033509,1.035641,1.034528,1.031159,1.029701,1.033969,1.021809,1.02614,
            1.017396,1.017012,1.013525,1.012976,1.007164,1.006868,1.00653,0.9983816,0.9962069,
            0.9987208,0.9958153,0.9902154,0.9837839,0.9805614,0.9825041,0.9821056,0.9785275,
            0.9793774,0.9739373,0.9722809,0.9728094,0.972367,0.9687113,0.96755,0.9635185,0.9605392,
            0.9610214,0.9614648,0.9591571,0.9603319,0.9610102,0.9675955,0.9609205,0.9605896,
            0.9625102,0.9589448,0.9624427,0.966783,0.9632197,0.9626284,0.9706073,0.9693101,
            0.9717702,0.9703041,0.9747158,0.9741852,0.9755416,0.9798203,0.9797912,0.9790047,0.9802287};
        return weightEP[wBin];
    }
    
    if (fWhichPeriod==2015 && iCent==0){
        Double_t weightEP[] = {
            0.973914,0.966993,0.965672,0.963712,0.960045,0.956476,0.950277,0.946806,0.943463,0.940825,
            0.93614,0.931359,0.9305,0.925571,0.92469,0.925216,0.919276,0.920457,0.916642,0.914368,
            0.920269,0.918602,0.922843,0.919984,0.921338,0.927387,0.924402,0.925999,0.932253,0.93329,
            0.938218,0.944999,0.948569,0.949461,0.959602,0.972347,0.973659,0.979171,0.983525,0.994905,
            0.99968,1.00558,1.0115,1.01978,1.02969,1.03895,1.04168,1.04536,1.05951,1.06122,
            1.06038,1.07297,1.07644,1.07012,1.07884,1.08196,1.0809,1.08975,1.08716,1.08548,
            1.08538,1.09225,1.08233,1.08725,1.08868,1.08359,1.087,1.07784,1.0757,1.07094,
            1.07284,1.0776,1.06449,1.07286,1.05787,1.05861,1.04914,1.05653,1.05081,1.03951,
            1.04044,1.03734,1.02986,1.01915,1.02697,1.02088,1.01577,1.01029,1.0081,1.0124,
            1.00109,0.997445,0.995353,0.996374,0.991111,0.988454,0.983032,0.983172,0.978019,0.972857};
        return weightEP[wBin];
    }
    
    if (fWhichPeriod==2015 && iCent==1){
        Double_t weightEP[] = {
            0.996076,0.990221,0.994366,0.987992,0.986511,0.983266,0.983224,0.977671,0.974447,0.971244,
            0.965005,0.964951,0.963008,0.956021,0.957552,0.952022,0.95382,0.950101,0.953375,0.950352,
            0.952027,0.949728,0.952035,0.952742,0.95481,0.956931,0.95628,0.957993,0.958721,0.963276,
            0.964477,0.969764,0.966897,0.974835,0.973854,0.978222,0.984634,0.987663,0.995971,1.00027,
            1.00719,1.01316,1.01295,1.02477,1.02878,1.03206,1.04044,1.0428,1.04487,1.04613,
            1.05293,1.05438,1.05678,1.05647,1.05818,1.05301,1.05522,1.05411,1.04824,1.04775,
            1.04576,1.04149,1.03924,1.03667,1.03695,1.03444,1.02967,1.03172,1.03047,1.02956,
            1.02521,1.0241,1.01965,1.02047,1.01926,1.0222,1.00993,1.01188,1.01158,1.01139,
            1.00735,1.00739,1.00276,1.00367,0.995799,0.997767,0.9945,0.998059,1.00061,0.995756,
            0.998437,0.995709,1.00109,0.994699,0.997776,0.99737,1.00331,0.998245,0.995871,0.99299};
        return weightEP[wBin];
    }
    
    if (fWhichPeriod==2015 && iCent==2){
        Double_t weightEP[] = {
            0.988976,0.98572,0.987085,0.981873,0.979138,0.980245,0.97458,0.969915,0.966464,0.958588,
            0.958306,0.958248,0.950522,0.953947,0.952742,0.950097,0.950357,0.947756,0.949941,0.952441,
            0.949711,0.951997,0.953755,0.955147,0.956948,0.957596,0.957573,0.959893,0.964721,0.968515,
            0.966962,0.968772,0.968664,0.975869,0.983702,0.981753,0.985076,0.989997,1.00021,1.0082,
            1.00865,1.01677,1.02724,1.0269,1.03676,1.04114,1.04834,1.04418,1.05015,1.05704,
            1.05682,1.06196,1.05924,1.06517,1.05863,1.05905,1.05259,1.05768,1.0523,1.04626,
            1.04595,1.03879,1.03954,1.03899,1.03678,1.03524,1.03161,1.03299,1.03235,1.02935,
            1.03097,1.02623,1.02722,1.02345,1.0232,1.02147,1.02105,1.01631,1.01202,1.00845,
            1.00767,1.0061,1.00757,1.00302,0.999704,0.995877,0.993068,0.992803,0.993291,0.991993,
            0.99276,0.99069,0.993547,0.990369,0.992041,0.996005,0.995567,0.992074,0.992973,0.989776};
        return weightEP[wBin];
    }
    
    else
        return 1;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::RejectEvent(Double_t cent, Int_t centbin)
{
    // Reject randomly event in 0-10% in order to flatten the centrality distribution in MB events
    Int_t wBin = centbin-1;
    if (wBin<0 || wBin>9) return kFALSE;
    
    if (fWhichPeriod==2011){
        Double_t weight[] = {0.858984,0.853393,0.835121,0.834615,0.835851,0.855657,0.849417,0.856341,0.877473,1};
        
        Double_t centDigits=cent-(Int_t)(cent*100.)/100.;
        
        if(centDigits*100.>weight[wBin]) return kTRUE;
    }
    
    return kFALSE;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::RejectEventPlane(Double_t EP, Int_t EPbin)
{
    // Reject randomly event plane in 0-10% in order to flatten the EP distribution in MB events
    Int_t wBin = EPbin-1;
    if (wBin<0 || wBin>99) return kFALSE;
    
    if (fWhichPeriod==2011){
        Double_t weight[] = {
            0.939203,0.944153,0.945827,0.947285,0.946505,0.951336,0.952845,0.956385,0.957765,0.962094,
            0.962358,0.963419,0.968211,0.966111,0.972971,0.974466,0.980815,0.977518,0.982518,0.981103,
            0.985491,0.984252,0.985769,0.984412,0.991329,0.988619,0.991017,0.993163,0.99141,0.99528,
            0.996371,0.991279,0.99946,1,0.995702,0.996693,0.993997,0.994626,0.997287,0.992805,0.992743,
            0.987471,0.989508,0.988445,0.985226,0.983833,0.987911,0.976292,0.98043,0.972075,0.971709,
            0.968377,0.967853,0.9623,0.962017,0.961694,0.953908,0.951831,0.954233,0.951456,0.946106,
            0.939961,0.936882,0.938738,0.938357,0.934939,0.935751,0.930553,0.92897,0.929475,0.929053,
            0.92556,0.92445,0.920598,0.917752,0.918213,0.918636,0.916431,0.917554,0.918202,0.924494,
            0.918116,0.9178,0.919635,0.916228,0.91957,0.923717,0.920313,0.919748,0.927371,0.926132,
            0.928482,0.927082,0.931297,0.93079,0.932086,0.936174,0.936146,0.935395,0.936564};
        
        Double_t centDigits=EP-(Int_t)(EP*100.)/100.;
        
        if(centDigits*100.>weight[wBin]) return kTRUE;
    }
    
    return kFALSE;
}

//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsFromHFdecay(AliAODMCParticle *particle)
{
    // Check if the mother comes from heavy-flavour decays
    Bool_t isHFdecay = kFALSE;
    Int_t partPDG = particle->GetPdgCode();
    
    Int_t idMother = particle->GetMother();
    if (idMother>0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
        Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
        
        // c decays
        if((motherPDG % 1000) / 100 == 4) isHFdecay = kTRUE;
        if(motherPDG / 1000 == 4) isHFdecay = kTRUE;
        
        // b decays
        if((motherPDG % 1000) / 100 == 5) isHFdecay = kTRUE;
        if(motherPDG / 1000 == 5) isHFdecay = kTRUE;
        
        // all particles related to  HF
        if(motherPDG==4 || motherPDG==5 || motherPDG == 211 || motherPDG ==13 || motherPDG ==2112 || motherPDG ==130 || motherPDG ==3122 ||
           motherPDG ==310 || motherPDG ==3222 || motherPDG ==2212 || motherPDG ==3112 || motherPDG ==321 ||
           motherPDG ==11 || motherPDG ==3212 || motherPDG ==311 || motherPDG ==20213 || motherPDG ==3312 ||
           motherPDG ==3334 || motherPDG ==3324 || motherPDG ==3322 || motherPDG ==1000010020 || motherPDG ==15
           || motherPDG ==10323 || motherPDG ==2114 || motherPDG ==1000010030 || motherPDG ==2214 || motherPDG ==2224
           || motherPDG ==1114 || motherPDG == 2214 || motherPDG == 3114 || motherPDG ==3224 || motherPDG ==3124
           || motherPDG ==3314 || motherPDG ==10323 || motherPDG == 3214) isHFdecay = kTRUE;
        
    }
    
    return isHFdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsFromLMdecay(AliAODMCParticle *particle)
{
    // Check if the mother comes from light-meson decays
    Bool_t isLMdecay = kFALSE;
    Int_t partPDG = particle->GetPdgCode();
    
    Int_t idMother = particle->GetMother();
    if (idMother>0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fMCarray->At(idMother);
        Int_t motherPDG = TMath::Abs(mother->GetPdgCode());
        
        if(motherPDG == 111 || motherPDG == 221 || motherPDG==223 || motherPDG==333 || motherPDG==331 ||
           motherPDG==113 || motherPDG==213 || motherPDG==313 || motherPDG==323) isLMdecay = kTRUE;
    }
    
    return isLMdecay;
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::IsPrimary(AliAODMCParticle *particle)
{
    // Check if the particle is primary
    Bool_t isprimary = kFALSE;
    
    Int_t idMother = particle->GetMother();
    if (idMother==-1) isprimary = kTRUE;
    //if (particle->IsPrimary()) isprimary = kTRUE;
    
    return isprimary;
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::SelectPhotonicElectron(Int_t iTracks,AliAODTrack *track,Bool_t &fFlagPhotonicElec,
                                                           Bool_t &fFlagPhotonicElecBCG,Double_t weight, Int_t iCent, Int_t iEnhance, Int_t iDecay, Double_t EovP, Double_t fTPCnSigma, Double_t evPlaneV0)
{
    //Identify non-heavy flavour electrons using Invariant mass method
    
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    
    Bool_t flagPhotonicElec = kFALSE;
    Bool_t flagPhotonicElecBCG = kFALSE;
    
    Int_t ntracks = -999;
    if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender)  ntracks = fTracks_tender->GetEntries();
    
    
    for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
        
        if(jTracks==iTracks) continue;
        
        AliVParticle* Vassotrack = 0x0;
        if(!fUseTender) Vassotrack  = fVevent->GetTrack(jTracks);
        if(fUseTender) Vassotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
        
        if (!Vassotrack) {
            printf("ERROR: Could not receive associated track %d\n", jTracks);
            continue;
        }
        AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(Vassotrack);
        if(!trackAsso) continue;
        if(!trackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(trackAsso->Pt() < fAssPtCut) continue;
        if(TMath::Abs(trackAsso->Eta())>0.9) continue;
        if(trackAsso->GetTPCNcls() < fAssTPCnCut) continue;
        if (fAssITSrefitCut && !(trackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
        if(!(trackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
        
        Double_t pt=-999., ptAsso=-999., nTPCsigmaAsso=-999.;
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t openingAngle = -999., mass=999., width = -999;
        Int_t chargeAsso = 0, charge = 0, pdgAsso = 0;
        
        nTPCsigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);        
        pt = track->Pt();
        ptAsso = trackAsso->Pt();
        chargeAsso = trackAsso->Charge();
        charge = track->Charge();
        
        if(TMath::Abs(nTPCsigmaAsso)>3) continue;

        
        Double_t dphiPhotElec = -9, phi=0;
        phi = trackAsso->Pt();
        dphiPhotElec = GetDeltaPhi(phi,evPlaneV0);
        
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1(*track, fPDGe1);
        AliKFParticle ge2(*trackAsso, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        openingAngle = ge1.GetAngle(ge2);
        
        if(fFlagLS) fOpeningAngleLS[iCent]->Fill(openingAngle,ptAsso);
        if(fFlagULS) fOpeningAngleULS[iCent]->Fill(openingAngle,ptAsso);
        
        //if(openingAngle > fOpeningAngleCut) continue;
        
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);
        
        Double_t elecMC[9]={(Double_t)iCent,pt,mass,(Double_t)fFlagLS,(Double_t)iEnhance,(Double_t)iDecay, EovP, fTPCnSigma,dphiPhotElec};
        fElecMC->Fill(elecMC,weight);
        
        if(fFlagLS) fInvmassLS[iCent]->Fill(mass,pt,weight);
        if(fFlagULS) fInvmassULS[iCent]->Fill(mass,pt,weight);
        
        if(mass<fInvmassCut) fElecPtInvmassCut[iCent]->Fill(pt,iDecay,weight);
        if(mass<fInvmassCut && fFlagULS) fElecPtULSInvmassCut[iCent]->Fill(pt,iDecay,weight);
        if(mass<fInvmassCut && fFlagLS) fElecPtLSInvmassCut[iCent]->Fill(pt,iDecay,weight);
    }
    fFlagPhotonicElec = flagPhotonicElec;
    fFlagPhotonicElecBCG = flagPhotonicElecBCG;
    
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::InclElecTrackCuts(AliAODTrack *ietrack)
{
    
    // quality track cuts for inclusive electron
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    
    if(!ietrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
    if(TMath::Abs(ietrack->Eta())>0.7) return kFALSE;
    if(ietrack->GetTPCNcls() < 100) return kFALSE;
    if(ietrack->GetITSNcls() < 5) return kFALSE;
    if(!ietrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;
    if(!ietrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
    if(!(ietrack->HasPointOnITSLayer(0) || ietrack->HasPointOnITSLayer(1))) return kFALSE;
    if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;
    
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[nVertices];
    Int_t nMotherKink = 0;
    for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
        AliAODVertex *vertex = fAOD->GetVertex(ivertex);
        if(!vertex) return kFALSE;
        if(vertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
            if(!mother) return kFALSE;
            Int_t idmother = mother->GetID();
            listofmotherkink[nMotherKink] = idmother;
            nMotherKink++;
        }
    }
    
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMotherKink; kinkmother++) {
        if(ietrack->GetID() == listofmotherkink[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    Double_t d0z0[2]={-999,-999}, cov[3];
    
    if(ietrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > 2.4 || TMath::Abs(d0z0[1]) > 3.2) return kFALSE;
    
    
    return kTRUE;
    
}
//_________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalEP::AssElecTrackCuts(AliAODTrack *aetrack)
{
    // quality track cuts for associate tracks
    
    // check these cuts:
    //    fAssTrackCuts->SetRequireSigmaToVertex(kTRUE);
    //    fAssTrackCuts->SetDCAToVertex2D(kTRUE);
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    
    
    if(aetrack->Pt() < fAssPtCut) return kFALSE;
    if(!aetrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
    if(TMath::Abs(aetrack->Eta())>0.9) return kFALSE;
    if(aetrack->GetTPCNcls() < fAssTPCnCut) return kFALSE;
    if(aetrack->GetITSNcls() < 2) return kFALSE;
    if (fAssITSrefitCut && !aetrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;
    if(!aetrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
    if(aetrack->GetTPCFoundFraction() < 0.6) return kFALSE;
    
    
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[nVertices];
    Int_t nMotherKink = 0;
    for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
        AliAODVertex *vertex = fAOD->GetVertex(ivertex);
        if(!vertex) return kFALSE;
        if(vertex->GetType()==AliAODVertex::kKink) {
            AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
            if(!mother) return kFALSE;
            Int_t idmother = mother->GetID();
            listofmotherkink[nMotherKink] = idmother;
            nMotherKink++;
        }
    }
    
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMotherKink; kinkmother++) {
        if(aetrack->GetID() == listofmotherkink[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    Double_t d0z0[2]={-999,-999}, cov[3];
    
    if(aetrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > 2.4 || TMath::Abs(d0z0[1]) > 3.2) return kFALSE;
    
    
    return kTRUE;
    
}
//_________________________________________
const AliQnCorrectionsQnVector *AliAnalysisTaskFlowTPCEMCalEP::GetQnVectorFromList(const TList *list,
                                                                                   const char *subdetector,
                                                                                   const char *expectedstep,
                                                                                   const char *altstep)
{
    AliQnCorrectionsQnVector *theQnVector = NULL;
    
    TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
    if (pQvecList != NULL) {
        /* the detector is present */
        if (TString(expectedstep).EqualTo("latest"))
            theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
        else
            theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expectedstep);
        
        if (theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) {
            /* the Qn vector for the expected step was not there */
            if (TString(altstep).EqualTo("latest"))
                theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
            else
                theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altstep);
        }
    }
    if (theQnVector != NULL) {
        /* check the Qn vector quality */
        if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
        /* not good quality, discarded */
            theQnVector = NULL;
    }
    return theQnVector;
}
//_________________________________________
void AliAnalysisTaskFlowTPCEMCalEP::InitParameters()
{
    // Init parameters
    //    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //    fTrackCuts->SetRequireTPCRefit(kTRUE);
    //    fTrackCuts->SetRequireITSRefit(kTRUE);
    //    fTrackCuts->SetEtaRange(-0.7,0.7);
    //    fTrackCuts->SetRequireSigmaToVertex(kTRUE);
    //    fTrackCuts->SetMaxChi2PerClusterTPC(3.5);
    //    fTrackCuts->SetMinNClustersTPC(100);
    //    fTrackCuts->SetPtRange(0.5,100);
}
