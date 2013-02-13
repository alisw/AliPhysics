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
//  Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation       //
//  Non-Photonic Electron identified with Invariant mass              //
//  analysis methos in function  SelectPhotonicElectron               //
//  DeltaPhi calculated in function  ElectronHadCorrel                // 
//                                                                    //
//  Author: Deepa Thomas (Utrecht University)                         //
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

#include "AliAnalysisTaskElecHadronCorrel.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"

#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliCentralitySelectionTask.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

//#include "AliEventPoolManager.h"

#include "AliCentrality.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

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
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

ClassImp(AliAnalysisTaskElecHadronCorrel)
//ClassImp(AliehDPhiBasicParticle)  
//________________________________________________________________________
  AliAnalysisTaskElecHadronCorrel::AliAnalysisTaskElecHadronCorrel(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fVevent(0)  
  ,fESD(0)
  ,fAOD(0)
  ,fGeom(0)
  ,fOutputList(0)
  ,fTrackCuts1(new AliESDtrackCuts)
  ,fTrackCuts2(new AliESDtrackCuts)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)
  ,fPIDqa(0)	       
  ,fInvmassCut(0.01)	
  ,fCentrality(0)
  ,fCentralityMin(0)
  ,fCentralityMax(0)
  ,fkCentralityMethod(0)  
//  ,fPoolMgr(0x0)  
  ,fNoEvents(0)
//  ,fTrkpt(0)
//  ,fTrkEovPBef(0)	 
//  ,fTrkEovPBefHad(0)	 
 // ,fdEdxBef(0)	 
  ,fSemiIncElecDphi(0) 	
  ,fSemiIncElecDphi1(0) 	
  ,fSemiIncElecDphi2(0) 	
  ,fSemiIncElecDphi3(0) 	
  ,fSemiIncElecDphi4(0) 	
  ,fPhotElecDphi(0)  	
  ,fPhotElecDphi1(0)  	
  ,fPhotElecDphi2(0)  	
  ,fPhotElecDphi3(0)  	
  ,fPhotElecDphi4(0)  	
  ,fInclusiveElecDphi(0)  	
  ,fInclusiveElecDphi1(0)  	
  ,fInclusiveElecDphi2(0)  	
  ,fInclusiveElecDphi3(0)  	
  ,fInclusiveElecDphi4(0)  	
  ,fInclusiveElecDphiEtaFS(0)
  ,fInclusiveElecDphiEtaFS1(0)
  ,fInclusiveElecDphiEtaFS2(0)
  ,fInclusiveElecDphiEtaFS3(0)
  ,fInclusiveElecDphiEtaFS4(0)
  ,fDphiULSMassLow(0)	
  ,fDphiULSMassLow1(0)	
  ,fDphiULSMassLow2(0)	
  ,fDphiULSMassLow3(0)	
  ,fDphiULSMassLow4(0)	
  ,fDphiLSMassLow(0)
  ,fDphiLSMassLow1(0)
  ,fDphiLSMassLow2(0)
  ,fDphiLSMassLow3(0)
  ,fDphiLSMassLow4(0)
  ,fDphiULSMassLowNoPartner(0)   
  ,fDphiULSMassLowNoPartner1(0)   
  ,fDphiULSMassLowNoPartner2(0)   
  ,fDphiULSMassLowNoPartner3(0)   
  ,fDphiULSMassLowNoPartner4(0)   
  ,fDphiLSMassLowNoPartner(0)
  ,fDphiLSMassLowNoPartner1(0)
  ,fDphiLSMassLowNoPartner2(0)
  ,fDphiLSMassLowNoPartner3(0)
  ,fDphiLSMassLowNoPartner4(0)
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fInclusiveElecPt(0)
  ,fULSElecPt(0)
  ,fLSElecPt(0)  
 // ,fTrackPtBefTrkCuts(0)	 
 // ,fTrackPtAftTrkCuts(0)
 // ,fTPCnsigma(0)
 // ,fNCellv1(0)
 // ,fClsEv1(0)
 // ,fNClusv1(0)
  ,fInvmassLS1(0)       
 // ,fInvmassLS2(0)       
 // ,fInvmassLS3(0)       
 // ,fInvmassLS4(0)       
 // ,fInvmassLS5(0)       
  ,fInvmassULS1(0)
 // ,fInvmassULS2(0)
 // ,fInvmassULS3(0)
 // ,fInvmassULS4(0)
 // ,fInvmassULS5(0)
  ,fcentrality(0)     
  ,fElecPhi(0)  
  ,fElecPhiTPChalf(0)
  ,fElecPhiPt(0)  
//  ,fElecPhiTPC(0)  
//  ,fElecPhiTPCEovP(0)  
  ,fHadronPhi(0)  
  ,fHadronPhiTPChalf(0)  
  ,fHadronPhiPt(0)  
/*  ,fTrackHFEcuts(0)
  ,fTrakPhiSPD1(0)
  ,fTrakPhiSPD2(0)
  ,fTrakPhiSPDOr(0)
  ,fTrakPhiSPDAnd(0)
  ,fTrackHFEcutsITS(0)  
*/
/*  ,fNoMixedEvents(0)
  ,fMixStat(0)       
  ,fMixStat1(0)        
  ,fMixedIncElecDphi(0)  
  ,fMixedIncElecDphi1(0)  
  ,fMixedIncElecDphi2(0)  
  ,fMixedPhotElecDphi(0)
  ,fMixedPhotElecDphi1(0)
  ,fMixedPhotElecDphi2(0)
  ,fMixedSemiIncElecDphi(0)  
  ,fMixedSemiIncElecDphi1(0)  
  ,fMixedSemiIncElecDphi2(0)  
  ,fMixedDphiULSMassLow(0)  
  ,fMixedDphiULSMassLow1(0)  
  ,fMixedDphiULSMassLow2(0)  
  ,fMixedDphiLSMassLow(0)  
  ,fMixedDphiLSMassLow1(0)  
  ,fMixedDphiLSMassLow2(0)  
*/   
  ,fHadronPt(0)  
  ,fCentralityPass(0)
  ,fCentralityNoPass(0)
  ,fHadronDphi(0)
  ,fHadronDphi1(0)
  ,fHadronDphi2(0)
  ,fHadronDphi3(0)
  ,fHadronDphi4(0)
  ,fPiPt(0)  
  ,fHadronDphiNoSS(0)
  ,fHadronDphiNoSS1(0)
  ,fHadronDphiNoSS2(0)
  ,fHadronDphiNoSS3(0)
  ,fHadronDphiNoSS4(0)
  ,fPiPtNoSS(0)
  ,fSparseElectron(0)  
    ,fvalueElectron(0)   
{
  //Named constructor

  fPID = new AliHFEpid("hfePid");
  fvalueElectron = new Double_t[6];

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskElecHadronCorrel::AliAnalysisTaskElecHadronCorrel() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
  ,fVevent(0)  
  ,fESD(0)
  ,fAOD(0)
  ,fGeom(0)  
  ,fOutputList(0)
  ,fTrackCuts1(new AliESDtrackCuts)
  ,fTrackCuts2(new AliESDtrackCuts)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)       
  ,fPIDqa(0)	       
  ,fInvmassCut(0.01)	
  ,fCentrality(0)
  ,fCentralityMin(0)
  ,fCentralityMax(0)
  ,fkCentralityMethod(0)  
//  ,fPoolMgr(0x0)    
  ,fNoEvents(0)
//  ,fTrkpt(0)
//  ,fTrkEovPBef(0)	 
//  ,fTrkEovPBefHad(0)	 
//  ,fdEdxBef(0)	 
  ,fSemiIncElecDphi(0) 	
  ,fSemiIncElecDphi1(0) 	
  ,fSemiIncElecDphi2(0) 	
  ,fSemiIncElecDphi3(0) 	
  ,fSemiIncElecDphi4(0) 	
  ,fPhotElecDphi(0)  	
  ,fPhotElecDphi1(0)  	
  ,fPhotElecDphi2(0)  	
  ,fPhotElecDphi3(0)  	
  ,fPhotElecDphi4(0)  	
  ,fInclusiveElecDphi(0)  	
  ,fInclusiveElecDphi1(0)  	
  ,fInclusiveElecDphi2(0)  	
  ,fInclusiveElecDphi3(0)  	
  ,fInclusiveElecDphi4(0)  	
  ,fInclusiveElecDphiEtaFS(0)
  ,fInclusiveElecDphiEtaFS1(0)
  ,fInclusiveElecDphiEtaFS2(0)
  ,fInclusiveElecDphiEtaFS3(0)
  ,fInclusiveElecDphiEtaFS4(0)
  ,fDphiULSMassLow(0)	
  ,fDphiULSMassLow1(0)	
  ,fDphiULSMassLow2(0)	
  ,fDphiULSMassLow3(0)	
  ,fDphiULSMassLow4(0)	
  ,fDphiLSMassLow(0)
  ,fDphiLSMassLow1(0)
  ,fDphiLSMassLow2(0)
  ,fDphiLSMassLow3(0)
  ,fDphiLSMassLow4(0)
  ,fDphiULSMassLowNoPartner(0)   
  ,fDphiULSMassLowNoPartner1(0)   
  ,fDphiULSMassLowNoPartner2(0)   
  ,fDphiULSMassLowNoPartner3(0)   
  ,fDphiULSMassLowNoPartner4(0)   
  ,fDphiLSMassLowNoPartner(0)
  ,fDphiLSMassLowNoPartner1(0)
  ,fDphiLSMassLowNoPartner2(0)
  ,fDphiLSMassLowNoPartner3(0)
  ,fDphiLSMassLowNoPartner4(0)
  ,fPhotoElecPt(0)
  ,fSemiInclElecPt(0)
  ,fInclusiveElecPt(0)
  ,fULSElecPt(0)
  ,fLSElecPt(0)  
//  ,fTrackPtBefTrkCuts(0)	 
//  ,fTrackPtAftTrkCuts(0)	 	  
//  ,fTPCnsigma(0)	
//  ,fNCellv1(0)  
//  ,fClsEv1(0)
//  ,fNClusv1(0)
  ,fInvmassLS1(0)   
//  ,fInvmassLS2(0)   
//  ,fInvmassLS3(0)   
//  ,fInvmassLS4(0)   
//  ,fInvmassLS5(0)   
  ,fInvmassULS1(0)  
//  ,fInvmassULS2(0)  
//  ,fInvmassULS3(0)  
//  ,fInvmassULS4(0)  
//  ,fInvmassULS5(0)  
  ,fcentrality(0)     
  ,fElecPhi(0)
  ,fElecPhiTPChalf(0)  
  ,fElecPhiPt(0)
//  ,fElecPhiTPC(0)
//  ,fElecPhiTPCEovP(0)  
  ,fHadronPhi(0)
  ,fHadronPhiTPChalf(0)
  ,fHadronPhiPt(0)
/*  ,fTrackHFEcuts(0)
  ,fTrakPhiSPD1(0)
  ,fTrakPhiSPD2(0)
  ,fTrakPhiSPDOr(0)
  ,fTrakPhiSPDAnd(0)
  ,fTrackHFEcutsITS(0)  
*/
/*  ,fNoMixedEvents(0)
  ,fMixStat(0)      
  ,fMixStat1(0)     
  ,fMixedIncElecDphi(0)  
  ,fMixedIncElecDphi1(0)  
  ,fMixedIncElecDphi2(0)  
  ,fMixedPhotElecDphi(0)
  ,fMixedPhotElecDphi1(0)
  ,fMixedPhotElecDphi2(0)
  ,fMixedSemiIncElecDphi(0)
  ,fMixedSemiIncElecDphi1(0)
  ,fMixedSemiIncElecDphi2(0)
  ,fMixedDphiULSMassLow(0) 
  ,fMixedDphiULSMassLow1(0) 
  ,fMixedDphiULSMassLow2(0) 
  ,fMixedDphiLSMassLow(0)      
  ,fMixedDphiLSMassLow1(0)      
  ,fMixedDphiLSMassLow2(0)      
*/
  ,fHadronPt(0)  
  ,fCentralityPass(0)
  ,fCentralityNoPass(0)
  ,fHadronDphi(0)
  ,fHadronDphi1(0)
  ,fHadronDphi2(0)
  ,fHadronDphi3(0)
  ,fHadronDphi4(0)
  ,fPiPt(0)
  ,fHadronDphiNoSS(0)
  ,fHadronDphiNoSS1(0)
  ,fHadronDphiNoSS2(0)
  ,fHadronDphiNoSS3(0)
  ,fHadronDphiNoSS4(0)
  ,fPiPtNoSS(0)
  ,fSparseElectron(0)  
    ,fvalueElectron(0)  
{
  //Default constructor
  fPID = new AliHFEpid("hfePid");
  fvalueElectron = new Double_t[6];

  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(3, TTree::Class());
}
//_________________________________________

AliAnalysisTaskElecHadronCorrel::~AliAnalysisTaskElecHadronCorrel()
{
  //Destructor 

  delete fOutputList;
  delete fGeom;
  delete fPID;
  delete fCFM;
  delete fPIDqa;
  delete fTrackCuts1;
  delete fTrackCuts2;
  delete fSparseElectron;
  delete []fvalueElectron;
}
//_________________________________________

void AliAnalysisTaskElecHadronCorrel::UserExec(Option_t*)
{
  //Main loop
  //Called for each event

  // create pointer to event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if(!(fESD || fAOD)){
    printf("ERROR: fESD & fAOD not available\n");
    return;
  }
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent) {
    printf("ERROR: fVEvent not available\n");
    return;
  }

  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }

  if(!fPID->IsInitialized()){ 
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
 
    if(IsAODanalysis())fPID->InitializePID(fAOD->GetRunNumber());
    else fPID->InitializePID(fESD->GetRunNumber());
  }

  // trigger selection
  if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kCentral))) return;

  // centrality selection 
  SetCentralityParameters(0., 7., "V0M");
  Bool_t pass = kFALSE; 
  CheckCentrality(fVevent,pass);
  if(!pass)return;

  Int_t fNOtrks =  fVevent->GetNumberOfTracks();
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();

  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();

  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(0);

  if(fNOtrks<2) return;

  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
     pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
   }

   fPID->SetPIDResponse(pidResponse);

   fCFM->SetRecEventInfo(fVevent);

   /*
   //Event mixing
   AliEventPool* pool = fPoolMgr->GetEventPool(centvalue, pVtxZ); // Get the buffer associated with the current centrality and z-vtx
   if (!pool)
     AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centvalue, pVtxZ));     
*/

   // Look for kink mother for AOD
   Double_t *listofmotherkink =0;
   Int_t numberofvertices = 0, numberofmotherkink = 0;
   if(IsAODanalysis()){
     numberofvertices = fAOD->GetNumberOfVertices();
     listofmotherkink = new Double_t[numberofvertices];
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
   }

   // Track loop 
   for (Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) {
     AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
     if (!Vtrack) {
       printf("ERROR: Could not receive track %d\n", iTracks);
       continue;
     }
     AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
     AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
     AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
     
     if(IsAODanalysis())
       if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

     if(track->Pt()<1) continue;

     // fTrackPtBefTrkCuts->Fill(track->Pt());		

     // RecKine: ITSTPC cuts  
     if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;

     // Reject kink mother
     if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
       if(IsAODanalysis()){
         Bool_t kinkmotherpass = kTRUE;
         for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
           if(track->GetID() == listofmotherkink[kinkmother]) {
             kinkmotherpass = kFALSE;
             continue;
           }
         }
         if(!kinkmotherpass) continue;
       }
       else{
         if(etrack->GetKinkIndex(0) != 0) continue;
       }
     }
     // RecPrim
     //     if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue; //gives warning for AOD, so not using

     // HFE cuts: TPC PID cleanup
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

     // fTrackHFEcuts->Fill(track->Phi());

     // HFEcuts: ITS layers cuts
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;

     //     fTrackHFEcutsITS->Fill(track->Phi());
     //fTrackPtAftTrkCuts->Fill(track->Pt());		

     Double_t fClsE = -999, p = -999, fEovP=-999, pt = -999, dEdx=-999, fTPCnSigma=0;
     pt = track->Pt();
     p = track->P();
     dEdx = track->GetTPCsignal();
     fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;

     //TPC electron phi
    // if(fTPCnSigma >= -2 && fTPCnSigma <= 2){
       //       fElecPhiTPC->Fill(track->Phi());
    // }

     //eta cut (-0.7,0.7)
     if(track->Eta() < -0.7 || track->Eta() > 0.7) continue;

     // Track extrapolation to EMCAL
     Int_t fClsId = track->GetEMCALcluster();
     if(fClsId <0) continue;
     AliVCluster *cluster = fVevent->GetCaloCluster(fClsId);
     if(!cluster->IsEMCAL()) continue;
     if(TMath::Abs(cluster->GetTrackDx())>0.05 || TMath::Abs(cluster->GetTrackDz())>0.05) continue;    
     //     fdEdxBef->Fill(p,dEdx);
     //fTPCnsigma->Fill(p,fTPCnSigma);

     //     fTrkpt->Fill(pt);
     fClsE = cluster->E();
     fEovP = fClsE/p;
        
     //--------THnsparse---------
     fvalueElectron[0] = pt;
     fvalueElectron[1] = fTPCnSigma;
     fvalueElectron[2] = fEovP;
     fvalueElectron[3] = cluster->GetM20();
     fvalueElectron[4] = cluster->GetM02();
     fvalueElectron[5] = cluster->GetDispersion();

     fSparseElectron->Fill(fvalueElectron);

     //----------------
     //hadron E/p and Dphi distribution with shower shape cuts
     if(((fTPCnSigma > -10) && (fTPCnSigma < -3.5)) && ((cluster->GetM20()>0.03) && (cluster->GetM20()<0.3)) &&  ((cluster->GetM02()>0.03) && (cluster->GetM02()<0.5)) && ((cluster->GetDispersion()<1))){
       ElectronHadCorrel(iTracks, track, fHadronDphi, fHadronDphi1,fHadronDphi2,fHadronDphi3,fHadronDphi4);
       fPiPt->Fill(pt);
     }

     //hadron E/p and Dphi distribution without shower shape cuts
     if((fTPCnSigma > -10) && (fTPCnSigma < -3.5)){
       ElectronHadCorrel(iTracks, track, fHadronDphiNoSS, fHadronDphiNoSS1,fHadronDphiNoSS2,fHadronDphiNoSS3,fHadronDphiNoSS4);
       fPiPtNoSS->Fill(pt);
     }

     //Electron id with TPC and E/p
     if(fTPCnSigma >= -2 && fTPCnSigma <= 2 && fEovP >= 0.8 && fEovP <=1.2) {
       //   fElecPhiTPCEovP->Fill(track->Phi());

       //Electron id with shower shape  
       if(cluster->GetM20()>0.03 && cluster->GetM20()<0.3 && cluster->GetM02()>0.03 && cluster->GetM02()< 0.5 && cluster->GetDispersion()<1){

         fElecPhi->Fill(track->Phi());
         fElecPhiPt->Fill(track->Phi(),track->Pt());
         if (track->Eta() >0 && track->Eta() <0.7) fElecPhiTPChalf->Fill(track->Phi());
         
         HadronInfo(iTracks);

         Bool_t fFlagPhotonicElec = kFALSE;
         // select photonic electron
         SelectPhotonicElectron(iTracks,track,fFlagPhotonicElec);

         //Inclusive electron-hadron correlation
         ElectronHadCorrel(iTracks, track, fInclusiveElecDphi, fInclusiveElecDphi1,fInclusiveElecDphi2,fInclusiveElecDphi3,fInclusiveElecDphi4);
         fInclusiveElecPt->Fill(pt);
         //       MixedEvent(track,fMixedIncElecDphi, fMixedIncElecDphi1,fMixedIncElecDphi2);

         //Inclusive electron-hadron correlation
         ElectronHadCorrelEtaFarSide(iTracks, track, fInclusiveElecDphiEtaFS, fInclusiveElecDphiEtaFS1,fInclusiveElecDphiEtaFS2,fInclusiveElecDphiEtaFS3,fInclusiveElecDphiEtaFS4);
         //       MixedEvent(track,fMixedIncElecDphi, fMixedIncElecDphi1,fMixedIncElecDphi2);

         // photonic electron
         if(fFlagPhotonicElec){
           //Electron hadron correlation
           ElectronHadCorrel(iTracks, track, fPhotElecDphi,fPhotElecDphi1,fPhotElecDphi2,fPhotElecDphi3,fPhotElecDphi4);
           fPhotoElecPt->Fill(pt);
           //       MixedEvent(track,fMixedPhotElecDphi, fMixedPhotElecDphi1,fMixedPhotElecDphi2);
         }

         // Semi inclusive electron 
         if(!fFlagPhotonicElec){
           //Electron hadron correlation
           ElectronHadCorrel(iTracks, track, fSemiIncElecDphi, fSemiIncElecDphi1,fSemiIncElecDphi2,fSemiIncElecDphi3,fSemiIncElecDphi4);
           fSemiInclElecPt->Fill(pt);
           //        MixedEvent(track,fMixedSemiIncElecDphi,fMixedSemiIncElecDphi1,fMixedSemiIncElecDphi2);
         }
         
       }
     }
   }
   
/*   //EMC clusters  
   Int_t clsNo = fVevent->GetNumberOfCaloClusters();
   fNClusv1->Fill(clsNo); 
   for(Int_t iclus=0; iclus<clsNo ; iclus++){ 
     AliVCluster* clus = fVevent->GetCaloCluster(iclus);
     if(!clus->IsEMCAL()) continue; 
     fNCellv1->Fill(clus->GetNCells());
     fClsEv1->Fill(clus->E());  
   }
*/
/*
   TObjArray* tracksClone = CloneAndReduceTrackList();
   tracksClone->SetOwner();
   pool->UpdatePool(tracksClone);
   */
   delete listofmotherkink;
   PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::UserCreateOutputObjects()
{
  //Create histograms

  AliDebug(3, "Creating Output Objects");
  // Automatic determination of the analysis mode
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis();
  } else {
    SetESDAnalysis();
  }
  printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");

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

  if(IsAODanalysis()) fCuts->SetAOD(); 
  fCuts->Initialize(fCFM);
  /*
  //Mixed event initialising
  Int_t trackDepth = 2000;
  Int_t poolsize   = 1000;

  Int_t nCentralityBins  = 5;
  Double_t CentralityBins[] = {0,2,4,6,8,10};

  Int_t nZvtxBins  = 4;
  Double_t vertexBins[] = {-10,-5,0,5,10};

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, (Double_t*) CentralityBins, nZvtxBins, (Double_t*) vertexBins);
*/
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->Add(fPIDqa->MakeList("PIDQA"));

  fNoEvents = new TH1F("fNoEvents","",1,0,1) ;
  fOutputList->Add(fNoEvents);

  fcentrality = new TH1F("fcentrality","centrality", 100,0,100);
  fOutputList->Add(fcentrality);

//  fTrkpt = new TH1F("fTrkpt","track pt",1000,0,50);
//  fOutputList->Add(fTrkpt);

//  fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",1000,0,50);
//  fOutputList->Add(fTrackPtBefTrkCuts);

//  fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",1000,0,50);
//  fOutputList->Add(fTrackPtAftTrkCuts);

//  fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma",1000,0,50,200,-10,10);
//  fOutputList->Add(fTPCnsigma);
  
//  fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",1000,0,50,100,0,2);
//  fOutputList->Add(fTrkEovPBef);

//  fTrkEovPBefHad = new TH2F("fTrkEovPBefHad","track E/p for TPCnsig < 3.5",1000,0,50,100,0,2);
//  fOutputList->Add(fTrkEovPBefHad);

//  fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",1000,0,50,150,0,150);
//  fOutputList->Add(fdEdxBef);

  fElecPhi = new TH1F("fElecPhi", "Electron phi",1000,0,6.28);
  fOutputList->Add(fElecPhi);

  fElecPhiPt = new TH2F("fElecPhiPt", "Electron phi vs pt; Electron phi; pt (GeV/c)",1000,0,6.28,1000,0,100);
  fOutputList->Add(fElecPhiPt);

  fElecPhiTPChalf = new TH1F("fElecPhiTPChalf", "Electron phi for 0<eta<0.7",1000,0,6.28);                            
  fOutputList->Add(fElecPhiTPChalf);

/*  fElecPhiTPC = new TH1F("fElecPhiTPC", "Electron phi after TPC cut",1000,0,6.28);
  fOutputList->Add(fElecPhiTPC);

  fElecPhiTPCEovP = new TH1F("fElecPhiTPCEovP", "Electron phi after TPC and E/p cut",1000,0,6.28);
  fOutputList->Add(fElecPhiTPCEovP);
*/
  fHadronPhi = new TH1F("fHadronPhi", "Hadron phi",1000,0,6.28);
  fOutputList->Add(fHadronPhi);

  fHadronPhiTPChalf = new TH1F("fHadronPhiTPChalf", "Hadron phi for 0<eta<0.9",1000,0,6.28);                          
  fOutputList->Add(fHadronPhiTPChalf);

  fHadronPhiPt = new TH2F("fHadronPhiPt", "Hadron phi vs pt; hadron phi; pt (GeV/c)",1000,0,6.28,1000,0,100);
  fOutputList->Add(fHadronPhiPt);

/*
  fTrackHFEcuts = new TH1F("fTrackHFEcuts","Track phi for HFE cuts",1000,0,6.28);
  fOutputList->Add(fTrackHFEcuts);

  fTrakPhiSPD1 = new TH1F("fTrakPhiSPD1","Track phi for hit in SPD layer 1",1000,0,6.28);
  fOutputList->Add(fTrakPhiSPD1);

  fTrakPhiSPD2 = new TH1F("fTrakPhiSPD2","Track phi for hit in SPD layer 2",1000,0,6.28);
  fOutputList->Add(fTrakPhiSPD2);

  fTrakPhiSPDOr = new TH1F("fTrakPhiSPDOr","Track phi for hit in any SPD layer",1000,0,6.28);
  fOutputList->Add(fTrakPhiSPDOr);

  fTrakPhiSPDAnd = new TH1F("fTrakPhiSPDAnd","Track phi for hit in both SPD layer",1000,0,6.28);
  fOutputList->Add(fTrakPhiSPDAnd);

  fTrackHFEcutsITS = new TH1F("fTrackHFEcutsITS","Track phi for HFE cuts + ITS HFE cuts",1000,0,6.28);
  fOutputList->Add(fTrackHFEcutsITS);
*/
  fSemiIncElecDphi = new TH2F("fSemiIncElecDphi", "Semi Inclusive elec-had Dphi correlation",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi);

  fSemiIncElecDphi1 = new TH2F("fSemiIncElecDphi1", "Semi Inclusive elec-had Dphi correlation for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi1);

  fSemiIncElecDphi2 = new TH2F("fSemiIncElecDphi2", "Semi Inclusive elec-had Dphi correlation for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi2);

  fSemiIncElecDphi3 = new TH2F("fSemiIncElecDphi3", "Semi Inclusive elec-had Dphi correlation for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi3);

  fSemiIncElecDphi4 = new TH2F("fSemiIncElecDphi4", "Semi Inclusive elec-had Dphi correlation for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi4);

  fPhotElecDphi = new TH2F("fPhotElecDphi", "Photon elec-had Dphi correlation",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi);

  fPhotElecDphi1 = new TH2F("fPhotElecDphi1", "Photon elec-had Dphi correlation for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi1);

  fPhotElecDphi2 = new TH2F("fPhotElecDphi2", "Photon elec-had Dphi correlation for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi2);

  fPhotElecDphi3 = new TH2F("fPhotElecDphi3", "Photon elec-had Dphi correlation for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi3);

  fPhotElecDphi4 = new TH2F("fPhotElecDphi4", "Photon elec-had Dphi correlation for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi4);

  fInclusiveElecDphi = new TH2F("fInclusiveElecDphi", "Inclusive elec-had Dphi correlation",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi);

  fInclusiveElecDphi1 = new TH2F("fInclusiveElecDphi1", "Inclusive elec-had Dphi correlation for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi1);

  fInclusiveElecDphi2 = new TH2F("fInclusiveElecDphi2", "Inclusive elec-had Dphi correlation for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi2);

  fInclusiveElecDphi3 = new TH2F("fInclusiveElecDphi3", "Inclusive elec-had Dphi correlation for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi3);

  fInclusiveElecDphi4 = new TH2F("fInclusiveElecDphi4", "Inclusive elec-had Dphi correlation for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi4);

  fInclusiveElecDphiEtaFS = new TH2F("fInclusiveElecDphiEtaFS", "Inclusive elec-had Dphi correlation (hadron 1<eta<1.6)",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphiEtaFS);

  fInclusiveElecDphiEtaFS1 = new TH2F("fInclusiveElecDphiEtaFS1", "Inclusive elec-had Dphi correlation for 1<pt^{asso}<3 (hadron 1<eta<1.6)",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphiEtaFS1);

  fInclusiveElecDphiEtaFS2 = new TH2F("fInclusiveElecDphiEtaFS2", "Inclusive elec-had Dphi correlation for 3<pt^{asso}<5 (hadron 1<eta<1.6)",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphiEtaFS2);

  fInclusiveElecDphiEtaFS3 = new TH2F("fInclusiveElecDphiEtaFS3", "Inclusive elec-had Dphi correlation for 5<pt^{asso}<7 (hadron 1<eta<1.6)",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphiEtaFS3);

  fInclusiveElecDphiEtaFS4 = new TH2F("fInclusiveElecDphiEtaFS4", "Inclusive elec-had Dphi correlation for 7<pt^{asso}<9 (hadron 1<eta<1.6)",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphiEtaFS4);

  fDphiULSMassLow = new TH2F("fDphiULSMassLow", "e-h Dphi ULS, mass<cut",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow);

  fDphiULSMassLow1 = new TH2F("fDphiULSMassLow1", "e-h Dphi ULS, mass<cut for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow1);

  fDphiULSMassLow2 = new TH2F("fDphiULSMassLow2", "e-h Dphi ULS, mass<cut for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow2);

  fDphiULSMassLow3 = new TH2F("fDphiULSMassLow3", "e-h Dphi ULS, mass<cut for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow3);

  fDphiULSMassLow4 = new TH2F("fDphiULSMassLow4", "e-h Dphi ULS, mass<cut for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow4);

  fDphiLSMassLow = new TH2F("fDphiLSMassLow", "e-h Dphi LS, mass<cut",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow);

  fDphiLSMassLow1 = new TH2F("fDphiLSMassLow1", "e-h Dphi LS, mass<cut for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow1);

  fDphiLSMassLow2 = new TH2F("fDphiLSMassLow2", "e-h Dphi LS, mass<cut for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow2);

  fDphiLSMassLow3 = new TH2F("fDphiLSMassLow3", "e-h Dphi LS, mass<cut for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow3);

  fDphiLSMassLow4 = new TH2F("fDphiLSMassLow4", "e-h Dphi LS, mass<cut for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow4);

  fDphiULSMassLowNoPartner = new TH2F("fDphiULSMassLowNoPartner", "e-h Dphi ULS with no partner, mass<mass cut,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner);

  fDphiULSMassLowNoPartner1 = new TH2F("fDphiULSMassLowNoPartner1", "e-h Dphi ULS with no partner, mass<mass cut for 1<pt^{asso}<3,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner1);

  fDphiULSMassLowNoPartner2 = new TH2F("fDphiULSMassLowNoPartner2", "e-h Dphi ULS with no partner, mass<mass cut for 3<pt^{asso}<5,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner2);

  fDphiULSMassLowNoPartner3 = new TH2F("fDphiULSMassLowNoPartner3", "e-h Dphi ULS with no partner, mass<mass cut for 5<pt^{asso}<7,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner3);

  fDphiULSMassLowNoPartner4 = new TH2F("fDphiULSMassLowNoPartner4", "e-h Dphi ULS with no partner, mass<mass cut for 7<pt^{asso}<9,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner4);

  fDphiLSMassLowNoPartner = new TH2F("fDphiLSMassLowNoPartner", "e-h Dphi LS with no partner, mass<mass cut",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner);

  fDphiLSMassLowNoPartner1 = new TH2F("fDphiLSMassLowNoPartner1", "e-h Dphi LS with no partner, mass<mass cut for 1<pt^{asso}<3,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner1);                                            

  fDphiLSMassLowNoPartner2 = new TH2F("fDphiLSMassLowNoPartner2", "e-h Dphi LS with no partner, mass<mass cut for 3<pt^{asso}<5,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner2);

  fDphiLSMassLowNoPartner3 = new TH2F("fDphiLSMassLowNoPartner3", "e-h Dphi LS with no partner, mass<mass cut for 5<pt^{asso}<7,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner3);

  fDphiLSMassLowNoPartner4 = new TH2F("fDphiLSMassLowNoPartner4", "e-h Dphi LS with no partner, mass<mass cut for 7<pt^{asso}<9,",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner4);

  fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",1000,0,100);
  fOutputList->Add(fPhotoElecPt);

  fSemiInclElecPt = new TH1F("fSemiInclElecPt", "Semi-inclusive electron pt",1000,0,100);
  fOutputList->Add(fSemiInclElecPt);

  fInclusiveElecPt = new TH1F("fInclElecPt", "Inclusive electron pt",1000,0,100);
  fOutputList->Add(fInclusiveElecPt);

  fULSElecPt = new TH1F("fULSElecPt", "ULS electron pt",1000,0,100);
  fOutputList->Add(fULSElecPt);

  fLSElecPt = new TH1F("fLSElecPt", "LS electron pt",1000,0,100);
  fOutputList->Add(fLSElecPt);

/*  fNCellv1 = new TH1F("fNCellv1","Ncell in clus (v1); NCell; count",100,0,100) ;
  fOutputList->Add(fNCellv1);

  fClsEv1 = new TH1F("fClsEv1", "Clus E(v1); Cls E; count",1000,0,100); 
  fOutputList->Add(fClsEv1); 

  fNClusv1 = new TH1F("fNClusv1","Nclus in event (v1); NClus; count",500,0,500) ; 
  fOutputList->Add(fNClusv1);
*/
  fInvmassLS1 = new TH1F("fInvmassLS1", "Inv mass of LS (e,e) for pt^{e}>2; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS1);

  fInvmassULS1 = new TH1F("fInvmassULS1", "Inv mass of ULS (e,e) for pt^{e}>2; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS1);
/*
  fInvmassLS2 = new TH1F("fInvmassLS2", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS2);

  fInvmassULS2 = new TH1F("fInvmassULS2", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS2);

  fInvmassLS3 = new TH1F("fInvmassLS3", "Inv mass of LS (e,e) for pt^{e}>2; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS3);                                                          

  fInvmassULS3 = new TH1F("fInvmassULS3", "Inv mass of ULS (e,e) for pt^{e}>2; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS3);

  fInvmassLS4 = new TH1F("fInvmassLS4", "Inv mass of LS (e,e) for pt^{e}>3; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS4);  

  fInvmassULS4 = new TH1F("fInvmassULS4", "Inv mass of ULS (e,e) for pt^{e}>3; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS4);

  fInvmassLS5 = new TH1F("fInvmassLS5", "Inv mass of LS (e,e) for pt^{e}>4; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS5);  

  fInvmassULS5 = new TH1F("fInvmassULS5", "Inv mass of ULS (e,e) for pt^{e}>4; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS5);
*/
/*  fNoMixedEvents = new TH1F("fNoMixedEvents","",1,0,1) ;
  fOutputList->Add(fNoMixedEvents);

  fMixStat = new TH2F("fMixStat","no of events in pool  vs Centrality;Nevent in pool;Centrality",200,0,200,5,0,10);
  fOutputList->Add(fMixStat);                                                             

  fMixStat1 = new TH2F("fMixStat1","no of events in pool  vs zvtx;Nevents in pool;zvtx",200,0,200,4,-10,10);
  fOutputList->Add(fMixStat1);

  fMixedIncElecDphi = new TH2F("fMixedIncElecDphi", "Mixed event - Inclusive elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedIncElecDphi);

  fMixedIncElecDphi1 = new TH2F("fMixedIncElecDphi1", "Mixed event - Inclusive elec-had Dphi correlation 1<pt<3",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedIncElecDphi1);

  fMixedIncElecDphi2 = new TH2F("fMixedIncElecDphi2", "Mixed event - Inclusive elec-had Dphi correlation 3<pt<5",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedIncElecDphi2);

  fMixedSemiIncElecDphi = new TH2F("fMixedSemiIncElecDphi", "Mixed event - Semi Inclusive elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedSemiIncElecDphi);

  fMixedSemiIncElecDphi1 = new TH2F("fMixedSemiIncElecDphi1", "Mixed event - Semi Inclusive elec-had Dphi correlation 1<pt<3",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedSemiIncElecDphi1);

  fMixedSemiIncElecDphi2 = new TH2F("fMixedSemiIncElecDphi2", "Mixed event - Semi Inclusive elec-had Dphi correlation 3<pt<5",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedSemiIncElecDphi2);

  fMixedPhotElecDphi = new TH2F("fMixedPhotElecDphi", "Mixed event - Photo elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedPhotElecDphi);

  fMixedPhotElecDphi1 = new TH2F("fMixedPhotElecDphi1", "Mixed event - Photo elec-had Dphi correlation 1<pt<3",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedPhotElecDphi1);

  fMixedPhotElecDphi2 = new TH2F("fMixedPhotElecDphi2", "Mixed event - Photo elec-had Dphi correlation 3<pt<5",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedPhotElecDphi2);

  fMixedDphiULSMassLow = new TH2F("fMixedDphiULSMassLow", "Mixed event - ULS mass < cut elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiULSMassLow);

  fMixedDphiULSMassLow1 = new TH2F("fMixedDphiULSMassLow1", "Mixed event - ULS mass < cut elec-had Dphi correlation 1<pt<3",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiULSMassLow1);

  fMixedDphiULSMassLow2 = new TH2F("fMixedDphiULSMassLow2", "Mixed event - ULS mass < cut elec-had Dphi correlation 3<pt<5",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiULSMassLow2);

  fMixedDphiLSMassLow = new TH2F("fMixedDphiLSMassLow", "Mixed event - LS mass < cut elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiLSMassLow);

  fMixedDphiLSMassLow1 = new TH2F("fMixedDphiLSMassLow1", "Mixed event - LS mass < cut elec-had Dphi correlation 1<pt<3",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiLSMassLow1);

  fMixedDphiLSMassLow2 = new TH2F("fMixedDphiLSMassLow2", "Mixed event - LS mass < cut elec-had Dphi correlation 3<pt<5",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiLSMassLow2);
*/
  fHadronPt = new TH1F("fHadronPt","hadron pt distribution",1000,0,100);
  fOutputList->Add(fHadronPt);

  fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass", 101, -1, 100);
  fOutputList->Add(fCentralityPass);

  fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass", 101, -1, 100);
  fOutputList->Add(fCentralityNoPass);

  fHadronDphi = new TH2F("fHadronDphi", "Hadron-had Dphi correlation",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphi);

  fHadronDphi1 = new TH2F("fHadronDphi1", "Hadron-had Dphi correlation for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphi1);

  fHadronDphi2 = new TH2F("fHadronDphi2", "Hadron-had Dphi correlation for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphi2);

  fHadronDphi3 = new TH2F("fHadronDphi3", "Hadron-had Dphi correlation for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphi3);

  fHadronDphi4 = new TH2F("fHadronDphi4", "Hadron-had Dphi correlation for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphi4);

  fPiPt = new TH1F("fPiPt","Pi (-10 <TPC nsig < -3.5) pt distribution",1000,0,100);
  fOutputList->Add(fPiPt);

  fHadronDphiNoSS = new TH2F("fHadronDphiNoSS", "Hadron-had Dphi correlation (NoSS cuts)",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphiNoSS);

  fHadronDphiNoSS1 = new TH2F("fHadronDphiNoSS1", "Hadron-had Dphi correlation (NoSS cuts) for 1<pt^{asso}<3",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphiNoSS1);

  fHadronDphiNoSS2 = new TH2F("fHadronDphiNoSS2", "Hadron-had Dphi correlation (NoSS cuts) for 3<pt^{asso}<5",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphiNoSS2);

  fHadronDphiNoSS3 = new TH2F("fHadronDphiNoSS3", "Hadron-had Dphi correlation (NoSS cuts) for 5<pt^{asso}<7",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphiNoSS3);

  fHadronDphiNoSS4 = new TH2F("fHadronDphiNoSS4", "Hadron-had Dphi correlation (NoSS cuts) for 7<pt^{asso}<9",200,0,20,64,-1.57,4.71);
  fOutputList->Add(fHadronDphiNoSS4);

  fPiPtNoSS = new TH1F("fPiPtNoSS","Pi (-10 <TPC nsig < -3.5) (NoSS cuts) pt distribution",1000,0,100);
  fOutputList->Add(fPiPtNoSS);

  Int_t binsv1[6]={500,200,100,100,100,100}; //pt, TPCnsig, E/p, M20, M02, dispersion 
  Double_t xminv1[6]={0,-10,0,0,0,0};
  Double_t xmaxv1[6]={25,10,2,2,2,2};
  fSparseElectron = new THnSparseD ("Electron","Electron",6,binsv1,xminv1,xmaxv1);
  fOutputList->Add(fSparseElectron);

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskElecHadronCorrel::Terminate(Option_t *)
{
  // Info("Terminate");
  AliAnalysisTaskSE::Terminate();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskElecHadronCorrel::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
  // Check single track cuts for a given cut step
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  return kTRUE;
}

//_________________________________________
void AliAnalysisTaskElecHadronCorrel::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec)
{
  //Identify non-heavy flavour electrons using Invariant mass method

  fTrackCuts1->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts1->SetRequireTPCRefit(kTRUE);
  fTrackCuts1->SetRequireITSRefit(kTRUE);
  fTrackCuts1->SetEtaRange(-0.9,0.9);
  fTrackCuts1->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts1->SetMaxChi2PerClusterTPC(4);
  fTrackCuts1->SetMinNClustersTPC(80);
  fTrackCuts1->SetMaxDCAToVertexZ(3.2);
  fTrackCuts1->SetMaxDCAToVertexXY(2.4);
  fTrackCuts1->SetDCAToVertex2D(kTRUE);

  Bool_t flagPhotonicElec = kFALSE;

  for(Int_t jTracks = 0; jTracks<fVevent->GetNumberOfTracks(); jTracks++){
    AliVParticle* VtrackAsso = fVevent->GetTrack(jTracks);
    if (!VtrackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }

    AliVTrack *trackAsso = dynamic_cast<AliVTrack*>(VtrackAsso);

    //track cuts applied
    if(IsAODanalysis()) { 
      AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
      if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if(atrackAsso->GetTPCNcls() < 80) continue;
      if((!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
    }
    else{
      AliESDtrack *etrackAsso = dynamic_cast<AliESDtrack*>(VtrackAsso);
      if(!fTrackCuts1->AcceptTrack(etrackAsso)) continue;
    }

    if(jTracks==itrack) continue;

    Double_t dEdxAsso = -999., ptAsso=-999., nsigma=-999.0;
    Double_t mass=-999., width = -999;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;

    dEdxAsso = trackAsso->GetTPCsignal();
    nsigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();

    if(ptAsso <0.3) continue;
    if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
    if(nsigma < -3 || nsigma > 3) continue;

    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;

    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;

    AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
    AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);

    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

    Int_t MassCorrect;
    MassCorrect = recg.GetMass(mass,width);

    if(fFlagLS) {
      if(track->Pt()>2)fInvmassLS1->Fill(mass);
     // if(track->Pt()> 1) fInvmassLS2->Fill(mass);
     // if(track->Pt()>2) fInvmassLS3->Fill(mass);
     // if(track->Pt()>3) fInvmassLS4->Fill(mass);
     // if(track->Pt()>4) fInvmassLS5->Fill(mass);
    }
    if(fFlagULS) {
      if(track->Pt()>2)fInvmassULS1->Fill(mass);
      //if(track->Pt() >1) fInvmassULS2->Fill(mass);
      //if(track->Pt() >2) fInvmassULS3->Fill(mass);
      //if(track->Pt() >3) fInvmassULS4->Fill(mass);
      //if(track->Pt() >4) fInvmassULS5->Fill(mass);
    }

    if(mass<fInvmassCut){
      if(fFlagULS)
      {
        ElectronHadCorrel(itrack,track,fDphiULSMassLow, fDphiULSMassLow1,fDphiULSMassLow2,fDphiULSMassLow3,fDphiULSMassLow4);
        fULSElecPt->Fill(track->Pt());
        //      MixedEvent(track,fMixedDphiULSMassLow,fMixedDphiULSMassLow1,fMixedDphiULSMassLow2);
      }
      if(fFlagLS)
      {
        ElectronHadCorrel(itrack,track,fDphiLSMassLow,fDphiLSMassLow1,fDphiLSMassLow2,fDphiLSMassLow3,fDphiLSMassLow4);
        fLSElecPt->Fill(track->Pt());
        //     MixedEvent(track,fMixedDphiLSMassLow,fMixedDphiLSMassLow1,fMixedDphiLSMassLow2);
      }
      if(fFlagLS) ElectronHadCorrelNoPartner(itrack,jTracks,track,fDphiLSMassLowNoPartner, fDphiLSMassLowNoPartner1,fDphiLSMassLowNoPartner2,fDphiLSMassLowNoPartner3,fDphiLSMassLowNoPartner4);
      if(fFlagULS) ElectronHadCorrelNoPartner(itrack,jTracks,track,fDphiULSMassLowNoPartner, fDphiULSMassLowNoPartner1,fDphiULSMassLowNoPartner2,fDphiULSMassLowNoPartner3,fDphiULSMassLowNoPartner4);
    }

    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
    //   }

}
fFlagPhotonicElec = flagPhotonicElec;
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrel(Int_t itrack, AliVTrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2,TH2F *DphiPt3,TH2F *DphiPt4)
{
  //Construct Delta Phi between electrons and hadrons

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts2->SetRequireTPCRefit(kTRUE);
  fTrackCuts2->SetRequireITSRefit(kTRUE);
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts2->SetMaxChi2PerClusterTPC(4);
  fTrackCuts2->SetMinNClustersTPC(80);
  fTrackCuts2->SetMaxDCAToVertexZ(3.2);
  fTrackCuts2->SetMaxDCAToVertexXY(2.4);
  fTrackCuts2->SetDCAToVertex2D(kTRUE);

  for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){
    AliVParticle* VtrackHad = fVevent->GetTrack(ktracks);
    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }

    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);

    if(IsAODanalysis()) {
      AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
      if(!atrackHad->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if((!(atrackHad->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrackHad->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(atrackHad->GetTPCNcls() < 80) continue; 
    }
    else{   
      AliESDtrack *etrackHad = dynamic_cast<AliESDtrack*>(VtrackHad); 
      if(!fTrackCuts2->AcceptTrack(etrackHad)) continue; 
    }
    
    if(ktracks == itrack) continue; //do not select the same electron

    Double_t ptHad= -999, pHad=-999., dEdxHad = -999;
    Double_t ptEle = -999;
    Double_t phiEle = -999, phiHad = -999, Dphi = -999;
    Double_t pi = 3.14;

    dEdxHad = trackHad->GetTPCsignal();
    ptHad = trackHad->Pt();
    pHad = trackHad->P();
    ptEle = track->Pt();

    //    if(ptHad <2) continue;
    if(ptHad > ptEle) continue;
    if(trackHad->Eta()<-0.9 || trackHad->Eta()>0.9) continue;

    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    if(ptHad>2) DphiPt->Fill(ptEle,Dphi);
    if(ptHad>2 && ptHad<4) DphiPt1->Fill(ptEle,Dphi);
    if(ptHad>4 && ptHad<6) DphiPt2->Fill(ptEle,Dphi);
    if(ptHad>6 && ptHad<8) DphiPt3->Fill(ptEle,Dphi);
    if(ptHad>8 && ptHad<10) DphiPt4->Fill(ptEle,Dphi);

  }
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrelNoPartner(Int_t itrack,Int_t jtrack, AliVTrack *track, TH2F *DphiPtNew, TH2F *DphiPtNew1,TH2F *DphiPtNew2,TH2F *DphiPtNew3,TH2F *DphiPtNew4)
{
  //Construct Delta Phi between electrons and hadrons for electrons from invariant mass calculation excluding associated track

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts2->SetRequireTPCRefit(kTRUE);
  fTrackCuts2->SetRequireITSRefit(kTRUE);
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts2->SetMaxChi2PerClusterTPC(4);
  fTrackCuts2->SetMinNClustersTPC(80);
  fTrackCuts2->SetMaxDCAToVertexZ(3.2);
  fTrackCuts2->SetMaxDCAToVertexXY(2.4);
  fTrackCuts2->SetDCAToVertex2D(kTRUE);

  for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){
    AliVParticle* VtrackHad = fVevent->GetTrack(ktracks);
    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);

    if(IsAODanalysis()) {
      AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
      if(!atrackHad->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if((!(atrackHad->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrackHad->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(atrackHad->GetTPCNcls() < 80) continue; 
    }
    else{   
      AliESDtrack *etrackHad = dynamic_cast<AliESDtrack*>(VtrackHad); 
      if(!fTrackCuts2->AcceptTrack(etrackHad)) continue; 
    }

    if(ktracks == itrack || ktracks == jtrack) continue; //do not select the same electron and associated track from inv mass cal

    Double_t ptHad= -999, pHad=-999., dEdxHad = -999;
    Double_t ptEle = -999;
    Double_t phiEle = -999, phiHad = -999, Dphi = -999;
    Double_t pi = 3.14;

    dEdxHad = trackHad->GetTPCsignal();
    ptHad = trackHad->Pt();
    pHad = trackHad->P();
    ptEle = track->Pt();

//    if(ptHad <2) continue;
    if(ptHad > ptEle) continue;
    if(trackHad->Eta()<-0.9 || trackHad->Eta()>0.9) continue;

    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    if(ptHad>2) DphiPtNew->Fill(ptEle,Dphi);
    if(ptHad>2 && ptHad<4) DphiPtNew1->Fill(ptEle,Dphi);
    if(ptHad>4 && ptHad<6) DphiPtNew2->Fill(ptEle,Dphi);
    if(ptHad>6 && ptHad<8) DphiPtNew3->Fill(ptEle,Dphi);
    if(ptHad>8 && ptHad<10) DphiPtNew4->Fill(ptEle,Dphi);
  }
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrelEtaFarSide(Int_t itrack, AliVTrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2,TH2F *DphiPt3,TH2F *DphiPt4)
{
  //Construct Delta Phi between electrons and hadrons for 1<eta(had)<1.6

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts2->SetRequireTPCRefit(kTRUE);
  fTrackCuts2->SetRequireITSRefit(kTRUE);
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts2->SetMaxChi2PerClusterTPC(4);
  fTrackCuts2->SetMinNClustersTPC(80);
  fTrackCuts2->SetMaxDCAToVertexZ(3.2);
  fTrackCuts2->SetMaxDCAToVertexXY(2.4);
  fTrackCuts2->SetDCAToVertex2D(kTRUE);

  for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){
    AliVParticle* VtrackHad = fVevent->GetTrack(ktracks);
    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);

    if(IsAODanalysis()) {
      AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
      if(!atrackHad->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if((!(atrackHad->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrackHad->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(atrackHad->GetTPCNcls() < 80) continue;
    }
    else{
      AliESDtrack *etrackHad = dynamic_cast<AliESDtrack*>(VtrackHad);
      if(!fTrackCuts2->AcceptTrack(etrackHad)) continue;
    }

    if(ktracks == itrack) continue; //do not select the same electron

    Double_t ptHad= -999, pHad=-999., dEdxHad = -999;
    Double_t ptEle = -999;
    Double_t phiEle = -999, phiHad = -999, Dphi = -999;
    Double_t pi = 3.14;

    dEdxHad = trackHad->GetTPCsignal();
    ptHad = trackHad->Pt();
    pHad = trackHad->P();
    ptEle = track->Pt();

    //    if(ptHad <2) continue;
    if(ptHad > ptEle) continue;
    if(trackHad->Eta()<1 || trackHad->Eta()>1.6) continue;

    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    if(ptHad>2) DphiPt->Fill(ptEle,Dphi);
    if(ptHad>2 && ptHad<4) DphiPt1->Fill(ptEle,Dphi);
    if(ptHad>4 && ptHad<6) DphiPt2->Fill(ptEle,Dphi);
    if(ptHad>6 && ptHad<8) DphiPt3->Fill(ptEle,Dphi);
    if(ptHad>8 && ptHad<10) DphiPt4->Fill(ptEle,Dphi);

  }
}
/*
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::MixedEvent(AliAODTrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2)
{

const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
Double_t zVtx;
zVtx = pVtx->GetZ();


AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality();
Double_t centvalue = fCentrality->GetCentralityPercentile("V0M");

AliEventPool* pool = fPoolMgr->GetEventPool(centvalue, zVtx); // Get the buffer associated with the current centrality and z-vtx
if (!pool)
AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centvalue, zVtx));

//  pool->PrintInfo();
if (pool->GetCurrentNEvents() >= 5) // start mixing when 5 events are in the buffer
{
Int_t nMix = pool->GetCurrentNEvents();
fNoMixedEvents->Fill(0);
fMixStat->Fill(pool->GetCurrentNEvents(),centvalue);
fMixStat1->Fill(pool->GetCurrentNEvents(),zVtx);

// cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
for (Int_t jMix=0; jMix<nMix; jMix++)  // mix with each event in the buffer
{
TObjArray* bgTracks = pool->GetEvent(jMix);
for (Int_t i=0;i<bgTracks->GetEntriesFast(); i++)
{
AliVParticle* mixtrk = (AliVParticle*) bgTracks->At(i);

Double_t mixtrkPhi = -999;
Double_t ptEle = -999;
Double_t phiEle = -999, Dphi = -999;
Double_t pi = 3.14;
Double_t ptmixtrk = -999;

ptEle = track->Pt();
ptmixtrk = mixtrk->Pt();
if(ptmixtrk > ptEle) continue;

mixtrkPhi = mixtrk->Phi();
phiEle = track->Phi();
Dphi = phiEle - mixtrkPhi;

if (Dphi > 3*pi/2)
Dphi = Dphi - 2*pi;
if (Dphi < -pi/2)
Dphi = Dphi + 2*pi;
if(ptmixtrk>2) DphiPt->Fill(ptEle,Dphi);
if(ptmixtrk>2 && ptmixtrk<4) DphiPt1->Fill(ptEle,Dphi);
if(ptmixtrk>4 && ptmixtrk<6) DphiPt2->Fill(ptEle,Dphi);
}
}

}

}
//___________________________________________
TObjArray*  AliAnalysisTaskElecHadronCorrel::CloneAndReduceTrackList()
{
// clones a track list by using AliehDPhiBasicParticle which uses much less memory (used for event mixing)

fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
fTrackCuts2->SetRequireTPCRefit(kTRUE);
fTrackCuts2->SetRequireITSRefit(kTRUE);
fTrackCuts2->SetEtaRange(-0.9,0.9);
fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
fTrackCuts2->SetMaxChi2PerClusterTPC(3.5);
fTrackCuts2->SetMinNClustersTPC(80);

TObjArray* tracksClone = new TObjArray;
tracksClone->SetOwner(kTRUE);

for(Int_t ktracks = 0; ktracks<fAOD->GetNumberOfTracks(); ktracks++){
  AliAODTrack* track = fAOD->GetTrack(ktracks);
  if (!track) {
    printf("ERROR: Could not receive track %d\n", ktracks);
    continue;
  }
  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

  //   if(ktracks == iTrack) continue;
  Double_t eta=-999,ptHad= -999, pHad=-999., phi=-999.0;
  Int_t label=-9999, id=-999;
  eta = track->Eta();
  ptHad = track->Pt();
  pHad = track->P();
  phi= track->Phi();
  label= track->GetLabel();
  id=track->GetID();

  if(track->Eta()<-0.9 || track->Eta()>0.9) continue; 
  if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue; 
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue; 
  if(ptHad <2) continue;
  //    if(!fTrackCuts2->AcceptTrack(track)) continue;

  AliVParticle* particle = (AliVParticle*) fAOD->GetTrack(ktracks);
  tracksClone->Add(new AliehDPhiBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge()));

}

return tracksClone;
}
*/
//___________________________________________
void AliAnalysisTaskElecHadronCorrel::HadronInfo(Int_t itrack)
{
  //Hadron information

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);                                           
  fTrackCuts2->SetRequireTPCRefit(kTRUE);                                                
  fTrackCuts2->SetRequireITSRefit(kTRUE);                                                
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);                                           
  fTrackCuts2->SetMaxChi2PerClusterTPC(4);                                             
  fTrackCuts2->SetMinNClustersTPC(80);                                                   
  fTrackCuts2->SetMaxDCAToVertexZ(3.2);
  fTrackCuts2->SetMaxDCAToVertexXY(2.4);
  fTrackCuts2->SetDCAToVertex2D(kTRUE);

  for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){                  
    AliVParticle* VtrackHad = fVevent->GetTrack(ktracks);                                     
    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);                            
      continue;                                                                          
    }

    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);

    if(IsAODanalysis()) {
      AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
      if(!atrackHad->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if((!(atrackHad->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrackHad->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(atrackHad->GetTPCNcls() < 80) continue; 
    }
    else{   
      AliESDtrack *etrackHad = dynamic_cast<AliESDtrack*>(VtrackHad); 
      if(!fTrackCuts2->AcceptTrack(etrackHad)) continue; 
    }

    if(ktracks == itrack) continue; //do not select the same electron
    
    Double_t ptHad= -999;
    ptHad = trackHad->Pt();
    
    if(trackHad->Eta()<-0.9 || trackHad->Eta()>0.9) continue;
   // cout << "pt had = " << ptHad <<endl;

    if(ptHad<2) continue;

    fHadronPhi->Fill(trackHad->Phi());
    fHadronPhiPt->Fill(trackHad->Phi(),ptHad);
    if (trackHad->Eta() >0 && trackHad->Eta() <0.9) fHadronPhiTPChalf->Fill(trackHad->Phi());

    fHadronPt->Fill(ptHad);
  }
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::CheckCentrality(AliVEvent* event, Bool_t &centralitypass)
{
  // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
  if (!fkCentralityMethod) AliFatal("No centrality method set! FATAL ERROR!");
  fCentrality = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethod);
 // cout << "Centrality evaluated-------------------------: " << fCentrality <<endl;

  if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
  {
    fCentralityNoPass->Fill(fCentrality);
  //  cout << "--------------Fill no pass-------------------------"<<endl;
    centralitypass = kFALSE;
  }else
  {
    fCentralityPass->Fill(fCentrality);
  //  cout << "--------------Fill pass-------------------------"<<endl;
    centralitypass = kTRUE;
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskElecHadronCorrel::SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod)
{
  // Set a centrality range ]min, max] and define the method to use for centrality selection
  fCentralityMin = CentralityMin;
  fCentralityMax = CentralityMax;
  fkCentralityMethod = CentralityMethod;
}
