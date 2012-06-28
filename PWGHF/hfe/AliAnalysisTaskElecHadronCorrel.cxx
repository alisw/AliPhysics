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
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
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
/*
#include "AliEventPoolManager.h"
#include "AliAnalysisTaskPhiCorrelations.h"
*/
#include "AliCentrality.h"
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
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

ClassImp(AliAnalysisTaskElecHadronCorrel)
//________________________________________________________________________
  AliAnalysisTaskElecHadronCorrel::AliAnalysisTaskElecHadronCorrel(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
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
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.01)	
//  ,fPoolMgr(0x0)  
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPBefHad(0)	 
  ,fTrkEovPAft(0)	
  ,fTrkEovPAftOwn(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fdEdxAftOwn(0)	 
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
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
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fTPCnsigmaAft(0)
  ,fTPCnsigmaAftOwn(0)
  ,fNCellv1(0)
  ,fClsEv1(0)
  ,fNClusv1(0)
  ,fKFParticleP(0)
  ,fKFParticleE(0)
  ,fInvmassLS1(0)       
  ,fInvmassLS2(0)       
  ,fInvmassLS3(0)       
  ,fInvmassLS4(0)       
  ,fInvmassLS5(0)       
  ,fInvmassULS1(0)
  ,fInvmassULS2(0)
  ,fInvmassULS3(0)
  ,fInvmassULS4(0)
  ,fInvmassULS5(0)
  ,fcentrality(0)     
  ,fElecPhi(0)  
  ,fElecPhiTPC(0)  
  ,fElecPhiTPCEovP(0)  
  ,fHadronPhi(0)  
  ,fTrackHFEcuts(0)
  ,fTrakPhiSPD1(0)
  ,fTrakPhiSPD2(0)
  ,fTrakPhiSPDOr(0)
  ,fTrakPhiSPDAnd(0)
  ,fTrackHFEcutsITS(0)  
/*  ,fNoMixedEvents(0)
  ,fMixStat(0)       
  ,fMixStat1(0)        
  ,fMixedIncElecDphi(0)  
  ,fMixedPhotElecDphi(0)
  ,fMixedSemiIncElecDphi(0)  
  ,fMixedDphiULSMassLow(0)  
  ,fMixedDphiLSMassLow(0)  
*/  ,fNLSminus(0)
  ,fNLSplus(0)
  ,fNULS(0)  
  ,fHadronIPxy(0)  
  ,fHadronIPz(0)  
//  ,fSparseElectron(0)  
//  ,fvalueElectron(0)   
{
  //Named constructor

  fPID = new AliHFEpid("hfePid");
//  fvalueElectron = new Double_t[8];

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
  ,fESD(0)
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
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.01)	
//  ,fPoolMgr(0x0)    
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPBefHad(0)	 
  ,fTrkEovPAft(0)	 
  ,fTrkEovPAftOwn(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fdEdxAftOwn(0)	 
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
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
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)	
  ,fTPCnsigmaAft(0)	
  ,fTPCnsigmaAftOwn(0)	
  ,fNCellv1(0)  
  ,fClsEv1(0)
  ,fNClusv1(0)
  ,fKFParticleP(0)
  ,fKFParticleE(0)
  ,fInvmassLS1(0)   
  ,fInvmassLS2(0)   
  ,fInvmassLS3(0)   
  ,fInvmassLS4(0)   
  ,fInvmassLS5(0)   
  ,fInvmassULS1(0)  
  ,fInvmassULS2(0)  
  ,fInvmassULS3(0)  
  ,fInvmassULS4(0)  
  ,fInvmassULS5(0)  
  ,fcentrality(0)     
  ,fElecPhi(0)
  ,fElecPhiTPC(0)
  ,fElecPhiTPCEovP(0)  
  ,fHadronPhi(0)
  ,fTrackHFEcuts(0)
  ,fTrakPhiSPD1(0)
  ,fTrakPhiSPD2(0)
  ,fTrakPhiSPDOr(0)
  ,fTrakPhiSPDAnd(0)
  ,fTrackHFEcutsITS(0)  
/*  ,fNoMixedEvents(0)
  ,fMixStat(0)      
  ,fMixStat1(0)     
  ,fMixedIncElecDphi(0)  
  ,fMixedPhotElecDphi(0)
  ,fMixedSemiIncElecDphi(0)
  ,fMixedDphiULSMassLow(0) 
  ,fMixedDphiLSMassLow(0)      
*/  ,fNLSminus(0)
  ,fNLSplus(0)    
  ,fNULS(0)       
  ,fHadronIPxy(0) 
  ,fHadronIPz(0)      
    //  ,fSparseElectron(0)  
    //    ,fvalueElectron(0)  
{
  //Default constructor
  fPID = new AliHFEpid("hfePid");
  //  fvalueElectron = new Double_t[8];

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
  //   delete fSparseElectron;
  //   delete []fvalueElectron;
}
//_________________________________________

void AliAnalysisTaskElecHadronCorrel::UserExec(Option_t*)
{
  //Main loop
  //Called for each event

  // create pointer to event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }

  if(!fPID->IsInitialized()){ 
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
    fPID->InitializePID(fESD->GetRunNumber());
  }

  //-------trigger selection
  UInt_t res = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if (res==0)
    return;

  //	if( (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kFastOnly) )
  //		return;

  if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kCentral))) return;

  AliCentrality *fCentrality = (AliCentrality*)fESD->GetCentrality();

  Float_t centvalue = fCentrality->GetCentralityPercentile("V0M");
  fcentrality->Fill(centvalue);    
 // cout << "cent val" << centvalue <<endl;
  if(centvalue<0 || centvalue>10) return;

 // cout << "event no : " <<fESD->GetRunNumber() <<endl;
  Int_t fNOtrks =  fESD->GetNumberOfTracks();
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();

  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();

  // Event cut
  //	if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) return;

  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(0);

  if(fNOtrks<2) return;

  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
     pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
   }

   fPID->SetPIDResponse(pidResponse);

   fCFM->SetRecEventInfo(fESD);
/*
   //Event mixing
   AliEventPool* pool = fPoolMgr->GetEventPool(centvalue, pVtxZ); // Get the buffer associated with the current centrality and z-vtx
   if (!pool)
     AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centvalue, pVtxZ));

   TObjArray* tracksClone = CloneAndReduceTrackList();
   tracksClone->SetOwner();
   pool->UpdatePool(tracksClone);
*/
   // Track loop 
   for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
     AliESDtrack* track = fESD->GetTrack(iTracks);
     if (!track) {
       printf("ERROR: Could not receive track %d\n", iTracks);
       continue;
     }

     if(track->Pt()<1) continue;

     fTrackPtBefTrkCuts->Fill(track->Pt());		

     // RecKine: ITSTPC cuts  
     if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;

     //RecKink
     if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
       if(track->GetKinkIndex(0) != 0) continue;
     } 

     // RecPrim
     if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

     // HFE cuts: TPC PID cleanup
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

     fTrackHFEcuts->Fill(track->Phi());

     //track phi distribution for diff ITS layer hit
     if(track->HasPointOnITSLayer(0)) fTrakPhiSPD1->Fill(track->Phi());
     if(track->HasPointOnITSLayer(1)) fTrakPhiSPD2->Fill(track->Phi());

     if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) fTrakPhiSPDOr->Fill(track->Phi());
     if(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1)) fTrakPhiSPDAnd->Fill(track->Phi());

     // HFEcuts: ITS layers cuts
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;

     fTrackHFEcutsITS->Fill(track->Phi());

     fTrackPtAftTrkCuts->Fill(track->Pt());		

     Double_t fClsE = -999, p = -999, fEovP=-999, pt = -999, dEdx=-999, fTPCnSigma=0;
     pt = track->Pt();
     p = track->P();
     dEdx = track->GetTPCsignal();
     fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;

     //TPC electron phi
     if(fTPCnSigma >= -2 && fTPCnSigma <= 2){
       fElecPhiTPC->Fill(track->Phi());
     }

     //eta cut (-0.7,0.7)
     if(track->Eta() < -0.7 || track->Eta() > 0.7) continue;

     // Track extrapolation to EMCAL
     Int_t fClsId = track->GetEMCALcluster();
     if(fClsId <0) continue;
     AliESDCaloCluster *cluster = fESD->GetCaloCluster(fClsId);
     if(TMath::Abs(cluster->GetTrackDx())>0.05 || TMath::Abs(cluster->GetTrackDz())>0.05) continue;    
     fdEdxBef->Fill(p,dEdx);
     fTPCnsigma->Fill(p,fTPCnSigma);

     fTrkpt->Fill(pt);
     fClsE = cluster->E();
     fEovP = fClsE/p;
     /*
        fvalueElectron[0] = pt;
        fvalueElectron[1] = p;
        fvalueElectron[2] = fTPCnSigma;
        fvalueElectron[3] = dEdx;
        fvalueElectron[4] = fEovP;
        fvalueElectron[5] = cluster->GetM20();
        fvalueElectron[6] = cluster->GetM02();
        fvalueElectron[7] = cluster->GetDispersion();

        fSparseElectron->Fill(fvalueElectron);
      */
     if(fTPCnSigma >= -2 && fTPCnSigma <= 2){
       fTrkEovPBef->Fill(pt,fEovP);
     }
     if(fTPCnSigma < -3.5)fTrkEovPBefHad->Fill(pt,fEovP);
     /*
        Int_t pidpassed = 0;
     //--- track accepted, do PID
     AliHFEpidObject hfetrack;
     hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
     hfetrack.SetRecTrack(track);
     hfetrack.SetPbPb();
     if(fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 1;

     if(pidpassed==1){
     cout << "dedx, E/p :  "<< dEdx << ", " << fEovP <<endl;
     fTrkEovPAft->Fill(pt,fEovP);
     fdEdxAft->Fill(p,dEdx);
     fTPCnsigmaAft->Fill(p,fTPCnSigma);
     }
      */

     //Electron id with TPC and E/p
     if(fTPCnSigma >= -2 && fTPCnSigma <= 2 && fEovP >= 0.8 && fEovP <=1.2) {
       fElecPhiTPCEovP->Fill(track->Phi());

       //Electron id with shower shape  
       if(cluster->GetM20()<0.2 && cluster->GetM02()< 0.5 && cluster->GetDispersion()<1){
         fElecPhi->Fill(track->Phi());
         fTrkEovPAftOwn->Fill(pt,fEovP);
         fdEdxAftOwn->Fill(p,dEdx);
         fTPCnsigmaAftOwn->Fill(p,fTPCnSigma);

         Bool_t fFlagPhotonicElec = kFALSE;
         // select photonic electron
         SelectPhotonicElectron(iTracks,track,fFlagPhotonicElec);
         //Inclusive electron-hadron correlation
         ElectronHadCorrel(iTracks, track, fInclusiveElecDphi, fInclusiveElecDphi1,fInclusiveElecDphi2,fInclusiveElecDphi3,fInclusiveElecDphi4);
         fInclusiveElecPt->Fill(pt);
//         MixedEvent(track,fMixedIncElecDphi);
         
         // photonic electron
         if(fFlagPhotonicElec){
         //Electron hadron correlation
         ElectronHadCorrel(iTracks, track, fPhotElecDphi,fPhotElecDphi1,fPhotElecDphi2,fPhotElecDphi3,fPhotElecDphi4);
         fPhotoElecPt->Fill(pt);
//         MixedEvent(track,fMixedPhotElecDphi);
         }

         // Semi inclusive electron 
         if(!fFlagPhotonicElec){
           //Electron hadron correlation
           ElectronHadCorrel(iTracks, track, fSemiIncElecDphi, fSemiIncElecDphi1,fSemiIncElecDphi2,fSemiIncElecDphi3,fSemiIncElecDphi4);
           fSemiInclElecPt->Fill(pt);
//           MixedEvent(track,fMixedSemiIncElecDphi);
         }
         
       }
     }
   }

   //EMC clusters  
   Int_t clsNo = fESD->GetNumberOfCaloClusters();
   fNClusv1->Fill(clsNo); 
   for(Int_t iclus=0; iclus<clsNo ; iclus++){ 
     AliESDCaloCluster* clus = fESD->GetCaloCluster(iclus);
     if(!clus->IsEMCAL()) continue; 
     fNCellv1->Fill(clus->GetNCells());
     fClsEv1->Fill(clus->E());  
   }


   PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::UserCreateOutputObjects()
{
  //Create histograms
  //  TGeoManager::Import("geometry.root");
  //  fGeom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

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

  fTrkpt = new TH1F("fTrkpt","track pt",1000,0,50);
  fOutputList->Add(fTrkpt);

  fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",1000,0,50);
  fOutputList->Add(fTrackPtBefTrkCuts);

  fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",1000,0,50);
  fOutputList->Add(fTrackPtAftTrkCuts);

  fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma",1000,0,50,200,-10,10);
  fOutputList->Add(fTPCnsigma);

  fTPCnsigmaAft = new TH2F("fTPCnsigmaAft", "TPC - n sigma after hfepid",1000,0,50,200,-10,10);
  fOutputList->Add(fTPCnsigmaAft);

  fTPCnsigmaAftOwn = new TH2F("fTPCnsigmaAftOwn", "TPC - n sigma after own pid",1000,0,50,200,-10,10);
  fOutputList->Add(fTPCnsigmaAftOwn);

  fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",1000,0,50,100,0,2);
  fOutputList->Add(fTrkEovPBef);

  fTrkEovPBefHad = new TH2F("fTrkEovPBefHad","track E/p for TPCnsig < 3.5",1000,0,50,100,0,2);
  fOutputList->Add(fTrkEovPBefHad);

  fTrkEovPAft = new TH2F("fTrkEovPAft","track E/p after HFE pid",1000,0,50,100,0,2);
  fOutputList->Add(fTrkEovPAft);

  fTrkEovPAftOwn = new TH2F("fTrkEovPAftOwn","track E/p after own pid",1000,0,50,100,0,2);
  fOutputList->Add(fTrkEovPAftOwn);

  fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",1000,0,50,150,0,150);
  fOutputList->Add(fdEdxBef);

  fdEdxAft = new TH2F("fdEdxAft","track dEdx vs p after HFE pid",1000,0,50,150,0,150);
  fOutputList->Add(fdEdxAft);

  fdEdxAftOwn = new TH2F("fdEdxAftOwn","track dEdx vs p own HFE pid",1000,0,50,150,0,150);
  fOutputList->Add(fdEdxAftOwn);

  fElecPhi = new TH1F("fElecPhi", "Electron phi",1000,0,6.28);
  fOutputList->Add(fElecPhi);

  fElecPhiTPC = new TH1F("fElecPhiTPC", "Electron phi after TPC cut",1000,0,6.28);
  fOutputList->Add(fElecPhiTPC);

  fElecPhiTPCEovP = new TH1F("fElecPhiTPCEovP", "Electron phi after TPC and E/p cut",1000,0,6.28);
  fOutputList->Add(fElecPhiTPCEovP);

  fHadronPhi = new TH1F("fHadronPhi", "Hadron phi",1000,0,6.28);
  fOutputList->Add(fHadronPhi);

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

  fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleLS);

  fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleULS);

  fSemiIncElecDphi = new TH2F("fSemiIncElecDphi", "Semi Inclusive elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi);

  fSemiIncElecDphi1 = new TH2F("fSemiIncElecDphi1", "Semi Inclusive elec-had Dphi correlation for 2<pt^{asso}<4",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi1);

  fSemiIncElecDphi2 = new TH2F("fSemiIncElecDphi2", "Semi Inclusive elec-had Dphi correlation for 4<pt^{asso}<6",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi2);

  fSemiIncElecDphi3 = new TH2F("fSemiIncElecDphi3", "Semi Inclusive elec-had Dphi correlation for 6<pt^{asso}<8",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi3);

  fSemiIncElecDphi4 = new TH2F("fSemiIncElecDphi4", "Semi Inclusive elec-had Dphi correlation for 8<pt^{asso}<10",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fSemiIncElecDphi4);

  fPhotElecDphi = new TH2F("fPhotElecDphi", "Photon elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi);

  fPhotElecDphi1 = new TH2F("fPhotElecDphi1", "Photon elec-had Dphi correlation for 2<pt^{asso}<4",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi1);

  fPhotElecDphi2 = new TH2F("fPhotElecDphi2", "Photon elec-had Dphi correlation for 4<pt^{asso}<6",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi2);

  fPhotElecDphi3 = new TH2F("fPhotElecDphi3", "Photon elec-had Dphi correlation for 6<pt^{asso}<8",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi3);

  fPhotElecDphi4 = new TH2F("fPhotElecDphi4", "Photon elec-had Dphi correlation for 8<pt^{asso}<10",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fPhotElecDphi4);

  fInclusiveElecDphi = new TH2F("fInclusiveElecDphi", "Inclusive elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi);

  fInclusiveElecDphi1 = new TH2F("fInclusiveElecDphi1", "Inclusive elec-had Dphi correlation for 2<pt^{asso}<4",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi1);

  fInclusiveElecDphi2 = new TH2F("fInclusiveElecDphi2", "Inclusive elec-had Dphi correlation for 4<pt^{asso}<6",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi2);

  fInclusiveElecDphi3 = new TH2F("fInclusiveElecDphi3", "Inclusive elec-had Dphi correlation for 6<pt^{asso}<8",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi3);

  fInclusiveElecDphi4 = new TH2F("fInclusiveElecDphi4", "Inclusive elec-had Dphi correlation for 8<pt^{asso}<10",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fInclusiveElecDphi4);

  fDphiULSMassLow = new TH2F("fDphiULSMassLow", "e-h Dphi ULS, mass<cut",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow);

  fDphiULSMassLow1 = new TH2F("fDphiULSMassLow1", "e-h Dphi ULS, mass<cut for 2<pt^{asso}<4",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow1);

  fDphiULSMassLow2 = new TH2F("fDphiULSMassLow2", "e-h Dphi ULS, mass<cut for 4<pt^{asso}<6",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow2);

  fDphiULSMassLow3 = new TH2F("fDphiULSMassLow3", "e-h Dphi ULS, mass<cut for 6<pt^{asso}<8",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow3);

  fDphiULSMassLow4 = new TH2F("fDphiULSMassLow4", "e-h Dphi ULS, mass<cut for 8<pt^{asso}<10",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLow4);

  fDphiLSMassLow = new TH2F("fDphiLSMassLow", "e-h Dphi LS, mass<cut",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow);

  fDphiLSMassLow1 = new TH2F("fDphiLSMassLow1", "e-h Dphi LS, mass<cut for 2<pt^{asso}<4",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow1);

  fDphiLSMassLow2 = new TH2F("fDphiLSMassLow2", "e-h Dphi LS, mass<cut for 4<pt^{asso}<6",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow2);

  fDphiLSMassLow3 = new TH2F("fDphiLSMassLow3", "e-h Dphi LS, mass<cut for 6<pt^{asso}<8",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow3);

  fDphiLSMassLow4 = new TH2F("fDphiLSMassLow4", "e-h Dphi LS, mass<cut for 8<pt^{asso}<10",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLow4);

  fDphiULSMassLowNoPartner = new TH2F("fDphiULSMassLowNoPartner", "e-h Dphi ULS with no partner, mass<mass cut,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner);

  fDphiULSMassLowNoPartner1 = new TH2F("fDphiULSMassLowNoPartner1", "e-h Dphi ULS with no partner, mass<mass cut for 2<pt^{asso}<4,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner1);

  fDphiULSMassLowNoPartner2 = new TH2F("fDphiULSMassLowNoPartner2", "e-h Dphi ULS with no partner, mass<mass cut for 4<pt^{asso}<6,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner2);

  fDphiULSMassLowNoPartner3 = new TH2F("fDphiULSMassLowNoPartner3", "e-h Dphi ULS with no partner, mass<mass cut for 6<pt^{asso}<8,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner3);

  fDphiULSMassLowNoPartner4 = new TH2F("fDphiULSMassLowNoPartner4", "e-h Dphi ULS with no partner, mass<mass cut for 8<pt^{asso}<10,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiULSMassLowNoPartner4);

  fDphiLSMassLowNoPartner = new TH2F("fDphiLSMassLowNoPartner", "e-h Dphi LS with no partner, mass<mass cut",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner);

  fDphiLSMassLowNoPartner1 = new TH2F("fDphiLSMassLowNoPartner1", "e-h Dphi LS with no partner, mass<mass cut for 2<pt^{asso}<4,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner1);                                            

  fDphiLSMassLowNoPartner2 = new TH2F("fDphiLSMassLowNoPartner2", "e-h Dphi LS with no partner, mass<mass cut for 4<pt^{asso}<6,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner2);

  fDphiLSMassLowNoPartner3 = new TH2F("fDphiLSMassLowNoPartner3", "e-h Dphi LS with no partner, mass<mass cut for 6<pt^{asso}<8,",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fDphiLSMassLowNoPartner3);

  fDphiLSMassLowNoPartner4 = new TH2F("fDphiLSMassLowNoPartner4", "e-h Dphi LS with no partner, mass<mass cut for 8<pt^{asso}<10,",200,0,20,100,-1.57,4.71);
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

  fNCellv1 = new TH1F("fNCellv1","Ncell in clus (v1); NCell; count",100,0,100) ;
  fOutputList->Add(fNCellv1);

  fClsEv1 = new TH1F("fClsEv1", "Clus E(v1); Cls E; count",1000,0,100); 
  fOutputList->Add(fClsEv1); 

  fNClusv1 = new TH1F("fNClusv1","Nclus in event (v1); NClus; count",500,0,500) ; 
  fOutputList->Add(fNClusv1);

  fKFParticleP = new TH1F("fKFParticleP","KFparticle rec P; P(GeV/c)",1000,0,50);
  fOutputList->Add(fKFParticleP);

  fKFParticleE = new TH1F("fKFParticleE", "KfParticle rec E; E; count",1000,0,100); 
  fOutputList->Add(fKFParticleE);

  fInvmassLS1 = new TH1F("fInvmassLS1", "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS1);

  fInvmassULS1 = new TH1F("fInvmassULS1", "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS1);

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
/*
  fNoMixedEvents = new TH1F("fNoMixedEvents","",1,0,1) ;
  fOutputList->Add(fNoMixedEvents);

  fMixStat = new TH2F("fMixStat","no of events in pool  vs multiplicity;Nevent in pool;N hadrons",200,0,200,500,0,500);
  fOutputList->Add(fMixStat);                                                             

  fMixStat1 = new TH2F("fMixStat1","no of events in pool  vs zvtx;Nevents in pool;zvtx",200,0,200,200,-11,11);
  fOutputList->Add(fMixStat1);

  fMixedIncElecDphi = new TH2F("fMixedIncElecDphi", "Mixed event - Inclusive elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedIncElecDphi);

  fMixedSemiIncElecDphi = new TH2F("fMixedSemiIncElecDphi", "Mixed event - Semi Inclusive elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedSemiIncElecDphi);

  fMixedPhotElecDphi = new TH2F("fMixedPhotElecDphi", "Mixed event - Photo elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedPhotElecDphi);

  fMixedDphiULSMassLow = new TH2F("fMixedDphiULSMassLow", "Mixed event - ULS mass < cut elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiULSMassLow);

  fMixedDphiLSMassLow = new TH2F("fMixedDphiLSMassLow", "Mixed event - LS mass < cut elec-had Dphi correlation",200,0,20,100,-1.57,4.71);
  fOutputList->Add(fMixedDphiLSMassLow);
*/
  fNLSminus = new TH1F("fNLSminus","No of LS negative pairs (e-,e-) ",1000,-0.5,999.5);
  fOutputList->Add(fNLSminus);

  fNLSplus = new TH1F("fNLSplus","No of LS positive pairs (e+,e+)",1000,-0.5,999.5);
  fOutputList->Add(fNLSplus);

  fNULS = new TH1F("fNULS","No of ULS pairs (e+,e-)",1000,-0.5,999.5);
  fOutputList->Add(fNULS);

  fHadronIPxy = new TH1F("fHadronIPxy", "hadron impact paramter XY",1000,-5,5);
  fOutputList->Add(fHadronIPxy);

  fHadronIPz = new TH1F("fHadronIPz", "hadron impact paramter Z",1000,-20,20);
  fOutputList->Add(fHadronIPz);

  /*
     Int_t binsv1[8]={1000,1000,200,150,100,100,100,100}; //pt, p, TPCnsig, dEdx, E/p, M20, M02, dispersion 
     Double_t xminv1[8]={0,0,-10,0,0,0,0,0};
     Double_t xmaxv1[8]={50,50,10,150,2,2,2,2};
     fSparseElectron = new THnSparseD ("Electron","Electron",8,binsv1,xminv1,xmaxv1);
     fOutputList->Add(fSparseElectron);
   */
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
void AliAnalysisTaskElecHadronCorrel::SelectPhotonicElectron(Int_t itrack, AliESDtrack *track, Bool_t &fFlagPhotonicElec)
{
  //Identify non-heavy flavour electrons using Invariant mass method

  fTrackCuts1->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts1->SetRequireTPCRefit(kTRUE);
  fTrackCuts1->SetEtaRange(-0.9,0.9);
  fTrackCuts1->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts1->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts1->SetMinNClustersTPC(80);

  //  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();

  Bool_t flagPhotonicElec = kFALSE;
  Int_t NLS_plus=0, NLS_minus=0, NULS=0;

  for(Int_t jTracks = itrack+1; jTracks<fESD->GetNumberOfTracks(); jTracks++){
    AliESDtrack* trackAsso = fESD->GetTrack(jTracks);
    if (!trackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }

    Double_t dEdxAsso = -999., ptAsso=-999., openingAngle = -999.,nsigma=-999.0;
    Double_t mass=-999., width = -999;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;

    dEdxAsso = trackAsso->GetTPCsignal();
    nsigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();

    if(ptAsso <0.3) continue;
    if(!fTrackCuts1->AcceptTrack(trackAsso)) continue;
//    if(dEdxAsso <70 || dEdxAsso>100) continue; //11a pass1
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

    openingAngle = ge1.GetAngle(ge2);
    if(fFlagLS) fOpeningAngleLS->Fill(openingAngle);
    if(fFlagULS) fOpeningAngleULS->Fill(openingAngle);

 //   if(openingAngle > fOpeningAngleCut) continue;

    recg.GetMass(mass,width);
  
    if(fFlagLS){
      if(track->Charge() > 0 ) NLS_plus++;
      if(track->Charge() < 0 ) NLS_minus++;
    }
    if(fFlagULS) NULS++;

    if(fFlagLS) {
      fInvmassLS1->Fill(mass);
      if(track->Pt()> 1) fInvmassLS2->Fill(mass);
      if(track->Pt()>2) fInvmassLS3->Fill(mass);
      if(track->Pt()>3) fInvmassLS4->Fill(mass);
      if(track->Pt()>4) fInvmassLS5->Fill(mass);
    }
    if(fFlagULS) {
      fInvmassULS1->Fill(mass);
      if(track->Pt() >1) fInvmassULS2->Fill(mass);
      if(track->Pt() >2) fInvmassULS3->Fill(mass);
      if(track->Pt() >3) fInvmassULS4->Fill(mass);
      if(track->Pt() >4) fInvmassULS5->Fill(mass);
    }

    if(mass<fInvmassCut){
      if(fFlagULS)
      {
        ElectronHadCorrel(itrack,track,fDphiULSMassLow, fDphiULSMassLow1,fDphiULSMassLow2,fDphiULSMassLow3,fDphiULSMassLow4);
        fULSElecPt->Fill(track->Pt());
//        MixedEvent(track,fMixedDphiULSMassLow);
      }
      if(fFlagLS)
      {
        ElectronHadCorrel(itrack,track,fDphiLSMassLow,fDphiLSMassLow1,fDphiLSMassLow2,fDphiLSMassLow3,fDphiLSMassLow4);
        fLSElecPt->Fill(track->Pt());
//        MixedEvent(track,fMixedDphiLSMassLow);
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

  fNLSminus->Fill(NLS_minus);
  fNLSplus->Fill(NLS_plus);
  fNULS->Fill(NULS);

}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrel(Int_t itrack, AliESDtrack *track, TH2F *DphiPt, TH2F *DphiPt1,TH2F *DphiPt2,TH2F *DphiPt3,TH2F *DphiPt4)
{
  //Construct Delta Phi between electrons and hadrons

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts2->SetRequireTPCRefit(kTRUE);
  fTrackCuts2->SetRequireITSRefit(kTRUE);
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts2->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts2->SetMinNClustersTPC(80);

  for(Int_t ktracks = 0; ktracks<fESD->GetNumberOfTracks(); ktracks++){
    AliESDtrack* trackHad = fESD->GetTrack(ktracks);
    if (!trackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
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
    if(!fTrackCuts2->AcceptTrack(trackHad)) continue;
  
    Float_t IPxy=-999.0, IPz=-999.0;
    trackHad->GetImpactParameters(IPxy,IPz);
   //cout << "xy,z: " << IPxy << ", " << IPz <<endl;
    fHadronIPxy->Fill(IPxy);
    fHadronIPz->Fill(IPz);

    fHadronPhi->Fill(trackHad->Phi());
    
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
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrelNoPartner(Int_t itrack,Int_t jtrack, AliESDtrack *track, TH2F *DphiPtNew, TH2F *DphiPtNew1,TH2F *DphiPtNew2,TH2F *DphiPtNew3,TH2F *DphiPtNew4)
{
  //Construct Delta Phi between electrons and hadrons for electrons from invariant mass calculation excluding associated track

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts2->SetRequireTPCRefit(kTRUE);
  fTrackCuts2->SetRequireITSRefit(kTRUE);
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts2->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts2->SetMinNClustersTPC(80);

  for(Int_t ktracks = 0; ktracks<fESD->GetNumberOfTracks(); ktracks++){
    AliESDtrack* trackHad = fESD->GetTrack(ktracks);
    if (!trackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
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

    if(ptHad <2) continue;
    if(ptHad > ptEle) continue;
    if(!fTrackCuts2->AcceptTrack(trackHad)) continue;

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
/*
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::MixedEvent(AliESDtrack *track, TH2F *DphiPt)
{

  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  Double_t zVtx;
  zVtx = pVtx->GetZ();


  AliCentrality *fCentrality = (AliCentrality*)fESD->GetCentrality();
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
    for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++)  // mix with each event in the buffer
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
        if(ptmixtrk > ptEle) continue;

        mixtrkPhi = mixtrk->Phi();
        phiEle = track->Phi();
        Dphi = phiEle - mixtrkPhi;

        if (Dphi > 3*pi/2)
          Dphi = Dphi - 2*pi;
        if (Dphi < -pi/2)
          Dphi = Dphi + 2*pi;
        DphiPt->Fill(ptEle,Dphi);
      }
    }

  }

}
//___________________________________________
TObjArray*  AliAnalysisTaskElecHadronCorrel::CloneAndReduceTrackList()
{
  // clones a track list by using AliDPhiBasicParticle which uses much less memory (used for event mixing)

  fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts2->SetRequireTPCRefit(kTRUE);
  fTrackCuts2->SetRequireITSRefit(kTRUE);
  fTrackCuts2->SetEtaRange(-0.9,0.9);
  fTrackCuts2->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts2->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts2->SetMinNClustersTPC(80);

  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for(Int_t ktracks = 0; ktracks<fESD->GetNumberOfTracks(); ktracks++){
    AliESDtrack* track = fESD->GetTrack(ktracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }

    //   if(ktracks == iTrack) continue;
    Double_t ptHad= -999, pHad=-999.;
    ptHad = track->Pt();
    pHad = track->P();
    
    if(ptHad <2) continue;
    if(!fTrackCuts2->AcceptTrack(track)) continue;

    AliVParticle* particle = (AliVParticle*) fESD->GetTrack(ktracks);
    tracksClone->Add(new AliDPhiBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge()));
  }

  return tracksClone;
}
*/
