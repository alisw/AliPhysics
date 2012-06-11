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
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPBefHad(0)	 
  ,fTrkEovPAft(0)	
  ,fTrkEovPAftOwn(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fdEdxAftOwn(0)	 
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fSemiIncElecDphi(0) 	
  ,fPhotElecDphi(0)  	
  ,fInclusiveElecDphi(0)  	
  ,fDphiMassHigh(0)		
  ,fDphiULSMassLow(0)	
  ,fDphiLSMassLow(0)
  ,fDphiULSMassLowNoPartner(0)   
  ,fDphiLSMassLowNoPartner(0)
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
  ,fInvmassULS1(0)
  ,fcentrality(0)     
  ,fElecPhi(0)  
  ,fHadronPhi(0)  
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
  ,fNoEvents(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPBefHad(0)	 
  ,fTrkEovPAft(0)	 
  ,fTrkEovPAftOwn(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fdEdxAftOwn(0)	 
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fSemiIncElecDphi(0) 	
  ,fPhotElecDphi(0)  	
  ,fInclusiveElecDphi(0)  	
  ,fDphiMassHigh(0)		
  ,fDphiULSMassLow(0)	
  ,fDphiLSMassLow(0)
  ,fDphiULSMassLowNoPartner(0)   
  ,fDphiLSMassLowNoPartner(0)
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
  ,fInvmassULS1(0)  
  ,fcentrality(0)     
  ,fElecPhi(0)
  ,fHadronPhi(0)
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
   cout << "cent val" << centvalue <<endl;
   if(centvalue<0 || centvalue>10) return;

   cout << "event no : " <<fESD->GetRunNumber() <<endl;
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

   // Track loop 
   for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
     AliESDtrack* track = fESD->GetTrack(iTracks);
     if (!track) {
       printf("ERROR: Could not receive track %d\n", iTracks);
       continue;
     }
    
     //eta cut (-0.7,0.7)
     if(track->Eta() < -0.7 || track->Eta() > 0.7) continue;

     fTrackPtBefTrkCuts->Fill(track->Pt());		
     // RecKine: ITSTPC cuts  
     if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;

     //RecKink
     if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
       if(track->GetKinkIndex(0) != 0) continue;
     } 

     // RecPrim
     if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;

     // HFEcuts: ITS layers cuts
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;

     // HFE cuts: TPC PID cleanup
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

     fTrackPtAftTrkCuts->Fill(track->Pt());		
       
     Double_t fClsE = -999, p = -999, fEovP=-999, pt = -999, dEdx=-999, fTPCnSigma=0;
     // Track extrapolation
     Int_t fClsId = track->GetEMCALcluster();
     if(fClsId <0) continue;
     AliESDCaloCluster *cluster = fESD->GetCaloCluster(fClsId);
     if(TMath::Abs(cluster->GetTrackDx())>0.05 || TMath::Abs(cluster->GetTrackDz())>0.05) continue;    

     pt = track->Pt();
     if(pt<1) continue;
     fTrkpt->Fill(pt);
     fClsE = cluster->E();
     cout << "cluster E = " << fClsE <<endl;
     p = track->P();
     dEdx = track->GetTPCsignal();
     fEovP = fClsE/p;
     fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
     fdEdxBef->Fill(p,dEdx);
     fTPCnsigma->Fill(p,fTPCnSigma);
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
     if(fTPCnSigma >= -2 && fTPCnSigma <= 2)fTrkEovPBef->Fill(pt,fEovP);
     if(fTPCnSigma < -4.5)fTrkEovPBefHad->Fill(pt,fEovP);
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
     if(fTPCnSigma >= -2 && fTPCnSigma <= 2 && fEovP >= 0.9 && fEovP <=1.3) {
       if(cluster->GetM20()<0.3 && cluster->GetM02()< 0.7 && cluster->GetDispersion()<1){
         fElecPhi->Fill(track->Phi());
         fTrkEovPAftOwn->Fill(pt,fEovP);
         fdEdxAftOwn->Fill(p,dEdx);
         fTPCnsigmaAftOwn->Fill(p,fTPCnSigma);

         Bool_t fFlagPhotonicElec = kFALSE;
         // select photonic electron
         SelectPhotonicElectron(iTracks,track,fFlagPhotonicElec);
         //Inclusive electron-hadron correlation
         ElectronHadCorrel(iTracks, track, fInclusiveElecDphi);
         fInclusiveElecPt->Fill(pt);

         // photonic electron
         if(fFlagPhotonicElec){
           //Electron hadron correlation
           ElectronHadCorrel(iTracks, track, fPhotElecDphi);
           fPhotoElecPt->Fill(pt);
         }

         // Semi inclusive electron 
         if(!fFlagPhotonicElec){
           //Electron hadron correlation
           ElectronHadCorrel(iTracks, track, fSemiIncElecDphi);
           fSemiInclElecPt->Fill(pt);
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

  fHadronPhi = new TH1F("fHadronPhi", "Hadron phi",1000,0,6.28);
  fOutputList->Add(fHadronPhi);

//  fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) if mass cal is correct; mass(GeV/c^2); counts;", 1000,0,1.0);
//  fOutputList->Add(fInvmassLS);

//  fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e); mass(GeV/c^2) if mass cal is correct; counts;", 1000,0,1.0);
//  fOutputList->Add(fInvmassULS);

  fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleLS);

  fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleULS);

  fSemiIncElecDphi = new TH2F("fSemiIncElecDphi", "Semi Inclusive elec-had Dphi correlation",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fSemiIncElecDphi);

  fPhotElecDphi = new TH2F("fPhotElecDphi", "Photon elec-had Dphi correlation",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fPhotElecDphi);

  fInclusiveElecDphi = new TH2F("fInclusiveElecDphi", "Inclusive elec-had Dphi correlation",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fInclusiveElecDphi);

  fDphiMassHigh = new TH2F("fDphiMassHigh", "e-h Dphi LS+ULS, mass>0.01",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fDphiMassHigh);

  fDphiULSMassLow = new TH2F("fDphiULSMassLow", "e-h Dphi ULS, mass<0.01",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fDphiULSMassLow);

  fDphiLSMassLow = new TH2F("fDphiLSMassLow", "e-h Dphi LS, mass<0.01",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fDphiLSMassLow);

  fDphiULSMassLowNoPartner = new TH2F("fDphiULSMassLowNoPartner", "e-h Dphi ULS with no partner, mass<mass cut,",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fDphiULSMassLowNoPartner);

  fDphiLSMassLowNoPartner = new TH2F("fDphiLSMassLowNoPartner", "e-h Dphi LS with no partner, mass<mass cut",200,0,20,100,-1.6,4.75);
  fOutputList->Add(fDphiLSMassLowNoPartner);

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

  for(Int_t jTracks = itrack+1; jTracks<fESD->GetNumberOfTracks(); jTracks++){
    AliESDtrack* trackAsso = fESD->GetTrack(jTracks);
    if (!trackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }

    Double_t dEdxAsso = -999., ptAsso=-999., openingAngle = -999.;
    Double_t mass=-999., width = -999;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;

    dEdxAsso = trackAsso->GetTPCsignal();
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();

    if(ptAsso <0.3) continue;
    if(!fTrackCuts1->AcceptTrack(trackAsso)) continue;
    if(dEdxAsso <70 || dEdxAsso>100) continue; //11a pass1

    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;

    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;

    AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
    AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);
/*
    Double_t recP2=-999.0, recP=-999.0, recE=-999.0, m2=-999.0, m=-999.0;
    recP2= (recg.GetPx()*recg.GetPx() + recg.GetPy()*recg.GetPy() + recg.GetPz()*recg.GetPz());
    recP = TMath::Sqrt(recP2);
    recE = recg.GetE();
    fKFParticleP->Fill(recP);
    fKFParticleE->Fill(recE);
    m2 = (recg.GetE()*recg.GetE() - recg.GetPx()*recg.GetPx() - recg.GetPy()*recg.GetPy() - recg.GetPz()*recg.GetPz());
    m = TMath::Sqrt(m2);
*/
    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

    openingAngle = ge1.GetAngle(ge2);
    if(fFlagLS) fOpeningAngleLS->Fill(openingAngle);
    if(fFlagULS) fOpeningAngleULS->Fill(openingAngle);

    if(openingAngle > fOpeningAngleCut) continue;

    Int_t MassCorrect;
    MassCorrect = recg.GetMass(mass,width);
    cout << "mass = " << mass <<endl;

    if(fFlagLS) fInvmassLS1->Fill(mass);
    if(fFlagULS) fInvmassULS1->Fill(mass);

    //    if(MassCorrect==1){
//    if(fFlagLS) fInvmassLS->Fill(m);
//    if(fFlagULS) fInvmassULS->Fill(m);

    if(mass>fInvmassCut){
      ElectronHadCorrel(itrack,track,fDphiMassHigh);
    }
    if(mass<fInvmassCut){
      if(fFlagULS)
      {
        ElectronHadCorrel(itrack,track,fDphiULSMassLow);
        fULSElecPt->Fill(track->Pt());
      }
      if(fFlagLS)
      {
        ElectronHadCorrel(itrack,track,fDphiLSMassLow);
        fLSElecPt->Fill(track->Pt());
      }
      if(fFlagLS) ElectronHadCorrelNoPartner(itrack,jTracks,track,fDphiLSMassLowNoPartner);
      if(fFlagULS) ElectronHadCorrelNoPartner(itrack,jTracks,track,fDphiULSMassLowNoPartner);
    }

    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
    //   }

  }
  fFlagPhotonicElec = flagPhotonicElec;

}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrel(Int_t itrack, AliESDtrack *track, TH2F *DphiPt)
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

    if(ptHad <0.3) continue;
    if(!fTrackCuts2->AcceptTrack(trackHad)) continue;

    fHadronPhi->Fill(trackHad->Phi());

    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    ptEle = track->Pt();

    DphiPt->Fill(ptEle,Dphi);
  }
}
//_________________________________________
void AliAnalysisTaskElecHadronCorrel::ElectronHadCorrelNoPartner(Int_t itrack,Int_t jtrack, AliESDtrack *track, TH2F *DphiPtNew)
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

    if(ptHad <0.3) continue;
    if(!fTrackCuts2->AcceptTrack(trackHad)) continue;

    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    ptEle = track->Pt();

    DphiPtNew->Fill(ptEle,Dphi);
  }
}
